#Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman
#PS 08.11.18 - added imp to load local python files common_formats, skippy.
#PS 08.16.18 - edited print statements, comments for clarity

#PS 08.18.18 - suppressing ignorable warnings
# https://github.com/numpy/numpy/pull/432
import warnings
warnings.simplefilter("ignore")

import imp
imp.load_source('common_formats','../analysis/common_formats.py')
imp.load_source('skippy','../analysis/skippy.py')

import string
import numpy
import sys
import os
import simplejson as json
import pysam



import common_formats
import skippy

SETTINGS = json.load(open('analysis.json','r'))
if not os.path.exists('bedgraph'):
	os.mkdir('bedgraph')

if not os.path.exists('curvegraph'):
	os.mkdir('curvegraph')

processed_dir = SETTINGS['processed_dir']
if not os.path.exists(processed_dir):
	os.mkdir(processed_dir)


POLSKI = {'positive':'+', 'negative':'-'}
boolPOLSKI = {'+':False, '-': True}


per_chromosome_TLS_data = {}
per_chromosome_features_masked_out = {}
unwanted_genes = set([])

common_formats.import_gff_file(SETTINGS['GFF_FILE'], per_chromosome_TLS_data, per_chromosome_features_masked_out, unwanted_genes)
# you calculate this by taking the log of the expected maximum number of items
MAX_ORF_SKIPLIST = SETTINGS['MAX_ORF_SKIPLIST']

#import skippy
orf_lookup = {}
position_to_orf_trees = {}
jj = 0
for achr, polar in per_chromosome_TLS_data.iterkeys():
    world = achr, boolPOLSKI[polar]
    me = skippy.IntervalSkiplist(MAX_ORF_SKIPLIST)
    # print world
    for orfname, furthest_TSS in per_chromosome_TLS_data[achr, polar].iteritems():
        theTIS = per_chromosome_features_masked_out[(achr, polar)][orfname][0]
        left, right = sorted([furthest_TSS, theTIS])
        if left!=right:
            orf_lookup[orfname] = jj
            me.insert_interval(orfname, left, right)
            jj += 1
    position_to_orf_trees[world] = me


just_T_or_N = frozenset(['T','N'])
just_A_or_N = frozenset(['A','N'])
crap = numpy.zeros((2,),dtype='int')

# we throw away STAR alignments that soft clip presumably meaningful components of the read
def read_quality_control(bb):
	if len(bb.cigar)==1:
		return True # assume it was just a match
	# check 5' end to see if we are only throwing away Ts
	if bb.cigar[0][0]==4: # soft clipping
		if frozenset(bb.query_sequence[:bb.cigar[0][1]]) > just_T_or_N:
			crap[0] += 1
			return False
	# check 3' end to see if we are only throwing away As
	if bb.cigar[-1][0]==4:
		if frozenset(bb.query_sequence[-bb.cigar[-1][1]:]) > just_A_or_N:
			crap[1] += 1
			return False
	return True

samfile_name = sys.argv[1]

# def quant_samfile(samfile_name):
maxseqlen = SETTINGS['rpf_maxseqlen']
minseqlen = SETTINGS['rpf_minseqlen']
sizeacc = {}
samfile = pysam.AlignmentFile(samfile_name, 'rb')
orf_quants = numpy.zeros((len(orf_lookup),))
print 'working on', samfile_name
chrdata = []
curvedata = []
permuted_chrdata = []
num_permutations = len(SETTINGS['permutation_seeds'])
for csize in samfile.lengths:
	# positive, then negative
	chrdata.append((numpy.zeros((csize,)), numpy.zeros((csize,))))
	curvedata.append((numpy.zeros((csize,), dtype='uint'), numpy.zeros((csize,), dtype='uint')))
	permuted_chrdata.append((numpy.zeros((csize,num_permutations)), numpy.zeros((csize,num_permutations))))

rand_sources = []
for i in SETTINGS['permutation_seeds']:
	rand_sources.append(numpy.random.RandomState(i))

for achrom, not_strand in sorted(position_to_orf_trees.iterkeys()):
	hitplace = achrom, not_strand
	if achrom not in samfile.references:
		continue
	for anorf, (aleft, aright) in position_to_orf_trees[achrom, not_strand].lookup.iteritems():
		if numpy.any(numpy.sum(permuted_chrdata[samfile.references.index(achrom)][not_strand][aleft:aright],axis=0))!=0:
			print 'error in sum of permuted data'
			import pdb; pdb.set_trace()
		for bb in samfile.fetch(reference=achrom, start=aleft-maxseqlen, end=aright+maxseqlen):
			if not bb.is_unmapped and bb.is_reverse==not_strand and 'N' not in bb.cigarstring and read_quality_control(bb):
				# this should work with soft clipping from STAR now
				seq_len = bb.reference_length
				if maxseqlen >= seq_len >= minseqlen:
					# all_reads += 1
					# reporting [a,b)
					if not bb.is_reverse:
						if seq_len >= 28:
							a,b = 12, seq_len-15
						else:
							a,b = seq_len-16, 12
					else:
						if seq_len >= 28:
							a,b = 15, seq_len-12
						else:
							a,b = seq_len-12, 16
					pre_a = a
					pre_b = b
					a += bb.reference_start
					b += bb.reference_start
					# ignore reads that don't land their p-site schmear entirely inside the UTR
					if ( aleft <= a <= aright ) and ( aleft <= b <= aright ):
						schmear = 1.0/(b-a)
						chrdata[bb.rname][bb.is_reverse][a:b] += schmear
						curvedata[bb.rname][bb.is_reverse][bb.reference_start:bb.reference_start+seq_len] += 1
						if seq_len not in sizeacc:
							sizeacc[seq_len] = 0
						sizeacc[seq_len] += 1
						for i in range(num_permutations):
							# this is to ensure that the shmear of the simulated read ALSO lands entirely inside the UTR
							if aright-pre_b+1 <= aleft-pre_a:
								print 'wtf', 1/0
							offset = rand_sources[i].randint(aleft-pre_a, aright-pre_b+1)
							a = pre_a + offset
							b = pre_b + offset
							permuted_chrdata[bb.rname][bb.is_reverse][a:b,i] += schmear
		# we've already adjusted for the reference_start!
		# for hitset, (begin, end) in position_to_orf_trees[hitplace].intervals_containing_node(a):
		orf_quants[orf_lookup[anorf]] = numpy.sum(chrdata[samfile.references.index(achrom)][not_strand][aleft:aright])
		# check to make sure permutations are close
		permuted_sums_here = numpy.sum(permuted_chrdata[samfile.references.index(achrom)][not_strand][aleft:aright],axis=0)
		try:
			numpy.testing.assert_allclose(permuted_sums_here, orf_quants[orf_lookup[anorf]])
			# print "I'm ok", permuted_sums_here
		except AssertionError:
			print anorf
			import pdb; pdb.set_trace()

all_reads = 0 #samfile.mapped
samfile.reset()
for bb in samfile.fetch(until_eof=True):
	if not bb.is_unmapped and 'N' not in bb.cigarstring and read_quality_control(bb):
		seq_len = bb.reference_length
		if maxseqlen >= seq_len >= minseqlen:
			all_reads += 1


print '...saving for', samfile_name
sn_short = samfile_name.split('/')[-1] 
for x,y in zip(samfile.references, chrdata):
	numpy.save(processed_dir+sn_short+'-'+x+'-positive_data.npy', y[0])
	# FLIP IT!!!
	numpy.save(processed_dir+sn_short+'-'+x+'-negative_data.npy', y[1][::-1])

for ii,x in enumerate(samfile.references):
	for jj in range(num_permutations):
		numpy.save(processed_dir+'p'+str(jj)+'-'+sn_short+'-'+x+'-positive_data.npy', permuted_chrdata[ii][0][:,jj])
		numpy.save(processed_dir+'p'+str(jj)+'-'+sn_short+'-'+x+'-negative_data.npy', permuted_chrdata[ii][1][:,jj][::-1])

fout_pos = open('bedgraph/psite-'+sn_short+'-watson.bedGraph','w')
fout_neg = open('bedgraph/psite-'+sn_short+'-crick.bedGraph','w')
for x,y in zip(samfile.references, chrdata):
	pos, = numpy.where(numpy.diff(y[0]) != 0)
	pos = numpy.concatenate(([0],pos+1,[len(y[0])]))
	rle = [(a,b,y[0][a]) for (a,b) in zip(pos[:-1],pos[1:])]
	for ss,se,vv in rle:
		fout_pos.write("\t".join([x, str(ss), str(se), str(vv)])+"\n")
	pos, = numpy.where(numpy.diff(y[1]) != 0)
	pos = numpy.concatenate(([0],pos+1,[len(y[1])]))
	rle = [(a,b,y[1][a]) for (a,b) in zip(pos[:-1],pos[1:])]
	for ss,se,vv in rle:
		fout_neg.write("\t".join([x, str(ss), str(se), str(vv)])+"\n")

fout_pos.close()
fout_neg.close()

for i in range(num_permutations):
	fout_pos = open('bedgraph/psite-'+sn_short+'-null'+str(i)+'-watson.bedGraph','w')
	fout_neg = open('bedgraph/psite-'+sn_short+'-null'+str(i)+'-crick.bedGraph','w')
	for x,yy in zip(samfile.references, permuted_chrdata):
		# unpack the permuted_chrdata:
		# permuted_chrdata[bb.rname][bb.is_reverse][a:b,i]
		y = yy[0][:,i], yy[1][:,i]
		pos, = numpy.where(numpy.diff(y[0]) != 0)
		pos = numpy.concatenate(([0],pos+1,[len(y[0])]))
		rle = [(a,b,y[0][a]) for (a,b) in zip(pos[:-1],pos[1:])]
		for ss,se,vv in rle:
			fout_pos.write("\t".join([x, str(ss), str(se), str(vv)])+"\n")
		pos, = numpy.where(numpy.diff(y[1]) != 0)
		pos = numpy.concatenate(([0],pos+1,[len(y[1])]))
		rle = [(a,b,y[1][a]) for (a,b) in zip(pos[:-1],pos[1:])]
		for ss,se,vv in rle:
			fout_neg.write("\t".join([x, str(ss), str(se), str(vv)])+"\n")

	fout_pos.close()
	fout_neg.close()



fout_pos = open('curvegraph/curve-'+sn_short+'-watson.bedGraph','w')
fout_neg = open('curvegraph/curve-'+sn_short+'-crick.bedGraph','w')
for x,y in zip(samfile.references, curvedata):
	pos, = numpy.where(numpy.diff(y[0]) != 0)
	pos = numpy.concatenate(([0],pos+1,[len(y[0])]))
	rle = [(a,b,y[0][a]) for (a,b) in zip(pos[:-1],pos[1:])]
	for ss,se,vv in rle:
		fout_pos.write("\t".join([x, str(ss), str(se), str(vv)])+"\n")
	pos, = numpy.where(numpy.diff(y[1]) != 0)
	pos = numpy.concatenate(([0],pos+1,[len(y[1])]))
	rle = [(a,b,y[1][a]) for (a,b) in zip(pos[:-1],pos[1:])]
	for ss,se,vv in rle:
		fout_neg.write("\t".join([x, str(ss), str(se), str(vv)])+"\n")

fout_pos.close()
fout_neg.close()

fout = open(samfile_name.split('/')[-1]+'.dump','w')
fout.write(' '.join([str(x) for x in crap])+"\n")
for x,y in sizeacc.iteritems():
	fout.write(str(x)+' '+str(y)+"\n")

fout.close()

fout = open(processed_dir+'/'+samfile_name.split('/')[-1]+'-quantified.csv','w')
fout.write('ALL_READS '+str(all_reads)+"\n")
for anorf,ii in orf_lookup.iteritems():
	fout.write(anorf+' '+str(orf_quants[ii])+"\n")
fout.close()
