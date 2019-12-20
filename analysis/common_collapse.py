#Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman
import skippy
from collections import namedtuple
import os
import numpy

MIN_RPF_CONSIDERED = 3
MAX_ORF_SKIPLIST = 4

RESULTS = namedtuple('RESULTS', 'tbl feats dfeats countdata freqdata NSAMPLES')

def apply_zscore(old_results, mf, sf):
	return RESULTS(old_results.tbl, (old_results.feats-mf)/sf, old_results.dfeats, old_results.countdata, old_results.freqdata, old_results.NSAMPLES)

def partition_results(aresult, xs):
	left = RESULTS([aresult.tbl[x] for x in xs], aresult.feats[xs], aresult.dfeats[xs], aresult.countdata[xs], aresult.freqdata[xs], aresult.NSAMPLES)
	notxs = [x for x in xrange(len(aresult.tbl)) if x not in set(xs)]
	right = RESULTS([aresult.tbl[x] for x in notxs], aresult.feats[notxs], aresult.dfeats[notxs], aresult.countdata[notxs], aresult.freqdata[notxs], aresult.NSAMPLES)
	return left, right


def within_sample_determination(PREFIX, feat1, dir1, feat2, dir2, dfeat__protection_sum):
	atbl = []
	for line in open(PREFIX+'.table'):
		z = tuple(line.rstrip().split())
		atbl.append((z[0], z[1], int(z[2]), int(z[3]), z[4]))

	feats = numpy.load(PREFIX+'-feats.npy')
	dfeats = numpy.load(PREFIX+'-design-feats.npy')

	pp = PREFIX.split('/')
	fout = open('/'.join(pp[:-1])+'/pre_collapse-'+pp[-1]+'.bed','w')
	pre_collapse = {}
	for ii,x in enumerate(atbl):
		fout.write("\t".join([x[0], str(x[2]), str(x[3]), x[1]+'.'+str(ii),'0',x[-1]])+"\n")
		if dfeats[ii,dfeat__protection_sum]>=MIN_RPF_CONSIDERED:
			if '+' in x[-1]:
				pre_collapse.setdefault((x[0],x[1],x[3],x[-1]),[]).append(ii)
			else:
				pre_collapse.setdefault((x[0],x[1],x[2],x[-1]),[]).append(ii)
	fout.close()

	# collapse all those sharing a stop codon by selecting the ones with the best in-phase of the max power frequency
	# and then breaking ties (if it happens) by the best PWM score
	nested_collapsed = {}
	for x,ys in pre_collapse.iteritems():
		if dir1:
			best_candidate_val = numpy.max(feats[ys,feat1])
		else:
			best_candidate_val = numpy.min(feats[ys,feat1])
		ppx = [ys[z] for z in numpy.where(feats[ys,feat1]==best_candidate_val)[0]]
		if dir2:
			nested_collapsed.setdefault((x[0],x[1],x[-1]),[]).append(ppx[numpy.argmax(feats[ppx,feat2])])
		else:
			nested_collapsed.setdefault((x[0],x[1],x[-1]),[]).append(ppx[numpy.argmin(feats[ppx,feat2])])


	# to resolve overlapping scores, simply greedily choose in order of best PWM scores, and ignore any conflicts
	acc = []
	for x,rawys in nested_collapsed.iteritems():
		# if 'YJL080C' in x:
		# 	import pdb; pdb.set_trace()

		rrys = [(atbl[y][2], atbl[y][3],y) for y in rawys]
		# ys = sorted(rrys, key=lambda x: feats[x[2],feat__break_ties], reverse=True)
		ys = sorted(rrys, key=lambda q: feats[q[2],feat2], reverse=dir2)
		# t = IntervalTree()
		t = skippy.IntervalSkiplist(3)
		for y in ys:
			if not t.is_interval_contained(y[0], y[1]):
				t.insert_interval(y[2], y[0], y[1])
		acc.extend(t.lookup.keys())

	return [atbl[x] for x in acc], feats[acc], dfeats[acc]

def combine_results(xs):
	return RESULTS([item for sublist in [x.tbl for x in xs] for item in sublist], \
		numpy.vstack([x.feats for x in xs]), \
		numpy.vstack([x.dfeats for x in xs]), \
		numpy.hstack([x.countdata for x in xs]), \
		numpy.hstack([x.freqdata for x in xs]), \
		xs[0].NSAMPLES)

def single_sample_determination(SETTINGS, wsample, feat1, dir1, feat2, dir2, dfeat__protection_sum):
	suffix = '-'.join([str(x) for x in [feat1, dir1, feat2, dir2]])
	tbl, feats, dfeats = within_sample_determination(SETTINGS['PREFIX']+'/'+wsample, feat1, dir1, feat2, dir2, dfeat__protection_sum)
	return RESULTS(tbl, feats, dfeats, numpy.ones((feats.shape[0],)), numpy.ones((feats.shape[0],)), 1)

def single_samples(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum):
	return [single_sample_determination(SETTINGS, w, feat1, dir1, feat2, dir2, dfeat__protection_sum) for w in SETTINGS['SAMPLES']]

# this will average them
def grouped_samples_determination(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, group_name, wsamples, extra_suffix='', force_replay=False):
	suffix = '-'+group_name+'-'.join([str(x) for x in [feat1, dir1, feat2, dir2]])+'-'+extra_suffix
	calc_dir = SETTINGS['PREFIX']+'/calculated'+extra_suffix+'/'
	if force_replay or not os.path.exists(calc_dir+'raw_averaged_dfeats'+suffix+'.npy'):
		tbls, feats, dfeats = zip(*[within_sample_determination(calc_dir+w, feat1, dir1, feat2, dir2, dfeat__protection_sum) for w in wsamples])
		tbls = [item for sublist in tbls for item in sublist]
		feats = numpy.vstack(feats)
		dfeats = numpy.vstack(dfeats)
		uorf_lookup = {}
		for ii,x in enumerate(tbls):
			uorf_lookup.setdefault(x,[]).append(ii)
		averaged_tbl = []
		averaged_feats = numpy.zeros((len(uorf_lookup),feats.shape[1]))
		averaged_dfeats = numpy.zeros((len(uorf_lookup),dfeats.shape[1]))
		averaged_countdata = numpy.zeros((len(uorf_lookup),), dtype='int')
		for ii,(x,ys) in enumerate(uorf_lookup.iteritems()):
			averaged_tbl.append(x)
			averaged_feats[ii,:] = numpy.mean(feats[ys,:],axis=0)
			averaged_dfeats[ii,:] = numpy.mean(dfeats[ys,:],axis=0)
			averaged_countdata[ii] = len(ys)
		numpy.save(calc_dir+'raw_averaged_feats'+suffix+'.npy', averaged_feats)
		numpy.save(calc_dir+'raw_averaged_dfeats'+suffix+'.npy', averaged_dfeats)
		numpy.save(calc_dir+'raw_averaged_countdata'+suffix+'.npy', averaged_countdata)
		with open(calc_dir+'raw_averaged_tbl'+suffix+'.tbl','w') as fout:
			for x in averaged_tbl:
				fout.write(' '.join([str(q) for q in x])+"\n")
		print group_name, len(averaged_tbl)
	else:
		averaged_feats = numpy.load(calc_dir+'raw_averaged_feats'+suffix+'.npy')
		averaged_dfeats = numpy.load(calc_dir+'raw_averaged_dfeats'+suffix+'.npy')
		averaged_countdata = numpy.load(calc_dir+'raw_averaged_countdata'+suffix+'.npy')
		averaged_tbl = []
		for line in open(calc_dir+'raw_averaged_tbl'+suffix+'.tbl'):
			z = tuple(line.rstrip().split())
			averaged_tbl.append((z[0], z[1], int(z[2]), int(z[3]), z[4]))
		print group_name, len(averaged_tbl)
	NSAMPLES = int(numpy.max(averaged_countdata))
	averaged_freqdata = averaged_countdata/float(NSAMPLES)
	results = RESULTS(averaged_tbl, averaged_feats, averaged_dfeats, averaged_countdata, averaged_freqdata, NSAMPLES)
	return results


def across_samples_determination(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, null_sample=None, force_replay=False):
	if null_sample is not None:
		extra_suffix = 'p'+str(null_sample)+'-'
	else:
		extra_suffix = ''
	acc = []
	all_samples = set([x['name'] for x in SETTINGS['SAMPLES']])
	if 'SAMPLE_GROUPS' in SETTINGS:
		for gg in SETTINGS['SAMPLE_GROUPS']:
			wsamples = gg['wsamples']
			# the nulls will get collapsed by name elsewhere
			acc.append((gg['groupname'], grouped_samples_determination(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, gg['groupname'], wsamples, extra_suffix, force_replay)))
			all_samples -= set(wsamples)
	for asamp in all_samples:
		acc.append((asamp, grouped_samples_determination(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, asamp, [asamp], extra_suffix)))

	return acc

def load_data_from_settings(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, force_replay=False):
	obsdata = across_samples_determination(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, None, force_replay=force_replay)
	nulldata = [item for sublist in [across_samples_determination(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, null_sample=i, force_replay=force_replay) for i in range(len(SETTINGS['permutation_seeds']))] for item in sublist]
	return obsdata, nulldata

def load_nonnull_data_from_settings(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, force_replay=False):
	obsdata = across_samples_determination(SETTINGS, feat1, dir1, feat2, dir2, dfeat__protection_sum, None, force_replay=force_replay)
	return obsdata



	# suffix = '-'.join([str(x) for x in [feat1, dir1, feat2, dir2]])+'-'+extra_suffix
	# calc_dir = SETTINGS['PREFIX']+'/calculated'+extra_suffix+'/'
	# if not os.path.exists(calc_dir+'raw_averaged_dfeats'+suffix+'.npy'):
	# 	tbls, feats, dfeats = zip(*[within_sample_determination(calc_dir+w['name'], feat1, dir1, feat2, dir2, dfeat__protection_sum, feat__break_ties) for w in SETTINGS['SAMPLES']])
	# 	tbls = [item for sublist in tbls for item in sublist]
	# 	feats = numpy.vstack(feats)
	# 	dfeats = numpy.vstack(dfeats)
	# 	uorf_lookup = {}
	# 	for ii,x in enumerate(tbls):
	# 		uorf_lookup.setdefault(x,[]).append(ii)
	# 	averaged_tbl = []
	# 	averaged_feats = numpy.zeros((len(uorf_lookup),feats.shape[1]))
	# 	averaged_dfeats = numpy.zeros((len(uorf_lookup),dfeats.shape[1]))
	# 	averaged_countdata = numpy.zeros((len(uorf_lookup),))
	# 	for ii,(x,ys) in enumerate(uorf_lookup.iteritems()):
	# 		averaged_tbl.append(x)
	# 		averaged_feats[ii,:] = numpy.mean(feats[ys,:],axis=0)
	# 		averaged_dfeats[ii,:] = numpy.mean(dfeats[ys,:],axis=0)
	# 		averaged_countdata[ii] = len(ys)
	# 	numpy.save(calc_dir+'raw_averaged_feats'+suffix+'.npy', averaged_feats)
	# 	numpy.save(calc_dir+'raw_averaged_dfeats'+suffix+'.npy', averaged_dfeats)
	# 	numpy.save(calc_dir+'raw_averaged_countdata'+suffix+'.npy', averaged_countdata)
	# 	with open(calc_dir+'raw_averaged_tbl'+suffix+'.tbl','w') as fout:
	# 		for x in averaged_tbl:
	# 			fout.write(' '.join([str(q) for q in x])+"\n")
	# else:
	# 	averaged_feats = numpy.load(calc_dir+'raw_averaged_feats'+suffix+'.npy')
	# 	averaged_dfeats = numpy.load(calc_dir+'raw_averaged_dfeats'+suffix+'.npy')
	# 	averaged_countdata = numpy.load(calc_dir+'raw_averaged_countdata'+suffix+'.npy')
	# 	averaged_tbl = []
	# 	for line in open(calc_dir+'raw_averaged_tbl'+suffix+'.tbl'):
	# 		z = tuple(line.rstrip().split())
	# 		averaged_tbl.append((z[0], z[1], int(z[2]), int(z[3]), z[4]))
	# NSAMPLES = int(numpy.max(averaged_countdata))
	# averaged_freqdata = averaged_countdata/NSAMPLES
	# results = RESULTS(averaged_tbl, averaged_feats, averaged_dfeats, averaged_countdata, averaged_freqdata, NSAMPLES)
	# return results

def partition_for_isolates(results, avoid={}):
	achromo = {}
	for ii,x in enumerate(results.tbl):
		achromo.setdefault((x[0],x[4]),skippy.IntervalSkiplist(MAX_ORF_SKIPLIST)).insert_interval(ii, x[2], x[3])
	# find the ones with no overlaps
	isolates = set([])
	for x,ys in achromo.iteritems():
		isolates |= ys.isolated_intervals()
		# for y in ys.lookup:
		# 	if len(ys.search(y))==1 and results.tbl[y.data] not in avoid:
		# 		isolates.append(y.data)
	#
	isolates = [x for x in isolates if x not in avoid]
	#
	iso_tbl = [results.tbl[x] for x in isolates]
	iso_feats = results.feats[isolates,:]
	iso_dfeats = results.dfeats[isolates,:]
	iso_freqdata = results.freqdata[isolates]
	iso_countdata = results.countdata[isolates]
	#
	nonisolated = [x for x in range(len(results.tbl)) if x not in set(isolates) and x not in avoid]
	#
	noniso_tbl = [results.tbl[x] for x in nonisolated]
	noniso_feats = results.feats[nonisolated,:]
	noniso_dfeats = results.dfeats[nonisolated,:]
	noniso_freqdata = results.freqdata[nonisolated]
	noniso_countdata = results.countdata[nonisolated]
	iso_results = RESULTS(iso_tbl, iso_feats, iso_dfeats, iso_countdata, iso_freqdata, results.NSAMPLES)
	noniso_results = RESULTS(noniso_tbl, noniso_feats, noniso_dfeats, noniso_countdata, noniso_freqdata, results.NSAMPLES)
	return iso_results, noniso_results

