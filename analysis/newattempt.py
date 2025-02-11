
#PS 08.11.18 - added imp to load local python files common_formats.
#PS 08.16.18 - edited print statements, comments for clarity
#PS 08.16.18: edited to handle log of zero errors on 'entropy_of_power' 

#PS 08.18.18 - suppressing ignorable warnings
# https://github.com/numpy/numpy/pull/432
import warnings
warnings.simplefilter("ignore")

import imp
imp.load_source('common_formats','../analysis/common_formats.py')

import numpy
import pylab
import scipy
import scipy.signal
import sys
import os
import string
import simplejson as json
import common_formats

SETTINGS = json.load(open('analysis.json','r'))

# this is the number of nucleotides we expect to have schmeared by the ribosome overhang
# even if ribosomes are mostly occupying at a P-site sitting on an AUG 
# via inspection of YEL009C
OVERHANG_BY_RIBOSOME = SETTINGS['OVERHANG_BY_RIBOSOME']

# A size limit to avoid absurdly small UTRs
TOO_SMALL_UTR = SETTINGS['TOO_SMALL_UTR']

# a lower bound bound on RPFs
MIN_RPF_CONSIDERED = SETTINGS['MIN_RPF_CONSIDERED']

PERMUTE_5pUTR = SETTINGS['PERMUTE_5pUTR']
PERMUTE_FRAME = SETTINGS['PERMUTE_FRAME']

in_phase_tolerance = SETTINGS['in_phase_tolerance']

predefined_kozak_sequences = SETTINGS['example_kozak_sequences']

FORCE_ONLY_LABELED = SETTINGS['FORCE_ONLY_LABELED']

POLSKI = {'positive':'+', 'negative':'-'}

# remove this
min_codon_length_filter = 0 # Rumi, Pieter had 3

# change
WTYPE = sys.argv[3]
SAMPLE_PREFIX = sys.argv[1]
results_dir = 'predictions.'+WTYPE+SAMPLE_PREFIX+'/'
predict_on = sys.argv[2]

if not os.path.exists(results_dir):
	os.mkdir(results_dir)

annotated = SETTINGS['labeled_data']


per_chromosome_TLS_data = {}
per_chromosome_features_masked_out = {}
unwanted_genes = set([])

common_formats.import_gff_file(SETTINGS['GFF_FILE'], per_chromosome_TLS_data, per_chromosome_features_masked_out, unwanted_genes)

per_chromosome_labeled = {}
for f in annotated:
	per_chromosome_labeled = common_formats.import_bed_smpl(f, per_chromosome_labeled)

base_index = {}
base_index['.'] = [0,0,0,0]
base_index['A'] = [1,0,0,0]
base_index['T'] = [0,1,0,0]
base_index['G'] = [0,0,1,0]
base_index['C'] = [0,0,0,1]
base_index['N'] = [1,1,1,1] 
base_index['R'] = [1,0,1,0] 
base_index['Y'] = [0,1,0,1] 
base_index['W'] = [1,1,0,0] 
base_index['S'] = [0,0,1,1] 
base_index['M'] = [1,0,0,1] 
base_index['K'] = [0,1,1,0] 
base_index['H'] = [1,1,0,1] 
base_index['B'] = [0,1,1,1] 
base_index['V'] = [1,0,1,1] 
base_index['D'] = [1,1,1,0] 

vals_to_base = dict([((x[0]<<3)+(x[1]<<2)+(x[2]<<1)+x[3],y) for (y,x) in base_index.items()])
vals_to_base[0] = '.'

def nicely_letters(xs):
	return ''.join([vals_to_base[x] for x in (xs[:,0]<<3)+(xs[:,1]<<2)+(xs[:,2]<<1)+xs[:,3]])

def lower_if_mismatch(aa,bb):
	acc = []
	for x in zip(aa,bb):
		if x[0]==x[1]:
			acc.append(x[1])
		else:
			acc.append(x[1].lower())
	return ''.join(acc)

start_windows = {}
for q in predefined_kozak_sequences:
	start_windows.setdefault(q[:-4]+q[-1:],set([])).add(q[-4:-1])

for the_chromosome in SETTINGS['chromosomes']:
	for polarity in ['positive', 'negative']:
		positions = per_chromosome_labeled.setdefault((the_chromosome, POLSKI[polarity]),[])

acc = []
for q in [numpy.vstack([base_index[z] for z in x]) for x in start_windows.keys()]:
	b = numpy.zeros((10,4), dtype='int64')
	b[:-4] = q[:-1]
	b[-1] = q[-1]
	acc.append(b)

start_windows = acc


THE_ATG = numpy.vstack(base_index[z] for z in 'ATG')

# Taken from:
# Miyasaka, H. (1999), The positive relationship between codon usage bias and translation initiation AUG context in Saccharomyces cerevisiae. Yeast, 15:633-637
PWM_FROM_STARTS = numpy.zeros((10,4))
# A C G T, so we need to convert to ATGC
PWM_FROM_STARTS[0,:] = [1.00, 0.167, 0.500, 0.833]
PWM_FROM_STARTS[1,:] = [1.00, 0.909, 0.091, 0.727]
PWM_FROM_STARTS[2,:] = [1.00, 0.333, 0.111, 0.222]
PWM_FROM_STARTS[3,:] = [1.00, 0.021, 0.250, 0.021]
PWM_FROM_STARTS[4,:] = [1.00, 0.217, 0.043, 0.043]
PWM_FROM_STARTS[5,:] = [1.00, 0.080, 0.040, 0.080]
PWM_FROM_STARTS[9,:] = [0.357, 0.143, 1.000, 0.643]

PWM_FROM_STARTS = numpy.hstack([PWM_FROM_STARTS[:,[0,3,2,1]]])
PWM_FROM_STARTS[6:9,:] = 1 #THE_ATG

# might want to emit sufficient statistics on start windows
print '... there are', len(start_windows), 'unique seed Kozak sequences'

# stops are TAG, TAA, TGA. find them and mark the nearest upstream in-frame stop codon for each base
stop_codons = [numpy.vstack([base_index[x] for x in zulu]) for zulu in ['TAG', 'TAA', 'TGA']]

def swizzle_seq(TSS, TIS, frame, seqdata):
	if TSS < TIS:
		right = TIS-frame
		ncodons = (right-TSS)/3
		left = right-ncodons*3
		wcodons = numpy.arange(ncodons)
		numpy.random.shuffle(wcodons)
		seqdata[left:right,:] = seqdata[left+numpy.hstack([numpy.arange(3)+3*x for x in wcodons]),:]

orf_characterization_mem = {}


# make predictions
the_chromosome = predict_on
this_chromosome_predictions = []
labeled_evidence = []
for polarity in ['positive', 'negative']:
	print 'forming predictions for chromosome', the_chromosome, '('+polarity+')'
	seqdata = numpy.load(SETTINGS['GENOME_CACHE']+the_chromosome+'-'+polarity+'_strand.npy').astype('int')
	pwm_scores = numpy.sum(numpy.vstack([numpy.convolve(seqdata[:,z],PWM_FROM_STARTS[:,z][::-1]) for z in range(0,4)]),axis=0)
	# # make the p-value lookup table for this chromosome
	# unique_pwm_scores = numpy.unique(pwm_scores)
	# vv = numpy.cumsum(numpy.bincount(numpy.digitize(pwm_scores, unique_pwm_scores))[::-1])[::-1][1:]/float(pwm_scores.shape[0])
	# pwm_pvalues = dict(zip(unique_pwm_scores,vv))
	sprefix = [WTYPE+x['RPF'].split('/')[-1] for x in SETTINGS['SAMPLES'] if x['name']==SAMPLE_PREFIX][0]
	ribdata_orig = numpy.load(SETTINGS['processed_dir']+sprefix+'-'+the_chromosome+'-'+polarity+'_data.npy')
	if ( (the_chromosome, POLSKI[polarity]) in per_chromosome_labeled ):
		if polarity == 'positive':
			labeled_positions = per_chromosome_labeled[(the_chromosome, POLSKI[polarity])]
			#print 'Change A'
			#print labeled_positions
		else:
			labeled_positions = [(seqdata.shape[0]-a[1], seqdata.shape[0]-a[0], a[2]) for a in per_chromosome_labeled[(the_chromosome, POLSKI[polarity])]]
	else:
		labeled_positions = []

	labeled_starts = [x[0] for x in labeled_positions]
	print 'change D'
	print labeled_positions
	print '...', len(labeled_starts), 'annotations available'

	# find possible starts
	seqdata_sans_orfs = numpy.zeros(seqdata.shape,dtype=seqdata.dtype)
	ribdata = numpy.zeros(ribdata_orig.shape,dtype=ribdata_orig.dtype)
	the_chromosome_mask = numpy.zeros(seqdata.shape[0], dtype='bool')
	orf_annotations = []
	if (the_chromosome, POLSKI[polarity]) in per_chromosome_TLS_data:
		# first, mask IN the 5'UTRs we have TSS/TIS data for
		if (the_chromosome, POLSKI[polarity]) in per_chromosome_TLS_data: 
			for orfname, furthest_TLS in per_chromosome_TLS_data[(the_chromosome, POLSKI[polarity])].iteritems():
				region_to_keep = furthest_TLS, per_chromosome_features_masked_out[(the_chromosome, POLSKI[polarity])][orfname][0]
				if abs(region_to_keep[0]-region_to_keep[1]) < TOO_SMALL_UTR:
					print '*** WARNING', orfname, "5'UTR is too small"
					continue
				if polarity == 'positive':
					a,b = region_to_keep
				else:
					a = seqdata.shape[0]-region_to_keep[0]
					b = seqdata.shape[0]-region_to_keep[1]
				# ensure that we have a CANONICAL start codon to get around misannotated cases
				if 'ATG' == nicely_letters(seqdata[b:b+3]):
					the_chromosome_mask[a:b-OVERHANG_BY_RIBOSOME] = True
					# add this to our available annotations
					# because we are going to next mask out other genome features, this ends up being an overapproximation
					# of what we can actually find hits to, but we're not going to optimize this yet
					orf_annotations.append((orfname, a, b))
					if PERMUTE_5pUTR:
						swizzle_seq(a, b, PERMUTE_FRAME, seqdata)
				else:
					print 'Change B'
					print nicely_letters(seqdata[b:b+3])
					the_chromosome_mask[a:b-OVERHANG_BY_RIBOSOME] = True
					print '*** PROBABLY MISANNOTATED TIS', orfname, region_to_keep, nicely_letters(seqdata[b:b+3]), nicely_letters(seqdata[b-2:b+1])
					# # this hack is to get around bayanus problems
					# if 'ATG' == nicely_letters(seqdata[b-2:b+1]):
					#     the_chromosome_mask[a:b-OVERHANG_BY_RIBOSOME] = True
					#     # add this to our available annotations
					#     # because we are going to next mask out other genome features, this ends up being an overapproximation
					#     # of what we can actually find hits to, but we're not going to optimize this yet
					#     orf_annotations.append((orfname, a, b-2))
					#     if PERMUTE_5pUTR:
					#         swizzle_seq(a, b, PERMUTE_FRAME, seqdata)
					# else:
					#     print '*** PROBABLY MISANNOTATED TIS', orfname, region_to_keep, nicely_letters(seqdata[b:b+3]), nicely_letters(seqdata[b-2:b+1])
		print 'chromosome mask-in', numpy.sum(the_chromosome_mask)
		# now mask OUT all the features of any main ORF coding region that might overlap
		if (the_chromosome, POLSKI[polarity]) in per_chromosome_features_masked_out: 
			for orfname, region_to_discard in per_chromosome_features_masked_out[(the_chromosome, POLSKI[polarity])].iteritems():
				if polarity == 'positive':
					a,b = region_to_discard
				else:
					a = seqdata.shape[0]-region_to_discard[0]
					b = seqdata.shape[0]-region_to_discard[1]
				if numpy.any(the_chromosome_mask[a:b-OVERHANG_BY_RIBOSOME]):
					print '** based on '+orfname+" we are masking out an overlapping 5'UTR"
				the_chromosome_mask[a:b-OVERHANG_BY_RIBOSOME] = False
		# now apply the mask directly for the RPF data
		ribdata[the_chromosome_mask] = ribdata_orig[the_chromosome_mask]
		seqdata_sans_orfs[the_chromosome_mask,:] = seqdata[the_chromosome_mask,:]


	matchvs = numpy.sum(numpy.vstack([[numpy.sum(numpy.vstack([numpy.convolve(seqdata[:,z], sw[:,z][::-1]) for z in range(4)]),axis=0)==3] for sw in stop_codons]),axis=0)
	# from -2 to +1 of the match
	stop_positions = numpy.where(matchvs)[0]-2
	stop_lookup = numpy.zeros((seqdata.shape[0],),dtype='uint')
	for relframe in range(0,3):
		# grab the stops that are in frame relframe, starting from...relframe
		stops_in_relative_frame = numpy.hstack([relframe,stop_positions[numpy.mod(stop_positions,3)==relframe]])
		for j in range(1,stops_in_relative_frame.shape[0]):
			stop_lookup[stops_in_relative_frame[j-1]:stops_in_relative_frame[j]:3] = stops_in_relative_frame[j]

	# look for ATGs in TLSs (strong) and one-offs (weak)
	matchvs = numpy.sum(numpy.vstack([numpy.convolve(seqdata_sans_orfs[:,z], THE_ATG[:,z][::-1], mode='same') for z in range(4)]),axis=0)
	# the -1 aligns it to the '1' position of ATG (the A)
	strong_possibles = numpy.where(matchvs==3)[0]-1
	weak_possibles = numpy.where(matchvs==2)[0]-1

	# LIMIT TO NUG
	possible_strong_kozaks = numpy.asarray([numpy.prod(numpy.sum(PWM_FROM_STARTS*seqdata[x-6:x+4],axis=1)) for x in strong_possibles])
	possible_weak_kozaks = numpy.asarray([numpy.prod(numpy.sum(PWM_FROM_STARTS*seqdata[x-6:x+4],axis=1)) for x in weak_possibles])
	#
	# idx = numpy.where(possible_weak_kozaks>=SETTINGS['weak_kozak_cutoff'])[0]
	idx = numpy.arange(possible_weak_kozaks.shape[0])
	weak_possibles = weak_possibles[idx]
	possible_weak_kozaks = possible_weak_kozaks[idx]
	# find starts
	putative_starts = sorted(zip(strong_possibles, possible_strong_kozaks)+zip(weak_possibles, possible_weak_kozaks), key=lambda x:x[0])
	print '...', len(putative_starts), 'possible starts for uORFs'

	# find possible in-frame stops
	unmeasured_orfs = [] 
	for ss,v in putative_starts:
		closest_in_frame_stop = stop_lookup[ss]
		if closest_in_frame_stop > ss:
			# the +3 is to include the stop codon in the uORF sequence
			unmeasured_orfs.append((int(ss), int(closest_in_frame_stop)+3, v))

	# mainorf, morfs, morfe = [x for x in orf_annotations if 'YJL183W' in x[0]][0] 
	def characterize_morf_things(mainorf, morfs, morfe):
		if mainorf not in orf_characterization_mem:
			# print mainorf, morfs, morfe
			summ = numpy.sum(ribdata[morfs:morfe])
			if summ <= MIN_RPF_CONSIDERED:
				return
			zulu_q = numpy.fft.rfft(ribdata[morfs:morfe]/summ)
			zulu_power_spectra = numpy.abs(zulu_q)
			zulu_freq = numpy.argmax(zulu_power_spectra[1:])+1
			zulu_phase = numpy.abs(numpy.angle(zulu_q[zulu_freq]))
			orf_characterization_mem[mainorf] = {}
			orf_characterization_mem[mainorf]['average_power'] = summ
			orf_characterization_mem[mainorf]['maxpower'] = zulu_power_spectra[zulu_freq]
			orf_characterization_mem[mainorf]['maxpower_angle'] = zulu_phase
			orf_characterization_mem[mainorf]['TLSlen'] = morfe-morfs                
			# identify a feature that is either an N-terminal extension or a misannotated start
			possible_extensions = [x for x in stop_positions if ( morfs <= x <= morfe ) and numpy.mod(x,3)==numpy.mod(morfe,3)]
			if len(possible_extensions)>0:
				# first T in the stop codon (should be in-frame to the stop)
				possible_extension = numpy.max(possible_extensions)
			else:
				# choose the position near the TSS that is in frame to the TIS
				possible_extension = numpy.floor(morfs/3)*3+numpy.mod(morfe,3)
			orf_characterization_mem[mainorf]['possible_extension'] = possible_extension


	def characterize_puorf(starti, endi, ribdata, mymap, seqdata, polarity, mainorf, morfs, morfe):
		characterize_morf_things(mainorf, morfs, morfe)
		if mainorf not in orf_characterization_mem:
			return False
		mymap['morf_5putr_region_power'] = orf_characterization_mem[mainorf]['average_power']
		summ = mymap['morf_5putr_region_power']
		within_summ = numpy.sum(ribdata[starti:endi])
		if within_summ <= MIN_RPF_CONSIDERED:
			return False
		puorf_fft = numpy.fft.rfft(ribdata[starti:endi]/summ)
		power_spectra = numpy.abs(puorf_fft)
		phase_spectra = numpy.abs(numpy.angle(puorf_fft))
		mymap['max_power_freq'] = numpy.argmax(power_spectra[1:])+1
		mymap['within_power_of_max_power_freq'] = power_spectra[mymap['max_power_freq']]
		mymap['within_phase_of_max_power_freq'] = abs(phase_spectra[mymap['max_power_freq']])
		if mymap['within_phase_of_max_power_freq'] < 10**-3:
			mymap['within_phase_of_max_power_freq'] = 0
		freq_closest_to_three = 3
		mymap['within_phase_of_in_frame'] = abs(phase_spectra[freq_closest_to_three])
		if mymap['within_phase_of_in_frame'] < 10**-3:
			mymap['within_phase_of_in_frame'] = 0
		mymap['within_power_of_in_frame'] = power_spectra[freq_closest_to_three]
		mymap['median_power'] = numpy.median(power_spectra[1:])
		mymap['median_phase'] = numpy.median(phase_spectra[1:])
		fraction_power = power_spectra/numpy.sum(power_spectra)
		mymap['weighted_avg_phase'] = numpy.nansum(fraction_power[1:]*phase_spectra[1:])
  
		# PS 08.16.18: edited to handle log of zero errors
		if fraction_power.all() != 0:
			mymap['entropy_of_power'] = -numpy.nansum(fraction_power*numpy.log(fraction_power))
		else:
			fraction_power = 1
			mymap['entropy_of_power'] = -numpy.nansum(fraction_power*numpy.log(fraction_power))
   
		# Brar et. al. concept...modifed so that it doesn't have divide by zero issues.
		mymap['relative_start_magnitude'] = (ribdata[starti]-numpy.mean(ribdata[(starti-2*3):starti]))
		# mymap['relmag'] = power_spectra[mymap['max_power_freq']]/(power_spectra[0]+1)
		mymap['protection_sum'] = within_summ
		if polarity == 'negative': 
			ss_start = seqdata.shape[0] - endi
			ss_end = seqdata.shape[0] - starti
		else:
			ss_start = starti
			ss_end = endi
		mymap['ss_start'] = ss_start
		mymap['ss_end'] = ss_end
		mymap['ss_sequence'] = nicely_letters(seqdata[starti-6:starti+4])
		mymap['pwm_score'] = pwm_scores[starti+3] # aligns to end of the "kozak" window
		# mymap['pwm_pvalue'] = pwm_pvalues[pwm_scores[starti+3]]
		mymap['polarity'] = POLSKI[polarity]
		mymap['uorf_sequence'] = nicely_letters(seqdata[starti:endi])
		mymap['morf_framediff'] = float(numpy.mod(starti,3)!=numpy.mod(morfe,3))
		mymap['morf_name'] = mainorf
		mymap['morf_dist_TIS'] = endi-morfe
		mymap['morf_dist_TSS'] = morfs-starti
		# this now measures the fraction of the uorf that doesn't overlap the main orf in the prefix
		mymap['morf_nonoverlap'] = (min(morfe,endi)-starti)/(endi-starti)
		# if mainorf not in orf_characterization_mem:
		#     zulu_q = numpy.fft.rfft(ribdata[morfs:morfe])
		#     zulu_power_spectra = numpy.abs(zulu_q)**2
		#     zulu_freq = numpy.argmax(zulu_power_spectra[1:])+1
		#     zulu_phase = numpy.angle(zulu_q[zulu_freq])
		#     orf_characterization_mem[mainorf] = (zulu_power_spectra[0],zulu_power_spectra[zulu_freq],numpy.angle(zulu_q[zulu_freq]))
		mymap['morf_len_TLS'] =  orf_characterization_mem[mainorf]['TLSlen'] 
		mymap['morf_5putr_peak_power'] = orf_characterization_mem[mainorf]['maxpower']
		mymap['morf_5putr_peak_phase'] = orf_characterization_mem[mainorf]['maxpower_angle']
		# wasted calculations
		# mymap['morf_5putr_region_power_nterm'] = orf_characterization_mem[mainorf]['nterm_average_power']
		# mymap['morf_5putr_peak_power_nterm'] = orf_characterization_mem[mainorf]['nterm_maxpower']
		# mymap['morf_5putr_peak_phase_nterm'] = orf_characterization_mem[mainorf]['nterm_maxpower_angle']
		mymap['morf_5putr_sum_upstream'] = numpy.sum(ribdata[morfs:starti])
		mymap['morf_5putr_sum_downstream'] = numpy.sum(ribdata[endi:morfe])
		# mymap['morf_nterm_region_framediff'] = numpy.abs(numpy.mod(starti,3)-numpy.mod(orf_characterization_mem[mainorf]['possible_extension'],3))
		mymap['morf_nterm_region_ends_before'] = int(orf_characterization_mem[mainorf]['possible_extension']>=endi)
		mymap['morf_nterm_region_start_after'] = int(orf_characterization_mem[mainorf]['possible_extension']<=starti)
		# 20160315 remove so we have an external variable to measure
		# mymap['rna_start_folding_energy'] = _RNA.fold(nicely_letters(seqdata[starti-28:starti+28]))[1]
		# print nicely_letters(seqdata[endi-3:endi])
		# if nicely_letters(seqdata[endi-3:endi]) not in ['TAG', 'TAA', 'TGA']:
		#     print 1/0
		# if 'YDR028C' in mainorf and ss_end==500910 or ss_end==500907:
		# 	import pdb
		# 	pdb.set_trace()
		# 	print [mymap[q] for q in ['ss_start', 'ss_end', 'within_phase_of_max_power_freq', 'within_phase_of_in_frame']]
		return True


	def gimme_the_main_orf_5pUTR(starti, endi):
		# a puorf MUST start after a TSS, and can only extend past the corresponding TIS 
		# iff it is in a different frame than the TIS
		starts_after_TSS = [x[1]<=starti for x in orf_annotations]
		starts_before_TIS = [starti<=x[2] for x in orf_annotations]
		end_before_TIS_or_different_frame = [(not (numpy.mod(x[2],3)==numpy.mod(endi,3))) or (endi < x[2]) for x in orf_annotations]
		mycases = [numpy.all(x) for x in zip(starts_after_TSS, starts_before_TIS, end_before_TIS_or_different_frame)]
		legal_mpointers = numpy.where(mycases)[0]
		if legal_mpointers.shape[0]==0:
			return None, None, None
		mpointer = legal_mpointers[numpy.argmin([starti-orf_annotations[x][1] for x in legal_mpointers])]
		# figure out what the ORF is for this uORF in the clumsiest fashion possible
		mainorf, morfs, morfe = orf_annotations[mpointer]
		return mainorf, morfs, morfe

	# compute scores for putative uORFs 
	for starti, endi, kozy in unmeasured_orfs:
		if (endi-starti)/3 >= min_codon_length_filter:
			mainorf, morfs, morfe = gimme_the_main_orf_5pUTR(starti, endi)
			# if mainorf is not None and 'YJL080C' in mainorf:
			# 	import pdb; pdb.set_trace()

			if mainorf is None:
				continue
			mymap = {}
			ok = characterize_puorf(starti, endi, ribdata, mymap, seqdata, polarity, mainorf, morfs, morfe)
			# rumi change 20141028: we changed from just relmag to also including the 5'UTR power
			# rumi change 20141206: added in protection sum minimum bound
			# rumi change 20150413: removed and (mymap['relmag'] > 0.0)
			# rumi change 20150919: changed the minimum bound to 3
			if ok and (mymap['protection_sum'] >= MIN_RPF_CONSIDERED) and (mymap['within_power_of_in_frame'] > 0) and (mymap['morf_5putr_region_power'] > 0.0): # and (numpy.abs(mymap['phase_of_max_power_freq']) <= in_phase_tolerance):
				was_I_labeled = False
				if starti in labeled_starts:
					was_I_labeled = True
				mymap['labeled'] = was_I_labeled
				mymap['pwm_score'] = kozy
				mymap['chr'] = the_chromosome
				this_chromosome_predictions.append(mymap)

	print '...', len(this_chromosome_predictions), 'putative uORFs identified for chromosome', the_chromosome, 'so far'

	# compute scores for any labeled data we have
	for starti,endi,lorf in labeled_positions:
		mymap = {}
		can_hits = [x for x in orf_annotations if lorf in x[0]]
		if len(can_hits) == 0:
			# try stripping off Ingolia silliness
			can_hits = [x for x in orf_annotations if '-'.join(lorf.split('-')[:-1]) in x[0]]
		if len(can_hits) == 0:
			print 'WARNING - missing data for labeled entry', starti, endi, lorf
			continue
		mainorf, morfs, morfe = can_hits[0]
		#mainorf, morfs, morfe = gimme_the_main_orf_5pUTR(starti, endi)
		if characterize_puorf(starti, endi, ribdata, mymap, seqdata, polarity, mainorf, morfs, morfe):
			mymap['labeled'] = True
			mymap['pwm_score'] = numpy.prod(numpy.sum(PWM_FROM_STARTS*seqdata[starti-6:starti+4],axis=1))
			mymap['chr'] = the_chromosome
			labeled_evidence.append(mymap)

#print(labeled_evidence)

print 'saving results...'
if not FORCE_ONLY_LABELED:
	outs = open(results_dir+the_chromosome+'.predictions','w')
	print outs
	json.dump(this_chromosome_predictions, outs)
	# for z in this_chromosome_predictions:
	#     outs.write(the_chromosome+' '+' '.join([f+'='+(feature_kind_map[f] % z[f]) for f in feature_ordering])+"\n")
	outs.close()

outs = open(results_dir+the_chromosome+'.evidence','w')
print outs
json.dump(labeled_evidence,outs)

outs.close()


