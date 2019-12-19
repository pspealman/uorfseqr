#Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman
#PS 08.18.18 - suppressing ignorable warnings
# https://github.com/numpy/numpy/pull/432
import warnings
warnings.simplefilter("ignore")

import string
import os
import numpy
import scipy
import scipy.stats
import sys
import simplejson as json

uorf_identifier = ['chr', 'morf_name', 'ss_start', 'ss_end', 'polarity']

SETTINGS = json.load(open('analysis.json','r'))

crappy_starts = set(SETTINGS['ignore_start_codons'])

def load_a_sample(kindof, aname):
	predictions = []
	evidence = []
	for achr in SETTINGS['chromosomes']:
		try:
			predictions.extend([x for x in json.load(open('predictions.'+kindof+aname+'/'+achr+'.predictions')) if x['uorf_sequence'][:3] not in crappy_starts])
		except json.decoder.JSONDecodeError:
			pass
		try:
			evidence.extend([x for x in json.load(open('predictions.'+kindof+aname+'/'+achr+'.evidence')) if x['uorf_sequence'][:3] not in crappy_starts])
		except json.decoder.JSONDecodeError:
			pass

	# load up abundance calculations
	orf_transcript_estimate = {}
	sname = [x['mRNA'] for x in SETTINGS['SAMPLES'] if x['name']==aname][0].split('/')[-1]
	for line in open(SETTINGS['processed_dir']+sname+'-quantified.csv'):
		Line = line.rstrip().split()
		orf_transcript_estimate[Line[0]] = float(Line[1])
	#
	orf_loading_estimates_lookup = {}
	sname = [x['RPF'] for x in SETTINGS['SAMPLES'] if x['name']==aname][0].split('/')[-1]
	for line in open(SETTINGS['processed_dir']+sname+'-quantified.csv'):
		Line = line.rstrip().split()
		orf_loading_estimates_lookup[Line[0]] = float(Line[1])
	return predictions, evidence, orf_transcript_estimate, orf_loading_estimates_lookup



feat_order = []
designfeats = ['ss_start', 'ss_end', 'max_power_freq', 'entropy_of_power', 'median_phase', 'median_power', 'morf_5putr_peak_phase', 'morf_5putr_peak_power', 'morf_5putr_region_power', 'morf_5putr_sum_downstream', 'morf_5putr_sum_upstream', 'morf_dist_TIS', 'morf_dist_TSS', 'within_phase_of_in_frame', 'within_phase_of_max_power_freq', 'within_power_of_in_frame', 'within_power_of_max_power_freq', 'protection_sum', 'pwm_score', 'relative_start_magnitude', 'weighted_avg_phase', 'morf_nterm_region_ends_before', 'morf_nterm_region_start_after', 'morf_nonoverlap','morf_len_TLS',  'within_power_of_max_power_freq', 'within_power_of_in_frame']
dfidx = {x[1]:x[0] for x in enumerate(designfeats)}
keepfeats = ['within_power_of_max_power_freq', 'within_power_of_in_frame',  'morf_dist_TIS', 'morf_dist_TSS', 'morf_nterm_region_ends_before', 'morf_nterm_region_start_after', 'morf_nonoverlap']

def compute_derived_features(rep_data, orf_transcript_ests, orf_loading_ests):
	n = len(rep_data)
	X = numpy.zeros((n,len(designfeats)))
	QQ = numpy.zeros((n,))
	total_rpf_depth = float(orf_loading_ests['ALL_READS'])
	ii = 0
	which_key = []
	for aval in rep_data:
		if aval['morf_name'] in orf_transcript_ests:
			X[ii,:] = [float(aval[z]) for z in designfeats]
			QQ[ii] = orf_transcript_ests[aval['morf_name']]
			if (QQ[ii] + X[ii,dfidx['morf_5putr_region_power']]) >= 20 and QQ[ii]>0:
				ii += 1
				which_key.append([aval[x] for x in uorf_identifier])

	X = X[:ii,:]
	QQ = QQ[:ii]
	n = ii

	# normalize mRNA relative abundances
	QQ /= float(orf_transcript_ests['ALL_READS'])

	synfeats = {}
	for k in keepfeats:
		synfeats[k] = X[:,dfidx[k]] 

	mylengths = numpy.abs(X[:,dfidx['ss_end']]-X[:,dfidx['ss_start']])
	synfeats['spacing_of_max_power'] = X[:,dfidx['max_power_freq']]/mylengths
	synfeats['relative_power_of_max_power_freq'] = X[:,dfidx['within_power_of_max_power_freq']]/QQ
	synfeats['relative_power_of_in_frame'] = X[:,dfidx['within_power_of_in_frame']]/QQ
	synfeats['relative_median_power'] = X[:,dfidx['median_power']]/QQ
	#
	synfeats['region_relative_morf_peak_power']= X[:,dfidx['morf_5putr_peak_power']]/QQ	
	#
	synfeats['abs_phase_of_max_power_freq'] = X[:,dfidx['within_phase_of_max_power_freq']]
	synfeats['abs_phase_of_in_frame'] = X[:,dfidx['within_phase_of_in_frame']]
	#
	synfeats['length_of_puorf'] = mylengths
	#
	#
	avg_upstream_reads = X[:,dfidx['morf_5putr_sum_upstream']]/X[:,dfidx['morf_dist_TSS']]/QQ
	avg_upstream_reads[numpy.isnan(avg_upstream_reads)] = 0
	synfeats['relative_avg_upstream_reads'] = avg_upstream_reads
	#
	avg_downstream_reads = X[:,dfidx['morf_5putr_sum_downstream']]/X[:,dfidx['morf_dist_TIS']]/QQ
	avg_downstream_reads[numpy.isnan(avg_downstream_reads)] = 0
	synfeats['relative_avg_downstream_reads'] = avg_downstream_reads

	synfeats['normalized_start_magnitude'] = X[:,dfidx['relative_start_magnitude']]/total_rpf_depth

	# print feat_order
	synfeatorder = sorted(synfeats.keys())
	if len(feat_order)==0:
		feat_order.extend(synfeatorder)
	else:
		synfeatorder = feat_order
	newX = numpy.vstack([synfeats[x] for x in synfeatorder]).T
	# import pdb
	# pdb.set_trace()
	return which_key, newX, X


for kindof in ['']+['p'+str(i)+'-' for i in range(len(SETTINGS['permutation_seeds']))]:
	outdir = 'calculated'+kindof
	if not os.path.exists(outdir):
		os.mkdir(outdir)

	for asample in SETTINGS['SAMPLES']:
		sname = asample['name']
		print asample, sname
		predictions, evidence, orf_transcript_estimate, orf_loading_estimates_lookup = load_a_sample(kindof, sname)
		which_key, newX, X = compute_derived_features(predictions, orf_transcript_estimate, orf_loading_estimates_lookup)
		fout = open(outdir+'/'+sname+'.table','w')
		for x in which_key:
			fout.write(' '.join([str(z) for z in x])+"\n")
		fout.close()
		numpy.save(outdir+'/'+sname+'-feats.npy', newX)
		numpy.save(outdir+'/'+sname+'-design-feats.npy', X)
		sname = asample['name']+'-evidence'
		which_key, newX, X = compute_derived_features(evidence, orf_transcript_estimate, orf_loading_estimates_lookup)
		fout = open(outdir+'/'+sname+'.table','w')
		for x in which_key:
			fout.write(' '.join([str(z) for z in x])+"\n")
		fout.close()
		numpy.save(outdir+'/'+sname+'-feats.npy', newX)
		numpy.save(outdir+'/'+sname+'-design-feats.npy', X)


	fout = open('here.forder','w')
	for x in feat_order:
		fout.write(x+"\n")
	fout.close()

	fout = open('here.dorder','w')
	for x in designfeats:
		fout.write(x+"\n")
	fout.close()

