# Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman

#PS 08.18.18 - edited to handle line commands for datasets, input settings:

#PS 08.18.18 - suppressing ignorable warnings
# https://github.com/numpy/numpy/pull/432
import warnings
warnings.simplefilter("ignore")

import matplotlib
# Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')

FDR_RECALL = 0.05
PERMUTATION_NULL = False

import string
import os
import numpy
import scipy
import scipy.stats
import sys
import pylab
# from intervaltree import Interval, IntervalTree
from collections import namedtuple
import itertools
import simplejson as json

# PS 08.18.18 - added argparse
import argparse

#PS 08.11.18 - added imp to load local python files common_formats.
import imp
imp.load_source('common_collapse','analysis/common_collapse.py')
imp.load_source('skippy','analysis/skippy.py')

from common_collapse import *
from skippy import *

import cPickle

# PS 08.18.18 - Added parser for argparse
parser = argparse.ArgumentParser()

parser.add_argument('-known',"--known_uorfs")
parser.add_argument('-i',"--input_dir")
parser.add_argument('-o',"--output_dir")

args = parser.parse_args()

if args.known_uorfs:
    known_uorfs = args.known_uorfs
else:
    known_uorfs = 'data/labelled_uorfs/STANDARD-golden.bed'

# PS 08.18.18 - added analysis_json_file_name variable
analysis_json_file_name = args.input_dir + '/analysis.json'

# PS 08.18.18 - added dataset_str variable
if '/' in args.input_dir:
    set_name_str = args.input_dir.rsplit('/',1)[1]
    dataset_str = set_name_str + '/analysis.json'
else:
    dataset_str = args.input_dir + '/analysis.json'

# PS 08.19.18 - Added output_dir variable
if args.output_dir:
    output_dir = args.output_dir
else:
    output_dir = 'results'
    
if not os.path.exists(output_dir):
	os.mkdir(output_dir)

####

gold_standard = []

# PS 08.18.18 - added known_uorfs varaible
for line in open(known_uorfs):
	Line = line.rstrip().split()
	gold_standard.append((Line[0], Line[3].split('.')[0], int(Line[1]), int(Line[2]), Line[5]))

golden = { x: ii for (ii,x) in enumerate(gold_standard) }

# PS 08.18.18 - added analysis_json_file_name variable
FUNC_SETTINGS = json.load(open(analysis_json_file_name,'r'))

feat_order = []
for line in open(FUNC_SETTINGS['PREFIX']+'/here.forder'):
	feat_order.append(line.rstrip())

dfeat_order = []
for line in open(FUNC_SETTINGS['PREFIX']+'/here.dorder'):
	dfeat_order.append(line.rstrip())

dfeat__protection_sum = dfeat_order.index('protection_sum')

(feat1, dir1, feat2, dir2) = (feat_order.index('abs_phase_of_max_power_freq'), False, feat_order.index('normalized_start_magnitude'), True)

# PS 08.18.18 - added dataset_str variable
datasets = [dataset_str]

import itertools
all_results = {}
all_null_results = {}
datasets_order = []
for adataset in datasets:
	x,y = load_data_from_settings(json.load(open(adataset)), feat1, dir1, feat2, dir2, dfeat__protection_sum, force_replay=False)
	for a,b in x:
		all_results[a] = b
		datasets_order.append(a)
	if not PERMUTATION_NULL:
		for a,b in y:
			all_null_results.setdefault(a,[]).append(b)

all_settings = { x : json.load(open(x)) for x in datasets }

standardized_per_sample = {}
for x,y in all_results.iteritems():
	mean_averaged_feats = numpy.mean(y.feats,axis=0)
	std_averaged_feats = numpy.std(y.feats,axis=0)
	standardized_per_sample[x] = (mean_averaged_feats, std_averaged_feats)

def extended_features(b):
	# znew = [z for z in itertools.combinations(range(b.feats.shape[1]),2)]
	# pairfeats = numpy.zeros((b.feats.shape[0], len(znew)))
	# for ii,(i,j) in enumerate(znew):
	# 	pairfeats[:,ii] = numpy.sqrt(numpy.abs(b.feats[:,i]*b.feats[:,j]))
	# return RESULTS(b.tbl, numpy.hstack([b.feats, pairfeats]), b.dfeats, b.countdata, b.freqdata, b.NSAMPLES)
	return b

all_standardized_results = { x : extended_features(apply_zscore(y, standardized_per_sample[x][0], standardized_per_sample[x][1])) for (x,y) in all_results.iteritems() }

def create_permutation_null(y):
	nfeats = numpy.zeros(y.feats.shape)
	idxn = numpy.arange(y.feats.shape[0])
	for j in xrange(y.feats.shape[1]):
		numpy.random.shuffle(idxn)
		nfeats[:,j] = y.feats[idxn,j]
	return RESULTS(y.tbl, nfeats, y.dfeats, y.countdata, y.freqdata, y.NSAMPLES)


all_permutation_null_results = {}
for x,y in all_results.iteritems():
	yy = create_permutation_null(apply_zscore(y, standardized_per_sample[x][0], standardized_per_sample[x][1]))
	all_permutation_null_results[x] = extended_features(yy)


if not PERMUTATION_NULL:
	all_null_results = { x : combine_results(ys) for (x,ys) in all_null_results.iteritems() }
	all_null_results = { x : extended_features(apply_zscore(y, standardized_per_sample[x][0], standardized_per_sample[x][1])) for (x,y) in all_null_results.iteritems() }
else:
	for d in all_standardized_results.keys():
		nfeats = numpy.zeros(all_standardized_results[d].feats.shape)
		for j in xrange(nfeats.shape[1]):
			nfeats[:,j] = all_standardized_results[d].feats[numpy.random.random_integers(0,nfeats.shape[0]-1, size=nfeats.shape[0]),j]
			all_null_results[d] = RESULTS(all_standardized_results[d].tbl, nfeats, all_standardized_results[d].dfeats, all_standardized_results[d].countdata, all_standardized_results[d].freqdata, all_standardized_results[d].NSAMPLES)


all_iso_results = {}
all_noniso_results = {}
for x,y in all_standardized_results.iteritems():
	print 'isolating in',x
	a,b = partition_for_isolates(y)
	all_iso_results[x] = a
	all_noniso_results[x] = b


import rpy2
import rpy2.robjects
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
ro.r("library(glmpath)")

import statsmodels.api as sm
# this only makes sense when you have more than one sample
per_sample_regressions = {}
REGRESSION = namedtuple('REGRESSION', 'fit2 mytau')
globals()[REGRESSION.__name__] = REGRESSION

def train_regression(iso_results):
	# themodel = sm.GLM(iso_results.freqdata, iso_results.feats, family=sm.families.Binomial())
	# themodel_results = themodel.fit()
	# the method used by glmnet and glmpath is to use the lambdas at a fixed percentage of the maximum lambda
	# which is a super clever way to do this
	fractions = numpy.linspace(0,1,100)
	gene_buckets = { x : ii for (ii,x) in enumerate(set([x[1] for x in iso_results.tbl])) }
	bucket_of = numpy.array([gene_buckets[x[1]] for x in iso_results.tbl])
	from sklearn.cross_validation import LabelKFold
	NFOLDS = 10
	mycv = LabelKFold(bucket_of, n_folds=NFOLDS)
	errs = numpy.zeros((NFOLDS, len(fractions)))
	for vv,(trainidxn, testidxn) in enumerate(mycv):
		fit = ro.r.glmpath(iso_results.feats[trainidxn,:],iso_results.freqdata[trainidxn],family="binomial")
		predictions = numpy.asarray(ro.r("predict.glmpath")(fit, iso_results.feats[testidxn], s=ro.r.seq(0,1,length=100), type="response", mode='norm.fraction', standardize=False))
		errs[vv,:] = numpy.sum((predictions-iso_results.freqdata[testidxn,numpy.newaxis])**2 ,axis=0)
	fit2 = ro.r.glmpath(iso_results.feats,iso_results.freqdata,family="binomial")
	mytau = numpy.argmin(numpy.mean(errs,axis=0))
	# return themodel_results
	return REGRESSION(fit2, mytau)


if not os.path.exists('regression'):
	os.mkdir('regression')

for sname, asampl in all_results.iteritems():
	if asampl.NSAMPLES > 1:
		if not os.path.exists('regression/'+sname+'.glmregression'):
			print 'working on', sname
			fit2 = train_regression(all_iso_results[sname])
			per_sample_regressions[sname] = fit2
			cPickle.dump(per_sample_regressions[sname], open('regression/'+sname+'.glmregression','w'))
		else:
			print '...recovering regression from', sname
			per_sample_regressions[sname] = cPickle.load(open('regression/'+sname+'.glmregression','r'))

fout = open('results/feature_weights_per_experiment.csv','w')
fout.write(', '.join(['experiment name', 'Intercept']+feat_order)+"\n")
for x,y in per_sample_regressions.iteritems():
	atau = y.mytau
	areg = y.fit2
	vv = numpy.array(ro.r("predict.glmpath")(areg, s=ro.r.seq(0,1,length=100), type="coefficients", mode='norm.fraction', standardize=False))[atau,:]
	fout.write(x+', '+', '.join([str(a) for a in vv])+"\n")

fout.close()


multiple_sample_datasets = sorted(per_sample_regressions.keys())

# BUILD THE SPREADS

datasets_ordered_by_size = sorted(datasets_order, key=lambda x: all_results[x].NSAMPLES, reverse=True)

nontrivial_datasets_ordered_by_size = [x for x in datasets_ordered_by_size if all_results[x].NSAMPLES>1]


from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.cross_validation import StratifiedKFold

SPREAD = namedtuple('SPREAD', 'nonnull score_real score_null score_perm probs_real probs_null probs_perm score_cutoff probs_cutoff NICE qvals')
globals()[SPREAD.__name__] = SPREAD

def bisect_right(a, x, lo=0, hi=None):
	if lo < 0:
		raise ValueError('lo must be non-negative')
	if hi is None:
		hi = len(a)
	while lo < hi:
		mid = (lo+hi)//2
		if numpy.any(x < a[mid]):
			lo = mid+1
		else: 
			hi = mid
		# print x, a[mid], lo, hi, mid
	return lo

def assess_result(regress, obs_RESULT, null_RESULT, perm_RESULT):
	y = numpy.hstack([numpy.ones((obs_RESULT.feats.shape[0],)), numpy.zeros((null_RESULT.feats.shape[0],))])
	gene_buckets = { x : ii for (ii,x) in enumerate(set([]).union(*[set([x[1] for x in z]) for z in [obs_RESULT.tbl, null_RESULT.tbl]])) }
	bucket_of = numpy.array([gene_buckets[x[1]] for x in obs_RESULT.tbl]+[gene_buckets[x[1]] for x in null_RESULT.tbl])
	working = numpy.vstack([obs_RESULT.feats, null_RESULT.feats])

	from sklearn.ensemble import RandomForestClassifier
	clf = RandomForestClassifier(n_estimators=101, class_weight={0 : 0.5, 1:0.5}, oob_score=False, n_jobs=-1)
	clf.fit(numpy.vstack([obs_RESULT.feats, null_RESULT.feats]), y)

	probs_real = clf.predict_proba(obs_RESULT.feats)[:,1]
	probs_null = clf.predict_proba(null_RESULT.feats)[:,1]
	probs_perm = clf.predict_proba(perm_RESULT.feats)[:,1]

	score_real = numpy.asarray(ro.r("predict.glmpath")(regress.fit2, obs_RESULT.feats, s=ro.r.seq(0,1,length=100), type="response", mode='norm.fraction', standardize=False))[:,regress.mytau]
	score_null = numpy.asarray(ro.r("predict.glmpath")(regress.fit2, null_RESULT.feats, s=ro.r.seq(0,1,length=100), type="response", mode='norm.fraction', standardize=False))[:,regress.mytau]
	score_perm = numpy.asarray(ro.r("predict.glmpath")(regress.fit2, perm_RESULT.feats, s=ro.r.seq(0,1,length=100), type="response", mode='norm.fraction', standardize=False))[:,regress.mytau]

	stacked = numpy.vstack([score_perm, probs_perm]).T
	thresholds = numpy.sort(stacked,axis=0)[::-1]
	hh = numpy.zeros((thresholds.shape[0]+1,))
	for x in stacked:
		hh[bisect_right(thresholds, x)] +=1
	hh = numpy.cumsum(hh)/stacked.shape[0]
	real_bundle = numpy.vstack([score_real, probs_real]).T
	pvals = numpy.ones((score_real.shape[0],))
	for ii in xrange(real_bundle.shape[0]):
		pvals[ii] = hh[bisect_right(thresholds, real_bundle[ii])]
	spvals = numpy.argsort(pvals)
	#G'Sell et al. ForwardStop
	unsorted_qvals = -numpy.cumsum(numpy.log(1-pvals[spvals]))/numpy.arange(1, pvals.shape[0]+1)
	final_qval = numpy.max(numpy.where(unsorted_qvals<=FDR_RECALL)[0])
	score_cutoff, probs_cutoff = thresholds[final_qval]
	qvals = numpy.zeros((pvals.shape[0],))
	qvals[spvals] = numpy.minimum(unsorted_qvals, 1.0)
	# resolve across replicates by choosing the very best score in overlapping predictions
	preNICE = numpy.where(qvals<=FDR_RECALL)[0]
	per_gene = {}
	for z,zz in [(obs_RESULT.tbl[x],x) for x in preNICE]:
		per_gene.setdefault(z[1],[]).append(zz)
	NICE = numpy.zeros((score_real.shape[0],),dtype='bool')
	for g,xs in per_gene.iteritems():
		if len(xs)>1:
			xs = sorted(xs, key=lambda x:score_real[x], reverse=True)
			t = skippy.IntervalSkiplist(3)
			for x in xs:
				y = obs_RESULT.tbl[x]
				if not t.is_interval_contained(y[2], y[3]):
					t.insert_interval(x, y[2], y[3])
			NICE[t.lookup.keys()] = True
		else:
			NICE[xs] = True
	return SPREAD(clf, score_real, score_null, score_perm, probs_real, probs_null, probs_perm, score_cutoff, probs_cutoff, NICE, qvals)


all_SPREADS = {}
for ii,d in enumerate(nontrivial_datasets_ordered_by_size):
	if not os.path.exists('regression/'+d+'.glmspread'):
		print 'working on spread for', d
		all_SPREADS[d] = assess_result(per_sample_regressions[d], all_standardized_results[d], all_null_results[d], all_permutation_null_results[d])
		cPickle.dump(all_SPREADS[d], open('regression/'+d+'.glmspread','w'))
	else:
		print '...recovering spread for', d
		all_SPREADS[d] = cPickle.load(open('regression/'+d+'.glmspread','r'))

##### Produce invidiual dataset scatter plots

obscolors = ['#e34a33', '#e34a33', '#e34a33']
nullcolors = ['#8856a7', '#8856a7', '#8856a7']


def plotagainst(dname, ax1, ax2, obs_RESULT, null_RESULT, the_SPREAD, i):
	if numpy.max(null_RESULT.countdata)<i:
		return
	ax1.set_title(dname)
	ax1.plot([the_SPREAD.probs_cutoff, the_SPREAD.probs_cutoff], [0,1], '-', linewidth=1, color='0.80')
	ax1.plot([0,1], [the_SPREAD.score_cutoff, the_SPREAD.score_cutoff], '-', linewidth=1, color='0.80')
	ax1.hexbin(the_SPREAD.probs_null[null_RESULT.countdata==i], the_SPREAD.score_null[null_RESULT.countdata==i], mincnt=1, marginals=False, cmap='Blues', bins='log', gridsize=50)
	#
	ax1.set_xticks([0.25, 0.5, 0.75, 1.0])
	ax1.set_yticks([0.25, 0.5, 0.75, 1.0])
	ax1.set_xlim([0,1])
	ax1.set_ylim([0,1])
	ax1.set_xlabel('non-null projection', fontsize=9)
	ax1.set_ylabel('regression score', fontsize=9)
	ax2.plot([the_SPREAD.probs_cutoff, the_SPREAD.probs_cutoff], [0,1], '-', linewidth=1, color='0.80')
	ax2.plot([0,1], [the_SPREAD.score_cutoff, the_SPREAD.score_cutoff], '-', linewidth=1, color='0.80')
	ax2.hexbin(the_SPREAD.probs_real[obs_RESULT.countdata==i], the_SPREAD.score_real[obs_RESULT.countdata==i], mincnt=1, marginals=False, cmap='Reds', bins='log', gridsize=50)
	Q = [(obs_RESULT.countdata[ii], the_SPREAD.probs_real[ii], the_SPREAD.score_real[ii]) for (ii,x) in enumerate(obs_RESULT.tbl) if x in gold_standard]
	if len(Q)>0:
		Q = numpy.vstack(Q)
		ax2.plot(Q[:,1], Q[:,2], 'k.')
	ax2.set_xticks([0.25, 0.5, 0.75, 1.0])
	ax2.set_yticks([0.25, 0.5, 0.75, 1.0])
	ax2.set_xlim([0,1])
	ax2.set_ylim([0,1])
	ax2.set_xlabel('non-null projection', fontsize=9)
	ax2.set_ylabel('regression score', fontsize=9)

for d in nontrivial_datasets_ordered_by_size:
	fout_all = open('results/'+d+'-all-qvals-with-pwm-hack.bed','w')
	for ii in range(len(all_results[d].tbl)):
		x = all_results[d].tbl[ii]
		if x[4]=='+':
			jj = str(x[2])
		else:
			jj = str(x[3])
		pre = "\t".join([x[0], str(x[2]), str(x[3]), x[1]+'.'+jj])
		bland = pre+'\t'+str(int(1000*all_SPREADS[d].qvals[ii]))+'\t'+str(x[4])+'\t'
		bland += str(all_results[d].dfeats[ii, dfeat_order.index('pwm_score')])+"\n"
		fout_all.write(bland)
	fout_all.close()

for d in nontrivial_datasets_ordered_by_size:
	fout_all = open('results/'+d+'-all-qvals.bed','w')
	fout_significants = open('results/'+d+'-significant-qvals.bed','w')
	fout_resolved = open('results/'+d+'-significant-resolved-qvals.bed','w')
	for ii in range(len(all_results[d].tbl)):
		x = all_results[d].tbl[ii]
		if x[4]=='+':
			jj = str(x[2])
		else:
			jj = str(x[3])
		pre = "\t".join([x[0], str(x[2]), str(x[3]), x[1]+'.'+jj])
		bland = pre+'\t'+str(int(1000*all_SPREADS[d].qvals[ii]))+'\t'+str(x[4])+"\n"
		fout_all.write(bland)
		if all_SPREADS[d].qvals[ii]<=FDR_RECALL:
			fout_significants.write(bland)
			if all_SPREADS[d].NICE[ii]:
				fout_resolved.write(bland)
	fout_all.close()
	fout_significants.close()
	fout_resolved.close()

for d in nontrivial_datasets_ordered_by_size:
	fout_all = open('results/'+d+'-all-scores.bed','w')
	fout_significants = open('results/'+d+'-significant-scores.bed','w')
	fout_resolved = open('results/'+d+'-significant-resolved-scores.bed','w')
	for ii in range(len(all_results[d].tbl)):
		x = all_results[d].tbl[ii]
		if x[4]=='+':
			jj = str(x[2])
		else:
			jj = str(x[3])
		pre = "\t".join([x[0], str(x[2]), str(x[3]), x[1]+'.'+jj])
		bland = pre+'\t'+str(int(1000*all_SPREADS[d].score_real[ii]))+'\t'+str(x[4])+"\n"
		fout_all.write(bland)
		if all_SPREADS[d].qvals[ii]<=FDR_RECALL:
			fout_significants.write(bland)
			if all_SPREADS[d].NICE[ii]:
				fout_resolved.write(bland)
	fout_all.close()
	fout_significants.close()
	fout_resolved.close()

for d in nontrivial_datasets_ordered_by_size:
	fout_all = open('results/'+d+'-all-specificity.bed','w')
	fout_significants = open('results/'+d+'-significant-specificity.bed','w')
	fout_resolved = open('results/'+d+'-significant-resolved-specificity.bed','w')
	for ii in range(len(all_results[d].tbl)):
		x = all_results[d].tbl[ii]
		if x[4]=='+':
			jj = str(x[2])
		else:
			jj = str(x[3])
		pre = "\t".join([x[0], str(x[2]), str(x[3]), x[1]+'.'+jj])
		bland = pre+'\t'+str(int(1000*all_SPREADS[d].probs_real[ii]))+'\t'+str(x[4])+"\n"
		fout_all.write(bland)
		if all_SPREADS[d].qvals[ii]<=FDR_RECALL:
			fout_significants.write(bland)
			if all_SPREADS[d].NICE[ii]:
				fout_resolved.write(bland)
	fout_all.close()
	fout_significants.close()
	fout_resolved.close()

all_NICE = {}
for d in nontrivial_datasets_ordered_by_size:
	nice_idxn = numpy.where(all_SPREADS[d].NICE)[0]
	newtbl = [all_results[d].tbl[x] for x in nice_idxn]
	newfeats = all_results[d].feats[nice_idxn]
	newdfeats = all_results[d].dfeats[nice_idxn]
	newcountdata = all_results[d].countdata[nice_idxn]
	newfreqdata = all_results[d].freqdata[nice_idxn]
	sigresult = RESULTS(newtbl, newfeats, newdfeats, newcountdata, newfreqdata, all_results[d].NSAMPLES)
	cPickle.dump(sigresult, open('results/NICE_OF-'+d+'.RESULT','w'))
	all_NICE[d] = sigresult



