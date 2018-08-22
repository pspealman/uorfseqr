# Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman
#PS 08.11.18 - added imp to load local python files common_formats, skippy.
#PS 08.16.18 - edited print statements, comments for clarity

#PS 08.18.18 - suppressing ignorable warnings
# https://github.com/numpy/numpy/pull/432
import warnings
warnings.simplefilter("ignore")

import imp
imp.load_source('common_formats','../analysis/common_formats.py')
imp.load_source('skippy','../analysis/skippy.py')

import numpy
import os
import sys
import simplejson as json

SETTINGS = json.load(open('analysis.json','r'))

keep_chromo = set(SETTINGS['chromosomes'])

POLSKI = {'positive':'+', 'negative':'-'}
# this is reverse of what you'd ordinarily think because we want to match "is_reverse" in pysam
boolPOLSKI = {'+':False, '-': True}

import common_formats

per_chromosome_TLS_data = {}
per_chromosome_features_masked_out = {}
unwanted_genes = set([])

common_formats.import_gff_file(SETTINGS['GFF_FILE'], per_chromosome_TLS_data, per_chromosome_features_masked_out, unwanted_genes)
print 'GFF KEYS', per_chromosome_TLS_data.keys()
# you calculate this by taking the log of the expected maximum number of items
MAX_ORF_SKIPLIST = SETTINGS['MAX_ORF_SKIPLIST']

import skippy
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

orf_quants = numpy.zeros((jj,),dtype='uint')

import pysam

samfile_name = sys.argv[1]

samfile = pysam.AlignmentFile(samfile_name, 'rb')
print 'SAMFILE references', samfile.references
print 'working on', samfile_name
for achrom, not_strand in sorted(position_to_orf_trees.iterkeys()):
    if achrom not in samfile.references:
        continue
    for anorf, (aleft, aright) in position_to_orf_trees[achrom, not_strand].lookup.iteritems():
        for bb in samfile.fetch(reference=achrom, start=aleft, end=aright):
            if bb.is_reverse==not_strand:
                if not_strand and aleft <= bb.reference_end <= aright:
                    orf_quants[orf_lookup[anorf]] += 1
                elif not not_strand and aleft <= bb.reference_start <= aright:
                    orf_quants[orf_lookup[anorf]] += 1

all_reads = samfile.mapped

fout = open(SETTINGS['processed_dir']+'/'+samfile_name.split('/')[-1]+'-quantified.csv','w')
fout.write('ALL_READS '+str(all_reads)+"\n")
for anorf,ii in orf_lookup.iteritems():
    fout.write(anorf+' '+str(orf_quants[ii])+"\n")

fout.close()
