# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 16:39:54 2018

nteseqr
Description:
Takes ribosome profiling data and annotated reference genomes to identify 
N-terminal extensions (NTE)

@author: Pieter Spealman

12.3.2019 (Concentrate Feed) Implementing improvements for MiMB publication
    _x_ nteseqr, identify high likelihood NTEs and filter them from uORFseqr
       
"""
import datetime

#
#import multiprocessing
import re
import numpy as np
import pickle
import scipy.stats as stats
import argparse
import subprocess
import os

def output_handler(output):
    if len(output.strip()) > 1:
        print(output)
        
#function handles help commands
def help_dialog():
    monolog=('Manual for nteseqr\n'+
    '#=============================================================#\n'+
    'Lead programmer: Pieter Spealman pspsealman@nyu.edu\n'+
    'Release version: 1.0 \n'+
    'Release date: 12.31.19 \n'+
    'Description:\n\t nteseqr identifies N-terminal extensions using ribosome \n'+
    '\t profiling and RNAseq data.\n'+
    '\t\t Briefly, NTE-seqr attempts to identify N-terminal extensions by first \n'+
    '\t finding all regions upstream of main ORF start codons and the nearest \n'+
    '\t in-frame upstream stop codon. These search regions are then scanned to \n'+
    '\t identify genes with large numbers of in-frame ribosomes. Search regions \n'+
    '\t are also scanned for AUG and NCC start codons. We presume that the start \n'+
    '\t codon most likely to function as the initiation site will have a confluence \n'+
    '\t of features: higher relative start magnitude, higher relative translational \n'+
    '\t efficiency, and a significant fraction of total in-frame ribosomes.'+
    'Citation:'+
    'Copyright MIT License - Pieter Spealman'
    '#=============================================================#\n'+
    'For demonstration use:\n\t python nteseqr.py -demo\n'+
    'To run a install test using defaults, use:\n\t python uorfseqr.py -test\n'+
    '')
    print(monolog)
   
# 
def demo():
    monolog = ('\tStep 1. -load command loads and assigns reads. This will need to be done for '+
               'every pair of RPF and RNA files. Here we use only two, the minimum number.\n')
    print(monolog)
    monolog = ('\t\tUsage:\n\tpython nteseqr.py -load -gff <path_to_gff_file> -fa <path to reference fasta file>\n'+
               '\t-sample <sample_name> <path_to_RPF_bam_file> <path_to_RNA_bam_file> -o <output_prefix>'+
    '\t\tExample:\n\tpython nteseqr.py -load -gff analysis/saccharomyces_cerevisiae.gff -fa data/reference_genomes/Scer_SacCer3.fa -samples Scer_A data/bam/Scer_A_RPF_10.bam data/bam/Scer_A_mRNA_10.bam -o Scer_A_nte\n')
    print(monolog)
    
    monolog = ('\tStep 2. -load command loads and assigns reads for replicate 2.\n')
    print(monolog)
    monolog = ('\t\tUsage:\n\tpython nteseqr.py -load -gff <path_to_gff_file> -fa <path to reference fasta file>\n'+
               '\t-sample <sample_name> <path_to_RPF_bam_file> <path_to_RNA_bam_file> -o <output_prefix>'+
    '\t\tExample:\n\tpython nteseqr.py -load -gff analysis/saccharomyces_cerevisiae.gff -fa data/reference_genomes/Scer_SacCer3.fa -samples Scer_B data/bam/Scer_B_RPF_10.bam data/bam/Scer_B_mRNA_10.bam -o Scer_B_nte\n')
    print(monolog)

    monolog = ('\tStep 3. -eval command generates candidate uORFs using the previouslt loaded samples.\n'+
               'Note that the -samples values here are the output (-o) from the previous two steps.\n')
    print(monolog)
    monolog = ('\t\tUsage:\n\tpython nteseqr.py -eval -samples <name_of_sample_1> <name_of_sample_2> -o <output_prefix>\n'+
    '\t\tExample:\n\tpython nteseqr.py -eval -samples Scer_A_nte Scer_B_nte -o scer.demo/combined\n')
    print(monolog)
    
    monolog = ('\tStep 4. NTE candidate file.\n'+
               '\tThe highest scoring alternative translation initiation site (aTIS) for each NTE event\n'+
               'is output in bed file format as )
    print(monolog)
# 
def test():
    monolog = ('=== Currently Testing nteseqr.py ===')
    print(monolog)

    monolog = ('\tTesting Step 1a. -load command loads and assigns reads for replicate 1.\n')
    print(monolog)    
    bashCommand = ('python nteseqr.py -load -gff analysis/saccharomyces_cerevisiae.gff -fa data/reference_genomes/Scer_SacCer3.fa -samples Scer_A ../data/bam/Scer_A_RPF_10.bam ../data/bam/Scer_A_mRNA_10.bam -o Scer_A_nte')
    print(bashCommand)
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    monolog = ('\tTesting Step 1b. -load command loads and assigns reads for replicate 2.\n')
    print(monolog)    
    bashCommand = ('python nteseqr.py -load -gff analysis/saccharomyces_cerevisiae.gff -fa data/reference_genomes/Scer_SacCer3.fa -samples Scer_B ../data/bam/Scer_B_RPF_10.bam ../data/bam/Scer_B_mRNA_10.bam -o Scer_B_nte')
    print(bashCommand)
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    monolog = ('\tTesting Step 3. -eval command to evaluate candidates based on loaded expression data.\n')
    print(monolog)    
    bashCommand = ('python nteseqr.py -eval -samples Scer_A_nte Scer_B_nte -o scer.demo/combined')
    print(bashCommand)
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))

### Argparser definitions        
parser = argparse.ArgumentParser()

##nteseqr
'''
python nteseqr.py -load -gff saccharomyces_cerevisiae.gff -fa Scer_SacCer3.fa -samples starved_r2 r2_RPF.bam r2_RNA.bam -o set_r2
python nteseqr.py -load -gff saccharomyces_cerevisiae.gff -fa Scer_SacCer3.fa -samples starved_r3 r3_RPF.bam r3_RNA.bam -o set_r3
python nteseqr.py -eval -samples set_r2 set_r3 -o set_combined
---
python nteseqr.py -load -gff saccharomyces_cerevisiae.gff -fa Scer_SacCer3.fa -samples starved_r2 ./Scer_r2Starved_RPF_genome_Aligned.out.bam ./Scer_r2Starved_mRNA_genome_Aligned.out.bam -o starved_r2
python nteseqr.py -load -gff saccharomyces_cerevisiae.gff -fa Scer_SacCer3.fa -samples starved_r3 ./Scer_r3Starved_RPF_genome_Aligned.out.bam ./Scer_r3Starved_mRNA_genome_Aligned.out.bam -o starved_r3
python nteseqr.py -eval -samples starved_r2 starved_r3 -o starved_combined

'''
# help dialog arguments
parser.add_argument('-man',"--manual", action='store_true')
parser.add_argument('-demo',"--demo",action='store_true')
parser.add_argument('-test',"--test",action='store_true')

# Load reads arguments
parser.add_argument('-load',"--load_reads", action='store_true')
parser.add_argument('-i',"--input_file")
parser.add_argument('-o',"--output_file")
parser.add_argument('-fa',"--fa_file")
parser.add_argument('-gff',"--gff_file")
parser.add_argument('-samples', '--sample_list', nargs='+')
parser.add_argument('-gt',"--gene_tag")
parser.add_argument('-tl',"--transcript_leader_tag")
parser.add_argument('-3p',"--three_prime_UTR_tag")
parser.add_argument('-min_tl',"--minimum_length_transcript_leader")
parser.add_argument('-min_3p',"--minimum_length_three_prime_UTR")
parser.add_argument('-mask_tl',"--mask_length_transcript_leader")
parser.add_argument('-mask_3p',"--mask_length_three_prime_UTR")
parser.add_argument('-defualt',"--default_search_region_length")

#evaluate arguments
parser.add_argument('-eval',"--evaluate", action='store_true') 
parser.add_argument('-min_reads',"--minimum_reads")
#
args = parser.parse_args()
###
#Common dictionaries
start_codons = ['ATG', 'TTG', 'GTG', 'CTG', 'ATC', 'ATA', 'ATT', 'ACG']
stop_codons = ['TAG', 'TAA', 'TGA']
strand_to_sign = {0:'+',1:'-'}
complement = {'A':'T','G':'C','T':'A','C':'G'}

###
if args.manual:
    help_dialog()
    
#if args.demo:
#    demo()
    
if args.test:
    test()

def parse_cigar(cigar, sequence):
    """This function calculates the offset for the read based on the match
    """
    # TODO - maybe improve to handle '28M1I4M', 'TCAGGGAAATATTGATTTACCCAAAAAAAGACG'
    if cigar.count('M') == 1:
        left_cut = 0
        right_cut = 0
        
        left_list = re.split('M|S|D|I|H|N', cigar.split('M')[0])[0:-1]
        M = re.split('M|S|D|I|H|N', cigar.split('M')[0])[-1]
        right_list = re.split('M|S|D|I|H|N', cigar.split('M')[1])
        
        for each in left_list:
            if each: 
                left_cut += int(each)
                
        for each in right_list:
            if each: 
                right_cut -= int(each)
        
        n_cigar = ('{}M').format(M)

        if right_cut:
            n_sequence = sequence[left_cut:right_cut]
        else:
            n_sequence = sequence[left_cut:]
                        
        #print (left_cut, right_cut,  n_cigar, n_sequence)
        return(True, n_cigar, n_sequence)
            
    else:
        return(False, '', '')

def unpackbits(x, num_bits=12):
    xshape = list(x.shape)
    x = x.reshape([-1,1])
    to_and = 2**np.arange(num_bits).reshape([1,num_bits])
    upb = (x & to_and).astype(bool).astype(int).reshape(xshape + [num_bits])

    #0  (rp)    read_paired
    #1  (rmp)    read_mapped_in_proper_pair
    #2  (ru)    read_unmapped
    #3  (mu)    mate_unmapped
    #4  (rrs)    read_reverse_strand
    #5  (mrs)    mate_reverse_strand
    #6  (fip)    first_in_pair
    #7  (sip)    second_in_pair
    #8  (npa)    not_primary_alignment
    #9  (rfp)    read_fails_platform
    #10 (pcr)    read_is_PCR_or_optical_duplicate
    #11 (sa)    supplementary_alignment
    
    """ DISCORDANT definition (from samblaster)
        Both side of the read pair are mapped (neither FLAG 0x4 or 0x8 is set).
        The properly paired FLAG (0x2) is not set.
        Note: We implemented an additional criteria to distinguish between strand re-orientations and distance issues
        Strand Discordant reads must be both on the same strand.
    """
        
    """ SPLIT READS
        Identify reads that have between two and --maxSplitCount [2] primary and supplemental alignments.
        Sort these alignments by their strand-normalized position along the read.
        Two alignments are output as splitters if they are adjacent on the read, and meet these criteria:
            each covers at least --minNonOverlap [20] base pairs of the read that the other does not.
            the two alignments map to different reference sequences and/or strands. 
            the two alignments map to the same sequence and strand, and represent a SV that is at least --minIndelSize [50] in length, 
            and have at most --maxUnmappedBases [50] of un-aligned base pairs between them.
        Split read alignments that are part of a duplicate read will be output unless the -e option is used.
    """
    
    return(upb)     

def os_mkdir(in_name):    
    if '/' in in_name:
        directory_name = in_name.rsplit('/',1)[0]
        
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
###

''' Handle inputs and defaults: 

 users can set their own gff file to define 5'UTRs (aka. transcript leaders), 3'UTRs,
 Transcription start sites, poly-a sites, and main orf coordinates.
 Otherwise the standard gff for S.cerevisiae (from Spealman and Naik, Genome Research, 2017) is loaded.

'''

if args.gff_file:
    gff_file = args.gff_file
else:
    gff_file = '../data/reference_genomes/saccharomyces_cerevisiae.gff'

''' Given a gff and fasta identity the upstream in frame stops, 
    then downstream inframe starts, 
    then ask if reads from a sam map to those starts. 
'''

if args.gene_tag:
    gene_tag = args.gene_tag
else:
    gene_tag = 'gene'
    
if args.transcript_leader_tag:
    tl_tag = args.transcript_leader_tag
else:
    tl_tag = 'five_prime_UTR'
    
if args.three_prime_UTR_tag:
    tp_tag = args.three_prime_UTR_tag
else:
    tp_tag = 'three_prime_UTR'

if args.minimum_length_transcript_leader:
    min_tl = int(args.minimum_length_transcript_leader)
else:
    min_tl = 15
    
# TODO - future version turns this on
#if args.minimum_length_three_prime_UTR:
#    min_3p = int(args.minimum_length_three_prime_UTR)
#else:
#    min_3p = 15
min_3p = 0
    
if args.mask_length_transcript_leader:
    mask_tl = int(args.mask_length_transcript_leader)
else:
    mask_tl = 3
    
if args.mask_length_three_prime_UTR:
    mask_3p = int(args.mask_length_three_prime_UTR)
else:
    mask_3p = 3
    
if args.default_search_region_length:
    dsrl = int(args.default_search_region_length)
else:
    dsrl = 0
    
if args.minimum_reads:
    minimum_reads = int(args.minimum_reads)
else:
    minimum_reads = 3
            
def parse_name(field):
    field = field.strip()
    
    if 'PARENT=' in field:
        field = field.split('PARENT=')[1]
    if ';' in field:
        field = field.split(';')[0]
    if '_mRNA' in field:
        field = field.split('_mRNA')[0]
        
    return(field)
    
def parse_gff(gff_name, gene_tag, tl_tag, tp_tag, min_utr, min_3p, mask_tl, mask_3p):
    global fasta_dict
    global filter_nt
    
    coord_dict = {}
    chromosome_set = set()
    gff_file = open(gff_name)
    
    #Build coord_dict
    for line in gff_file:
        if line[0] != '#':
            name = parse_name(line.split('\t')[8])
            chromo = line.split('\t')[0]
            region = line.split('\t')[2]
            start = int(line.split('\t')[3])-1
            stop = int(line.split('\t')[4])
            sign = line.split('\t')[6]   
            
            if chromo not in chromosome_set:
                chromosome_set.add(chromo)
            
            if name not in coord_dict:
                coord_dict[name] = {'chromo': chromo, 'sign': line.split('\t')[6], 'tl':'', 'gene':'', 'tp': '', 'tl_mask':'', 'tp_mask':''}
            
            if chromo not in filter_nt:
                filter_nt[chromo] = set()
            
            if region == gene_tag:
                coord_dict[name]['gene'] = (start, stop)
            
            
            if region == tl_tag:
                if abs(stop - start) > min_utr:
                    coord_dict[name]['tl'] = (start, stop)
                    
                    if sign == '+':
                        for nt in range(stop-mask_tl, stop+1):
                            filter_nt[chromo].add(nt)
                    else:
                        for nt in range(start, start+mask_tl+1):
                            filter_nt[chromo].add(nt)
                else:
                    coord_dict[name]['tl'] = 'too_small'
                    
# TODO: Future version - C-terminal extensions should be a similar work flow
#            if region == tp_tag:
#                if abs(stop-start) > min_3p:
#                    coord_dict[name]['tp'] = (start, stop)
#                    
#                    if sign == '+':
#                        for nt in range(start, start+mask_3p+1):
#                            filter_nt[chromo].add(nt)
#                    else:
#                        for nt in range(stop-mask_3p, stop+1):
#                            filter_nt[chromo].add(nt)
#                else:
#                    coord_dict[name]['tp'] = 'too_small'
                
    gff_file.close()
    
    # cycle through coord_dict, remove anything with 'too_small'
    remove_set = set()
    
    for name, region_dict in coord_dict.items():
        if not region_dict['gene']:
            #TODO improve readout
            remove_set.add(name)
        else:
            for each_region in ['tl','tp']:
                if region_dict[each_region] == 'too_small':
                    remove_set.add(name)
    
    for remove in remove_set:
        pop = coord_dict.pop(remove)
        if '@' not in remove:
            outline = ('Removing {} for to short of a UTR: {}').format(remove, pop)
            print(outline)
    
    return(coord_dict, chromosome_set)
    
def parse_fasta(fasta_name):
    fasta_file = open(fasta_name)
    
    fasta_dict = {}
    
    for line in fasta_file:
        line = line.strip()
        if line[0] == '>':
            name = line.split('>')[1].split(' ')[0]
            name = name.strip()
            
            if name not in fasta_dict:
                fasta_dict[name]=''
            else:
                print('Error in FASTQ chromosome name, duplicate names identified.\n Names are after the carrot ">" and before the space " " - Make sure each name is unique. ')
                quit()
        else:
            fasta_dict[name]+=line
            
    return(fasta_dict)
            
def rev_comp(seq):
    seq = seq.upper()
        
    seq = seq[::-1]
    rev_seq = ''

    for each in seq:
        rev_seq+= complement[each]
        
    return(rev_seq)
        
def use_dsrl(name, region, first_pass_dict, coord_dict):
    chromo = first_pass_dict[name]['chromo']
    sign = first_pass_dict[name]['sign']
    gleast = coord_dict[name]['gene'][0]
    gmost = coord_dict[name]['gene'][1]
    
    if region == 'tl': 
        if sign == '+':
            start = int(gleast)-dsrl
            stop = int(gleast)-1
            first_pass_dict[name][region] = fasta_dict[chromo][start:stop]
            coord_dict[name][region]=(start, stop)
            
        if sign == '-':
            start = int(gmost)+1
            stop = int(gmost)+dsrl
            first_pass_dict[name][region] = fasta_dict[chromo][start:stop]
            coord_dict[name][region]=(start, stop)
            
    if region == 'tp': 
        if sign == '-':
            start = int(gleast)-dsrl
            stop = int(gleast)-1
            first_pass_dict[name][region] = fasta_dict[chromo][start:stop]
            coord_dict[name][region]=(start, stop)
            
        if sign == '+':
            start = int(gmost)+1
            stop = int(gmost)+dsrl
            first_pass_dict[name][region] = fasta_dict[chromo][start:stop]
            coord_dict[name][region]=(start, stop)
        
    return(first_pass_dict, coord_dict)
    
#def derive_coordinates(start, stop, triplet_step, sign, runmode):
#    if (runmode == 'tl' and sign == '+') or (runmode == 'tp' and sign == '-'):
#        sr_start = stop - (triplet_step*3)-3
#        sr_stop = stop
#    if (runmode == 'tl' and sign == '-') or (runmode == 'tp' and sign == '+'):        
#        sr_start = start
#        sr_stop = start + (triplet_step*3)+3
#    
#    return(sr_start, sr_stop)
#    
#def recover_sequence(triplet_step, tl_list, sign, runmode):
#    if (runmode == 'tl' and sign == '+') or (runmode == 'tp' and sign == '-'):
#        tl_list = tl_list[::-1]
#        tl_list = tl_list[1:]
        
def find_stop(name, start, stop, seq, sign, runmode):
    tl_list = []
    sr_seq = []

    #if (runmode == 'tl' and sign == '+') or (runmode == 'tp' and sign == '-'):
    if (runmode == 'tl' and sign == '+'):
        seq = seq[::-1]
        
        for triplet_step in range(len(seq)//3):
            triplet = seq[(3*triplet_step):(3*triplet_step)+3]
            tl_list.append(triplet[::-1])         
        
        triplet_step = 0
        for triplet in tl_list:
            sr_seq.append(triplet)
            if triplet in stop_codons:
                sr_start = stop - (3*triplet_step)-3
                sr_stop = stop
                return(sr_start, sr_stop, sr_seq)
                
            triplet_step += 1
            
        return(start, stop, sr_seq)

    #if (runmode == 'tl' and sign == '-') or (runmode == 'tp' and sign == '+'):      
    if (runmode == 'tl' and sign == '-'):  
        #reverse so the first codon is the one next to the start...
        seq = seq[::-1]
        #step out and reverse each codon to the original (negative) orientation, add to list
        for triplet_step in range((len(seq))//3):
            triplet = seq[(3*triplet_step):(3*triplet_step)+3]
            tl_list.append(triplet[::-1])
        
        # scan for stop, first stop send the modified sequence out for detection
        for triplet in tl_list:
            sr_seq.append(triplet)
            
            if triplet in stop_codons:
                sr_start = start
                sr_stop = start + (3*triplet_step)+3
                return(sr_start, sr_stop, sr_seq)

            triplet_step += 1
      
        return(start, stop, sr_seq)
        
    return(start, stop, sr_seq)

def find_starts(name, chromo, sign, start, stop, tl_list, runmode):
    global atis_id_dict
    global coord_dict
    
    tl_list = tl_list[::-1]
    
    seq = ''
    if sign == '-':
        for each_codon in tl_list:
            each_codon = each_codon[::-1]
            seq += each_codon
    else:
        for each_codon in tl_list:
            seq += each_codon
            
    start_coords = {'full':(start+1, stop)}
    #start_coords['full']=(start+1, stop)
    
    #if (runmode == 'tl' and sign == '+') or (runmode == 'tp' and sign == '-'):
    if (runmode == 'tl' and sign == '+'):
        triplet_step = 0       
        for triplet in tl_list:
            #sr_seq += triplet
            
            if triplet in start_codons:
                sr_seq = ''
                for codon in tl_list[triplet_step:]:
                    sr_seq += codon
                                    
                hash_line = ('{}_{}_{}_{}').format(name, triplet, triplet_step, runmode)
                atis_id = hash(hash_line)
                start_coords = {'atis':{},'sr':{}, 'up':{}, 'gene':[coord_dict[name]['gene'][0],coord_dict[name]['gene'][1]], 'meta':{'name':name, 'chromo':chromo, 'sign': sign, 'region':runmode, 'seq':sr_seq, 'triplet':triplet}, 'full': (start+1, stop)}
                
                #calc aTIS
                #+1 for gff format
                sr_start = start + (3*triplet_step) + 1
                sr_stop = sr_start + 3
                start_coords['atis'] = (sr_start, sr_stop)

                #calc whole region
                sr_start = start + (3*triplet_step) + 1
                sr_stop = stop
                start_coords['sr'] = (sr_start, sr_stop)

                #calc upstream
                sr_start = start + 1
                sr_stop = start + (3*triplet_step)
                start_coords['up'] = (sr_start, sr_stop)

                atis_id_dict[atis_id] = start_coords

            triplet_step += 1
    
    #if (runmode == 'tl' and sign == '-') or (runmode == 'tp' and sign == '+'):
    if (runmode == 'tl' and sign == '-'):
        #tl_list = tl_list[::-1]
        triplet_step = 0
        for triplet in tl_list:
            
            if triplet in start_codons:
                rt_step = len(tl_list)-triplet_step
                sr_seq = ''
                for codon in tl_list[triplet_step:]:
                    sr_seq += codon
                    
                hash_line = ('{}_{}_{}_{}').format(name, triplet, triplet_step, runmode)
                atis_id = hash(hash_line)
                start_coords = {'atis':{},'sr':{}, 'up':{}, 'gene':[coord_dict[name]['gene'][0],coord_dict[name]['gene'][1]], 'meta':{'name':name, 'chromo':chromo, 'sign': sign, 'region':runmode, 'seq':sr_seq, 'triplet':triplet}, 'full': (start+1, stop)}
                
                #calc atis
                sr_start = start + (3*rt_step) - 2
                sr_stop = start + (3*rt_step)
                start_coords['atis'] = (sr_start, sr_stop)
                
                #calc whole region
                sr_start = start + 1 
                sr_stop = start + (3*rt_step)
                start_coords['sr'] = (sr_start, sr_stop)
                
                #calc upstream
                sr_start = start + (3*rt_step) + 1
                sr_stop = stop
                start_coords['up'] = (sr_start, sr_stop)

                atis_id_dict[atis_id] = start_coords                    
           
            triplet_step += 1
        
    return(start_coords)
                
def build_search_region(coord_dict, fasta_dict, dsrl):
    first_pass_dict = {}

    #make each, if possible
    for name, deets in coord_dict.items():
        chromo = deets['chromo']
        if chromo in fasta_dict:
            sign = deets['sign']
            first_pass_dict[name] = {'chromo': chromo, 'sign': sign, 'tl':'', 'gene':'', 'tp':''}
            
            for region in ['tl', 'gene', 'tp']:
                if deets[region]:
                    start = int(deets[region][0])
                    stop = int(deets[region][1])+1
                            
                    if sign == '+':
                        first_pass_dict[name][region] = fasta_dict[chromo][start:stop-1]
                        
                    if sign == '-':
                        first_pass_dict[name][region] = rev_comp(fasta_dict[chromo][start:stop])
                
    #use 'gene' and default search region lenght to fill in those that are absent
    for name, deets in first_pass_dict.items():
        chromo = deets['chromo']
        sign = deets['sign']
        
        for region in ['tl', 'tp']:
            if not deets[region]:
                first_pass_dict, coord_dict = use_dsrl(name, region, first_pass_dict, coord_dict) 
    
    #scan for stops and define regions       
    search_region_dict = {}
    assign_region_dict = {'tl':{}, 'tp':{}}
    flanking_region_dict = {'tl':{}, 'tp':{}}
    
    for name, deets in first_pass_dict.items():
        chromo = deets['chromo']
        sign = deets['sign']
        
        search_region_dict[name] = {'chromo': chromo, 'sign': sign, 'tl':'', 'tp':'', 'starts':{'tl':{},'tp':{}}}
        
        for region in ['tl', 'tp']:
            r_start = coord_dict[name][region][0]
            r_stop = coord_dict[name][region][1]
            
            sr_start, sr_stop, sr_seq = find_stop(name, r_start, r_stop, deets[region], sign, region)

            search_region_dict[name]['starts'][region] = find_starts(name, chromo, sign, sr_start, sr_stop, sr_seq, region) 
            
            if chromo not in assign_region_dict[region]:
                assign_region_dict[region][chromo]={}
                
            if chromo not in flanking_region_dict[region]:
                flanking_region_dict[region][chromo]={}
                
            for nt in range(sr_start, sr_stop+1):
                if nt not in assign_region_dict[region][chromo]:
                    gene_set = set()
                    gene_set.add(name)
                    assign_region_dict[region][chromo][nt] = gene_set                
                else:
                    assign_region_dict[region][chromo][nt].add(name)
                
            for f_nt in range(r_start, r_stop+1):
                if f_nt not in range(sr_start, sr_stop+1):
                    if f_nt not in flanking_region_dict[region][chromo]:
                        gene_set = set()
                        gene_set.add(name)
                        flanking_region_dict[region][chromo][f_nt] = gene_set
                    else:
                        flanking_region_dict[region][chromo][f_nt].add(name)
                        
    return(search_region_dict, assign_region_dict, flanking_region_dict)
    
def convert_to_sam(each_sample):
    if each_sample[-4:] == '.bam':	
        new_name = each_sample.split('.bam')[0]+'.sam'
        
        monolog = ('\tLoading bam file {}.\n').format(str(each_sample))	
        print(monolog)    	
        bashCommand = ('samtools view -h -o {} {}').format(new_name, each_sample)	
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))	
        
        return(new_name)
    
    if each_sample[-4:] == '.sam':	
        return(each_sample)

def assign_reads(uid, chromo, start, stop, sign, runmode):   
    global sample_dict
    global search_region_dict
    global assign_region_dict
    global flanking_region_dict
    global filter_nt
    
    hit_ct = 0

    filter_region = filter_nt[chromo]

    for region in ['tl','tp']:        
        for nt in range(start, stop+1):
            if nt not in filter_region:
                if nt in assign_region_dict[region][chromo]:
                    gene_set = assign_region_dict[region][chromo][nt]
    
                    for gene in gene_set:
                        if sign == search_region_dict[gene]['sign']:
                            
                            hit_ct += 1 
                            if gene not in sample_dict[runmode]['sr']:
                                uid_set = set()
                                uid_set.add(uid)
                                sample_dict[runmode]['sr'][gene] = {}
                                sample_dict[runmode]['sr'][gene][nt] = uid_set
                            
                            else:
                                if nt not in sample_dict[runmode]['sr'][gene]:
                                    uid_set = set()
                                    uid_set.add(uid)
                                    sample_dict[runmode]['sr'][gene][nt] = uid_set
                                else:
                                    sample_dict[runmode]['sr'][gene][nt].add(uid)
    
                else:
                    if nt in flanking_region_dict[region][chromo]:
                        gene_set = flanking_region_dict[region][chromo][nt]
                        for gene in gene_set:
                            if sign == search_region_dict[gene]['sign']:
                                hit_ct += 1 
                                if gene not in sample_dict[runmode]['fl']:
                                    uid_set = set()
                                    uid_set.add(uid)
                                    sample_dict[runmode]['fl'][gene] = {}
                                    sample_dict[runmode]['fl'][gene][nt] = uid_set     
                                    
                                else:
                                    if nt not in sample_dict[runmode]['fl'][gene]:
                                        uid_set = set()
                                        uid_set.add(uid)
                                        sample_dict[runmode]['fl'][gene][nt] = uid_set
                                    else:
                                        sample_dict[runmode]['fl'][gene][nt].add(uid)
                     
    return(hit_ct)  
    
def load_reads(convert_name, chromosome_set, search_region_dict, assign_region_dict, flanking_region_dict, runmode):
    global sample_bam_dict
    #TODO: Future version - Basic version only extracts 28M reads, make this better

    sam_file = open(convert_name)
    
    ct = 0
    hit_ct = 0        
    header_list = []

    for line in sam_file:
        if line[0] == '@':
            header_list.append(line)
            
        else:
            ct += 1
            cigar = line.split('\t')[5]
            chromo = line.split('\t')[2]
            
            if chromo in chromosome_set:
                uid = line.split('\t')[0]+'~'+str(ct)
                flag = line.split('\t')[1]
        
                start = int(line.split('\t')[3])
                
                sign = strand_to_sign[unpackbits(np.array([int(flag)]))[0][4]]
                sequence = line.split('\t')[9]
                
                process, n_cigar, n_sequence = parse_cigar(cigar, sequence)
                
                stop = start + len(n_sequence)

                new_hits = assign_reads(uid, chromo, start, stop, sign, runmode)
                hit_ct += new_hits
                if new_hits > 0:    
                    sample_bam_dict[runmode][uid] = line
                
                if '28M' in cigar and runmode == 'RPF':                    
                    if chromo in chromosome_set:
                        if process:                           
                            if sign == '-':
                                psite = stop - 13
                                p_seq = n_sequence[-13]
                            else:
                                psite = start + 12
                                p_seq = n_sequence[12]
                                
                            mapq = line.split('\t')[4]
                            mid = str(line.split('\t')[6:9]).replace('[','').replace(']','').replace(',','\t').replace("'",'').replace(' ','')
                            qual = str(line.split('\t')[10:]).replace('[','').replace(']','').replace(',','\t').replace("'",'').replace(' ','')

                            new_hits = assign_reads(uid, chromo, psite, psite, sign, 'psites')
                            hit_ct += new_hits

                            new_line = ('{uid}\t{flag}\t{chromo}\t{psite}\t{mapq}\t1M\t{mid}\t{p_seq}\t{qual}\n').format(uid=uid, flag=flag, chromo=chromo,
                                       psite=psite, mapq=mapq, mid=mid, p_seq=p_seq, qual=qual)

                            sample_bam_dict['psites'][uid]= new_line
                    
    sam_file.close()
    sample_bam_dict['header'] = header_list       
        
    return()
        
def output_bam(output_dir, header_list, new_bam_dict):
    file_name = ('{}.sam').format(output_dir)
    new_sam_file = open(file_name, 'w')
    
    for header in header_list:
        new_sam_file.write(header)
        
    for uid, line in new_bam_dict.items():
        if uid != 'header':
            new_sam_file.write(line)
    new_sam_file.close()
    
    print('\tConverting to bam file\n') 
    bashCommand = ('samtools view -Sb {sample_name}.sam > {sample_name}_unsorted.bam').format(sample_name=args.output_file)
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    bashCommand = ('samtools sort -o {sample_name}.bam {sample_name}_unsorted.bam').format(sample_name=args.output_file)
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    bashCommand = ('samtools index {sample_name}.bam').format(sample_name=args.output_file)
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
def output_bed(name, atis_id, tis_scores, atis_sample_dict, outputfile):
        
    for each_sample, atis_id_dict in atis_sample_dict.items():       
        score = np.median(tis_scores)
        
        if atis_id in atis_id_dict:
            check_name = atis_id_dict[atis_id]['meta']['name']
            chromo = atis_id_dict[atis_id]['meta']['chromo']
            sign =  atis_id_dict[atis_id]['meta']['sign']
            seq = atis_id_dict[atis_id]['meta']['seq']
            triplet = atis_id_dict[atis_id]['meta']['triplet']
            start, stop = atis_id_dict[atis_id]['sr']
                        
            if name != check_name:
                print('atis name disagreement', name, check_name, atis_id)
                1/0
            
            else:
                outline = ('{chromo}\t{start}\t{stop}\t{name}_{triplet}_{seq}_nte\t{score}\t{sign}\n').format(chromo=chromo, start=start, stop=stop, name=name, seq=seq, triplet=triplet, score=score, sign=sign)
                outputfile.write(outline)
                return()
                
    print('not found', name, atis_id, score)    
    
def parse_samples(sample_list):
    name_dict = {}
    
    if sample_list:
        if len(sample_list)%3!= 0:
            print('Each sample requires a Name, RPF bam file, and RNA bam file.')
        else:
            for i in range(int(len(sample_list)/3)):
                name_dict[sample_list[i*3]]=[sample_list[(i*3)+1], sample_list[(i*3)+2]]
                
    return(name_dict)
    
def eval_atis(atis_id):
    global atis_id_dict
    global score_dict
    global eval_atis_dict
    global quantified_search_regions_dict
    
    process_ct = 0
    
    for each_sample in args.sample_list:
        if atis_id in atis_id_dict[each_sample]:
            if atis_id_dict[each_sample][atis_id]['meta']['region'] == 'tl':
                name = atis_id_dict[each_sample][atis_id]['meta']['name']

                if name in quantified_search_regions_dict[each_sample]:
                    if atis_id in quantified_search_regions_dict[each_sample][name]:
                        if quantified_search_regions_dict[each_sample][name][atis_id]['sr']['psites'][0] >= minimum_reads:
                            process_ct += 1
                            
                        if (process_ct/float(len(args.sample_list)) > 0.5) and process_ct >= 2:
                                                    
                            for each_sample in args.sample_list:                                
                                atis_details = quantified_search_regions_dict[each_sample][name][atis_id]
                                
                                if atis_id not in eval_atis_dict:
                                    eval_atis_dict[atis_id] = {'rte':[], 'rsm':[], 'tis_score':[], 'pval':[], 'atis_psites':[], 'sr_psites':[], 'up_psites':[]}
                                
                                rte = (atis_details['sr']['RPF'])/float(max(1, atis_details['sr']['RNA']))
                                rsm = (atis_details['atis']['psites'][0])/float(max(1, sum(atis_details['up']['psites'])))
                                weight = (atis_details['atis']['psites'][0]/float(atis_details['full']['psites'][0]))
                                eval_atis_dict[atis_id]['tis_score'].append(rte * rsm * weight)                      
                                eval_atis_dict[atis_id]['pval'].append(stats.binom_test([atis_details['sr']['psites'][0], atis_details['sr']['psites'][1] + atis_details['sr']['psites'][2]], p=athird))

                            if np.median(eval_atis_dict[atis_id]['pval']) <= 0.05 and ((atis_details['sr']['psites'][0] >= atis_details['sr']['psites'][1]) and (atis_details['sr']['psites'][0] >= atis_details['sr']['psites'][2])):
                                
                                output_bed(name, atis_id, eval_atis_dict[atis_id]['tis_score'], atis_id_dict, nte_potential_file)
                                
                                if name not in score_dict:
                                    score_dict[name] = {}
                                
                                if atis_id not in score_dict[name]:
                                    score_dict[name][atis_id] = (eval_atis_dict[atis_id]['tis_score'])
                                else:
                                    print('atis_id error', atis_id)
                                    1/0
                
    
if __name__ == '__main__':
    atis_id_dict = {}
    filter_nt = {}
    
    os_mkdir(args.output_file)
    
    if args.load_reads:
        print('Starting nteseq ... ')
        name_dict = parse_samples(args.sample_list)
        
        #2 parse fasta
        print('Parsing fasta file... ')
        fasta_dict = parse_fasta(args.fa_file)

        #1 parse_gff
        print('Parsing GFF file for genes with transcript leaders...')
        coord_dict, chromosome_set = parse_gff(gff_file, gene_tag, tl_tag, tp_tag, min_tl, min_3p, mask_tl, mask_3p)
                
        #3 for each gene get fasta of upstream, scan for first stop codon, all start codons
        print('Defining NTE search regions ...')
        search_region_dict, assign_region_dict, flanking_region_dict = build_search_region(coord_dict, fasta_dict, dsrl) 

        #4 load bam, sam file
        print('Loading reads from bam/sam files into candidate regions... ')
        sample_dict = {'RNA':{'sr':{}, 'fl':{}, 'uid':{}}, 'RPF':{'sr':{}, 'fl':{}, 'uid':{}}, 'psites':{'sr':{}, 'fl':{}, 'uid':{}}}
        sample_bam_dict = {'RNA':{}, 'RPF':{}, 'psites':{}, 'header':[]}
        
        #for i in range(len(name_dict)+1):
        jobs = []
        for each_sample, RPF_RNA_pair in name_dict.items():    
            print('4... ', each_sample, datetime.datetime.now())
            
            RPF_name, RNA_name = RPF_RNA_pair
            
            convert_rpf_name = convert_to_sam(RPF_name)
            convert_rna_name = convert_to_sam(RNA_name)
        
            load_reads(convert_rpf_name, chromosome_set, search_region_dict, assign_region_dict, flanking_region_dict, 'RPF')
            load_reads(convert_rna_name, chromosome_set, search_region_dict, assign_region_dict, flanking_region_dict, 'RNA')
            
            #TODO: future versions multiprocess read loading
            #p = multiprocessing.Process(target=load_reads, args=(convert_rpf_name, chromosome_set, search_region_dict, assign_region_dict, flanking_region_dict, 'RPF',))
            #jobs.append(p)
            #p = multiprocessing.Process(target=load_reads, args=(convert_rna_name, chromosome_set, search_region_dict, assign_region_dict, flanking_region_dict, 'RNA',))
            #jobs.append(p)
            
            #p.start()
            #p.join()
                
        #  
        print('generating P-site bam files ... ')
        output_bam(args.output_file, sample_bam_dict['header'], sample_bam_dict['psites'])
        
        resource_pickle_name = ('{}_sample_dict.p').format(args.output_file)    
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(sample_dict, file)
        
        #
        print('loading reads and p-sites into candidate regions ... ')
        quantified_search_regions_dict = {}
        
        rep_psite_dict = sample_dict['psites']['sr']
        rep_RPF_dict = sample_dict['RPF']['sr']
        rep_RNA_dict = sample_dict['RNA']['sr']
        
        for atis_id, etc in atis_id_dict.items():
            name = etc['meta']['name']
            sign = etc['meta']['sign']
            
            if sign == '+':
                gene_start_codon = etc['gene'][0]
            else:
                gene_start_codon = etc['gene'][1]
                
            if name in rep_psite_dict:
                sr_psites_dict = rep_psite_dict[name]
                sr_RPF_dict = rep_psite_dict[name]
                sr_RNA_dict = rep_psite_dict[name]
                                        
                if name not in quantified_search_regions_dict:
                    quantified_search_regions_dict[name] = {}
                    
                if atis_id not in quantified_search_regions_dict[name]:
                    quantified_search_regions_dict[name][atis_id] = {'atis':{}, 'sr': {}, 'up':{}, 'full':{}} 
                    for atis_region in ['atis', 'sr', 'up', 'full']:
                        quantified_search_regions_dict[name][atis_id][atis_region] = {'psites':{}, 'RNA':0, 'RPF':0}
                        quantified_search_regions_dict[name][atis_id][atis_region]['psites'] = [0,0,0]

                        for nt in range(etc[atis_region][0], etc[atis_region][1]+1):
                            if nt in sr_psites_dict:
                                psite_ct = len(sr_psites_dict[nt])
                                rpf_ct = len(sr_RPF_dict[nt])
                                rna_ct = len(sr_RNA_dict[nt])
                                
                                if sign == '+':
                                    reading_frame = (nt - gene_start_codon - 1 ) % 3
                                else:
                                    reading_frame = (gene_start_codon - nt) % 3

                                quantified_search_regions_dict[name][atis_id][atis_region]['psites'][reading_frame] += psite_ct
                                quantified_search_regions_dict[name][atis_id][atis_region]['RPF'] += rpf_ct
                                quantified_search_regions_dict[name][atis_id][atis_region]['RNA'] += rna_ct
                                        
        resource_pickle_name = ('{}_quantified_search_regions_dict.p').format(args.output_file)
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(quantified_search_regions_dict, file)
            
        resource_pickle_name = ('{}_atis_id_dict.p').format(args.output_file)    
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(atis_id_dict, file)
        
    if args.evaluate:
        union_set = set()
        quantified_search_regions_dict = {'union':{}}
        sample_dict = {}
        atis_id_dict = {}
        
        for each_sample in args.sample_list:             

            pickle_out = ('{}_quantified_search_regions_dict.p').format(each_sample)
            quantified_search_regions_dict[each_sample] = pickle.load(open(pickle_out))

            pickle_out = ('{}_sample_dict.p').format(each_sample)
            sample_dict[each_sample] = pickle.load(open(pickle_out))
        
            pickle_out = ('{}_atis_id_dict.p').format(each_sample)
            atis_id_dict[each_sample] = pickle.load(open(pickle_out))
             
            for atis_id, _etc in atis_id_dict[each_sample].items():
                union_set.add(atis_id)
                                                                  
        athird = float(1/3)
                
        eval_atis_dict = {}
        score_dict = {}
        
        nte_candidate_file_name = ('{}.bed').format(args.output_file)
        nte_candidate_file = open(nte_candidate_file_name, 'w')
        
        nte_potential_file_name = ('{}_potential.bed').format(args.output_file)
        nte_potential_file = open(nte_potential_file_name, 'w')
        
        print('Evaluating candidate NTE events... ')
        
        for atis_id in union_set:
            eval_atis(atis_id)
                        
        for name, atis_scores in score_dict.items():
            best_score = 0
            best_set = [0,0,0]
            best_atis = ''
                
            for atis_id, atis_score in atis_scores.items():
                if len(atis_score) > 0:
                    calc_score = sum(atis_score)
                    if best_atis != atis_id:
                        if calc_score > best_score:
                            best_score = calc_score
                            best_atis = atis_id
                            best_set = atis_score
                            
                        if calc_score == best_score:
                             if np.median(atis_score) > np.median(best_set):
                                 best_score = calc_score
                                 best_atis = atis_id
                                 best_set = atis_score
                                                            
            if best_atis != '':                      
                output_bed(name, best_atis, best_score, atis_id_dict, nte_candidate_file)
                                    
        nte_potential_file.close()
        nte_candidate_file.close()