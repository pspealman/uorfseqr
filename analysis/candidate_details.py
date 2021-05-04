# -*- coding: utf-8 -*-
"""
Created on Mon May  3 17:11:52 2021

Small script to extract full

@author: pspealman
"""
import os
import argparse
import subprocess

### Argparser definitions        
parser = argparse.ArgumentParser()

'''
    python candidate_details.py -f test.fa \
        -i test.bed \
        -o test.output
'''
# basic arguments
parser.add_argument('-f',"--fasta")
parser.add_argument('-i',"--input_candidate_uORF_bed")
parser.add_argument('-o',"--output_file")
parser.add_argument('-manual',"--manual",action='store_true')
parser.add_argument('-full',"--full_sequence",action='store_true')

args = parser.parse_args()

fasta_file_name = args.fasta
input_file_name = args.input_candidate_uORF_bed
output_file_name = args.output_file

def output_handler(output):
    if len(output.strip()) > 1:
        print(output)

def help_dialog():
    monolog = ('\tThis script uses bedtools getfasta command to collect seqeunces '
               'from provided fasta file based on the coordinates found in the '
               'input_uorf_bed file.\n')
    print(monolog)
    
    monolog = ('\t-full, --full_sequence [optional]. Using this command will also '
               'output the full candidate uORF sequence instead of just the start codon.\n')
    print(monolog)
    
    monolog = ('\t\tUsage:\n\tpython candidate_details.py -f <path_to_genome_fasta_file> ' 
               '-i <path_to_*candidate_uORF.bed_file> -o <path_for_output>\n')
    print(monolog)
    
    quit()
    
if args.manual:
    help_dialog()
    quit()
            
def os_mkdir(in_name):    
    if '/' in in_name:
        directory_name = in_name.rsplit('/',1)[0]
        
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)    

if __name__ == '__main__':    
    cuorf_dict = {}
    
    if args.output_file:
        os_mkdir(args.output_file)
    else:
        quit()
    
#    input_file = open(input_file_name)
#    
#    for line in input_file:
#        #chrXVI	593068	593110	YPR016C.593110	29     -
#        if line[0]!='#':
#            chromo, start, stop, cuorf_id, score, sign
#            if cuorf_id not in cuorf_dict:
#                cuorf_dict[cuorf_id] = {}

    
    print('\tConverting to bam file\n') 
    bashCommand = ('bedtools getfasta -fi {fasta} -bed {input_f} -s -name -tab -fo {output_f}_raw_results.log').format(
            fasta = fasta_file_name, input_f = input_file_name, output_f = output_file_name)
    
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    rawfile_name = ('{output_f}_raw_results.log').format(output_f = output_file_name)
    rawfile = open(rawfile_name)
    
    for line in rawfile:
        #forward::chr1:20-25(+)  CGCTA
        line = line.strip()
        deets, sequence = line.split('\t')
        
        cuorf_id, locus = deets.split('::')
        
        if cuorf_id not in cuorf_dict:
            cuorf_dict[cuorf_id] = {'locus': locus, 'sequence': sequence, 'TIS': sequence[0:3]}
        else:
            print('Error: Duplicate cuorf_id present in candidate_uORF.bed file: ', cuorf_id)

    rawfile.close()

    outfile = open(output_file_name, 'w')

    for cuorf_id in cuorf_dict:
        locus = cuorf_dict[cuorf_id]['locus']
        sequence = cuorf_dict[cuorf_id]['sequence']
        TIS = cuorf_dict[cuorf_id]['TIS']
        
        if args.full_sequence:
            outline = ('>{cuorf_id}, {locus}, full_sequence\n{sequence}\n').format(
                    cuorf_id = cuorf_id, locus = locus, sequence = sequence)
            outfile.write(outline)
        else:
            outline = ('>{cuorf_id}, {locus}, start_codon\n{TIS}\n').format(
                    cuorf_id = cuorf_id, locus = locus, TIS = TIS)
            outfile.write(outline)
        
    outfile.close()
    
    print('Summary completed.')
        
        
        
        
        
        
        
        
        

