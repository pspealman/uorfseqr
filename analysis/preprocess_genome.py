#Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman
#PS - 08.07.18: updated syntax for print commands
import string
import numpy
import sys
import os

genome = []
cache = []
commands = open(sys.argv[1],'r')
for rline in commands:
    line = rline.rstrip()
    Line = string.split(line)
    if '#' in Line[0]:
        continue
    if 'genome' in Line[0]:
        genome.append(Line[1])
    if 'cachedir' in Line[0]:
        cache.append(Line[1])

if len(genome) < 1:
    print("you must define a genome to process")
    exit(1)

if len(cache) == 0:
    cache = ''
else:
    cache = cache[0]
    if cache[-1] != '/':
        cache = cache + '/'

if not os.path.exists(cache):
    os.mkdir(cache)

base_index = {}
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



for which_genome in range(len(genome)):
    z = open(genome[which_genome], 'r')

    chrom = ([], [], [], [])
    chrom_size = {}
    def flush_chromosome(chromosome_name, chrom):
        outv = numpy.zeros((len(chrom[0]),4), dtype='bool')
        for i in range(4):
            outv[:,i] = chrom[i]
        numpy.save(cache+chromosome_name+'-positive_strand', outv)
        for i in range(4):
            chrom[i].reverse()
        outv[:,0] = chrom[1]
        outv[:,1] = chrom[0]
        outv[:,2] = chrom[3]
        outv[:,3] = chrom[2]
        numpy.save(cache+chromosome_name+'-negative_strand', outv)

    for zline in z:
        line = zline.rstrip()
        if line[0] == '>':
            if len(chrom[0]) > 0:
                flush_chromosome(chromosome_name, chrom)
                chrom_size[chromosome_name] = len(chrom[0])
            Line = string.split(line[1:])
            chromosome_name = Line[0]
            print('handling chromosome', chromosome_name)
            chrom = ([], [], [], [])
        else:
            for ii,i in enumerate(line):
                if i not in base_index:
                    print(ii, line[ii-20:ii+1])
                updates = base_index[i]
                for j in range(4):
                    chrom[j].append(updates[j])

    flush_chromosome(chromosome_name, chrom)
    chrom_size[chromosome_name] = len(chrom[0])

