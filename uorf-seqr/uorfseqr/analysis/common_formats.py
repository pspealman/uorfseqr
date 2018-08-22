#Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman

import urllib


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

def import_bed_smpl(fname, per_chromosome):
    z = open(fname, 'r')
    print 'importing bed file', fname
    for zline in z:
        line = zline.rstrip()
        if ( len(line) == 0 ):
            continue
        if ( ( 'track_name' in line ) or ( '#' in line[0] ) ):
            continue
        Line = line.split()
        if ( Line[0], Line[5] ) not in per_chromosome:
            per_chromosome[( Line[0], Line[5] )] = set([])
        per_chromosome[( Line[0], Line[5] )].add((int(Line[1]), int(Line[2]), Line[3]))
    return per_chromosome


def import_gff_file(fname, per_chromosome_TLS_data, per_chromosome_features_masked_out, unwanted_genes):
	print "importing GFF data", fname
	for line in open(fname, 'r'):
		if '#' in line[0]:
			continue
		if '##FASTA' in line:
			break
		Line = line.rstrip().split("\t")
		if len(Line) < 8:
			break
		seqname = Line[0]
		source = Line[1]
		feature = Line[2]
		start = int(Line[3])-1
		end = int(Line[4])
		score = Line[5]
		strand = Line[6]
		frame = Line[7]
		attribute = Line[8]
		attribute_components = {x:y for (x,y) in [z.split('=') for z in attribute.split(';')]} 
		if 'gene' in feature:
			myname = urllib.unquote(attribute_components['ID'])
			if 'orf_classification' in attribute_components:
				if 'Dubious' in attribute_components['orf_classification']:
					unwanted_genes.add(myname)
			if 'gene' != feature:
				unwanted_genes.add(myname)
			if (seqname, strand) not in per_chromosome_features_masked_out:
				per_chromosome_features_masked_out[(seqname, strand)] = {}
			if '+' in strand:
				per_chromosome_features_masked_out[(seqname, strand)][myname] = ( start, end )
			else:
				per_chromosome_features_masked_out[(seqname, strand)][myname] = ( end, start )
		# elif 'intron' in feature:
		# 	unwanted_genes.add(urllib.unquote('_'.join(attribute_components['Parent'].split('_')[:-1])))
		elif 'blocked_reading_frame' in feature:
			myname = urllib.unquote(attribute_components['ID'])
			unwanted_genes.add(myname)
		elif 'transcription_start_site' in feature:
			myname = urllib.unquote(attribute_components['Parent'])
			if myname not in unwanted_genes:
				if (seqname, strand) not in per_chromosome_TLS_data:
					per_chromosome_TLS_data[(seqname, strand)] = {}
				me = per_chromosome_TLS_data[(seqname, strand)].setdefault(myname,start)
				if '+' in strand:
					per_chromosome_TLS_data[(seqname, strand)][myname] = min(per_chromosome_TLS_data[(seqname, strand)][myname], start)
				else:
					per_chromosome_TLS_data[(seqname, strand)][myname] = max(per_chromosome_TLS_data[(seqname, strand)][myname], start)





class gff_annotation:
	def __init__(self):
		self.chromosome = None
		self.strand = None
		self.frame = None
		self.TSSs = []
		self.TIS = None
		self.STOP = None
		self.pAs = []
		self.introns = None

def import_gff_file_extended(fname, genome_annotations, unwanted_genes):
	print "importing GFF data", fname
	for line in open(fname, 'r'):
		if '#' in line[0]:
			continue
		if '##FASTA' in line:
			break
		Line = line.rstrip().split("\t")
		if len(Line) < 8:
			break
		seqname = Line[0]
		source = Line[1]
		feature = Line[2]
		start = int(Line[3])-1
		end = int(Line[4])
		score = Line[5]
		strand = Line[6]
		frame = Line[7]
		attribute = Line[8]
		attribute_components = {x:y for (x,y) in [z.split('=') for z in attribute.split(';')]} 
		if 'gene' in feature:
			myname = urllib.unquote(attribute_components['ID'])
			if 'orf_classification' in attribute_components:
				if 'Dubious' in attribute_components['orf_classification']:
					unwanted_genes.add(myname)
			if 'gene' != feature:
				unwanted_genes.add(myname)
			if myname not in genome_annotations:
				genome_annotations[myname] = gff_annotation()
			if '+' == strand:
				genome_annotations[myname].TIS = start
				genome_annotations[myname].STOP = end
			else:
				genome_annotations[myname].TIS = end
				genome_annotations[myname].STOP = start
			genome_annotations[myname].strand = strand
			genome_annotations[myname].chromosome = seqname
		elif 'blocked_reading_frame' in feature:
			myname = urllib.unquote(attribute_components['ID'])
			unwanted_genes.add(myname)
		elif 'transcription_start_site' in feature:
			myname = urllib.unquote(attribute_components['Parent'])
			if myname not in unwanted_genes:
				if myname not in genome_annotations:
					genome_annotations[myname] = gff_annotation()
				try:
					genome_annotations[myname].TSSs.append((start,float(score)))
				except ValueError:
					genome_annotations[myname].TSSs.append((start,0))
		elif 'polyA_site' in feature:
			myname = urllib.unquote(attribute_components['Parent'])
			if myname not in unwanted_genes:
				if myname not in genome_annotations:
					genome_annotations[myname] = gff_annotation()
				try:
					genome_annotations[myname].pAs.append((start,float(score)))
				except ValueError:
					genome_annotations[myname].pAs.append((start,0))
