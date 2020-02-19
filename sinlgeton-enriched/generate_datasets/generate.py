#!/usr/bin/python3

import random
import math
import copy
import sys
from statistics import *
from Bio import pairwise2


#from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
#import plotly.plotly as py
#import plotly.graph_objs as go
#init_notebook_mode(connected=True)


#nof_genes = 1000
#nof_genomes = 2000
nof_genomes = 1000

#min_gene_len = 50
#avg_gene_len = 350
#max_gene_len = 1200

gene_variation_prob = 0.8	#probability of variation when ancestor gene is aquired
#gene_variation_percentage = 0.005	#variation ammount as percentage of gene length

locus_variation_prob = 0.05

#locus_variation_add = 0.33
#locus_variation_remove = 0.66
#probability to variate a locus wihtin a gene
gene_duplication_prob = 0.001	#probability of duplicate a gene
#gene_duplication_prob = 0.0
geneset_variation = 0.01	#percentage of variation in gene sets, it includes creation of new genes and removal of inherited ones
geneset_variation_add = 1.0
geneset_variation_remove = 0.1

#geneset_variation = 0.05	#percentage of variation in gene sets, it includes creation of new genes and removal of inherited ones
#geneset_variation_add = 0.05
#geneset_variation_remove = 1.0


alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
def generate_random_gene(length):
	s = ""
	for j in range(0,length):
		s += random.choice(alphabet)
	return s

def variate_gene_byperc(seq, variation_percentage):
	nof_variations = math.ceil(seq * variation_percentage)
	for v in range(nof_variation):
		position = random.randrange( len(seq)+1 )
		operation = random.choice( [0,1,2] )
		if operation == 0: #snp
			seq = seq[0:position] + random.choice(alphabet) + seq[position+1:]
		elif operation == 2: #delete
			seq = seq[0:position] + seq[position+1:]
		else: #insert
			seq = seq[0:position] + random.choice(alphabet) + seq[position:]
	return seq

def variate_gene_byprob(seq, variation_probability):
	for position in range(len(seq)):
		if random.random() <= variation_probability:
			operation = random.choice( [0,1,2] )
			if operation == 0: #snp
				seq = seq[0:position] + random.choice(alphabet) + seq[position+1:]
			elif operation == 2: #delete
				seq = seq[0:position] + seq[position+1:]
			else: #insert
				seq = seq[0:position] + random.choice(alphabet) + seq[position:]
	return seq


def print_ancestors_d3js(genome_ancestors):
	for g in sorted(genome_ancestors.keys()):
		#print(g)
		s = str(g) + ",1"
		a = g
		while( genome_ancestors[a] != 0 ):
			a = genome_ancestors[a]
			s = str(a) +"."+s
		s = "0."+s
		print(s)


genome_ancestors = {0:0}
genomes = dict()
genomes[0] = list()
gene_count = 0


ifile = sys.argv[1]

locus_variation_prob = float(sys.argv[3])

i = 0
for line in open(ifile, 'r'):
	if i%2 == 0:
		pass
	else:
		seq = line.strip()
		genomes[0].append( (gene_count, seq) )
		gene_count += 1
	i += 1

min_gene_len = min( [ len(s[1]) for s in genomes[0] ] )
avg_gene_len = int(mean( [ len(s[1]) for s in genomes[0] ] ))
max_gene_len = max( [ len(s[1]) for s in genomes[0] ] )

print("nof genes", len(genomes[0]))
print("min length", min_gene_len)
print("max length", max_gene_len)
print("mean length", avg_gene_len)



for g in range(1, nof_genomes):
	#select ancestor genome from previous created ones
	ag = random.randrange(g) #ancestor genome index
	genome_ancestors[g] = ag
	print('parent',g,ag)
	#copy genes from ancestor genome and variate them
	genomes[g] = list()
	for ggi in range(len(genomes[ag])):
		#inheritance
		gg_id = genomes[ag][ggi][0]
		gg_seq = copy.deepcopy( genomes[ag][ggi][1] )
		#variation
		if random.random() <= gene_variation_prob:
			gg_seq = variate_gene_byprob(gg_seq, locus_variation_prob)
		genomes[g].append(  (gg_id, gg_seq) )
		#duplication
		if random.random() <= gene_duplication_prob:
			genomes[g].append(  (gg_id, gg_seq) )
	#add/remove genes
	geneset_size = len(genomes[g])
	for i in range( math.ceil(geneset_size * geneset_variation)  ):
		if random.random() <= geneset_variation_remove:
			#remove a gene
			p = random.randrange( geneset_size )
			genomes[g] = genomes[g][:p] + genomes[g][p+1:]
			geneset_size -= 1
		#else:
		if random.random() <= geneset_variation_add:
			#add a totally new gene
			gene_count += 1
			l = math.ceil(random.triangular(min_gene_len, max_gene_len, avg_gene_len))
			genomes[g].append( (gene_count,generate_random_gene(l))  )
	#shuffle genes
	random.shuffle(genomes[g])


seqs = dict()
seq_count = dict()

with open(sys.argv[2], 'w') as ofile:
	for g,genes in genomes.items():
		para_count = dict()
		for gene in genes:
			seqs[ gene[0] ] = gene[1]
			seq_count[ gene[0] ] = seq_count.get(gene[0],0)+1
			para_count[gene[0]] = para_count.get(gene[0], 0) + 1
			pc = para_count[gene[0]]
			ofile.write('genome_'+ str(g)+'\tsequence_'+ str(gene[0])+'_'+ str(g) +'_'+str(pc) +'\tsequence_'+ str(gene[0]) +'\n')
			ofile.write(gene[1] +'\n')
		
