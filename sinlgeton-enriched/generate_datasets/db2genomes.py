#!/usr/bin/pythpon3

import sys

ifile = sys.argv[1]
ofdir = sys.argv[2]


genome_sequences = dict()



geneacc = ""
genome = ""

lineno = 0
for line in open(ifile,'r'):
	if lineno % 2 == 0:
		cc = line.strip().split('\t')
		family = cc[2]
		gene = cc[2].replace('sequence_','')
		rep = cc[1].split('_')[3]
		genome = cc[0].replace('genome_','')
		geneacc = "G"+gene+":"+rep+"_SE"+genome
	else:
		if genome not in genome_sequences:
			genome_sequences[genome] = ""
		genome_sequences[genome] += ">"+geneacc+"\n"
		genome_sequences[genome] += line
	lineno += 1


for k,v in genome_sequences.items():
	with open(ofdir+"/SE"+k+"_aa.fasta",'w') as of:
		of.write(v)
