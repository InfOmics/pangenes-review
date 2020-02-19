#!/usr/bin/pythpon3

import sys

ifile = sys.argv[1]

families = dict()



lineno = 0
for line in open(ifile,'r'):
	if lineno % 2 == 0:
		#print(line)
		cc = line.strip().split('\t')
		#print(cc)
		family = cc[2]
		gene = cc[2].replace('sequence_','')
		rep = cc[1].split('_')[3]
		genome = cc[0].replace('genome_','')
		#G108_SE043
		geneacc = "G"+gene+":"+rep+"_SE"+genome
		#print(geneacc)
		if family not in families:
			families[family] = set()
		families[family].add(geneacc)
	lineno += 1


for k,v in families.items():
	print(' '.join(v))
