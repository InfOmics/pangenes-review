#!/usr/bin/python3


import sys
from os import listdir
from os.path import isfile, join
import re
from Bio import SeqIO


idir = sys.argv[1]
odir = sys.argv[2]
print("reading gbk files from", idir)

fastafiles = [f for f in listdir(idir) if isfile(join(idir, f)) and re.match('^.+\.fasta$', f)]
print(fastafiles)

for fastafile in fastafiles:
	newname = fastafile.replace('_aa.fasta','@').replace('_','').replace('@','_aa.fasta')
	print(newname)
	suff = newname.replace('_aa.fasta','')
	iname = -1
	with open(odir+'/'+newname, 'w') as off:
		for line in open(idir+'/'+fastafile):
			if line[0] == '>':
				#line = line.strip().replace('_','').replìace(':','').replace('.','').replace('-','') + ":1_"+suff+"\n"
				iname += 1
				line = ">"+ str(iname)+ ":1_"+suff+"\n"
			off.write(line)
