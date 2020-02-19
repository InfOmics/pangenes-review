#!/usr/bin/ptyhon3

import sys
import os
import random

ifile = sys.argv[1]
odir = sys.argv[2]


rna2dna_maps = {
'A' : ['GCT', 'GCC', 'GCA', 'GCG'],
'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
'N' : ['AAT', 'AAC'],
'D' : ['GAT', 'GAC'],
'C' : ['TGT', 'TGC'],
'Q' : ['CAA', 'CAG'],
'E' : ['GAA', 'GAG'],
'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
'H' : ['CAT', 'CAC'],
'I' : ['ATT', 'ATC', 'ATA'],
'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
'K' : ['AAA', 'AAG'],
'M' : ['ATG'],
'F' : ['TTT', 'TTC'],
'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
'W' : ['TGG'],
'Y' : ['TAT', 'TAC'],
'V' : ['GTT', 'GTC', 'GTA', 'GTG'],
'START' : ['ATG'],
'STOP' : ['TAA', 'TGA', 'TAG']
}

def rna2dna(s):
	d = ''+ random.choice(rna2dna_maps['START'])
	for c in s:
		d += random.choice(rna2dna_maps[c])
	d += random.choice(rna2dna_maps['STOP'])
	return d


igenomes = dict()
with open(ifile,'r') as iff:
	lines = iff.readlines()
	for i in range(0,len(lines), 2):
		cc = lines[i].strip().split('\t')
		gname = cc[0]
		sname = cc[1]
		sdescr = cc[2]
		print(i,cc)
		if gname not in igenomes:
			igenomes[gname] = list()
		igenomes[gname].append( (sname, sdescr, lines[i+1].strip() ) )


def write_80cols(s, off):
	for i in range(0,len(s),80):
		off.write( s[i: i+80 if i+80<len(s) else len(s)] +'\n')

try:
	os.stat(odir)
except:
	os.mkdir(odir)

for gname,seqs in sorted(igenomes.items()):
	print(gname)
	ofile = odir+"/"+gname+".gff"
	print(ofile)

	gseq = ''
	scoords = list()
	for seqi in seqs:
		scoords.append(  (len(gseq)+3, len(gseq)+3 + (len(seqi[2])*3 ))  )
		gseq += rna2dna(seqi[2])
		
		#print(seqi[2])
		#print(rna2dna(seqi[2]))

	
	off = open(ofile,'w')
	off.write('##gff-version 3\n')
	off.write('##sequence-region '+gname+' 1 '+ str(len(gseq))+'\n')
	off.write('# conversion-by db2gff3.py\n')
	off.write(gname+'	Dicupan	region	1	'+str(len(gseq))+'	.	+	1	ID='+gname+'\n')

	ind = 0
	for seqi in seqs:
		coords = scoords[ind]
		off.write(gname+'	Dicupan	gene	'+str(coords[0])+'	'+str(coords[1])+'	.	+	1	ID='+seqi[0]+';Name='+seqi[0]+'\n')
		off.write(gname+'	Dicupan	mRNA	'+str(coords[0])+'	'+str(coords[1])+'	.	+	1	ID='+seqi[0]+'.t01;Parent='+seqi[0]+'\n')
		off.write(gname+'	Dicupan	CDS	'+str(coords[0])+'	'+str(coords[1])+'	.	+	1	ID='+seqi[0]+'.p01;Parent='+seqi[0]+'.t01;Name='+seqi[0]+';product='+seqi[1]+'\n')
		off.write(gname+'	Dicupan	exon	'+str(coords[0])+'	'+str(coords[1])+'	.	+	1	Parent='+seqi[0]+'.t01\n')
		ind += 1

	off.write('##FASTA\n')
	off.write('>' + gname +'\n')
	write_80cols(gseq, off)

	for seqi in seqs:
		off.write('>' + seqi[0] +'.p01\n')
		write_80cols(seqi[2],off)

	off.flush()
	off.close()


