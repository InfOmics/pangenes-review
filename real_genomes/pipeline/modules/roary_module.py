#!/usr/bin/ptyhon3
from pathlib import Path
import sys
import os
import random
from subprocess import call
import re
from Bio import SeqIO
from Bio.Seq import Seq
from tempfile import NamedTemporaryFile

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

def write_80cols(s, off):
	for i in range(0,len(s),80):
		off.write( s[i: i+80 if i+80<len(s) else len(s)] +'\n')

def rna2dna(s):
	d = ''+ random.choice(rna2dna_maps['START'])
	for c in s:
		d += random.choice(rna2dna_maps[c])
	d += random.choice(rna2dna_maps['STOP'])
	return d


def createRoaryInput(in_path, out_path):
	ifile = in_path
	print(ifile)
	tmp_file = NamedTemporaryFile()
	aux = open(tmp_file.name,'w')
	for i,line in enumerate(open(ifile,'r')):
		if i%2 ==0:
			aux.write(line.replace(':','_'))
		else:
			aux.write(line)
	aux.close()
	""" for i,l in enumerate(open(tmp_file.name,'r')):
		print(i,l)
		if i == 5:
			quit()
	tmp_file.close()"""
	""" for l in open(tmp_file.name,'r'):
		print(l) """
	ifile = tmp_file.name
	odir = Path(out_path,'input')
	odir.mkdir(parents=True,exist_ok=True)
	#odir = str(odir)

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


	try:
		os.stat(odir)
	except:
		os.mkdir(odir)

	for gname,seqs in sorted(igenomes.items()):
		print(gname)
		ofile = str(odir)+"/"+gname+".gff"
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
	#	off.write(gname+'	Dicupan	region	1	'+str(len(gseq))+'	.	+	1	ID='+gname+'\n')

		ind = 0
		for seqi in seqs:
			coords = scoords[ind]
	#		off.write(gname+'	Dicupan	gene	'+str(coords[0])+'	'+str(coords[1])+'	.	+	1	ID='+seqi[0]+';Name='+seqi[0]+'\n')
	#		off.write(gname+'	Dicupan	mRNA	'+str(coords[0])+'	'+str(coords[1])+'	.	+	1	ID='+seqi[0]+'.t01;Parent='+seqi[0]+'\n')
	#		off.write(gname+'	Dicupan	CDS	'+str(coords[0])+'	'+str(coords[1])+'	.	+	1	ID='+seqi[0]+'.p01;Parent='+seqi[0]+'.t01;Name='+seqi[0]+';product='+seqi[1]+'\n')
	#		off.write(gname+'	Dicupan	exon	'+str(coords[0])+'	'+str(coords[1])+'	.	+	1	Parent='+seqi[0]+'.t01\n')
			#off.write(gname+'	Dicupan	CDS	'+str(coords[0] + 1)+'	'+str(coords[1])+'	.	+	0	ID='+seqi[1]+'_'+str(ind)+';gene='+seqi[0]+';inference=synthetic;locus_tag='+seqi[1]+'_'+str(ind)+';product='+seqi[1]+'\n')
			off.write(gname+'	Dicupan	CDS	'+str(coords[0] + 1)+'	'+str(coords[1])+'	.	+	0	ID='+seqi[0]+'_'+str(ind)+';gene='+seqi[0]+';inference=synthetic;locus_tag='+seqi[0]+'_'+str(ind)+';product='+seqi[1]+'\n')


			ind += 1

		off.write('##FASTA\n')
		off.write('>' + gname +'\n')
		write_80cols(gseq, off)

	#	for seqi in seqs:
	#		off.write('>' + seqi[0] +'.p01\n')
	#		write_80cols(seqi[2],off)

		off.flush()
		off.close()
	tmp_file.close()

	
from modules.resource_control import call_program
def callRoary(path_roary_data):
	
	path_output = Path(path_roary_data,'output')
	#cores = 6
	#roary_software = Path('pangenome_softwares','roary','bin','roary')
	
	#-p 6 ( 6 threads)
	#-f output directory
	#-b blast executable [blastp]
	#-c mcl executable [mcl]
	#-m makeblastdb executable
	#-v verbose
	
	gffs = [i for i in Path(path_roary_data,'input').glob('*')]
	#forse posso anche definire di utilizzare blast,mcl e makeblastdb dalla cartella in comune
	#cmd = ['roary', '-v', '-p', str(cores), '-f', path_output ]+gffs
	parameters = [ path_output ] + gffs
	stat = call_program(parameters,'roary')
	return stat 
	
def roary2families(path_roary_data, path_Gfamily):
	roary_file = Path(path_roary_data,'output','clustered_proteins')
	
	if roary_file.exists():
		with open(Path(path_Gfamily,'roary_families.clus'),'w') as out:
			for line in open(roary_file,'r'):
				#G6341_1_SE713: G6341_1_SE713_188

				cc = (line.strip().split(' ')[1]).split('\t')
				buffer = list()
				for v in cc:
					vc = v.split('_')
					buffer.append( vc[0]+":"+vc[1]+"_"+vc[2] )
				#aux = re.sub(r'_[0-9]+\t',' ',re.sub(r'^.+\:','',line).strip())+'\n'
				#buffer = list()
				#for gene in aux.rstrip().split(' '):
				#	aux_gene = gene.split('_')
				#	buffer.append(aux_gene[0]+':'+aux_gene[1]+'_'+aux_gene[2])
				out.write(' '.join(buffer)+'\n')
	else:
		print('MESSAGE:\n File "'+str(roary_file)+'" does not exist, families will not be computed (file.clus)')
	
