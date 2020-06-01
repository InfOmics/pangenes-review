from pathlib import Path
import re
from shutil import copy, move
from Bio import SeqIO
from subprocess import call
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature , FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import datetime
from modules.resource_control import call_program

def createPanocInput(path_roary_data,path_panoct_data, in_path):
	
	#path_genomes = Path(path_alfdata,'genome','genome','DB')
	out = Path(path_panoct_data,'input') #the resulting files will be the input for panoct
	out.mkdir(exist_ok=True)
	in_gffsDir = Path(path_roary_data,'input')

	tags_file = Path(out,'exp.tags')
	gene_att_file = Path(out,'exp.gene_att')
	peps_file = Path(out,'exp.pep')
	dbname='exp_prot_db'
	db_path=Path(out,dbname,dbname)
	blast_output = Path(out,'blastall.out')

	#gffs files list
	gffs = sorted([i for i in in_gffsDir.glob('*.gff')])
	
	#output file of genomes names (exp.tags)
	with open(tags_file,'w') as t:
		for i in gffs:
			t.write(i.stem+'\n')

	
	#output file of gene attributes
	with open(gene_att_file,'w') as gattr:
		for gff in gffs:
			with open(gff,'r') as f:
				a = f.readlines()

				del a[0:3] #delete header gff
				for record in a:
					if re.match('^##FASTA',record):
						break
					aux = record.split('\t')
					
					contig_id = '1'
					genome_id = aux[0]
		
					start = aux[3]
					end = aux[4]

					aux2 = aux[-1].rstrip().split(';')
					
					annotation = ''#';'.join(aux2[2:]) #leave empty the annotation part because it's not needed
					
					aux_protein_id = aux2[1].split('=')[-1].split('_')
					protein_id = aux_protein_id[0]+':'+aux_protein_id[1]+'_'+aux_protein_id[2]
					
					gene_attr_record = '\t'.join([contig_id,protein_id,start,end,annotation,genome_id+'\n'])
					gattr.write(gene_attr_record)
	
	#output exp.pep contenente i fasta delle proteine
	# >gene_id
	# VMLAK....

	#copying genomes fasta
	Path(out,'fasta').mkdir(exist_ok=True)
	for i in Path(in_path,'protein').glob('*_aa.fasta'):
		copy(i,Path(out,'fasta'))
	
	fastas = sorted([i for i in Path(out,'fasta').glob('*_aa.fasta')])
	
	pep_list = list()
	for f in fastas:
		for protein in SeqIO.parse(f,'fasta'):
			header = protein.id.rstrip(',')
			sequence = Seq(str(protein.seq), IUPAC.protein)
			sR = SeqRecord(seq=sequence, id=header, description='')
			pep_list.append(sR)
			
	SeqIO.write(pep_list,peps_file,'fasta')


def callPanoct(in_path,path_data):
	
	blast_file = Path(path_data,'blastDB','blastall.out')
	
	tags = Path(in_path,'input','exp.tags')
	gattr = Path(in_path,'input','exp.gene_att')
	peps = Path(in_path,'input','exp.pep')

	print('Starting Panoct :',datetime.datetime.now())
	#panoct_path = Path('pangenome_softwares','panoct','panoct_v3.23','bin','panoct.pl')

	""" cmd = ['perl', str(panoct_path), 
	'-t', str(blast_file), 
	'-f', str(tags), 
	'-g', str(gattr), 
	'-P', str(peps), 
	'-S', 'Y', '-L', '1', '-M', 'Y', '-H', 'Y', '-V', 'Y', '-N', 'Y', 
	'-F', '1.33', 
	'-G', 'y', 
	'-c', '0,25,50,75,100', 
	'-T']
	call(cmd)
	"""
	parameters = [ str(blast_file), str(tags), str(gattr), str(peps) ]
	#parameters = [ str(blast_file), str(tags), str(peps) ]	
	stat = call_program(parameters,'panoct')
	print('Panoct finished at :',datetime.datetime.now())

	result_folder = Path('pangenome_softwares','panoct','results')
	
	#moving result files to result folder
	out_path = Path(in_path,'output')
	out_path.mkdir(exist_ok=True)
	
	for i in Path().glob('*.txt'):
		move(str(i), str(out_path))
	
	move('centroids.fasta',str(out_path))
	
	return stat
	
def panoct2families(in_path, out_path):
	
	panoct_clus = Path(out_path,'panoct_families.clus')
	hits_panoct = Path(in_path,'output','hits.txt')

	with open(panoct_clus,'w') as out:
		for line in open(hits_panoct,'r'):
			buffer = list()
			aux = re.split('\t+',line.strip())
			del aux[0]

			for gene in aux:
				buffer.append(gene.split(' ')[0])
			out.write(' '.join(buffer)+'\n')