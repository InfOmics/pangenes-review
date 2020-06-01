from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature , FeatureLocation
from Bio.Seq import Seq
from collections import OrderedDict
from subprocess import call


def createPanxInput(in_path, out_path):
	
	for f in (i for i in Path(in_path,'dna').glob('*_dna.fasta')):
		genome_name = f.stem.split('_')[0]
		
		features = list()
		genome_seq = list()
		position=0 #il FeatureLocation viene aumentato di 1 nello start ma non nell'end perchè gbk ha base 1 ma FL ha base 0
		
		for gene in SeqIO.parse(f,'fasta', IUPAC.ambiguous_dna):
			gene_len= len(gene.seq)	
			locus_tag = gene.id.rstrip(',')

			sf_gene = SeqFeature( id=gene.id, type='gene', location=FeatureLocation(position, position+gene_len-1, strand=1) )
			sf_gene.qualifiers = {'locus_tag':[locus_tag]}

			sf_cds = SeqFeature( id=gene.id, type='CDS', location=FeatureLocation(position, position+gene_len-1, strand=1) )
			sf_cds.qualifiers = {'locus_tag':[locus_tag]}
			sf_cds.qualifiers['product'] = 'hypothetical protein'
			sf_cds.qualifiers['translation'] = gene.seq.translate(table=11)
			
			features.append(sf_gene)
			features.append(sf_cds)

			genome_seq.append(str(gene.seq))
			position += gene_len

		annotations = {'source': 'ALF (A Simulation Framework for Genome Evolution)', 'organism': 'Synthetic genome'}
		""" source = SeqFeature(type='source', location=FeatureLocation(0,len(genome_seq)-1))
		source.qualifier = {'organism': ['Synthetic mycoplasma genitalium'] , 'mol_type': ['genomic DNA'], 'strain',  } """
		s = SeqRecord(id=genome_name, seq=Seq(''.join(genome_seq),IUPAC.ambiguous_dna), features=features, annotations=annotations)
		SeqIO.write(s,Path(out_path,genome_name+'.gbk'),'genbank')

from modules.resource_control import call_program
def callPanX(in_path):
	panx = Path('pangenome_softwares','panX','panX.py')
	diamond = Path('pangenome_softwares','panX','diamond-linux64','diamond')
	# -dmp : directory diamond ( MSA )
	# -fn : folder name
	# -sl : species name for temporary folders
	# -t : number of threads
	#cmd = ['python', panx, '-dmp', diamond,'-fn', in_path, '-sl','synthetic','-t', '6','-ct']
	parameters = [in_path]
	print(in_path)
	stat = call_program(parameters,'panx')
	return stat


import os
def panX2families(in_path, out_path):
	print("@vb@panx families")	
	i = 0
	out = Path(out_path,'panx_families.clus')
	with open(out,'w') as off:
		for file in os.listdir(str(in_path) + os.sep + 'geneCluster'):
			#print("@vb@"+str(file))
			if  file.startswith('GC') and file.endswith('.fna'):
				print('@vb@'+file)
				for line in open(str(in_path) + os.sep + 'geneCluster'+os.sep+file):
					if len(line) >0 and line[0]=='>':
						cc = line.strip().split('|')[1].split('-')[0]
						off.write(cc+" ")
				off.write('\n')
	#panx_file = Path(in_path,'allclusters_final.tsv')
	#if panx_file.exists():
	#	with open(out,'w') as out:
	#		for line in open(panx_file,'r'):
	#			i+=1
	#			buffer = list()
	#			for gene in line.rstrip().split('\t'):
	#				buffer.append(gene.split('|')[1])
	#			
	#			out.write(' '.join(buffer)+'\n')
	#else:
	#	print('MESSAGE:\n File "'+str(panx_file)+'" does not exist, families will not be computed (file.clus)')

