#create input for panseq
import re
from pathlib import Path
from subprocess import call
import networkx as nx

#fixes di headers of the fasta files
def createPanseqInput(in_path, out_path):
	files = Path(in_path,'dna').glob('*_dna.fasta')
	
	Path(out_path,'input').mkdir(exist_ok=True)
	
	for in_file in files:
		f_out = open(Path(out_path,'input',in_file.stem+'.fasta'),'w')
		i=0
		for line in open(in_file,'r'):
			i+=1

			#remove spaces and :, symbols with _ (to have only 1 character for spacing)
			a=re.sub(r'(:|,| )','_',line)

			#insert a label that identifies which genes belogs to a genome (required, otherwise all sequences are considered genomes)
			a=re.sub(r'>','>lcl|'+in_file.stem+'|',a)
			f_out.write(a)

		f_out.flush()
		f_out.close()

def createPanseqSettings(in_path):
	ncore = 6
	path_software_panseq = Path('pangenome_softwares','Panseq')
	with open(Path(in_path,'settings.txt'),'w') as f:

		settings_file = ['queryDirectory\t'+str(Path(in_path,'input')),
		'baseDirectory\t'+str(Path(in_path,'output')),
		'numberOfCores\t'+str(ncore),
		'mummerDirectory\t'+str(Path(path_software_panseq,'MUMmer3.23')),
		'blastDirectory\t'+str(Path(path_software_panseq,'blast2.9')),
		'muscleExecutable\t'+str(Path(path_software_panseq,'muscle3.8.31','muscle3.8.31')),
		'cdhitDirectory\t'+str(Path(path_software_panseq,'cd-hit-v4.8.1')),
		'fragmentationSize\t500',
		'percentIdentityCutoff\t85',
		'runMode\tpan',
		'overwrite\t1']
		
		f.write('\n'.join(settings_file))

from modules.resource_control import call_program
def callPanseq(settings_file):
	#cmd = ['perl','pangenome_softwares/Panseq/lib/panseq.pl', settings_file]
	parameters = [ settings_file ]
	stat = call_program(parameters, 'panseq')
	return stat
	
def panseq2families(in_path, out_path):
	panseq_file = Path(in_path,'pan_genome.txt')
	if panseq_file.exists():
		G = nx.Graph()
		for line in open(panseq_file,'r').readlines()[1:]:
			try: #because panseq has a bug where it doesn't print all the record till the end, leaving half a line of results.
				aux_gene_a = re.split(r'_{2}',line.rstrip().split('\t')[1].split('|')[2])[0].split('_')

				gene_a = aux_gene_a[0]+':'+aux_gene_a[1]+'_'+aux_gene_a[2]
				
				if line.rstrip().split('\t')[6] != 'NA':
					aux_gene_b = re.split(r'_{2}',line.rstrip().split('\t')[6].split('|')[2])[0].split('_')
					gene_b = aux_gene_b[0]+':'+aux_gene_b[1]+'_'+aux_gene_b[2]
					if gene_a != gene_b:
						G.add_node(gene_a)
						G.add_node(gene_b)
						G.add_edge(gene_a,gene_b)
			except:
				pass

		with open(Path(out_path,'panseq_families.clus'),'w') as handle:
			for f in list(nx.connected_components(G)):
				handle.write(' '.join(f)+'\n')
	else:
		print('MESSAGE:\n File "'+str(panseq_file)+'" does not exist, families will not be computed (file.clus)')
