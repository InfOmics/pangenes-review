from Bio import SeqIO
import re
from pathlib import Path
from shutil import copy, move
from subprocess import call

#creates files .nuc (fna) , .pep (faa) and .function
def createPgapInput(in_path, out_path):
	#copy and rename _nn.fa into .nuc
	out_path = Path(out_path,'input')
	out_path.mkdir(parents=True,exist_ok=True)

	#copy fasta files and change their names
	# *_aa.fa -> .pep
	for f_aa in Path(in_path,'protein').glob('*_aa.fasta'):
		nickname = f_aa.stem.split('_')[0]+'.pep'
		copy(f_aa,Path(out_path, nickname))
	
	# *_dna.fa -> .nuc
	for f_nn in Path(in_path,'dna').glob('*_dna.fasta'):
		nickname = f_nn.stem.split('_')[0]+'.nuc'
		copy(f_nn,Path(out_path, nickname))

	# .function files
	nucs = Path(out_path).glob('*.nuc')
	strains = list()
	for nuc in nucs:
		strains.append(nuc.stem)
		functions = list()
		for seq in SeqIO.parse(nuc,'fasta'):
			aux = seq.description.split(',')
			line = aux[0]+'\t-\thypothetical protein'
			functions.append(line)
		with open(Path(out_path,nuc.stem+'.function'),'w') as out:
			out.write('\n'.join(functions))
	with open(Path(out_path,'strains.txt'),'w') as st:
		st.write('+'.join(sorted(strains)))	

from modules.resource_control import call_program
def callPgap(in_path):
	#threads = 6
	#create the output folder
	output = Path(in_path,'output')
	output.mkdir(parents=True,exist_ok=True)
	
	strains = open(Path(in_path,'input','strains.txt'),'r').readline()
	#print(repr(strains))

	""" cmd = ['perl',Path('pangenome_softwares','pgap','PGAP.pl'),
	'--strains',strains,
	'--input',Path(in_path,'input'),
	'--output', output,
	'--cluster',
	#'--pangenome',
	'--method', 'GF',
	'--thread', str(threads)] """
	parameters = [strains, Path(in_path,'input'), output]
	stat = call_program(parameters,'pgap')
	
	#the orthologues files (All.cluster) is saved in the current directory instead of the 
	# output folder
	# idem for other files
	try:
		for f in Path().glob('All.*'):
			move(str(f),output)
		move('error.log',output)
		move('formatdb.log',output)
		move('genelist',output)

		#delete useless .pep files
		for d in Path().glob('*.pep'):
			d.unlink()
	except:
		pass
	return stat 
	
def pgap2families(in_path, out_path):
	
	#Note: i think the genes listes with double ,, means that they are paralogues
	cluster_file = Path(in_path,'output','1.Orthologs_Cluster.txt')
	i = 0
	if cluster_file.exists():
		out = open(Path(out_path,'pgap_families.clus'),'w')
		for cl in open(cluster_file,'r'):
			i+=1
			
			if cl.split('\t')[0] == 'ClutserID':
				continue
			
			cl = re.sub(r'^[0-9]+\t','',cl)
			
			aa = re.sub(r'(\t|,|-)+',' ',cl)
			out.write(aa)
			""" if i == 3:
				quit() """
		out.flush()
		out.close()
	else:
		print('MESSAGE:\n File "'+str(cluster_file)+'" does not exist, families will not be computed (file.clus)')
	