from pathlib import Path
from shutil import copy, move
from subprocess import call
import re
from Bio import SeqIO
from Bio.Seq import Seq

def createPangetInput(in_path, out_path):
	
	files = Path(in_path,'protein').glob('*_aa.fasta')
	Path(out_path,'input').mkdir(parents=True, exist_ok=True)
	
	for f in files:
		buffer = list()
		for seq in SeqIO.parse(f,'fasta'):
			name = seq.name.replace(':','_')
			seq.name = name
			seq.id = name
			seq.description = name
			buffer.append(seq)
		
		SeqIO.write(buffer,Path(out_path,'input',f.stem+'.faa'),'fasta')
	
	""" for f in files:
		copy(f, Path(out_path,'input',f.stem+'.faa')) """

	
from modules.resource_control import call_program
def callPanget(in_path):
	#path_parameters = Path(in_path,'parameters.txt')
	#cmd = ['perl', str(Path('pangenome_softwares','panget','PanGet_blastp_modified.pl')), str(path_parameters)]
	parameters = [str(in_path)]
	stat = call_program(parameters,'panget')

	#move output directory on the same level of the input folder
	#nota: da sistemare Destination path 'softwares_data/species10/indels02/panget/output' already exists
	try:
		move(str(Path(in_path,'input','output')),str(in_path))
	except:
		pass
	return stat
	
def panget2families(in_path,out_path):
	panget_file = Path(in_path,'output','process_input','conseverd')
	if panget_file.exists():
		with open(Path(out_path,'panget_families.clus'),'w') as out:
			for line in open(panget_file,'r'):
				aux = list()
				for g in re.split( r'\*+',line.rstrip()):
					if g != '':
						aux_gene = g.split('_')
						gene = aux_gene[0]+':'+aux_gene[1]+'_'+aux_gene[2]
						aux.append(gene)
				
				out.write(' '.join(aux) + '\n' )
	else:
		print('MESSAGE:\n File "'+str(panget_file)+'" does not exist, families will not be computed (file.clus)')
	
