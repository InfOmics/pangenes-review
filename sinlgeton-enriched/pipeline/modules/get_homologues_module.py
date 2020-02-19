from shutil import copy, move
from pathlib import Path
from subprocess import call
import re

def createGetHomologuesInput(in_path, out_path):
	
	files_path = Path(in_path,'protein')
	out_path = Path(out_path,'input')
	out_path.mkdir(parents=True, exist_ok=True)
	
	files = files_path.glob('*_aa.fasta')
	for f in files :
		print('copying file', f.stem, 'to', out_path)
		copy(f,out_path)

from modules.resource_control import call_program	
def callGetHomologues(in_path):
	path = Path('pangenome_softwares','get_homologues','get_homologues.pl')
	# -n number of cores ( default = 2 )
	#cmd = ['perl', path, '-d', Path(in_path,'input'), '-M', '-n', '6']
	parameters = [path, str(Path(in_path,'input'))]
	
	stat = call_program(parameters, 'get_homologues')
	try: #nel caso non venissero generati files
		move('input_homologues',in_path)
	except:
		pass

	return stat

def get_homologues2families(in_path, families_folder):
	path = Path(in_path,'input_homologues','tmp')
	ids = Path(path,'all.p2o.csv')
	clusters = Path(path,'all_ortho.mcl')

	id_dict = dict()
	if ids.exists() and clusters.exists():
		for line in open(ids):
			aux = line.split(',')
			id = aux[0]
			gene = aux[2]
			id_dict[id]=gene.rstrip()
		
		family_out = open(Path(families_folder,'get_homologues.clus'),'w')
		
		data = False
		families = list()
		genes = list()
	

		for line in open(clusters,'r'):
			
			if re.match('^begin',line):
				data = True
				continue

			if data == True:
				if re.match('^[0-9]',line):
					
					aux= re.split('\s{2,}',line.strip()) #splitta dove trova almeno 2 spazi consecutivi
					genes += aux[1].split(' ')
					
				else:
					genes += line.strip().split(' ')
				
				if genes[-1] == '$' or genes[-1] == ')':
					
					families.append(genes[:-1])
					genes = list()	

		del families[-2:]

		for family in families:
			
			if '0' in family: #the end of file is signed with a 0, this 0 has no meaning for homology
				#print(len(family), family)
				continue
			family_out.write( ' '.join([ id_dict[g_id] for g_id in family ])+'\n' )
							

		family_out.flush()
		family_out.close()
	else:
		print('MESSAGE:\n File "'+str(ids)+'" does not exist, families will not be computed (file.clus)')
