from pathlib import Path
from subprocess import call
import re

def createMicropanInput(in_path,path_micropan_data):
	
	path_blastall = Path(in_path,'blastDB','blastall.out')
	
	input_micropan = Path(path_micropan_data,'input')
	input_micropan.mkdir(exist_ok=True, parents=True)
	
	files_dict = dict()
	for line in open(path_blastall,'r'):
		if re.match('^#', line):
			continue
		
		aux = line.split('\t')
		g1 = aux[0].split('_')[1]
		g2 = aux[1].split('_')[1]
		#print(g1,g2)

		if g1 not in files_dict:
			files_dict[g1] = dict()
		
		if g2 not in files_dict[g1]:
			files_dict[g1][g2]=list()

		files_dict[g1][g2].append(line)

	for gene1 in files_dict:
		for gene2 in files_dict[gene1]:

			filename = (gene1+'_vs_'+gene2+'.txt').replace('SE','GID')
			#print(filename)
			
			with open(Path(input_micropan,filename),'w') as out:
				out.write(''.join(files_dict[gene1][gene2]))
	
from modules.resource_control import call_program
def callMicropan(path_micropan_data):
	
	in_path = Path(path_micropan_data,'input')
	out_path = Path(path_micropan_data,'output')
	out_path.mkdir(exist_ok=True)

	#micropan_script = Path('pangenome_softwares','micropan','run_micropan.R')
	#cmd = ['Rscript', micropan_script, in_path, out_path]
	parameters = [in_path, out_path]
	stat = call_program(parameters, 'micropan')
	return stat
	
def micropan2families(in_path, out_path):

	with open(Path(out_path,'micropan_families.clus'),'w') as out:
		clus = dict()
		with open(Path(in_path,'output','micropan.out'),'r') as handle:
			next(handle) #skip the first line
			for line in handle:
				aux = line.rstrip().split('\t')
		
				if int(aux[1]) not in clus:
					clus[int(aux[1])] = list()
				clus[int(aux[1])].append(aux[0].strip('"'))
			for i in clus:
				out.write(' '.join(clus[i])+'\n')



