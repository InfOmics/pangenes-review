from pathlib import Path
from subprocess import call
import re

def createMicropanInput(in_path,path_micropan_data):
	
#	path_blastall = Path(in_path,'blastDB','blastall.out')
#	input_micropan = Path(path_micropan_data,'input')
#	input_micropan.mkdir(exist_ok=True, parents=True)
#	files_dict = dict()
#	for line in open(path_blastall,'r'):
#		if re.match('^#', line):
#			continue
#		aux = line.split('\t')
#		g1 = aux[0].split('_')[1]
#		g2 = aux[1].split('_')[1]
#		#print(g1,g2)
#		if g1 not in files_dict:
#			files_dict[g1] = dict()
#		if g2 not in files_dict[g1]:
#			files_dict[g1][g2]=list()
#		files_dict[g1][g2].append(line)
#	for gene1 in files_dict:
#		for gene2 in files_dict[gene1]:
#			#filename = (gene1+'_vs_'+gene2+'.txt').replace('SE','GID')
#			filename = ("GID"+gene1+'_vs_GID'+gene2+'.txt')
#			#print(filename)
#			with open(Path(input_micropan,filename),'w') as out:
#				out.write(''.join(files_dict[gene1][gene2]))

	path_blastall = Path(in_path,'blastDB','blastall.out')

	input_micropan = Path(path_micropan_data,'input')
	input_micropan.mkdir(exist_ok=True, parents=True)

	gidmap = dict()

	files_dict = dict()
	for line in open(path_blastall,'r'):
		if re.match('^#', line):
		        continue

		aux = line.split('\t')
		#print(aux)
		g1 = aux[0].split('_')[1]
		g2 = aux[1].split('_')[1]
		#print(g1,g2)

		if g1 not in files_dict:
		        files_dict[g1] = dict()

		if g2 not in files_dict[g1]:
		        files_dict[g1][g2]=list()

		files_dict[g1][g2].append(line)

	for gene1 in files_dict:
		if gene1 not in gidmap:
		        gidmap[gene1] = str(len(gidmap)+1)
		for gene2 in files_dict[gene1]:
		        if gene2 not in gidmap:
		                gidmap[gene2] = str(len(gidmap)+1)

		        gid1 = gidmap[gene1]
		        gid2 = gidmap[gene2]
		        #filename = (gene1+'_vs_'+gene2+'.txt').replace('SE','GID')
		        filename = ("GID"+gid1+'_vs_GID'+gid2+'.txt')
		        #print(filename)

		        with open(Path(input_micropan,filename),'w') as out:
		                olines = files_dict[gene1][gene2]
		                for i in range(len(olines)):
		                        #print(olines[i])
		                        cc = olines[i].split('\t')
		                        gid ='GID'+ gidmap[cc[0].split('_')[1]]
		                        geneid ='seq'+ cc[0].split(':')[0]
		                        cc[0] = gid+"_"+geneid
		                        gid ='GID'+ gidmap[cc[1].split('_')[1]]
		                        geneid ='seq'+ cc[1].split(':')[0]
		                        cc[1] = gid+"_"+geneid
		                        cc = cc[0:2] + cc[10:]
		                        olines[i] = '\t'.join(cc)+"\n"
		                        #print(olines[i])
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



