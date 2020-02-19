from ete3 import Tree, TreeStyle, TreeNode, TextFace, NodeStyle
import random
random.seed(4)
from pathlib import Path
from shutil import copy, move
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from subprocess import call
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import json
import tarfile
import networkx as nx
import argparse

#UTILITY FUNCTIONS
#used to sort tuples based on the third field "leaves distance"
def getKey(item):
	return item[2]

#method to recursively remove a folder
def delete_folder(pth) :
    for sub in pth.iterdir() :
        if sub.is_dir() :
            delete_folder(sub)
        else :
            sub.unlink()
    pth.rmdir()


#METHODS
#backtranslate protein genes of the genomes into DNA
def backtranslate(path):
	aa_codon_table = { #from translation table 11 NCBI
		'F':'TTT' ,
		'Y':'TAT' ,
		'C':'TGT' ,
		'W':'TGG' ,
		'H':'CAT' ,
		'L':'CTC' ,
		'P':'CCG' ,
		'Q':'CAG' ,
		'I':'ATC' ,
		'T':'ACC' ,
		'N':'AAC' ,
		'S':'AGC' , 
		'M':'ATG' , 
		'K':'AAG' ,
		'R':'AGG' ,
		'D':'GAC' ,
		'V':'GTG' , 
		'A':'GCG' ,
		'E':'GAG' , 
		'G':'GGG' ,
		'*':'TAG' 
	}

	
	#path = Path(path,'genome','genome','DB')
	Path(path,'dna').mkdir(parents=True, exist_ok=True)
	Path(path,'protein').mkdir(parents=True, exist_ok=True)

	for filepath in sorted(path.glob('*_aa.fa')):
			
		to_fasta = list()
		for sequence in SeqIO.parse(filepath,'fasta'):
					
			aa_seq = list()
			for char in sequence.seq:
				aa_seq.append(aa_codon_table[char])
			
			sequence.seq = Seq(''.join(aa_seq), IUPAC.ambiguous_dna )
			to_fasta.append(sequence)
	
		outfile = filepath.stem.replace('_aa','_dna')+'.fasta'
		SeqIO.write(to_fasta,Path(path,'dna',outfile),'fasta')

		new_path = Path(path,'protein',filepath.name.replace('.fa','.fasta'))
		#move(filepath, new_path) #per ora copio per poter rilanciare il programma
		copy(filepath, new_path)
		#print(outfile)

#CVTree (composition vector tree)
from subprocess import PIPE
def runCVTree(in_path, out_path):
	#cmd = ['./parallel.out', str(in_path)+'/']
	cmd = ['./cvtree.sh', str(in_path)+'/']
	if Path('results-parallel.csv').exists():
		Path('results-parallel.csv').unlink()
	print("calling cvtree repeatedly...")
	while not Path('results-parallel.csv').exists():
		call(cmd)
		
	#call(cmd)
	if Path(out_path,'results-parallel.csv').exists():
		Path(out_path,'results-parallel.csv').unlink()
	move('results-parallel.csv', out_path)
		

#CVTree to matrix
def csv2matrix(in_path):

	data = pd.read_csv(in_path)
	m_size = max(max(data['indexA']),max(data['indexB']))+1
	matrix = np.zeros(shape=(m_size,m_size))

	labels = list()
	
	for index, row in data.iterrows():
		
		idxA = int(row['indexA'])
		idxB = int(row['indexB'])
		correlation = row['correlation']
		
		#distance = 1-correlation
		matrix[idxA,idxB] = 1.0 - correlation
		matrix[idxB,idxA] = 1.0 - correlation
		
		bactA = Path(row['bacteriaA']).stem
		if bactA not in labels:
			labels.append(bactA)
		
	labels.append(Path(data['bacteriaB'].iloc[-1]).stem)

	df = pd.DataFrame(matrix, index=labels, columns=labels)
	
	#i'm removing the '_dna' from indexes and columns in because they are kept
	df.index = [i.replace('_dna','') for i in df.index]
	df.columns = [c.replace('_dna','') for c in df.columns]
	df.to_csv(Path(in_path.parent,'CVTree_matrix.csv'))
	
	return df



#input path of the tar file, output directory where you save the family clus
def alf2families(in_path, out_dir): 

	#extract VP.tgz ( it contains data about homologues, orthologues, paralogues and xenologues )
	tar = tarfile.open(Path(in_path,'VP.tgz'), "r:gz")
	tar.extractall(path=in_path)
	tar.close()

	print('\nGenerating graph...')
	G = nx.Graph() #inizializzo il grafo
	orthologues_filepath = Path(in_path,'VP','OP.drw')
	
	#ortologhi
	orth = open(orthologues_filepath) #non prendo il file HP con tutti gli omologhi perchè è composto ortologhi, paraloghi e xenologhi(LGT). Usare file separati permette di fare statistiche su 3 tipi di omologhi in modo piu chiaro senza impattare sui tempi di analisi del codice
	a = orth.readlines()
	b = a[1].lstrip('OP := ').rstrip(':\n')

	array = json.loads(b) #dati del file di omologhi ( ortologhi, paraloghi, non dovrebbero essere inclusi gli xenologhi perchè è genoma i vs genoma j e le LGT sono tra 1 genoma e molti genomi paralleli )

	for i in range(len(array)): #indice genoma di riferimento --> "genoma i" vs
		
		id_i = '{0}'.format(str(i+1).zfill(3)) # i+1 perchè i genomi partono da SE001 mentre gli indici partono da 0
		genome_i = 'SE'+id_i

		for j in range(len(array[i])):#indice genoma a cui punta il genoma i --> genoma i vs "genoma j"
			
			id_j = '{0}'.format(str(j+1).zfill(3))
			genome_j = 'SE'+id_j
			#print(genome_i,'vs',genome_j)

			for p in range(len(array[i][j])): #vado a vedere il contenuto delle posizioni
				
				if array[i][j][p]: #se non è una posizione vuota
					gene_i = 'G'+str(p+1)
					gene_j = 'G'+str(array[i][j][p][0])

					node_i = (genome_i,gene_i)
					node_j = (genome_j,gene_j)
					#print( node_i, node_j )
					
					if node_i not in G:
						G.add_node(node_i)
					if node_j not in G:
						G.add_node(node_j)
					
					G.add_edge(node_i,node_j)

	orth.flush()
	orth.close()


	#paraloghi ( da aggiungere al grafo ) [i paraloghi sono segnati sul file correttamente]
	paralogues_filepath = Path(in_path,'VP','PP.drw')
	para = open(paralogues_filepath) #non prendo il file HP con tutti gli omologhi perchè è composto ortologhi, paraloghi e xenologhi(LGT). Usare file separati permette di fare statistiche su 3 tipi di omologhi in modo piu chiaro senza impattare sui tempi di analisi del codice
	a = para.readlines()
	b = a[1].lstrip('PP := ').rstrip(':\n')

	array = json.loads(b)

	for d in range(len(array)):
		id_i = '{0}'.format(str(d+1).zfill(3)) # prendo solo le diagonali quindi i == j, uso d come indice
		genome_i = 'SE'+id_i

		id_j = '{0}'.format(str(d+1).zfill(3))
		genome_j = 'SE'+id_j

		for p in range(len(array[d][d])): #vado a vedere il contenuto delle posizioni
			gene_i = 'G'+str(p+1)
			node_i = (genome_i,gene_i)
			if node_i not in G:
				G.add_node(node_i)

			if array[d][d][p]: #se non è una posizione vuota
				gene_i = 'G'+str(p+1)
				gene_j = 'G'+str(array[d][d][p][0])

				node_i = (genome_i,gene_i)
				node_j = (genome_j,gene_j)
				#print( node_i, node_j )
				
				if node_i not in G:
					G.add_node(node_i)
				if node_j not in G:
					G.add_node(node_j)
				
				G.add_edge(node_i,node_j)
	para.flush()
	para.close()
	
	#print('Gene Families found',len(list(nx.connected_components(G))),'\n')
	
	#file.clus
	cluspath = Path(out_dir,'alf_families.clus')
	print('Saving families.clus in :',cluspath)
	with open(cluspath,'w') as handle:
		for f in list(nx.connected_components(G)):
			#print(type(f), sorted(f)) #f è un set
			aux = list()
			for s in sorted(f):
				aux.append(s[1]+'_'+s[0])

			#print(repr('\t'.join(aux)+'\n'))			
			handle.write(' '.join(aux)+'\n')


def getGenomeDistribution(path_clus):
	#qui calcolo l'istogramma della "genomes per class distribution"

	histogram = dict() #numero di famiglie che toccano X genomi
	for line in open(path_clus,'r'):
		genes = line.strip().split(' ')
		
		aux = set() #set perchè non ha ripetizioni, paraloghi appartengono allo stesso genoma
		for g in genes:
			genome = g.split('_')[1]
			aux.add(genome)

		histogram[int(len(aux))] = histogram.get(len(aux),0)+1

	return histogram
 

def get_datasetFASTA(tree, pth_mintree_dataset, pth_alf_dataset):
	fasta_path = Path(pth_alf_dataset,'DB','dna')
	
	dataset_path = Path(pth_mintree_dataset,'dna')
	dataset_path.mkdir(parents=True, exist_ok=True)

	leaves = [leaf.name for leaf in tree]
	for f in [leaf+'_dna.fasta' for leaf in leaves]:
		copy(Path(fasta_path,f),Path(dataset_path,f))
	
	fasta_path = Path(pth_alf_dataset,'DB','protein')
	dataset_path = Path(pth_mintree_dataset,'protein')
	dataset_path.mkdir(parents=True, exist_ok=True)

	for f in [leaf+'_aa.fasta' for leaf in leaves]:
		copy(Path(fasta_path,f),Path(dataset_path,f))
	
	return dataset_path.parent


#create the gene family file from the complete family .clus file of ALF
def get_familiesSelection(tree,dataset_name,out_path, pth_alf_dataset):
	with open(Path(out_path,dataset_name+'.clus'),'w') as out:
		leaves = [leaf.name for leaf in tree]
		for line in open(Path(pth_alf_dataset,'alf_families.clus'),'r'):
			buffer = list()
			family = line.strip().split(' ')
			#rimuovo i geni che non appartengono ai genomi delle foglie
			for gene in family:
				genome = gene.split('_')[1]
				if genome in leaves:
					buffer.append(gene)
				else:
					pass
			#print(buffer)
			#print(len(buffer))
			if buffer: #se il buffer non è vuoto aggiunge la famiglia sul file
				out.write(' '.join(buffer)+'\n')


#shows the selected genomes on the complete phylogenetic tree
def highlight_selection(tree, selection, out_name, seed=True):

	highlight_tree = tree.copy() 
	# Draws nodes as small red spheres of diameter equal to 10 pixels
	nstyle = NodeStyle()
	nstyle["shape"] = "sphere"
	nstyle["size"] = 10
	nstyle["fgcolor"] = "darkred"

	#evidenzio diversamente il genoma seed
	nstyle_seed = NodeStyle()
	nstyle_seed["shape"] = "square"
	nstyle_seed["size"] = 10
	nstyle_seed["fgcolor"] = "yellow"

	# Applies the same static style to all nodes in the tree. Note that,
	# if "nstyle" is modified, changes will affect to all nodes
	for genome in selection:
		node = highlight_tree.get_leaves_by_name(genome)[0]
		node.set_style(nstyle)

	#nodo centrale
	if seed:
		seed = highlight_tree.get_leaves_by_name(selection[0])[0]
		seed.set_style(nstyle_seed)
	
	highlight_tree.render(out_name,tree_style=ts)

def histogramDistribution(in_path, hist_name):
	histogram = getGenomeDistribution(in_path)
	nof_genomes = max(histogram)
	
	hist = [0]*nof_genomes
	for k in histogram:
		hist[k-1] = histogram[k] 
		
	x = range(1,len(hist)+1)

	plt.bar(x,hist,color='blue')
	#plt.xticks([],[]) #no x-axis ticks
	plt.xlabel('Genomes')
	plt.ylabel('Gene Families')
	plt.title('Genomes per family distribution')
	
	plt.savefig(Path(in_path.parent,hist_name))
	plt.clf()
	
def createALFinput(nof_species, seed, loss_rate=0, indel_rate=0.003, dupl_level=0):
	input_file = ["SetRand("+str(seed)+"): # use this with any number, if you want reproducable results",
	"webRequest := false;",
	"uuid := 's0-uuid';",
	"# name of simulation - you may want to change this",
	"mname := dataset;",
	"# directories for file storage - you may want to change these",
	"wdir := 'alf_results/datasets/'.mname.'/';",
	"dbdir := 'DB/';",
	"#dbAncdir := 'DBancestral/';",
	"# time scale for simulation (PAM is default)",
	"unitIsPam := false:",
	"# parameters concerning the root genome",
	"realorganism := './input/mycoplasma.g-1.db';",
	"# parameters concerning the species tree",
	"treeType := 'BDTree';",
	"birthRate := 0.2;",
	"deathRate := 0.15;",
	"mutRate := 5000;",
	"NSpecies := "+str(nof_species)+";",
	"ultrametric := false;",
	"#scaleTree := false;# parameters concerning the substitution models",
	"substModels := [SubstitutionModel('WAG')];",
	"#indelModels := [IndelModel(0.002,GEOM,[0.333],50)]; #qui ho una maxLen di indels messa a 50, quindi posso anche avere delle indels lunghe fino a 50",
	"indelModels := [IndelModel("+str(indel_rate)+", ZIPF, [1.821])];",
	"rateVarModels := [RateVarModel()];",
	"# parameters concerning gene duplication",
	"geneDuplRate := "+str(dupl_level)+";",
	"numberDupl := 10;",
	"fissionDupl := 0.0;",
	"fusionDupl := 0.0;",
	"# parameters concerning gene loss",
	"geneLossRate := "+str(loss_rate)+";",
	"numberLoss := 10;",
	"# parameters concerning LGT",
	"lgtRate := 0;",
	"#orthRep := 0;",
	"#lgtGRate := 0;",
	"#lgtGSize := 10;",
	"# parameters concerning rate heterogeneity among genes",
	"#amongGeneDistr := 'Gamma';",
	"#aGAlpha := 1;",
	"amongGeneDistr := 'None';",
	"# select the output you want (apart from the species tree and genomes)",
	"simOutput := { 'GeneTrees' , 'VP', 'Fasta', NULL }:"]

	with open('input_alf.drw','w') as f:
		f.write('\n'.join(input_file))

#MAIN


parser = argparse.ArgumentParser(description='Input parameters')
parser.add_argument('-n', nargs=1 , type=int, dest="subset_dimension",required=True, 
help="dimension of the minimum and random datasets")
parser.add_argument('-s', nargs=1 , type=int, dest="nof_species",required=True, 
help="number of species")
parser.add_argument('-r', nargs=1 , type=int, dest="repetitions",required=False, default=[1], 
help="generate multiple dataset")
parser.add_argument('-loss', nargs=1 , type=float, dest="loss_rate",required=False, default=[0],
help="loss rate for ALF")
parser.add_argument('-indel', nargs=1 , type=float, dest="indel_rate",required=False, default=[0.003],
help="indel rate for ALF")
parser.add_argument('-dupl', nargs=1 , type=float, dest="dupl_level",required=False, default=[0.0],
help="duplication level for ALF")
parser.add_argument('--analysis_only', action='store_false', dest='analysis_only', 
help="use this flag to run only the analysis of the dataset obtained with alf")

args = parser.parse_args()

genome_group_dimension = args.subset_dimension[0]
nof_species = args.nof_species[0]
repetitions = args.repetitions[0]
loss_rate = args.loss_rate[0]
indel_rate = args.indel_rate[0]
dupl_level = args.dupl_level[0]

print('nof_species:', nof_species)
print('repetitions dataset:', repetitions)
print('subset dimension:',genome_group_dimension)
print('loss rate:', loss_rate)
print('indel rate:', indel_rate)
print('duplication level', dupl_level)
#print('analysis only:', 'no' if args.analysis_only else 'yes', end='\n\n')

input_datasets = Path('input_datasets')

#tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c" #circular tree
ts.show_branch_length = True
#ts.scale = 20 #per scalare(ridurre o aumentare) la dimensione dell'albero
ts.arc_start = 0 # 0 degrees = 3 o'clock
ts.arc_span = 360



#ALF dataset#

#convert mycoplasma genome into darwin format for ALF

call(['alf/bin/fasta2darwin','input/mycoplasma.g-1.fasta'])
if Path('alf_results').exists():
		delete_folder(Path('alf_results'))

if Path('input_datasets').exists():
	delete_folder(Path('input_datasets'))


startseed = 123
startseed = abs(int(repetitions * startseed) + startseed)
endseed = startseed + repetitions


for seed in range(startseed, endseed):
	createALFinput(nof_species, seed, loss_rate=loss_rate, indel_rate=indel_rate, dupl_level=dupl_level)
	##run alf
	call(['./alf/bin/alfsim','input_alf.drw'])

for dataset in Path('alf_results','datasets','dataset').glob('*'):

	#backtranlslate protein files into dna because alf has a bug
	backtranslate(Path(dataset,'DB'))

	#retrieve genomes families (alf_families.clus)
	print('Writing gene families file : alf_families.clus')
	alf2families(dataset,dataset)
	
	"""
	#CVTree on dna fasta files
	print('Running CVTree on :',dataset)
	runCVTree(Path(dataset,'DB','dna'),dataset)

	df = csv2matrix(Path(dataset,'results-parallel.csv'))

	#heatmap with cluster
	sns.clustermap(df, cmap='RdBu')
	plt.savefig(Path(dataset,'heatmap_alfFull.png'))
	plt.clf() #clear figure
	"""

	#computing genomes per family distribution
	histogramDistribution(Path(dataset,'alf_families.clus'),'histogram_alfFull.png')

	#RETRIEVE PHYLOGENETIC TREES FOR ALF DATASETS
	nwk = Path(dataset,'RealTree.nwk')
	t = Tree(str(nwk))
	t.render(str(dataset)+'/TreeFull.png',tree_style=ts)
	#t.show(tree_style=ts)
	

	###CREATE CLOSE GENOMES SUBSETS : Tree -> selection -> retrieval from alf data###
	#dataset folder mintree (the number corresponds to the relative alf dataset)
	if len(dataset.stem.split('_'))>1:
		dataset_mintree = Path(input_datasets,'dataset_mintree','mintree_'+dataset.stem.split('_')[1])
	else:
		dataset_mintree = Path(input_datasets,'dataset_mintree','mintree')
	
	dataset_mintree.mkdir(parents=True, exist_ok=True)	

	mindist_list = ([],float('inf'))
	leaves = [leaf.name for leaf in t]
	
	for i in leaves:
		
		buffer = list()
		for j in leaves:		
			if i != j:
				buffer.append((i,j, t.get_distance(i,j)))

		# meno 1 perchè calcola i 49 genomi piu vicini +1 che è il genoma iniziale a fare un gruppo da 50		
		selection = sorted(buffer,key=getKey)[0:genome_group_dimension-1]

		distance_sum = 0
		for q in selection:
			distance_sum += q[2]
		
		if distance_sum < mindist_list[1]:
			mindist_list = (selection,distance_sum)


	#minimum distance list
	min_list = [mindist_list[0][0][0]] #il genoma confrontato con tutti gli altri (seed da acui è calcolata la distanza dagli altri nodi)
	#print(min_list)
	for g in mindist_list[0]: #per ogni tupla
		min_list.append(g[1])
	
	#pruning tree with min
	min_tree = t.copy()
	min_tree.prune(min_list)
	min_tree.render(str(dataset_mintree)+'/mintree.png',tree_style=ts)

	#Retrieve fasta files for the genomes selected for the mintree dataset
	#Parameter passed are the phylogenetic tree, the path to the mintree dataset, the relative path to alf dataset
	get_datasetFASTA(min_tree, dataset_mintree, dataset)
	get_familiesSelection(min_tree, dataset_mintree.stem, dataset_mintree, dataset)
	highlight_selection(t,min_list,str(dataset_mintree)+'/highlight_mintree.png')
		
	#histogram genomes per family distribution
	print('computing histogram distribution of close distance genomes...')
	histogramDistribution(Path(dataset_mintree,dataset_mintree.stem+'.clus'),'histogram_mintree.png')
	print('...done!')
	
	####
	"""
	#CVTree on randomtree dataset
	print('Running CVTree on random genomes...')
	runCVTree(Path(dataset_mintree,'dna'), dataset_mintree)
	df = csv2matrix(Path(dataset_mintree,'results-parallel.csv'))
	print('...done!')

	#heatmap with cluster
	print('computing heatmap of random genomes...')
	sns.clustermap(df, cmap='RdBu')
	plt.savefig(Path(dataset_mintree,'heatmap_mintree.png'))
	plt.clf()
	print('...done!') 
	"""
	

	##CREATE RANDOM SELECTION GENOMES SUBSETS##
	#dataset folder randomtree (the number corresponds to the relative alf dataset)
	if len(dataset.stem.split('_'))>1:
		dataset_randomtree = Path(input_datasets,'dataset_randomtree','randomtree_'+dataset.stem.split('_')[1])
	else:
		dataset_randomtree = Path(input_datasets,'dataset_randomtree','randomtree')
	
	dataset_randomtree.mkdir(parents=True, exist_ok=True)
	
	#leaves = [leaf.name for leaf in t]


	random_list = random.sample(leaves, k=genome_group_dimension)

	random_tree = t.copy()

	random_tree.prune(random_list)
	#for l in random_tree:
	#	if l.name not in random_list:
	#		l.delete()

	
	random_tree.render(str(dataset_randomtree)+'/randomtree.png',tree_style=ts)

	get_datasetFASTA(random_tree,dataset_randomtree, dataset)
	get_familiesSelection(random_tree, dataset_randomtree.stem, dataset_randomtree,dataset)
	highlight_selection(t,random_list,str(dataset_randomtree)+'/highlight_randomtree.png',seed=False)
	
	"""
	#CVTree on randomtree dataset
	print('Running CVTree on random genomes...')
	runCVTree(Path(dataset_randomtree,'dna'), dataset_randomtree)
	df = csv2matrix(Path(dataset_randomtree,'results-parallel.csv'))
	print('...done!')

	#heatmap with cluster
	print('computing heatmap of random genomes...')
	sns.clustermap(df, cmap='RdBu')
	plt.savefig(Path(dataset_randomtree,'heatmap_randomtree.png'))
	plt.clf()
	print('...done!')
	"""

	#histogram genomes per family distribution
	print('computing histogram distribution of random genomes...')
	histogramDistribution(Path(dataset_randomtree,dataset_randomtree.stem+'.clus'),'histogram_randomtree.png')
	print('...done!')
	
