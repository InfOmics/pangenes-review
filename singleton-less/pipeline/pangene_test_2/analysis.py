from pathlib import Path
from subprocess import call
import numpy as np
import pandas as pd


def getGenomeDistribution(path_clus):
	#qui calcolo l'istogramma della "genomes per class distribution"

	histogram = dict() #numero di famiglie che toccano X genomi
	for line in open(path_clus,'r'):
		genes = line.strip().split(' ')
		
		aux = set() #set perchè non ha ripetizioni, paraloghi appartengono allo stesso genoma
		for g in genes:
			genome = g.split('_')[1]
			aux.add(genome)

		histogram[len(aux)] = histogram.get(len(aux),0)+1

	return histogram


def getClasses(path_alf):
	family_number = -1
	nof_seqs = 0
	classes = dict()

	for line in open(path_alf,'r'):
		family_number +=1
		alf_cluster = [gene for gene in line.strip().split(' ')]
		nof_seqs += len(alf_cluster)
		classes['family_'+str(family_number)] = set(alf_cluster)	

	return classes, nof_seqs

#performance analysis with True/False Positives/Negatives
def performance_analysis(path_clus, dbclasses):

	name = path_clus.stem
	
	clusters = list()
	for line in open(path_clus,'r'):
		cc = line.strip().split(' ')
		clusters.append( cc )

	#print('number of clusters:', len(clusters))

	e_classes = set()
	for k,vv in dbclasses.items():
		for v in [ (i,j) for i in vv for j in vv if i != j] :
			e_classes.add(v)

	e_clusters = set()
	for vv in clusters:
		for v in [ (i,j) for i in vv for j in vv if i != j] :
			e_clusters.add(v)

	
	tp = len( e_classes & e_clusters ) #True positives
	fp = len( e_clusters - e_classes) #False positives
	fn = len( e_classes - e_clusters) #False negatives
	tn = ((nof_seqs * nof_seqs) - nof_seqs) - len(e_classes) #True negatives
	precision = tp/(tp + fp)
	recall = tp/(tp + fn)
	result = {'classes_links':len(e_classes),
	'cluster_links':len(e_clusters),
	'TP':tp,
	'FP':fp,
	'FN':fn,
	'TN':tn,
	'precision':tp/(tp + fp),
	'recall':tp/(tp + fn),
	'true_negative_rate':tn/(tn + fp),
	'c-diff':None,
	'f1-score':2*(precision*recall)/(precision+recall) if (precision+recall) != 0 else 0}
	
	#labels = ['classes_links','cluster_links','TP','FP','FN','TN','precision','recall','true_negative_rate','c-diff']
	#result = [len(e_classes),len(e_clusters),tp,fp,fn,tn,tp/(tp + fp),tp/(tp + fn),tn/(tn + fp)]
	return result
	



##MAIN NEW

#mintree dataset analysis
n_genes_mintree = dict()
mintree_hist_list = dict()
mintree_hist_diff = dict()

analysis_result = Path('analysis')
analysis_result.mkdir(exist_ok=True)

#dataset mintree analysis
mintree_result = Path(analysis_result,'analysis_mintree')
mintree_result.mkdir(exist_ok=True)

clus_mintree = Path('input_datasets', 'dataset_mintree','mintree','mintree.clus')

#analyze dataset_mintree.clus, used as reference to compare other pangenome softwares

#classes are the families of genes, nof_seqs is the number of sequences/genes
classes, nof_seqs = getClasses(clus_mintree)

#histogram in the form of dictionary (only values that occour at least once)
mintree_gen_distr = getGenomeDistribution(clus_mintree)

#number of genomes
genomes_list = set()
for key, value in classes.items():
	for gene in value:
		genomes_list.add(gene.split('_')[1])

nof_genomes = len(genomes_list)

#create the histogram
hist_mintree = [0]*nof_genomes
for k in mintree_gen_distr:
	hist_mintree[k-1] = mintree_gen_distr[k]

#saving data
n_genes_mintree[clus_mintree.stem] = nof_seqs
mintree_hist_list[clus_mintree.stem] = hist_mintree

#calculate histograms for all the clus files obtained with pangenome softwares
for clus_file in Path('gene_families','mintree','mintree').glob('*'):
	print(clus_file)

	f_classes, f_nof_seqs = getClasses(clus_file)

	file_gen_distr = getGenomeDistribution(clus_file)
	
	#number of genomes
	f_genomes_list = set()
	for key, value in classes.items():
		for gene in value:
			f_genomes_list.add(gene.split('_')[1])

	f_nof_genomes = len(f_genomes_list)
	print(f_nof_genomes)
	print(file_gen_distr)
		
	#create the histogram
	f_hist_file = [0]*f_nof_genomes
	for k in file_gen_distr:
		f_hist_file[k-1] = file_gen_distr[k]
	#print(hist_file)

	#saving
	n_genes_mintree[clus_file.stem] = f_nof_seqs
	mintree_hist_list[clus_file.stem] = f_hist_file

	hist_aux = list()
	for i in range(len(hist_mintree)):
		hist_aux.append(abs(hist_mintree[i] - f_hist_file[i]))
	""" print(hist_mintree)
	print(f_hist_file)
	print(hist_aux) """
	mintree_hist_diff[clus_file.stem] = hist_aux
	

#converting into dataframe to have a clear layout of the data
mintree_hist_df = pd.DataFrame(mintree_hist_list)
mintree_hist_df.index +=1

mintree_hist_diff_df = pd.DataFrame(mintree_hist_diff)
mintree_hist_diff_df.index +=1

c_diff = dict()

for k,v in mintree_hist_diff.items():
	c_diff[k]=sum(v)
#print(c_diff)


#saving histograms as csv file
mintree_hist_df.to_csv(Path(mintree_result,'histograms.csv'),sep='\t')
mintree_hist_diff_df.to_csv(Path(mintree_result,'histograms_difference.csv'),sep='\t')

with open(Path(mintree_result,'nof_genes.csv'),'w') as f:
	f.write(''+'\t'+'nof_genes'+'\t'+'difference\n')
	for k,v in n_genes_mintree.items():
		if k == clus_mintree.stem:
			f.write(k+'\t'+str(v)+'\n')
		else:
			f.write(k+'\t'+str(v)+'\t'+str(abs(n_genes_mintree[clus_mintree.stem] - v))+'\n')


#PARAMETERS MINSTREE DATASET

params_analysis = dict()
#params_analysis['dataset_mintree']=analysis(clus_mintree) #used only to check (false negative/positive is 0, precision is 1.0)
for clus_file in Path('gene_families','mintree','mintree').glob('*'):
	print(clus_file)
	params_analysis[clus_file.stem] = performance_analysis(clus_file,classes)
	params_analysis[clus_file.stem]['c-diff']= c_diff[clus_file.stem]
params_df = pd.DataFrame(params_analysis)
params_df.to_csv(Path(mintree_result,'parameters.csv'),sep='\t')	
print(params_df)


#----------------------#

#randomtree dataset analysis
n_genes_randomtree = dict()
randomtree_hist_list = dict()
randomtree_hist_diff = dict()

analysis_result = Path('analysis')
analysis_result.mkdir(exist_ok=True)

#dataset mintree analysis
randomtree_result = Path(analysis_result,'analysis_randomtree')
randomtree_result.mkdir(exist_ok=True)

clus_randomtree = Path('input_datasets', 'dataset_randomtree','randomtree','randomtree.clus')

#analyze dataset_mintree.clus, used as reference to compare other pangenome softwares

#classes are the families of genes, nof_seqs is the number of sequences/genes
classes, nof_seqs = getClasses(clus_randomtree)

#histogram in the form of dictionary (only values that occour at least once)
randomtree_gen_distr = getGenomeDistribution(clus_randomtree)

#number of genomes
genomes_list = set()
for key, value in classes.items():
	for gene in value:
		genomes_list.add(gene.split('_')[1])

nof_genomes = len(genomes_list)
print('@',nof_genomes)

#create the histogram
hist_randomtree = [0]*nof_genomes
for k in randomtree_gen_distr:
	hist_randomtree[k-1] = randomtree_gen_distr[k]

#saving data
n_genes_randomtree[clus_randomtree.stem] = nof_seqs
randomtree_hist_list[clus_randomtree.stem] = hist_randomtree

#calculate histograms for all the clus files obtained with pangenome softwares
for clus_file in Path('gene_families','randomtree','randomtree').glob('*'):
	print('rand',clus_file)

	f_classes, f_nof_seqs = getClasses(clus_file)

	file_gen_distr = getGenomeDistribution(clus_file)
	
	#number of genomes
	f_genomes_list = set()
	for key, value in classes.items():
		for gene in value:
			f_genomes_list.add(gene.split('_')[1])

	f_nof_genomes = len(f_genomes_list)
		
	#create the histogram
	f_hist_file = [0]*f_nof_genomes
	#print(sorted(f_hist_file))
	#print(f_nof_genomes)
	for k in sorted(file_gen_distr):
		f_hist_file[k-1] = file_gen_distr[k]
	#print(hist_file)

	#saving
	n_genes_randomtree[clus_file.stem] = f_nof_seqs
	randomtree_hist_list[clus_file.stem] = f_hist_file

	hist_aux = list()
	for i in range(len(hist_randomtree)):
		hist_aux.append(abs(hist_randomtree[i] - f_hist_file[i]))
	""" print(hist_randomtree)
	print(f_hist_file)
	print(hist_aux) """
	randomtree_hist_diff[clus_file.stem] = hist_aux
	

#converting into dataframe to have a clear layout of the data
randomtree_hist_df = pd.DataFrame(randomtree_hist_list)
randomtree_hist_df.index +=1

randomtree_hist_diff_df = pd.DataFrame(randomtree_hist_diff)
randomtree_hist_diff_df.index +=1

c_diff = dict()

for k,v in randomtree_hist_diff.items():
	c_diff[k]=sum(v)
print(c_diff)


#saving histograms as csv file
randomtree_hist_df.to_csv(Path(randomtree_result,'histograms.csv'),sep='\t')
randomtree_hist_diff_df.to_csv(Path(randomtree_result,'histograms_difference.csv'),sep='\t')

with open(Path(randomtree_result,'nof_genes.csv'),'w') as f:
	f.write(''+'\t'+'nof_genes'+'\t'+'difference\n')
	for k,v in n_genes_randomtree.items():
		if k == clus_randomtree.stem:
			f.write(k+'\t'+str(v)+'\n')
		else:
			f.write(k+'\t'+str(v)+'\t'+str(abs(n_genes_randomtree[clus_randomtree.stem] - v))+'\n')


#PARAMETERS MINSTREE DATASET

params_analysis = dict()
#params_analysis['dataset_randomtree']=analysis(clus_randomtree) #used only to check (false negative/positive is 0, precision is 1.0)
for clus_file in Path('gene_families','randomtree','randomtree').glob('*'):
	params_analysis[clus_file.stem] = performance_analysis(clus_file,classes)
	params_analysis[clus_file.stem]['c-diff']= c_diff[clus_file.stem]
params_df = pd.DataFrame(params_analysis)
params_df.to_csv(Path(randomtree_result,'parameters.csv'),sep='\t')	
print(params_df)




			
		
