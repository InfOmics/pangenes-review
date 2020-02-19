from pathlib import Path
from subprocess import call
import numpy as np
import pandas as pd
import sys

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

	#print(clusters[0])
	#print(dbclasses[   next(iter(dbclasses.keys()))  ])

	e_classes = set()
	for k,vv in dbclasses.items():
		for v in ( (i,j) for i in vv for j in vv if i != j) :
			e_classes.add(v)

	e_clusters = set()
	for vv in clusters:
		for v in ( (i,j) for i in vv for j in vv if i != j ) :
			e_clusters.add(v)

	
	tp = len( e_classes & e_clusters ) #True positives
	fp = len( e_clusters - e_classes) #False positives
	fn = len( e_classes - e_clusters) #False negatives
	tn = ((nof_seqs * nof_seqs) - nof_seqs) - len(e_classes) #True negatives
	precision = tp/(tp + fp) if (tp + fp) != 0 else 0
	recall = tp/(tp + fn) if (tp + fn) != 0 else 0
	result = {'classes_links':len(e_classes),
	'cluster_links':len(e_clusters),
	'TP':tp,
	'FP':fp,
	'FN':fn,
	'TN':tn,
	'precision':tp/(tp + fp) if (tp + fp) != 0 else '#DIV0',
	'recall':tp/(tp + fn) if (tp + fn) != 0 else '#DIV0',
	'true_negative_rate':tn/(tn + fp) if (tn + fp) != 0 else '#DIV0',
	'c-diff':None,
	'f1-score':2*(precision*recall)/(precision+recall) if (precision+recall) != 0 else 0}
	
	#labels = ['classes_links','cluster_links','TP','FP','FN','TN','precision','recall','true_negative_rate','c-diff']
	#result = [len(e_classes),len(e_clusters),tp,fp,fn,tn,tp/(tp + fp),tp/(tp + fn),tn/(tn + fp)]
	return result
	


##MAIN NEW


#mintree dataset analysis
n_genes_dataset = dict()
dataset_hist_list = dict()
dataset_hist_diff = dict()

analysis_result = Path('analysis')
analysis_result.mkdir(exist_ok=True)


dataset = Path(sys.argv[1])
gene_families = Path(sys.argv[2])

#clus_dataset = Path('datasets', 'root','mycoplasma.g-3.synth-2000.lvp-005.root-50','mycoplasma.g-3.synth-2000.lvp-005.root-50.fa.clus')
clus_dataset = list(dataset.glob('*.clus'))[0]

#classes are the families of genes, nof_seqs is the number of sequences/genes
dbclasses, nof_seqs = getClasses(clus_dataset)
""" for k,v in dbclasses.items():
	print(k,':',v)
quit() """

#histogram in the form of dictionary (only values that occour at least once)
dataset_gen_distr = getGenomeDistribution(clus_dataset)

#number of genomes
genomes_list = set()
for key, value in dbclasses.items():
	for gene in value:
		genomes_list.add(gene.split('_')[1])

nof_genomes = len(genomes_list)

#create the histogram
hist_dataset = [0]*nof_genomes
for k in dataset_gen_distr:
	hist_dataset[k-1] = dataset_gen_distr[k]

#saving data
n_genes_dataset[clus_dataset.stem] = nof_seqs
dataset_hist_list[clus_dataset.stem] = hist_dataset

#calculate histograms for all the clus files obtained with pangenome softwares
#for clus_file in Path('gene_families','root','mycoplasma.g-3.synth-2000.lvp-005').glob('*'):
for clus_file in gene_families.glob('*'):
	print(clus_file)
	
	f_classes, f_nof_seqs = getClasses(clus_file)
	file_gen_distr = getGenomeDistribution(clus_file)
	
	#number of genomes
	f_genomes_list = set()
	for key, value in dbclasses.items():
		for gene in value:
			f_genomes_list.add(gene.split('_')[1])

	f_nof_genomes = len(f_genomes_list)
		
	#create the histogram
	f_hist_file = [0]*f_nof_genomes
	for k in file_gen_distr:
		f_hist_file[k-1] = file_gen_distr[k]
	#print(hist_file)

	#saving
	n_genes_dataset[clus_file.stem] = f_nof_seqs
	dataset_hist_list[clus_file.stem] = f_hist_file

	hist_aux = list()
	for i in range(len(hist_dataset)):
		hist_aux.append(abs(hist_dataset[i] - f_hist_file[i]))
	""" print(hist_dataset)
	print(f_hist_file)
	print(hist_aux) """
	dataset_hist_diff[clus_file.stem] = hist_aux

""" for k,v in dataset_hist_diff.items():
	print(k,':',v) """



#converting into dataframe to have a clear layout of the data
dataset_hist_df = pd.DataFrame(dataset_hist_list)
dataset_hist_df.index +=1

dataset_hist_diff_df = pd.DataFrame(dataset_hist_diff)
dataset_hist_diff_df.index +=1

c_diff = dict()

for k,v in dataset_hist_diff.items():
	c_diff[k]=sum(v)
#print(c_diff)


#saving histograms as csv file
dataset_hist_df.to_csv(Path(analysis_result,'histograms.csv'),sep='\t')
dataset_hist_diff_df.to_csv(Path(analysis_result,'histograms_difference.csv'),sep='\t')

""" with open(Path(analysis_result,'nof_genes.csv'),'w') as f:
	f.write(''+'\t'+'nof_genes'+'\t'+'difference\n')
	for k,v in n_genes_dataset.items():
		if k == clus_dataset.stem:
			f.write(k+'\t'+str(v)+'\n')
		else:
			f.write(k+'\t'+str(v)+'\t'+str(abs(n_genes_dataset[clus_dataset.stem] - v))+'\n') """


#PARAMETERS DATASET

params_analysis = dict()

for clus_file in gene_families.glob('*'):
	print(clus_file.stem)
	if clus_file.stem != 'micropan_families':
		params_analysis[clus_file.stem] = performance_analysis(clus_file,dbclasses)
		params_analysis[clus_file.stem]['c-diff']= c_diff[clus_file.stem]
params_df = pd.DataFrame(params_analysis)
params_df.to_csv(Path(analysis_result,'parameters.csv'),sep='\t')	
print(params_df)





			
		
