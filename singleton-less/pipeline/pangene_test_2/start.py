
from pathlib import Path
from modules.pipelines import *
from shutil import move
import datetime
import pandas as pd 



def convert_time(timestring):
	if len(timestring.split(':'))>2:
		if "." in timestring:
			pt =datetime.datetime.strptime(timestring,'%H:%M:%S.%f')
		else:
			pt =datetime.datetime.strptime(timestring,'%H:%M:%S')
	else:
		if "." in timestring:
			pt =datetime.datetime.strptime(timestring,'%M:%S.%f')
		else:
			pt =datetime.datetime.strptime(timestring,'%M:%S')
	return pt

#method to recursively remove a folder
def delete_folder(pth) :
    for sub in pth.iterdir() :
        if sub.is_dir() :
            delete_folder(sub)
        else :
            sub.unlink()
    pth.rmdir()


if Path('log_running').exists():
	Path('log_running').unlink()
if Path('gene_families').exists():
	delete_folder(Path('gene_families'))
if Path('softwares_data').exists():
	delete_folder(Path('softwares_data'))
if Path('execution_stats').exists():
	delete_folder(Path('execution_stats'))


#the two datasets we want to use
pth_dataset_mintree = Path('datasets','dataset_mintree')
pth_dataset_randomtree = Path('datasets','dataset_randomtree')

#foler containing the data for the execution of the pangenome softwares
pth_software_data = Path('softwares_data')
pth_software_data.mkdir(exist_ok=True)

#gene families folder where we have the .clus files
pth_gene_families = Path('gene_families')
pth_gene_families.mkdir(exist_ok=True)

#temporary folder
pth_tmp = Path('tmp')
pth_tmp.mkdir(exist_ok=True)

#execution stats folder
stats_folder = Path('execution_stats')
stats_folder.mkdir(exist_ok=True)


datasets = {
	'mintree':Path('input_datasets','dataset_mintree'),
	'randomtree':Path('input_datasets','dataset_randomtree')
}

for dtset_type, pth in datasets.items():
	for data in pth.glob('*'):

		#software execution data ( memory, elapsed_time )
		stats= list()

		software_data = Path(pth_software_data, dtset_type, data.stem)
		software_data.mkdir(parents=True,exist_ok=True)

		gene_families = Path(pth_gene_families, dtset_type, data.stem)
		gene_families.mkdir(parents=True, exist_ok=True)

		Path(stats_folder,dtset_type,data.stem).mkdir(parents=True, exist_ok=True)

		#le statistiche di esecuzione poi vengono salvate su un file quindi non devo cambiare nulla,
		#mi basta definire dove salvare il file nell' open()
		
		#PANDELOS
		pandelos_stat = pandelos(data, gene_families, software_data)
		if pandelos_stat != None:
			stats += pandelos_stat
		else: #need to remove pandelos empty clus file when it gets halted before completing (only good result files are kept)
			p = Path(gene_families,'pandelos_families.clus')
			if p.exists():
				p.unlink()
		print()

		#PANX
		panx_stat = panx(data,gene_families,software_data)
		if panx_stat != None:
			stats += panx_stat
		print()

		#PANSEQ
		panseq_stat = panseq(data,gene_families,software_data)	
		if panseq_stat != None:
			stats += panseq_stat
		print()

		gethomologues_stat = gethomologues(data,gene_families,software_data)	
		if gethomologues_stat != None:
			stats += gethomologues_stat
		print()

		pgap_stat = pgap(data,gene_families,software_data)
		if pgap_stat != None:
			stats += pgap_stat
		print()

		panget_stat = panget(data,gene_families,software_data)	
		if panget_stat != None:
			stats += panget_stat
		print()

		roary_stat = roary(data,gene_families,software_data)	
		if roary_stat != None:
			stats += roary_stat
		print()

		#BLAST all vs all , input for panoct and micropan
		print('Running BLAST...')
		blast_stat = call_program([data,software_data], 'blast')
		print('...BLAST done!')
		path_blastall = Path(software_data, 'blastDB','blastall.out')

		if blast_stat != None:
			stats += blast_stat
			panoct_stat = panoct(data,gene_families,software_data)
			if panoct_stat != None:
				stats += panoct_stat

			micropan_stat = micropan(data,gene_families,software_data)
			if micropan_stat != None:
				stats += micropan_stat
		else:
			print('MESSAGE: BLAST was terminated, panoct and micropan will not be executed')
		





		#cambiare il tempo in numero di secondi
		if len(stats)>0: #if i have at least one result (no results if all softwares take more than 2h to compute)

			print('saving resources used for:')
			print('dataset type:',dtset_type)
			print('dataset name:',data.stem)

			###SISTEMARE L'OUTPUT IN MODO CHE I FILE VENGANO GENERATI DIRETTAMENTE NELLA CARTELLA CORRETTA
			#idea -> salvare come file csv in modo da poter utilizzare i dataframe per fare la media di
			#tutti i risultati

			stat_dict = dict()
			blast_ram = int()
			blast_time = int()
			
			for s in stats:
				d = s.split(' ')[1:]

				if d[0] != 'blast': #ram and time will be added to panoct and micropan ram and time
					stat_dict[d[0]]= dict()

				if d[0] == 'blast':
					blast_ram = int(d[1])
					pt = convert_time(d[2])
					blast_time = pt.second+pt.minute*60+pt.hour*3600
				
				elif d[0] == 'panoct' or d[0] == 'micropan':
					stat_dict[d[0]]['ram'] = max(blast_ram, int(d[1]))

					pt = convert_time(d[2])
					software_time = pt.second+pt.minute*60+pt.hour*3600
					stat_dict[d[0]]['time'] = blast_time + software_time
				
				else:
					stat_dict[d[0]]['ram'] = int(d[1])
					pt = convert_time(d[2])
					software_time = pt.second+pt.minute*60+pt.hour*3600
					stat_dict[d[0]]['time'] = software_time
			
			df_stats = pd.DataFrame(stat_dict)
			df_stats.to_csv(Path(stats_folder,dtset_type,data.stem,'running_data.csv'), sep='\t')
			break #test con solo 1 dataset mintree e 1 randomtree

