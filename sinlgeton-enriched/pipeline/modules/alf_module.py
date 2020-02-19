from pathlib import Path
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

#creates the alf parameters file
def alf_params_file(Nspecies, indels, dupl, loss, lgt, path):
	
	max_gap_len = 50 #maximum of indel length

	zipf_param = 1.821 #zipf distribution
	zipf = 'ZIPF,['+str(zipf_param)+']'
	#geometric distribution
	geom_param = 0.333 # #non ha spikes per gap lenght elevati
	geom = 'GEOM,['+str(geom_param)+']'
	#modello utilizzato
	model = geom

	file = [
		'SetRand(2345): # use this with any number, if you want reproducable results', 

		'webRequest := false;',
		"uuid := 's0-uuid';",

		'# name of simulation - you may want to change this',
		'mname := genome;',

		'# directories for file storage - you may want to change these',
		
		"wdir := './"+str(path)+"/'.mname.'/';",
		"dbdir := 'DB/';",
		"dbAncdir := 'DBancestral/';",

		'# time scale for simulation (PAM is default)',
		'unitIsPam := false:',

		'# parameters concerning the root genome',
		"realorganism := 'input/mycoplasma.g-1.db';", ##da sistemare perchè se no non legge nessun input

		'# parameters concerning the species tree',
		"treeType := 'BDTree';",
		'birthRate := 0.5;',
		'deathRate := 0.01;',
		'mutRate := 5000;',
		'NSpecies := '+str(Nspecies)+';',
		'ultrametric := false;',	#la parte tratteggiata dell'albero è per spaziare e rendere l'albero leggibile
		'scaleTree := false;'

		'# parameters concerning the substitution models',
		"substModels := [SubstitutionModel('WAG')];",

		#IndelModel(rate:non negative , model:string, parameters:list, maxLen:positive integer)
		'indelModels := [IndelModel('+str(indels)+','+model+','+str(max_gap_len)+')];',
		
		'rateVarModels := [RateVarModel()];', # impostato a zero , rappresenta il fatto che alcuni geni mutano piu di altri ( es. proteine di membrana )
		'modelAssignments := [1]:',
		'modelSwitchS := [[1]]:',
		'modelSwitchD := [[1]]:',

		'# parameters concerning gene duplication',
		'geneDuplRate := '+str(dupl)+';',
		'numberDupl := '+str(10)+';', #Maximum number of consecutive genes involved in one duplication event
		'transDupl := 0;', #Probability of a tranlocation after duplication (?) 

		#non usiamo fission e fusion
		'fissionDupl := 0;',
		'fusionDupl := 0;',

		'# parameters concerning gene loss',
		'geneLossRate := '+str(loss)+';',#Rate of gene losses (relative to substitutions)
		'numberLoss := '+str(10)+';',##Maximum number of consecutive genes involved in one loss event.

		'# parameters concerning LGT',
		'lgtRate := '+str(lgt)+';', #Rate of single lateral gene transfers (relative to substitutions)
		'orthRep := 0;', #Proportion of lateral gene transfers that are orthologous replacements
		'lgtGRate := 0.00004;', #Rate of lateral transfers of groups of genes
		'lgtGSize := 10;', #Maximum number of genes which can be transferred in one event

		'# parameters concerning rate heterogeneity among genes',
		"amongGeneDistr := 'Gamma';",
		'aGAlpha := 1;',

		'# select the output you want (apart from the species tree and genomes)',
		"simOutput := { 'GeneTrees' , 'VP', 'Fasta' , NULL }:"
	]

	print('\n'.join(file) , file=open(str(path)+'/alfparams.drw','w')) #mi conviene tenere il nome in alf-params.drw perchè così so che una volta raggiunta la cartella il file è sempre quello. Idem per i dati contenuti nella cartella genome


#create the dna sequences of the genomes from the aminoacidic sequences of genomes created with alf
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

	#questo serve nel caso ci fossero già dei risultati
	#lastfolder = sorted([i for i in path.glob('*')])[-1]
	
	
	path = Path(path,'genome','genome','DB')
	for filepath in path.glob('*'):
		if filepath.match('*_aa.fa'):
			
			to_fasta = list()
			for sequence in SeqIO.parse(filepath,'fasta'):
						
				aa_seq = list()
				for char in sequence.seq:
					aa_seq.append(aa_codon_table[char])
				
				sequence.seq = Seq(''.join(aa_seq), IUPAC.ambiguous_dna )
				to_fasta.append(sequence)
				
			outfile = str(filepath).replace('_aa','_dna')
			SeqIO.write(to_fasta,outfile,'fasta')
			print(outfile)

#metodo per ricavare il file familes.clus da ALF
import json
import tarfile
import networkx as nx
def alf2families(path, out_dir): #input path of the tar file, output directory where you save the family clus

	#extract VP.tgz ( it contains data about homologues, orthologues, paralogues and xenologues )
	tar = tarfile.open(Path(path,'VP.tgz'), "r:gz")
	tar.extractall(path=path)
	tar.close()

	print('\nGenerating graph...')
	G = nx.Graph() #inizializzo il grafo
	orthologues_filepath = Path(path,'VP','OP.drw')
	
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
	paralogues_filepath = Path(path,'VP','PP.drw')
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
	
	print('Gene Families found',len(list(nx.connected_components(G))),'\n')
	
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

