from pathlib import Path
from Bio import SeqIO
#create input file for pandelos, same folder of the genomes
def alf2idb(in_path, out_path):
	
	print('creating synthetic_idb.faa ( for Pandelos )...')
	
	files = sorted([i for i in Path(in_path,'protein').glob('*_aa.fasta')])
	with open(Path(out_path,'synthetic_idb.faa'),'w') as s :
		for fa in files:
			for genome in SeqIO.parse(fa,'fasta'):
				
				gene_id = genome.id.rstrip(',')
				genome_id = gene_id.split('_')[-1]
				seq = str(genome.seq)

				s.write(genome_id+'\t'+gene_id+'\t'+'hypothetical protein'+'\n')
				s.write(seq+'\n')
	print('...Done!')

	print('creating synthetic_idb.fna ( for Pandelos )...')
	files = sorted([i for i in Path(in_path,'dna').glob('*_dna.fasta')])
	with open(Path(out_path,'synthetic_idb.fna'),'w') as s :
		for fa in files:
			for genome in SeqIO.parse(fa,'fasta'):
				
				gene_id = genome.id.rstrip(',')
				genome_id = gene_id.split('_')[-1]
				seq = str(genome.seq)

				s.write(genome_id+'\t'+gene_id+'\t'+'hypothetical protein'+'\n')
				s.write(seq+'\n')
	print('...Done!')


#method to run pandelos
from modules.resource_control import call_program
def callPandelos(path_idb, gene_families_folder):
	#cmd = ['bash', Path('pangenome_softwares','Pandelos_java10','pandelos.sh') , path_idb, Path(gene_families_folder,'pandelos_families')]
	parameters = [ path_idb, Path(gene_families_folder,'pandelos_families') ]
	stat = call_program(parameters, 'pandelos')
	return stat
