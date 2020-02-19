from pathlib import Path
from subprocess import call
import datetime
from Bio.Blast.Applications import NcbiblastpCommandline

#in questo script faccio avvenire 2 cose, la creazione del database e la ricerca con blast, i tempi
#di esecuzione sono calcolati come tempi di questo script che comprende entrambe le azioni

# creo un file temporaneo con le sequenze peptidiche ( copio la parte iniziale di panoct )
import sys

#peps_file = sys.argv[1]
#dbname = sys.argv[2]
#db_path = sys.argv[3]

path_seqs = sys.argv[1]
path_softwares_data = sys.argv[2]
cores = sys.argv[3]

print('path seqs:', path_seqs)
print('path softwares data:', path_softwares_data)
print('cores:',cores)

#creazione del file delle sequenze
from tempfile import TemporaryFile, NamedTemporaryFile
from Bio import SeqIO

peps_file = NamedTemporaryFile(dir='tmp')
with open(peps_file.name,'a') as tmp:
	for p in Path(path_seqs,'protein').glob('*_aa.fasta'):
		#print(p)
		#for line in open(p,'r'):
		#	tmp.write(line)
		buffer = list()
		for seq in SeqIO.parse(p,'fasta'):
			seq.id = seq.id.strip(',')
			seq.name = seq.name.strip(',')
			seq.description = ''
			buffer.append(seq)
		SeqIO.write(buffer,tmp,'fasta')


#parameter = path_seqs.parents[2].stem
blast_db = Path(path_softwares_data,'blastDB')
blast_db.mkdir(parents=True,exist_ok=True)
dbname = 'blastDB'
db_path = Path(blast_db,dbname)


#MAKEBLASTDB
makeblastdb = Path('pangenome_softwares','panoct','blast','makeblastdb')
#makeblastdb -in $peps -parse_seqids -dbtype prot -title exp_prot_db -out exp_prot_db
print('running makeblastdb...')
call([makeblastdb,'-in', str(peps_file.name), '-parse_seqids','-dbtype', 'prot', '-title', dbname, '-out', str(db_path)])
print('...makeblastdb done!')

#BLASTP
#blastall -m 9 -d exp_prot_db -i $peps -p blastp -o $bout
blastp = Path('pangenome_softwares','panoct','blast','blastp')

print('blastp STARTED at :',datetime.datetime.now())

out=Path(blast_db,'blastall.out')

cline = NcbiblastpCommandline(cmd=str(blastp), num_threads=cores, max_hsps=1, query=str(peps_file.name), db=str(db_path), remote=False, out=str(out), outfmt='7')

print('command :', cline)

cline()

peps_file.close()
print('blastp COMPLETED at :',datetime.datetime.now())
