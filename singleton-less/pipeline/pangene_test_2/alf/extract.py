#from gbk file, extract a subset or all the CDS and save it as fasta file. Fasta file can be converted into darwin format with fasta2darwin
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature , FeatureLocation
from Bio.SeqRecord import SeqRecord


for seq in SeqIO.parse('NC_00098.gbk','genbank'):
	to_fasta = list()
	i = 0
	featurelist = seq.features
	cdsfeatures = list()
	#shuffle(featurelist)
	print(len(featurelist))
	for feature in featurelist:
		if feature.type == 'CDS':
		
			if i%1==0: #module 1 = all CDS , 5 = 1 CDS every five, etc...
				seqRec = feature.extract(seq)
				
				seqRec.seq = seqRec.seq.translate(stop_symbol='@')
				
				product = feature.qualifiers['product'][0] if 'product' in feature.qualifiers else 'hypothetical protein'
				seqRec.id = seq.name+'\t'+ feature.qualifiers['locus_tag'][0] +'\t'+ product
				seqRec.description = ''
				
				to_fasta.append(seqRec)
				#print(i)
				
			
			i += 1
	
	print('features estratte:',str(len(to_fasta)))
	seq.features = to_fasta
	file_out = seq.name+'extr'+str(len(to_fasta))+'_aa.fasta'
	SeqIO.write( to_fasta, file_out, 'fasta')
