#!/usr/bin/python3

import os


ifolder = 'gene_families/root'
ispecies = ['proch_cdb', 'mycoplasma_cdb','salmonella_cdb', 'xanthomonas_cdb', 'escherichia_cdb']
#ispecies = ['proch_cdb']

itools = ['get_homologues', 'pandelos', 'panoct', 'pgap', 'micropan', 'panget', 'panseq', 'panx', 'roary']
itools = ['pandelos']

for ispecie in ispecies:
	table = dict()
	for itool in itools:
		if itool == 'get_homologues':
			ifile = ifolder + os.sep + ispecie + os.sep + itool + '.clus'
		else:
			ifile = ifolder + os.sep + ispecie + os.sep + itool + '_families.clus'
		print(ifile)
		if os.path.isfile(ifile):
			table[itool] = dict()
			for line in open(ifile):
				cc = line.strip().split(' ')
				genomes = set()
				for c in cc:
					if itool == 'micropan':
						genome = c.split('_')[0]
					else:
						try:
							genome = c.split(':')[1].split('_')[1]
						except:
							genome = c.split(':')[1]
					genomes.add(genome)
				table[itool][len(genomes)] = table[itool].get(len(genomes),0)+1
	print('#\t', end='')
	print( '\t'.join( [i for i in sorted(itools) if i in table] ) )
	keys = set()
	for  d in table.values():
		keys |= d.keys()
	for key in sorted(keys):
		print(str(key)+'\t', end='')
		for tool in [i for i in sorted(itools) if i in table] :
			if key in table[tool]:
				print(str(table[tool][key])+'\t', end='')
			else:
				print('0\t', end='')
		print()
	print('-'*80)
