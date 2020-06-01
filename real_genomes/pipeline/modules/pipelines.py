#modules
from modules.alf_module import *
from modules.pandelos_module import *
from modules.panX_module import *
from modules.panseq_module import *
from modules.get_homologues_module import *
from modules.pgap_module import *
from modules.panget_module import *
from modules.roary_module import createRoaryInput, callRoary, roary2families
from modules.panoct_module import *
from modules.micropan_module import *
from modules.resource_control import call_program

#packages
from pathlib import Path
from subprocess import call

def pandelos(in_path, path_Gfamily, path_data):
	
	pth_pandelos_data = Path(path_data,'pandelos')
	pth_pandelos_data.mkdir(exist_ok=True)

	#creating idb input file
	alf2idb(in_path, pth_pandelos_data)
	
	#run Pandelos
	print('Running Pandelos...')
	stat = callPandelos(Path(pth_pandelos_data,'synthetic_idb.faa'), path_Gfamily)	
	print('...Pandelos done!')
	return stat

def panx(in_path, path_Gfamily, path_data):

	path_panx_data = Path(path_data,'panx','data','input_GenBank')
	path_panx_data.mkdir(parents=True, exist_ok=True)

	#create PanX input files
	createPanxInput(in_path, path_panx_data)

	#run PanX
	print('Running Panx...')
	stat = callPanX(path_panx_data.parent)
	print('...PanX done!')

	print('@vb@panxfamilies', path_panx_data.parent, path_Gfamily)
	panX2families(path_panx_data.parent, path_Gfamily)
	return stat

def panseq(in_path, path_Gfamily, path_data):
	print('Running Panseq...')
	path_panseq_data = Path(path_data,'panseq')
	path_panseq_data.mkdir(parents=True,exist_ok=True)
	
	#create input files and settings file
	createPanseqInput(in_path, path_panseq_data)
	createPanseqSettings(path_panseq_data)

	#run panseq		
	stat = callPanseq(Path(path_panseq_data,'settings.txt'))#,path_GfamilyIndels)
	
	#ricavare il file clus
	print('panseq2families')
	panseq2families(Path(path_panseq_data,'output'),path_Gfamily)
	print('...Panseq done!')
	
	return stat

def gethomologues(in_path, path_Gfamily, path_data):
	print('Running GET_HOMOLOGUES...')
	path_gethomologues_data = Path(path_data,'get_homologues')
	path_gethomologues_data.mkdir(parents=True,exist_ok=True)

	
	createGetHomologuesInput(in_path,path_gethomologues_data)
	stat = callGetHomologues(path_gethomologues_data)

	#extract gene families
	print('get_homologues2families')
	get_homologues2families(path_gethomologues_data, path_Gfamily)
	print('...GET_HOMOLOGUES done!')
	
	return stat

def pgap(in_path, path_Gfamily, path_data):
	print('Running PGAP...')
	#set output folder
	path_pgap_data = Path(path_data,'pgap')
	path_pgap_data.mkdir(parents=True,exist_ok=True)
	
	#prepare input
	createPgapInput(in_path, path_pgap_data)
	#run pgap
	stat = callPgap(path_pgap_data)
	
	#extract families from pgap data
	pgap2families(path_pgap_data,path_Gfamily)
	print('...PGAP done!')
	return stat

def panget(in_path, path_Gfamily, path_data):
	print('Running PANGET...')
	#panget data path
	path_panget_data = Path(path_data,'panget')
	
	#prepare input files
	createPangetInput(in_path, path_panget_data)
	#run panget
	
	stat = callPanget(path_panget_data)
	#extract families from panget data
	print('panget2families')
	panget2families(path_panget_data,path_Gfamily)
	print('...PANGET done!') 
	
	return stat

def roary(in_path, path_Gfamily, path_data):
	
	print('Running Roary...')
	#set output folder
	path_roary_data = Path(path_data,'roary')
	path_roary_data.mkdir(parents=True,exist_ok=True)

	#prepare input
	#input = pandelos idb files
	#out = folder where gff files are saved
	idb_file = Path(path_data,'pandelos','synthetic_idb.faa')
	print('creating roary input (.gff)')
	createRoaryInput(idb_file , path_roary_data)
	
	#run roary
	stat = callRoary(path_roary_data)

	#extract families from roary data
	roary2families(path_roary_data, path_Gfamily)
	print('...Roary done!')
	
	return stat

def panoct(in_path, path_Gfamily, path_data):
	print('Running PANOCT...')
	#set output folder
	path_roary_data = Path(path_data,'roary')
	path_panoct_data = Path(path_data,'panoct')
	path_panoct_data.mkdir(parents=True,exist_ok=True)

	#prepare input files from the input file (gff) used for roary
	createPanocInput(path_roary_data,path_panoct_data,in_path)
	
	#run panoct
	stat = callPanoct(path_panoct_data,path_data)
	#extract families from panoct data
	panoct2families(path_panoct_data,path_Gfamily)
	print('...PANOCT done!')
	
	return stat

def micropan(in_path, path_Gfamily, path_data):
	print('Running Micropan...')

	#set output folder
	path_panoct_data = Path(path_data,'panoct')
	path_micropan_data = Path(path_data,'micropan')
	path_micropan_data.mkdir(parents=True,exist_ok=True)

	#prepare input files from panoct blastall.out, micropan wants genome pairs blastp results
	#we can get them from the blastall.out file that is a blast all vs all
	createMicropanInput(path_data, path_micropan_data)
	
	#run micropan ( R library )
	stat = callMicropan(path_micropan_data)
	#extract families from micropan data
	micropan2families(path_micropan_data, path_Gfamily)
	print('...Micropan done!')
	return stat
