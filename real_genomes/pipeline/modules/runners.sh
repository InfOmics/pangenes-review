#!/bin/bash

#$1 is the first argument and specifies the name of the software we want to run
#the other argv are the values we use to run pangenome softwares
#echo $1 #type of program to call ( for the case selection )
#echo "${@:2}" #all the others are parameters to use to run the program
cores=6 #number of cores to use for the pangenome programs
ulimit -m 60000000 #limit the ram to 60GB = 60 000 000 kbytes

case $1 in
	#alf) #$2 is the path to the alf-param.drw file containing the parameters to generate the genome
	#	/usr/bin/time -f"#memory# %M" bash alf/bin/alfsim $2
	#	echo $1: bash alf/bin/alfsim $2 >> log_running	
	#;;
	get_homologues)
		echo $1: perl $2 -d $3 -M -n $cores >> log_running #-M
		/usr/bin/time -f"#panreview# $1 %M %E" perl $2 -d $3 -M -n $cores
		
	;;
	micropan)
	echo $1: Rscript pangenome_softwares/micropan/run_micropan.R $2 $3 >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" Rscript pangenome_softwares/micropan/run_micropan.R $2 $3
		
	;;
	pandelos)
		echo $1: bash pangenome_softwares/Pandelos_java10/pandelos.sh $2 $3 >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" bash pangenome_softwares/Pandelos_java10/pandelos.sh $2 $3
		
	;;
	panget)
		PATH_input=$2/input
		FILE=$2/parameters.txt
		echo Input_directory	$PATH_input >> $FILE
		echo Output_directory	output >> $FILE
		echo Hvalue	0.42 >> $FILE
		echo Evalue	1e-5 >> $FILE
		echo Threads $cores >> $FILE
		
		echo $1: perl pangenome_softwares/panget/PanGet_blastp_modified.pl $FILE >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" perl pangenome_softwares/panget/PanGet_blastp_modified.pl $FILE
			
	;;
	panoct)
		panoct=pangenome_softwares/panoct/panoct_v3.23/bin/panoct.pl
		echo $1: $panoct -t $2 -f $3 -g $4 -P $5 -M Y -H Y -V Y -N Y -G y -T -c 0,50,75,100 >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" $panoct -t $2 -f $3 -g $4 -P $5 -M Y -H Y -V Y -N Y -G y -T
		
	;;
	panseq)
		echo $1: perl pangenome_softwares/Panseq/lib/panseq.pl $2 >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" perl pangenome_softwares/Panseq/lib/panseq.pl $2
		
	;;
	panx) #uses python2.7
		panx=pangenome_softwares/panX/panX.py
		diamond=pangenome_softwares/panX/diamond-linux64/diamond
		echo $1: python $panx -dmp $diamond -fn $2 -sl synthetic -t $cores -ct >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" python $panx -dmp $diamond -fn $2 -sl synthetic -t $cores -ct
		
	;;
	pgap)
		echo $1: perl pangenome_softwares/pgap/PGAP.pl --strains $2 --input $3 --output $4 --cluster --method GF --thread $cores >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" perl pangenome_softwares/pgap/PGAP.pl --strains $2 --input $3 --output $4 --cluster --method GF --thread $cores  
		
	;;
	roary)
		echo $1: roary -v -p $cores -f $2 "${@:3}" >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" roary -v -p $cores -f $2 "${@:3}"
		
	;;
	blast)
		echo $1: python3 modules/blast_module.py $2 $3 $cores >> log_running
		/usr/bin/time -f"#panreview# $1 %M %E" python3 modules/blast_module.py $2 $3 $cores
		
	
esac
