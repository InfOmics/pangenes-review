#!/bin/bash


l_s="100"
l_s=""

l_n="5 10 20 50 100"
l_n="5 10 50 100"
l_n="75 100"
l_n="35 50 75"
l_n=""

l_loss="0.003 0.03 0.1 0.01"
l_loss=""

l_indel="0.003 0.01 0.02"
l_indel="0.003 0.01 0.021 0.03"
l_indel="0.021"
l_indel="0.03"
l_indel=""

l_dupl="0.001 0.01 0.05 0.1"
l_dupl=""



d_s="100"
d_n="20"
d_loss="0.003"
d_indel="0.003"
d_dupl="0"



mkdir results


echo "################################################################################"
echo "LOSS"
echo "################################################################################"
for p_loss in ${l_loss}
do
	echo "########################################"
	echo ${p_loss}
	date
	dout="results/${d_s}_${d_n}_${d_indel}_${p_loss}_${d_dupl}" 
	rm $dout
	mkdir $dout
	echo "########################################"
	cd pangene_test_2
	echo "db creation"
	date
	python3 -u create_dataset.py -s ${d_s} -n ${d_n} -indel ${d_indel} -loss ${p_loss} -dupl ${d_dupl} &> create_dataset.run.log
	echo "running tools"
	date
	python3 -u start.py &> start.run.log
	echo "statistics"
	date
	python3 -u analysis.py &> analysis.run.log
	echo "data copy"
	date
	cd ..
	cp -r pangene_test_2/alf_results $dout
	cp -r pangene_test_2/input_datasets $dout
	cp -r pangene_test_2/gene_families $dout
	cp -r pangene_test_2/execution_stats $dout
	cp -r pangene_test_2/analysis $dout
	cp pangene_test_2/input_alf.drw $dout
	cp pangene_test_2/log_running $dout
	cp pangene_test_2/*.run.log $dout
	date

done


echo "################################################################################"
echo "INDEL"
echo "################################################################################"
for p_indel in ${l_indel}
do
	echo "########################################"
	echo ${p_indel}
	date
	dout="results/${d_s}_${d_n}_${p_indel}_${d_loss}_${d_dupl}" 
	rm $dout
	mkdir $dout
	echo "########################################"
	cd pangene_test_2
	echo "db creation"
	date
	python3 -u create_dataset.py -s ${d_s} -n ${d_n} -indel ${p_indel} -loss ${d_loss} -dupl ${d_dupl} &> create_dataset.run.log
	echo "running tools"
	date
	python3 -u start.py &> start.run.log
	echo "statistics"
	date
	python3 -u analysis.py &> analysis.run.log
	echo "data copy"
	date
	cd ..
	cp -r pangene_test_2/alf_results $dout
	cp -r pangene_test_2/input_datasets $dout
	cp -r pangene_test_2/gene_families $dout
	cp -r pangene_test_2/execution_stats $dout
	cp -r pangene_test_2/analysis $dout
	cp pangene_test_2/input_alf.drw $dout
	cp pangene_test_2/log_running $dout
	cp pangene_test_2/*.run.log $dout
	date

done




echo "################################################################################"
echo "DUPL"
echo "################################################################################"
for p_dupl in ${l_dupl}
do
	echo "########################################"
	echo ${p_dupl}
	date
	dout="results/${d_s}_${d_n}_${d_indel}_${d_loss}_${p_dupl}" 
	rm $dout
	mkdir $dout
	echo "########################################"
	cd pangene_test_2
	echo "db creation"
	date
	python3 -u create_dataset.py -s ${d_s} -n ${d_n} -indel ${d_indel} -loss ${d_loss} -dupl ${p_dupl} &> create_dataset.run.log
	echo "running tools"
	date
	python3 -u start.py &> start.run.log
	echo "statistics"
	date
	python3 -u analysis.py &> analysis.run.log
	echo "data copy"
	date
	cd ..
	cp -r pangene_test_2/alf_results $dout
	cp -r pangene_test_2/input_datasets $dout
	cp -r pangene_test_2/gene_families $dout
	cp -r pangene_test_2/execution_stats $dout
	cp -r pangene_test_2/analysis $dout
	cp pangene_test_2/input_alf.drw $dout
	cp pangene_test_2/log_running $dout
	cp pangene_test_2/*.run.log $dout
	date

done












echo "################################################################################"
echo "SUBPOPULATION"
echo "################################################################################"
for p_n in ${l_n}
do
	echo "########################################"
	echo ${p_n}
	date
	dout="results/${d_s}_${p_n}_${d_indel}_${d_loss}_${d_dupl}" 
	rm $dout
	mkdir $dout
	echo "########################################"
	cd pangene_test_2
	echo "db creation"
	date
	python3 -u create_dataset.py -s ${d_s} -n ${p_n} -indel ${d_indel} -loss ${d_loss} -dupl ${d_dupl} &> create_dataset.run.log
	echo "running tools"
	date
	python3 -u start.py &> start.run.log
	echo "statistics"
	date
	python3 -u analysis.py &> analysis.run.log
	echo "data copy"
	date
	cd ..
	cp -r pangene_test_2/alf_results $dout
	cp -r pangene_test_2/input_datasets $dout
	cp -r pangene_test_2/gene_families $dout
	cp -r pangene_test_2/execution_stats $dout
	cp -r pangene_test_2/analysis $dout
	cp pangene_test_2/input_alf.drw $dout
	cp pangene_test_2/log_running $dout
	cp pangene_test_2/*.run.log $dout
	date

done
