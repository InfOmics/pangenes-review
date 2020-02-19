#!/bin/bash


configurations="_-n_100_-s_20_-indel_0.003_-loss_0.003_-dupl_0 \
_-n_100_-s_20_-indel_0.01_-loss_0.003_-dupl_0 \
_-n_100_-s_20_-indel_0.02_-loss_0.003_-dupl_0 \
_-n_100_-s_20_-indel_0.03_-loss_0.003_-dupl_0 \
_-n_100_-s_20_-indel_0.003_-loss_0.01_-dupl_0 \
_-n_100_-s_20_-indel_0.003_-loss_0.03_-dupl_0 \
_-n_100_-s_20_-indel_0.003_-loss_0.10_-dupl_0 \
_-n_100_-s_20_-indel_0.003_-loss_0.003_-dupl_0.001 \
_-n_100_-s_20_-indel_0.003_-loss_0.003_-dupl_0.01 \
_-n_100_-s_20_-indel_0.003_-loss_0.003_-dupl_0.05 \
_-n_100_-s_5_-indel_0.003_-loss_0.003_-dupl_0 \
_-n_100_-s_10_-indel_0.003_-loss_0.003_-dupl_0 \
_-n_100_-s_35_-indel_0.003_-loss_0.003_-dupl_0 \
_-n_100_-s_50_-indel_0.003_-loss_0.003_-dupl_0 
"


for configuration in $configurations
do
echo "-> "$configuration

mkdir $configuration
rm -r $configuration/*

pars=`echo "$configuration" | sed s/\_/\ /g`
echo $pars

python3 create_dataset.py $pars
mv alf_results $configuration/
mv input_datasets $configuration/

done
