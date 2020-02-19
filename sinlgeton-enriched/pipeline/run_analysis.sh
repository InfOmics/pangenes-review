#!/bin/bash

genomes="2 3"
types="root leaf"
percs="001 0005"


mkdir myco_results


for g in $genomes
do
for t in $types
do
for p in $percs
do
echo $g $t $p

odir="myco_results/${t}/${p}/${g}" 
mkdir -p $odir

cp gene_families_1/${t}/mycoplasma.g-${g}.synth-2000.lvp-${p}/* gene_families/${t}/mycoplasma.g-${g}.synth-2000.lvp-${p}/ 
cp gene_families_pandelos/${t}/mycoplasma.g-${g}.synth-2000.lvp-${p}/* gene_families/${t}/mycoplasma.g-${g}.synth-2000.lvp-${p}/ 

python3 myco_analysis.py datasets_complete/${t}/mycoplasma.g-${g}.synth-2000.lvp-${p}.${t}-50/ gene_families/${t}/mycoplasma.g-${g}.synth-2000.lvp-${p}/  > ${odir}/measures.csv

cp -r analysis ${odir}/


done
done
done

