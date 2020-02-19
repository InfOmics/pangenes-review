#!/bin/bash

#gid="NC_000912"
gid="NC_000908"
#gid="NC_002771"

idir="input/${gid}/"
odir="output/${gid}/"
mkdir -p ${odir}
python3 gbk2ig.py ${idir} ${odir}/${gid}.db


ifile="${odir}/${gid}.db"

vperc="0.01"
ofile="${odir}/${gid}_${vperc}.db"
parents="${odir}/${gid}_${vperc}.parent"

python3 generate.py $ifile $ofile $vperc | grep -P "^parent " > $parents

iseqs="$ofile"
nof="50"
ofile="${odir}/${gid}_${vperc}_roots_${nof}.db"

python3 get_roots.py $parents $iseqs 50 $ofile

ogff="${odir}/${gid}_${vperc}_roots_${nof}.gff3"
mkdir -p $ogff
python3 db2gff3.py $ofile $ogff



ofile="${odir}/${gid}_${vperc}_leafs_${nof}.db"

python3 get_leafs.py $parents $iseqs 50 $ofile

ogff="${odir}/${gid}_${vperc}_leafs_${nof}.gff3"
mkdir -p $ogff
python3 db2gff3.py $ofile $ogff
