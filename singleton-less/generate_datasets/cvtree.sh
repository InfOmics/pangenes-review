#!/bin/bash

mkdir dnatmp
rm -r dnatmp/*
cp $1/*.fasta dnatmp/

for f in `ls dnatmp/*.fasta`
do
echo $f
sed -i s/\>.*$/A/g $f
cat $f | tr -d '\n' > tmp
cat tmp | tr -d '\r' > $f
rm tmp
done

./parallel.out ./dnatmp/
