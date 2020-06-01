#!/bin/bash
cd datasets
bash run_examples.sh
cd ..
cd pipeline
mkdir datasets
mkdir datasets/root/
cp ../datasets/*_cdb datasets/root/


if [ -f "softwares_data" ]; then
echo "softwares already downloaded"
else 
echo "downloading softwares"
wget ncrnadb.scienze.univr.it/pangenes-software/pangenome_softwares_r.zip
unzip pangenome_softwares_r.zip
rm pangenome_softwares_r.zip
fi

python3 start.py
python3 get_curves.py

