#!/bin/bash

list="escherichia salmonella xanthomonas mycoplasma"


for s in $list
do
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "running $s 's example ..."
bash download.sh ${s}.list.txt ${s}
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#exit
done

unzip proch.zip
mv pro2 proch

bash gbk2files.sh

