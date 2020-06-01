list="mycoplasma escherichia salmonella proch xanthomonas"

for l in $list
do
mkdir ${l}_db
mkdir ${l}_cdb
python3 gbk2dbfiles.py ${l}/ ${l}_db/
python3 naming_convertion.py ${l}_db/ ${l}_cdb/
done
