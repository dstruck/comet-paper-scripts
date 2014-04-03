cat *.csv > all.csv
sed -i "/name./d" all.csv
sed -i "/Sequence error/d" all.csv

