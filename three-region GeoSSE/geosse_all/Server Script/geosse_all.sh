#!/bin/bash
for ((i=0; i<10; i++))
do
   sed -e "s/"REPLACE"/"${i}"/g" geosse_all.R > geosse_all${i}.R
   Rscript geosse_all${i}.R &
done