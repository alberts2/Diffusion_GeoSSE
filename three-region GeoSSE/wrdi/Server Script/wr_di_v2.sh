#!/bin/bash
for ((i=0; i<10; i++))
do
   sed -e "s/"REPLACE"/"${i}"/g" wr_di_v2.R > wr_di_v2${i}.R
   Rscript wr_di_v2${i}.R &
done