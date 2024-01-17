#!/bin/bash
for ((i=0; i<10; i++))
do
   sed -e "s/"REPLACE"/"${i}"/g" wr_ext.R > wr_ext${i}.R
   Rscript wr_ext${i}.R &
done