#!/bin/bash
for ((i=0; i<10; i++))
do
   sed -e "s/"REPLACE"/"${i}"/g" compute_1repli_2reg.R > compute_1repli_2reg${i}.R
   Rscript compute_1repli_2reg${i}.R &
done