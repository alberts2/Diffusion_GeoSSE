#!/bin/bash
for ((i=0; i<10; i++))
do
   sed -e "s/"REPLACE"/"${i}"/g" wr_bet.R > wr_bet${i}.R
   Rscript wr_bet${i}.R &
done