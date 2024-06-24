#!/bin/bash
for ((i=0; i<7; i++))
do
   if [[ "$i" = "0" ]]
   then
         sed -e "s/"REPLACE"/"${i}"/g" geosse_all_A.R > geosse_all_A${i}.R
         Rscript geosse_all_A${i}.R &
   fi
   #
   if [[ "$i" = "1" ]]
   then
         sed -e "s/"REPLACE"/"${i}"/g" geosse_all_B.R > geosse_all_B${i}.R
         Rscript geosse_all_B${i}.R &
   fi
   #
   if [[ "$i" = "2" ]]
   then
         sed -e "s/"REPLACE"/"${i}"/g" geosse_all_C.R > geosse_all_C${i}.R
         Rscript geosse_all_C${i}.R &
   fi
   #
   if [[ "$i" = "3" ]]
   then
         sed -e "s/"REPLACE"/"${i}"/g" geosse_all_AB.R > geosse_all_AB${i}.R
         Rscript geosse_all_AB${i}.R &
   fi
   #
   if [[ "$i" = "4" ]]
   then
         sed -e "s/"REPLACE"/"${i}"/g" geosse_all_AC.R > geosse_all_AC${i}.R
         Rscript geosse_all_AC${i}.R &
   fi
   #
   if [[ "$i" = "5" ]]
   then
         sed -e "s/"REPLACE"/"${i}"/g" geosse_all_BC.R > geosse_all_BC${i}.R
         Rscript geosse_all_BC${i}.R &
   fi
   #
   if [[ "$i" = "6" ]]
   then
         sed -e "s/"REPLACE"/"${i}"/g" geosse_all_ABC.R > geosse_all_ABC${i}.R
         Rscript geosse_all_ABC${i}.R &
   fi
done