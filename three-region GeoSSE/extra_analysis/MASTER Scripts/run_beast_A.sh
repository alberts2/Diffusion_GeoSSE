#!/bin/sh

for ((i=0; i<150; i++))
do
   sed -e "s/REPLACE/geosse${i}/g" geosse_A.xml > geosse_A${i}.xml
   beast geosse_A${i}.xml
done