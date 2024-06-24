#!/bin/sh

for ((i=300; i<450; i++))
do
   sed -e "s/REPLACE/geosse${i}/g" geosse_C.xml > geosse_C${i}.xml
   beast geosse_C${i}.xml
done