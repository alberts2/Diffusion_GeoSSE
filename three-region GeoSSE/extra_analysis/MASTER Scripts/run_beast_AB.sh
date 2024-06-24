#!/bin/sh

for ((i=450; i<600; i++))
do
   sed -e "s/REPLACE/geosse${i}/g" geosse_AB.xml > geosse_AB${i}.xml
   beast geosse_AB${i}.xml
done