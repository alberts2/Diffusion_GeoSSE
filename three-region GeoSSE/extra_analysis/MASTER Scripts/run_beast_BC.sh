#!/bin/sh

for ((i=750; i<900; i++))
do
   sed -e "s/REPLACE/geosse${i}/g" geosse_BC.xml > geosse_BC${i}.xml
   beast geosse_BC${i}.xml
done