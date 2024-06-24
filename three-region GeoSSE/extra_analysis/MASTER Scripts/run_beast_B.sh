#!/bin/sh

for ((i=150; i<300; i++))
do
   sed -e "s/REPLACE/geosse${i}/g" geosse_B.xml > geosse_B${i}.xml
   beast geosse_B${i}.xml
done