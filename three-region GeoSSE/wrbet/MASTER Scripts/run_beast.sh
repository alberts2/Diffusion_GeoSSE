#!/bin/sh

for ((i=0; i<1000; i++))
do
   sed -e "s/REPLACE/geosse${i}/g" geosse.xml > geosse${i}.xml
   beast geosse${i}.xml
done