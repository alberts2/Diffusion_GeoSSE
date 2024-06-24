#!/bin/sh

for ((i=600; i<750; i++))
do
   sed -e "s/REPLACE/geosse${i}/g" geosse_AC.xml > geosse_AC${i}.xml
   beast geosse_AC${i}.xml
done