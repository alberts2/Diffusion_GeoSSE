#!/bin/sh

for ((i=900; i<1050; i++))
do
   sed -e "s/REPLACE/geosse${i}/g" geosse_ABC.xml > geosse_ABC${i}.xml
   beast geosse_ABC${i}.xml
done