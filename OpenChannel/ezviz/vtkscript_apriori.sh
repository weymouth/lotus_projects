#!/bin/bash

N=$2
name=$1
scene=$3
ext=png
if [ $scene == 'apriori_iso_scene.ini' ]; then
	ext=vtk
fi
if [ $scene == 'apriorilam2_iso_scene.ini' ]; then
	ext=vtk
fi


files=`seq -f %04g 0 $N`

for s in $files
do
	echo "$s"
	/usr/local/DAAC/bin/ezVizGeneric $scene -input $name.t.$s -output frame-$s.$ext -v
done
echo "$name had $N files"


