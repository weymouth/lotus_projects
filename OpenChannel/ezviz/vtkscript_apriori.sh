#!/bin/bash

N=$2
name=$1
files=`seq -f %04g 0 $N`

for s in $files
do
	echo "$s"
	aprun -N 1 -n 1 /usr/local/DAAC/bin/ezVizGeneric apriori_scene.ini -input $name.t.$s -output frame-$s.png -v
done
echo "$name had $N files"


