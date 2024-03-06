#!/bin/bash
# Description: Loop through all original files that have XYZ and FOD positions. 
# Read them and create a file with the UP/DOWN FODs from this XYZ files.
# Author: Angel-Emilio Villegas Sanchez

for file in $(ls *_og.xyz)
do
	# Write a file with the all positions
	echo $file
	target=$(echo $file | awk -F '_' '{print $1}')
	echo "Creating target ${target}_Target from $file"

	# Get UP count and positions
	sed -n '/X/p' $file | awk '{printf "%-15.10f %-15.10f %-15.10f \n", $2,$3,$4}' > tmp1
	up=$(wc -l < tmp1)

	# Get DOWN count and positions
	sed -n '/He/p' $file | awk '{printf "%-15.10f %-15.10f %-15.10f \n", $2,$3,$4}' > tmp2
	down=$(wc -l < tmp2)
	
	#Info
	echo up:$up - down:$down 

	#Create file 
	targetfile="${target}_Target"
	echo "$up $down" > $targetfile
	cat tmp1 tmp2 >> $targetfile
done
rm tmp1 tmp2
