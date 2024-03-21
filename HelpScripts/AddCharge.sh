#!/bin/bash
# Description: This adds a neutral charge to all XYZ files
# Author: Angel-Emilio Villegas Sanchez

for file in $(ls *_in.xyz)
do
	sed -i '2s/$/ 0/' $file
done
