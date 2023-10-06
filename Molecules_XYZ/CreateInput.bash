#!/bin/bash

#Loop through all xyz files
for xyz in $(ls *.xyz); do
    echo $xyz
    newname=$(echo $xyz | awk -F '_' '{print $1}')
    echo $newname
done