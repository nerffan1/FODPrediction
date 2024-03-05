#!/bin/bash

for file in $(ls *_in)
do
	new=$(echo $file | sed 's/_in/_in.xyz/')
	mv $file $new
done
