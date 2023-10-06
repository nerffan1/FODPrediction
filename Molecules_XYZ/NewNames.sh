#!/bin/bash

#Loop through all xyz files
for in in $(ls *in); do
    sed -i "2i $in " $in
done