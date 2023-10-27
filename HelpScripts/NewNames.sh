#!/bin/bash
# Description: This script went several revisions. The original files did not have the '_og' postfix. This was added for quick reference.
# The '_in' postix was added to the same files with their FODs removed, as to compare outputs. As of 10/24/23, the data is not fully exemplary 
# of other systems. More systems with double and triple bonds are needed.
#Author: Angel-Emilio Villegas S.

#Loop through all xyz files
for in in $(ls *in); do
    sed -i "2i $in " $in
done
