#!/usr/bin/python3
#Atom Class
#Version 0.1.6
from  globaldata import GlobalData
from Molecule import *
import sys

dat = GlobalData()
if len(sys.argv) == 1:
    print("No arguments were given, please provide XYZ file name")
    exit()
elif len(sys.argv) == 2:
    print("One argument passed. Creating FOD Prediction.")
    mol = Molecule(sys.argv[1])
elif len(sys.argv) == 3:
    print("Two arguments passed. Reverse Determination of Relaxed FODs.")
    mol = Molecule(sys.argv[1], sys.argv[2])
    mol._debug_printBFODsXYZ()
