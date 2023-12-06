#!/usr/bin/python3
#Atom Class
#Version 0.1.6
from FODHeuristic import *
from  globaldata import GlobalData
import numpy as np
import sys

if len(sys.argv) == 1:
    print("No arguments were given, please provide XYZ file name")
    exit()
dat = GlobalData()
mol = Molecule(sys.argv[1])
#mol._debug_printAtoms()
mol.CreateXYZ()