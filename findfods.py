#!/usr/bin/python3
#Atom Class
#Version 0.1.6
from  globaldata import GlobalData
from FODHeuristic import *
import numpy as np
import sys

if len(sys.argv) == 1:
    print("No arguments were given, please provide XYZ file name")
    exit()
dat = GlobalData()
mol = Molecule(sys.argv[1])
mol._debug_printBFODs()
#mol.debug_printAtoms()
mol._debug_printBFODsXYZ()
#mol.CreateXYZ()
