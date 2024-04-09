#!/usr/bin/python3
#Atom Class
# Author: Angel-Emilio Villegas Sanchez

from  globaldata import GlobalData
from Molecule import *
from Analysis import *
from graphing import *
import sys

dat = GlobalData()
if len(sys.argv) == 1:
    print("No arguments were given, please provide XYZ file name")
    exit(1)

# Read a file and get parameters for multiple molecules/systems
elif sys.argv[1] == "list":
    assert len(sys.argv) == 3, "You did not provide a list of files to analyze."
    print("You provided the 'list' flag. The following file is expected to have several filenames")
    mols = CreateMolecules(sys.argv[2])
    #GridPropRatios(mols)
    #Histogram_Radii(mols)
    Histogram_Deviation(mols)
    #Angles_Hist(mols)

# Determine parameters just for one molecule
elif len(sys.argv) == 2:
    print("One argument passed. Creating FOD Prediction.")
    mol = Molecule(sys.argv[1])
    mol.CreateCLUSTER()
    mol.CreateFRMORB()
    mol._debug_printBFODsXYZ()

# Determine parameters and compare with a target file of FRMORBs
elif len(sys.argv) == 3:
    print("Two arguments passed. Reverse Determination of Relaxed FODs.")
    mol = Molecule(sys.argv[1], sys.argv[2])
    mol._debug_printBFODsXYZ()
