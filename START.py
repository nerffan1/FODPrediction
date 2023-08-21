#!/usr/bin/python3
#Atom Class
from FODHeuristic import *
from  globaldata import GlobalData

dat = GlobalData()
GlobalData._debug_samplenames()
mol = Molecule("Molecules_XYZ/test3.xyz")
mol._debug_printAtoms()
mol.CreateXYZ()
mol.DetermineFODs()