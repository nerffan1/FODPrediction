#Description: This file contains as set of classes that will implement the FOD Heuristic Solution 
#  following the paradigm of Object-Oriented Programming (OOP). Encapsulation for FOD_Structure, FODs,
#  FOD_Orbital, among other things will be included in this file.
# Roadmap: Use polymorphism for Closed Hybrid Shells (e.g. sp3).
#  -  Add Tetrahedra class and several attributes/methods to manipulate them
#  - In far future, somehow implement the triaugmented triangular prism that corresponds to sp3d5 ( 9 FODs, 18 electrons) 
#Author: Angel-Emilio Villegas S.
from  globaldata import GlobalData
import scipy.spatial.transform 
from numpy.linalg import inv, det
import csv

#FOD STRUCTURE
class FODStructure:
    def __init__(self, atom):
        self.mCore = []
        self.mValence = []
        self.mfods = []
        self.mGeometry = []

    def DetermineCore(self,atom):
        """
        This function will determine the creation of Core FODs, those that are not 
        related to bonding. 
        Comment: The scheme is easy in the first 3 periods of the Periodic
        Table, but it will become trickier ahead if Hybridization heuristics don't work.
        Currently it only works for closed shell calculations (V 0.1.0)
        """
        #Begin with atoms preceding the transition metals
        #Set 1s FOD, assume every atom will have it in the current iteration of code
        #TODO: Add 1S here
        electrons = atom.mZ + atom.mValenceELec 
        for shellelecs in GlobalData.mLadder_3p:
            if (electrons-shellelecs) == 0:
                #TODO: Here Initialize the geometries of the closed shells
                pass

class Tetrahedron:
    def __init__(self) -> None:
        pass

class FOD:
    def __init__(self, parent, mPos = [0.0, 0.0, 0.0] ) -> None:
        self.mAtom = parent
        self.mPos = mPos
