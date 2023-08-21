#Description: This file contains as set of classes that will implement the FOD Heuristic Solution 
#  following the paradigm of Object-Oriented Programming (OOP). Encapsulation for FOD_Structure, FODs,
#  FOD_Orbital, among other things will be included in this file.
# Roadmap: Use polymorphism for Closed Hybrid Shells (e.g. sp3).
#  -  Add Tetrahedra class and several attributes/methods to manipulate them
#  - In far future, somehow implement the triaugmented triangular prism that corresponds to sp3d5 ( 9 FODs, 18 electrons) 
#Author: Angel-Emilio Villegas S.
from  globaldata import GlobalData
import math3d
import scipy.spatial.transform 
from numpy.linalg import inv, det
import csv

#FOD STRUCTURE

class Bond:
    def __init__(self,start,end,order):
        self.mAtoms = (start,end)
        self.mOrder = order
    def __str__(self) -> str:
        return f"From {self.mAtoms[0]} to {self.mAtoms[1]}. Order: {self.mOrder}"
    
class Atom:
    def __init__(self, index: int, Name: str, Pos) -> None:
        #Known Attributes
        self.mName = Name
        self.mPos = Pos 
        self.mI = index
        self.mZ = GlobalData.GetElementAtt(self.mName, "AtomicNumber")
        self.mPeriod = int(GlobalData.GetZAtt(self.mZ, "Period" ))
        self.mGroup = int(GlobalData.GetZAtt(self.mZ, "Group" ))
        self.mValenceELec = self._FindValence()
        #Undetermined Attributes
        self.mSteric = 0
        self.mCharge = 0 #Temporarily zero for all atoms. Ionic atoms can be dealt with in the future
        self.mBondTo = []
        self.mFODStruct = FODStructure(self)
        self.mCompleteValence = False
        
    def _AddBond(self, atom2: int, order: int):
        self.mBondTo.append(Bond(self.mI, atom2,order))
        self.mSteric += 1

    def _FindValence(self):
        """
        This method finds the number of valence electrons by finding the 
        difference between the current Group and the next FullShell Group.
        """
        for shell in GlobalData.mClosedGroups: 
            if self.mGroup <= shell:
                if self.mGroup == shell:
                    return 0
                else:
                    return  (shell - self.mGroup)
                    
    def _CheckFullShell(self):
        """
        Check that the atom has a closed shell.
        Future: Add a variable that saves the info so that looping every time
        this is called (if called more than once) is unnecesary
        """
        checkshell = self.mValenceELec
        for bond in self.mBondTo:
            checkshell -= bond.mOrder
        if checkshell == 0:
            return True
        else:
            return False

class FODShell:
    def __init__(self):
        self.mShape = ''
        self.mFOD = []
        self.mFree = 0
        self.mBonded = 0
        pass

class Tetrahedron(FODShell):
    """
    Tetrahedron Class: FOD Geometry corresponding to sp3 "hybridization' geometry.
    Roadmap: There will  be different functions to create compound transformations of FODs (e.g. the base, or peak
    of the tetrahedron), and to rotate them in the proper direction as well.
    """
    def __init__(self):
        super().__init__()
        self.mShape = 'tetra'
    
    #Class Methods
    def CreateTetra(self):
        """
        This method arranges 4 FOD points in a tetrahedral form. Depending on the bonding, and free electrons,
        the direction vector will be different. 
        """
        pass

    def RotateTetra(self):
        pass

    def Direction(self): #You can base this direction
        pass

class FODStructure:
    def __init__(self, parent: Atom):
        self.mAtom = parent
        self.mfods = [] #All FODs, but are they necessary, more likley than yes, for debugging or future purposes
        self.mCore = [] #A list of FODShells, since there can be more than 1
                        # Should we create just the an xyz list, ordered for the shells
        self.mValence = [] #Only one shell really, a list of FODs
        self.mGeometry = [] # A list of strings that accounts for mCore (m shells) and mValence shell (1) shapes, m + 1 elements

    def DetermineShells(self):
        """
        This function will determine the creation of Core FODs, those that are not 
        related to bonding. 
        Comment: The scheme is easy in the first 3 periods of the Periodic
        Table, but it will become trickier ahead if Hybridization heuristics don't work.
        Currently it only works for closed shell calculations (V 0.1.0)
        """
        #Begin with atoms preceding the transition metals
        #Set 1s FOD, assume every atom will have it in the current iteration of code
        
        # 1S Orbitals
        self.oneS()
        #Following orbitals

        #electrons = atom.mZ + atom.mValenceELec 
        #for shellelecs in GlobalData.mLadder_3p:
        #   if (electrons-shellelecs) == 0:
        #        #TODO: Here Initialize the geometries of the closed shells
        #        pass
    
    def oneS(self):
        self.mfods.append(self.mAtom.mPos) #mfods is a list of positions
        self.mCore.append([0]) #This is the index of where the coordinates found in mfods
        self.mGeometry.append("point")


    def CreateShell(self) -> FODShell:
        pass