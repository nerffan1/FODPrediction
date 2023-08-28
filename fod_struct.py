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
import numpy as np 
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
    def __init__(self, index: int, Name: str, Pos):
        #Known Attributes
        self.mName = Name
        self.mPos = Pos 
        self.mI = index
        self.mZ = GlobalData.GetElementAtt(self.mName, "AtomicNumber")
        self.mPeriod = int(GlobalData.GetZAtt(self.mZ, "Period" ))
        self.mGroup = int(GlobalData.GetZAtt(self.mZ, "Group" ))
        self.mValCount = self._FindValence()
        #Undetermined Attributes
        self.mSteric = 0
        self.mCharge = 0 #Temporarily zero for all atoms. Ionic atoms can be dealt with in the future
        self.mBonds = []
        self.mFODStruct = FODStructure(self)
        self.mCompleteVal = False
        
    def _AddBond(self, atom2: int, order: int):
        self.mBonds.append(Bond(self.mI, atom2,order))

    def CalcSteric(self) -> None:
        self.mSteric = self.mValCount
        print(self.mName)
        for bond in self.mBonds:
            self.mSteric += bond.mOrder
        self.mSteric /= 2
    
    def _FindValence(self):
        """
        This method finds the number of electrons in the valence shell by finding the 
        difference between the current Group and the last ClosedShell Group. Only works up
        to 5th period.
        TODO: This can be saved in data
        """
        if self.mGroup < 4:
            return self.mGroup
        else:
            if self.mPeriod < 4:
                return (2 + (self.mGroup - 12))
            elif self.mPeriod < 6:
                return (self.mGroup)



                    
    def _CheckFullShell(self):
        """
        Check that the atom has a closed shell.
        Future: Add a variable that saves the info so that looping every time
        this is called (if called more than once) is unnecesary
        """
        checkshell = self.mValCount
        for bond in self.mBonds:
            checkshell -= bond.mOrder
        if checkshell == 0:
            return True
        else:
            return False
    
    # Additional Functions
    def __str__(self):
        pass    

######################## FOD Structure Class  ########################
class FODStructure:
    def __init__(self, parent: Atom):
        self.mAtom = parent
        self.mCore = [] #A list of FODShells
        self.mValence = [] #A list of FODs
        self.mfods = [] #All Finalized FODs

    def PrepareShells(self, atoms):
        """
        This function will determine the creation of Core FODs, those that are not 
        related to bonding. 
        Comment: The scheme is easy in the first 3 periods of the Periodic
        Table, but it will become trickier ahead if Hybridization heuristics don't work.
        Currently it only works for closed shell calculations (V 0.1.0)
        """
        #Prepare the valence shell first, since it will help determine the 
        for bond in self.mAtom.mBonds:
            if bond.mOrder == 1:
                #Add FOD in-between
                pass
            elif bond.mOrder == 2:
                #Add 2 FODs, perpendicular to other 2.
                #Probably need to check how many 
                #Hardest one
                pass
            elif bond.mOrder == 3:
                pass

        #Count core electrons and
        for shell_count in GlobalData.mGeo_Ladder:
            if self.mAtom.mZ < shell_count:
                core_elec = shell_count
        for shell in GlobalData.mGeo_Ladder[core_elec]:
            if shell == 'point':
                self.mCore.append(self.Point)
            elif shell == 'tetra':
                self.mCore.append(self.Tetrahedron)

    def FinalizeFODs(self):
        for shell in self.mCore:
            self.mfods.append(shell.mfods)


    def AddFOD(self):
        """
        Create a function that adds FOD information to the FOD Structure
        """

        pass
    
    class FODShell:
        def __init__(self, shape, fods):
            self.mShape = shape 
            self.mCenter = [0,0,0]
            self.mfods = fods

    class Point(FODShell):
        def __init__(self):
            super().__init__('point', [0,0,0])

    class Tetrahedron(FODShell):
        """
        Tetrahedron Class: FOD Geometry corresponding to sp3 "hybridization' geometry.
        Roadmap: There will  be different functions to create compound transformations of FODs (e.g. the base, or peak
        of the tetrahedron), and to rotate them in the proper direction as well.
        """
        def __init__(self):
            super().__init__('tetra', GlobalData.mTetraGeo)

        #Class Methods
        def CreateTetra(self):
            """
            This method arranges 4 FOD points in a tetrahedral form. Depending on the bonding, and free electrons,
            the direction vector will be different. 
            """
            pass

        def RotateTetra(self):
            pass
