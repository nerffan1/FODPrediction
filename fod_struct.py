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
    def __init__(self, index: int, Name: str, Pos) -> None:
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
        self.mBondTo = []
        self.mFODStruct = FODStructure(self)
        self.mCompleteValence = False
        
    def _AddBond(self, atom2: int, order: int):
        self.mBondTo.append(Bond(self.mI, atom2,order))
        #self.mSteric += 1

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
        checkshell = self.mValCount
        for bond in self.mBondTo:
            checkshell -= bond.mOrder
        if checkshell == 0:
            return True
        else:
            return False
    
    def _CalculateStericity(self):
        """
        Calculates the steric number of the atom
        """
        if self.mSteric == 0:
            self.mSteric = (sum([bond.mOrder for bond in self.mBondTo]) + self.mValCount)/2
        else:
            print("The stericity has already been determined")

        #TODO: THIS SECTION NEED REVISION ASAP

class FODStructure:
    def __init__(self, parent: Atom):
        self.mAtom = parent
        self.mCore = [] #A list of Shells, since there can be more than 1
        self.mValence = [] #Only one shell really, a list of FODs
        self.mfods = [] #All FINALIZED FODs

    def PrepareShells(self, atoms):
        """
        This function will determine the creation of Core FODs, those that are not 
        related to bonding. 
        Comment: The scheme is easy in the first 3 periods of the Periodic
        Table, but it will become trickier ahead if Hybridization heuristics don't work.
        Currently it only works for closed shell calculations (V 0.1.0)
        """
        #Prepare the valence shell first, since it will help determine the 
        for bond in self.mAtom.mBondTo:
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
