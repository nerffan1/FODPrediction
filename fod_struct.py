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
        self.mPos = np.array(Pos) 
        self.mI = index
        self.mZ = GlobalData.GetElementAtt(self.mName, "AtomicNumber")
        self.mPeriod = int(GlobalData.GetZAtt(self.mZ, "Period" ))
        self.mGroup = int(GlobalData.GetZAtt(self.mZ, "Group" ))
        self.mValCount = self._FindValence()
        #Undetermined Attributes
        self.mSteric = 0
        self.mFreePairs = 0
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
        TODO: This can be saved in GlobalData
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
        TODO: This will assume that we are doing up to the n=3 shell, with sp3 hybridization
        TODO: Take into account free pairs of electrons
        """
        #Prepare the valence shell first, since it will help determine the
        # orientation of the inner shells 
        for bond in self.mAtom.mBonds:
            if bond.mOrder == 1:
                #TODO: Add repulsion logic of the form a/r^k_{ij}
                At1 = bond.mAtoms[0]
                At2 = bond.mAtoms[1]
                dx = (atoms[At2].mPos - atoms[At1].mPos)/2
                self.mValence.append(atoms[At1].mPos + dx)
            elif bond.mOrder == 2:
                #Add 2 FODs, perpendicular to other 2.
                #Probably need to check how many 
                #Hardest one
                pass
            elif bond.mOrder == 3:
                pass

        #Count core electrons and
        core_elec = self.mAtom.mZ - self.mAtom.mValCount
        if core_elec != 0:
            for shell in GlobalData.mGeo_Ladder[core_elec]:
                if shell == 'point':
                    self.mCore.append(self.Point())
                elif shell == 'tetra':
                    self.mCore.append(self.Tetrahedron())

    def FinalizeFODs(self):
        """
        Add all FODs in the FODStructure, to the atom, so that they can be
        printed when creating the XYZ output file
        """
        #Add the core FODs
        print(self.mAtom.mName, f'({self.mAtom.mI}) at {self.mAtom.mPos}')
        print('line 136', [x.mfods for x in self.mCore])
        print('line 137', [x.mShape for x in self.mCore])
        for shell in self.mCore:
            print("Core Shells (", shell.mShape, '): ', shell.mfods)
            shell.mfods *= 0.2
            shell.mfods += self.mAtom.mPos   
            if self.mfods == []:
                self.mfods = [shell.mfods]
            else:
                self.mfods = np.vstack((self.mfods,shell.mfods))  ###HOW TO concatenate FODs, easily
        #Add the Valence FODs
        print('Valence FODs: ', self.mValence)
        print('Before CAT: ',self.mfods)
        if self.mfods == []:
            self.mfods = self.mValence  
        else:
            self.mfods = np.vstack((self.mfods,self.mValence))
        print('All FODs(?): ',self.mfods)

    def AddFOD(self):
        """
        Create a function that adds FOD information to the FOD Structure
        """

        pass
    
    class FODShell:
        def __init__(self, shape, fods):
            self.mShape = shape 
            self.mfods = np.array(fods)
        def __str__(self):
            return self.mShape

    class Point(FODShell):
        def __init__(self):
            super().__init__('point', [0.0,0.0,0.0])

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
