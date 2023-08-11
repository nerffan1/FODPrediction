#!/usr/bin/python3
#Atom Class
#Version 0.1.0

import numpy as np
import scipy.spatial.transform 
from numpy.linalg import inv, det
import csv

#RDKit for BondPrediction
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import AllChem

class GlobalData:
    
    def __init__(self):
        GlobalData.mElementInfo = self.LoadElements()
        GlobalData.mElemNames = self.LoadNames()
    
    #Public Functions
    def LoadElements(self):
        with open('elements2.0', mode ='r') as file:
            return np.genfromtxt(file, delimiter=',', encoding=None, dtype=None)

    def LoadNames(self):
        """
        Names are loaded in their abbreviated form, since XYZ files receive them this way
        """
        return GlobalData.mElementInfo[1:,2]

    def GetElementAtt(name: str, attr: str):
        """
        Gets the attribute found in the elements2.0 file

        Parameters:
        name (str): The name of the atom, in abbreviated form.
        attr (str): The attribute you want of chosen atom, found in the first row of elements2.0

        """

        if attr == "AtomicNumber":
            attrib_i = np.where(GlobalData.mElemNames == name)  
            # We must offset index by +1 because the array starts at zero instead of 1
            return attrib_i[0].item() + 1      
        else:
            element_i = np.where(GlobalData.mElemNames == name)[0]
            attrib_i = np.where(GlobalData.mElementInfo[0] == attr)[0]
            #We offset by 1 because the Info list has attribute names on row 1
            return GlobalData.mElementInfo[element_i + 1,attrib_i]
        
    def GetZAtt(Z: int, attr: str):
        attrib_i = np.where(GlobalData.mElementInfo[0] == attr)
        #We offset by 1 because the Info list has attribute names on row 1
        return GlobalData.mElementInfo[Z,attrib_i[0].item()]

    #Degugging tests
    def _debug_samplenames():
        for att in ["AtomicNumber", "Group", "Period", "Element","Metal", "NumberofShells","NumberofValence"]:
            print(f'{att}: {GlobalData.GetElementAtt("Ga", att)}')
            print(GlobalData.GetZAtt(31, att))

    #Class Variables    
    mElementInfo = []
    mElemNames = []
    mClosedGroups = [2,12,18]
    mLadder_3p = [4,10,12,18]
    #This scheme below assumes that a 1s Core FOD will already have been placed in the FODStruct 
    mLadder_3p_geometry = { 4: ['offcenter'], 10: ['tetra'], 12: ['tetra', 'offcenter'], 18: ['tetra', 'tetra'] }

dat = GlobalData()

class Molecule:
    def __init__(self, xyzfile) -> None:
        self.mFile = xyzfile
        self.mComment = ""
        self.mAtoms = [Atom]
        self.mAtoms.clear()
        self.mBonds = []
        self.mfods = [] #Replace this - exclude it to Atom classes, and loop over them
        self.__LoadXYZ()
        self.__RD_Bonds()
        self.CalculateFODs() 
    
    #String Output
    def __str__(self) -> str:
        str1 = self.mComment + "\n"
        str2 = f"Atoms: {len(self.mAtoms)}\n"
        str3 = f"Bonds: {len(self.mBonds)}\n"
        str4 = f"FODs: {len(self.mfods)}\n"
        return str1+str2+str3+str4 

    #Public Methods
    def CalculateFODs(self):
        """
        Predict FOD positions
        """
        self.__CoreFODs()
        pass
    
    def CreateXYZ(self):
        """
        Create an XYZ file with
        """
        with open("output",'w') as output:
            #First 2 lines
            output.write(str(len(self.mAtoms) + len(self.mfods)))
            output.write(self.mComment + "(with calculated FODs)\n")
            #Write all atoms
            for atom in self.mAtoms:
                output.write(' '.join([atom.mName,*atom.mPos]) + '\n')
            #Write all FODs
            for fod in self.mfods:
                output.write(' '.join(["X", *[str(f) for f in fod.mPos]]) + '\n')

    def ClosedMol(self):
        """
        Checks that all atoms in the molecule have a 
        """
        for atom in self.mAtoms:
            if atom._CheckFullShell() == False:
                return False
        return True

    #Private Methods
    def __LoadXYZ(self):
        """
        Load the XYZ file
        """
        XYZ = open(self.mFile, "r")
        count = int(XYZ.readline()) #Read Size
        self.mComment = XYZ.readline() #Comment
        for i in range(count):
            coor = XYZ.readline().split()
            self.mAtoms.append(Atom(coor[0],coor[1:4])) #Name and Position
        XYZ.close()
    
    def __RD_Bonds(self):
        """
        Calculate Bonds using the RDKit library.
        This will be used for prototyping  
        """
        mol = Chem.MolFromXYZFile(self.mFile) 
        rdDetermineBonds.DetermineConnectivity(mol)
        rdDetermineBonds.DetermineBondOrders(mol, charge=0)
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtomIdx()   
            atom2 = bond.GetEndAtomIdx() 
            self.mBonds.append(Bond(atom1, atom2,
                                   bond.GetBondTypeAsDouble()))
            self.mAtoms[atom1]._AddBond(atom2, bond.GetBondTypeAsDouble())
            self.mAtoms[atom2]._AddBond(atom1, bond.GetBondTypeAsDouble())
    
    def __CoreFODs(self):
        """Populate Core FODs for each atom
        Create 1s FODs first.
        """    
        #Create 1s for all atoms
        for atom in self.mAtoms:
            self.mfods.append(FOD(atom, atom.mPos))
    
    def DetermineFODs(self):
        """Create additional FODs"""
        for bond in self.mBonds:
            print(bond)

    #Debugging Methods 
    def _debug_printAtoms(self):
        """Print atom names and positions"""
        for atom in self.mAtoms:  
            print("---------------------------------")
            print(atom.mName, "at", atom.mPos)
            print("Group:", atom.mGroup)
            print(f'Valency: {atom.mValenceELec}')
            print(f'BondedAtoms: {atom.mBondTo}' )
            closedshell = atom._CheckFullShell()
            print (f'Shell Full: {closedshell}')
            if (closedshell == False): print ("###NONCLOSED SHELL SYSTEM###")

    def _debug_printBondMatrix(self):
        print("##BOND MATRIX##")
        str4 = [[0] * len(self.mAtoms) for _ in range(len(self.mAtoms))]
        for b in self.mBonds:
            str4[b.mAtoms[0]][b.mAtoms[1]] = b.mOrder
        for atom in str4:
            print(atom)
    
class Atom:
    def __init__(self, mName, mPos) -> None:
        #Atom Attributes
        self.mName = mName
        self.mPos = mPos 
        self.mZ = GlobalData.GetElementAtt(self.mName, "AtomicNumber")
        self.mPeriod = GlobalData.GetZAtt(self.mZ, "Period" )
        self.mGroup = int(GlobalData.GetZAtt(self.mZ, "Group" ))
        self.mValenceELec = self._FindValence()
        self.mBondTo = []
        self.mFODS = self.FOD_Structure(self)
        
    def _AddBond(self, atom2: int, order: int):
        self.mBondTo.append((atom2,order))

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
            checkshell -= bond[1]
        if checkshell == 0:
            return True
        else:
            return False

    #Attributes
    
    #FOD STRUCTURE
    class FOD_Structure:
        def __init__(self, atom):
            self.mCore = self.DetermineCore(atom) 
            self.mValence = []
            self.mfods = []

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

class FOD:
    def __init__(self, parent, mPos = [0.0, 0.0, 0.0] ) -> None:
        self.mAtom = parent
        self.mPos = mPos

class Bond:
    def __init__(self,start,end,order):
        self.mAtoms = (start,end)
        self.mOrder = order
    def __str__(self) -> str:
        return f"From {self.mAtoms[0]} to {self.mAtoms[1]}. Order: {self.mOrder}"

GlobalData._debug_samplenames()
mol = Molecule("Molecules_XYZ/test3.xyz")
mol._debug_printAtoms()
mol.CreateXYZ()
mol.DetermineFODs()