#!/usr/bin/python3
#Atom Class
#Version 0.1.0
from  fod_struct import *
from  globaldata import GlobalData
import numpy as np

#RDKit for BondPrediction
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import AllChem

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
        self.mFODS = FODStructure(self)
        
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
    

class Bond:
    def __init__(self,start,end,order):
        self.mAtoms = (start,end)
        self.mOrder = order
    def __str__(self) -> str:
        return f"From {self.mAtoms[0]} to {self.mAtoms[1]}. Order: {self.mOrder}"

dat = GlobalData()
GlobalData._debug_samplenames()
mol = Molecule("Molecules_XYZ/test3.xyz")
mol._debug_printAtoms()
mol.CreateXYZ()
mol.DetermineFODs()