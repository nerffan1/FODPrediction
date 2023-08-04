#!/usr/bin/python3
#Atom Class

from numpy import array, dot, cross
from numpy.linalg import inv, det
from scipy.spatial.transform import Rotation

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
        self.mfods = []
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
    """Predict FOD positions"""
    def CalculateFODs(self):
        self.__CoreFODs()
        pass
    
    """Create an XYZ file with"""
    def CreateXYZ(self):
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

    #Private Methods
    """Load the """
    def __LoadXYZ(self):
        XYZ = open(self.mFile, "r")
        count = int(XYZ.readline()) #Read Size
        self.mComment = XYZ.readline() #Comment
        for i in range(count):
            coor = XYZ.readline().split()
            self.mAtoms.append(Atom(coor[0],coor[1:4])) #Name and Position
        XYZ.close()
    
    """Calculate Bonds using the RDKit library.
    This will be used for prototyping  
    """
    def __RD_Bonds(self):
        mol = Chem.MolFromXYZFile(self.mFile) 
        rdDetermineBonds.DetermineConnectivity(mol)
        rdDetermineBonds.DetermineBondOrders(mol, charge=0)
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtomIdx()   
            atom2 = bond.GetEndAtomIdx() 
            self.mBonds.append(Bond(atom1, atom2,
                                   bond.GetBondTypeAsDouble()))
    
    """Populate Core FODs for each atom
    Create 1s FODs first.
    """    
    def __CoreFODs(self):
        #Create 1s for all atoms
        for atom in self.mAtoms:
            self.mfods.append(FOD(atom, atom.mPos))
    
    """Create additional FODs"""
    def DetermineFODs(self):
        for bond in self.mBonds:
            print(bond)

    #Debugging Methods 
    """Print atom names and positions"""
    def _debug_printAtoms(self):
        for atom in self.mAtoms:  
            print(atom.mName, "at", atom.mPos)
    
    def _debug_printBondMatrix(self):
        print("##BOND MATRIX##")
        str4 = [[0] * len(self.mAtoms) for _ in range(len(self.mAtoms))]
        for b in self.mBonds:
            str4[b.mAtoms[0]][b.mAtoms[1]] = b.mOrder
        for atom in str4:
            print(atom)
    
class Atom:
    def __init__(self, mName, mPos) -> None:
        self.mName = mName
        self.mPos = mPos 
        self.mPeriod = self.__DeterminePeriod()
        self.mCoreFod = []
        self.mValenceFod = []
        self.mfods = []
        self.mBondTo = []
    
    def SetBond(self, atom2):
        self.mBondTo = atom2
    
    def __DeterminePeriod():
        pass

    def __DetermineGroup():
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



mol = Molecule("Molecules_XYZ/test3.xyz")
mol._debug_printAtoms()
mol._debug_printBondMatrix()
mol.CreateXYZ()
mol.DetermineFODs()