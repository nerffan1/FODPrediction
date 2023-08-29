#RDKit for BondPrediction
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import AllChem

#Custom Made library
from globaldata import GlobalData
from fod_struct import Atom #For intellisense
from fod_struct import *

class Molecule:
    def __init__(self, xyzfile) -> None:
        self.mFile = xyzfile
        self.mComment = ""
        self.mAtoms = []
        self.mBonds = []
        self.__LoadXYZ()
        self.__RD_Bonds()
        self.CheckStericity()
        GlobalData.mAtoms = self.mAtoms
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
        Loop through the atoms and call respective methods to 
        calculate the FOD shells
        """
        for atom in self.mAtoms:
            atom.mFODStruct.PrepareShells(self.mAtoms)
            #Will change this feature soon
            atom.mFODStruct.FinalizeFODs()
    
    def CreateXYZ(self):
        """
        Create an XYZ file with
        """
        with open("output",'w') as output:
            #First 2 lines
            output.write(str(len(self.mAtoms) + self.CountFODs()))
            output.write(self.mComment + "(with calculated FODs)\n")
            #Write all atoms
            for atom in self.mAtoms:
                output.write(' '.join([atom.mName,*[str(x) for x in atom.mPos]]) + '\n')
            #Write all FODs
            for atom in self.mAtoms:
                for fod in atom.mFODStruct.mfods:
                    print([str(x) for x in fod])
                    xyz = " ".join([str(x) for x in fod])   
                    output.write(f"X {xyz}\n")

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
            atom_xyz = [float(x) for x in coor[1:4]]
            self.mAtoms.append(Atom(i, coor[0],atom_xyz)) #Name and Position
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
            order = bond.GetBondTypeAsDouble()
            self.mBonds.append(Bond(atom1, atom2,order))
            self.mAtoms[atom1]._AddBond(atom2, order)
            self.mAtoms[atom2]._AddBond(atom1, order)
    
    def CheckStericity(self):
        """
        Determine Steric number of all atoms. Currently assumes that the atoms are closed-shell
        TODO: Implement some way to easily add an open-shell calculation, in which there might be 
        open-shells
        """
        for atom in self.mAtoms:
            atom.CalcSteric()

    def CountFODs(self):
        """
        Returns the amout of FODs in the molecule
        """
        count = 0
        for atom in self.mAtoms:
            count += len(atom.mFODStruct.mfods)
        return count

    #Debugging Methods 
    def _debug_printAtoms(self):
        """Print atom names and positions"""
        for atom in self.mAtoms:  
            print("---------------------------------")
            print(atom.mName, "at", atom.mPos)
            print(f'Valency: {atom.mValCount}')
            print(f'Steric Number: {atom.mSteric}')
            print('BondedAtoms: ')
            for b in atom.mBonds:
                bonded = self.mAtoms[b.mAtoms[0]].mName
                print(f'-- Bonded to {bonded}({b.mAtoms[1]}). Order {b.mOrder}')
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
