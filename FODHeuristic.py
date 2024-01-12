#RDKit for BondPrediction
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem

#Custom Made library
from globaldata import GlobalData
from ElementaryClasses import *

class Molecule:
    def __init__(self, xyzfile) -> None:
        self.mFile = xyzfile
        self.rdmol = Chem.MolFromXYZFile(xyzfile)
        self.mComment = ""
        self.mAtoms: List[Atom] = []
        GlobalData.mAtoms = self.mAtoms
        self.mBonds = []
        self.__LoadXYZ()
        self.__RD_Bonds()
        self.CheckStericity()
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
        with open("out.xyz",'w') as output:
            #First 2 lines
            output.write(str(len(self.mAtoms) + self.CountFODs()) + '\n')
            output.write(self.mComment)
            #Write all atoms
            for atom in self.mAtoms:
                output.write(' '.join([atom.mName,*[str(x) for x in atom.mPos]]) + '\n')
            #Write all FODs
            for atom in self.mAtoms:
                for fod in atom.mFODStruct.mfods:
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
        rdDetermineBonds.DetermineConnectivity(self.rdmol)
        print(self.mComment)
        rdDetermineBonds.DetermineBondOrders(self.rdmol, charge=0)
        rdmolops.Kekulize(self.rdmol)
        for bond in self.rdmol.GetBonds():
            
            atom1 = bond.GetBeginAtomIdx()   
            atom2 = bond.GetEndAtomIdx() 
            order = bond.GetBondTypeAsDouble()
            self.mBonds.append(Bond(atom1, atom2,order))
            self.mAtoms[atom1].AddBond(atom2, order)
            self.mAtoms[atom2].AddBond(atom1, order)
            # Add the Bonding FODs
            self.mAtoms[atom1].mFODStruct.mBFODs.append
            boldMeek = BoldMeek(self.mAtoms[atom1],self.mAtoms[atom2])
            if order == 1:
                self.mAtoms[atom1].mFODStruct.mBFODs.append(SBFOD())
                self.mBFOD = SBFOD(*boldMeek)
            elif order == 2:
                self.mBFOD = DBFOD(*boldMeek)
            elif order == 3:
                self.mBFOD = TBFOD(*boldMeek)
 
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
    def debug_printAtoms(self):
        """Print atom names and positions"""
        for atom in self.mAtoms:  
            print("---------------------------------")
            print(atom.mName, "at", atom.mPos)
            at = self.rdmol.GetAtomWithIdx(atom.mI)
            print(f'RDValency: {at.GetTotalValence()}')
            print(f'RDImplicitVal: {at.GetImplicitValence()}')
            print(f'RDExplicitVal: {at.GetExplicitValence()}')
            print(f'RDFormalQ: {at.GetFormalCharge()}')
            print(f'RDNeighbors: {[x.GetSymbol() for x in at.GetNeighbors()]}')
            print(f'RDNeighbors: {at.GetHybridization()}')
            print(f'Valency: {atom.mValCount}')
            print(f'Steric Number: {atom.mSteric}')
            print(f'Free Pairs: {atom.mFreePairs}')
            print("Shell (Core) Structure:",*atom.mFODStruct.mCore)
            
            print('BondedAtoms: ')
            for b in atom.mBonds:
                bonded = self.mAtoms[b.mAtoms[1]].mName
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
