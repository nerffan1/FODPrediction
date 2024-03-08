# RDKit for BondPrediction
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
# Numpy
from scipy.spatial import distance
# Custom Made library
from globaldata import GlobalData
from ElementaryClasses import *
from Bond import *
from FOD import FOD
from BFOD import *

class Molecule:
    def __init__(self, xyzfile, RelaxedFODs = None) -> None:
        self.mAtoms: List[Atom] = []
        self.mFile = xyzfile
        self.mBonds = []
        # Load File
        if xyzfile[-3:] == "pdb":  # WiP. This is a test
            self.rdmol = Chem.MolFromPDBFile(xyzfile)
            self.__LoadPDB()
            exit
        elif xyzfile[-3:] == "xyz":
            self.rdmol = Chem.MolFromXYZFile(xyzfile)
            self.__LoadXYZ()
        GlobalData.mAtoms = self.mAtoms
        self.__RD_Bonds()
        self.CheckStericity()
        self.CalculateFODs()

        # Reverse Determine Parameters with a target file
        if RelaxedFODs != None:
            self.mRelaxPos = []
            self.__LoadTargetFODs(RelaxedFODs)
            self.ReverseDetermination()

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

    def ReverseDetermination(self):
        """
        This function executes the reversedetermination of paramters for all Target FODs that have
        been associated with the predicted FODs.
        """
        self.__AssociateTargets()
        # Bonding?
        for bfod in GlobalData.mBFODs:
            bfod.mAssocFOD.RevDet()
            print(50*'-')
            print(type(bfod))
            bfod.PrintParams()
            print(30*'*')
            bfod.mAssocFOD.PrintParams()

        # Free?
        # Core?
        from Shells import FODShell
        for atom in self.mAtoms:
            for shell in atom.mFODStruct.mCoreShells:
                # For a tetra, get average radii
                if shell.mShape == 'point':
                    core_dist = norm(atom.mPos - shell.mfods[0].mPos)
                    print(core_dist)
                elif shell.mShape == 'tetra':
                    targets = np.vstack([x.mAssocFOD.mPos for x in shell.mfods])
                    target_R = distance.cdist([atom.mPos], targets, 'euclidean')
                    # Pred. Stats. of Target Shell
                    shell.mTarget_u_R = np.mean(target_R)
                    shell.mTarget_s2_R = np.var(target_R)
                    # Pred. Stats. of Predicted Shell
                    pred_radii = [x.mR for x in shell.mfods]
                    shell.mPred_u_R = np.mean(pred_radii)
                    shell.mTarget_s2_R = np.var(pred_radii)

    #Private Methods
    def __AssociateTargets(self):
        """
        This function associates the target FOD that is nearest to the predicted FOD (i.e. the output of this program).
        """
        from FFOD import FFOD

        # Create a vertical vector
        c = np.vstack([x for x in self.mRelaxPos])
        for pfod in GlobalData.mFODs:
            # Get the minimum distance to Target FOD
            distances = distance.cdist([pfod.mPos],c, 'sqeuclidean')
            index = np.argmin(distances[0])
            # Print general information
            # Create appropriate associate fod
            if isinstance(pfod, BFOD):
                assoc = BFOD(
                    pfod.mBold,
                    pfod.mMeek,
                    c[index])
            elif isinstance(pfod, FFOD):
                assoc = FFOD(
                    pfod.mAtom,
                    target=c[index])
            elif isinstance(pfod, CFOD):
                assoc = CFOD(pfod.mAtom, c[index])
            else:
                print("Invalid classification for FOD")
            # Create the associate!
            pfod.mAssocFOD = assoc

            import logging
            logger = logging.getLogger('Logger')
            logger.setLevel(logging.DEBUG)
            logger.debug('Logger test')

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
        print(self.mComment) # Print the comment from the XYZ file 

    def __LoadTargetFODs(self, fods):
        """
        Load a file of FOD position to reversedetermine their parameters.
        """
        Target = open(fods, "r")
        # Read size of UP/DOWN FODs
        fods = Target.readline().split()
        upcount = int(fods[0])
        downcount = int(fods[1])
        # Load positions as ndarrays
        for i in range(upcount):
            coor = Target.readline().split()
            atom_xyz = np.array([float(x) for x in coor[0:3]])  # May have to convert to angstrom?
            self.mRelaxPos.append(atom_xyz)  # Name and Position
        Target.close()

    def __LoadPDB(self):
        """
        Credit to Betelgeuse in stackoverflow describing how to get this:
        https://stackoverflow.com/questions/71915443/rdkit-coordinates-for-atoms-in-a-molecule
        """
        for i, atom in enumerate(self.rdmol.GetAtoms()):
            position = self.rdmol.GetConformer().GetAtomPosition(i)
            pos = [position.x, position.y, position.z]
            self.mAtoms.append(Atom(i, atom.GetSymbol(), pos))

    def __RD_Bonds(self):
        """
        Calculate Bonds using the RDKit library.
        This will be used for prototyping  
        """ 
        #RDKit Functionality to determine the bonding.
        rdDetermineBonds.DetermineConnectivity(self.rdmol)
        rdDetermineBonds.DetermineBondOrders(self.rdmol, charge=0)
        rdmolops.Kekulize(self.rdmol)
        #Loop through the
        for rdbond in self.rdmol.GetBonds():
            
            atom1 = rdbond.GetBeginAtomIdx()   
            atom2 = rdbond.GetEndAtomIdx() 
            order = rdbond.GetBondTypeAsDouble()
            self.mBonds.append(Bond(atom1, atom2,order))
            self.mAtoms[atom1].AddBond(atom2, order)
            self.mAtoms[atom2].AddBond(atom1, order)
            # Add the Bonding FODs
            self.mAtoms[atom1].mFODStruct.mBFODs.append
            boldMeek = BoldMeek(self.mAtoms[atom1],self.mAtoms[atom2])
 
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
    def debug_printTargPred():
        c = np.vstack([x for x in self.mRelaxPos])
        for pfod in GlobalData.mFODs:
            # Get the minimum distance to Target FOD
            distances = distance.cdist([pfod.mPos],c, 'sqeuclidean')
            index = np.argmin(distances[0])
            if __debug__:
                print('-'*50)
                print(f"Index: {index}")
                print(f'Predicted: {pfod.mPos}')
                print(f'{type(pfod)}')
                print(f'Target: {c[index]}')
                print(f'Distance: {sqrt(distances[0][index])}')

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

    def _debug_printBFODs(self):
        from ElementaryClasses import Atom
        for atom in GlobalData.mAtoms:
            print(f'In atom {atom.mI}:')
            for bfod in atom.mFODStruct.mBFODs:
                print(bfod)

    def _debug_printBFODsXYZ(self):
        with open("out.xyz",'w') as output:
            #First 2 lines
            output.write(str(len(self.mAtoms) + len(GlobalData.mFODs)) + '\n')
            output.write(self.mComment)
            #Write all atoms
            for atom in self.mAtoms:
                output.write(' '.join([atom.mName,*[str(x) for x in atom.mPos]]) + '\n')
            #Write all FODs
            for bfod in GlobalData.mFODs:
                xyz = " ".join([str(x) for x in bfod.mPos])   
                output.write(f"X {xyz}\n")
        pass

    def _debug_CompareTargetFODs(self):
        """
        This function enumerates the number of predicted FODs (generated by FODLego)
        and the number of FODs from your target file.
        """
        print("-"*30)
        print(f'You have {len(GlobalData.mFODs)} Predicted FODs')
        print(f'You have {len(self.mRelaxPos)} Target FODs')
