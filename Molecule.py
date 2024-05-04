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
    def __init__(self, source, RelaxedFODs = None) -> None:
        self.mAtoms: List[Atom] = []
        self.mFile = source
        self.mQ = 0
        self.mBonds = []
        self.mBFODs = set()
        self.mFFODs = set()
        self.mCFODs = set()

        # Load File
        if source[-3:] == "pdb":  # WiP. This is a test
            self.rdmol = Chem.MolFromPDBFile(source)
            self.__LoadPDB()
            exit
        elif source[-3:] == "xyz":
            self.rdmol = Chem.MolFromXYZFile(source)
            self.rdmol
            self.__LoadXYZ()
        else:    # SMILES
            self.mComment = source + '\n'
            m = Chem.MolFromSmiles(source)
            self.rdmol = Chem.AddHs(m)
            self.__LoadSMILES()
        Chem.rdMolDescriptors.CalcOxidationNumbers(self.rdmol)
        GlobalData.mAtoms = self.mAtoms
        self.__RD_Bonds()
        self.CheckStericity()
        self.CalculateFODs()
        self._QCFODs()

        # Reverse Determine Parameters with a target file
        if RelaxedFODs != None:
            self.mTargetFile = RelaxedFODs
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
            # Add the calculated FODs to the molecule
            for bfod in atom.mFODStruct.mBFODs:
                self.mBFODs.add(bfod)
            for ffod in atom.mFODStruct.mFFODs:
                self.mFFODs.add(ffod)
            for cfod in atom.mFODStruct.mCore:
                self.mCFODs.add(cfod)

    def _QCFODs(self):
        # Create a list with all FODs. Make into List for index use.
        allFODs = set.union(self.mBFODs, self.mFFODs)
        allFODs = list(allFODs)
        n = len(allFODs)
        # Double loop without double counting
        for i in range(n):
            fod1 = allFODs[i]
            for j in range(i+1,n):
                fod2 = allFODs[j]
                if dist(fod1.mPos, fod2.mPos) < 1:
                    print(f'Valence FOD at {fod1.mPos} is too close to FOD at {fod2.mPos}')

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
    
    def CreateCLUSTER(self):
        """
        Creates a CLUSTER file that will serve as an input file for FLOSIC to begin
        """
        cluster = open("CLUSTER", "w")
        # CLUSTER Preamble
        cluster.write("LDA-PW91*LDA-PW91\n")
        cluster.write("NONE\n")
        cluster.write(f"{len(self.mAtoms)}\n")
        # Loop thorugh each atom for coordinates
        for atom in self.mAtoms:
            coordinate = " ".join(str(x*1.88973) for x in atom.mPos) + " " + str(atom.mZ) + " ALL" + '\n'
            cluster.write(coordinate)
        cluster.write(f"{self.mQ} {0.0}")  # TODO: Make a variable that contains sum of all spins
        cluster.close()
    
    def CreateFRMORB(self):
        cluster = open("FRMORB", "w")
        # CLUSTER Preamble
        cluster.write(f"{len(GlobalData.mFODs)} 0\n")
        # Loop thorugh each atom for coordinates
        for fod in GlobalData.mFODs:
            coordinate = " ".join(str(x*1.88973) for x in fod.mPos) + '\n'
            cluster.write(coordinate)
        cluster.close()

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
        import FOD_Print
        self.__AssociateTargets()
        # Bonding?
        #for bfod in self.mBFODs:
            #FOD_Print.PrintSidebySide(bfod,bfod.mAssocFOD)

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
        First it checks for CFODs and removes them from the list (assuming we have the most accuracy on them), and then we start finding the FFODs and BFODs.
        Note: This is a bit messy.
        """
        from FFOD import FFOD

        # Create a vertical vector
        c = np.vstack([x for x in self.mRelaxPos])
        fods = self.mCFODs
        for pfod in fods:
            # Get the minimum distance to Target FOD
            distances = distance.cdist([pfod.mPos],c, 'sqeuclidean')
            index = np.argmin(distances[0])

            # Print general information
            # Create appropriate associate fod
            if isinstance(pfod, CFOD):
                pfod.mAssocFOD = CFOD(pfod.mAtom, c[index])
                # Exclude the FOD that has been associated from the relaxed list.
                c = np.delete(c,index,axis=0)
                print(c)
            else:
                print("Invalid classification for associated FOD")

        # Now do the bonding FODs
        fods = set.union(self.mBFODs, self.mFFODs)
        for pfod in fods:
            # Get the minimum distance to Target FOD
            distances = distance.cdist([pfod.mPos],c)
            index = np.argmin(distances[0])

            # Print general information
            # Create appropriate associate fod
            if isinstance(pfod, BFOD):
                pfod.mAssocFOD = BFOD(
                    pfod.mBold,
                    pfod.mMeek,
                    c[index])
            elif isinstance(pfod, FFOD):
                pfod.mAssocFOD = FFOD(
                    pfod.mAtom,
                    target=c[index])
            else:
                print("Invalid classification for associated FOD")

    def __LoadXYZ(self):
        """
        Load the XYZ file
        """
        XYZ = open(self.mFile, "r")
        count = int(XYZ.readline()) #Read Size
        # Extract comment and charge
        self.mComment = XYZ.readline()

        # If blank
        if self.mComment != '\n':
            self.mQ = int(self.mComment.split()[-1])
        
        # Fill atom information 
        for i in range(count):
            coor = XYZ.readline().split()
            atom_xyz = [float(x) for x in coor[1:4]]
            self.mAtoms.append(Atom(i, coor[0],atom_xyz)) #Name and Position
        XYZ.close()
        print(f"Comment in xyz file: {self.mComment}") # Print the comment from the XYZ file 

    def __LoadSMILES(self):
        """
        Creates a more appropriate molecule according to the "working with 3D Molecules section of the RDKit documentation.
        """
        # Prepare SMILES Molecule with rdkit
        AllChem.EmbedMolecule(self.rdmol, maxAttempts=6000)
        AllChem.MMFFOptimizeMolecule(self.rdmol)

        # Load onto FODLego scheme
        for i,atom in enumerate(self.rdmol.GetAtoms()):
            coor = self.rdmol.GetConformer().GetAtomPosition(i)
            name = atom.GetSymbol()
            self.mAtoms.append(Atom(i, name,coor)) #Name and Position

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
            if self.mTargetFile[-6:] == 'FRMORB':
                atom_xyz = np.array([float(x)*(1/1.88973) for x in coor[0:3]])
            else:
                atom_xyz = np.array([float(x) for x in coor[0:3]])
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
        rdDetermineBonds.DetermineBondOrders(self.rdmol, charge=self.mQ)
        rdmolops.Kekulize(self.rdmol)
        #Loop through the
        for rdbond in self.rdmol.GetBonds():
            at1 = rdbond.GetBeginAtomIdx()   
            at2 = rdbond.GetEndAtomIdx() 
            order = rdbond.GetBondTypeAsDouble()
            self.mBonds.append(Bond(at1, at2,order))
            self.mAtoms[at1].AddBond(at2, order)
            self.mAtoms[at2].AddBond(at1, order)
            # Add the Bonding FODs
            self.mAtoms[at1].mFODStruct.mBFODs.append
            boldMeek = BoldMeek(self.mAtoms[at1],self.mAtoms[at2])

    def __Checks(self):
        """
        Check the valency of the atoms.
        """
        for atom in self.rdmol.GetAtoms():
            if atom.GetTotalValence() > 4:
                pass
 
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

    def _debug_printAtoms(self):
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
            print(f'RD Bonds: {at.GetBonds()}')
            print(f'RD Total Degree of atom: {at.GetTotalDegree()}')
            print(f'RD Oxidation: {at.GetProp("OxidationNumber")}')
            
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