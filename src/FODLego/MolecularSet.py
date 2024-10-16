from FODLego.Molecule import *
from typing import List
from FODLego.Bond import Bond
from FODLego.FFOD import DFFOD, SFFOD, TFFOD
import pickle

class MolecularSet:
    def __init__(self, mols: List[Molecule]):
       self.mMols = mols
    
    def GetBonds_alternate(self, order=None) -> List[Bond]:
        """
        Returns a Python list with all bonds across all object molecules.
        If order is specified, only bonds with that order are returned.
        If order is None, bonds of all orders are returned
        """
        bonds = []
        for mol in self.mMols:
            if mol.mValidStruct == True:
                if order is None:
                    bonds.extend(mol.mBonds)
                else:
                    for bond in mol.mBonds:
                        if bond.mOrder == order:
                            bonds.append(bond)
        return bonds

    def GetBonds(self, order=None) -> List[Bond]:
        """
        Returns a python list of all bonds across all molecules. Because they are ordered by index, then there the 'greater than' condition assures no repeats
        """
        bonds = []
        for mol in self.mMols:
            if mol.mValidStruct == True:
                for atom in mol.mAtoms:
                    for bond in atom.mBonds:
                        bonded_at = bond.mAtoms[1]
                        if bonded_at.mI  > atom.mI:
                            if order is None:
                                bonds.append(bond)
                            else:
                                if bond.mOrder == order:
                                    bonds.append(bond)
        return bonds

    def tally_ffod_edges(self, typ=FOD):
        """
        Tallies the distances of the FFODs and BFODs in all atoms.
        The atom is the key, while the observations are the value.
        """
        # Useful Variables
        tally = {}
        validmols = 0

        # Loop through molecules and create an array for each atom in the dictionary
        for index, mol in enumerate(self.mMols):
            if mol.mValidStruct == True:
                # Keep track of succesful molecules
                validmols += 1
                for at in mol.mAtoms:
                    edges = at.GetAssocEdges_B_F_FOD(typ)
                    if len(edges) > 0:
                        if 0.0 in edges:
                            print(f"Found 0.0 in pairdist for molecule {mol.mSrc}, atom {at.mI}")
                            print(edges)
                            exit()
                        if at.mZ not in tally:
                            tally[at.mZ] = edges
                        else:
                            tally[at.mZ] = np.concatenate((tally[at.mZ], edges))
            else:
                print(f"Molecule {index} not valid ")
                pass

        return tally, validmols

    def tally_ffod_bfod_radii(self, typ=FOD):
        """
        Tallies the distances of the FFODs and BFODs in all atoms.
        The atom is the key, while the observations are the value.
        """
        # Useful Variables
        ffod_tally = {}
        bfod_tally = {}

        # Loop through molecules and create an array for each atom in the dictionary
        for index, mol in enumerate(self.mMols):
            if mol.mValidStruct == True:
                for at in mol.mAtoms:
                    bfod_radii, ffod_radii = at.get_assoc_radii_b_f_fod(typ)

                   # BFODs
                    if len(bfod_radii) > 0:
                        if at.mZ not in bfod_tally:
                            bfod_tally[at.mZ] = bfod_radii
                        else:
                            bfod_tally[at.mZ] = np.concatenate((bfod_tally[at.mZ], bfod_radii))
                   # FFODs
                    if len(ffod_radii) > 0:
                        if at.mZ not in ffod_tally:
                            ffod_tally[at.mZ] = ffod_radii
                        else:
                            ffod_tally[at.mZ] = np.concatenate((ffod_tally[at.mZ], ffod_radii))
            else:
                print(f"Molecule {index} not valid ")

        return ffod_tally, bfod_tally

    def get_atoms_w_FFOD_Type(self, typ):
        """
        Returns a set of atoms that contains a specific type of FOD.
        """
        atoms = []
        for mol in self.mMols:
            if mol.mValidStruct:
                for at in mol.mAtoms:
                    if exists(typ, at.GetFFODs()):
                        atoms.append(at)
        return atoms

    def get_fod_deviations(self,type='all'):
        """
        Returns an array with the distance deviation of all FODs in the set.
        """
        devi = []

        # Loop through the molecules and get valence electrons
        for mol in self.mMols:
            if type == 'ffod':
                vfods = mol.mFFODs
            elif type == 'bfod':
                vfods = mol.mBFODs
            elif type == 'sbfod':
                vfods = mol.get_bfod(SBFOD)
            elif type == 'dbfod':
                vfods = mol.get_bfod(DBFOD)
            elif type == 'tbfod':
                vfods = mol.get_bfod(TBFOD)
            elif type == 'sffod':
                vfods = mol.get_ffod(SFFOD)
            elif type == 'dffod':
                vfods = mol.get_ffod(DFFOD)
            elif type == 'tffod':
                vfods = mol.get_ffod(TFFOD)
            else:
                vfods = set.union(mol.mBFODs, mol.mFFODs)

            if mol.mValidStruct:
                for vfod in vfods:
                    d = dist(vfod.mPos, vfod.mAssocFOD.mPos)
                    if not np.isnan(d):
                        devi.append(d)
        return devi

    def get_atoms_Z(self, Z:int):
        """
        Returns the set of atoms across the set with a certain Z.
        """
        atoms_z = []
        for mol in self.mMols:
            if mol.mValidStruct == True:
                atoms_z += [at for at in mol.mAtoms if at.mZ == Z]
        return atoms_z

    ##############  Printing ##############

    def print_names(self):
       for mol in self.mMols:
           print(mol.mSrc)

    ############## Pickle Variables ##############
    def saveset(self, filename='mols.pickle'):
        with open(filename, 'wb') as file:
            pickle.dump(self, file)
        print('You have saved', filename)

    @staticmethod
    def loadpickle(filename='mols.pickle'):
        with open(filename, 'rb') as file:
                return pickle.load(file)

    ############## Has ? Methods ################
    def hasSFFOD(self):
        print('SFFOD at:')
        for mol in self.mMols:
            mol.has_fodtype(SFFOD)
