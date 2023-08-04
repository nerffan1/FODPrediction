#!/usr/bin/python3

import rdkit
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

#Create an Atom class

        

mol = Chem.MolFromXYZFile("test3.xyz")
rdDetermineBonds.DetermineConnectivity(mol)
rdDetermineBonds.DetermineBondOrders(mol, charge=0)

#TEST
num_atoms = mol.GetNumAtoms()
bond_matrix = [[0] * num_atoms for _ in range(num_atoms)]
for bond in mol.GetBonds():
    atom1_idx = bond.GetBeginAtomIdx()
    atom2_idx = bond.GetEndAtomIdx()
    print(mol.GetAtomWithIdx(atom1_idx).GetSymbol(), 
          mol.GetAtomWithIdx(atom2_idx).GetSymbol(),
          atom1_idx,
          atom2_idx)
    print(mol.GetAtomWithIdx(atom1_idx).GetBonds()[0].GetEndAtom().GetSymbol())
    print(mol.GetBondBetweenAtoms(atom1_idx,atom2_idx).GetBondType())
    bond_order = bond.GetBondTypeAsDouble()
    bond_matrix[atom1_idx][atom2_idx] = bond_order
    bond_matrix[atom2_idx][atom1_idx] = bond_order

for row in bond_matrix:
    print(row)

size = mol.GetNumAtoms()
print(mol.GetBondBetweenAtoms(0,1).GetBondType())
print(mol.GetAtomWithIdx(0).GetSymbol())
print(mol.GetAtomWithIdx(1).GetSymbol())

#Testing
AllChem.EmbedMolecule(mol)
#AllChem.Compute2DCoords(mol)
print(Chem.MolToMolBlock(mol))
Chem.AssignStereochemistry(mol)    
d = Draw.rdMolDraw2D.MolDraw2DCairo(250, 200)
d.drawOptions().addStereoAnnotation = True
Chem.rdmolops.FindPotentialStereo(mol)
d.drawOptions().addAtomIndices = True
d.DrawMolecule(mol)
d.FinishDrawing()
d.WriteDrawingText('atom_annotation_1.png')