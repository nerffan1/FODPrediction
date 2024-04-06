#!/bin/python3
from  globaldata import GlobalData
from matplotlib import pyplot as pp
from typing import Type

def graph2sp3():
    x = [x for x in GlobalData.mRadii[10]]
    y = [GlobalData.mRadii[10][z] for z in x ]
    fig, ax = pp.subplots()
    ax.plot(x,y)
    ax.scatter(x,y)
    ax.set_xlabel('Z (Atomic Number)',fontsize=20)
    ax.set_ylabel('Radii (Angstrom)', fontsize=20)
    ax.set_title('2SP3 Radius vs. Atomic Number',fontsize=20)
    ax.tick_params(axis='x', which='major', labelsize=18)  # Adjusting tick mark size on x-axis
    ax.tick_params(axis='y', which='major', labelsize=18)  
    #Label
    for i, j in zip(x, y):
        ax.annotate(str(i), (i, j), textcoords="offset points", xytext=(0, 10), ha='center',fontsize=10)
    #Save
    fig.set_figwidth(16)
    fig.set_figheight(9) 
    pp.savefig('testfig.svg', dpi=400)

def Histogram_Radii(molecule, atom_name: str):
    """
    Create a histogram for the radii of a certain atom with the name and list of molecules passed as arguments.
    """
    import Molecule
    # Initialize Data
    radii = []
    specified_atoms = []

    # Obtain all atoms from the molecules
    for mol in molecule:
        specified_atoms += [entry for entry in mol.mAtoms if entry.mName == atom_name]

    # Loop through atoms and extract distances
    for at in specified_atoms:
        for bfod in at.mFODStruct.:



