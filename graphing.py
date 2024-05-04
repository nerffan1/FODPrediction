#!/bin/python3
from  globaldata import GlobalData
from matplotlib import pyplot as pp
from Funcs import dist, AngleBetween
import numpy as np
from scipy.spatial.distance import cdist

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
    fig.set_figwidth(14)
    fig.set_figheight(5.5) 
    pp.savefig('testfig.svg', dpi=400)

def Histogram_Radii(molecules):
    """
    Create a histogram for the radii of a certain atom with the name and list of molecules passed as arguments.
    """
    import Molecule
    atomZ = input("What atom (Z) would you like to generate Radii Histogram for?")
    atomZ = int(atomZ)
    monoatomicR = GlobalData.mRadii[10][atomZ]
    
    # Initialize Data
    radii = []
    specified_atoms = []
    fig, ax = pp.subplots()

    # Obtain all atoms from the molecules
    for mol in molecules:
        specified_atoms += [entry for entry in mol.mAtoms if entry.mZ == atomZ]

    # Loop through atoms and extract distances
    for at in specified_atoms:
        for bfod in at.mFODStruct.mValence:
            radii += [dist(bfod.mAssocFOD.mPos,at.mPos)]
    
    print(len(radii))
    
    # Histogram
    # Create Histogram
    ax.hist(radii, bins=15, edgecolor='black', alpha=0.7, weights=np.ones(len(radii)) / len(radii))
    
    # Calculate average and plot a line
    average_radius = np.mean(radii)
    ax.axvline(average_radius, color='r', linestyle='dashed', linewidth=4, ymax=1, label='Average Radius')
    ax.axvline(monoatomicR, color='g', linestyle='dashed', linewidth=4, ymax=1, label='Monoatomic Radius')
    
    shift = 0.01
    # Markings
    ax.text(average_radius - 6*shift, 0.15, f'Average\nRadii: {average_radius:.2f} Å', color='r', verticalalignment='center', fontsize=18)
    ax.text(monoatomicR + shift, 0.15, f'Monoatomic\nRadii: {monoatomicR:.2f} Å', color='g', verticalalignment='center', fontsize=18)

    # Title
    ax.set_title('Carbon (Z=6) Radii', fontsize=25)
    
    # Labels
    ax.set_xlabel('Radius (Å)', fontsize=20)
    ax.set_ylabel('Density', fontsize=20)

    # Ticks
    ax.tick_params(axis='x', which='major', labelsize=15)
    ax.tick_params(axis='y', which='major', labelsize=15)  
    #ax.xaxis.set_ticks(np.arange(0, max, .1))  # Set tick positions

    #Size
    fig.set_figwidth(5)
    fig.set_figheight(4) 

    pp.show() 

def Histogram_Deviation(molecules):
    """
    Visualizes the deviation from prediction and optimized FODs. 
    """
    # Initialize Data
    devi = []
    fig, ax = pp.subplots()

    # Loop through atoms in all your molecules and find load the distances
    for mol in molecules:
        for vfod in mol.mBFODs:
            d = dist(vfod.mPos, vfod.mAssocFOD.mPos)
            if not np.isnan(d):
                devi.append(d)
        for vfod in mol.mFFODs:
            d = dist(vfod.mPos, vfod.mAssocFOD.mPos)
            if not np.isnan(d):
                devi.append(d)

    # Useful data
    min = np.min(devi)
    max = np.max(devi)
    mean = np.mean(devi)
    print(devi)
    
    # Create histogram
    size = len(devi) 
    bins = np.linspace(0, .6, 25)
    ax.hist(devi, bins=bins, edgecolor='black', alpha=0.7, weights=np.ones(size)/size)

    # Title
    ax.set_title(f' Distance deviation ({size} Observations)',fontsize=20)

    # Labels
    ax.set_xlabel('Deviation',fontsize=20)
    ax.set_ylabel('Density', fontsize=20)

    # Ticks
    ax.tick_params(axis='x', which='major', labelsize=15)
    ax.tick_params(axis='y', which='major', labelsize=15)  
    ax.xaxis.set_ticks(np.arange(0, 0.6, .05))  # Set tick positions

    # Limits
    ax.set_xlim(min,0.6)

    # Mean line and text
    ax.axvline(mean, color='red', linestyle='dashed', linewidth=3)
    pp.text(mean + .02, 0.13, f'Average: {mean:.2f} Å', color='r', verticalalignment='center', fontsize='17')

    #Size
    fig.set_figwidth(16)
    fig.set_figheight(9) 

    pp.show()

def Angles_Hist(molecules):
    # Initialize Data
    from BFOD import DBFOD
    from FFOD import DFFOD
    angles = []
    sumofangles = []
    fig, ax = pp.subplots()

    # Loop through atoms in all your molecules and find load the distances
    for mol in molecules:
        for atom in mol.mAtoms:
            countff = 0
            countbf = 0
            for vfod in atom.mFODStruct.mValence:
                if isinstance(vfod, DFFOD):
                    a = vfod.mAssocFOD.mPos - atom.mPos
                    b = vfod.mSiblings[0].mAssocFOD.mPos - atom.mPos
                    angles.append(np.rad2deg(AngleBetween(a,b)))
                    countff += 1
            if countbf + countff == 4:
                sumofangles.append(angles[-1] + angles[-2])
    print(len(angles))
    uniq = list(set(angles))
    print(len(uniq))
    print(uniq, np.mean(uniq))

def EdgeD_FFODs(molecules):
    # Loop through your various molecules
    tally = {}
    for mol in molecules:
        for at in mol.mAtoms:
            if len(at.mFODStruct.mFFODs) > 0:
                ffods = [x.mAssocFOD.mPos for x in at.mFODStruct.mFFODs]
                bfods = [x.mAssocFOD.mPos for x in at.mFODStruct.mBFODs]
                pairdist = cdist(ffods, bfods)
                if at.mZ not in tally:
                    tally[at.mZ] = pairdist.ravel()
                else:
                    tally[at.mZ] = np.concatenate((tally[at.mZ], pairdist.ravel()))
    tally = dict(sorted(tally.items(), key=lambda item: item[0]))
    print(tally)
    # Calculate mean and standard deviation for each atom
    means = []
    stds = []
    atoms = []
    for atom, distances in tally.items():
        atoms.append(atom)
        means.append(np.mean(distances))
        stds.append(np.std(distances))

    # Get edges for monoatomic atoms
    edges_monoatomic = []
    for at in atoms:
        if at < 11:
            edges_monoatomic.append(GlobalData.mVert[10][at])
        else:
            edges_monoatomic.append(GlobalData.mVert[18][at])

    # Create bar chart with error bars for atoms found in the data
    categories = [symbols[x] for x in atoms]
    fig, ax = pp.subplots()
    ax.bar(categories, means, yerr=stds, align='center', alpha=0.5, ecolor='black', capsize=10)

    # Plot empirical data as scatter points
    ax.scatter(categories, edges_monoatomic, label='Empirical', color='red', s=140)  # Adjust scatter marker size here

    # Customize plot
    ax.set_xlabel('Atomic Number')
    ax.set_ylabel('Mean Distance')
    ax.set_title('Mean Distance with Error Bars for Each Atom')
    ax.yaxis.grid(True)

    # Show plot
    pp.show()

symbols = {
1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar'
}