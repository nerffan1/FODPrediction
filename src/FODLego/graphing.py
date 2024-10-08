#!/bin/python3
from FODLego.MolecularSet import MolecularSet
from FODLego.globaldata import GlobalData
from matplotlib import pyplot as pp
from FODLego.Funcs import *
import numpy as np
from scipy.spatial.distance import cdist
from typing import List
from rdkit.Chem import PeriodicTable as PT
from rdkit.Chem import GetPeriodicTable
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
import matplotlib.transforms as mtransforms
from FODLego.MolecularSet import MolecularSet

symbols = {
1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar'
}
savepath='/home/angelemv/thesis_markdown/source/figures/'

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
    pp.show()

def Histogram_Radii(molecules):
    """
    Create a histogram for the radii of a certain atom with the name and list of molecules passed as arguments.
    """
    import Molecule
    import numpy as np
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(10,7))
    atoms = [[6,7],[8,16]]
    x_max = [[1.1, 1.2], [.9, 1]]
    bins = [[75,30],[100,20]]
    # Plot the histogram on each subplot
    for i in range(2):
        for j in range(2):
            ax = axs[i, j]
    
            def dist(pos1, pos2):
                return np.linalg.norm(np.array(pos1) - np.array(pos2))

            atomZ = atoms[i][j]
            vale = GetPeriodicTable().GetDefaultValence(atomZ)
            shell = vale + atomZ
            monoatomicR = GlobalData.mRadii[shell][atomZ]

            # Initialize Data
            radii = []
            specified_atoms = []

            # Obtain all atoms from the molecules
            for mol in molecules:
                if mol.mValidStruct == True:
                    specified_atoms += [entry for entry in mol.mAtoms if entry.mZ == atomZ]

            # Loop through atoms and extract distances
            for at in specified_atoms:
                for bfod in at.mFODStruct.mFFODs:
                    radii += [dist(bfod.mAssocFOD.mPos, at.mPos)]

            # Create Histogram
            ax.hist(radii, bins=np.arange(0,2,.02), edgecolor='black', alpha=0.7, weights=np.ones(len(radii)) / len(radii))

            # Calculate average and plot a line
            average_radius = np.mean(radii)
            ax.axvline(average_radius, color='r', linestyle='dashed', linewidth=3, ymax=1, label=f'Average R: {average_radius:.3f} Å')
            ax.axvline(monoatomicR, color='g', linestyle='dashed', linewidth=3, ymax=1,label=f'Monoatomic R: {monoatomicR:.3f} Å')

            shift = 0.03
            ax.legend(markerscale=0.5, fontsize=8)

            # Title
            atomName = GetPeriodicTable().GetElementName(atomZ)
            ax.set_title(f'{atomName} Valence FOD Radii ({len(radii)} Observations)', fontsize=9)

            # Labels
            ax.set_xlabel('Radius (Å)', fontsize=9)
            ax.set_ylabel('Density', fontsize=9)

            # Ticks
            ax.tick_params(axis='x', which='major', labelsize=9)
            ax.tick_params(axis='y', which='major', labelsize=9)
            ax.set_xlim(np.min(radii) - .1, x_max[i][j])

            # Printing info:
            print('Atom', f'Average Radius', 'MonoatomicR', 'Observations')
            print(atomName, f'{average_radius:.3f}', f'{monoatomicR:.3f}', f'{len(radii)}')

    # Show the plot
    pp.tight_layout()
    #fig.savefig(savepath + 'radii_select.svg')
    plt.show()

def get_mono_r(molset:MolecularSet):
    """
    Create a histogram for the radii of a certain atom with the name and list of molecules passed as arguments.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    atoms = [6,7,8,9,14,15,16,17]
    for atomZ in atoms:
        vale = GetPeriodicTable().GetDefaultValence(atomZ)
        shell = vale + atomZ
        monoatomicR = GlobalData.mRadii[shell][atomZ]

        # Initialize Data
        radii = []
        specified_atoms = molset.get_atoms_Z(atomZ)

        # Loop through atoms and extract distances
        for at in specified_atoms:
            for bfod in at.mFODStruct.mFFODs:
                radii += [dist(bfod.mAssocFOD.mPos, at.mPos)]

        average_radius = np.mean(radii)
        # Printing info:
        atomName = GetPeriodicTable().GetElementName(atomZ)
        print('Atom', f'Average Radius', 'MonoatomicR', 'Observations')
        print(atomName, f'{average_radius:.3f}', f'{monoatomicR:.3f}', f'{len(radii)}')

def Histogram_Deviation(mols: MolecularSet, ax):
    """
    Visualizes the deviation from prediction and optimized FODs. 
    """
    # Initialize Data
    devi = mols.get_fod_deviations()
    ffod_dev = mols.get_fod_deviations('ffod')
    bfod_dev = mols.get_fod_deviations('bfod')
 
    # Useful data
    min = np.min(devi)
    print(np.sort(devi)[-56:])
    max = np.max(devi)
    mean = np.mean(devi)
   
    # Create histogram
    size = len(devi)
    size2 = len(ffod_dev)
    size3 = len(bfod_dev)
    bins = np.linspace(0, 2, 160)
    #ax.hist(devi, bins=bins, edgecolor='black', alpha=1, weights=np.ones(size)/size, linewidth=0.5)
    ax.hist(bfod_dev, bins=bins, edgecolor='black', alpha=0.7, weights=np.ones(size3)/size, linewidth=0.5, label='BFODs')
    ax.hist(ffod_dev, bins=bins, edgecolor='black', alpha=0.6, weights=np.ones(size2)/size, linewidth=0.5, label='FFODs')

    # Title
    ax.set_title(f' Deviation of Valence FODs ({size} FODs)')

    # Labels
    ax.set_ylabel('Density', fontsize=12)

    # Ticks
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.tick_params(axis='y', which='major', labelsize=12)
    ax.xaxis.set_ticks(np.arange(0, max, .05))  # Set tick positions

    # Limits
    ax.set_xlim(-0.025,.8)

    # Mean line and text
    ffod_mean = np.mean(ffod_dev)
    ax.axvline(ffod_mean, color='darkorange', linestyle='dashed', linewidth=2, label=f'{ffod_mean:.2f} Å')
    bfod_mean = np.mean(bfod_dev)
    ax.axvline(bfod_mean, color='blue', linestyle='dashed', linewidth=2, label=f'{bfod_mean:.2f} Å')
    # Labels
    ax.legend(fontsize=12)

    pp.tight_layout()
    #fig.savefig(savepath + 'deviation.svg')
    return ax

def histogram_bfod_deviation(mols: MolecularSet, ax):
    """
    Visualizes the deviation from prediction and optimized FODs.
    """
    # Initialize Data
    devi = mols.get_fod_deviations()
    bfod1 = mols.get_fod_deviations('sbfod')
    bfod2 = mols.get_fod_deviations('dbfod')
    bfod3 = mols.get_fod_deviations('tbfod')

    # Useful data
    min = np.min(devi)
    max = np.max(devi)*.6
    mean = np.mean(devi)

    # Create histogram
    size1 = len(bfod1)
    size2 = len(bfod2)
    size3 = len(bfod3)
    size = size1 + size2 + size3
    bins = np.linspace(0, 0.5, 80)
    ax.hist(bfod1, bins=bins, edgecolor='black', alpha=1, weights=np.ones(size1)/size,
            linewidth=0.5, label='SBFODs')
    ax.hist(bfod2, bins=bins, edgecolor='black', alpha=1, weights=np.ones(size2)/size,
            linewidth=0.5, label='DBFODs', color='mediumblue')
    ax.hist(bfod3, bins=bins, edgecolor='black', alpha=0.9, weights=np.ones(size3)/size,
            linewidth=0.5, label='TBFODs', color='turquoise')

    # Title
    ax.set_title(f' Deviation of BFODs ({size} FODs)')

    # Labels
    ax.set_ylabel('Density',fontsize=12)

    # Ticks
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.tick_params(axis='y', which='major', labelsize=12)
    ax.xaxis.set_ticks(np.arange(0, max, .05))  # Set tick positions

    # Limits
    ax.set_xlim(-0.025,0.5)

    # Mean line and text
    ffod_mean = np.mean(bfod1)
    ax.axvline(ffod_mean, linestyle='dashed', linewidth=2, label=f'{ffod_mean:.2f} Å')
    bfod_mean = np.mean(bfod2)
    ax.axvline(bfod_mean, color='mediumblue', linestyle='dashed', linewidth=2, label=f'{bfod_mean:.2f} Å')
    bfod_mean = np.mean(bfod3)
    ax.axvline(bfod_mean, color='turquoise', linestyle='dashed', linewidth=2, label=f'{bfod_mean:.2f} Å')

    # Legend
    ax.legend(fontsize=12)

    pp.tight_layout()
    #fig.savefig(savepath + 'deviation.svg')
    return ax

def histogram_ffod_deviation(mols: MolecularSet, ax):
    """
    Visualizes the deviation from prediction and optimized FODs.
    """
    # Initialize Data
    devi = mols.get_fod_deviations()
    ffod1 = mols.get_fod_deviations('sffod')
    ffod2 = mols.get_fod_deviations('dffod')
    ffod3 = mols.get_fod_deviations('tffod')

    # Useful data
    min = np.min(devi)
    max = np.max(devi)*.6
    mean = np.mean(devi)

    # Create histogram
    size1 = len(ffod1)
    size2 = len(ffod2)
    size3 = len(ffod3)
    print('tffods: ', size3)
    size = size1 + size2 + size3
    bins = np.linspace(0, 0.5, 80)
    ax.hist(ffod2, bins=bins, edgecolor='black', alpha=0.9,
            weights=np.ones(size2)/size, linewidth=0.5,
            label='DFFODs', color='darkorange')
    ax.hist(ffod3, bins=bins, edgecolor='black', alpha=0.9,
            weights=np.ones(size3)/size, linewidth=0.5,
            label='TFFODs', color='saddlebrown')
    ax.hist(ffod1, bins=bins, edgecolor='black', alpha=0.8,
            weights=np.ones(size1)/size, linewidth=0.5,
            label='SFFODs', color='tan')

    # Title
    ax.set_title(f' Deviation of FFODs ({size} FODs)')

    # Labels
    ax.set_xlabel('Deviation (Angstrom)')
    ax.set_ylabel('Density',fontsize=12)

    # Ticks
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.tick_params(axis='y', which='major', labelsize=12)
    ax.xaxis.set_ticks(np.arange(0, max, .05))  # Set tick positions

    # Limits
    ax.set_xlim(-0.025,0.5)

    mean = np.mean(ffod1)
    ax.axvline(mean, color='darkorange', linestyle='dashed', linewidth=2, label=f'{mean:.2f} Å')
    mean = np.mean(ffod2)
    ax.axvline(mean, color='saddlebrown', linestyle='dashed', linewidth=2, label=f'{mean:.2f} Å')
    mean = np.mean(ffod3)
    ax.axvline(mean, color='tan', linestyle='dashed', linewidth=2, label=f'{mean:.2f} Å')

    # Legend
    ax.legend(fontsize=12)

    pp.tight_layout()
    #fig.savefig(savepath + 'deviation.svg')
    return ax

def deviation(mols):
    fig, ax = pp.subplots(3,1, figsize=(16,9))
    Histogram_Deviation(mols, ax[0])
    histogram_bfod_deviation(mols, ax[1])
    histogram_ffod_deviation(mols, ax[2])
    fig.savefig(savepath + 'deviation_all_pp.svg')
    #pp.show()

def Angles_Hist(mols: MolecularSet, z: int):
    """
    Create a graph with the addition of the DFFOD angles and the other FODs.
    """
    # Initialize Data
    from BFOD import DBFOD
    from FFOD import DFFOD
    bfod_angles = []
    ffod_angles = []
    fig, ax = pp.subplots(figsize=(12,6))

    bfod_tally = {}
    ffod_tally = {}
    addition = {}
    observed_atoms = set()
    
    # Loop through atoms in all your molecules and find load the distances
    atoms_dffod = mols.get_atoms_w_FFOD_Type(DFFOD)

    for atom in atoms_dffod:
        # FFODs
        observed_atoms.add(atom.mZ)
        ffod = atom.GetFFODs()[0]
        a = ffod.mAssocFOD.mPos - atom.mPos
        b = ffod.mSiblings[0].mAssocFOD.mPos - atom.mPos
        free_ang = np.rad2deg(AngleBetween(a,b))
        #ffod_angles.append(free_ang)

        # Free Tally
        if atom.mZ not in ffod_tally:
            ffod_tally[atom.mZ] = [free_ang]
        else:
            ffod_tally[atom.mZ].append(free_ang)

        # Find the angle between the other BFODs. THere should only be 2
        bfods = [x for x in atom.GetBFODs()]
        a = bfods[0].mAssocFOD.mPos - atom.mPos
        b = bfods[1].mAssocFOD.mPos - atom.mPos
        bond_ang = np.rad2deg(AngleBetween(a,b))
        #bfod_angles.append(bond_ang)
        
        # Bond Tally
        if atom.mZ not in bfod_tally:
            bfod_tally[atom.mZ] = [bond_ang]
        else:
            bfod_tally[atom.mZ].append(bond_ang)
            if atom.mZ not in addition:
                addition[atom.mZ] = [bond_ang + free_ang]
            else:
                addition[atom.mZ].append(bond_ang + free_ang)
        
    # Barchart information to create sequential bar charts for different atoms
    lines = 0
    for at in observed_atoms:
        # Data for specific atom
        bfod_angles += bfod_tally[at]
        ffod_angles += ffod_tally[at]
        lines += len(bfod_tally[at])
        startrect = 0
        colors = np.random.rand(len(observed_atoms), 3)  # Generates random RGB colors
        for i,at in enumerate(observed_atoms):
            frac = len(bfod_tally[at])/lines
            ax.add_patch(pp.Rectangle((startrect, 0), frac, ax.get_ylim()[1], transform=ax.transAxes, alpha=0.25, color=colors[i]))
            startrect += len(bfod_tally[at])/lines

    # Create a bar chart
    ctgs = range(lines)
    ax.bar(ctgs, bfod_angles, color='blue',label='BFOD Angle', align='edge')
    ax.bar(ctgs, ffod_angles, color='red',label='FFOD Angle', bottom=bfod_angles, align='edge')

    # fig2, ax2 = pp.subplots(2,1)
    # ax2[0].hist(addition[8], bins=30, linewidth=1, edgecolor='black')
    # ax2[1].hist(addition[16], bins=30, edgecolor='black')


    # Line for Average
    # ax.axhline(np.mean(addition), color='g', linestyle='dashed', linewidth=4, label=f'Mean: {np.mean(addition):3.3f}')
    
    # Add legend for data counts per atom
    # str_datums = ["Observations"] + [f'{symbols[atom]}: {len(bfod_tally[atom])} Observations' for atom in observed_atoms]
    # str_datums = '\n'.join(str_datums)
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # ax.text(1.06, 0.5, str_datums, transform=ax.transAxes, fontsize=14,
    #     verticalalignment='center', horizontalalignment='center', bbox=props)

    #Limits
    ax.set_xlim([0, lines])

    #Labels
    ax.set_title('BFOD and DFFOD Angle Addition', fontsize=15)
    ax.set_xlabel('Observation',fontsize=12)
    ax.set_ylabel('Sum of Angles', fontsize=12)
    ax.legend(loc='best')

    #Print stats
    for item, stuff in addition.items():
        print(f'{item}: {np.mean(stuff)} pm {np.std(stuff)} with {len(stuff)} observations')

    pp.tight_layout()
    pp.show()
    fig.savefig(savepath + 'anglesum.svg')

def EdgeDist_FFODs(mols:MolecularSet):
    from FFOD import SFFOD, DFFOD, TFFOD
    # Loop through your various molecules
    tally, validmols = mols.tally_ffod_edges(SFFOD)

    # Sort the values based on the first
    tally = dict(sorted(tally.items(), key=lambda item: item[0]))

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
    fig, ax = pp.subplots(figsize=(8,4.5))

    # Try to split background
    ax.yaxis.grid(True)
    ax.set_axisbelow(True)
    ax.add_patch(pp.Rectangle((0, 0), .5, ax.get_ylim()[1], transform=ax.transAxes, color='lime', alpha=0.07))
    ax.add_patch(pp.Rectangle((0.5, 0), .5, ax.get_ylim()[1], transform=ax.transAxes, color='coral', alpha=0.07))
    ax.bar(categories, means, yerr=stds, align='center', alpha=1, ecolor='black', capsize=10, color='slategray')

    # Plot empirical data as scatter points
    ax.scatter(categories, edges_monoatomic, label='Empirical', color='red', s=30)  # Adjust scatter marker size here
    for n, y in enumerate(edges_monoatomic):
        ax.hlines(y, n-0.4, n+0.4, color='r', alpha=0.8, linestyles='dashed')

    # Customize plot
    ax.set_xlabel('Atomic Number')
    ax.set_ylabel('Mean Distance')
    #ax.set_title('Mean Distance of BFOD and FFODs')

    # Add legend for data counts per atom
    str_datums = [f'{symbols[atom]}: {len(tally[atom])}' for atom in tally]
    str_datums = '\n'.join(str_datums)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.90, 0.95, str_datums, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

    # Statistics for each bar
    for i, (mean, std) in enumerate(zip(means, stds)):
            ax.text(i, 0.05, f"{mean:.2f} ± {std:.2f}",
                    horizontalalignment='center',
                    transform=ax.get_xaxis_transform(),
                    fontsize=9,
                    color='black',
                    rotation='vertical'
                    )

    # Show plot
    pp.tight_layout()
    fig.savefig(savepath + 'Edges.svg')
    pp.show()

def ffod_radii(molSet):
    from FFOD import SFFOD, DFFOD, TFFOD

    def tally_descriptive_stats(dictionary: dict):
        """
        Return the mean and standard of deviation for different elements
        in a dictionary.
        """
        means = {}
        stds = {}
        for atom, rad in dictionary.items():
            means[atom] = np.mean(rad)
            stds[atom] = np.std(rad)
        return means, stds

    # Helper Functions
    def CreateBarChart(ffod_tally, bfod_tally, colorbond='darkorange', colorfree='royalblue'):
        # Calculate mean and standard deviation for each atom
        atoms = []

        ffod_means, ffod_stds = tally_descriptive_stats(ffod_tally)
        bfod_means, bfod_stds = tally_descriptive_stats(bfod_tally)

        for atom in ffod_tally:
            atoms.append(atom)

        # Get edges for monoatomic atoms
        monoatomic_radii = []
        for at in atoms:
            if at < 11:
                monoatomic_radii.append(GlobalData.mRadii[10][at])
            else:
                monoatomic_radii.append(GlobalData.mRadii[18][at])

        # Create bar chart with error bars for atoms found in the data
        categories = [symbols[x] for x in atoms]
        fig, ax = pp.subplots()

        # Try to split background
        # ax.add_patch(pp.Rectangle((0, 0), .5, ax.get_ylim()[1], transform=ax.transAxes, color='lime', alpha=0.07))
        # ax.add_patch(pp.Rectangle((0.5, 0), .5, ax.get_ylim()[1], transform=ax.transAxes, color='coral', alpha=0.07))
        range = np.arange(0,len(categories))
        width = .4
        ffods = [f[1] for f in ffod_means.items()]
        ff_std = [f[1] for f in ffod_stds.items()]
        bfods = [f[1] for f in bfod_means.items()]
        bf_std = [f[1] for f in bfod_stds.items()]
        ax.bar(range - width/2, ffods, width, yerr=ff_std, align='center', alpha=1, ecolor='black', capsize=10, color=colorfree)
        ax.bar(range + width/2, bfods, width, yerr=bf_std, align='center', alpha=1, ecolor='black', capsize=10, color=colorbond)

        # Plot empirical data as scatter points
        ax.scatter(categories, monoatomic_radii, label='Empirical', color='red', s=140)  # Adjust scatter marker size here

        # Customize plot 1
        ax.set_xlabel('Atomic Number', fontsize=20)
        ax.set_ylabel('Mean Distance', fontsize=20)
        ax.yaxis.grid(True)
        # Customize plot
        ax.yaxis.grid(True)

        # Add legend for data counts per atom
        ffod_datums = [f'{symbols[atom]}: {len(sf[atom])}' for atom in sf]
        ffod_datums = '\n'.join(ffod_datums)
        bfod_datums = [f'{symbols[atom]}: {len(sb[atom])}' for atom in sb]
        bfod_datums = '\n'.join(bfod_datums)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, bfod_datums, transform=ax.transAxes, fontsize=20,
            verticalalignment='top', bbox=props)

        # Statistics for each bar
        for i,j in enumerate(atoms):
            mean = ffod_means[j]
            std = ffod_stds[j]
            ax.text(i - width/2, 0.15, f"{mean:.2f}\n± {std:.2f}", horizontalalignment='center', transform=ax.get_xaxis_transform(), fontsize=18, color='black')
            mean = bfod_means[j]
            std = bfod_stds[j]
            ax.text(i + width/2, 0.15, f"{mean:.2f}\n± {std:.2f}", horizontalalignment='center', transform=ax.get_xaxis_transform(), fontsize=18, color='black')

        return fig, ax

    def count_bars(*args: dict):
        """
        Count how many bars will there be for each atom Z.
        """
        tally = {}
        for arg in args:
            for atom in arg:
                if atom in tally:
                    tally[atom] += 1
                else:
                    tally[atom] = 1
        return tally

    def CreateBarChart2(a, b, c, d, e, f):
        # Calculate mean and standard deviation for each atom
        sf = tally_descriptive_stats(a)
        sb = tally_descriptive_stats(b)
        df = tally_descriptive_stats(c)
        db = tally_descriptive_stats(d)
        tf = tally_descriptive_stats(e)
        tb = tally_descriptive_stats(f)
        countbars = count_bars(a,b,c,d,e,f)

        atoms = [at for at in countbars]
        categories = [symbols[x] for x in atoms]

        # Get edges for monoatomic atoms
        monat_r = []
        for at in sorted(atoms):
            if at < 11:
                monat_r.append(GlobalData.mRadii[10][at])
            else:
                monat_r.append(GlobalData.mRadii[18][at])

        # Create bar chart with error bars for atoms found in the data
        categories = [symbols[x] for x in atoms]
        fig, ax = pp.subplots(figsize=(16,9))


        # Colors and labels for the different bars
        colors = ['bisque', 'darkorange', 'skyblue', 'royalblue', 'limegreen', 'forestgreen']
        labels = ['SFFOD', 'valence BFOD', 'DFFOD', 'valence BFOD', 'TFFOD', 'valence BFOD']

        # Try to split background
        for i,atom in enumerate(sorted(countbars)):
            barnumber = countbars[atom]
            width = 0.8/barnumber
            if barnumber != 1:
                barnumber /= -2
            for j, bars in enumerate([sf, sb, df, db, tf, tb]):
                if atom in bars[0]:
                    mean = bars[0][atom]
                    std = bars[1][atom]
                    ax.bar(i+barnumber*width, mean, width, yerr=std, align='edge', ecolor='black', capsize=5, color=colors[j])
                    ax.text(i+barnumber*width + width/2, 0.05, f"{mean:.2f} ± {std:.2f}", horizontalalignment='center',
                            transform=ax.get_xaxis_transform(), fontsize=13, color='black', rotation='vertical')
                    ax.hlines(monat_r[i], i+barnumber*width, i+barnumber*width + width, color='r', alpha=0.8, linestyles='dashed')
                    barnumber += 1

        # Plot empirical data as scatter points
        #ax.scatter(categories, monoatomic_radii, label='Empirical', color='red', s=140)  # Adjust scatter marker size here
        # Create the legend
        from matplotlib.patches import Patch
        legend_elements = []
        for col,label in zip(colors, labels):
            legend_elements.append(Patch(facecolor=col, edgecolor=col, label=label))
        ax.legend(handles=legend_elements, fontsize=13, loc='center right')

        # X Labels
        ax.set_xlabel('Atomic Number', fontsize=20)
        ax.set_xticks(np.arange(0,len(countbars)))
        ax.set_xticklabels([symbols[i] for i in sorted(countbars)])

        # Y Label
        ax.set_ylabel('Mean Radius (angstrom)', fontsize=20)

        # Grid Options
        ax.yaxis.grid(True)
        ax.set_axisbelow(True)

        # # Add legend for data counts per atom
        # ffod_datums = [f'{symbols[atom]}: {len(sf[atom])}' for atom in sf]
        # ffod_datums = '\n'.join(ffod_datums)
        # bfod_datums = [f'{symbols[atom]}: {len(sb[atom])}' for atom in sb]
        # bfod_datums = '\n'.join(bfod_datums)
        # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # ax.text(0.05, 0.95, bfod_datums, transform=ax.transAxes, fontsize=20,
        #     verticalalignment='top', bbox=props)

        # # Statistics for each bar
        # for i, (mean, std) in enumerate(zip(sf_u, sf_std)):
        #         ax.text(i - width/2, 0.15, f"{mean:.2f}\n± {std:.2f}", horizontalalignment='center', transform=ax.get_xaxis_transform(), fontsize=18, color='black')
        # for i, (mean, std) in enumerate(zip(sb_u, bf_std)):
        #         ax.text(i + width/2, 0.15, f"{mean:.2f}\n± {std:.2f}", horizontalalignment='center', transform=ax.get_xaxis_transform(), fontsize=18, color='black')

        return fig, ax

    # Sort the values based on the first
    sf, sb = molSet.tally_ffod_bfod_radii(SFFOD)
    sf = dict(sorted(sf.items(), key=lambda item: item[0]))
    sb = dict(sorted(sb.items(), key=lambda item: item[0]))

    # DFFOD
    df, db = molSet.tally_ffod_bfod_radii(DFFOD)
    df = dict(sorted(df.items(), key=lambda item: item[0]))
    db = dict(sorted(db.items(), key=lambda item: item[0]))

    # TFFOD
    tf, tb = molSet.tally_ffod_bfod_radii(TFFOD)
    tf = dict(sorted(tf.items(), key=lambda item: item[0]))
    tb = dict(sorted(tb.items(), key=lambda item: item[0]))

    colors = ['bisque', 'darkorange', 'skyblue', 'royalblue', 'limegreen', 'forestgreen']
    # fig1, ax1 = CreateBarChart(sf, sb, colors[1], colors[0])
    # fig2, ax2 = CreateBarChart(df, db, colors[3], colors[2])
    # fig3, ax3 = CreateBarChart(tf, tb, colors[5], colors[4])
    fig, ax = CreateBarChart2(sf, sb, df, db, tf, tb)
    ax.set_title('Radii Comparison of FFODs and BFODs', fontsize=20)

    # Show plot
    pp.tight_layout()
    fig.savefig(savepath + 'free_bfod_radii_pp.svg')
    pp.show()

def r_ratio(Z1, Z2):
    from numpy import sqrt

    if Z1 == Z2:
        return 0.5
    else:
        Z1 = 0.4 if Z1 == 1 else Z1
        Z2 = 0.4 if Z2 == 1 else Z2
        r = sqrt(Z1/Z2)
        g = r/(1+r)
        if g < 0.5:
            return g
        else:
            return (1-g)

def Intra_DBFOD_Dist(mols: MolecularSet):
    """
    Graphs a Heatmap of bond distance for different forms of bonds. The purpose is
    to identify the ways in which edge distance might correlate with the dominance of the bond. 
    """
    import numpy as np
    from MolecularSet import MolecularSet

    # Dictionary to store bond classes and their associated bond distance
    dbfod_dist = {}
    dbfod_mMeekR = {}
    dbfod_mBoldR = {}
    dbfod_mMeekAng = {}
    dbfod_mBoldAng = {}
    dbfod_port = {}
    # Intra-Distance Stats
    avg_bond_dist = {}
    std_dev_bond_dist = {}
    # BFOD BoldPortion
    avg_port = {}
    std_dev_port = {}
    # BFOD Distances Stats
    avg_MeekR = {}
    std_dev_MeekR = {}
    average_BoldR = {}
    std_dev_BoldR = {}
    # BFOD Angles Stats
    avg_MeekAng = {}
    std_dev_MeekAng = {}
    average_BoldAng = {}
    std_dev_BoldAng = {}

    # Iterate through each bfod. Use the MolecularSet Class for less clutter
    for bond in mols.GetBonds(2):
        # Create a unique key-value pair for each bonding pair
        key = frozenset({bond.mAtoms[0].mZ,bond.mAtoms[1].mZ})

        # Add key to dictionary if nonexistent.
        # TODO: Can probably add this to the <MolecularSet> class 
        if key not in dbfod_dist:
            dbfod_port[key] = []
            dbfod_dist[key] = []
            dbfod_mMeekR[key] = []
            dbfod_mBoldR[key] = []
            dbfod_mMeekAng[key] = []
            dbfod_mBoldAng[key] = []
        dbfod_port[key].append(bond.GetPort())
        dbfod_dist[key].append(bond.GetDist())
        dbfod_mMeekR[key].extend(bond.GetMeekR())
        dbfod_mBoldR[key].extend(bond.GetBoldR())
        dbfod_mMeekAng[key].extend(bond.GetMeekAng())
        dbfod_mBoldAng[key].extend(bond.GetBoldAng())

    # Calculate the average and std. deviation distances for each diatomic
    for bond_pair, distances in dbfod_dist.items():
        avg_bond_dist[bond_pair] = np.mean(distances)
        std_dev_bond_dist[bond_pair] = np.std(distances)
    for bond_pair, r in dbfod_mMeekR.items():
        avg_MeekR[bond_pair] = np.mean(r)
        std_dev_MeekR[bond_pair] = np.std(r)
    for bond_pair, r in dbfod_mBoldR.items():
        average_BoldR[bond_pair] = np.mean(r)
        std_dev_BoldR[bond_pair] = np.std(r)
    for bond_pair, r in dbfod_port.items():
        avg_port[bond_pair] = np.mean(r)
        std_dev_port[bond_pair] = np.std(r)
    for bond_pair, r in dbfod_mMeekAng.items():
        avg_MeekAng[bond_pair] = np.mean(r)
        std_dev_MeekAng[bond_pair] = np.std(r)
    for bond_pair, r in dbfod_mBoldAng.items():
        average_BoldAng[bond_pair] = np.mean(r)
        std_dev_BoldAng[bond_pair] = np.std(r)

    # Extract unique atomic numbers from bond pairs
    unique_atomic_numbers = set()
    for bond_pair in avg_port.keys():
        unique_atomic_numbers.update(bond_pair)

    # Create a sorted list of unique atomic numbers
    atomic_numbers_sorted = sorted(unique_atomic_numbers)

    # Create a matrix to store the average bond proportions
    matrix_size = len(atomic_numbers_sorted)
    bond_matrix = np.zeros((matrix_size, matrix_size))

    # Figure 2
    matrix2 = np.zeros((matrix_size, matrix_size))
    # Fill in the matrix with the average bond proportions
    for bond_pair, avg_dbfod_dist in avg_port.items():
        atomic_numbers_list = sorted(list(bond_pair))
        if len(atomic_numbers_list) == 1:
            atomic_number = atomic_numbers_list[0]
            index = atomic_numbers_sorted.index(atomic_number)
            matrix2[index, index] = avg_dbfod_dist
        else:
            atomic_number_1, atomic_number_2 = atomic_numbers_list
            index_1 = atomic_numbers_sorted.index(atomic_number_1)
            index_2 = atomic_numbers_sorted.index(atomic_number_2)
            matrix2[index_1, index_2] = avg_dbfod_dist
            matrix2[index_2, index_1] = avg_dbfod_dist

    fig2, ax2 = pp.subplots(figsize=(8,8))
    heatmap2 = ax2.imshow(matrix2, cmap='summer', interpolation='nearest', vmin=matrix2.min(), vmax=matrix2.max())

    for i in range(matrix_size):
        m = atomic_numbers_sorted[i]
        for j in range(matrix_size):
            n = atomic_numbers_sorted[j]
            avg_port_grid = avg_port.get(frozenset({m,n}), 0)
            s2_port_grid = std_dev_port.get(frozenset({m,n}), 0)
            
            # Distances
            s2_port_grid = std_dev_port.get(frozenset({m,n}), 0)
            num_dbfods = len(dbfod_dist.get(frozenset({m,n}), []))
            # MeekR
            avg_MeekAngle_grid = np.rad2deg(avg_MeekAng.get(frozenset({m,n}), 0))
            s2_MeekAngle_grid = np.rad2deg(std_dev_MeekAng.get(frozenset({m,n}), 0))
            # BoldAngle
            avg_BoldAngle_grid = np.rad2deg(average_BoldAng.get(frozenset({m,n}), 0))
            s2_BoldAngle_grid = np.rad2deg(std_dev_BoldAng.get(frozenset({m,n}), 0))
            # MeekR
            avg_MeekR_grid = avg_MeekR.get(frozenset({m,n}), -1)
            s2_MeekR_grid = std_dev_MeekR.get(frozenset({m,n}), 0)
            # BoldR
            avg_BoldR_grid = average_BoldR.get(frozenset({m,n}), -1)
            s2_BoldR_grid = std_dev_BoldR.get(frozenset({m,n}), 0)
            # Notes
            # note_title = f'{num_dbfods} Observations\n\n\n\n\n\n\n\n\n'
            # note_IntraDist = f'\n\n\nBold Proportion\n{avg_port_grid:.4f} ± {s2_port_grid:.4f}\n\n'
            # note_MeekAngle = f'MeekAngle\n{avg_MeekAngle_grid:.3f} ± {s2_MeekAngle_grid:.3f}\n\n'
            # note_BoldAngle = f'BoldAngle\n{avg_BoldAngle_grid:.3f} ± {s2_BoldAngle_grid:.3f}\n\n'
            # notation = note_IntraDist + note_MeekAngle + note_BoldAngle
            note_observs = f'{num_dbfods} Observations\n'
            note_port = f'Bold Portion\n{avg_port_grid:.3f} $\pm$ {s2_port_grid:.3f}\n'
            note_MeekR = f'MeekR\nμ = {avg_MeekR_grid:.2f} $\pm$ {s2_MeekR_grid:.2f}\n'
            note_BoldR = f'BoldR\nμ = {avg_BoldR_grid:.2f} $\pm$ {s2_BoldR_grid:.2f}'
            notation = note_observs + note_port + note_MeekR + note_BoldR
            ax2.annotate(notation, (j, i), color='black', ha='center', va='center',fontsize=7)

    # Labels
    pp.colorbar(heatmap2, ax=ax2, label='Bold Proportion on Bonding Axis', fraction=0.04, pad=0.04)
    pp.xticks(np.arange(matrix_size), atomic_numbers_sorted)
    pp.yticks(np.arange(matrix_size), atomic_numbers_sorted)
    pp.ylabel('Atomic Number')
    pp.xlabel('Atomic Number')
    pp.title('Average Bold Proportion Heatmap in DBFODs with Angle Statistics')
    
    #Size
    # fig2.subplots_adjust(left=0.1,
    #                 bottom=0.05,
    #                 right=1,
    #                 top=0.95,
    #                 wspace=0.4,
    #                 hspace=0.4)
    pp.tight_layout()
    fig2.savefig(savepath + '/heatmap2.svg')
    return fig2

def Intra_SBFOD_Dist(mols: MolecularSet):
    """
    Analyze SBFOD of whole list of molecules
    """
    from MolecularSet import MolecularSet

    # Dictionary to store bond classes and their associated bond distance
    sbfod_dist = {}
    sbfod_mMeekR = {}
    sbfod_mBoldR = {}
    sbfod_mMeekAng = {}
    sbfod_mBoldAng = {}
    sbfod_port = {}
    # BFOD BoldPortion
    avg_port = {}
    std_dev_port = {}
    # BFOD Distances Stats
    avg_MeekR = {}
    std_dev_MeekR = {}
    average_BoldR = {}
    std_dev_BoldR = {}
    # BFOD Angles Stats
    avg_MeekAng = {}
    std_dev_MeekAng = {}
    average_BoldAng = {}
    std_dev_BoldAng = {}

    # Iterate through each bfod. Use the MolecularSet Class for less clutter
    for bond in mols.GetBonds(1):
        # Create a unique key-value pair for each bonding pair
        key = frozenset({bond.mAtoms[0].mZ,bond.mAtoms[1].mZ})

        # Add key to dictionary if nonexistent.
        # TODO: Can probably add this to the <MolecularSet> class
        if key not in sbfod_dist:
            sbfod_port[key] = []
            sbfod_dist[key] = []
            sbfod_mMeekR[key] = []
            sbfod_mBoldR[key] = []
            sbfod_mMeekAng[key] = []
            sbfod_mBoldAng[key] = []
        sbfod_port[key].append(bond.GetPort())
        sbfod_dist[key].append(bond.GetDist())
        sbfod_mMeekR[key].extend(bond.GetMeekR())
        sbfod_mBoldR[key].extend(bond.GetBoldR())
        sbfod_mMeekAng[key].extend(bond.GetMeekAng())
        sbfod_mBoldAng[key].extend(bond.GetBoldAng())

    # Calculate the average and std. deviation distances for each diatomic
    for bond_pair, distances in sbfod_dist.items():
        avg_port[bond_pair] = np.mean(distances)
        std_dev_port[bond_pair] = np.std(distances)
    for bond_pair, r in sbfod_mMeekR.items():
        avg_MeekR[bond_pair] = np.mean(r)
        std_dev_MeekR[bond_pair] = np.std(r)
    for bond_pair, r in sbfod_mBoldR.items():
        average_BoldR[bond_pair] = np.mean(r)
        std_dev_BoldR[bond_pair] = np.std(r)
    for bond_pair, r in sbfod_port.items():
        avg_port[bond_pair] = np.mean(r)
        std_dev_port[bond_pair] = np.std(r)
    for bond_pair, r in sbfod_mMeekAng.items():
        avg_MeekAng[bond_pair] = np.mean(r)
        std_dev_MeekAng[bond_pair] = np.std(r)
    for bond_pair, r in sbfod_mBoldAng.items():
        average_BoldAng[bond_pair] = np.mean(r)
        std_dev_BoldAng[bond_pair] = np.std(r)

    # Extract unique atomic numbers from bond pairs
    unique_atomic_numbers = set()
    for bond_pair in avg_port.keys():
        unique_atomic_numbers.update(bond_pair)

    # Create a sorted list of unique atomic numbers
    atomic_numbers_sorted = sorted(unique_atomic_numbers)

    # Create a matrix to store the average bond proportions
    matrix_size = len(atomic_numbers_sorted)
    bond_matrix = np.zeros((matrix_size, matrix_size))

    # Fill in the matrix with the average bond proportions
    for bond_pair, avg_sbfod_dist in avg_port.items():
        atomic_numbers_list = sorted(list(bond_pair))
        if len(atomic_numbers_list) == 1:
            atomic_number = atomic_numbers_list[0]
            index = atomic_numbers_sorted.index(atomic_number)
            bond_matrix[index, index] = avg_sbfod_dist
        else:
            atomic_number_1, atomic_number_2 = atomic_numbers_list
            index_1 = atomic_numbers_sorted.index(atomic_number_1)
            index_2 = atomic_numbers_sorted.index(atomic_number_2)
            bond_matrix[index_1, index_2] = avg_sbfod_dist
            bond_matrix[index_2, index_1] = avg_sbfod_dist

    # Plot the heatmap with modified colormap and scaled data
    fig, ax = pp.subplots(figsize=(9,8))
    heatmap = ax.imshow(bond_matrix, cmap='summer', interpolation='nearest', vmin=bond_matrix.min(), vmax=bond_matrix.max())

    # Fig1: Annotate the averages and standard deviations in the grid boxes
    for i in range(matrix_size):
        m = atomic_numbers_sorted[i]
        for j in range(matrix_size):
            n = atomic_numbers_sorted[j]
            # Distances
            num_sbfods = len(sbfod_dist.get(frozenset({m,n}), []))
            # MeekR
            avg_MeekR_grid = avg_MeekR.get(frozenset({m,n}), -1)
            s2_MeekR_grid = std_dev_MeekR.get(frozenset({m,n}), 0)
            # BoldR
            avg_BoldR_grid = average_BoldR.get(frozenset({m,n}), -1)
            s2_BoldR_grid = std_dev_BoldR.get(frozenset({m,n}), 0)
            # Portions
            avg_port_grid = avg_port.get(frozenset({m,n}), 0)
            s2_port_grid = std_dev_port.get(frozenset({m,n}), 0)
            # Notes
            note_observs = f'{num_sbfods} Obs.\n'
            # note_port = f'Bold Portion\n{avg_port_grid:.3f} $\pm$ {s2_port_grid:.3f}\n'
            note_port = f'r = {avg_port_grid:.2f}$\pm${s2_port_grid:.2f}\n'
            note_MeekR = f'$\\overline{{m}}$ = {avg_MeekR_grid:.2f}$\pm${s2_MeekR_grid:.2f}\n'
            note_BoldR = f'$\\overline{{b}}$ = {avg_BoldR_grid:.2f}$\pm${s2_BoldR_grid:.2f}'
            annotation = note_observs + note_port + note_MeekR + note_BoldR
            if bond_matrix[i, j] != 0:
                ax.annotate(annotation, (j, i), color='black', ha='center', va='center', size=8.5)
            ax.axvline(0.5, c='black')
            ax.axvline(1.5, c='black')
            ax.axvline(2.5, c='black')
            ax.axvline(3.5, c='black')
            ax.axvline(4.5, c='black')
            ax.axhline(0.5, c='black')
            ax.axhline(1.5, c='black')
            ax.axhline(2.5, c='black')
            ax.axhline(3.5, c='black')
            ax.axhline(4.5, c='black')

    # Labels
    pp.colorbar(heatmap, ax=ax, label='Axial ratio of $b$', fraction=0.03, pad=0.04, shrink=0.6)
    pp.xticks(np.arange(matrix_size), [symbols[i] for i in atomic_numbers_sorted])
    pp.yticks(np.arange(matrix_size), [symbols[i] for i in atomic_numbers_sorted])
    pp.xlabel('Atom')
    pp.ylabel('Atom')
    pp.title('SBFOD Bold Portions and Distances')

    # Labels
    pp.tight_layout()
    fig.savefig(savepath + 'heatmap.png')
    return fig

def Intra_TBFOD_Dist(mols: MolecularSet):
    """
    Graphs a Heatmap of bond distance for different forms of bonds. The purpose is
    to identify the ways in which edge distance might correlate with the dominance of the bond.
    """
    import numpy as np
    from MolecularSet import MolecularSet

    # Dictionary to store bond classes and their associated bond distance
    tbfod_dist = {}
    tbfod_mMeekR = {}
    tbfod_mBoldR = {}
    tbfod_mMeekAng = {}
    tbfod_mBoldAng = {}
    tbfod_port = {}
    # Intra-Distance Stats
    avg_bond_dist = {}
    std_dev_bond_dist = {}
    # BFOD BoldPortion
    avg_port = {}
    std_dev_port = {}
    # BFOD Distances Stats
    avg_MeekR = {}
    std_dev_MeekR = {}
    average_BoldR = {}
    std_dev_BoldR = {}
    # BFOD Angles Stats
    avg_MeekAng = {}
    std_dev_MeekAng = {}
    average_BoldAng = {}
    std_dev_BoldAng = {}

    # Iterate through each bfod. Use the MolecularSet Class for less clutter
    for bond in mols.GetBonds(3):
        # Create a unique key-value pair for each bonding pair
        key = frozenset({bond.mAtoms[0].mZ,bond.mAtoms[1].mZ})

        # Add key to dictionary if nonexistent.
        # TODO: Can probably add this to the <MolecularSet> class
        if key not in tbfod_dist:
            tbfod_port[key] = []
            tbfod_dist[key] = []
            tbfod_mMeekR[key] = []
            tbfod_mBoldR[key] = []
            tbfod_mMeekAng[key] = []
            tbfod_mBoldAng[key] = []
        tbfod_port[key].append(bond.GetPort())
        tbfod_dist[key].append(bond.GetDist())
        tbfod_mMeekR[key].extend(bond.GetMeekR())
        tbfod_mBoldR[key].extend(bond.GetBoldR())
        tbfod_mMeekAng[key].extend(bond.GetMeekAng())
        tbfod_mBoldAng[key].extend(bond.GetBoldAng())

    # Calculate the average and std. deviation distances for each diatomic
    for bond_pair, distances in tbfod_dist.items():
        avg_bond_dist[bond_pair] = np.mean(distances)
        std_dev_bond_dist[bond_pair] = np.std(distances)
    for bond_pair, r in tbfod_mMeekR.items():
        avg_MeekR[bond_pair] = np.mean(r)
        std_dev_MeekR[bond_pair] = np.std(r)
    for bond_pair, r in tbfod_mBoldR.items():
        average_BoldR[bond_pair] = np.mean(r)
        std_dev_BoldR[bond_pair] = np.std(r)
    for bond_pair, r in tbfod_port.items():
        avg_port[bond_pair] = np.mean(r)
        std_dev_port[bond_pair] = np.std(r)
    for bond_pair, r in tbfod_mMeekAng.items():
        avg_MeekAng[bond_pair] = np.mean(r)
        std_dev_MeekAng[bond_pair] = np.std(r)
    for bond_pair, r in tbfod_mBoldAng.items():
        average_BoldAng[bond_pair] = np.mean(r)
        std_dev_BoldAng[bond_pair] = np.std(r)

    # Extract unique atomic numbers from bond pairs
    unique_atomic_numbers = set()
    for bond_pair in avg_port.keys():
        unique_atomic_numbers.update(bond_pair)

    # Create a sorted list of unique atomic numbers
    atomic_numbers_sorted = sorted(unique_atomic_numbers)

    # Create a matrix to store the average bond proportions
    matrix_size = len(atomic_numbers_sorted)
    bond_matrix = np.zeros((matrix_size, matrix_size))

    # Figure 2
    matrix2 = np.zeros((matrix_size, matrix_size))
    # Fill in the matrix with the average bond proportions
    for bond_pair, avg_tbfod_dist in avg_port.items():
        atomic_numbers_list = sorted(list(bond_pair))
        if len(atomic_numbers_list) == 1:
            atomic_number = atomic_numbers_list[0]
            index = atomic_numbers_sorted.index(atomic_number)
            matrix2[index, index] = avg_tbfod_dist
        else:
            atomic_number_1, atomic_number_2 = atomic_numbers_list
            index_1 = atomic_numbers_sorted.index(atomic_number_1)
            index_2 = atomic_numbers_sorted.index(atomic_number_2)
            matrix2[index_1, index_2] = avg_tbfod_dist
            matrix2[index_2, index_1] = avg_tbfod_dist

    fig, ax = pp.subplots(figsize=(8,8))
    heatmap2 = ax.imshow(matrix2, cmap='summer', interpolation='nearest', vmin=matrix2.min(), vmax=matrix2.max())

    for i in range(matrix_size):
        m = atomic_numbers_sorted[i]
        for j in range(matrix_size):
            n = atomic_numbers_sorted[j]
            avg_port_grid = avg_port.get(frozenset({m,n}), 0)
            s2_port_grid = std_dev_port.get(frozenset({m,n}), 0)

            # Distances
            s2_port_grid = std_dev_port.get(frozenset({m,n}), 0)
            num_tbfods = len(tbfod_dist.get(frozenset({m,n}), []))
            # MeekR
            avg_MeekAngle_grid = np.rad2deg(avg_MeekAng.get(frozenset({m,n}), 0))
            s2_MeekAngle_grid = np.rad2deg(std_dev_MeekAng.get(frozenset({m,n}), 0))
            # BoldAngle
            avg_BoldAngle_grid = np.rad2deg(average_BoldAng.get(frozenset({m,n}), 0))
            s2_BoldAngle_grid = np.rad2deg(std_dev_BoldAng.get(frozenset({m,n}), 0))
            # Notes
            note_title = f'{num_tbfods} Observations\n\n\n\n\n\n\n\n\n'
            note_IntraDist = f'\n\n\nBold Proportion\n{avg_port_grid:.4f} ± {s2_port_grid:.4f}\n\n'
            note_MeekAngle = f'MeekAngle\n{avg_MeekAngle_grid:.3f} ± {s2_MeekAngle_grid:.3f}\n\n'
            note_BoldAngle = f'BoldAngle\n{avg_BoldAngle_grid:.3f} ± {s2_BoldAngle_grid:.3f}\n\n'
            notation = note_IntraDist + note_MeekAngle + note_BoldAngle
            ax.annotate(note_title, (j, i), color='black', ha='center', va='center', weight='semibold',fontsize=7)
            ax.annotate(notation, (j, i), color='black', ha='center', va='center',fontsize=7)

    # Labels
    pp.colorbar(heatmap2, ax=ax, label='Bold Proportion on Bonding Axis', fraction=0.03, pad=0.04, )
    pp.xticks(np.arange(matrix_size), atomic_numbers_sorted)
    pp.yticks(np.arange(matrix_size), atomic_numbers_sorted)
    pp.ylabel('Atomic Number')
    pp.xlabel('Atomic Number')
    pp.title('Average Bold Proportion Heatmap in TBFODs with Angle Statistics')

    #Size
    # fig2.subplots_adjust(left=0.1,
    #                 bottom=0.05,
    #                 right=1,
    #                 top=0.95,
    #                 wspace=0.4,
    #                 hspace=0.4)
    pp.tight_layout()
    #fig2.savefig(savepath + '/heatmap3.svg')
    #pp.show()
    return fig

def tables(mols):
    fig1 = Intra_SBFOD_Dist(mols)
    fig2 = Intra_DBFOD_Dist(mols)
    fig3 = Intra_TBFOD_Dist(mols)
