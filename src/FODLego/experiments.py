#!/usr/bin/env python3

import numpy as np



def readXYZ(input_file):
    """
    Create a relaxed position in which the FODs have a minum of energy,
    based off an electrostatic potential.
    """
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        atoms = int(lines[0])
        coords = np.empty((atoms,3))

        # Process the remaining lines
        for i, line in enumerate(lines[2:]):
            parts = line.split()
            strnum = np.array(parts[1:4])  # First three are the coordinates
            coords[i] = strnum.astype(float)

    # Create the shells
    magnitude = np.linalg.norm(coords, axis=1)
    sorted_i = np.argsort(magnitude)
    print(coords[sorted_i])
    print (np.linalg.norm(coords[sorted_i], axis=1))



readXYZ("test2shellsxyz.xyz")
    # Load the shells
    #Shell1, Shell2 = 0,2
