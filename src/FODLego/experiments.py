import numpy as np
from scipy.spatial.distance import cdist


def CreateXYZ_Shells(sh1, sh2) -> None:
    """
    Create an XYZ file with
    """
    with open("shells.xyz",'w') as output:
        #First 2 lines
        output.write(str(len(sh1) + len(sh2)) + '\n')
        output.write("These are 2 shells\n")
        # Write all atoms
        for fod in sh1:
            atom_coords = ' '.join([f"{x:7.4f}" for x in fod])
            output.write(f"X {atom_coords}\n")

        for fod in sh2:
            atom_coords = ' '.join([f"{x:7.4f}" for x in fod])
            output.write(f"He {atom_coords}\n")

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
    sp3 = coords[sorted_i[:4]]
    sp3d5 = coords[sorted_i[4:]]

    #Return both shells
    return sp3, sp3d5

def rotate(coords, angle, axis):
    rotation_matrix = np.eye(3)
    if axis == 'x':
        rotation_matrix[1, 1] = np.cos(angle)
        rotation_matrix[1, 2] = -np.sin(angle)
        rotation_matrix[2, 1] = np.sin(angle)
        rotation_matrix[2, 2] = np.cos(angle)
    elif axis == 'y':
        rotation_matrix[0, 0] = np.cos(angle)
        rotation_matrix[0, 2] = np.sin(angle)
        rotation_matrix[2, 0] = -np.sin(angle)
        rotation_matrix[2, 2] = np.cos(angle)
    elif axis == 'z':
        rotation_matrix[0, 0] = np.cos(angle)
        rotation_matrix[0, 1] = -np.sin(angle)
        rotation_matrix[1, 0] = np.sin(angle)
        rotation_matrix[1, 1] = np.cos(angle)
    return np.dot(coords, rotation_matrix)

def electrostatic_potential_energy(coords1, coords2):
    # CDIST gets the pairwise distances, and we want the inverse of all those values
    inverse_values = 1/cdist(coords1, coords2)
    # We flatten the array and get the sum of all of these contributions, ignoring the
    # factors/constants
    return np.sum(inverse_values)

def find_minimum_energy_configuration(coords1, coords2):
    min_energy = float('inf')
    best_coords = None
    angles = np.linspace(0, 2 * np.pi, 500)
    for angle in angles:
        x_rot = rotate(coords2, angle, 'x')
        for angle2 in angles:
            xy_rot = rotate(x_rot, angle2, 'y')
            for angle3 in angles:
                rot_coor = rotate(xy_rot, angle3, 'z')
                energy = electrostatic_potential_energy(coords1, rot_coor)
                if energy < min_energy:
                    print(angle,angle2,angle3)
                    min_energy = energy
                    best_coords = rot_coor
    return best_coords, min_energy

sp3, sp3d5 = readXYZ("test2shellsxyz.xyz")

best_coords, min_energy = find_minimum_energy_configuration(sp3, sp3d5)
CreateXYZ_Shells(sp3, best_coords)
