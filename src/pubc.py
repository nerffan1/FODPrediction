#!/usr/bin/python3
import pubchempy as pcp

def create_xyz_file(compound, filename):
    # Retrieve compound information
    compound_info = pcp.get_compounds(compound, 'name', record_type='3d')[0]

    # Write XYZ file
    with open(filename, 'w') as f:
        f.write(str(len(compound_info.atoms)) + '\n')
        f.write('\n')

        for atom in compound_info.atoms:
            f.write(f"{atom.element} {atom.x} {atom.y} {atom.z}\n")

    print(f"XYZ file '{filename}' created successfully.")

# Example usage
if __name__ == "__main__":
    compound_name = input("Enter the name of the compound: ")
    filename = input("Enter the filename for the XYZ file: ")

    create_xyz_file(compound_name, filename)
