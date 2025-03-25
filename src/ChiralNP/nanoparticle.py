import numpy as np
from crystal import FCCLattice
from shapes import Cube
from point_in_mesh import *

class Nanoparticle:
    def __init__(self, xyz_file=None):
        self.positions = None
        self.atom_types = None

        if xyz_file is not None:
            self.load_atoms(xyz_file)
    
    def set_atoms(self, positions, atom_types):
        self.positions = np.array(positions) # Nx3 array of coordinates
        self.atom_types = np.array(atom_types) # N array of atom types
    
    def load_atoms(self, xyz_file):
        try:
            # Use trimesh to load the mesh (supports many formats)
            positions = []
            atom_types = []
            with open(xyz_file, 'r') as f:
                natoms = int(f.readline())
                comment = f.readline()
                for _ in range(natoms):
                    line = f.readline().split()
                    atom_types.append(line[0])
                    positions.append([float(x) for x in line[1:]])
        except Exception as e:
            print(f"Error loading xyz: {e}")

    def translate(self, vector):
        self.positions += np.array(vector)

    def rotate(self, axis, angle):
        # Implement rotation using a rotation matrix
        pass

    def get_center_of_mass(self):
        # Calculate center of mass
        pass

    def write_xyz(self, filename):
        with open(filename, 'w') as f:
            # Number of atoms
            f.write(f"{len(self.positions)}\n")
        
            # Comment line
            f.write('Nanoparticle\n')
            
            # Atom positions
            for element, position in zip(self.atom_types, self.positions):
                x, y, z = position
                f.write(f"{element} {x:.6f} {y:.6f} {z:.6f}\n")

# class ChiralNanoparticle(Nanoparticle):
#     def is_chiral(self, method='hausdorff'):
#         # Implement chirality check using ChiralityChecker
#         if method == 'hausdorff':
#             pass #Use a Hausdorff distance based check
#         elif method == 'index':
#             pass #Use a chirality index
#         else:
#             raise ValueError("Invalid chirality check method")

#     def induce_chirality(self, method='removal', **kwargs):
#         # Implement chirality modification using ChiralityModifier
#         modifier = chirality.ChiralityModifier()
#         return modifier.modify_chirality(self, method=method, **kwargs)
    