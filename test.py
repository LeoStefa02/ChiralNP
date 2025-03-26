from ChiralNP import nanoparticle, shapes, crystal, chirality

# Example Usage
# Create an XYZ nanoparticle
xyz_nanoparticle = nanoparticle.Nanoparticle('../seeds/Cu4631.xyz')

sphere = shapes.Sphere(radius=50)                      # Sphere shape
cube = shapes.Cube(edge=50)                            # Cube Shape
squirrel = shapes.Stl('../seeds/stl/squirrel.stl', 50) # Squirrel Shape

# Create FCC lattice
fcc_lattice = crystal.FCCLattice(a=4.08)

# Generate FCC points within the sphere
fcc_points_sphere = fcc_lattice.generate_points(sphere)
fcc_points_cube = fcc_lattice.generate_points(cube)
fcc_points_squirrel = fcc_lattice.generate_points(squirrel)


nanoparticle = nanoparticle.Nanoparticle()
# Create a Nanoparticle from the FCC points
# (You'd need to determine atom types based on the desired material)
# For example for gold:
atom_types = ['Au'] * len(fcc_points_sphere)
nanoparticle.set_atoms(fcc_points_sphere, atom_types)
nanoparticle.write_xyz('../seeds/gold_sphere.xyz')

atom_types = ['Au'] * len(fcc_points_squirrel)
nanoparticle.set_atoms(fcc_points_squirrel, atom_types)
nanoparticle.write_xyz('../seeds/gold_squirrel.xyz')

# And for Cu
atom_types = ['Cu'] * len(fcc_points_cube)
nanoparticle.set_atoms(fcc_points_cube, atom_types)
nanoparticle.write_xyz('../seeds/Cu_cube.xyz')

# Check chirality
# checker = ChiralityChecker()
# is_chiral = checker.check_chirality(gold_nanoparticle)
# print(f'Is initially chiral: {is_chiral}')

# Induce chirality (example using removal)
rules = [
    {'type': 'plane', 'normal': [0, 0, 1], 'threshold': 1.5},   # Remove near plane
    {'type': 'cylinder', 'axis': [0, 0, 1], 'radius': 3.0},     # Remove inside cylinder
    {'type': 'sphere', 'center': [0, 0, 0], 'radius': 30.0}      # Remove inside sphere
]
modifier = chirality.ChiralityModifier()
modifier.modify_chirality(nanoparticle, rules)

# Check chirality again
# is_chiral_modified = checker.check_chirality(modified_nanoparticle)
# print(f"Is chiral after modification: {is_chiral_modified}")

# Save the modified nanoparticle
nanoparticle.write_xyz('../seeds/modified.xyz')