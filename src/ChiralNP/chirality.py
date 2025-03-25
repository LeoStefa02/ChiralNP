from nanoparticle import *

class ChiralityChecker:
    def check_chirality(self, nanoparticle, method='hausdorff'):
        # Implement different chirality checking algorithms
        pass

class ChiralityModifier:
    def modify_chirality(self, nanoparticle, rules):
        for rule in rules:
            rule_type = rule['type']
            if rule_type == 'plane':
                self.remove_near_plane(
                    nanoparticle, rule['normal'], rule['threshold']
                )
            elif rule_type == 'cylinder':
                self.remove_inside_cylinder(
                    nanoparticle, rule['axis'], rule['radius']
                )
            elif rule_type == 'sphere':
                self.remove_inside_sphere(
                    nanoparticle, rule['center'], rule['radius']
                )
            else:
                raise ValueError(f'Unknow rule type: {rule_type}')
        
        if not nanoparticle.positions.any():
            print('Warning: All atoms has been removed.')
                
    def remove_near_plane(self, nanoparticle, normal, threshold):
        new_positions = []
        new_atom_types = []
        normal = np.array(normal) / np.linalg.norm(normal) # Normalize
        for i, pos in enumerate(nanoparticle.positions):
            distance = np.dot(pos, normal) # Simplified for plane through
            if abs(distance) > threshold:
                new_positions.append(pos)
                new_atom_types.append(nanoparticle.atom_types[i])
        nanoparticle.set_atoms(new_positions, new_atom_types)
    
    def remove_inside_cylinder(self, nanoparticle, axis, radius):
        new_positions = []
        new_atom_types = []
        axis = np.array(axis) / np.linalg.norm(axis) # Normalize
        for i, pos in enumerate(nanoparticle.positions):
            proj_onto_axis = np.dot(pos, axis) * axis
            dist_to_axis = np.linalg.norm(pos - proj_onto_axis)
            if dist_to_axis > radius: #Outside
                new_positions.append(pos)
                new_atom_types.append(nanoparticle.atom_types[i])
        nanoparticle.set_atoms(new_positions, new_atom_types)
    
    def remove_inside_sphere(self, nanoparticle, center, radius):
        new_positions = []
        new_atom_types = []
        center = np.array(center)
        for i, pos in enumerate(nanoparticle.positions):
            dist_to_center = np.linalg.norm(pos - center)
            if dist_to_center > radius: #Outside
                new_positions.append(pos)
                new_atom_types.append(nanoparticle.atom_types[i])
        nanoparticle.set_atoms(new_positions, new_atom_types)