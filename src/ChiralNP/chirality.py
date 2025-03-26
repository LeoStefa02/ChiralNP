from .nanoparticle import *
from scipy.spatial import distance
import tensorflow as tf

class ChiralityChecker:
    def check_chirality(self, nanoparticle, method='hausdorff'):
        if method == 'hausdorff':
            # 1. Mirror image
            mirrored_positions = nanoparticle.positions.copy()
            mirrored_positions[:, 0] *= -1  # Invert x-coordinates

            R, t, rmsd = self._kabsch(nanoparticle.positions, mirrored_positions)

            # 2.d Apply Rotation
            rotated_mirrored_positions = np.dot(R, mirrored_positions.T)

            # 3. Calculate Hausdorff distance
            hausdorff_distance = self._calculate_hausdorff_distance(
                nanoparticle.positions, rotated_mirrored_positions.T
            )
            
            # 4. Normalize by the diameter of the nanoparticle
            diameter = nanoparticle.get_diameter()
            chirality_measure = hausdorff_distance / diameter
            
            return chirality_measure
        else:
            raise ValueError(f"Method '{method}' not supported. Use 'hausdorff'.")
    
    def _kabsch(self, P, Q):
        assert P.shape == Q.shape, "Matrix dimensions must match"

        # Compute centroids
        centroid_P = tf.reduce_mean(P, axis=0)
        centroid_Q = tf.reduce_mean(Q, axis=0)

        # Optimal translation
        t = centroid_Q - centroid_P

        # Center the points
        p = P - centroid_P
        q = Q - centroid_Q

        # Compute the covariance matrix
        H = tf.matmul(tf.transpose(p), q)

        # SVD
        S, U, V = tf.linalg.svd(H)

        # Validate right-handed coordinate system
        d = tf.linalg.det(tf.matmul(V, tf.transpose(U)))

        # Create two possible rotation matrices and choose the right one
        R_positive = tf.matmul(V, tf.transpose(U))


        # Create a modified V with the last row flipped
        V_modified = tf.concat([V[:-1], -V[-1:]], axis=0)
        R_negative = tf.matmul(V_modified, tf.transpose(U))

        # Choose the right rotation matrix based on determinant
        R = tf.cond(d < 0.0, lambda: R_negative, lambda: R_positive)

        # RMSD
        rmsd = tf.sqrt(tf.reduce_sum(tf.square(tf.matmul(p, tf.transpose(R)) - q)) / P.shape[0])

        return R, t, rmsd
    
    def _calculate_hausdorff_distance(self, set_a, set_b):
        # Calculate distances from each point in set_a to the closest point in set_b
        distances_a_to_b = np.array([
            np.min(distance.cdist([point_a], set_b)) for point_a in set_a
        ])
        
        # Calculate distances from each point in set_b to the closest point in set_a
        distances_b_to_a = np.array([
            np.min(distance.cdist([point_b], set_a)) for point_b in set_b
        ])
        
        # The Hausdorff distance is the maximum of these minimum distances
        forward_hausdorff = np.max(distances_a_to_b)
        backward_hausdorff = np.max(distances_b_to_a)
        
        return max(forward_hausdorff, backward_hausdorff)


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