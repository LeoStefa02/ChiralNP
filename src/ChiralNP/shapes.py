import numpy as np
from point_in_mesh import *

class Shape:
    def is_inside(self, point):
        # Check if a point is inside the shape
        raise NotImplementedError
    
    def boundary_box(self):
        # Returns dimension of the box containing the shape as [xEdge, yEdge, zEdge]
        raise NotImplementedError

class Sphere(Shape):
    def __init__(self, radius):
        self.radius = radius

    def is_inside(self, point):
        return np.linalg.norm(point) <= self.radius
    
    def boundary_box(self):
        minimum = int(np.floor(-self.radius))
        maximum = int(np.ceil(self.radius))
        m = max(-minimum, maximum) * 2
        # Box containing the sphere
        return [m, m, m]

class Cube(Shape):
    def __init__(self, edge):
        self.edge = edge

    def is_inside(self, point):
        x, y, z = point
        half_edge = self.edge / 2
        if (((x < half_edge) and (x > - half_edge))
            and ((y < half_edge) and (y > - half_edge))
            and ((z < half_edge) and (z > - half_edge))):
            return True
        else:
            return False
    
    def boundary_box(self):
        return [self.edge, self.edge, self.edge]
        
class Stl(Shape):
    def __init__(self, filename, size):
        self.mesh = trimesh.load(filename)
        self.mesh = self.mesh.apply_scale((1.0/max(self.mesh.extents))*size)
        bounds = self.mesh.bounds
        bbox_center = (bounds[0] + bounds[1]) / 2
        self.mesh.apply_translation(-bbox_center)

        self.fast_mesh = FastPointInMesh()
        vertices = np.array(self.mesh.vertices, dtype=np.float64)
        faces = np.array(self.mesh.faces, dtype=np.int32)
        self.fast_mesh.set_mesh(vertices, faces)
    
    def is_inside(self, point):
        return self.fast_mesh.contains_single_point(point)
    
    def boundary_box(self):
        bounds = self.mesh.bounds
        edge_lengths = bounds[1] - bounds[0]
        return edge_lengths
