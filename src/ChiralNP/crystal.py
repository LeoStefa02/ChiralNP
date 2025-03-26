import numpy as np
from math import pi

class Lattice:
    def __init__(self, a, b, c, alpha, beta, gamma):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.lattice_vectors = self._calculate_lattice_vectors()

    def _calculate_lattice_vectors(self):
        # Calculate lattice vectors from lattice parameters
        # Chosen a1 direction along x axis
        a1 = np.array([self.a, 0, 0])
        a2 = np.array([self.b * np.cos(self.gamma), self.b * np.sin(self.gamma), 0])
        a3 = np.array([
            self.c * np.cos(self.beta),
            self.c * (np.cos(self.alpha) - np.cos(self.b) * np.cos(self.gamma))**0.5,
            self.c * (1 - np.cos(self.beta)**2 - (np.cos(self.alpha) - np.cos(self.beta) * np.cos(self.gamma))**2)**0.5
        ])

    def generate_points(self, shape):
        # Generate lattice points within the given shape
        pass

class FCCLattice(Lattice):
    def __init__(self, a):
        super().__init__(a, a, a, pi/2, pi/2, pi/2)
        # FCC has cubic symmetry
    def generate_points(self, shape):
        # Generate FCC lattice points within the shape
        points = []
        
        # Basis vectors for the four atoms in the FCC unit cell.
        basis = self.a * np.array([[0.0, 0.0, 0.0],
                        [0.5, 0.5, 0.0],
                        [0.5, 0.0, 0.5],
                        [0.0, 0.5, 0.5]])
        
        edgeX, edgeY, edgeZ = shape.boundary_box()
        nx = int(np.floor((edgeX / 2) / self.a)) + 1
        ny = int(np.floor((edgeY / 2) / self.a)) + 1
        nz = int(np.floor((edgeZ / 2) / self.a)) + 1

        initial_points = []

        for i in range(-nx, nx + 1, 1):
            for j in range(-ny, ny + 1, 1):
                for k in range(-nz, nz + 1, 1):
                    # Translate the basis to each unit cell.
                    for b in basis:
                        point = b + np.array([i, j, k]) * self.a
                        initial_points.append(point)
        
        try:
            points = np.array(initial_points)[shape.is_inside_multi(initial_points)]
        except:
            for point in initial_points:
                if shape.is_inside(point):
                    points.append(point)

        return np.array(points)