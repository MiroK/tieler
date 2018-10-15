from tieler.periodicity import compute_vertex_periodicity, compute_entity_periodicity
from dolfin import UnitCubeMesh, UnitSquareMesh, CompiledSubDomain
import numpy as np
import unittest


class TestPeriodicity(unittest.TestCase):
    
    def test_vertex_2d(self):
        mesh = UnitSquareMesh(8, 8)

        x = mesh.coordinates()
        min_ = np.min(x, axis=0)
        max_ = np.max(x, axis=0)

        tol = 1E-9  # The precision in gmsh isn't great, probably why DOLFIN's
        # periodic boundary computation is not working

        # Check x periodicity
        master = CompiledSubDomain('near(x[0], A, tol)', A=min_[0], tol=tol)
        slave = CompiledSubDomain('near(x[0], A, tol)', A=max_[0], tol=tol)

        shift_x = np.array([max_[0]-min_[0], 0])
        to_master = lambda x, shift=shift_x: x - shift

        error, mapping = compute_vertex_periodicity(mesh, master, slave, to_master)
        self.assertTrue(error < 10*tol)

        _, mapping = compute_entity_periodicity(1, mesh, master, slave, to_master)
        self.assertTrue(len(list(mapping.keys())) == 8)
        

    def test_vertex_3d(self):
        mesh = UnitCubeMesh(8, 8, 8)

        x = mesh.coordinates()
        min_ = np.min(x, axis=0)
        max_ = np.max(x, axis=0)

        tol = 1E-9  # The precision in gmsh isn't great, probably why DOLFIN's
        # periodic boundary computation is not working

        # Check x periodicity
        master = CompiledSubDomain('near(x[1], A, tol)', A=min_[1], tol=tol)
        slave = CompiledSubDomain('near(x[1], A, tol)', A=max_[1], tol=tol)

        shift_x = np.array([0, max_[1]-min_[1], 0])
        to_master = lambda x, shift=shift_x: x - shift

        error, mapping = compute_vertex_periodicity(mesh, master, slave, to_master)
        self.assertTrue(error < 10*tol)

        _, mapping = compute_entity_periodicity(2, mesh, master, slave, to_master)
        self.assertTrue(len(list(mapping.keys())) == 128)
