import unittest
from tieler.tile import powers2, TileMesh
from tieler.tile_data import mf_from_data, load_data
from dolfin import (UnitSquareMesh, MeshFunction, CompiledSubDomain,
                    SubsetIterator)
import numpy as np


class TestTile(unittest.TestCase):

    def test_powers(self):
        check = all(sum(2**p for p in powers2(n)) == n for n in (9, 13, 425, 123))

        self.assertTrue(check)


    def test_tiling(self):

        tile = UnitSquareMesh(2, 2)

        mf = MeshFunction('size_t', tile, tile.topology().dim()-1, 0)
        CompiledSubDomain('near(x[0], 0.5) || near(x[1], 0.5)').mark(mf, 1)

        mesh_data = {}
        mesh_data = load_data(tile, mf, dim=1, data=mesh_data)

        mesh, mesh_data = TileMesh(tile, shape=(23, 13), mesh_data=mesh_data)
        f = mf_from_data(mesh, mesh_data)[1]

        self.assertTrue(np.linalg.norm(
            mesh.coordinates().min(axis=0) - np.zeros(2)) < 1E-13)
        
        self.assertTrue(np.linalg.norm(
            mesh.coordinates().max(axis=0) - np.array([23., 13.])) < 1E-13)
        
        self.assertTrue(23*13*4 == sum(1 for _ in SubsetIterator(f, 1)))
