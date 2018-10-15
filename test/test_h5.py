from tieler import TileMesh, mf_from_data, load_data
from mpi4py import MPI as pyMPI
import numpy as np
from dolfin import *

import unittest

class TestH5(unittest.TestCase):
    '''Check that the tiling works if the tile comes from H5 file'''
    # We have 2 physical volumes (0, 1)
    # And one physical surface tag 1
    # The idea is to check that in the tiled domain integrals over the
    # marked domains grow as expected
    def __volume(self, n, domain_id):
        # Relative to root
        path = './test/tile_1_narrow_GMSH306.h5'
        
        comm = pyMPI.COMM_WORLD
        h5 = HDF5File(comm, path, 'r')
        tile = Mesh()
        h5.read(tile, 'mesh', False)

        # Read in mesh data
        collection = MeshFunction('size_t', tile, tile.topology().dim(), 0)
        h5.read(collection, 'volumes')
        # Encode for evolution
        data = {}
        load_data(tile, collection, dim=3, data=data)

        # Eval volume integral with domain m and tags provided by f
        get_value = lambda m, f, domain_id=domain_id: assemble(
            CellVolume(m)*dx(domain=m, subdomain_data=f, subdomain_id=domain_id)
        )
                
        mesh, mesh_data = TileMesh(tile, (n, n+1), mesh_data=data)
        # Function on evolve tile, i.e. on mesh
        f = mf_from_data(mesh, mesh_data)[3]

        # On base we have
        tile_value = get_value(tile, collection)
        evolved_value = get_value(mesh, f)

        return (tile_value, evolved_value)

    def test_volume0(self):
        # Extracellular
        for n in (2, 3):
            tile_value, evolved_value = self.__volume(n, 0)
            # Sensible
            self.assertTrue(tile_value > 1E-13 and evolved_value > 1E-13)
            # Relative
            error = abs(n*(n+1)*tile_value - evolved_value)
            error /= n*(n+1)*tile_value

            self.assertTrue(error < 1E-13)

    def test_volume1(self):
        # Intracellular
        for n in (2, 3):
            tile_value, evolved_value = self.__volume(n, 1)
            # Sensible
            self.assertTrue(tile_value > 1E-13 and evolved_value > 1E-13)
            # Relative
            error = abs(n*(n+1)*tile_value - evolved_value)
            error /= n*(n+1)*tile_value

            self.assertTrue(error < 1E-13)


    def __area(self, n):
        # Relative to root
        path = './test/tile_1_narrow_GMSH306.h5'
        
        comm = pyMPI.COMM_WORLD
        h5 = HDF5File(comm, path, 'r')
        tile = Mesh()
        h5.read(tile, 'mesh', False)

        # Read in mesh data
        collection = MeshFunction('size_t', tile, tile.topology().dim()-1, 0)
        h5.read(collection, 'surfaces')
        # NOTE: the marked facets (cell surfaces) of the tile extend to
        # boundary. When evolving ths boundary one would be tricky to
        # account for in test
        DomainBoundary().mark(collection, 300)
        # Encode for evolution
        data = {}
        load_data(tile, collection, dim=2, data=data)

        # Eval volume integral with domain m and tags provided by f
        get_value = lambda m, f: assemble(
            avg(FacetArea(m))*dS(domain=m, subdomain_data=f, subdomain_id=1)
        )
                
        mesh, mesh_data = TileMesh(tile, (n, n+1), mesh_data=data)
        # Function on evolve tile, i.e. on mesh
        f = mf_from_data(mesh, mesh_data)[2]  # For facets

        # On base we have
        tile_value = get_value(tile, collection)
        evolved_value = get_value(mesh, f)

        return (tile_value, evolved_value)

    def test_area(self):
        for n in (2, 3):
            tile_value, evolved_value = self.__area(n)
            # Sensible
            self.assertTrue(tile_value > 1E-13 and evolved_value > 1E-13)
            # Relative
            error = abs(n*(n+1)*tile_value - evolved_value)
            error /= n*(n+1)*tile_value

            self.assertTrue(error < 1E-13)
