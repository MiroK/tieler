from tieler.pattern import pattern_mesh
import dolfin as df
import numpy as np


h5 = df.HDF5File(df.mpi_comm_world(), '../tiles/tile_1_narrow_GMSH307.h5', 'r')
tile1 = df.Mesh()
h5.read(tile1, 'mesh', False)

mf1 = df.MeshFunction('size_t', tile1, 3, 0)
h5.read(mf1, 'volumes')

tile2 = df.Mesh(tile1)
tile2.coordinates()[:, 0] *= 2.
mf2 = df.MeshFunction('size_t', tile2, 3, 0)
mf2.array()[:] = mf1.array()

mf3 = df.MeshFunction('size_t', tile1, 3, 0)

tiles = {1: mf3, 0: mf1, 2: mf2}
pattern = np.array([[[1, 0, 0, 2, 0, 0, 1],
                     [0, 1, 0, 2, 0, 1, 0],
                     [1, 1, 1, 2, 1, 1, 1],
                     [0, 1, 0, 2, 0, 1, 0],
                     [1, 0, 0, 2, 0, 0, 1]],
                    ##
                    [[0, 1, 1, 2, 1, 1, 0],
                     [1, 0, 1, 2, 1, 0, 1],
                     [0, 0, 0, 2, 0, 0, 0],
                     [1, 0, 1, 2, 1, 0, 1],
                     [0, 1, 1, 2, 1, 1, 0]]]                    
)

foo = pattern_mesh(pattern, tiles)
print(foo.mesh().num_cells())
df.File('x.pvd') << foo

# FIXME: shifts computing (np.take?)
#        use shifts in mesh making
#        making tiles with periodicity
#
#        tiled_mesh using gmsh API
