from tieler.pattern import pattern_mesh, get_pattern_shifts
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

xx = get_pattern_shifts(pattern, tiles, axis=0)


def make_consistent(shift_patterns, axis):
    if shift_patterns.ndim == 1:
        return shift_patterns

    new_pattern = []
    for i in range(shift_patterns.shape[axis]):
        slice_ = np.take(shift_patterns, i, axis)
        values = np.unique(slice_.flat)

        if np.any(np.isnan(values)):
            values = np.delete(values, np.where(np.isnan(values)))
        # Fail if there were more - inconsisten pattern            
        v, = values
        v, = values

        slice_[np.where(np.isnan(slice_))] = v
        new_pattern.append(slice_)
    return np.stack(new_pattern, axis=axis)

        
foo = pattern_mesh(pattern, tiles)
print(foo.mesh().num_cells())
df.File('x.pvd') << foo

# FIXME: shifts computing (np.take?)
#        use shifts in mesh making
#        making tiles with periodicity
#
#        tiled_mesh using gmsh API
