from tieler.pattern import pattern_mesh
import dolfin as df
import numpy as np
from mpi4py import MPI

if __name__ == '__main__':
    from dolfin import HDF5File, Timer, File, Mesh, info
    import argparse

    parser = argparse.ArgumentParser(description='Put n tiles in x axis, m in y axis.')
    parser.add_argument('tile', type=str, help='H5 file that is the file')
    parser.add_argument('-m', type=int, default=1, help='Number of tiles in x dir')
    parser.add_argument('-n', type=int, default=1, help='Number of tiles in y dir')
    parser.add_argument('-k', type=int, default=1, help='Number of tiles in z dir')
    parser.add_argument('--save_pvd', dest='save_pvd', action='store_true')
    args = parser.parse_args()
    
    tile = args.tile
    m, n, k = args.m, args.n, args.k

    h5 = df.HDF5File(MPI.COMM_WORLD, tile, 'r')
    tile1 = df.Mesh()
    h5.read(tile1, 'mesh', False)

    mf1 = df.MeshFunction('size_t', tile1, 3, 0)
    h5.read(mf1, 'volumes')

    tiles = {1: mf1}
    shape = (k, m, n)
    pattern = np.ones(shape)

    foo = pattern_mesh(pattern, tiles)
    
    if args.save_pvd:
        df.File('result/x.pvd') << foo
    
    root = tile.split(".h5")[0]
    shape_str = ('_%d'*len(shape) % shape)
    h5_file = '%s%s.h5' % (root, shape_str)

    fout = df.HDF5File(MPI.COMM_WORLD, h5_file, 'w')
    fout.write(foo.mesh(), 'mesh')
    fout.write(foo, 'volumes')
    fout.close()
