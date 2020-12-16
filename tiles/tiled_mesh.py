from tieler import TileMesh, load_data, mf_from_data
from mpi4py import MPI as pyMPI


def get_comm_world():
    # FEniCS 2017
    try:
        from dolfin import mpi_comm_world 
        return mpi_comm_world()
    # FEniCS 2018
    except ImportError:
        return pyMPI.COMM_WORLD

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import HDF5File, Timer, File, Mesh, info
    import argparse, os

    parser = argparse.ArgumentParser(description='Put n tiles in x axis, m in y axis.')
    parser.add_argument('tile', type=str, help='H5 file that is the file')
    parser.add_argument('-m', type=int, default=1, help='Number of tiles in x dir')
    parser.add_argument('-n', type=int, default=1, help='Number of tiles in y dir')
    parser.add_argument('-k', type=int, default=1, help='Number of tiles in z dir')    
    # NOTE: give option to scale the tile rather than the final mesh
    parser.add_argument('-scale_x', type=float, default=1.,
                        help='Scale factor for coordinates')
    
    parser.add_argument('-facet_tags', type=str, default='surfaces',
                        help='name under which H5 stores facet tags')
    parser.add_argument('-cell_tags', type=str, default='volumes',
                        help='name under which H5 stores volume tags')

    save_pvd_parser = parser.add_mutually_exclusive_group(required=False)
    save_pvd_parser.add_argument('--save_pvd', dest='save_pvd', action='store_true')
    save_pvd_parser.add_argument('--no_save_pvd', dest='save_pvd', action='store_false')
    parser.set_defaults(save_pvd=False)

    # NOTE: give option to shift the tile rather than the final mesh
    shift_origin_parser = parser.add_mutually_exclusive_group(required=False)
    shift_origin_parser.add_argument('--shift_origin', dest='shift_origin', action='store_true')
    shift_origin_parser.add_argument('--no_shift_origin', dest='shift_origin', action='store_false')
    parser.set_defaults(shift_origin=True)
    
    args = parser.parse_args()

    # Some sanity
    root, ext = os.path.splitext(args.tile)
    assert args.m > 0 and args.n > 0
    assert ext == '.h5'

    if args.k == 1:
        shape = (args.m, args.n)
    else:
        shape = (args.m, args.n, args.k)
    
    # Load the tile mesh
    h5 = HDF5File(get_comm_world(), args.tile, 'r')
    tile = Mesh()
    h5.read(tile, 'mesh', False)

    # Shift and so
    x = tile.coordinates()
    xmin = x.min(axis=0)
    dx = x.max(axis=0) - xmin

    if args.shift_origin: x[:] -= xmin
    
    x[:] *= args.scale_x

    # Did it work?
    print(dx, 'vs', x.max(axis=0)-x.min(axis=0), xmin, x.min(axis=0))
    
    
    data = {}
    cell_dim = tile.topology().dim()
    facet_dim = cell_dim - 1

    if args.facet_tags: 
        data = load_data(tile, (h5, args.facet_tags), facet_dim, data)
    
    if args.cell_tags: 
        data = load_data(tile, (h5, args.cell_tags), cell_dim, data)

    t = Timer('tile')
    mesh, mesh_data = TileMesh(tile, shape, mesh_data=data)
    info('\nTiling took %g s; nvertices %d, nfacets %d, ncells %d' % (t.stop(),
                                                                      mesh.num_vertices(),
                                                                      mesh.init(2),
                                                                      mesh.num_cells()))

    x_ = mesh.coordinates()
    print('Final mesh size', x_.min(axis=0), x_.max(axis=0)-x_.min(axis=0))

    # Saving
    t = Timer('save')
    shape_str = ('%d'*len(shape) % shape)
    h5_file = '%s_%s.h5' % (root, shape_str)
        
    out = HDF5File(mesh.mpi_comm(), h5_file, 'w')
    out.write(mesh, 'mesh')
    
    tt = Timer('data')
    # To mesh functions
    if mesh_data:
        mfs = mf_from_data(mesh, mesh_data)

        for dim, name in (zip((facet_dim, cell_dim), (args.facet_tags, args.cell_tags))):
            if name:
                out.write(mfs[dim], name)
            
                if args.save_pvd:
                    File('%s_%s_%s.pvd' % (root, shape_str, name)) << mfs[dim]
                
    info('\t\tGetting data as MeshFoo took %g s' % tt.stop())
    
    info('\tSaving took %g' % t.stop())
