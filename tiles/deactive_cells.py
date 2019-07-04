from tieler.deactivation import (deactivate_volumes, interfaces_between,
                                 all_facets_of)
from dolfin import CompiledSubDomain


def deactivate_cells(surfaces, volumes, ncells_x, ncells_y):
    '''The tile dimension is 0.1 x 0.025 x 0.025
    If the grid is
      oooooo
      oxxxxo
      oxxxxo
      oooooo
    I want to deactive o
    '''

    front = 0.025
    back = (ncells_y-1)*0.025
    left = 0.1
    right = (ncells_x-1)*0.1
    # Marking the boundary guys
    in_layer = CompiledSubDomain('x[0] > right || x[0] < left || x[1] > back || x[1] < front',
                                 left=left,
                                 right=right,
                                 front=front,
                                 back=back)
        
    predicate = lambda c: in_layer.inside(c.midpoint().array(), False)
    # Mark by only midpoint check
    bdry_cells = deactivate_volumes(volumes, tag=1, predicate=predicate)

    # Set the layer to 2 
    volumes.array()[bdry_cells] = 2
    # Now we have created a difference between layer and interior
    # Collect exterior ports of vol1 as surfaces between (1, 2)
    ext_ports = interfaces_between(volumes,
                                   ext_tag=1,
                                   tagged_cells=bdry_cells,
                                   int_tag=2,
                                   with_true_bdry=False)
    # Get all the facets in 2
    two_points = all_facets_of(mesh, bdry_cells, volumes, 2)
    # If I remove from them the ports I can safely kill of the remaining
    # facet; i.e. make them exterior
    surfaces.array()[list(set(two_points) - set(ext_ports))] = 0
    
    # Do that for volumes as well
    volumes.array()[bdry_cells] = 2

    return surfaces, volumes

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import HDF5File, Timer, File, Mesh, info, MeshFunction
    from tiled_mesh import get_comm_world
    import argparse, os

    parser = argparse.ArgumentParser(description='Remove boundary cells from mesh')
    parser.add_argument('mesh', type=str, help='H5 file that is the file')
    parser.add_argument('-m', type=int, default=1, help='Number of cells in x dir')
    parser.add_argument('-n', type=int, default=1, help='Number of cells in y dir')

    save_pvd_parser = parser.add_mutually_exclusive_group(required=False)
    save_pvd_parser.add_argument('--save_pvd', dest='save_pvd', action='store_true')
    save_pvd_parser.add_argument('--no_save_pvd', dest='save_pvd', action='store_false')
    parser.set_defaults(save_pvd=False)

    inplace_parser = parser.add_mutually_exclusive_group(required=False)
    inplace_parser.add_argument('--in_place', dest='in_place', action='store_true')
    inplace_parser.add_argument('--no_in_place', dest='in_place', action='store_false')
    parser.set_defaults(in_place=False)

    args = parser.parse_args()

    # Some sanity
    root, ext = os.path.splitext(args.mesh)
    assert args.m > 2 and args.n > 2
    assert ext == '.h5'

    # Typically we have files named foo_ncellsX_ncellsY.h5
    assert tuple(map(int, root.split('_')[-2:])) == (args.m, args.n)

    # Load the original 
    h5 = HDF5File(get_comm_world(), args.mesh, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)

    surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    h5.read(surfaces, 'surfaces')

    volumes = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    h5.read(volumes, 'volumes')
    
    h5.close()

    # Remove layer
    tt = Timer('cleanup')
    ncells_x, ncells_y = args.m, args.n
    surfaces, volumes = deactivate_cells(surfaces, volumes, ncells_x, ncells_y)

    info('Removing boundary took %g s' % tt.stop())

    # Write 
    tt = Timer('write')
    # Completely new
    if not args.in_place:
        h5_file = ''.join(['_'.join([root, 'noBdry']), ext])

        h5 = HDF5File(get_comm_world(), h5_file, 'w')
        # FIXME: replacing
        h5.write(mesh, 'mesh')
        h5.write(surfaces, 'surfaces')
        h5.write(volumes, 'volumes')
        h5.close()
    else:
        import h5py
        
        assert False
            
    info('Writing new entity functions took %g s' % tt.stop())

    # Optional visualization
    if args.save_pvd:
        File('%s_nobdry_volumes.pvd' % root) << volumes
        File('%s_nobdry_surfaces.pvd' % root)  << surfaces
