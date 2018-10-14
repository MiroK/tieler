from dolfin import Mesh, MeshFunction, HDF5File, info
import subprocess, os


def convert(msh_file, h5_file):
    '''Temporary version of convert from msh to h5'''
    root, _ = os.path.splitext(msh_file)
    assert os.path.splitext(msh_file)[1] == '.msh'
    assert os.path.splitext(h5_file)[1] == '.h5'

    # Get the xml mesh
    xml_file = '.'.join([root, 'xml'])
    subprocess.call(['dolfin-convert %s %s' % (msh_file, xml_file)], shell=True)
    # Success?
    assert os.path.exists(xml_file)

    mesh = Mesh(xml_file)
    out = HDF5File(mesh.mpi_comm(), h5_file, 'w')
    out.write(mesh, 'mesh')

    info('Mesh has %d cells' % mesh.num_cells())
    info('Mesh size %g %g' % (mesh.hmin(), mesh.hmax()))
    
    # Save ALL data as facet_functions
    names = ('surfaces', 'volumes')
    for name, region in zip(names, ('facet_region.xml', 'physical_region.xml')):
        r_xml_file = '_'.join([root, region])

        f = MeshFunction('size_t', mesh, r_xml_file)
        out.write(f, name)

    return True
    

def cleanup(files=None, exts=()):
    '''Get rid of xml'''
    if files is not None:
        return map(os.remove, files)
    else:
        files = list(filter(lambda f: any(map(f.endswith, exts)), os.listdir('.')))
        info('Removing %r' % files) 
        return cleanup(files)
                    
# --------------------------------------------------------------------

if __name__ == '__main__':
    from mpi4py import MPI
    import argparse

    parser = argparse.ArgumentParser(description='Convert msh file to h5')
    parser.add_argument('io', type=str, nargs='+', help='input [output]')
    parser.add_argument('--cleanup', type=str, nargs='+',
                        help='extensions to delete', default=('.xml'))

    # Save the mesh markers for visualizing
    save_pvd_parser = parser.add_mutually_exclusive_group(required=False)
    save_pvd_parser.add_argument('--save_pvd', dest='save_pvd', action='store_true')
    save_pvd_parser.add_argument('--no_save_pvd', dest='save_pvd', action='store_false')
    parser.set_defaults(save_pvd=False)

    args = parser.parse_args()

    # Protecting self
    assert not(set(('geo', '.geo')) & set(args.cleanup))
    
    try:
        msh_file, h5_file = args.io[:2]
    except ValueError:
        msh_file = args.io[0]

        root, ext = os.path.splitext(msh_file)
        h5_file = '.'.join([root, 'h5'])

    assert convert(msh_file, h5_file)
    
    # VTK visualize tags
    if args.save_pvd:
        from dolfin import File 

        h5 = HDF5File(MPI.COMM_WORLD, h5_file, 'r')
        mesh = Mesh()
        h5.read(mesh, 'mesh', False)

        surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
        volumes = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)

        h5.read(surfaces, 'surfaces')
        h5.read(volumes, 'volumes')
        
        File('results/%s_surf.pvd' % root) << surfaces
        File('results/%s_vols.pvd' % root) << volumes    

    cleanup(exts=args.cleanup)
