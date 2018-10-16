from tieler.tile_cpp import fill_mesh_function
    
from dolfin import MeshFunction, info, SubsetIterator
from collections import defaultdict
import numpy as np
import six


def mf_from_data(mesh, data):
    '''Build tdim mesh function from the data of TileMesh'''
    return _mx_from_data(mesh, data,
                         fill=fill_mesh_function,
                         init_container=lambda m, t: MeshFunction('size_t', m, t, 0))


def groupby(pairs, index):
    '''Organize pairs by pairs[index]'''
    groups = defaultdict(list)
    for pair in pairs: groups[pair[index]].append(pair)

    for item in six.iteritems(groups):
        yield item

        
def _mx_from_data(mesh, data, fill, init_container):
    '''Fill the container over mesh by data. Get back dict{tdim -> MeshFoo}'''
    try: 
        assert mesh.mpi_comm().tompi4py().size == 1
    # FEniCS 2018
    except AttributeError:
        assert mesh.mpi_comm().size == 1

    containers = {}
    # We have define entities in terms of vertex numbering
    # Order keys such by tdim (the first key)
    for tdim, keys in groupby(data.keys(), 0):
        # So we'll be getting the entity index by lookup
        mesh.init(tdim)
        mesh.init(0, tdim)
        # Build the meshfunction from data
        f = init_container(mesh, tdim)
        for key in keys:
            indices = data[key]
            indices.shape = (np.prod(indices.shape), )
            # These entity indices get the 'color'
            fill(mesh, indices, tdim, key[1], f)
        containers[tdim] = f

    return containers


def load_data(mesh, mesh_f, dim, data):
    '''
    Represent mesh_f over dim entities of mesh as collection of vertices.
    Can have mesh as mesh function or (h5_file, data_set)
    '''
    try:
        h5_file, data_set = mesh_f
        mf = MeshFunction('size_t', mesh, dim, 0)
        h5_file.read(mf, data_set)
    except ValueError:
        mf = mesh_f
        data_set = '%d dim entites' % mf.dim()
            
    # Mapf for encoding entities as vertices
    mesh.init(dim, 0)
    e2v = mesh.topology()(dim, 0)

    tags = set(mf.array())
    # Don't encode zero - we initialize to it
    if 0 in tags: tags.remove(0)
    info('%s evolves tags %r' % (data_set, tags))

    for tag in tags:
        data[(dim, tag)] = np.array([e2v(e.index()) for e in SubsetIterator(mf, tag)],
                                    dtype='uintp')
    return data
