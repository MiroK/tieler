from dolfin import MeshFunction, Mesh, info, HDF5File
import numpy as np


class TileData(object):
    '''
    Tile represents a mesh of box which is meshed periodically in all its axis.
    '''
    def __init__(self, vertices, cells, mesh_data={}, periodic_maps=[]):
        self.vertices = vertices
        self.cells = cells
        self.mesh_data = mesh_data
        self.periodic_maps = periodic_maps

        self.dims = np.max(vertices, axis=0) - np.min(vertices, axis=0)

        
def Tile(mesh_path, mesh_data, TOL=1E-8):
    '''
    Construct a tile based on mesh in HDF5File
    
    INPUT:
      mesh_path = string to HDF5FILE
      mesh_data = dict data_set_name -> [dim] + list of tags to keep track of
                  
    OUTPUT:
      TileData instance
    '''
    tile = Mesh()
    comm = tile.mpi_comm()
    h5 = HDF5File(comm, mesh_path, 'r')
    h5.read(tile, 'mesh', False)

    mesh_data = {}
    for data_set in mesh_data:
        dim, tags = mesh_data[data_set][0], mesh_data[data_set][1:]
        # Output as last
        load_data(mesh, mesh_path, data_set, dim, tags, data)

    x = mesh.coordinates()
    min_x = np.min(x, axis=0)
    max_x = np.max(x, axis=0)
    shifts = max_x - min_x
    
    shifts_x = []  # Geometry
    vertex_mappings = []  # Periodicity
    # Compute geometrical shift for EVERY direction:
    for axis in range(len(shifts)):
        shift = shifts[axis]
        # Vector to shift cell vertices
        shift_x = np.zeros_like(shifts); shift_x[axis] = shift
        shifts_x.append(shift_x)

        # Compute periodicity in the vertices
        to_master = lambda x, shift=shift_x: x - shift
        # Mapping facets
        master = lambda x, axis=axis, A=min_x[axis], TOL=TOL: np.abs(x[:, axis]-A) < TOL
        slave = lambda x, axis=axis, A=max_x[axis], TOL=TOL: np.abs(x[:, axis]-A) < TOL

        error, vertex_mapping = compute_vertex_periodicity(tile, master, slave, to_master)
        # Fail when exended direction is no periodic
        assert error < 10*TOL, error
        
        vertex_mappings.append(vertex_mapping)

    return TileData(x, mesh.cells(), mesh_data, vertex_mappings)

# ---
        
def load_data(mesh, h5_file, data_set, dim, tags, data):
    '''
    Fill the data dictionary with data_set representing mesh function with 
    dim over mesh read from h5_file according to key spec expected by tiling 
    algorithm.
    '''
    mf = MeshFunction('size_t', mesh, dim, 0)
    h5_file.read(mf, data_set)

    tags = set(mf.array()) & set(tags)
    # Don't evolve zero - we initialize to it
    if 0 in tags: tags.remove(0)
    info('%s evolves tags %r' % (data_set, tags))
        
    # Data to evolve
    if dim > 0:
        mesh.init(dim, 0)
        e2v = mesh.topology()(dim, 0)

        for tag in tags:
            data[(dim, tag)] = np.array([e2v(e.index()) for e in SubsetIterator(mf, tag)],
                                        dtype='uintp')
    else:
        for tag in tags:
            data[(dim, tag)] = np.array([e.index() for e in SubsetIterator(mf, tag)],
                                        dtype='uintp')

    return data


def compute_vertex_periodicity(mesh, master, slave, to_master):
    '''Compute mapping from slave vertices to master vertices'''
    x = mesh.coordinates()

    master_vertices = (np.where(master(x))[0]).tolist()
    slave_vertices = (np.where(slave(x))[0]).tolist()

    assert len(master_vertices) == len(slave_vertices), (len(master_vertices), len(slave_vertices))

    error, mapping = 0., {}
    while slave_vertices:
        s = slave_vertices.pop()
        xs = x[s]
        mapped = to_master(xs)
        # Local to master_vertex_x
        master_vertex_x = x[master_vertices]

        dist = np.sqrt(np.sum((master_vertex_x - mapped)**2, axis=1))
        mapped_index = np.argmin(dist)
        # Wrt to vertex numbering
        m = master_vertices[mapped_index]
        mapping[s] = m
        error = max(error, dist[mapped_index])

        master_vertices.remove(m)
    assert not slave_vertices and not master_vertices
    
    return error, mapping

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import UnitSquareMesh, mpi_comm_world

    mesh = UnitSquareMesh(10, 10)
    f = HDF5File(mpi_comm_world(), 'test.h5', 'w')
    f.write(mesh, 'mesh')

    tile = Tile('test.h5', mesh_data={})

    print tile.periodic_maps
