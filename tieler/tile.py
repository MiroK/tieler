from tieler.tile_cpp import fill_mesh

from tieler.periodicity import compute_vertex_periodicity

from dolfin import Mesh, info, Timer, CompiledSubDomain

from collections import namedtuple
from itertools import izip
from copy import deepcopy
import numpy as np

# MAIN IDEA:
# We have x, cells, ... representing a single cell. To get the tile
# repeated ntimes I don't simply stack next to each other.
# Consider 1+1+1+1+1+1+1+1 has 7 additions. But (4+4) where 4 = 2+2 where
# 2 = 1+1 is three! So providided that `+` scales reasonably with input
# size this is more efficient. Here + is really unary. For input other
# than power of 2 it is necessary to extend it to binary, e.g. 5 = (2+2)+1

def powers2(num):
    '''List (descending) of powers of two in the number'''
    powers = [int(power)
              for power, value in enumerate(reversed(format(num, 'b')))
              if value != '0']
    info('Number of binary + %d' % (len(powers) - 1))
    info('Largest tile does %d unary +' % (max(powers) - 1))
    
    return powers[::-1]


# Basic data structure representing tile consists of coordinates and
# cells as indices in to the coord array. master/slave vertices define a
# map for gluing tiles in certain direction. The periodic maps for the
# remaining dirs are in vertex_mappings. Marked entities (encode by their
# vertex numbering) are in data
Tile = namedtuple('Tile', ('coords', 'cells',
                           'master_vertices', 'slave_vertices', 'mappings',
                           'data'))


def make_tile(x, cells, master_vertices, slave_vertices, vertex_mappings, data):
    '''Freeze the tile from input data.'''
    # As the data further evolves I need to make copy
    return Tile(deepcopy(x), deepcopy(cells),
                np.copy(master_vertices), np.copy(slave_vertices),
                [vm.copy() for vm in vertex_mappings],
                data.copy())


def merge(x, x_cells, shift, master_vertices, slave_vertices, y=None, y_cells=None):
    '''
    The + operation for tiles (x, x_cells) [(y, y_cells)]. Tiles are glued 
    by shifting y by shift. Return x, cells for the resulting merged mesh 
    and a map from tile y local vertex numbering to the global numbering 
    of the mesh
    '''
    if y is None and y_cells is None: y, y_cells = x, x_cells

    n = len(y)
    # To make the tile piece we add all but the master vertices (as these
    # are already (slaves) in x-tile
    new_vertices = set(range(n))
    new_vertices.difference_update(master_vertices)
    new_vertices = np.fromiter(new_vertices, dtype=int)
    assert np.all(np.diff(new_vertices) > 0)
        
    # Offset the free
    translate = np.empty(n, dtype=int)
    translate[new_vertices] = len(x) + np.arange(len(new_vertices))
    # Those at master positions take slave values
    translate[master_vertices] = slave_vertices

    # Verices of the glued tiles
    x = np.vstack([x, y[new_vertices] + shift])

    # Cells of the glued tiles
    new_cells = np.empty_like(y_cells)
    new_cells.ravel()[:] = translate[y_cells.flat]

    x_cells = np.vstack([x_cells, new_cells])

    return x, x_cells, translate


def TileMesh(tile, shape, mesh_data={}, TOL=1E-9):
    '''
    [tile tile;
     tile tile;
     tile tile;
     tile tile]

    The shape is an ntuple describing the number of pieces put next 
    to each other in the i-th axis. mesh_data : (tdim, tag) -> [entities] 
    is the way to encode mesh data of the tile.
    '''
    # Sanity for glueing
    gdim = tile.geometry().dim()
    assert len(shape) <= gdim, (shape, gdim)
    # While evolve is general mesh writing is limited to simplices only (FIXME)
    # so we bail out early
    assert str(tile.ufl_cell()) in ('interval', 'triangle', 'tetrahedron')

    t = Timer('evolve')
    # Do nothing
    if all(v == 1 for v in shape): return tile, mesh_data

    # We want to evolve cells, vertices of the mesh using geometry information
    # and periodicity info
    x = tile.coordinates()
    min_x = np.min(x, axis=0)
    max_x = np.max(x, axis=0)
    shifts = max_x - min_x
    
    shifts_x = []  # Geometry
    vertex_mappings = []  # Periodicity
    # Compute geometrical shift for EVERY direction:
    for axis in range(len(shape)):
        shift = shifts[axis]
        # Vector to shift cell vertices
        shift_x = np.zeros(gdim); shift_x[axis] = shift
        shifts_x.append(shift_x)

        # Compute periodicity in the vertices
        to_master = lambda x, shift=shift_x: x - shift
        # Mapping facets
        master = CompiledSubDomain('near(x[i], A, tol)', i=axis, A=min_x[axis], tol=TOL)
        slave = CompiledSubDomain('near(x[i], A, tol)', i=axis, A=max_x[axis], tol=TOL)

        error, vertex_mapping = compute_vertex_periodicity(tile, master, slave, to_master)
        # Fail when exended direction is no periodic
        assert error < 10*TOL, error
        
        vertex_mappings.append(vertex_mapping)
    # The final piece of information is cells
    cells = np.empty((tile.num_cells(), tile.ufl_cell().num_vertices()), dtype='uintp')
    cells.ravel()[:] = tile.cells().flat
        
    # Evolve the mesh by tiling along the last direction in shape
    while shape:
        x, cells, vertex_mappings, shape = \
            evolve(x, cells, vertex_mappings, shape, shifts_x, mesh_data=mesh_data)
        
    info('\tEvolve took %g s ' % t.stop())

    # We evolve data but the mesh function is made out of it outside
    mesh = make_mesh(x, cells, tdim=tile.topology().dim(), gdim=gdim)

    return mesh, mesh_data

        
def evolve(x, cells, vertex_mappings, shape, shifts_x, mesh_data={}):
    '''Evolve tile [and its data]along the last axis'''
    axis, gdim = len(shape) - 1, x.shape[1]
    assert gdim > axis >= 0

    # Do nothing for 1
    if shape[axis] == 1:
        vertex_mappings.pop()  # No longer needed
        shifts_x.pop()  # Do not touch x and cells
        return x, cells, vertex_mappings, shape[:-1]

    # Glue how many times
    refine = shape[axis]
    # Get the sizes of tiles (the real sizes are powers of two of these)
    # The idea is to do unary + on each base tile and power. Then binary
    # + on the resulting tile
    tile_sizes = powers2(refine)
    # Iff refine was a power of two then there is no need to merge as
    # at the end of tile evolution is the final mesh data. Otherwise
    # we develop the tiles in their LOCAL numbering => they each start
    # from their copy
    vertex_mapping, shift_x = vertex_mappings.pop(), shifts_x.pop()
    master_vertices = vertex_mapping.values()
    slave_vertices = vertex_mapping.keys()

    tiles = []
    # Are we even or odd (to keep the initial tile)
    if 0 in tile_sizes:
        tiles = [make_tile(x, cells,
                           master_vertices, slave_vertices, vertex_mappings,
                           mesh_data)]
        # Done with this size
        tile_sizes.pop()

    # We need to evolve tiles of sizes which are power of 2. The idea
    # is to get them while evolving to the largest one
    size = 1
    max_size = max(tile_sizes)

    target_size = tile_sizes.pop()
    # The body of the loop corresponds to unary + 
    while size <= max_size:
        x, cells, translate = merge(x, cells, shift_x, master_vertices, slave_vertices)

        # For the directions that do not evolve we add the new periodic pairs
        for vm in vertex_mappings:
            keys, values = np.array(vm.items()).T
            vm.update(dict(izip(translate[keys], translate[values])))

        # Update the periodicty mapping - slaves are new
        slave_vertices = translate[slave_vertices]
        
        # Add the entities defined in terms of the vertices
        if mesh_data: mesh_data = evolve_data(mesh_data, translate)

        # We might be at a needed size
        if size == target_size:
            tiles.append(make_tile(
                x, cells, master_vertices, slave_vertices, vertex_mappings, mesh_data))
            
            if tile_sizes: target_size = tile_sizes.pop()
      
        # Iterate
        size += 1  # Remember powers of 2  (like divide)
        shift_x *= 2  # so we double
    assert not tile_sizes

    # Now that the tiles are prepare we want to glue them
    tile0 = tiles.pop()
    # No glueing of one tile
    if not tiles:
        return tile0.coords, tile0.cells, tile0.mappings, shape[:-1]

    # Otherwise extend data of the largest tile
    x, cells, vertex_mappings = tile0.coords, tile0.cells, tile0.mappings
    slave_vertices = tile0.slave_vertices

    # Now the binary +
    dx = shift_x / 2**max_size  # Recover length of the base tile
    shift_x = 0

    shifts = (2**p for p in powers2(refine))
    while tiles:
        next_tile = tiles.pop()
        # Shift for next tile to be added to x
        shift_x = shift_x + dx*next(shifts)        
        x, cells, translate = merge(x, cells, shift_x,
                                    next_tile.master_vertices, slave_vertices,
                                    next_tile.coords, next_tile.cells)

        # Updata periodicity mappings of next tile using new global map
        for vm, next_vm in zip(vertex_mappings, next_tile.mappings):
            keys, values = np.array(next_vm.items()).T
            vm.update(dict(izip(translate[keys], translate[values])))

        # Data evolve
        if mesh_data: mesh_data = evolve_data(mesh_data, translate, next_tile.data)
        
        # New global slabes
        slave_vertices = translate[next_tile.slave_vertices]

    return x, cells, vertex_mappings, shape[:-1]


def evolve_data(data, mapping, other=None):
    '''
    Add data of new tile. A new tile is this one or other. Mapping has 
    local to global from new tile to the summed tile.
    '''
    if other is None: other = data
    
    for key in data.keys():
        old = data[key]
            
        new = np.empty_like(other[key])
        new.ravel()[:] = mapping[other[key].flat]
        data[key] = np.vstack([old, new])
    return data


def make_mesh(coordinates, cells, tdim, gdim):
    '''Mesh by MeshEditor from vertices and cells'''
    mesh = Mesh()
    assert mesh.mpi_comm().tompi4py().size == 1

    fill_mesh(coordinates.flatten(), cells.flatten(), tdim, gdim, mesh)
    
    return mesh
