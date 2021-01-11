from tieler.periodicity import vertex_translate_stitch_mapping
from tieler.tile import merge, make_mesh
from collections import deque
import numpy as np



def pattern_mesh(pattern, tiles, axis=None):
    '''Build a mesh of tiles translate-stitched according to a pattern.'''
    # Tiles are cell functions and we return a cell function where the nonzero
    # subdomains in tiles are numbered as they appear in the pattern
    if axis is None:
        axis = np.arange(pattern.ndim)[::-1]
    assert len(axis) == pattern.ndim
    assert all(a >= 0 for a in axis)

    # An nd-dim pattern is patten_mesh of the stitched tiles - so we
    # reduce to base case
    if pattern.ndim > 1:
        # Glue along the first axis
        stitch_pattern, stitch_tiles = [], {}
        for tile_id, axis_pattern in enumerate(pattern):
            print(axis, tile_id, pattern)
            stitch_tiles[tile_id] = pattern_mesh(axis_pattern, tiles, axis[1:])
            stitch_pattern.append(tile_id)
        # We have now reduced the dim of pattern and built entirely new tiles
        # to be glued along the remainging axis
        return pattern_mesh(np.array(stitch_pattern), stitch_tiles, axis[:1])
    # The base case
    assert len(pattern) > 0
    return translate_stitch(pattern, tiles, axis[0])


def translate_stitch(pattern, tiles, axis, tagged_cells=[]):
    '''
    Union of tiles according to pattern along axis. Tagged_cells are the cells 
    indices (in the final mesh) of the tiles as they appear in the pattern. 
    NOTE: the method is ! on tiles and tagged_cells
    '''
    assert len(pattern) > 0
    pattern = deque(pattern)
    # In general case we will create a new tile to be inserted into the pattern
    key = max(pattern) + 1

    first = pattern.popleft()
    # For uninitialized computation the mappings of the (A, B) pattern
    # are the non-zero cells of A
    if not tagged_cells:
        values = tiles[first].array()
        for tag in sorted(np.unique(values)):
            tag and tagged_cells.append(np.where(values == tag)[0])
    # There was only one cell 
    if not pattern: return tiles[first]
    
    second = pattern.popleft()
    # Build a new piece to be inserted into a pattern. This new piaces
    # is made of firs 2 tiles
    first, second = tiles[first], tiles[second]
    # Now we glue
    mesh1, mesh2 = first.mesh(), second.mesh()
    gdim, tdim = mesh1.geometry().dim(), mesh1.topology().dim()
    assert gdim == mesh2.geometry().dim()
    assert tdim == mesh2.topology().dim()
    
    # Using the port vertices
    error, mapping = vertex_translate_stitch_mapping(mesh1, mesh2, axis)

    x1, x1_cells = mesh1.coordinates(), mesh1.cells()
    x2, x2_cells = mesh2.coordinates(), mesh2.cells()
    # And by translating 2 to the end of 1
    shift_size = x1.max(axis=0)[axis] - x2.min(axis=0)[axis]
    shift = np.zeros(gdim)
    shift[axis] = shift_size

    master_vertices, slave_vertices = map(list, zip(*mapping.items()))
    # Get vertices that make up mesh2 in the union
    x12, x12_cells, translated = merge(x1, x1_cells,
                                       shift,
                                       master_vertices, slave_vertices,
                                       x2, x2_cells)
    mesh12 = make_mesh(x12.ravel(),
                       np.asarray(x12_cells.ravel(), dtype='uintp'), tdim, gdim)
    # For consistent return type (value type of the dictionary) we construct
    # the mesh function ...
    mesh12_f = df.MeshFunction('size_t', mesh12, tdim, 0)
    # and we can already set tiles in it from the first one
    n = len(x1_cells)
    mesh12_array = mesh12_f.array()
    mesh12_array[:n] = first.array()
    
    # Since we always append x2_cells we just offset nonzero cells
    global_tag = len(tagged_cells)
    nz_tags = np.unique(second.array())
    for tag in sorted(nz_tags):
        if tag:
            cells = n + np.where(second.array() == tag)[0]
            tagged_cells.append(cells)
            mesh12_array[cells] = global_tag

            global_tag += 1

    pattern.appendleft(key)
    tiles[key] = mesh12_f
    # With the updated pattern and tiles we continue to reduce to len(pattern) == 1
    # case
    return translate_stitch(pattern, tiles, axis, tagged_cells)
    
# Merge cell functions and return one

# Shifts can be precomputed and stitch mappings reused
# Allow for None as pattern entry
# Speed up stitching
# Do not build more tiles than needed

# --------------------------------------------------------------------

if __name__ == '__main__':
    import dolfin as df

    mesh1 = df.UnitCubeMesh(4, 4, 4)
    mesh1_f = df.MeshFunction('size_t', mesh1, 3, 0)
    df.CompiledSubDomain('((0.25 - tol < x[0]) && (x[0] < 0.75 + tol)) && ((0.25 - tol < x[1]) && (x[1] < 0.75 + tol))',
                         tol=df.DOLFIN_EPS).mark(mesh1_f, 1)
    
    mesh2 = df.UnitCubeMesh(4, 4, 4)
    mesh2_f = df.MeshFunction('size_t', mesh2, 3, 0)
    df.CompiledSubDomain('((0.25 - tol < x[0]) && (x[0] < 0.75 + tol)) && ((0.25 - tol < x[1]) && (x[1] < 0.75 + tol))',
                         tol=df.DOLFIN_EPS).mark(mesh2_f, 1)

    mesh3_f = df.MeshFunction('size_t', mesh2, 3, 0)
    
    tiles = {0: mesh1_f, 1: mesh2_f, 2: mesh3_f}
    pattern = np.array([[[0, 1, 2, 0, 1, 1],
                         [0, 2, 1, 0, 2, 2],
                         [0, 1, 2, 0, 1, 1],
                         [0, 1, 2, 0, 2, 2]],
                        [[0, 2, 1, 0, 1, 1],
                        [0, 2, 1, 0, 2, 2],
                        [0, 1, 2, 0, 1, 1],
                        [0, 1, 2, 0, 2, 2]]])

    #pattern = np.array([[0, 1, 2, 0, 1, 1],
    #                     [0, 2, 1, 0, 2, 2],
    #                     [0, 1, 2, 0, 1, 1],
    #                     [0, 1, 2, 0, 2, 2]])
    
    #axis = 0
    #mappings = []
    #foo = translate_stitch(pattern, tiles, axis, mappings)    
    # shifts = pattern_shifts(pattern, tiles)

    foo = pattern_mesh(pattern, tiles)
    #axis = 0
    #aa = stitch(pattern, tiles, axis)

    #df.File('fpp.pvd') << aa
