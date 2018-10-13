from dolfin import MeshFunction, SubsetIterator
import numpy as np


def compute_vertex_periodicity(mesh, master, slave, to_master):
    '''
    Compute mapping from slave vertices to master vertices

    INPUT:
      mesh = dolfin.Mesh
      master = dolfin.SubDomain instance (to mark master vertices)
      slave = dolfin.SubDomain instance (to mark slave vertices)
      to_mastar = x (a slave vertex coord) -> a master vertex coord

    OUTPUT:
      error = largest distance between master and slave vertex coordinates
      mapping = dict{slave vertex index -> master vertex index}
    '''
    f = MeshFunction('size_t', mesh, 0, 0)  # As a vertex function
    master.mark(f, 2)
    slave.mark(f, 3)

    master_vertices = [v.index() for v in SubsetIterator(f, 2)]
    slave_vertices = set(v.index() for v in SubsetIterator(f, 3))

    assert len(master_vertices) == len(slave_vertices), (len(master_vertices), len(slave_vertices))

    x = mesh.coordinates()

    error, mapping = 0., {}
    while slave_vertices:
        s = slave_vertices.pop()
        xs = x[s]  # Its coord
        mapped = to_master(xs)  # Coord in master domain
        # Local to master_vertex_x
        master_vertex_x = x[master_vertices]
        # Perform search to find closest master vertex
        dist = np.sqrt(np.sum((master_vertex_x - mapped)**2, axis=1))
        mapped_index = np.argmin(dist)
        # Wrt to vertex numbering
        m = master_vertices[mapped_index]
        mapping[s] = m
        # Do we have a match ?
        error = max(error, dist[mapped_index])
        # Assume we do so no more match possible
        master_vertices.remove(m)
    assert not slave_vertices and not master_vertices
    
    return error, mapping


def compute_entity_periodicity(tdim, mesh, master, slave, to_master):
    '''
    Compute mapping from slave tdim entities to master tdim entities

    INPUT:
      mesh = dolfin.Mesh
      master = dolfin.SubDomain instance (to mark master vertices)
      slave = dolfin.SubDomain instance (to mark slave vertices)
      to_mastar = x (a slave vertex coord) -> a master vertex coord

    OUTPUT:
      error = largest distance between master and slave vertex coordinates
      mapping = dict{slave tdim entity index -> master tdim entity index}
    '''
    assert 0 <= tdim < mesh.topology().dim()

    error, vertex_mapping = compute_vertex_periodicity(mesh, master, slave, to_master)
    # Done for vertices
    if tdim == 0:
        return error, vertex_mapping

    # Other entities are established using their definition in terms
    # of vertex indices
    f = MeshFunction('size_t', mesh, tdim, 0)
    master.mark(f, 2)
    slave.mark(f, 3)

    master_entities = (e.index() for e in SubsetIterator(f, 2))
    slave_entities = (e.index() for e in SubsetIterator(f, 3))
        
    mesh.init(tdim, 0)
    e2v = mesh.topology()(tdim, 0)
    # Define in terms of vertices. Invert for look up
    master_vertices = {tuple(sorted(e2v(e))): e for e in master_entities}
    # For slave we define via mapped vertices
    slave_entities = {e: tuple(sorted(vertex_mapping[v] for v in e2v(e))) for e in slave_entities}

    assert len(master_vertices) == len(slave_entities)

    mapping = {}
    while slave_entities:
        s, vertices = slave_entities.popitem()
        # Look up the master entity
        m = master_vertices[vertices]
        mapping[s] = m
        
        master_vertices.pop(vertices)
    assert not slave_entities and not master_vertices
    
    return error, mapping
