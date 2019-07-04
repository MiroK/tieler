from dolfin import SubsetIterator, parameters, facets
import itertools, operator
import numpy as np


def facet_to_dof_map(V):
    '''facet index -> DLT0 dofs'''
    assert V.ufl_element().family() == 'HDiv Trace'
    assert V.ufl_element().degree() == 0
    assert V.ufl_element().value_shape() == ()

    dm = V.dofmap()
    
    mesh = V.mesh()
    cdim = mesh.topology().dim()
    fdim = cdim - 1
    mesh.init(fdim, cdim)
    mesh.init(cdim, fdim)
    
    c2f = mesh.topology()(cdim, fdim)
    mapping = np.zeros(mesh.num_entities(fdim), dtype='uintp')
    for f in facets(mesh):
        c = f.entities(cdim)[0]
        cell_dofs = dm.cell_dofs(c)

        fid = f.index()
        dof = cell_dofs[c2f(c).tolist().index(fid)]
        mapping[fid] = dof

    return mapping


def cells_of(cell_f, tag):
    '''Cells of tagged volume'''
    return list(map(operator.methodcaller('index'), SubsetIterator(cell_f, tag)))


def interfaces_between(cell_f, ext_tag, tagged_cells, int_tag, with_true_bdry):
    '''
    Which are the facets between which are shared by ext_tag and int_tag
    cells?
    '''
    mesh = cell_f.mesh()
    assert mesh.topology().dim() == cell_f.dim()
    assert parameters['ghost_mode'] == 'none'

    if tagged_cells is None:
        taggled_cells = cells_of(cell_f, int_tag)

    # They are facets that have at most one tagged cell
    cdim = mesh.topology().dim()
    fdim = cdim - 1

    mesh.init(fdim)
    mesh.init_global(fdim)
    mesh.init(cdim, fdim)
    mesh.init(fdim, cdim)

    f2c, c2f = mesh.topology()(fdim, cdim), mesh.topology()(cdim, fdim)
    cell_tags = cell_f.array()
    
    # The difficult part is handling facets that sit on cputs
    global_facet_idx = mesh.topology().global_indices(fdim)
    # These are all of them with local id -> {cpu}. 
    facet_cpus = mesh.topology().shared_entities(fdim)
    # Let me simplify
    facet_cpus = {k: next(iter(v)) for k, v in facet_cpus.items()}
    # Global-facet-id -> (local-facet-id, my conclusions about it)
    my_cpu_bdry = {}

    # Even if I dont have the marked volome I should still check my
    # cpu interfaces!
    facets = itertools.chain(itertools.chain(*map(c2f, tagged_cells)),
                             facet_cpus.keys())

    tags = set((int_tag, ext_tag))
        
    interface_facets = set()
    for facet in facets:
        if facet in interface_facets: continue
    
        values = list(cell_tags[f2c(facet)])
        # A facet is connected to at most 2 cells
        assert 0 < len(values) < 3, f2c(facet)

        # I own both
        if len(values) == 2:
            # Same tags mean that this is facet internal to the volume
            set(values) == tags and interface_facets.add(facet)
        # Maybe this is a true boundary or cpu boundary
        else:
            tag, = values
            # True bdry:
            if facet not in facet_cpus:
                (with_true_bdry and tag == int_tag) and interface_facets.add(facet)
            # Cpu
            else:
                # Store local facet index for syncing with other
                my_cpu_bdry[global_facet_idx[facet]] = (facet, tag)
    comm = mesh.mpi_comm()

    # Collect the conclusions
    all_cpu_bdries = comm.allgather(my_cpu_bdry)

    # Now overybody can work locally
    for facet in my_cpu_bdry:
        local_f, my_tag = my_cpu_bdry[facet]
        shared_cpu = all_cpu_bdries[facet_cpus[local_f]]
        if facet in shared_cpu:
            _, other_tag = shared_cpu[facet]

            set((my_tag, other_tag)) == tags and interface_facets.add(local_f)
    
    return np.fromiter(interface_facets, dtype='uintp')


def all_facets_of(mesh, tagged_cells, cell_f, tag):
    '''
    All facets of tagged cells
    '''
    # They are facets that have at most one tagged cell
    cdim = mesh.topology().dim()
    fdim = cdim - 1

    mesh.init(fdim)
    mesh.init(cdim, fdim)
    mesh.init(fdim, cdim)
    f2c, c2f = mesh.topology()(fdim, cdim), mesh.topology()(cdim, fdim)
    
    tagged_facets = list(map(c2f, tagged_cells))
    if tagged_facets:
        tagged_facets = np.hstack(tagged_facets).tolist()
    
    # NOTE: this must work also if the cpu would not have the volumes
    # to check; that is we check the CPU bdry facets and add those
    # that at some cpu are connected to tag cell
    facet_cpus = mesh.topology().shared_entities(fdim)
    # Let me simplify
    facet_cpus = {k: next(iter(v)) for k, v in facet_cpus.items()}

    global_facet_idx = mesh.topology().global_indices(fdim)
    
    my_cpu_bdry = {}
    for facet in facet_cpus:
        cell, = f2c(facet)
        cell_tag = cell_f[cell]

        # Store local facet index for syncing with other
        my_cpu_bdry[global_facet_idx[facet]] = (facet, cell_tag == tag)

    # Sync
    comm = mesh.mpi_comm()
    # Collect the conclusions
    all_cpu_bdries = comm.allgather(my_cpu_bdry)

    interface_facets = set()
    # Now overybody can work locally
    for facet in my_cpu_bdry:
        local_f, my_is = my_cpu_bdry[facet]
    
        shared_cpu = all_cpu_bdries[facet_cpus[local_f]]
        if facet in shared_cpu:
            _, other_is = shared_cpu[facet]

            my_is or other_is and interface_facets.add(local_f)
    
    tagged_facets = np.r_[tagged_facets, np.fromiter(interface_facets, dtype='uintp')]

    return np.fromiter(np.unique(tagged_facets), dtype='uintp')


def deactivate_volumes(cell_f, tag, predicate):
    '''Tagged cells of cell_f that satisfy the predicate'''
    return [cell.index()
            for cell in SubsetIterator(cell_f, tag)
            if predicate(cell)]

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *
    import os

    # parameters['ghost_mode'] = 'shared_facet'
    
    if False:
        mesh = UnitSquareMesh(32, 32)

        cell_f = MeshFunction('size_t', mesh, 2, 0)
        tag = 1

        CompiledSubDomain('x[0] < 0.5+DOLFIN_EPS || (x[1] > 0.5-DOLFIN_EPS)').mark(cell_f, tag)

        cpu_domains = MeshFunction('size_t', mesh, 2, 0)
        cpu_domains.array()[:] = mesh.mpi_comm().rank
        File('cpu_subs.pvd') << cpu_domains
    
        tagged_cells = cells_of(cell_f, tag)
        tagged_facets = interfaces_between(cell_f, ext_tag=1, tagged_cells=tagged_cells, int_tag=2)
        # tagged_facets = all_facets_of(mesh, tagged_cells)

        f = MeshFunction('size_t', mesh, 1, 0)
        f.array()[:] *= 0
        f.array()[tagged_facets] = 1

        g = MeshFunction('size_t', mesh, 2, 0)
        g.array()[tagged_cells] = 1

        File('bar.pvd') << g
        File('foo.pvd') << f

        dS_ = Measure('dS', domain=mesh, subdomain_data=f, subdomain_id=1)
        ds_ = Measure('ds', domain=mesh, subdomain_data=f, subdomain_id=1)

        print(assemble(Constant(1)*ds_+Constant(1)*dS_))

    # The tile dimension is 0.1 x 0.025 x 0.025
    # If the grid is
    #
    # oooooo
    # oxxxxo
    # oxxxxo
    # oooooo
    #
    # I want to deactive o

    mesh_file = '../tiles/tile_1_hein_GMSH307_12_6.h5'

    ncells_x, ncells_y = map(int, os.path.splitext(mesh_file)[0].split('_')[-2:])

    front = 0.025
    back = (ncells_y-1)*0.025
    left = 0.1
    right = (ncells_x-1)*0.1
    
    in_layer = CompiledSubDomain('x[0] > right || x[0] < left || x[1] > back || x[1] < front',
                                 left=left,
                                 right=right,
                                 front=front,
                                 back=back)
        
    # print((left, right), (front, back))
    comm = MPI.comm_world  # FIXME!
    h5 = HDF5File(comm, mesh_file, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)

    # print(mesh.coordinates().min(axis=0), mesh.coordinates().max(axis=0))
    surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    h5.read(surfaces, 'surfaces')

    volumes = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    h5.read(volumes, 'volumes')

    # print(mesh.coordinates().max(axis=0) - mesh.coordinates().min(axis=0))

    predicate = lambda c: in_layer.inside(c.midpoint().array(), False)
    # in_layer.mark(volumes, 2)
    bdry_cells = deactivate_volumes(volumes, tag=1, predicate=predicate)

    # Set the layer to 2 
    volumes.array()[bdry_cells] = 2
    # Now we have created a difference between layer and interior
    # Collect exterior ports of vol1 as surfaces of vol2 tagged with 1
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

    mesh_ranks = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    mesh_ranks.array()[:] = comm.rank

    
    h5 = HDF5File(comm, 'foo.h5', 'w')
    h5.write(mesh, 'mesh')
    h5.write(surfaces, 'surfaces')
    h5.write(volumes, 'volumes')
    
    File('bdry_volumes.pvd') << volumes
    File('bdry_surfs.pvd') << surfaces
    File('mesh_rankd.pvd') << mesh_ranks

    #dS = Measure('dS', domain=mesh, subdomain_data=surfaces)
    #print(assemble(Constant(1)*dS(1)))
