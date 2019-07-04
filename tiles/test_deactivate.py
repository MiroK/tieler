from dolfin import *
import sys


mesh_file = 'tile_1_hein_GMSH307_4_4_noBdry.h5'

parameters['ghost_mode'] = 'shared_facet'

comm = MPI.comm_world  # FIXME!
h5 = HDF5File(comm, mesh_file, 'r')
mesh = Mesh()
h5.read(mesh, 'mesh', False)

surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
h5.read(surfaces, 'surfaces')

volumes = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
h5.read(volumes, 'volumes')

# Volume of the remaining cells
print(assemble(Constant(1)*dx(domain=mesh, subdomain_data=volumes, subdomain_id=1)))

dS = dS(domain=mesh, subdomain_data=surfaces)
ds = ds(domain=mesh, subdomain_data=surfaces)

# Surface area of all the surfaces
print(1, assemble(Constant(1)*dS(1) + Constant(1)*ds(1)))
