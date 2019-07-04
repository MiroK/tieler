# This is an auxiliary for testing deactivation:
# mpirun -np X python3 deactivation.py
#
# deactivates the boundary cells and dums a new mesh foo.h5
# Here I test that no matter what X is the cells that should be there
# are indeed present.

from dolfin import *

comm = MPI.comm_world  # FIXME!
h5 = HDF5File(comm, 'foo.h5', 'r')
mesh = Mesh()
h5.read(mesh, 'mesh', False)


surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
h5.read(surfaces, 'surfaces')

volumes = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
h5.read(volumes, 'volumes')

print(assemble(Constant(1)*dx(domain=mesh, subdomain_data=volumes, subdomain_id=1)))

dS = dS(domain=mesh, subdomain_data=surfaces)
ds = ds(domain=mesh, subdomain_data=surfaces)

for i in (1,):
    print(i, assemble(Constant(1)*dS(i) + Constant(1)*ds(i)))
