from tile import TileMesh
from tile_data import mf_from_data, load_data
from dolfin import *
import numpy as np


tile = UnitSquareMesh(2, 2)

mf = MeshFunction('size_t', tile, tile.topology().dim()-1, 0)
CompiledSubDomain('near(x[0], 0.5) || near(x[1], 0.5)').mark(mf, 1)

mesh_data = {}
mesh_data = load_data(tile, mf, 1, mesh_data)

mesh, mesh_data = TileMesh(tile, shape=(23, 13), mesh_data=mesh_data)
fs = mf_from_data(mesh, mesh_data)

print mesh.coordinates().min(axis=0)
print mesh.coordinates().max(axis=0)
print 23*13*4, sum(1 for _ in SubsetIterator(fs[1], 1))

File('mesh_f.pvd') << fs[1]

# tile = UnitSquareMesh(1, 1)
# mesh, mesh_data = TileMesh(tile, shape=(13, 8), mesh_data={})

# if True:
#     ns = [128, 256, 1024, 2048, 4096]
#     dts = []
#     for n in ns:
#         shape = (n+1, n-1)
    
#         t = Timer('s')
#         mesh, mesh_data = TileMesh(tile, shape, mesh_data={})
#         dts.append(t.stop())
#         print mesh.num_cells()
        
#     import matplotlib.pyplot as plt
#     import numpy as np

#     a, b = np.polyfit(np.log(ns), np.log(dts), 1)
#     print a, b

#     ns, dts = map(np.array, (ns, dts))
    
#     plt.figure()
#     plt.loglog(ns, dts, basex=2., basey=2., marker='x')
#     plt.loglog(ns, np.exp(b)*ns**a, basex=2., basey=2., linestyle='dashed')

#     plt.show()
