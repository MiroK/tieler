from dolfin import UnitSquareMesh, Timer
import matplotlib.pyplot as plt
import numpy as np
from tieler import TileMesh

tile = UnitSquareMesh(1, 1)

ns = [128, 256, 1024, 2048, 4096]
dts = []
for n in ns:
    shape = (n+1, n-1)  # To get odd as well
    
    t = Timer('s')
    mesh, mesh_data = TileMesh(tile, shape, mesh_data={})
    # Get tile for mesh write as well
    dts.append(t.stop())
    print mesh.num_cells()
        
a, b = np.polyfit(np.log(ns), np.log(dts), 1)

ns, dts = map(np.array, (ns, dts))
    
plt.figure()
plt.loglog(ns, dts, basex=2., basey=2., marker='x')
plt.loglog(ns, np.exp(b)*ns**a, basex=2., basey=2., linestyle='dashed',
           label='O(N^%.2f)' % a)
plt.ylabel('T')
plt.xlabel('N')
plt.legend(loc='best')
plt.show()
