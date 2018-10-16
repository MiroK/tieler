from dolfin import UnitSquareMesh, Timer
import matplotlib.pyplot as plt
import numpy as np
from tieler import TileMesh
import os, psutil
# For psutil see https://github.com/giampaolo/psutil
process = psutil.Process(os.getpid())

tile = UnitSquareMesh(1, 1)

ns = [128, 256, 1024, 2048, 4096]
dts = []
memory = []
for n in ns:
    shape = (n+1, n-1)  # To get odd as well
    
    t = Timer('s')
    mesh, mesh_data = TileMesh(tile, shape, mesh_data={})
    # Get tile for mesh write as well
    dts.append(t.stop())
    memory.append(process.memory_full_info().rss)
    print mesh.num_cells()
    print process.memory_full_info().rss
        
a, b = np.polyfit(np.log(ns), np.log(dts), 1)

ns, dts, memory = map(np.array, (ns, dts, memory))
    
plt.figure()
plt.loglog(ns, dts, basex=2., basey=2., marker='x')
plt.loglog(ns, np.exp(b)*ns**a, basex=2., basey=2., linestyle='dashed',
           label='O(N^%.2f)' % a)
plt.ylabel('T')
plt.xlabel('N')
plt.legend(loc='best')

memory /= 2**20
a, b = np.polyfit(np.log(ns), np.log(memory), 1)

plt.figure()
plt.loglog(ns, dts, basex=2., basey=2., marker='x')
plt.loglog(ns, np.exp(b)*ns**a, basex=2., basey=2., linestyle='dashed',
           label='O(N^%.2f)' % a)
plt.ylabel('T')
plt.xlabel('N')
plt.legend(loc='best')


plt.show()
