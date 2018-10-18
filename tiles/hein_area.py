# This script illustrates how tweaking the Hein (superellipse) parameter
# changes the relative volume of the extracellular space
# 2.00 0.33
# 3.00 0.24
# 3.50 0.24
# 4.00 0.20
# 4.50 0.20
# 5.00 0.18
# 5.50 0.18

from msh_convert import convert
from dolfin import (Measure, assemble, File, SubsetIterator, sqrt, avg, TestFunction,
                    FunctionSpace)
import matplotlib.pyplot as plt
import subprocess, os
import numpy as np

ratios = []
exponents = [2, 3, 3.5, 4.0, 4.5, 5.0, 5.5]
for hein in exponents[-1:]:
    subprocess.call(['gmsh -3 -clscale 0.075 -setnumber Hein %d tile_1_hein.geo' % hein],
                    shell=True)

    
    outputs = convert('tile_1_hein.msh', 'tile_1_hein.h5')
    mesh, surfaces, volumes = outputs

    dx = Measure('dx', domain=mesh, subdomain_data=volumes)

    intra = assemble(1*dx(1))
    extra = assemble(1*dx(0))
    ratios.append(extra/(intra + extra))
    print('Extracellular volume %.2f' % ratios[-1])

    # Now want to get ratio of cell surface area to its volume
    # Want to keep track of port areas - these are boundary integrals
    ds = Measure('ds', domain=mesh, subdomain_data=surfaces)
    dS = Measure('dS', domain=mesh, subdomain_data=surfaces)

    area_no_port = assemble(1*dS(1))
    area_port = assemble(1*ds(1))

    volume = intra 
    print('Area/Volume %.3f, Area(!port)/Volume %.3f' % ((area_port+area_no_port)/volume,
                                                         area_no_port/volume))

    volume = intra + extra
    print('Area/Volume %.3f, Area(!port)/Volume %.3f' % ((area_port+area_no_port)/volume,
                                                         area_no_port/volume))

    # Mesh size
    mesh.init(2)
    mesh.init(2, 0)
    f2v = mesh.topology()(2, 0)

    tri_area = lambda x: np.linalg.norm(np.cross(x[0] - x[1], x[0] - x[2]), 2)/2
    x = mesh.coordinates()
    areas = [tri_area(x[f2v(f.index())]) for f in SubsetIterator(surfaces, 1)]
    
    min_l = sqrt(4*np.min(areas)/sqrt(3))
    max_l = sqrt(4*np.max(areas)/sqrt(3))

    print('Mesh size %g %g' % (min_l, max_l))
    print('Domain sizes %g %g %g' % tuple(np.max(x, axis=0)-np.min(x, axis=0)))

    File('foo.pvd') << volumes
    File('bar.pvd') << surfaces

    map(os.remove,
        filter(lambda f: any(f.endswith(ext) for ext in ('.xml', '.msh')), os.listdir('.')))

    
for er in zip(exponents, ratios): print '%.2f %.2f' % er
    
plt.figure()
plt.plot(exponents, ratios, '-bx')
plt.show()
