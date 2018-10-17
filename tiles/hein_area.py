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
from dolfin import Measure, assemble, File
import matplotlib.pyplot as plt
import subprocess, os

ratios = []
exponents = [2, 3, 3.5, 4.0, 4.5, 5.0, 5.5]
for hein in exponents:
    subprocess.call(['gmsh -3 -clscale 0.1 -setnumber Hein %d tile_1_hein.geo' % hein],
                    shell=True)

    
    outputs = convert('tile_1_hein.msh', 'tile_1_hein.h5')
    mesh, _, volumes = outputs

    dx = Measure('dx', domain=mesh, subdomain_data=volumes)

    intra = assemble(1*dx(1))
    extra = assemble(1*dx(0))
    ratios.append(extra/(intra + extra))
    print('Extracellular form %.2f' % ratios[-1])

    File('foo.pvd') << volumes


    map(os.remove,
        filter(lambda f: any(f.endswith(ext) for ext in ('.xml', '.msh')), os.listdir('.')))

    
for er in zip(exponents, ratios): print '%.2f %.2f' % er
    
plt.figure()
plt.plot(exponents, ratios, '-bx')
plt.show()
