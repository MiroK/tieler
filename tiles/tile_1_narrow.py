from math import sqrt
import numpy as np
import gmsh, sys

gmsh.initialize(sys.argv)

model = gmsh.model
occ = model.occ

tol_ = 1E-10  # Our tolerance for center of mass checking
gmsh.option.setNumber('Print.X3dPrecision', 1E-15)
gmsh.option.setNumber('Geometry.Tolerance', 1E-12)

# Cell dimensions
radius = 10
radius_x = 8
radius_y = 6
length = 100
length_x = 4
length_y = 4
padz = 4

depth = sqrt(radius*radius-radius_y*radius_y)
x_shift = length + 2*length_x
y_shift = 2*depth + 2*length_y

min_x = -length/2 - length_x
min_y = -depth - length_y
min_z = -radius - padz

max_x = length/2 + length_x
max_y = depth + length_y
max_z = radius + padz

#  The cell
v0 = occ.addCylinder(0-length/2, 0, 0, length, 0, 0, radius)
v1 = occ.addCylinder(0-length/2-length_x, 0, 0, length+2*length_x, 0, 0, radius_x)
v2 = occ.addCylinder(0, -depth-length_y, 0, 0, y_shift, 0, radius_y)
cell, _ = occ.fuse([(3, v0)], [(3, v1), (3, v2)])

# Bounding box
box = occ.addBox(min_x, min_y, min_z, max_x-min_x, max_y-min_y, max_z-min_z)

vols, _ = occ.fragment(cell, [(3, box)])
cell, box = vols

occ.synchronize()

# gmsh.fltk.initialize()
# gmsh.fltk.run()

lower_bounds, upper_bounds = (min_x, min_y, min_z), (max_x, max_y, max_z)
# Now for periodicity we want x
for obj in (cell, box):
    boundary = model.getBoundary([obj])
    print(obj)
    for axis in (0, 1, 2):
        lower = []
        for (dim, tag) in boundary:
            com = occ.getCenterOfMass(dim, abs(tag))
            abs(com[axis]-lower_bounds[axis]) < tol_ and lower.append((dim, abs(tag)))

        upper = []
        for (dim, tag) in boundary:
            com = occ.getCenterOfMass(dim, abs(tag))
            abs(com[axis]-upper_bounds[axis]) < tol_ and upper.append((dim, abs(tag)))

        if upper and lower:
            if len(upper) == len(lower) == 1:
                dimU, tagU = upper[0]
                dimL, tagL = lower[0]
                assert dimU == dimL == 2

                shift = np.zeros(3)
                shift[axis] = upper_bounds[axis] - lower_bounds[axis]
                transform = np.eye(4)
                transform[[0, 1, 2], -1] = shift
                transform = transform.flatten().tolist()
                
                gmsh.model.mesh.setPeriodic(dimU, [tagU], [tagL], transform)

occ.synchronize()

model.addPhysicalGroup(3, [cell[1]], 1)
model.addPhysicalGroup(3, [box[1]], 0)

cell_boundary = {abs(tag) for dim, tag in model.getBoundary([cell])}
box_boundary = {abs(tag) for dim, tag in model.getBoundary([box])}
interface = cell_boundary & box_boundary

model.addPhysicalGroup(2, list(interface), 1)

occ.synchronize()

model.mesh.generate(3)

gmsh.write('tile_1_narrow.msh')
