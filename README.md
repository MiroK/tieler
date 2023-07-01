# TIELER [tyler]

Constructing FEniCS meshes for structured geometries by `tiling` the domain with 
some representative tile.

## Dependencies

The `master` branch is compatible with FEniCS 2017.2.0 (python2, no pybind11)
and FEniCS 2018.+ (python3, pybind11). 

## Installation

Put this directory on python path, e.g. by running `source setup.rc`. For 
checking the installation execute `python -m unittest discover test` or
`python3 -m unittest discover test`

## EMI usage

In order to get the mesh for EMI experiments navigate to `./tiles`. Make 
sure that your Gmsh version is at least 3.0.6 but less than 4.+. In the following 
the tile `tile_1_narrow_GMSH306.geo` is repeated 3 times in x and 4 times in y 
direction to create the mesh.

1. Get the mesh for the tile: `gmsh -3 -clscale 0.3 tile_1_narrow_GMSH306.geo`
2. Wrap the mesh with data to H5: `python msh_convert.py tile_1_narrow_GMSH306.msh`
3. Tile the domain: `python tiled_mesh.py tile_1_narrow_GMSH306.h5 -m 3 -n 4 -scale_x 1E-3`

Where `scale_x` argument makes the units of mesh to be milimiters. Further, `python` is 
either python2 or python3. Note that Gmsh version should match that of the tile. 
(The .geo files have hardcoded mapping for periodic entities and their numbering seems 
to change between versions.). Also note that 1.-3. can only run in serial.

There is an optional step 4 which removes those (heart) cells of the mesh that touch
the mesh boundary 

4. `mpirun -np 3 python3 deactive_cells.py tile_1_narrow_GMSH306_3_3.h5 -m 3 -n 3`

## EMI usage alternatively
Instead of specifying the cell in the geo file we can use Gmsh Python API directly 

1. Get the mesh for the tile: `python3 tile_1_narrow.py -clscale 0.3 -format msh2`
2. Wrap the mesh with data to H5: `python msh_convert.py tile_1_narrow.msh`
3. Tile the domain: `python tiled_mesh.py tile_1_narrow.h5 -m 3 -n 4 -scale_x 1E-3`

Note that `-format msh2` flag is necessary for step 2 to work. For conversion betweem
msh and dolfin compatible formats you could instead use [meshio](https://github.com/nschloe/meshio).

## CI testing
This package is CI tested against FEniCS packages for `ubuntu 16.04 LTS` [![Build Status](https://travis-ci.org/MiroK/tieler.svg?branch=master)](https://travis-ci.org/MiroK/tieler)
