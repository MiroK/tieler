* TIELER [tyler]

Constructing FEniCS meshes for structured geometries by `tiling` the domain with 
some representative tile.

** Dependencies

The `master` branch is compatible with FEniCS 2017.2.0 (python2, no pybind11). The 
`pybind` branch has support for FEniCS 2018.+ (python3, pybind11). 

** Installation

Put this directory on python path, e.g. by running `source setup.rc`. For 
checking the installation execute `python -m unittest discover test`

** EMI usage

In order to get the mesh for EMI experiments navigate to `./tiles`. Make 
sure that your Gmsh version is at least 3.0.5 but less than 4.+. In the following 
the tile `tile_1_narrow_GMSH306.geo` is repeated 3 times in x and 4 times in y 
direction to create the mesh.

1. Get the mesh for the tile `gmsh -3 -clscale 0.3 tile_1_narrow_GMSH306.geo`
2. Wrap the mesh with data to H5 `python msh_convert.py tile_1_narrow_GMSH306.geo`
3. Tile the domain python tiled_mesh.py tile_1_narrow_GMSH306.h5 -n 3 -m 4
