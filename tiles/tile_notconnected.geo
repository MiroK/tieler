Print.X3dPrecision = 1E-15;
Geometry.Tolerance = 1E-12;
Mesh.PreserveNumberingMsh2 = 1;

DefineConstant[
radius = {11.5, Name "radius of cell"}
radius_x = {8, Name "radius of connection in x direction"}
radius_y = {6, Name "radius of connection in y direction"}
length = {96, Name "length of cell (body)"}
length_x = {5.0, Name "length of connection in x direction"}
length_y = {5.0, Name "length of connection in y direction"}
padz = {5.0, Name "bounding box padding in z direction"}
Hein = {5.0, Name "Superellipse exponent"}
];

// The length units here micro meters
SetFactory("OpenCASCADE");

// Create the main body of the cell. Well quarete and then rotate for rest
r = Hein;
npts = 20;  // Control points for the ellipse
size = 10;

pointFirst = newp;
Point(pointFirst) = {0, radius, -length/2, size};
P[] = {pointFirst};

P0 = newp;
dtheta = Pi/2/(npts-1);
For i In {1:npts-2}
    theta = i*dtheta;
    xi = radius*Exp((2./r)*Log(Sin(theta)));
    yi = radius*Exp((2./r)*Log(Cos(theta)));

    Point(P0 + i - 1) = {xi, yi, -length/2, size};
    P[] += {P0 + i - 1};
EndFor
pointLast = newp;
Point(pointLast) = {radius, 0, -length/2, size};
P[] += {pointLast};

// The arch
L = newl;
Bezier(L) = {P[]};

// The quarter by extrusion
origin = newp;
Point(origin) = {0, 0, -length/2, size};
//+

toFirst = newl;
Line(toFirst) = {pointFirst, origin};
toLast = newl;
Line(toLast) = {pointLast, origin};

loop = newl;
Line Loop(loop) = {L, toFirst, -toLast};

base = news;
Plane Surface(base) = {loop};

quarter = Extrude {0, 0, length} { Surface{base}; };
// Finish the volume by symmetry operations
quarterY = Symmetry {0, 1, 0, 0} {Duplicata { Volume{quarter[1]}; } };

quarterX = Symmetry {1, 0, 0, 0} {Duplicata { Volume{quarter[1]}; } };
quarterYX = Symmetry {1, 0, 0, 0} {Duplicata { Volume{quarterY[]}; } };

// Now the cell is
cell() = BooleanUnion{ Volume{quarter[1]}; Delete; }{ Volume{quarterY, quarterX, quarterYX}; Delete; };

// Enclose in a bounding box
// Bounding box
box = newv;

min_x = -radius - padz;
max_x = radius + padz;
min_y = -radius - length_y;
max_y = radius + length_y;
min_z = -length/2-length_x;
max_z = length/2+length_x;

Box(box) = {min_x, min_y, min_z, max_x-min_x, max_y-min_y, max_z-min_z}; 

v() = BooleanFragments {Volume{box}; Delete; }{Volume{cell}; Delete; };
// Rotate to get ports in x, y
v[] = Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} {Volume{v[]}; };

// Tags
cell = v[0];
box = v[1];

Physical Volume(1) = {cell};
Physical Volume(0) = {box};

// X (this is gmsh 4.4.1 numbering, might need adjusting for other versions)
surfMaster = 35;
surfSlave = 33;
boundMaster[] = {64, 57, 56, 62};
boundSlave[] = {61, 59, 54, 60};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

// Y
surfMaster = 34;
surfSlave = 32;
boundMaster[] = {62, 63, 60, 55};
boundSlave[] = {57, 58, 59, 53};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

// Z
surfMaster = 36;
surfSlave = 31;
boundMaster[] = {58, 61, 63, 64};
boundSlave[] = {53, 54, 55, 56};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

cell_surf[] = Unique(Abs(Boundary{ Volume{cell}; }));  
Physical Surface(1) = {cell_surf[]};
