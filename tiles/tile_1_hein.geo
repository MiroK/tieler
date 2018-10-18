Print.X3dPrecision = 1E-15;
Geometry.Tolerance = 1E-12;

DefineConstant[
radius = {11.5, Name "radius of cell"}
radius_x = {8, Name "radius of connection in x direction"}
radius_y = {6, Name "radius of connection in y direction"}
length = {96, Name "length of cell (body)"}
length_x = {2, Name "length of connection in x direction"}
length_y = {1.0, Name "length of connection in y direction"}
padz = {1.0, Name "bounding box padding in z direction"}
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
body() = BooleanUnion{ Volume{quarter[1]}; Delete; }{ Volume{quarterY, quarterX, quarterYX}; Delete; };

// Let's add ports in z direction. Just a cylinder
port_x = newv;
Cylinder(port_x) = {0, 0, -length/2-length_x, 0, 0, length+2*length_x, radius_x, 2*Pi};

port_y = newv;
Cylinder(port_y) = {0, -radius-length_y, 0, 0, 2*radius+2*length_y, 0, radius_y, 2*Pi};

// And the cell
cell() = BooleanUnion{ Volume{body}; Delete; }{ Volume{port_x, port_y}; Delete; };

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

interfaces[] = Unique(Abs(Boundary{ Volume{cell}; }));  
boundary[] = Unique(Abs(Boundary{ Volume{box}; }));
boundary[] -= {interfaces[]};

Physical Surface(1) = {interfaces[]};

// Periodicity
// Ports in y, circle parts
surfMaster = 23;
surfSlave = 41;
boundMaster[] = {55};
boundSlave[] = {91};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

// Ports in y
surfMaster = 44;
surfSlave = 46;
boundMaster[] = {55, 97, 103, 102, 101};
boundSlave[] = {91, 99, 104, 107, 106};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

// Ports in x, circle parts
surfMaster = 42;
surfSlave = 40;
boundMaster[] = {96};
boundSlave[] = {86};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

// Ports in x
surfMaster = 45;
surfSlave = 47;
boundMaster[] = {98, 104, 105, 103, 96};
boundSlave[] = {100, 106, 108, 101, 86};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };
