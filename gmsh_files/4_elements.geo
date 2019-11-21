//+ points
Point(1) = {0.0, 0.0, 0, 1.0};
//+
Point(2) = {1.0, 0.0, 0, 1.0};
//+
Point(3) = {1.0, 1.0, 0, 1.0};
//+
Point(4) = {0.0, 1.0, 0, 1.0};

//+lines
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};

//+surface
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};

//+Boundary
Physical Curve("Left") = {1};
//+
Physical Curve("Right") = {3};
//+
Physical Curve("Top") = {4};
//+
Physical Curve("Bottom") = {2};

//+
Physical Surface("plane") = {1};

//+Transfinite
Transfinite Curve {1, 3} = 17 Using Progression 1;
//+
Transfinite Curve {4, 2} = 17 Using Progression 1;
//+
Transfinite Surface {1};

//+Recombine
Recombine Surface {1};

