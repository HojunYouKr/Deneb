

Point(1) = {-1,-1,0,100};
Point(2) = {1,-1,0,100};
Point(3) = {1,1,0,100};
Point(4) = {-1,1,0,100};

Line(7) = {1,2};
Line(8) = {2,3};
Line(9) = {3,4};
Line(10) = {4,1};

Line Loop(14) = {7, 8, 9, 10};

Plane Surface(16) = {14};

Physical Surface(1) = {16};

// 2D algorithm : Delaunay for quad
// 2D recombination algorithm : Standard
