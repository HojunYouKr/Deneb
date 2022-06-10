

Point(1) = {-1,-1,0,10};
Point(2) = {1,-1,0,10};
Point(3) = {-1,1,0,10};

Line(7) = {1,2};
Line(8) = {2,3};
Line(9) = {3,1};

Line Loop(14) = {7, 8, 9};

Plane Surface(16) = {14};

Physical Surface(1) = {16};
