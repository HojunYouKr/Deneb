

Point(1) = {-1,-1,-1,10};
Point(2) = {1,-1,-1,10};
Point(3) = {1,1,-1,10};
Point(4) = {-1,1,-1,10};
Point(5) = {-1,-1,1,10};
Point(6) = {1,-1,1,10};
Point(7) = {1,1,1,10};
Point(8) = {-1,1,1,10};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line Loop(13) = {1, 2, 3, 4};
Line Loop(14) = {5, 6, 7, 8};
Line Loop(15) = {1, 10, -5, -9};
Line Loop(16) = {2, 11, -6, -10};
Line Loop(17) = {3, 12, -7, -11};
Line Loop(18) = {4, 9, -8, -12};

Plane Surface(17) = {13};
Plane Surface(18) = {14};
Plane Surface(19) = {15};
Plane Surface(20) = {16};
Plane Surface(21) = {17};
Plane Surface(22) = {18};

Surface Loop(24) = {17, 18, 19, 20, 21, 22};
Volume(25) = {24};

Transfinite Line "*" = 1 Using Bump 0.25;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
