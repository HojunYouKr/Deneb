

Point(1) = {-1,-1,-1,10};
Point(2) = {1,-1,-1,10};
Point(3) = {-1,1,-1,10};
Point(4) = {-1,-1,1,10};

Line(7) = {1,2};
Line(8) = {1,3};
Line(9) = {1,4};
Line(10) = {2,3};
Line(11) = {2,4};
Line(12) = {3,4};

Line Loop(13) = {7, 10, -8};
Line Loop(14) = {7, 11, -9};
Line Loop(15) = {8, 12, -9};
Line Loop(16) = {10, 12, -11};

Plane Surface(17) = {13};
Plane Surface(18) = {14};
Plane Surface(19) = {15};
Plane Surface(20) = {16};



//Physical Surface(1) = {16};
//+
Surface Loop(21) = {17, 18, 20, 19};
//+
Volume(22) = {21};
//+
Volume(23) = {21};
