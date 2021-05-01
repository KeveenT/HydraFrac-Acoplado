Point(1) = {0.0, 0.0, 0.0};

Point(2) = {10.0, 0.0, 0.0};

Point(3) = {30.0, 0.0, 0.0};

Point(4) = {30.0, 15.0, 0.0};

Point(5) = {30.0, 30.0, 0.0};

Point(6) = {10.0, 30.0, 0.0};

Point(7) = {0.0, 30.0, 0.0};

Point(8) = {0.0, 15.0005, 0.0};

Point(9) = {10.0, 15.0, 0.0};

Point(10) = {0.0, 14.9995, 0.0};

Line(1) = {1, 2};

Line(2) = {2, 3};

Line(3) = {3, 4};

Line(4) = {4, 5};

Line(5) = {5, 6};

Line(6) = {6, 7};

Line(7) = {7, 8};

Line(8) = {8, 9};

Line(9) = {9, 10};

Line(10) = {10, 1};

Line(11) = {2, 9};

Line(12) = {9, 6};

Line(13) = {9, 4};

Curve Loop(1) = {1, 11, 9, 10};

Curve Loop(2) = {2, 3, -13, -11};

Curve Loop(3) = {13, 4, 5, -12};

Curve Loop(4) = {8, 12, 6, 7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Physical Curve("South") = {1, 2};

Physical Curve("East") = {3, 4};

Physical Curve("North") = {5, 6};

Physical Curve("West") = {7, 10};

Physical Curve("Sup") = {8};

Physical Curve("Inf") = {9};

Physical Surface("Body") =  {1, 2, 3, 4};

Transfinite Curve 1 = 41 Using Progression 1;

Transfinite Curve 2 = 21 Using Progression 1;

Transfinite Curve 3 = 15 Using Progression 1;

Transfinite Curve 4 = 15 Using Progression 1;

Transfinite Curve 5 = 21 Using Progression 1;

Transfinite Curve 6 = 41 Using Progression 1;

Transfinite Curve 7 = 15 Using Progression 1;

Transfinite Curve 10 = 15 Using Progression 1;

Transfinite Curve 8 = 41 Using Progression 1;

Transfinite Curve 9 = 41 Using Progression 1;

Transfinite Curve 11 = 15 Using Progression 1;

Transfinite Curve 12 = 15 Using Progression 1;

Transfinite Curve 13 = 21 Using Progression 1;

Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
