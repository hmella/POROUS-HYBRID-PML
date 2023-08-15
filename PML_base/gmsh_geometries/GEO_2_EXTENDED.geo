// Copyright (C) 2019 Hernan Mella

Include "dimensions";

// Depth of the second strata
h1 = 10.0;
h2 = 50.0;
h3 = 100.0;

// Get maximun characteristic length
cl = cl_base;

// Points
theta_1 = Asin(h1/R_ext);
theta_2 = Asin(h2/R_ext);
theta_3 = Asin(h3/R_ext);
theta_4 = Asin(Ly/R_ext);
Point(2) = {0.5*Lx, 0, 0, cl};
Point(3) = {0.5*Lx, -Ly, 0, cl};
Point(6) = {0.5*b+xc, 0, 0, cl};
Point(10) = {0.5*Lx, -h1, 0, cl};
Point(12) = {-0.5*Lx, 0, 0, cl};
Point(13) = {-0.5*Lx, -Ly, 0, cl};
Point(16) = {-0.5*b+xc, 0, 0, cl};
Point(17) = {-0.5*Lx, -h1, 0, cl};
Point(19) = {0.5*Lx, -h2, 0, cl};
Point(21) = {-0.5*Lx, -h2, 0, cl};
Point(23) = {0.5*Lx, -h3, 0, cl};
Point(25) = {-0.5*Lx, -h3, 0, cl};
Point(27) = {-R_ext, 0, 0, cl};
Point(28) = {R_ext, 0, 0, cl};
Point(29) = {xc, yc, 0, cl}; 
Point(30) = {-R_ext*Cos(theta_1), -R_ext*Sin(theta_1), 0, cl}; 
Point(31) = {R_ext*Cos(theta_1), -R_ext*Sin(theta_1), 0, cl};
Point(32) = {-R_ext*Cos(theta_2), -R_ext*Sin(theta_2), 0, cl}; 
Point(33) = {R_ext*Cos(theta_2), -R_ext*Sin(theta_2), 0, cl};
Point(34) = {-R_ext*Cos(theta_3), -R_ext*Sin(theta_3), 0, cl}; 
Point(35) = {R_ext*Cos(theta_3), -R_ext*Sin(theta_3), 0, cl};
Point(36) = {-R_ext*Cos(theta_4), -R_ext*Sin(theta_4), 0, cl}; 
Point(37) = {R_ext*Cos(theta_4), -R_ext*Sin(theta_4), 0, cl};
Point(38) = {0.5*Lx+LPML, 0, 0, cl};
Point(39) = {0.5*Lx+LPML, -Ly-LPML, 0, cl};
Point(40) = {0.5*Lx+LPML, -h1, 0, cl};
Point(41) = {-0.5*Lx-LPML, 0, 0, cl};
Point(42) = {-0.5*Lx-LPML, -Ly-LPML, 0, cl};
Point(43) = {-0.5*Lx-LPML, -h1, 0, cl};
Point(44) = {0.5*Lx+LPML, -h2, 0, cl};
Point(45) = {-0.5*Lx-LPML, -h2, 0, cl};
Point(46) = {0.5*Lx+LPML, -h3, 0, cl};
Point(47) = {-0.5*Lx-LPML, -h3, 0, cl};
Point(48) = {0, 0, 0, cl};


// Generate interface lines to achieve a fixed resolution
// Add points to the interface to enforce the same interior mesh
// Left vertical line
nPoints1 = Ceil(h1/cl); dx = h1/nPoints1;
For i In {1 : nPoints1-1}
  pList1[i-1] = newp;
  Point(pList1[i-1]) = {-0.5*Lx, -i*dx, 0, cl};
EndFor

nPoints2 = Ceil((h2-h1)/cl); dx = (h2-h1)/nPoints2;
For i In {1 : nPoints2-1}
  pList2[i-1] = newp;
  Point(pList2[i-1]) = {-0.5*Lx, -h1-i*dx, 0, cl};
EndFor

nPoints3 = Ceil((h3-h2)/cl); dx = (h3-h2)/nPoints3;
For i In {1 : nPoints3-1}
  pList3[i-1] = newp;
  Point(pList3[i-1]) = {-0.5*Lx, -h2-i*dx, 0, cl};
EndFor

nPoints4 = Ceil((Ly-h3)/cl); dx = (Ly-h3)/nPoints4;
For i In {1 : nPoints4-1}
  pList4[i-1] = newp;
  Point(pList4[i-1]) = {-0.5*Lx, -h3-i*dx, 0, cl};
EndFor

// Bottom horizontal line
nPoints5 = Ceil(Lx/cl); dx = Lx/nPoints5;
For i In {1 : nPoints5-1}
  pList5[i-1] = newp;
  Point(pList5[i-1]) = {-0.5*Lx+i*dx, -Ly, 0, cl};
EndFor

// Right vertical line
nPoints6 = Ceil((Ly-h3)/cl); dx = (Ly-h3)/nPoints6;
For i In {1 : nPoints6-1}
  pList6[i-1] = newp;
  Point(pList6[i-1]) = {0.5*Lx, -Ly+i*dx, 0, cl};
EndFor

nPoints7 = Ceil((h3-h2)/cl); dx = (h3-h2)/nPoints7;
For i In {1 : nPoints7-1}
  pList7[i-1] = newp;
  Point(pList7[i-1]) = {0.5*Lx, -h3+i*dx, 0, cl};
EndFor

nPoints8 = Ceil((h2-h1)/cl); dx = (h2-h1)/nPoints8;
For i In {1 : nPoints8-1}
  pList8[i-1] = newp;
  Point(pList8[i-1]) = {0.5*Lx, -h2+i*dx, 0, cl};
EndFor

nPoints9 = Ceil(h1/cl); dx = h1/nPoints9;
For i In {1 : nPoints9-1}
  pList9[i-1] = newp;
  Point(pList9[i-1]) = {0.5*Lx, -h1+i*dx, 0, cl};
EndFor

// Generate interface lines in the exterior PML layer
// to achieve a fixed resolution
// Add points to the interface to enforce the same interior mesh
// Left vertical line
nPoints10 = Ceil(h1/cl); dx = h1/nPoints10;
For i In {1 : nPoints10-1}
  pList10[i-1] = newp;
  Point(pList10[i-1]) = {-0.5*Lx-LPML, -i*dx, 0, cl};
EndFor

nPoints11 = Ceil((h2-h1)/cl); dx = (h2-h1)/nPoints11;
For i In {1 : nPoints11-1}
  pList11[i-1] = newp;
  Point(pList11[i-1]) = {-0.5*Lx-LPML, -h1-i*dx, 0, cl};
EndFor

nPoints12 = Ceil((h3-h2)/cl); dx = (h3-h2)/nPoints12;
For i In {1 : nPoints12-1}
  pList12[i-1] = newp;
  Point(pList12[i-1]) = {-0.5*Lx-LPML, -h2-i*dx, 0, cl};
EndFor

nPoints13 = Ceil((Ly+LPML-h3)/cl); dx = (Ly+LPML-h3)/nPoints13;
For i In {1 : nPoints13-1}
  pList13[i-1] = newp;
  Point(pList13[i-1]) = {-0.5*Lx-LPML, -h3-i*dx, 0, cl};
EndFor

// Bottom horizontal line
nPoints14 = Ceil((Lx+2*LPML)/cl); dx = (Lx+2*LPML)/nPoints14;
For i In {1 : nPoints14-1}
  pList14[i-1] = newp;
  Point(pList14[i-1]) = {-0.5*Lx-LPML+i*dx, -Ly-LPML, 0, cl};
EndFor

// Right vertical line
nPoints15 = Ceil((Ly+LPML-h3)/cl); dx = (Ly+LPML-h3)/nPoints15;
For i In {1 : nPoints15-1}
  pList15[i-1] = newp;
  Point(pList15[i-1]) = {0.5*Lx+LPML, -Ly-LPML+i*dx, 0, cl};
EndFor

nPoints16 = Ceil((h3-h2)/cl); dx = (h3-h2)/nPoints16;
For i In {1 : nPoints16-1}
  pList16[i-1] = newp;
  Point(pList16[i-1]) = {0.5*Lx+LPML, -h3+i*dx, 0, cl};
EndFor

nPoints17 = Ceil((h2-h1)/cl); dx = (h2-h1)/nPoints17;
For i In {1 : nPoints17-1}
  pList17[i-1] = newp;
  Point(pList17[i-1]) = {0.5*Lx+LPML, -h2+i*dx, 0, cl};
EndFor

nPoints18 = Ceil(h1/cl); dx = h1/nPoints18;
For i In {1 : nPoints18-1}
  pList18[i-1] = newp;
  Point(pList18[i-1]) = {0.5*Lx+LPML, -h1+i*dx, 0, cl};
EndFor


// Lines (interface boundary)
Line(1) = {12, pList1[0]};
Spline(2) = pList1[];
Line(3) = {pList1[nPoints1-2],17};

Line(4) = {17,pList2[0]};
Spline(5) = pList2[];
Line(6) = {pList2[nPoints2-2],21};

Line(7) = {21,pList3[0]};
Spline(8) = pList3[];
Line(9) = {pList3[nPoints3-2],25};

Line(10) = {25,pList4[0]};
Spline(11) = pList4[];
Line(12) = {pList4[nPoints4-2],13};

Line(13) = {13,pList5[0]};
Spline(14) = pList5[];
Line(15) = {pList5[nPoints5-2],3};

Line(16) = {3,pList6[0]};
Spline(17) = pList6[];
Line(18) = {pList6[nPoints6-2],23};

Line(19) = {23,pList7[0]};
Spline(20) = pList7[];
Line(21) = {pList7[nPoints7-2],19};

Line(22) = {19,pList8[0]};
Spline(23) = pList8[];
Line(24) = {pList8[nPoints8-2],10};

Line(25) = {10,pList9[0]};
Spline(26) = pList9[];
Line(27) = {pList9[nPoints9-2],2};

// Lines (PML exterior boundary)
Line(28) = {41, pList10[0]};
Spline(29) = pList10[];
Line(30) = {pList10[nPoints10-2],43};

Line(31) = {43,pList11[0]};
Spline(32) = pList11[];
Line(33) = {pList11[nPoints11-2],45};

Line(34) = {45,pList12[0]};
Spline(35) = pList12[];
Line(36) = {pList12[nPoints12-2],47};

Line(37) = {47,pList13[0]};
Spline(38) = pList13[];
Line(39) = {pList13[nPoints13-2],42};

Line(40) = {42,pList14[0]};
Spline(41) = pList14[];
Line(42) = {pList14[nPoints14-2],39};

Line(43) = {39,pList15[0]};
Spline(44) = pList15[];
Line(45) = {pList15[nPoints15-2],46};

Line(46) = {46,pList16[0]};
Spline(47) = pList16[];
Line(48) = {pList16[nPoints16-2],44};

Line(49) = {44,pList17[0]};
Spline(50) = pList17[];
Line(51) = {pList17[nPoints17-2],40};

Line(52) = {40,pList18[0]};
Spline(53) = pList18[];
Line(54) = {pList18[nPoints18-2],38};

// Lines (remaining)
Circle(55) = {27, 48, 30};
Circle(56) = {30, 48, 32};
Circle(57) = {32, 48, 34};
Circle(58) = {34, 48, 35};
Circle(59) = {35, 48, 33};
Circle(60) = {33, 48, 31};
Circle(61) = {31, 48, 28};
Line(62) = {27, 41};
Line(63) = {41, 12};
Line(64) = {12, 16};
Line(65) = {16, 6}; // Neumann
Line(66) = {6, 2};
Line(67) = {2, 38};
Line(68) = {38, 28};
Line(69) = {31, 40};
Line(70) = {40, 10};
Line(71) = {10, 17};
Line(72) = {17, 43};
Line(73) = {43, 30};
Line(74) = {32, 45};
Line(75) = {45, 21};
Line(76) = {21, 19};
Line(77) = {19, 44};
Line(78) = {44, 33};
Line(79) = {35, 46};
Line(80) = {46, 23};
Line(81) = {23, 25};
Line(82) = {25, 47};
Line(83) = {47, 34};

// Plane surfaces
Curve Loop(1) = {55, -73, -30, -29, -28, -62};
Plane Surface(1) = {1};
Curve Loop(2) = {28, 29, 30, -72, -3, -2, -1, -63};
Plane Surface(2) = {2};
Curve Loop(3) = {64, 65, 66, -27, -26, -25, 71, -3, -2, -1};
Plane Surface(3) = {3};
Curve Loop(4) = {25, 26, 27, 67, -54, -53, -52, 70};
Plane Surface(4) = {4};
Curve Loop(5) = {68, -61, 69, 52, 53, 54};
Plane Surface(5) = {5};
Curve Loop(6) = {56, 74, -33, -32, -31, 73};
Plane Surface(6) = {6};
Curve Loop(7) = {72, 31, 32, 33, 75, -6, -5, -4};
Plane Surface(7) = {7};
Curve Loop(8) = {71, 4, 5, 6, 76, 22, 23, 24};
Plane Surface(8) = {8};
Curve Loop(9) = {70, -24, -23, -22, 77, 49, 50, 51};
Plane Surface(9) = {9};
Curve Loop(10) = {69, -51, -50, -49, 78, 60};
Plane Surface(10) = {10};
Curve Loop(11) = {57, -83, -36, -35, -34, -74};
Plane Surface(11) = {11};
Curve Loop(12) = {75, 7, 8, 9, 82, -36, -35, -34};
Plane Surface(12) = {12};
Curve Loop(13) = {8, 9, -81, 19, 20, 21, -76, 7};
Plane Surface(13) = {13};
Curve Loop(14) = {77, -48, -47, -46, 80, 19, 20, 21};
Plane Surface(14) = {14};
Curve Loop(15) = {78, -59, 79, 46, 47, 48};
Plane Surface(15) = {15};
Curve Loop(16) = {58, 79, -45, -44, -43, -42, -41, -40, -39, -38, -37, 83};
Plane Surface(16) = {16};
Curve Loop(17) = {38, 39, 40, 41, 42, 43, 44, 45, 80, -18, -17, -16, -15, -14, -13, -12, -11, -10, 82, 37};
Plane Surface(17) = {17};
Curve Loop(18) = {11, 12, 13, 14, 15, 16, 17, 18, 81, 10};
Plane Surface(18) = {18};

// Physical surfaces
If (flag==0) // Extended domain
  For i In {1 : 18}
    Physical Surface(i) = {i};
  EndFor
ElseIf (flag==1) // PML domain
  Physical Surface(1) = {2};
  Physical Surface(2) = {3}; // RD
  Physical Surface(3) = {4};
  Physical Surface(4) = {7};
  Physical Surface(5) = {8}; // RD
  Physical Surface(6) = {9};
  Physical Surface(7) = {12};
  Physical Surface(8) = {13}; // RD
  Physical Surface(9) = {14};
  Physical Surface(10) = {17};
  Physical Surface(11) = {18};  // RD
ElseIf (flag==2) // Paraxial domain
  Physical Surface(1) = {3};
  Physical Surface(2) = {8};
  Physical Surface(3) = {13};
  Physical Surface(4) = {18};
EndIf

// Physical lines
Physical Line(3) = {65};                     // Neumann boundary (load)
If (flag==0) // Extended domain
  Physical Line(1) = {55,56,57,58,59,60,61}; // Dirichlet boundary
  Physical Line(4) = {62,63,64,66,67,68};    // Neumann boundary (remaining)
ElseIf (flag==1) // PML domain
  Physical Line(1) = {28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54}; // Dirichlet boundary
  Physical Line(2) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};          // Interface boundary  
  Physical Line(4) = {63,64,66,67};  // Neumann boundary (remaining)
ElseIf (flag==2) // Paraxial domain
  Physical Line(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};          // Interface boundary  
  Physical Line(4) = {64,66};  // Neumann boundary (remaining)
EndIf

Coherence;

// Refinement of the surroundings of the source
// Field[2] = Box;
// Field[2].VIn = cl_source;
// Field[2].VOut = cl_base;
// Field[2].XMin = -5.0*cl_source+xc;
// Field[2].XMax = +5.0*cl_source+xc;
// Field[2].YMin = -5.0*cl_source;
// Field[2].YMax = 0.0;
// Field[2].Thickness = 7.0*cl_source;
Field[1] = Distance;
Field[1]. PointsList = {29};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = cl_source;
Field[2].SizeMax = cl_base;
Field[2].DistMin = 5.0*cl_source;
Field[2].DistMax = 5.0*cl_base;
// Field[2].StopAtDistMax = 5.0*cl_base;

// Refinement of the surroundings of air water-transitions
Field[3] = Distance;
Field[3].CurvesList = {1,2,3,9,12,13,57,59};
Field[3].NumPointsPerCurve = 1000;

Field[4] = Threshold;
Field[4].InField = 3;
// Field[4].SizeMin = cl_base/2;
Field[4].SizeMin = cl_base;
Field[4].SizeMax = cl_base;
Field[4].DistMin = 2.0*cl_base;
Field[4].DistMax = 3.0*cl_base;
Field[4].StopAtDistMax = 3.0*cl_base;

// Finally, let's use the minimum of all the fields as the background mesh size
// field
Field[5] = Min;
Field[5].FieldsList = {2};
Background Field = 5;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms
MeshAlgorithm Surface {3,8,13,18} = 6; //5;         // RD surfaces
MeshAlgorithm Surface {2,4,7,9,12,14,17} = 6; //9;  // PML surfaces
MeshAlgorithm Surface {1,5,6,10,11,15,16} = 6; //5; // Remaining
Mesh.Smoothing = 1;