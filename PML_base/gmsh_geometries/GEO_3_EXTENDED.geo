// Copyright (C) 2019 Hernan Mella

Include "dimensions";

// Get maximun characteristic length
cl = cl_base;

// Depth of each stratum
y0 = 0.0;
y1 = -10;
y2 = -50;
y3 = -100;

// Atan amplitude and intersection
A = 0.3;
B = y3 - A*Ly*Atan(15*(-0.5*Lx-0.15*Lx)/Lx);

// Points
Point(1) = {-0.5*Lx-LPML, 0, 0, cl};
Point(2) = {-0.5*Lx, 0, 0, cl};
Point(3) = {0.5*Lx, -Ly, 0, cl};
Point(4) = {-0.5*Lx, -Ly, 0, cl};
Point(7) = {0.5*Lx+LPML, -Ly-LPML, 0, cl};  // exterior
Point(8) = {-0.5*Lx-LPML, -Ly-LPML, 0, cl}; // exterior

// Intersection points between stratums
x_ref_1 = Tan((y0 - B)/(A*Ly))/(15/Lx) + Lx*0.15;
x_ref_2 = Tan((y1 - B)/(A*Ly))/(15/Lx) + Lx*0.15;
x_ref_3 = Tan((y2 - B)/(A*Ly))/(15/Lx) + Lx*0.15;
Point(9) = {x_ref_1, y0, 0, cl};
Point(10) = {x_ref_2, y1, 0, cl};
Point(11) = {x_ref_3, y2, 0, cl};

// Exterior points to extend the stratum layers
Point(12) = {-0.5*Lx-LPML, y1, 0.0, cl};
Point(13) = {-0.5*Lx, y1, 0.0, cl};
Point(14) = {-0.5*Lx-LPML, y2, 0.0, cl};
Point(15) = {-0.5*Lx, y2, 0.0, cl};
y4 = A*Ly*Atan(15*(0.5*Lx-Lx*0.15)/Lx) + B;
Point(16) = {0.5*Lx + LPML, y4, 0.0, cl};
y5 = A*Ly*Atan(15*(-0.5*Lx-Lx*0.15)/Lx) + B;
Point(17) = {-0.5*Lx - LPML, y5, 0.0, cl};

// Points for Neumannn boundary
Point(18) = {xc - 0.5*b, 0, 0, cl};
Point(10000) = {xc, 0, 0, cl};
Point(19) = {xc + 0.5*b, 0, 0, cl};

// Center point for circles
Point(20) = {0.0, y4, 0.0, cl};


// Exterior points to extend the stratum layers
y4 = A*Ly*Atan(15*(0.5*Lx-Lx*0.15)/Lx) + B;
y5 = A*Ly*Atan(15*(-0.5*Lx-Lx*0.15)/Lx) + B;
theta_1 = Asin((y1-y4)/R_ext);
theta_2 = Asin((y2-y4)/R_ext);
theta_4 = Asin((y4)/R_ext);
theta_5 = Asin((y5-y4)/R_ext);
Point(21) = {-R_ext*Cos(theta_1), y1, 0.0, cl};
Point(22) = {-R_ext*Cos(theta_2), y2, 0.0, cl};
Point(23) = {-R_ext*Cos(theta_5), y5, 0.0, cl};
Point(24) = {-R_ext*Cos(theta_4), 0, 0, cl}; // OK
Point(25) = {R_ext, y4, 0, cl};  // OK
Point(26) = {0.0, -(R_ext-y4), 0, cl}; // OK

////////////////////////////////////////////////////
//  Curved stratum curves
///////////////////////////////////////////////////

// Surface curve to the left of the load boundary
nPoints1 = 5; // Number of discretization points
For i In {1 : nPoints1+1}
  x = x_ref_2 + (x_ref_1 - x_ref_2)*(i-1)/(nPoints1);
  If (i == 1)
    pList1[i-1] = 10;
  ElseIf (i == nPoints1 + 1)
    pList1[i-1] = 9;
  Else
    pList1[i-1] = newp;
    Point(pList1[i-1]) = {x,
                      A*Ly*Atan(15*(x - 0.15*Lx)/Lx) + B,
                      0,
                      cl};
  EndIf
EndFor

// Surface curve to the left of the load boundary
nPoints2 = 25; // Number of discretization points
For i In {1 : nPoints2+1}
  x = -0.5*Lx + (0.5*Lx + x_ref_3)*(i-1)/(nPoints2);
  If (i == nPoints2 + 1)
    pList2[i-1] = 11;
  Else
    pList2[i-1] = newp;
    Point(pList2[i-1]) = {x,
                  A*Ly*Atan(15*(x - 0.15*Lx)/Lx) + B,
                  0,
                  cl};
  EndIf
EndFor

// Surface curve to the left of the load boundary
nPoints3 = 15; // Number of discretization points
For i In {1 : nPoints3+1}
  x = x_ref_3 + (x_ref_2 - x_ref_3)*(i-1)/(nPoints3);
  If (i == 1)
    pList3[i-1] = 11;
  ElseIf (i == nPoints3 + 1)
    pList3[i-1] = 10;
  Else
    pList3[i-1] = newp;
    Point(pList3[i-1]) = {x,
                      A*Ly*Atan(15*(x - 0.15*Lx)/Lx) + B,
                      0,
                      cl};
  EndIf
EndFor

// Surface curve to the left of the load boundary
nPoints4 = 25; // Number of discretization points
For i In {1 : nPoints4+1}
  x = x_ref_1 + (0.5*Lx - x_ref_1)*(i-1)/(nPoints4);
  If (i == 1)
    pList4[i-1] = 9;
  Else
    pList4[i-1] = newp;
    Point(pList4[i-1]) = {x,
                      A*Ly*Atan(15*(x - 0.15*Lx)/Lx) + B,
                      0,
                      cl};
  EndIf
EndFor

// Generate interface lines to achieve a fixed resolution
// Add points to the interface to enforce the same interior mesh
// Left vertical line
nPoints5 = Ceil(Abs(y1)/cl); dx = Abs(y1)/nPoints5;
For i In {1 : nPoints5-1}
  pList5[i-1] = newp;
  Point(pList5[i-1]) = {-0.5*Lx, -i*dx, 0, cl};
EndFor

nPoints6 = Ceil(Abs(y2-y1)/cl); dx = Abs(y2-y1)/nPoints6;
For i In {1 : nPoints6-1}
  pList6[i-1] = newp;
  Point(pList6[i-1]) = {-0.5*Lx, y1-i*dx, 0, cl};
EndFor

nPoints7 = Ceil(Abs(y3-y2)/cl); dx = Abs(y3-y2)/nPoints7;
For i In {1 : nPoints7-1}
  pList7[i-1] = newp;
  Point(pList7[i-1]) = {-0.5*Lx, y2-i*dx, 0, cl};
EndFor

nPoints8 = Ceil(Abs(-Ly-y3)/cl); dx = Abs(-Ly-y3)/nPoints8;
For i In {1 : nPoints8-1}
  pList8[i-1] = newp;
  Point(pList8[i-1]) = {-0.5*Lx, y3-i*dx, 0, cl};
EndFor

// Bottom horizontal line
nPoints9 = Ceil(Lx/cl); dx = Lx/nPoints9;
For i In {1 : nPoints9-1}
  pList9[i-1] = newp;
  Point(pList9[i-1]) = {-0.5*Lx+i*dx, -Ly, 0, cl};
EndFor

// Right vertical line
nPoints10 = Ceil(Abs(Ly + y4)/cl); dx = Abs(Ly+y4)/nPoints10;
For i In {1 : nPoints10-1}
  pList10[i-1] = newp;
  Point(pList10[i-1]) = {0.5*Lx, -Ly+i*dx, 0, cl};
EndFor

// Generate interface lines in the exterior PML layer
// to achieve a fixed resolution
// Add points to the interface to enforce the same interior mesh
// Left vertical line
nPoints11 = Ceil(Abs(y1)/cl); dx = Abs(y1)/nPoints11;
For i In {1 : nPoints11-1}
  pList11[i-1] = newp;
  Point(pList11[i-1]) = {-0.5*Lx-LPML, -i*dx, 0, cl};
EndFor

nPoints12 = Ceil(Abs(y2-y1)/cl); dx = Abs(y2-y1)/nPoints12;
For i In {1 : nPoints12-1}
  pList12[i-1] = newp;
  Point(pList12[i-1]) = {-0.5*Lx-LPML, y1-i*dx, 0, cl};
EndFor

nPoints13 = Ceil(Abs(y3-y2)/cl); dx = Abs(y3-y2)/nPoints13;
For i In {1 : nPoints13-1}
  pList13[i-1] = newp;
  Point(pList13[i-1]) = {-0.5*Lx-LPML, y2-i*dx, 0, cl};
EndFor

nPoints14 = Ceil(Abs(-Ly-LPML-y3)/cl); dx = Abs(-Ly-LPML-y3)/nPoints14;
For i In {1 : nPoints14-1}
  pList14[i-1] = newp;
  Point(pList14[i-1]) = {-0.5*Lx-LPML, y3-i*dx, 0, cl};
EndFor

// Bottom horizontal line
nPoints15 = Ceil((Lx+2*LPML)/cl); dx = (Lx+2*LPML)/nPoints15;
For i In {1 : nPoints15-1}
  pList15[i-1] = newp;
  Point(pList15[i-1]) = {-0.5*Lx-LPML+i*dx, -Ly-LPML, 0, cl};
EndFor

// Right vertical line
nPoints16 = Ceil(Abs(Ly+y4+LPML)/cl); dx = Abs(Ly+y4+LPML)/nPoints16;
For i In {1 : nPoints16-1}
  pList16[i-1] = newp;
  Point(pList16[i-1]) = {0.5*Lx+LPML, -Ly-LPML+i*dx, 0, cl};
EndFor


// Splines for stratum curves
Spline(1) = pList1[];
Spline(2) = pList2[];
Spline(3) = pList3[];
Spline(4) = pList4[];

// Lines
Line(5) = {2, 1};
Line(6) = {16, pList4(nPoints4)};
Line(7) = {2, 18};
Line(8) = {19, pList1[nPoints1]};
Line(9) = {13, 12};
Line(10) = {14, 15};
Line(11) = {15, pList2[nPoints2]};
Line(12) = {17, pList2[0]};
Line(13) = {pList3[nPoints3], 13};
Circle(14) = {24, 20, 21};
Circle(15) = {21, 20, 22};
Circle(16) = {22, 20, 23};
Circle(17) = {23, 20, 26};
Circle(18) = {26, 20, 25};

// Lines (interface boundary)
Line(19) = {2, pList5[0]};
Spline(20) = pList5[];
Line(21) = {pList5[nPoints5-2],13};

Line(22) = {13,pList6[0]};
Spline(23) = pList6[];
Line(24) = {pList6[nPoints6-2],15};

Line(25) = {15,pList7[0]};
Spline(26) = pList7[];
Line(27) = {pList7[nPoints7-2],pList2[0]};

Line(28) = {pList2[0],pList8[0]};
Spline(29) = pList8[];
Line(30) = {pList8[nPoints8-2],4};

Line(31) = {4,pList9[0]};
Spline(32) = pList9[];
Line(33) = {pList9[nPoints9-2],3};

Line(34) = {3,pList10[0]};
Spline(35) = pList10[];
Line(36) = {pList10[nPoints10-2],pList4[nPoints4]};

// Lines (PML exterior boundary)
Line(37) = {1, pList11[0]};
Spline(38) = pList11[];
Line(39) = {pList11[nPoints11-2],12};

Line(40) = {12,pList12[0]};
Spline(41) = pList12[];
Line(42) = {pList12[nPoints12-2],14};

Line(43) = {14,pList13[0]};
Spline(44) = pList13[];
Line(45) = {pList13[nPoints13-2],17};

Line(46) = {17,pList14[0]};
Spline(47) = pList14[];
Line(48) = {pList14[nPoints14-2],8};

Line(49) = {8,pList15[0]};
Spline(50) = pList15[];
Line(51) = {pList15[nPoints15-2],7};

Line(52) = {7,pList16[0]};
Spline(53) = pList16[];
Line(54) = {pList16[nPoints16-2],16};

Line(55) = {24, 1};
Line(56) = {22, 14};
Line(57) = {23, 17};
Line(58) = {16, 25};
Line(59) = {21, 12};

// Neumann boundary line (load)
Line(60) = {18, 19};


// Plane surfaces
Curve Loop(1) = {14, 59, -39, -38, -37, -55};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 37, 38, 39, -9, -21, -20, -19};
Plane Surface(2) = {2};
Curve Loop(3) = {20, 21, -13, 1, -8, -60, -7, 19};
Plane Surface(3) = {3};
Curve Loop(4) = {15, 56, -42, -41, -40, -59};
Plane Surface(4) = {4};
Curve Loop(5) = {9, 40, 41, 42, 10, -24, -23, -22};
Plane Surface(5) = {5};
Curve Loop(6) = {23, 24, 11, 3, 13, 22};
Plane Surface(6) = {6};
Curve Loop(7) = {16, 57, -45, -44, -43, -56};
Plane Surface(7) = {7};
Curve Loop(8) = {10, 25, 26, 27, -12, -45, -44, -43};
Plane Surface(8) = {8};
Curve Loop(9) = {26, 27, 2, -11, 25};
Plane Surface(9) = {9};
Curve Loop(10) = {17, 18, -58, -54, -53, -52, -51, -50, -49, -48, -47, -46, -57};
Plane Surface(10) = {10};
Curve Loop(11) = {32, 33, 34, 35, 36, -6, -54, -53, -52, -51, -50, -49, -48, -47, -46, 12, 28, 29, 30, 31};
Plane Surface(11) = {11};
Curve Loop(12) = {2, 3, 1, 4, -36, -35, -34, -33, -32, -31, -30, -29, -28};
Plane Surface(12) = {12};

// Physical surfaces
If (flag==0) // Extended domain
  For i In {1 : 12}
    Physical Surface(i) = {i};
  EndFor
ElseIf (flag==1) // PML domain
  Physical Surface(1) = {2};
  Physical Surface(2) = {3}; // RD
  Physical Surface(3) = {5};
  Physical Surface(4) = {6}; // RD
  Physical Surface(5) = {8};
  Physical Surface(6) = {9}; // RD
  Physical Surface(7) = {11};
  Physical Surface(8) = {12}; // RD
ElseIf (flag==2) // Paraxial domain
  Physical Surface(1) = {3};
  Physical Surface(2) = {6};
  Physical Surface(3) = {9};
  Physical Surface(4) = {12};
EndIf

// Physical lines
Physical Line(3) = {60};                     // Neumann boundary (load)
If (flag==0) // Extended domain
  Physical Line(1) = {14,15,16,17,18};  // Dirichlet boundary
  Physical Line(4) = {4,5,6,7,8,55,58}; // Neumann boundary (remaining)
ElseIf (flag==1) // PML domain
  Physical Line(1) = {37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54}; // Dirichlet boundary
  Physical Line(2) = {19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36}; // Interface boundary  
  Physical Line(4) = {4,5,6,7,8};  // Neumann boundary (remaining)
ElseIf (flag==2) // Paraxial domain
  Physical Line(1) = {19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36};          // Interface boundary  
  Physical Line(4) = {4,7,8};  // Neumann boundary (remaining)
EndIf

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
Field[1]. PointsList = {10000};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = cl_source;
Field[2].SizeMax = cl_base;
Field[2].DistMin = 5.0*cl_source;
Field[2].DistMax = 5.0*cl_base;
// Field[2].StopAtDistMax = 5.0*cl_base;

// Refinement of the surroundings of air water-transitions
Field[3] = Distance;
Field[3].CurvesList = {69,70,71,72,73,79,80,81,82,83};
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
MeshAlgorithm Surface {3,6,9,12} = 6; //5;  // RD surfaces
MeshAlgorithm Surface {2,5,8,11} = 6; //9;  // PML surfaces
MeshAlgorithm Surface {1,4,7,10} = 6; //5;  // Remaining
Mesh.Smoothing = 1;