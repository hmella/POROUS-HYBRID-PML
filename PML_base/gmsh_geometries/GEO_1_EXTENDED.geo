// Copyright (C) 2019 Hernan Mella

Include "dimensions";

// Get maximun characteristic length
cl = cl_base;

// Points
Point(1) = {-R_ext, 0, 0, cl};
Point(2) = {0, -R_ext, 0, cl};
Point(3) = {R_ext, 0, 0, cl};
Point(4) = {0, 0, 0, cl}; 
Point(5) = {xc+0.5*b, 0, 0, cl};
Point(6) = {xc, yc, 0, cl}; 
Point(7) = {xc-0.5*b, 0, 0, cl};
Point(8) = {0.5*Lx, 0, 0, cl};
Point(9) = {0.5*Lx, -Ly, 0, cl};
Point(10) = {-0.5*Lx, 0, 0, cl};
Point(11) = {-0.5*Lx, -Ly, 0, cl};
Point(12) = {0.5*Lx+LPML, 0, 0, cl};
Point(13) = {0.5*Lx+LPML, -Ly-LPML, 0, cl};
Point(14) = {-0.5*Lx-LPML, 0, 0, cl};
Point(15) = {-0.5*Lx-LPML, -Ly-LPML, 0, cl};
Point(16) = {0, 0, 0, cl};

// Generate interface lines to achieve a fixed resolution
// Add points to the interface to enforce the same interior mesh
// Left vertical line
nPoints1 = Ceil(Ly/cl); dx = Ly/nPoints1;
For i In {1 : nPoints1-1}
  pList1[i-1] = newp;
  Point(pList1[i-1]) = {-0.5*Lx, -i*dx, 0, cl};
EndFor

// Bottom horizontal line
nPoints2 = Ceil(Lx/cl); dx = Lx/nPoints2;
For i In {1 : nPoints2-1}
  pList2[i-1] = newp;
  Point(pList2[i-1]) = {-0.5*Lx+i*dx, -Ly, 0, cl};
EndFor

// Right vertical line
nPoints3 = Ceil(Ly/cl); dx = Ly/nPoints3;
For i In {1 : nPoints3-1}
  pList3[i-1] = newp;
  Point(pList3[i-1]) = {0.5*Lx, -Ly+i*dx, 0, cl};
EndFor


// Generate interface lines in the exterior PML layer
// to achieve a fixed resolution
// Add points to the interface to enforce the same interior mesh
// Left vertical line
nPoints4 = Ceil((Ly+LPML)/cl); dx = (Ly+LPML)/nPoints4;
For i In {1 : nPoints4-1}
  pList4[i-1] = newp;
  Point(pList4[i-1]) = {-0.5*Lx-LPML, -i*dx, 0, cl};
EndFor

// Bottom horizontal line
nPoints5 = Ceil((Lx+2*LPML)/cl); dx = (Lx+2*LPML)/nPoints5;
For i In {1 : nPoints5-1}
  pList5[i-1] = newp;
  Point(pList5[i-1]) = {-0.5*Lx-LPML+i*dx, -Ly-LPML, 0, cl};
EndFor

// Right vertical line
nPoints6 = Ceil((Ly+LPML)/cl); dx = (Ly+LPML)/nPoints6;
For i In {1 : nPoints6-1}
  pList6[i-1] = newp;
  Point(pList6[i-1]) = {0.5*Lx+LPML, -Ly-LPML+i*dx, 0, cl};
EndFor

// Lines (interface boundary)
Line(1) = {10, pList1[0]};
Spline(2) = pList1[];
Line(3) = {pList1[nPoints1-2],11};

Line(4) = {11,pList2[0]};
Spline(5) = pList2[];
Line(6) = {pList2[nPoints2-2],9};

Line(7) = {9,pList3[0]};
Spline(8) = pList3[];
Line(9) = {pList3[nPoints3-2],8};

// Lines (PML exterior boundary)
Line(10) = {14, pList4[0]};
Spline(11) = pList4[];
Line(12) = {pList4[nPoints4-2],15};

Line(13) = {15,pList5[0]};
Spline(14) = pList5[];
Line(15) = {pList5[nPoints5-2],13};

Line(16) = {13,pList6[0]};
Spline(17) = pList6[];
Line(18) = {pList6[nPoints6-2],12};

// Lines
Circle(19) = {1, 16, 2};
Circle(20) = {2, 16, 3};
Line(21) = {1, 14};
Line(22) = {14, 10};
Line(23) = {10, 7};
Line(24) = {7, 5};
Line(25) = {5, 8};
Line(26) = {8, 12};
Line(27) = {12, 3};

// Surfaces
Curve Loop(1) = {21, 10, 11, 12, 13, 14, 15, 16, 17, 18, 27, -20, -19};
Plane Surface(1) = {1};
Curve Loop(2) = {22, 1, 2, 3, 4, 5, 6, 7, 8, 9, 26, -18, -17, -16, -15, -14, -13, -12, -11, -10};
Plane Surface(2) = {2};
Curve Loop(3) = {23, 24, 25, -9, -8, -7, -6, -5, -4, -3, -2, -1};
Plane Surface(3) = {3};

// Physical surfaces
If (flag==0) // Extended domain
  Physical Surface(1) = {1};
  Physical Surface(2) = {2};
  Physical Surface(3) = {3};
ElseIf (flag==1) // PML domain
  Physical Surface(1) = {2}; // PML
  Physical Surface(2) = {3}; // RD
ElseIf (flag==2) // Paraxial domain
  Physical Surface(1) = {3};
EndIf

// Physical lines
Physical Line(3) = {24};                  // Neumann boundary (load)
If (flag==0) // Extended domain
  Physical Line(1) = {19,20};             // Dirichlet boundary
  Physical Line(4) = {21,22,23,25,26,27}; // Neumann boundary (remaining)
ElseIf (flag==1) // PML domain
  Physical Line(1) = {10,11,12,13,14,15,16,17,18}; // Dirichlet boundary
  Physical Line(2) = {1,2,3,4,5,6,7,8,9};          // Interface boundary  
  Physical Line(4) = {22,23,25,26};    // Neumann boundary (remaining)
ElseIf (flag==2) // Paraxial domain
  Physical Line(1) = {1,2,3,4,5,6,7,8,9}; // Interface boundary  
  Physical Line(4) = {23,25};             // Neumann boundary (remaining)
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
Field[1]. PointsList = {6};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = cl_source;
Field[2].SizeMax = cl_base;
Field[2].DistMin = 5.0*cl_source;
Field[2].DistMax = 5.0*cl_base;
Field[2].StopAtDistMax = 5.0*cl_base;

// Coarsening
Field[3] = Distance;
Field[3].CurvesList = {19,20};
Field[3].NumPointsPerCurve = 1000;

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = cl_coarse;
Field[4].SizeMax = cl_base;
Field[4].DistMin = 0.9*(R_ext-R_int);
Field[4].DistMax = R_ext-R_int;
// Field[4].StopAtDistMax = (R_ext-R_int);

// Define global mesh size
Field[7] = Min;
Field[7].FieldsList = {2, 4};
Background Field = 7;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms
MeshAlgorithm Surface {3} = 6; //5; // RD surfaces
MeshAlgorithm Surface {2} = 6; //9; // PML surfaces
MeshAlgorithm Surface {1} = 6; //5; // Remaining
Mesh.Smoothing = 1;