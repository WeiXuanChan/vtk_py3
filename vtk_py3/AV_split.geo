//Geometry.HideCompounds = 1;

Mesh.CharacteristicLengthFactor = <<mesh_d>>;  // Coarse

Mesh.Algorithm    = 5; // (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) (Default=2)
Mesh.RecombineAll = 0;

//Mesh.RemeshAlgorithm = 1; // (0=no split, 1=automatic, 2=automatic only with metis) (Default=0)

//Mesh.RemeshParametrization = 7; // (0=harmonic_circle, 1=conformal_spectral, 2=rbf, 3=harmonic_plane, 4=convex_circle, 5=convex_plane, 6=harmonic square, 7=conformal_fe) (Default=4)

Mesh.Algorithm3D    = 4; // (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree) (Default=1)
Mesh.Recombine3DAll = 0;

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh.Smoothing = 0;

Merge "<<Endofilename>>";
endo_atr_id = <<endo_atr_id>>;
endo_ven_id = <<endo_ven_id>>;
epi_atr_id = <<epi_atr_id>>;
epi_ven_id = <<epi_ven_id>>;
avj_id = <<avj_id>>;
Coherence Mesh;

DefineConstant[
  // Angle between two triangles above which an edge is considered as sharp
  angle = {40, Min 20, Max 120, Step 1,
    Name "Parameters/Angle for surface detection"},
  // For complex geometries, patches can be too complex, too elongated or too
  // large to be parametrized; setting the following option will force the
  // creation of patches that are amenable to reparametrization:
  forceParametrizablePatches = {0, Choices{0,1},
    Name "Parameters/Create surfaces guaranteed to be parametrizable"},
  // For open surfaces include the boundary edges in the classification process:
  includeBoundary = 1,
  // Force curves to be split on given angle:
  curveAngle = 180
];
ClassifySurfaces{angle * Pi/180, includeBoundary, forceParametrizablePatches,
                 curveAngle * Pi / 180};
CreateGeometry;

CreateTopology;

ss[] = Surface "*";
endo_atr_bnd() = Abs(Boundary{ Surface {ss[endo_atr_id]}; });
endo_ven_bnd() = Abs(Boundary{ Surface {ss[endo_ven_id]}; });
epi_atr_bnd() = Abs(Boundary{ Surface {ss[epi_atr_id]}; });
epi_ven_bnd() = Abs(Boundary{ Surface {ss[epi_ven_id]}; });

bbox_endo_atr_bnd0() = BoundingBox Curve { endo_atr_bnd(0) };
bbox_endo_atr_bnd1() = BoundingBox Curve { endo_atr_bnd(1) };
bbox_endo_ven_bnd0() = BoundingBox Curve { endo_ven_bnd(0) };
bbox_endo_ven_bnd1() = BoundingBox Curve { endo_ven_bnd(1) };
If ( Fabs(bbox_endo_atr_bnd0(1)-bbox_endo_ven_bnd0(1))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_endo_atr_bnd0(3)-bbox_endo_ven_bnd0(3))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_endo_atr_bnd0(5)-bbox_endo_ven_bnd0(5))<(0.01*Mesh.CharacteristicLengthFactor))
  endo_atr_bnd_use = 1;
  endo_ven_bnd_use = 1;
ElseIf ( Fabs(bbox_endo_atr_bnd0(1)-bbox_endo_ven_bnd1(1))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_endo_atr_bnd0(3)-bbox_endo_ven_bnd1(3))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_endo_atr_bnd0(5)-bbox_endo_ven_bnd1(5))<(0.01*Mesh.CharacteristicLengthFactor))
  endo_atr_bnd_use = 1;
  endo_ven_bnd_use = 0;
ElseIf ( Fabs(bbox_endo_atr_bnd1(1)-bbox_endo_ven_bnd0(1))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_endo_atr_bnd1(3)-bbox_endo_ven_bnd0(3))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_endo_atr_bnd1(5)-bbox_endo_ven_bnd0(5))<(0.01*Mesh.CharacteristicLengthFactor))
  endo_atr_bnd_use = 0;
  endo_ven_bnd_use = 1;
Else
  endo_atr_bnd_use = 0;
  endo_ven_bnd_use = 0;
EndIf

bbox_epi_atr_bnd0() = BoundingBox Curve { epi_atr_bnd(0) };
bbox_epi_atr_bnd1() = BoundingBox Curve { epi_atr_bnd(1) };
bbox_epi_ven_bnd0() = BoundingBox Curve { epi_ven_bnd(0) };
bbox_epi_ven_bnd1() = BoundingBox Curve { epi_ven_bnd(1) };
If ( Fabs(bbox_epi_atr_bnd0(1)-bbox_epi_ven_bnd0(1))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_epi_atr_bnd0(3)-bbox_epi_ven_bnd0(3))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_epi_atr_bnd0(5)-bbox_epi_ven_bnd0(5))<(0.01*Mesh.CharacteristicLengthFactor))
  epi_atr_bnd_use = 1;
  epi_ven_bnd_use = 1;
ElseIf ( Fabs(bbox_epi_atr_bnd0(1)-bbox_epi_ven_bnd1(1))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_epi_atr_bnd0(3)-bbox_epi_ven_bnd1(3))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_epi_atr_bnd0(5)-bbox_epi_ven_bnd1(5))<(0.01*Mesh.CharacteristicLengthFactor))
  epi_atr_bnd_use = 1;
  epi_ven_bnd_use = 0;
ElseIf ( Fabs(bbox_epi_atr_bnd1(1)-bbox_epi_ven_bnd0(1))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_epi_atr_bnd1(3)-bbox_epi_ven_bnd0(3))<(0.01*Mesh.CharacteristicLengthFactor) && Fabs(bbox_epi_atr_bnd1(5)-bbox_epi_ven_bnd0(5))<(0.01*Mesh.CharacteristicLengthFactor))
  epi_atr_bnd_use = 0;
  epi_ven_bnd_use = 1;
Else
  epi_atr_bnd_use = 0;
  epi_ven_bnd_use = 0;
EndIf

L_endo_inlet = newll; Curve Loop(L_endo_inlet) = endo_atr_bnd(endo_atr_bnd_use);
L_epi_inlet = newll; Curve Loop(L_epi_inlet) = epi_atr_bnd(epi_atr_bnd_use);
L_endo_outlet = newll; Curve Loop(L_endo_outlet) = endo_ven_bnd(endo_ven_bnd_use);
L_epi_outlet = newll; Curve Loop(L_epi_outlet) = epi_ven_bnd(epi_ven_bnd_use);



S_inlet = news; Plane Surface(S_inlet) = { L_epi_inlet};
S_outlet = news; Plane Surface(S_outlet) = { L_epi_outlet};


S_inlet_close = news; Plane Surface(S_inlet_close) = { L_endo_inlet };
S_outlet_close = news; Plane Surface(S_outlet_close) = { L_endo_outlet };

Physical Surface("endoatr") = {ss[endo_atr_id],S_inlet_close};
Physical Surface("endoven") = {ss[endo_ven_id],S_outlet_close};
Physical Surface("epiatr") = {ss[epi_atr_id]};
Physical Surface("epiven") = {ss[epi_ven_id]};
Physical Surface("avj") = {ss[avj_id]};

Physical Surface("INLET") = {S_inlet };
Physical Surface("OUTLET") = {S_outlet};

SL_atr_wall = newsl; Surface Loop(SL_atr_wall) = {S_inlet_close,ss[endo_atr_id],ss[avj_id],ss[epi_atr_id],S_inlet };
SL_ven_wall = newsl; Surface Loop(SL_ven_wall) = {S_outlet_close,ss[endo_ven_id],ss[avj_id],ss[epi_ven_id],S_outlet};

V_atr_wall = newv; Volume(V_atr_wall) = {SL_atr_wall };
V_ven_wall = newv; Volume(V_ven_wall) = {SL_ven_wall };
Physical Volume("ATRIAL") = {V_atr_wall };
Physical Volume("VENTRICLE") = {V_ven_wall };

//BooleanFragments{ Volume{V_wall }; Delete; }{ Surface{ss[epi_atr_id]}; Delete; }

avj_bnd() = Abs(Boundary{ Surface {ss[avj_id]}; });
bbox_avj_bnd0() = BoundingBox Curve { avj_bnd(0) };
bbox_avj_bnd1() = BoundingBox Curve { avj_bnd(1) };
avj_bnd_use = 0;
If ( bbox_avj_bnd1(3) < bbox_avj_bnd0(3))
  avj_bnd_use = 1;
EndIf

Field[1] = Distance;
Field[1].CurvesList = {avj_bnd(avj_bnd_use)};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = <<mesh_d>> * <<rsizemin>>;
Field[2].SizeMax = <<mesh_d>> * <<rsizemax>>;
Field[2].DistMin = <<mesh_d>> * <<rdistmin>>;
Field[2].DistMax = <<mesh_d>> * <<rdistmax>>;
Background Field = 2;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;