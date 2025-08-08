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

wall = newsl; Surface Loop(wall) = {ss};

V_wall = newv; Volume(V_wall) = {wall };
