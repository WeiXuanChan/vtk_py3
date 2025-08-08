'''
################################################
MIT License
Copyright (c) 2021 L. C. Lee
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
################################################
File: __init__.py
Description: vtk support modules 
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: L.C. LEE                   08MAR2019           - Created
  Author: w.x.chan@gmail.com         06Dec2022           - v4.1.0
                                                            - added activeFiberStrength
                                                            - added activeFiberDelay
  Author: w.x.chan@gmail.com         07Jun2024           - v4.4.0
                                                            - added geometry: AVsplit to extractFeNiCsBiVFacet
                                                            - added createAVmesh
'''
_version='4.4.0'

import os as os
import sys

from .addFieldPrincipalDirections         import *
from .addLocalCylindricalDirections       import *
from .addLocalFiberOrientation            import *
from .addLocalFiberOrientation2           import *
from .addLocalProlateSpheroidalDirections import *
from .addMappingFromPointsToCells         import *
from .addSystolicStrains                  import *
from .clipDomainForCutLVMesh              import *
from .clipSurfacesForCutLVMesh            import *
from .clipSurfacesForFullLVMesh           import *
from .createFloatArray                    import *
from .createIntArray                      import *
from .exportDynaDeformationGradients      import *
from .findPointsInCell                    import *
from .getABPointsFromBoundsAndCenter      import *
from .getABPointsFromTTTSectors           import *
from .getCellCenters                      import *
from .mapCellDataToCellData               import *
from .mapPointDataToCellData              import *
from .readAbaqusDeformationGradients      import *
from .readAbaqusMesh                      import *
from .readAbaqusStresses                  import *
from .readDynaDeformationGradients        import *
from .readDynaDeformationGradients_vtk    import *
from .readDynaMesh                        import *
from .readDynaMeshStructured              import *
from .readFiberOrientation                import *
from .readPData                           import *
from .readSTL                             import *
from .readUGrid                           import *
from .readXMLImage                        import *
from .readXMLPData                        import *
from .readXMLUGrid                        import *
from .readXMLPUGrid                       import *
from .rotateSymmetricMatrix               import *
from .splitDomainBetweenEndoAndEpi        import *
from .writeFiberOrientationFileForAbaqus  import *
from .writeFiberOrientationFileForGNUPlot import *
from .writePData                          import *
from .writeSTL                            import *
from .writeUGrid                          import *
from .writeXMLImage                       import *
from .writeXMLPData                       import *
from .writeXMLUGrid                       import *
from .transform_mesh_w_4x4mat             import *
from .transform_pdata                     import *
from .scale_pdata                         import *
from .sliceheart				 import *
from .clean_pdata			 import *
from .extract_slice_thickness		 import *
from .clipheart				 import *
from .clipSurfacesForCutLVMesh_PLY        import *
from .readPLY                             import *
from .extract_slice_thickness_LVonly	 import *
from .computeVolume	 		 import *
from .computeRegionsForBiV 		 import *
from .getCellLocator			 import *
from .CreateVertexFromPoint		 import *
from .convertUGridtoPdata		 import *
from .getcentroid			 import *
from .createLVmesh		         import *
from .createAVmesh		         import *
from .computeLVthickness			 import *
from .translate_mesh			 import *
#from .translate_pdata 			 import *
from .convertCellDataToXML		 import *
from .readMSHGrid                     	 import *
from .addLocalFiberOrientation_infarct    import *
from .convertXMLMeshToUGrid		 import *
from .extractFeNiCsBiVFacet               import *
#from .SetBiVFiber			 import *
#from .SetBiVFiber_Quad			 import *
from .SetBiVFiber_Quad_PyQ                import *
from .extractCellFromPData 		 import *
from .readVTI				 import *
from .create_ellipsoidal_LV		 import *
from .addLVfiber				 import *
from .addBiVfiber				 import *
from .convertQuadDataToVTK		 import *
from .convertUGridToXMLMesh		 import *
from .Set4ChamberDirection		 import *
from .create_BiVmesh			 import *
from .rotateUGrid			 import *
#from .convertXMLMeshToUGrid2D             import *
#from .SetBiVFiber_Ce			 import *
from .SetBiVFiber_J 			 import *
#from .SlerpTestJ 			 import *
from .rotatePData_w_axis			 import *
from .convertPDataToXMLMesh		 import *
from .extractUGridBasedOnThreshold	 import *
