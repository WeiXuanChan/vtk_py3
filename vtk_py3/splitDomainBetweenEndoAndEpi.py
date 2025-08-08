########################################################################

import sys
import vtk

from .mat_vec_tools import *

########################################################################

def splitDomainBetweenEndoAndEpi(domain, verbose=True,RV=False):

    if (verbose): print ('*** splitDomainBetweenEndoAndEpi ***')

    connectivity0 = vtk.vtkConnectivityFilter()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        connectivity0.SetInputData(domain)
    else:
       connectivity0.SetInput(domain)
    connectivity0.SetExtractionModeToSpecifiedRegions()
    connectivity0.AddSpecifiedRegion(0)
    connectivity0.Update()
    ugrid0_temp = connectivity0.GetOutput()

    geom0 = vtk.vtkGeometryFilter()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
    	geom0.SetInputData(ugrid0_temp)
    else:
    	geom0.SetInput(ugrid0_temp)
    geom0.Update()
    pdata0_temp = geom0.GetOutput()

    tfilter0 = vtk.vtkTriangleFilter()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
    	tfilter0.SetInputData(pdata0_temp)
    else:
    	tfilter0.SetInput(pdata0_temp)
    tfilter0.Update()

    connectivity1 = vtk.vtkConnectivityFilter()
    connectivity1.SetExtractionModeToSpecifiedRegions()
    connectivity1.AddSpecifiedRegion(1)
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        connectivity1.SetInputData(domain)
    else:
        connectivity1.SetInput(domain)
    connectivity1.Update()
    ugrid1_temp = connectivity1.GetOutput()
    geom1 = vtk.vtkGeometryFilter()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
    	geom1.SetInputData(ugrid1_temp)
    else:
    	geom1.SetInput(ugrid1_temp)
    geom1.Update()
    pdata1_temp = geom1.GetOutput()


    tfilter1 = vtk.vtkTriangleFilter()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
    	tfilter1.SetInputData(pdata1_temp)
    else:
    	tfilter1.SetInput(pdata1_temp)
    tfilter1.Update()

    pdata1 = tfilter1.GetOutput()
    pdata0 = tfilter0.GetOutput()    
    
    p0bd = pdata0.GetBounds()
    p1bd = pdata1.GetBounds()
    
    if RV:
        connectivity2 = vtk.vtkConnectivityFilter()
        connectivity2.SetExtractionModeToSpecifiedRegions()
        connectivity2.AddSpecifiedRegion(2)
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            connectivity2.SetInputData(domain)
        else:
            connectivity2.SetInput(domain)
        connectivity2.Update()
        ugrid2_temp = connectivity2.GetOutput()
        geom2 = vtk.vtkGeometryFilter()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        	geom2.SetInputData(ugrid2_temp)
        else:
        	geom2.SetInput(ugrid2_temp)
        geom2.Update()
        pdata2_temp = geom2.GetOutput()
    
    
        tfilter2 = vtk.vtkTriangleFilter()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        	tfilter2.SetInputData(pdata2_temp)
        else:
        	tfilter2.SetInput(pdata2_temp)
        tfilter2.Update()
    
        pdata2 = tfilter2.GetOutput()
        p2bd = pdata2.GetBounds()
        
        if (abs(p0bd[0] - p0bd[1]))>=(abs(p1bd[0] - p1bd[1])) and (abs(p0bd[0] - p0bd[1]))>=(abs(p2bd[0] - p2bd[1])):
            pdata_epi = pdata0
            if p1bd[1] > p2bd[1]:
                pdata_endo1 = pdata2
                pdata_endo2 = pdata1
            else:
                pdata_endo1 = pdata1
                pdata_endo2 = pdata2
        elif (abs(p1bd[0] - p1bd[1]))>(abs(p2bd[0] - p2bd[1])) and (abs(p1bd[0] - p1bd[1]))>(abs(p0bd[0] - p0bd[1])):
            pdata_epi = pdata1
            if p0bd[1] > p2bd[1]:
                pdata_endo1 = pdata2
                pdata_endo2 = pdata0
            else:
                pdata_endo1 = pdata2
                pdata_endo2 = pdata0
        elif (abs(p2bd[0] - p2bd[1]))>(abs(p0bd[0] - p0bd[1])) and (abs(p2bd[0] - p2bd[1]))>(abs(p1bd[0] - p1bd[1])):
            pdata_epi = pdata2
            if p0bd[1] > p1bd[1]:
                pdata_endo1 = pdata1
                pdata_endo2 = pdata0
            else:
                pdata_endo1 = pdata1
                pdata_endo2 = pdata0
        return pdata_epi, pdata_endo1 ,pdata_endo2
    else:
        if (abs(p1bd[0] - p1bd[1]) < abs(p0bd[0] - p0bd[1])):
            pdata_epi = pdata0
            pdata_endo = pdata1
        else:
            pdata_epi = pdata1
            pdata_endo = pdata0
    
        return pdata_epi, pdata_endo



