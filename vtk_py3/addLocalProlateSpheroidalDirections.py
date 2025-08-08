########################################################################

import sys
import math
import numpy
import vtk

from .createFloatArray               import *
from .getABPointsFromBoundsAndCenter import *
from .getCellCenters                 import *
from .getPDataNormals                import *
from .writePData			    import *

########################################################################

def addLocalProlateSpheroidalDirections(ugrid_wall,
                                        pdata_end,
                                        pdata_epi,
                                        type_of_support="cell",
                                        points_AB=None,
					eCCname="eCC",
					eLLname="eLL",
					eRRname="eRR",
                                        verbose=True):

    if (verbose): print ('*** addLocalProlateSpheroidalDirections ***')
    if (verbose): print ("Initializing cell locators...")

    cell_locator_end = vtk.vtkCellLocator()
    cell_locator_end.SetDataSet(pdata_end)
    cell_locator_end.Update()

    cell_locator_epi = vtk.vtkCellLocator()
    cell_locator_epi.SetDataSet(pdata_epi)
    cell_locator_epi.Update()

    closest_point_end = [0.]*3
    closest_point_epi = [0.]*3
    generic_cell = vtk.vtkGenericCell()
    cellId_end = vtk.mutable(0)
    cellId_epi = vtk.mutable(0)
    subId = vtk.mutable(0)
    dist_end = vtk.mutable(0.)
    dist_epi = vtk.mutable(0.)
    
    
    if (points_AB == None):
        points_AB = getABPointsFromBoundsAndCenter(pdata_epi, verbose)
    assert (points_AB.GetNumberOfPoints() == 2), "points_AB must have two points. Aborting."
    point_A = numpy.array([0.]*3)
    point_B = numpy.array([0.]*3)
    points_AB.GetPoint(0, point_A)
    points_AB.GetPoint(1, point_B)
    print("DISPLAY point_A,point_B",point_A,point_B)
    eL  = point_B - point_A #isapexflip=Flase
    eL /= numpy.linalg.norm(eL)

    if (type_of_support == "cell"):
        pdata_cell_centers = getCellCenters(ugrid_wall)

    if (verbose): print ("Computing cell normals...")
    
    pdata_epi_temp = getPDataNormals(pdata_epi, flip=1)
    pdata_epi_centers = getCellCenters(pdata_epi_temp)
    pdata_pointoutwards_sum=0.
    for num_cell in range(pdata_epi_centers.GetPoints().GetNumberOfPoints()):
        cell_center = numpy.array(pdata_epi_centers.GetPoints().GetPoint(num_cell))
        cell_locator_epi.FindClosestPoint(cell_center, closest_point_epi, generic_cell, cellId_epi, subId, dist_epi)
        normal_epi = numpy.reshape(pdata_epi_temp.GetCellData().GetNormals().GetTuple(cellId_epi), (3))
        if cell_center[2]>(point_A[2]/2.):
            pdata_pointoutwards_sum+=numpy.sum(cell_center[:2]*normal_epi[:2])/numpy.linalg.norm(cell_center[:2])/numpy.linalg.norm(normal_epi[:2])
    if(pdata_pointoutwards_sum<0):
    	pdata_epi = getPDataNormals(pdata_epi, flip=0)
    else:
    	pdata_epi = getPDataNormals(pdata_epi, flip=1)
    
    pdata_end_temp = getPDataNormals(pdata_end, flip=0)
    pdata_end_centers = getCellCenters(pdata_end_temp)
    pdata_pointoutwards_sum=0.
    for num_cell in range(pdata_end_centers.GetPoints().GetNumberOfPoints()):
        cell_center = numpy.array(pdata_end_centers.GetPoints().GetPoint(num_cell))
        cell_locator_end.FindClosestPoint(cell_center, closest_point_end, generic_cell, cellId_end, subId, dist_end)
        normal_end = numpy.reshape(pdata_end_temp.GetCellData().GetNormals().GetTuple(cellId_end), (3))
        if cell_center[2]>(point_A[2]/2.):
            pdata_pointoutwards_sum+=numpy.sum(cell_center[:2]*normal_end[:2])/numpy.linalg.norm(cell_center[:2])/numpy.linalg.norm(normal_end[:2])
    if(pdata_pointoutwards_sum>0):
    	pdata_end = getPDataNormals(pdata_end, flip=1)
    else:
    	pdata_end = getPDataNormals(pdata_end, flip=0)
   
    if (verbose): print ("Computing surface bounds...")

    bounds_end = pdata_end.GetBounds()
    bounds_epi = pdata_epi.GetBounds()
    z_min_end = bounds_end[4]
    z_min_epi = bounds_epi[4]
    z_max_end = bounds_end[5]
    z_max_epi = bounds_epi[5]
    L_end = z_max_end-z_min_end
    L_epi = z_max_epi-z_min_epi

    

    if (verbose): print ("Computing local prolate spheroidal directions...")

    if (type_of_support == "cell"):
        nb_cells = ugrid_wall.GetNumberOfCells()
    elif (type_of_support == "point"):
        nb_cells = ugrid_wall.GetNumberOfPoints()

    farray_norm_dist_end = createFloatArray("norm_dist_end", 1, nb_cells)
    farray_norm_dist_epi = createFloatArray("norm_dist_epi", 1, nb_cells)

    farray_norm_z_end = createFloatArray("norm_z_end", 1, nb_cells)
    farray_norm_z_epi = createFloatArray("norm_z_epi", 1, nb_cells)

    farray_eRR = createFloatArray(eRRname, 3, nb_cells)
    farray_eCC = createFloatArray(eCCname, 3, nb_cells)
    farray_eLL = createFloatArray(eLLname, 3, nb_cells)

    for num_cell in range(nb_cells):
        if (type_of_support == "cell"):
            cell_center = numpy.array(pdata_cell_centers.GetPoints().GetPoint(num_cell))
        elif (type_of_support == "point"):
            cell_center = numpy.array(ugrid_wall.GetPoints().GetPoint(num_cell))
        cell_locator_end.FindClosestPoint(cell_center, closest_point_end, generic_cell, cellId_end, subId, dist_end)
        cell_locator_epi.FindClosestPoint(cell_center, closest_point_epi, generic_cell, cellId_epi, subId, dist_epi)

        norm_dist_end = dist_end/(dist_end+dist_epi)
        norm_dist_epi = dist_epi/(dist_end+dist_epi)
        farray_norm_dist_end.InsertTuple(num_cell, [norm_dist_end])
        farray_norm_dist_epi.InsertTuple(num_cell, [norm_dist_epi])

        norm_z_end = (closest_point_end[2]-z_min_end)/L_end
        norm_z_epi = (closest_point_epi[2]-z_min_epi)/L_epi
        farray_norm_z_end.InsertTuple(num_cell, [norm_z_end])
        farray_norm_z_epi.InsertTuple(num_cell, [norm_z_epi])

        normal_end = numpy.reshape(pdata_end.GetCellData().GetNormals().GetTuple(cellId_end), (3))
        normal_epi = numpy.reshape(pdata_epi.GetCellData().GetNormals().GetTuple(cellId_epi), (3))
        eRR  = -1*(1.-norm_dist_end) * normal_end + (1.-norm_dist_epi) * normal_epi
        eRR /= numpy.linalg.norm(eRR)
        eCC  = numpy.cross(eL, eRR)
        eCC /= numpy.linalg.norm(eCC)
        eLL  = numpy.cross(eRR, eCC)
        farray_eRR.InsertTuple(num_cell, eRR)
        farray_eCC.InsertTuple(num_cell, eCC)
        farray_eLL.InsertTuple(num_cell, eLL)

    if (verbose): print ("Filling mesh...")

    if (type_of_support == "cell"):
        ugrid_wall.GetCellData().AddArray(farray_norm_dist_end)
        ugrid_wall.GetCellData().AddArray(farray_norm_dist_epi)
        ugrid_wall.GetCellData().AddArray(farray_norm_z_end)
        ugrid_wall.GetCellData().AddArray(farray_norm_z_epi)
        ugrid_wall.GetCellData().AddArray(farray_eRR)
        ugrid_wall.GetCellData().AddArray(farray_eCC)
        ugrid_wall.GetCellData().AddArray(farray_eLL)
    elif (type_of_support == "point"):
        ugrid_wall.GetPointData().AddArray(farray_norm_dist_end)
        ugrid_wall.GetPointData().AddArray(farray_norm_dist_epi)
        ugrid_wall.GetPointData().AddArray(farray_norm_z_end)
        ugrid_wall.GetPointData().AddArray(farray_norm_z_epi)
        ugrid_wall.GetPointData().AddArray(farray_eRR)
        ugrid_wall.GetPointData().AddArray(farray_eCC)
        ugrid_wall.GetPointData().AddArray(farray_eLL)

def addLocalProlateSpheroidalDirections_AVsplit(ugrid_wall,
                                        pdata_endAtr,
                                        pdata_endVen,
                                        pdata_epiAtr,
                                        pdata_epiVen,
                                        pdata_inlet,
                                        pdata_avj,
                                        pdata_outlet,
                                        type_of_support="cell",
					eCCname="eCC",
					eLLname="eLL",
					eRRname="eRR",
                    atrial_delay=0.,
                                        verbose=True):

    if (verbose): print ('*** addLocalProlateSpheroidalDirections ***')
    if (verbose): print ("Initializing cell locators...")

    cell_locator_endAtr = vtk.vtkCellLocator()
    cell_locator_endAtr.SetDataSet(pdata_endAtr)
    cell_locator_endAtr.Update()
    
    cell_locator_endVen = vtk.vtkCellLocator()
    cell_locator_endVen.SetDataSet(pdata_endVen)
    cell_locator_endVen.Update()

    cell_locator_epiAtr = vtk.vtkCellLocator()
    cell_locator_epiAtr.SetDataSet(pdata_epiAtr)
    cell_locator_epiAtr.Update()
    
    cell_locator_epiVen = vtk.vtkCellLocator()
    cell_locator_epiVen.SetDataSet(pdata_epiVen)
    cell_locator_epiVen.Update()

    closest_point_endAtr = [0.]*3
    closest_point_epiAtr = [0.]*3
    closest_point_endVen = [0.]*3
    closest_point_epiVen = [0.]*3
    generic_cell = vtk.vtkGenericCell()
    cellId_endAtr = vtk.mutable(0)
    cellId_epiAtr = vtk.mutable(0)
    cellId_endVen = vtk.mutable(0)
    cellId_epiVen = vtk.mutable(0)
    subId = vtk.mutable(0)
    dist_endAtr = vtk.mutable(0.)
    dist_epiAtr = vtk.mutable(0.)
    dist_endVen = vtk.mutable(0.)
    dist_epiVen = vtk.mutable(0.)
    
    pointsAtr_AB = getABPointsFromTwoCenters(pdata_inlet,pdata_avj, verbose)
    pointsVen_AB = getABPointsFromTwoCenters(pdata_avj,pdata_outlet, verbose)
    assert (pointsAtr_AB.GetNumberOfPoints() == 2), "points_AB must have two points. Aborting."
    assert (pointsVen_AB.GetNumberOfPoints() == 2), "points_AB must have two points. Aborting."

    

    if (verbose): print ("Computing local prolate spheroidal directions...")

    if (type_of_support == "cell"):
        nb_cells = ugrid_wall.GetNumberOfCells()
    elif (type_of_support == "point"):
        nb_cells = ugrid_wall.GetNumberOfPoints()

    farray_norm_dist_end = createFloatArray("norm_dist_end", 1, nb_cells)
    farray_norm_dist_epi = createFloatArray("norm_dist_epi", 1, nb_cells)

    farray_norm_z_end = createFloatArray("norm_z_end", 1, nb_cells)
    farray_norm_z_epi = createFloatArray("norm_z_epi", 1, nb_cells)

    farray_eRR = createFloatArray(eRRname, 3, nb_cells)
    farray_eCC = createFloatArray(eCCname, 3, nb_cells)
    farray_eLL = createFloatArray(eLLname, 3, nb_cells)
    
    pointAtr_A = numpy.array([0.]*3)
    pointAtr_B = numpy.array([0.]*3)
    pointsAtr_AB.GetPoint(0, pointAtr_A)
    pointsAtr_AB.GetPoint(1, pointAtr_B)
    print("DISPLAY Atrial point_A,point_B",pointAtr_A,pointAtr_B)
    eL_Atr  = pointAtr_B - pointAtr_A #isapexflip=Flase
    eL_Atr /= numpy.linalg.norm(eL_Atr)
    
    pointVen_A = numpy.array([0.]*3)
    pointVen_B = numpy.array([0.]*3)
    pointsVen_AB.GetPoint(0, pointVen_A)
    pointsVen_AB.GetPoint(1, pointVen_B)
    print("DISPLAY Ventricle point_A,point_B",pointVen_A,pointVen_B)
    eL_Ven  = pointVen_B - pointVen_A #isapexflip=Flase
    eL_Ven /= numpy.linalg.norm(eL_Ven)
    
    eL_Z=numpy.array([0.,0.,-1.])

    if (type_of_support == "cell"):
        pdata_cell_centers = getCellCenters(ugrid_wall)

    if (verbose): print ("Computing cell normals...")
    
    pdata_epiAtr_temp = getPDataNormals(pdata_epiAtr, flip=1)
    pdata_epiAtr_centers = getCellCenters(pdata_epiAtr_temp)
    pdata_pointoutwards_Atrsum=0.
    for num_cell in range(pdata_epiAtr_centers.GetPoints().GetNumberOfPoints()):
        cell_center = numpy.array(pdata_epiAtr_centers.GetPoints().GetPoint(num_cell))
        cell_locator_epiAtr.FindClosestPoint(cell_center, closest_point_epiAtr, generic_cell, cellId_epiAtr, subId, dist_epiAtr)
        normal_epi = numpy.reshape(pdata_epiAtr_temp.GetCellData().GetNormals().GetTuple(cellId_epiAtr), (3))
        if cell_center[2]<(pointAtr_A[2]/2.):
            pdata_pointoutwards_Atrsum+=numpy.sum(cell_center[:2]*normal_epi[:2])/numpy.linalg.norm(cell_center[:2])/numpy.linalg.norm(normal_epi[:2])
    if(pdata_pointoutwards_Atrsum<0):
    	pdata_epiAtr = getPDataNormals(pdata_epiAtr, flip=0)
    else:
    	pdata_epiAtr = getPDataNormals(pdata_epiAtr, flip=1)
    
    pdata_epiVen_temp = getPDataNormals(pdata_epiVen, flip=1)
    pdata_epiVen_centers = getCellCenters(pdata_epiVen_temp)
    pdata_pointoutwards_Vensum=0.
    for num_cell in range(pdata_epiVen_centers.GetPoints().GetNumberOfPoints()):
        cell_center = numpy.array(pdata_epiVen_centers.GetPoints().GetPoint(num_cell))
        cell_locator_epiVen.FindClosestPoint(cell_center, closest_point_epiVen, generic_cell, cellId_epiVen, subId, dist_epiVen)
        normal_epi = numpy.reshape(pdata_epiVen_temp.GetCellData().GetNormals().GetTuple(cellId_epiVen), (3))
        if cell_center[2]>(pointVen_B[2]/2.):
            pdata_pointoutwards_Vensum+=numpy.sum(cell_center[:2]*normal_epi[:2])/numpy.linalg.norm(cell_center[:2])/numpy.linalg.norm(normal_epi[:2])
    if(pdata_pointoutwards_Vensum<0):
    	pdata_epiVen = getPDataNormals(pdata_epiVen, flip=0)
    else:
    	pdata_epiVen = getPDataNormals(pdata_epiVen, flip=1)
        
        
    pdata_endAtr_temp = getPDataNormals(pdata_endAtr, flip=0)
    pdata_endAtr_centers = getCellCenters(pdata_endAtr_temp)
    pdata_pointoutwards_Atrsum=0.
    for num_cell in range(pdata_endAtr_centers.GetPoints().GetNumberOfPoints()):
        cell_center = numpy.array(pdata_endAtr_centers.GetPoints().GetPoint(num_cell))
        cell_locator_endAtr.FindClosestPoint(cell_center, closest_point_endAtr, generic_cell, cellId_endAtr, subId, dist_endAtr)
        normal_end = numpy.reshape(pdata_endAtr_temp.GetCellData().GetNormals().GetTuple(cellId_endAtr), (3))
        if cell_center[2]<(pointAtr_A[2]/2.):
            pdata_pointoutwards_Atrsum+=numpy.sum(cell_center[:2]*normal_end[:2])/numpy.linalg.norm(cell_center[:2])/numpy.linalg.norm(normal_end[:2])
    if(pdata_pointoutwards_Atrsum>0):
    	pdata_endAtr = getPDataNormals(pdata_endAtr, flip=1)
    else:
    	pdata_endAtr = getPDataNormals(pdata_endAtr, flip=0)
        
    pdata_endVen_temp = getPDataNormals(pdata_endVen, flip=0)
    pdata_endVen_centers = getCellCenters(pdata_endVen_temp)
    pdata_pointoutwards_Vensum=0.
    for num_cell in range(pdata_endVen_centers.GetPoints().GetNumberOfPoints()):
        cell_center = numpy.array(pdata_endVen_centers.GetPoints().GetPoint(num_cell))
        cell_locator_endVen.FindClosestPoint(cell_center, closest_point_endVen, generic_cell, cellId_endVen, subId, dist_endVen)
        normal_end = numpy.reshape(pdata_endVen_temp.GetCellData().GetNormals().GetTuple(cellId_endVen), (3))
        if cell_center[2]>(pointVen_B[2]/2.):
            pdata_pointoutwards_Vensum+=numpy.sum(cell_center[:2]*normal_end[:2])/numpy.linalg.norm(cell_center[:2])/numpy.linalg.norm(normal_end[:2])
    if(pdata_pointoutwards_Vensum>0):
    	pdata_endVen = getPDataNormals(pdata_endVen, flip=1)
    else:
    	pdata_endVen = getPDataNormals(pdata_endVen, flip=0)
   
    if (verbose): print ("Computing surface bounds...")

    bounds_endAtr = pdata_endAtr.GetBounds()
    bounds_epiAtr = pdata_epiAtr.GetBounds()
    bounds_endVen = pdata_endVen.GetBounds()
    bounds_epiVen = pdata_epiVen.GetBounds()
    z_min_end = bounds_endVen[4]
    z_min_epi = bounds_epiVen[4]
    z_max_end = bounds_endAtr[5]
    z_max_epi = bounds_epiAtr[5]
    L_end = z_max_end-z_min_end
    L_epi = z_max_epi-z_min_epi
    
    pdata_inlet_temp = getPDataNormals(pdata_inlet, flip=0)
    pdata_inlet_normvec=numpy.array(pdata_inlet_temp.GetCellData().GetNormals().GetTuple(0))
    if numpy.sum(pdata_inlet_normvec*eL_Atr)<0:
        pdata_inlet_normvec=-pdata_inlet_normvec
    pdata_inlet_normvec /= numpy.linalg.norm(pdata_inlet_normvec)
    
    pdata_outlet_temp = getPDataNormals(pdata_outlet, flip=0)
    pdata_outlet_normvec=numpy.array(pdata_outlet_temp.GetCellData().GetNormals().GetTuple(0))
    if numpy.sum(pdata_outlet_normvec*eL_Ven)<0:
        pdata_outlet_normvec=-pdata_outlet_normvec
    pdata_outlet_normvec /= numpy.linalg.norm(pdata_outlet_normvec)
    
    for num_cell in range(nb_cells):
        if (type_of_support == "cell"):
            cell_center = numpy.array(pdata_cell_centers.GetPoints().GetPoint(num_cell))
        elif (type_of_support == "point"):
            cell_center = numpy.array(ugrid_wall.GetPoints().GetPoint(num_cell))
        cell_locator_endAtr.FindClosestPoint(cell_center, closest_point_endAtr, generic_cell, cellId_endAtr, subId, dist_endAtr)
        cell_locator_epiAtr.FindClosestPoint(cell_center, closest_point_epiAtr, generic_cell, cellId_epiAtr, subId, dist_epiAtr)
        cell_locator_endVen.FindClosestPoint(cell_center, closest_point_endVen, generic_cell, cellId_endVen, subId, dist_endVen)
        cell_locator_epiVen.FindClosestPoint(cell_center, closest_point_epiVen, generic_cell, cellId_epiVen, subId, dist_epiVen)
        
        if (dist_endAtr+dist_epiAtr)<(dist_endVen+dist_epiVen):
            norm_dist_end = dist_endAtr/(dist_endAtr+dist_epiAtr)
            norm_dist_epi = dist_epiAtr/(dist_endAtr+dist_epiAtr)
        else:
            norm_dist_end = dist_endVen/(dist_endVen+dist_epiVen)
            norm_dist_epi = dist_epiVen/(dist_endVen+dist_epiVen)
        farray_norm_dist_end.InsertTuple(num_cell, [norm_dist_end])
        farray_norm_dist_epi.InsertTuple(num_cell, [norm_dist_epi])
        
        if (dist_endAtr+dist_epiAtr)<(dist_endVen+dist_epiVen):
            norm_z_end = (closest_point_endAtr[2]-z_min_end)/L_end
            norm_z_epi = (closest_point_epiAtr[2]-z_min_epi)/L_epi
        else:
            norm_z_end = (closest_point_endVen[2]-z_min_end)/L_end
            norm_z_epi = (closest_point_epiVen[2]-z_min_epi)/L_epi
        farray_norm_z_end.InsertTuple(num_cell, [norm_z_end])
        farray_norm_z_epi.InsertTuple(num_cell, [norm_z_epi])
        
        if (dist_endAtr+dist_epiAtr)<(dist_endVen+dist_epiVen):
            normal_end = numpy.reshape(pdata_endAtr.GetCellData().GetNormals().GetTuple(cellId_endAtr), (3))
            normal_epi = numpy.reshape(pdata_epiAtr.GetCellData().GetNormals().GetTuple(cellId_epiAtr), (3))
        else:
            normal_end = numpy.reshape(pdata_endVen.GetCellData().GetNormals().GetTuple(cellId_endVen), (3))
            normal_epi = numpy.reshape(pdata_epiVen.GetCellData().GetNormals().GetTuple(cellId_epiVen), (3))
        eRR  = -1*(1.-norm_dist_end) * normal_end + (1.-norm_dist_epi) * normal_epi
        eRR /= numpy.linalg.norm(eRR)
        if (dist_endAtr+dist_epiAtr)<(dist_endVen+dist_epiVen):
            disttoA=numpy.abs(numpy.sum((cell_center-pointAtr_A)*pdata_inlet_normvec))
            disttoB=numpy.abs(numpy.sum((cell_center-pointAtr_B)*eL_Z))
            tol=numpy.sum((pointAtr_B-pointAtr_A)**2.)*1e-3
            if abs(disttoA-disttoB)<tol:
                eL_temp=pdata_inlet_normvec+eL_Z
            else:
                eL_temp=pdata_inlet_normvec*disttoB+eL_Z*disttoA
            eCC  = numpy.cross(eL_temp, eRR)
        else:
            disttoA=numpy.abs(numpy.sum((cell_center-pointVen_A)*eL_Z))
            disttoB=numpy.abs(numpy.sum((cell_center-pointVen_B)*pdata_outlet_normvec))
            tol=numpy.sum((pointVen_B-pointVen_A)**2.)*1e-3
            if abs(disttoA-disttoB)<tol:
                eL_temp=pdata_outlet_normvec+eL_Z
            else:
                eL_temp=pdata_outlet_normvec*disttoA+eL_Z*disttoB
            eCC  = numpy.cross(eL_temp, eRR)
        eCC /= numpy.linalg.norm(eCC)
        eLL  = numpy.cross(eRR, eCC)
        farray_eRR.InsertTuple(num_cell, eRR)
        farray_eCC.InsertTuple(num_cell, eCC)
        farray_eLL.InsertTuple(num_cell, eLL)

    if (verbose): print ("Filling mesh...")

    if (type_of_support == "cell"):
        ugrid_wall.GetCellData().AddArray(farray_norm_dist_end)
        ugrid_wall.GetCellData().AddArray(farray_norm_dist_epi)
        ugrid_wall.GetCellData().AddArray(farray_norm_z_end)
        ugrid_wall.GetCellData().AddArray(farray_norm_z_epi)
        ugrid_wall.GetCellData().AddArray(farray_eRR)
        ugrid_wall.GetCellData().AddArray(farray_eCC)
        ugrid_wall.GetCellData().AddArray(farray_eLL)
    elif (type_of_support == "point"):
        ugrid_wall.GetPointData().AddArray(farray_norm_dist_end)
        ugrid_wall.GetPointData().AddArray(farray_norm_dist_epi)
        ugrid_wall.GetPointData().AddArray(farray_norm_z_end)
        ugrid_wall.GetPointData().AddArray(farray_norm_z_epi)
        ugrid_wall.GetPointData().AddArray(farray_eRR)
        ugrid_wall.GetPointData().AddArray(farray_eCC)
        ugrid_wall.GetPointData().AddArray(farray_eLL)

if (__name__ == "__main__"):
    assert (len(sys.argv) in [2]), "Number of arguments must be 1. Aborting."
    basename = sys.argv[1]
    ugrid_wall = readUGrid(basename + "-Mesh.vtk")
    pdata_end = readSTL(basename + "-End.stl")
    pdata_epi = readSTL(basename + "-Epi.stl")
    addLocalProlateSpheroidalDirections(ugrid_wall, pdata_end, pdata_epi)
    writeUGrid(ugrid_wall, basename + "-Mesh.vtk")
