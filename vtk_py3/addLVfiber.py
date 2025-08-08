from dolfin import *
import vtk as vtk
from . import *

def addLVfiber(mesh, S,V, casename, endo_angle, epi_angle,  casedir, clipheight=0.05,fiberSheetletAngle=0.,fiberSheetletWidth=0.,radialFiberAngle=0.,fiberLength=0.,cstrainFactor=1.,activeFiberStrength=1.,activeFiberDelay=0.,fiberDirectionVec=None):


	fiberV = Function(V)
	sheetV = Function(V)
	sheetnormV = Function(V)
	cV = Function(V)
	lV = Function(V)
	rV = Function(V)
	cstrainFactorS=Function(V)
	activefiberStrengthS = Function(V)
	activefiberDelayS = Function(V)

	ugrid=vtk_py.convertXMLMeshToUGrid(mesh)
	pdata = vtk_py.convertUGridtoPdata(ugrid)
	C = vtk_py.getcentroid(pdata)

    #for inverted laxis
	#ztop = pdata.GetBounds()[4]
	#C = [C[0], C[1], ztop+clipheight]
	#clippedheart = vtk_py.clipheart(pdata, C, [0,0,-1], True)

	#for laxis apex towards -ve z
	ztop = pdata.GetBounds()[5]
	C = [C[0], C[1], ztop-clipheight]
	clippedheart = vtk_py.clipheart(pdata, C, [0,0,1], True)

	epi, endo= vtk_py.splitDomainBetweenEndoAndEpi(clippedheart)

	cleanepipdata = vtk.vtkCleanPolyData()
	if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
			cleanepipdata.SetInputData(epi)
	else:
			cleanepipdata.SetInput(epi)
	cleanepipdata.Update()
	pdata_epi = cleanepipdata.GetOutput()


	cleanendopdata = vtk.vtkCleanPolyData()
	if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
			cleanendopdata.SetInputData(endo)
	else:
			cleanendopdata.SetInput(endo)
	cleanendopdata.Update()
	pdata_endo = cleanendopdata.GetOutput()

	L_epi = pdata_epi.GetBounds()[5]  -  pdata_epi.GetBounds()[4]
	L_endo = pdata_endo.GetBounds()[5] - pdata_endo.GetBounds()[4]



	if(L_endo > L_epi):
		temp = pdata_epi
		pdata_epi = pdata_endo
		pdata_endo = temp
		


	# Quad points
	gdim = mesh.geometry().dim()
	xdofmap = V.sub(0).dofmap().dofs()
	ydofmap = V.sub(1).dofmap().dofs()
	zdofmap = V.sub(2).dofmap().dofs()
	sdofmap = S.dofmap().dofs()

	xq = V.tabulate_dof_coordinates().reshape((-1, gdim))
	xq0 = xq[xdofmap]  
	#if(dolfin.dolfin_version() == '1.6.0'):
	#	xq = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, gdim))
	#	xq0 = xq[xdofmap]  

	# Create an unstructured grid of Gauss Points
	points = vtk.vtkPoints()
	vertices = vtk.vtkCellArray()
	ugrid = vtk.vtkUnstructuredGrid()
	cnt = 0;
	for pt in xq0:
		points.InsertNextPoint([pt[0], pt[1], pt[2]])
		vertex = vtk.vtkVertex()
		vertex.GetPointIds().SetId(0, cnt)
		vertices.InsertNextCell(vertex)
		cnt += 1    

	ugrid.SetPoints(points)
	ugrid.SetCells(0, vertices)
    
	CreateVertexFromPoint(ugrid)
	addLocalProlateSpheroidalDirections(ugrid, pdata_endo, pdata_epi, type_of_support="cell")
	addLocalFiberOrientation(ugrid, endo_angle, epi_angle,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,cstrainFactor=cstrainFactor,activeFiberStrength=activeFiberStrength,activeFiberDelay=activeFiberDelay,fiberDirectionVec=fiberDirectionVec, type_of_support="cell")

	fiber_vector =  ugrid.GetCellData().GetArray("fiber vectors")
	sheet_vector =  ugrid.GetCellData().GetArray("sheet vectors")
	sheetnorm_vector =  ugrid.GetCellData().GetArray("sheet normal vectors")
	
	eCC_vector =  ugrid.GetCellData().GetArray("eCC")
	eLL_vector =  ugrid.GetCellData().GetArray("eLL")
	eRR_vector =  ugrid.GetCellData().GetArray("eRR")
    
	csf_vector =  ugrid.GetCellData().GetArray("cstrain_factor")
	afs_vector =  ugrid.GetCellData().GetArray("fiber_active_fiber_strength")
	afd_vector =  ugrid.GetCellData().GetArray("fiber_active_fiber_delay")

 
	cnt = 0
	for pt in xq0:

		fvec = fiber_vector.GetTuple(cnt)
		svec = sheet_vector.GetTuple(cnt)
		nvec = sheetnorm_vector.GetTuple(cnt)

		cvec = eCC_vector.GetTuple(cnt)
		lvec = eLL_vector.GetTuple(cnt)
		rvec = eRR_vector.GetTuple(cnt)

		csf = csf_vector.GetTuple(cnt)
		afs = afs_vector.GetTuple(cnt)
		afd = afd_vector.GetTuple(cnt)
        
		fvecnorm = sqrt(fvec[0]**2 + fvec[1]**2 + fvec[2]**2)
		svecnorm = sqrt(svec[0]**2 + svec[1]**2 + svec[2]**2)
		nvecnorm = sqrt(nvec[0]**2 + nvec[1]**2 + nvec[2]**2)


		if(abs(fvecnorm - 1.0) > 1e-7 or  abs(svecnorm - 1.0) > 1e-6 or abs(nvecnorm - 1.0) > 1e-7):
			print (fvecnorm)
			print (svecnorm)
			print (nvecnorm)

		#print xdofmap[cnt], ydofmap[cnt], zdofmap[cnt]
		fiberV.vector()[xdofmap[cnt]] = fvec[0]; fiberV.vector()[ydofmap[cnt]] = fvec[1]; fiberV.vector()[zdofmap[cnt]] = fvec[2];
		sheetV.vector()[xdofmap[cnt]] = svec[0]; sheetV.vector()[ydofmap[cnt]] = svec[1]; sheetV.vector()[zdofmap[cnt]] = svec[2];
		sheetnormV.vector()[xdofmap[cnt]] = nvec[0]; sheetnormV.vector()[ydofmap[cnt]] = nvec[1]; sheetnormV.vector()[zdofmap[cnt]] = nvec[2];

		cV.vector()[xdofmap[cnt]] = cvec[0];  cV.vector()[ydofmap[cnt]] = cvec[1]; cV.vector()[zdofmap[cnt]] = cvec[2]; 
		lV.vector()[xdofmap[cnt]] = lvec[0];  lV.vector()[ydofmap[cnt]] = lvec[1]; lV.vector()[zdofmap[cnt]] = lvec[2]; 
		rV.vector()[xdofmap[cnt]] = rvec[0];  rV.vector()[ydofmap[cnt]] = rvec[1]; rV.vector()[zdofmap[cnt]] = rvec[2]; 
        
		cstrainFactorS.vector()[xdofmap[cnt]] = csf[0]; cstrainFactorS.vector()[ydofmap[cnt]] = csf[1]; cstrainFactorS.vector()[zdofmap[cnt]] = csf[2];
		activefiberStrengthS.vector()[xdofmap[cnt]] = afs[0]; activefiberStrengthS.vector()[ydofmap[cnt]] = afs[1]; activefiberStrengthS.vector()[zdofmap[cnt]] = afs[2];
		activefiberDelayS.vector()[xdofmap[cnt]] = afd[0]; activefiberDelayS.vector()[ydofmap[cnt]] = afd[1]; activefiberDelayS.vector()[zdofmap[cnt]] = afd[2];
		#cstrainFactorS.vector()[sdofmap[cnt]] = csf[0]
		#activefiberStrengthS.vector()[sdofmap[cnt]] = afs[0]
		#activefiberDelayS.vector()[sdofmap[cnt]] = afd[0]


		cnt += 1


	writeXMLUGrid(ugrid, casedir+"fiber.vtu")

	nodal_points = vtk.vtkPoints()
	nodal_vertices = vtk.vtkCellArray()
	nodal_ugrid = vtk.vtkUnstructuredGrid()
	cnt = 0
	for pt in mesh.coordinates():
		nodal_points.InsertNextPoint([pt[0], pt[1], pt[2]])
		vertex = vtk.vtkVertex()
		vertex.GetPointIds().SetId(0, cnt)
		nodal_vertices.InsertNextCell(vertex)
		cnt += 1  
	nodal_ugrid.SetPoints(nodal_points)
	nodal_ugrid.SetCells(0, nodal_vertices)
	addLocalProlateSpheroidalDirections(nodal_ugrid, pdata_endo, pdata_epi, type_of_support="point")
	addLocalFiberOrientation(nodal_ugrid, endo_angle, epi_angle,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,cstrainFactor=cstrainFactor,activeFiberStrength=activeFiberStrength,activeFiberDelay=activeFiberDelay,fiberDirectionVec=fiberDirectionVec, type_of_support="point")
	writeXMLUGrid(nodal_ugrid, casedir+"fiber_nodal.vtu")
	
	return fiberV, sheetV, sheetnormV, cV, lV, rV,cstrainFactorS,activefiberStrengthS,activefiberDelayS


def addLVfiber_AVsplit(mesh, S,V, casename, endo_angle, epi_angle,  casedir, shifted_endoAtr_pdata,shifted_endoVen_pdata, shifted_epiAtr_pdata, shifted_epiVen_pdata,shifted_inlet_pdata, shifted_avj_pdata, shifted_outlet_pdata,clipheight=0.05,fiberSheetletAngle=0.,fiberSheetletWidth=0.,radialFiberAngle=0.,fiberLength=0.,cstrainFactor=1.,activeFiberStrength=1.,activeFiberDelay=0.,fiberDirectionVec=None):


	fiberV = Function(V)
	sheetV = Function(V)
	sheetnormV = Function(V)
	cV = Function(V)
	lV = Function(V)
	rV = Function(V)
	cstrainFactorS=Function(V)
	activefiberStrengthS = Function(V)
	activefiberDelayS = Function(V)

	ugrid=vtk_py.convertXMLMeshToUGrid(mesh)


	# Quad points
	gdim = mesh.geometry().dim()
	xdofmap = V.sub(0).dofmap().dofs()
	ydofmap = V.sub(1).dofmap().dofs()
	zdofmap = V.sub(2).dofmap().dofs()
	sdofmap = S.dofmap().dofs()

	xq = V.tabulate_dof_coordinates().reshape((-1, gdim))
	xq0 = xq[xdofmap]  
	#if(dolfin.dolfin_version() == '1.6.0'):
	#	xq = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, gdim))
	#	xq0 = xq[xdofmap]  

	# Create an unstructured grid of Gauss Points
	points = vtk.vtkPoints()
	vertices = vtk.vtkCellArray()
	ugrid = vtk.vtkUnstructuredGrid()
	cnt = 0;
	for pt in xq0:
		points.InsertNextPoint([pt[0], pt[1], pt[2]])
		vertex = vtk.vtkVertex()
		vertex.GetPointIds().SetId(0, cnt)
		vertices.InsertNextCell(vertex)
		cnt += 1

	ugrid.SetPoints(points)
	ugrid.SetCells(0, vertices)

	CreateVertexFromPoint(ugrid)
	addLocalProlateSpheroidalDirections_AVsplit(ugrid, shifted_endoAtr_pdata,shifted_endoVen_pdata, shifted_epiAtr_pdata, shifted_epiVen_pdata,shifted_inlet_pdata, shifted_avj_pdata, shifted_outlet_pdata, type_of_support="cell")
	addLocalFiberOrientation(ugrid, endo_angle, epi_angle,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,cstrainFactor=cstrainFactor,activeFiberStrength=activeFiberStrength,activeFiberDelay=activeFiberDelay,fiberDirectionVec=fiberDirectionVec, type_of_support="cell")

	fiber_vector =  ugrid.GetCellData().GetArray("fiber vectors")
	sheet_vector =  ugrid.GetCellData().GetArray("sheet vectors")
	sheetnorm_vector =  ugrid.GetCellData().GetArray("sheet normal vectors")
	
	eCC_vector =  ugrid.GetCellData().GetArray("eCC")
	eLL_vector =  ugrid.GetCellData().GetArray("eLL")
	eRR_vector =  ugrid.GetCellData().GetArray("eRR")
    
	csf_vector =  ugrid.GetCellData().GetArray("cstrain_factor")
	afs_vector =  ugrid.GetCellData().GetArray("fiber_active_fiber_strength")
	afd_vector =  ugrid.GetCellData().GetArray("fiber_active_fiber_delay")

 
	cnt = 0
	for pt in xq0:

		fvec = fiber_vector.GetTuple(cnt)
		svec = sheet_vector.GetTuple(cnt)
		nvec = sheetnorm_vector.GetTuple(cnt)

		cvec = eCC_vector.GetTuple(cnt)
		lvec = eLL_vector.GetTuple(cnt)
		rvec = eRR_vector.GetTuple(cnt)

		csf = csf_vector.GetTuple(cnt)
		afs = afs_vector.GetTuple(cnt)
		afd = afd_vector.GetTuple(cnt)
        
		fvecnorm = sqrt(fvec[0]**2 + fvec[1]**2 + fvec[2]**2)
		svecnorm = sqrt(svec[0]**2 + svec[1]**2 + svec[2]**2)
		nvecnorm = sqrt(nvec[0]**2 + nvec[1]**2 + nvec[2]**2)


		if(abs(fvecnorm - 1.0) > 1e-7 or  abs(svecnorm - 1.0) > 1e-6 or abs(nvecnorm - 1.0) > 1e-7):
			print (fvecnorm)
			print (svecnorm)
			print (nvecnorm)

		#print xdofmap[cnt], ydofmap[cnt], zdofmap[cnt]
		fiberV.vector()[xdofmap[cnt]] = fvec[0]; fiberV.vector()[ydofmap[cnt]] = fvec[1]; fiberV.vector()[zdofmap[cnt]] = fvec[2];
		sheetV.vector()[xdofmap[cnt]] = svec[0]; sheetV.vector()[ydofmap[cnt]] = svec[1]; sheetV.vector()[zdofmap[cnt]] = svec[2];
		sheetnormV.vector()[xdofmap[cnt]] = nvec[0]; sheetnormV.vector()[ydofmap[cnt]] = nvec[1]; sheetnormV.vector()[zdofmap[cnt]] = nvec[2];

		cV.vector()[xdofmap[cnt]] = cvec[0];  cV.vector()[ydofmap[cnt]] = cvec[1]; cV.vector()[zdofmap[cnt]] = cvec[2]; 
		lV.vector()[xdofmap[cnt]] = lvec[0];  lV.vector()[ydofmap[cnt]] = lvec[1]; lV.vector()[zdofmap[cnt]] = lvec[2]; 
		rV.vector()[xdofmap[cnt]] = rvec[0];  rV.vector()[ydofmap[cnt]] = rvec[1]; rV.vector()[zdofmap[cnt]] = rvec[2]; 
        
		cstrainFactorS.vector()[xdofmap[cnt]] = csf[0]; cstrainFactorS.vector()[ydofmap[cnt]] = csf[1]; cstrainFactorS.vector()[zdofmap[cnt]] = csf[2];
		activefiberStrengthS.vector()[xdofmap[cnt]] = afs[0]; activefiberStrengthS.vector()[ydofmap[cnt]] = afs[1]; activefiberStrengthS.vector()[zdofmap[cnt]] = afs[2];
		activefiberDelayS.vector()[xdofmap[cnt]] = afd[0]; activefiberDelayS.vector()[ydofmap[cnt]] = afd[1]; activefiberDelayS.vector()[zdofmap[cnt]] = afd[2];
		#cstrainFactorS.vector()[sdofmap[cnt]] = csf[0]
		#activefiberStrengthS.vector()[sdofmap[cnt]] = afs[0]
		#activefiberDelayS.vector()[sdofmap[cnt]] = afd[0]
		


		cnt += 1


	writeXMLUGrid(ugrid, casedir+"fiber.vtu")
    
	nodal_points = vtk.vtkPoints()
	nodal_vertices = vtk.vtkCellArray()
	nodal_ugrid = vtk.vtkUnstructuredGrid()
	cnt = 0
	for pt in mesh.coordinates():
		nodal_points.InsertNextPoint([pt[0], pt[1], pt[2]])
		vertex = vtk.vtkVertex()
		vertex.GetPointIds().SetId(0, cnt)
		nodal_vertices.InsertNextCell(vertex)
		cnt += 1  
	nodal_ugrid.SetPoints(nodal_points)
	nodal_ugrid.SetCells(0, nodal_vertices)
	addLocalProlateSpheroidalDirections_AVsplit(nodal_ugrid, shifted_endoAtr_pdata,shifted_endoVen_pdata, shifted_epiAtr_pdata, shifted_epiVen_pdata,shifted_inlet_pdata, shifted_avj_pdata, shifted_outlet_pdata, type_of_support="point")
	addLocalFiberOrientation(nodal_ugrid, endo_angle, epi_angle,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,cstrainFactor=cstrainFactor,activeFiberStrength=activeFiberStrength,activeFiberDelay=activeFiberDelay,fiberDirectionVec=fiberDirectionVec, type_of_support="point")
	writeXMLUGrid(nodal_ugrid, casedir+"fiber_nodal.vtu")
	
	return fiberV, sheetV, sheetnormV, cV, lV, rV,cstrainFactorS,activefiberStrengthS,activefiberDelayS