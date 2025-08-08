import dolfin as fenics
#import vtk as vtk
#from heartFEM.lcleeHeart.vtk_py import *
import os
try:
    import ldrb
except:
    pass
import vtk_py3 as vtk_py
#from .extractUGridBasedOnThreshold import extractUGridBasedOnThreshold as extractUGridBasedOnThreshold
#from .convertUGridtoPdata import convertUGridtoPdata as convertUGridtoPdata
#from .readXMLPUGrid import readXMLPUGrid as readXMLPUGrid
#from .readXMLUGrid import readXMLUGrid as readXMLUGrid 
#from .computeRegionsForBiV import addRegionsToBiV as addRegionsToBiV
from vtk.util import numpy_support

def GetSurfaces(directory, filebasename, fieldvariable, isparallel):
    import glob
    filenames = glob.glob(directory + "/" + filebasename + "*")
    filenames.sort()
    filenames = [filename for filename in filenames if filename[-3:] != "pvd"]

    if(filenames[0][-4:] == "pvtu" and isparallel):
        ugrid = vtk_py.readXMLPUGrid(filenames[0])
	
    elif(filenames[0][-4:] == ".vtu" and (not isparallel)):
        ugrid = vtk_py.readXMLUGrid(filenames[0])

    Epi = vtk_py.extractUGridBasedOnThreshold(ugrid, fieldvariable, 1)
    Epi  = vtk_py.convertUGridtoPdata(Epi)
    LVendo = vtk_py.extractUGridBasedOnThreshold(ugrid, fieldvariable, 2)
    LVendo  = vtk_py.convertUGridtoPdata(LVendo)
    RVendo = vtk_py.extractUGridBasedOnThreshold(ugrid, fieldvariable, 3)
    RVendo  = vtk_py.convertUGridtoPdata(RVendo)

    return LVendo, RVendo, Epi

def addBiVfiber(mesh, facetboundaries,endo_angle, epi_angle, fiberSheetletAngle=0.,saveDir=None):#,fiberSheetletWidth=0.,radialFiberAngle=0.,fiberLength=0.):
    
    angles={}
    if isinstance(endo_angle,(int,float)):
        angles['alpha_endo_lv']=endo_angle
        angles['alpha_endo_sept']=endo_angle
        angles['alpha_endo_rv']=endo_angle
    elif isinstance(endo_angle,dict):
        angles['alpha_endo_lv']=endo_angle['LV']
        angles['alpha_endo_sept']=endo_angle['septum']
        angles['alpha_endo_rv']=endo_angle['RV']
    elif len(endo_angle)==1:
        angles['alpha_endo_lv']=endo_angle[0]
        angles['alpha_endo_sept']=endo_angle[0]
        angles['alpha_endo_rv']=endo_angle[0]
    elif len(endo_angle)==2:
        angles['alpha_endo_lv']=endo_angle[0]
        angles['alpha_endo_sept']=0.5*endo_angle[0]+0.5*endo_angle[1]
        angles['alpha_endo_rv']=endo_angle[1]
    else:
        angles['alpha_endo_lv']=endo_angle[0]
        angles['alpha_endo_sept']=endo_angle[2]
        angles['alpha_endo_rv']=endo_angle[1]
        
    if isinstance(epi_angle,(int,float)):
        angles['alpha_epi_lv']=epi_angle
        angles['alpha_epi_sept']=epi_angle
        angles['alpha_epi_rv']=epi_angle
    elif isinstance(epi_angle,dict):
        angles['alpha_epi_lv']=epi_angle['LV']
        angles['alpha_epi_sept']=epi_angle['septum']
        angles['alpha_epi_rv']=epi_angle['RV']
    elif len(epi_angle)==1:
        angles['alpha_epi_lv']=epi_angle[0]
        angles['alpha_epi_sept']=epi_angle[0]
        angles['alpha_epi_rv']=epi_angle[0]
    elif len(epi_angle)==2:
        angles['alpha_epi_lv']=epi_angle[0]
        angles['alpha_epi_sept']=0.5*epi_angle[0]+0.5*epi_angle[1]
        angles['alpha_epi_rv']=epi_angle[1]
    else:
        angles['alpha_epi_lv']=epi_angle[0]
        angles['alpha_epi_sept']=epi_angle[2]
        angles['alpha_epi_rv']=epi_angle[1]
        
    if isinstance(fiberSheetletAngle,(int,float)):
        angles['beta_endo_lv']=fiberSheetletAngle
        angles['beta_endo_sept']=fiberSheetletAngle
        angles['beta_endo_rv']=fiberSheetletAngle
        angles['beta_epi_lv']=fiberSheetletAngle
        angles['beta_epi_sept']=fiberSheetletAngle
        angles['beta_epi_rv']=fiberSheetletAngle
    elif isinstance(fiberSheetletAngle,dict):
        if 'LV' in fiberSheetletAngle:
            angles['beta_endo_lv']=fiberSheetletAngle['LV']
            angles['beta_epi_lv']=fiberSheetletAngle['LV']
        else:
            angles['beta_endo_lv']=fiberSheetletAngle['LV_endo']
            angles['beta_epi_lv']=fiberSheetletAngle['LV_epi']
        if 'septum' in fiberSheetletAngle:
            angles['beta_endo_sept']=fiberSheetletAngle['septum']
            angles['beta_epi_sept']=fiberSheetletAngle['septum']
        else:
            angles['beta_endo_sept']=fiberSheetletAngle['Septum_endo']
            angles['beta_epi_sept']=fiberSheetletAngle['Septum_epi']
        if 'RV' in fiberSheetletAngle:
            angles['beta_endo_rv']=fiberSheetletAngle['RV']
            angles['beta_epi_rv']=fiberSheetletAngle['RV']
        else:
            angles['beta_endo_rv']=fiberSheetletAngle['RV_endo']
            angles['beta_epi_rv']=fiberSheetletAngle['RV_epi']
    elif len(fiberSheetletAngle)==1:
        if isinstance(fiberSheetletAngle[0],(int,float)):
            angles['beta_endo_lv']=fiberSheetletAngle[0]
            angles['beta_endo_sept']=fiberSheetletAngle[0]
            angles['beta_endo_rv']=fiberSheetletAngle[0]
            angles['beta_epi_lv']=fiberSheetletAngle[0]
            angles['beta_epi_sept']=fiberSheetletAngle[0]
            angles['beta_epi_rv']=fiberSheetletAngle[0]
        else:
            angles['beta_endo_lv']=fiberSheetletAngle[0][0]
            angles['beta_endo_sept']=fiberSheetletAngle[0][0]
            angles['beta_endo_rv']=fiberSheetletAngle[0][0]
            angles['beta_epi_lv']=fiberSheetletAngle[0][1]
            angles['beta_epi_sept']=fiberSheetletAngle[0][1]
            angles['beta_epi_rv']=fiberSheetletAngle[0][1]
    elif len(fiberSheetletAngle)==2:
        if isinstance(fiberSheetletAngle[0],(int,float)):
            angles['beta_endo_lv']=fiberSheetletAngle[0]
            angles['beta_endo_sept']=0.5*fiberSheetletAngle[0]+0.5*fiberSheetletAngle[1]
            angles['beta_endo_rv']=fiberSheetletAngle[1]
            angles['beta_epi_lv']=fiberSheetletAngle[0]
            angles['beta_epi_sept']=0.5*fiberSheetletAngle[0]+0.5*fiberSheetletAngle[1]
            angles['beta_epi_rv']=fiberSheetletAngle[1]
        else:
            angles['beta_endo_lv']=fiberSheetletAngle[0][0]
            angles['beta_endo_sept']=0.5*fiberSheetletAngle[0][0]+0.5*fiberSheetletAngle[1][0]
            angles['beta_endo_rv']=fiberSheetletAngle[1][0]
            angles['beta_epi_lv']=fiberSheetletAngle[0][1]
            angles['beta_epi_sept']=0.5*fiberSheetletAngle[0][1]+0.5*fiberSheetletAngle[1][1]
            angles['beta_epi_rv']=fiberSheetletAngle[1][1]
    else:
        if isinstance(fiberSheetletAngle[0],(int,float)):
            angles['beta_endo_lv']=fiberSheetletAngle[0]
            angles['beta_endo_sept']=fiberSheetletAngle[2]
            angles['beta_endo_rv']=fiberSheetletAngle[1]
            angles['beta_epi_lv']=fiberSheetletAngle[0]
            angles['beta_epi_sept']=fiberSheetletAngle[2]
            angles['beta_epi_rv']=fiberSheetletAngle[1]
        else:
            angles['beta_endo_lv']=fiberSheetletAngle[0][0]
            angles['beta_endo_sept']=fiberSheetletAngle[2][0]
            angles['beta_endo_rv']=fiberSheetletAngle[1][0]
            angles['beta_epi_lv']=fiberSheetletAngle[0][1]
            angles['beta_epi_sept']=fiberSheetletAngle[2][1]
            angles['beta_epi_rv']=fiberSheetletAngle[1][1]
            
    try:
        # Choose space for the fiber fields
        # This is a string on the form {family}_{degree}
        quad_deg = 4
        fiber_space = "Quadrature_"+str(quad_deg)
        ffun = fenics.MeshFunction("size_t", mesh, 2)
        ffun.set_all(0)
        markers={'base': 4, 'rv': 3, 'lv': 2, 'epi': 1}
        # Compute the microstructure
        fiberV, sheetV, sheetnormV = ldrb.dolfin_ldrb(mesh=mesh, fiber_space=fiber_space, ffun=ffun, markers=markers, **angles)
        if saveDir is not None:
            ldrb.fiber_to_xdmf(fiberV, saveDir+"fiber")
    except:
        from . import SetBiVFiber_Quad_PyQ as SetBiVFiber_Quad_PyQ
        import dolfin
        directory,meshname = os.path.split(saveDir)
        directory+='/'
        matid = dolfin.MeshFunction('size_t', mesh, 3, mesh.domains())
        LVendo, RVendo, Epi = GetSurfaces(directory, "facetboundaries000000.vtu", "f", False)
        ugrid = vtk_py.readXMLUGrid(directory + "mesh000000.vtu")
        vtk_py.addRegionsToBiV(ugrid, LVendo, RVendo, Epi)
        matid_vtk = numpy_support.vtk_to_numpy(ugrid.GetCellData().GetArray("region_id"))
        matid.array()[:] = matid_vtk
        dolfin.File(directory + "matid.pvd") << matid;
        fiber_angle_param = {"mesh": mesh,\
    	 "facetboundaries": facetboundaries,\
    	 "LV_fiber_angle": [angles['alpha_endo_lv'],angles['alpha_epi_lv']], \
    	 "LV_sheet_angle": [angles['beta_endo_lv'],angles['beta_epi_lv']], \
    	 "Septum_fiber_angle": [angles['alpha_endo_sept'],angles['alpha_epi_sept']],\
    	 "Septum_sheet_angle": [angles['beta_endo_sept'],angles['beta_epi_sept']],\
    	 "RV_fiber_angle": [angles['alpha_endo_rv'],angles['alpha_epi_rv']],\
    	 "RV_sheet_angle": [angles['beta_endo_rv'],angles['beta_epi_rv']],\
    	 "LV_matid": 0,\
    	 "Septum_matid": 1,\
    	 "RV_matid": 2,\
    	 "matid": matid,\
    	 "isrotatept": False,\
    	 "isreturn": True,\
    	 "outfilename":  meshname,\
    	 "outdirectory": directory,\
    	 "epiid": 1,\
    	 "rvid": 3,\
    	 "lvid": 2,\
    	 "degree": 4}
        fiberV, sheetV, sheetnormV = SetBiVFiber_Quad_PyQ(fiber_angle_param)
    return fiberV, sheetV, sheetnormV


