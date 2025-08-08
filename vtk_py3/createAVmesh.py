########################################################################

import sys
import vtk
import os
import inspect
import vtk_py3 as vtk_py

########################################################################

def createAVmesh(casename, meshsize, endofilename,split_avj_refinement,endo_atr_id=0,endo_ven_id=1,epi_atr_id=2,epi_ven_id=3,avj_id=4, verbose=True):

    if (verbose): print ('*** createAVmesh ***')

    cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
    savePath=os.path.dirname(os.path.abspath(casename))
    LVgeofile = cur_dir + "/AV_split.geo"
    LVtempgeofile = savePath+"/AV_splittemp.geo"
    mshfilename = casename + ".msh"
    vtkfilename = casename + ".vtk"
    print('cur_dir',cur_dir)
    cmd = "cp " + LVgeofile + " " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<mesh_d>>'/'" + str(meshsize) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<Endofilename>>'/'" + endofilename.replace('/','\/') + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<endo_atr_id>>'/'" + str(endo_atr_id) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<endo_ven_id>>'/'" + str(endo_ven_id) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<epi_atr_id>>'/'" + str(epi_atr_id) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<epi_ven_id>>'/'" + str(epi_ven_id) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<avj_id>>'/'" + str(avj_id) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<rsizemin>>'/'" + '{0:.1f}'.format(split_avj_refinement[0]) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<rsizemax>>'/'" + '{0:.1f}'.format(split_avj_refinement[1]) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<rdistmin>>'/'" + '{0:.1f}'.format(split_avj_refinement[2]) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<rdistmax>>'/'" + '{0:.1f}'.format(split_avj_refinement[3]) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "gmsh -3 "+LVtempgeofile+" -o " + mshfilename
    os.system(cmd)
    cmd = "gmsh -3 "+LVtempgeofile+" -o " + vtkfilename
    os.system(cmd)
    #cmd = "rm "+LVtempgeofile
    #os.system(cmd)split_avj_refinement



 
