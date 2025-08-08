import vtk
import vtk_py3 as vtk_py
import dolfin as dolfin
import numpy as np

def extractFeNiCsBiVFacet(ugrid,savePath='', geometry="BiV", tol=1e-2):

	#tol = 1e-2
	
	#ugrid = vtk_py.readUGrid(meshfilename)
	
	# Extract surface
    geom = vtk.vtkGeometryFilter()
    if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        geom.SetInput(ugrid)
    else:
        geom.SetInputData(ugrid)
    geom.Update()
    surf = geom.GetOutput()
    
    bc_pts_locator = []
    bc_pts = []
    bc_pts_range = []
    bc_pts_min = []
    bc_pts_max = []
    bc_pts_map = []
    	
    	# Extract Surface Normal
    normal = vtk.vtkPolyDataNormals()
    if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        normal.SetInput(surf)
    else:
        normal.SetInputData(surf)
    normal.ComputeCellNormalsOn()
    normal.Update()
    surf_w_norm = normal.GetOutput()
    
    	#vtk_py.writePData(normal.GetOutput(), "normal.vtk")
    	
    zmax = surf_w_norm.GetBounds()[5]
    	
    surf_w_norm.BuildLinks()
    idlist = vtk.vtkIdList()
    basecellidlist = vtk.vtkIdTypeArray()
    basesurf = vtk.vtkPolyData()
    for p in range(0, surf_w_norm.GetNumberOfCells()):
        zvec = surf_w_norm.GetCellData().GetNormals().GetTuple3(p)[2]
        
        surf_w_norm.GetCellPoints(p, idlist)
        zpos = surf_w_norm.GetPoints().GetPoint(idlist.GetId(0))[2]
        
        if((abs(zvec - 1.0) < tol or abs(zvec + 1.0) < tol) and (abs(zmax - zpos) < tol)):
            surf_w_norm.DeleteCell(p)
            basecellidlist.InsertNextValue(p)
        
    basesurf = vtk_py.extractCellFromPData(basecellidlist, surf)
    baseptlocator = vtk.vtkPointLocator()
    baseptlocator.SetDataSet(basesurf)
    baseptlocator.BuildLocator()
    
    	#######################################################################
    surf_w_norm.RemoveDeletedCells()
    
    cleanpdata = vtk.vtkCleanPolyData()
    if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        cleanpdata.SetInput(surf_w_norm)
    else:
        cleanpdata.SetInputData(surf_w_norm)
    cleanpdata.Update()
    
    connfilter = vtk.vtkPolyDataConnectivityFilter()
    if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        connfilter.SetInput(cleanpdata.GetOutput())
    else:
        connfilter.SetInputData(cleanpdata.GetOutput())
    connfilter.Update()
    
    print ("Total_num_points = ",  cleanpdata.GetOutput().GetNumberOfPoints())
    tpt = 0
    
    if(geometry=="BiV"):
        nsurf = 3
    elif (geometry=="AVsplit"):
        nsurf = 5
    else:
        nsurf = 2

	
    for p in range(0,nsurf):
        
        pts = vtk.vtkPolyData()
    	
        connfilter.SetExtractionModeToSpecifiedRegions()
        [connfilter.DeleteSpecifiedRegion(k) for k in range(0,nsurf)]
        connfilter.AddSpecifiedRegion(p)
        connfilter.ScalarConnectivityOff()
        connfilter.FullScalarConnectivityOff()
        connfilter.Update()
    	
        cleanpdata2 = vtk.vtkCleanPolyData()
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
            cleanpdata2.SetInput(connfilter.GetOutput())
        else:
            cleanpdata2.SetInputData(connfilter.GetOutput())
        cleanpdata2.Update()
    	
        pts.DeepCopy(cleanpdata2.GetOutput())
    	
        tpt = tpt + cleanpdata2.GetOutput().GetNumberOfPoints()
    	
        ptlocator = vtk.vtkPointLocator()
        ptlocator.SetDataSet(pts)
        ptlocator.BuildLocator()
    	
        bc_pts_locator.append(ptlocator)
        bc_pts.append(pts)
        bc_pts_range.append([abs(pts.GetBounds()[k+1]-pts.GetBounds()[k]) for k in range(0, 6, 2)])
        bc_pts_min.append([pts.GetBounds()[k] for k in range(0, 6, 2)])
        bc_pts_max.append([pts.GetBounds()[k+1] for k in range(0, 6, 2)])


	#vtk_py.writePData(connfilter.GetOutput(), "/home/likchuan/Research/fenicsheartmesh/ellipsoidal/Geometry/test.vtk")
	
    print ("Total_num_points = ",  tpt)

    Epiid = np.argmax(np.array([pts[2] for pts in bc_pts_range]))
    maxzrank =  np.array([pts[0] for pts in bc_pts_max]).argsort()
    if(geometry=="BiV"):
        LVid = maxzrank[0] 
        RVid = 3 - (LVid + Epiid)
        bc_pts_map = [4, 4, 4, 4]
        bc_pts_map[Epiid] = 1; bc_pts_map[LVid] = 2; bc_pts_map[RVid] = 3
        baseid  = 3;
    elif(geometry=="AVsplit"):
        
        zrange=np.array([pts[2] for pts in bc_pts_range]).argsort()
        AVJid=zrange[0]
        Inletid=zrange[1]
        Outletid=zrange[2]
        if bc_pts_max[Inletid][2] < bc_pts_max[Outletid][2]:
            Inletid=zrange[2]
            Outletid=zrange[1]
        zmax=np.array([pts[2] for pts in bc_pts_max]).argsort()
        zmax.remove(AVJid)
        zmax.remove(Inletid)
        zmax.remove(Outletid)
        EpiAtrid=zmax[-1]
        EndoAtrid=zmax[-2]
        zmin=np.array([pts[2] for pts in bc_pts_min]).argsort()
        zmin.remove(AVJid)
        zmin.remove(Inletid)
        zmin.remove(Outletid)
        zmin.remove(EpiAtrid)
        zmin.remove(EndoAtrid)
        EndoVenid=zmin[1]
        EpiVenid=zmin[0]

        bc_pts_map = [7,7,7,7,7,7,7]
        bc_pts_map[EpiAtrid] = 1; bc_pts_map[EpiVenid] = 2; bc_pts_map[EndoAtrid] = 3; bc_pts_map[EndoVenid] = 4; bc_pts_map[Inletid] = 5; bc_pts_map[Outletid] = 6; bc_pts_map[AVJid] = 7
        baseid=6
    else:
        LVid = maxzrank[0]
        bc_pts_map = [4, 4, 4]
        bc_pts_map[Epiid] = 1; bc_pts_map[LVid] = 2
        baseid  = 2;
    
    bc_pts_locator.append(baseptlocator)
    bc_pts.append(basesurf)

	
    dolfin_mesh = vtk_py.convertUGridToXMLMesh(ugrid)
	#dolfin_facets = dolfin.FacetFunction('size_t', dolfin_mesh)
    dolfin_facets = dolfin.MeshFunction('size_t', dolfin_mesh,dolfin_mesh.topology().dim()-1, dolfin_mesh.domains())
    dolfin_facets.set_all(0)

    for facet in dolfin.SubsetIterator(dolfin_facets, 0):
        for locator in range(0,nsurf+1):
            cnt = 0
            for p in range(0,3):
                v0 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(0)
                v1 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(1)
                v2 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(2)
                ptid = bc_pts_locator[locator].FindClosestPoint(v0, v1, v2)
                x0 =  bc_pts[locator].GetPoints().GetPoint(ptid)
                dist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], x0)
                if(dist < 1e-5*tol):
                    cnt = cnt + 1
            if(cnt == 3):
                dolfin_facets[facet] = bc_pts_map[locator]
					

	#dolfin_edges = dolfin.EdgeFunction('size_t', dolfin_mesh)
    dolfin_edges = dolfin.MeshFunction('size_t', dolfin_mesh,1, dolfin_mesh.domains())
    dolfin_edges.set_all(0)

    #For AVsplit, avj:1,2 inlet:3,4 ; outlet:5,6 , lower number is epi, 
    if geometry=="AVsplit":
        epilocator = EpiVenid
        lvendolocator = EndoVenid
    else:
        epilocator = Epiid
        lvendolocator = LVid
    for edge in dolfin.SubsetIterator(dolfin_edges, 0):
        cnt_epi = 0; cnt_lvendo = 0;
        for p in range(0,2):
            v0 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(0)
            v1 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(1)
            v2 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(2)

            epiptid = bc_pts_locator[epilocator].FindClosestPoint(v0, v1, v2)
            epix0 =  bc_pts[epilocator].GetPoints().GetPoint(epiptid)
            epidist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], epix0)

            topptid = bc_pts_locator[baseid].FindClosestPoint(v0, v1, v2)
            topx0 =  bc_pts[baseid].GetPoints().GetPoint(topptid)
            topdist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], topx0)

            lvendoptid = bc_pts_locator[lvendolocator].FindClosestPoint(v0, v1, v2)
            lvendox0 =  bc_pts[lvendolocator].GetPoints().GetPoint(lvendoptid)
            lvendodist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], lvendox0)

            if(topdist < 1e-5*tol and epidist < 1e-5*tol):
                cnt_epi = cnt_epi + 1

            if(topdist < 1e-5*tol and lvendodist < 1e-5*tol):
                cnt_lvendo = cnt_lvendo + 1

            if(cnt_epi == 2):
                dolfin_edges[edge] = 1

            if(cnt_lvendo == 2):
                dolfin_edges[edge] = 2

    dolfin.File(savePath+"temp.pvd") << dolfin_facets

    return dolfin_mesh, dolfin_facets, dolfin_edges	

def extractFeNiCsBiVFacet_AVsplit(ugrid,savePath='', tol=1e-2,split_avj_extend=None):
	
    bc_pts_locator = []
    bc_pts = []
    bc_pts_range = []
    bc_pts_min = []
    bc_pts_max = []
    bc_pts_map = []
    thresd=[]
    nsurf=7
    #ext_cc=[]
    #ext_centroid=[]
    #cc_locator=[]
    geom = vtk.vtkGeometryFilter()
    if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        geom.SetInput(ugrid)
    else:
        geom.SetInputData(ugrid)
    geom.Update()
    surf = geom.GetOutput()
    ext_cc=vtk.vtkCellCenters()
    if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        ext_cc.SetInput(surf)
    else:
        ext_cc.SetInputData(surf)
    ext_cc.Update()
    ext_centroid=ext_cc.GetOutput()
    cc_locator=vtk.vtkPointLocator()
    cc_locator.SetDataSet(ext_centroid)
    cc_locator.BuildLocator()
    for p in range(0,nsurf+2):
        #thresd = vtk.vtkPolyData()
        pts_vtkThreshold=vtk.vtkThreshold()
        pts_vtkThreshold.SetInputData(ugrid)
        pts_vtkThreshold.ThresholdBetween(p+0.5,p+1.5)
        pts_vtkThreshold.SetInputArrayToProcess(0, 0, 0, 1, 'CellEntityIds')
        pts_vtkThreshold.Update()
        #cleanpdata2 = vtk.vtkCleanUnstructuredGrid()
        #if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        #    cleanpdata2.SetInput(pts.GetOutput())
        #else:
        #    cleanpdata2.SetInputData(pts.GetOutput())
        #cleanpdata2.Update()
        #thresd.DeepCopy(pts.GetOutput())
        thresd.append(pts_vtkThreshold.GetOutput())
        #vtk_py.writeUGrid(thresd,savePath+"/splitsurf"+str(p)+".vtk")
        bc_pts_locator.append(vtk.vtkPointLocator())
        bc_pts_locator[-1].SetDataSet(thresd[-1])
        bc_pts_locator[-1].BuildLocator()##C+++ cant get it to overwite !!!!!
    	
        bc_pts.append(thresd[-1])
        bc_pts_range.append([abs(thresd[-1].GetBounds()[k+1]-thresd[-1].GetBounds()[k]) for k in range(0, 6, 2)])
        bc_pts_min.append([thresd[-1].GetBounds()[k] for k in range(0, 6, 2)])
        bc_pts_max.append([thresd[-1].GetBounds()[k+1] for k in range(0, 6, 2)])
        '''
        ext_cc.append(vtk.vtkCellCenters())
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
            ext_cc[-1].SetInput(thresd[-1])
        else:
            ext_cc[-1].SetInputData(thresd[-1])
        ext_cc[-1].Update()
        ext_centroid.append(ext_cc[-1].GetOutput())
        cc_locator.append(vtk.vtkPointLocator())
        cc_locator[-1].SetDataSet(ext_centroid[-1])
        cc_locator[-1].BuildLocator()
        '''
    zrange=list(np.array([pts[2] for pts in bc_pts_range[:-2]]).argsort())
    AVJid=zrange[0]
    zrange.remove(AVJid)
    zmax=list(np.array([pts[2] for pts in bc_pts_max[:-2]]).argsort())
    zmax.remove(AVJid)
    EpiAtrId=zmax[-1]
    Inletid=zmax[-2]
    EndoAtrId=zmax[-3]
    if bc_pts_range[EpiAtrId][2]<bc_pts_range[Inletid][2]:
        temp=EpiAtrId
        EpiAtrId=Inletid
        Inletid=temp
    if (bc_pts_range[EndoAtrId][0]**2.+bc_pts_range[EndoAtrId][1]**2.)<(bc_pts_range[Inletid][0]**2.+bc_pts_range[Inletid][1]**2.):
        temp=EndoAtrId
        EndoAtrId=Inletid
        Inletid=temp
    zmin=list(np.array([pts[2] for pts in bc_pts_min[:-2]]).argsort())
    zmin.remove(AVJid)
    EpiVenId=zmin[0]
    Outletid=zmin[1]
    EndoVenId=zmin[2]
    if bc_pts_range[EpiVenId][2]<bc_pts_range[Outletid][2]:
        temp=EpiVenId
        EpiVenId=Outletid
        Outletid=temp
    if (bc_pts_range[EndoVenId][0]**2.+bc_pts_range[EndoVenId][1]**2.)<(bc_pts_range[Outletid][0]**2.+bc_pts_range[Outletid][1]**2.):
        temp=EndoVenId
        EndoVenId=Outletid
        Outletid=temp
    zmax.remove(EpiAtrId)
    zmax.remove(EndoAtrId)
    zmax.remove(Inletid)
    zmax.remove(EpiVenId)
    zmax.remove(EndoVenId)
    zmax.remove(Outletid)
    zmin.remove(EpiAtrId)
    zmin.remove(EndoAtrId)
    zmin.remove(Inletid)
    zmin.remove(EpiVenId)
    zmin.remove(EndoVenId)
    zmin.remove(Outletid)
    #closeInletid=zmax[-1]
    #closeOutletid=zmin[0]

    bc_pts_map = [7]*7
    bc_pts_map[EpiAtrId] = 5
    bc_pts_map[EpiVenId] = 1
    bc_pts_map[EndoAtrId] = 6
    bc_pts_map[EndoVenId] = 2
    bc_pts_map[Inletid] = 7
    bc_pts_map[Outletid] = 8
    bc_pts_map[AVJid] = 4
    #bc_pts_map[closeInletid] = 8
    #bc_pts_map[closeOutletid] = 9
    baseid=4
    print('bc_pts_map',bc_pts_map)
    
    #pts = vtk.vtkThreshold()
    #pts.SetInputData(ugrid)
    #pts.ThresholdBetween(5.5,7.5)
    #pts.SetInputArrayToProcess(0, 0, 0, 1, 'CellEntityIds')
    #pts.Update()
    #thresd=pts.GetOutput()
    dolfin_mesh = vtk_py.convertUGridToXMLMesh(ugrid)
	#dolfin_facets = dolfin.FacetFunction('size_t', dolfin_mesh)
    dolfin_facets = dolfin.MeshFunction('size_t', dolfin_mesh,dolfin_mesh.topology().dim()-1, dolfin_mesh.domains())
    dolfin_facets.set_all(0)
    rangeorder=list(range(0,nsurf))
    rangeorder.pop(AVJid)
    rangeorder.insert(0,AVJid)
    store_avj_points=np.zeros((0,3))
    for facet in dolfin.SubsetIterator(dolfin_facets, 0):
        avj_cnt=0
        for locator in rangeorder:
            cnt = 0
            c0=0.
            c1=0.
            c2=0.
            for p in range(0,3):
                v0 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(0)
                v1 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(1)
                v2 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(2)
                c0+=v0/3.
                c1+=v1/3.
                c2+=v2/3.
                ptid = bc_pts_locator[locator].FindClosestPoint(v0, v1, v2)
                x0 =  bc_pts[locator].GetPoints().GetPoint(ptid)
                dist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], x0)
                if(dist < 1e-2*tol):
                    cnt = cnt + 1
                    if locator==AVJid:
                        avj_cnt = avj_cnt + 1
            if(cnt >= 3):
                if locator==EndoAtrId and avj_cnt >=1:
                    dolfin_facets[facet] = 9
                    for p in range(0,3):
                        v0 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(0)
                        v1 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(1)
                        v2 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(2)
                        store_avj_points=np.concatenate((store_avj_points,np.array([v0,v1,v2]).reshape((1,-1))),axis=0)
                elif locator==EndoVenId and avj_cnt >=1:
                    dolfin_facets[facet] = 10
                    for p in range(0,3):
                        v0 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(0)
                        v1 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(1)
                        v2 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(2)
                        store_avj_points=np.concatenate((store_avj_points,np.array([[v0,v1,v2]])),axis=0)
                else:
                    ptid_cc = cc_locator.FindClosestPoint(c0, c1, c2)
                    x0_cc =  ext_centroid.GetPoints().GetPoint(ptid_cc)
                    dist = vtk.vtkMath.Distance2BetweenPoints([c0,c1,c2], x0_cc)
                    if(dist < 1e-2*tol):
                        dolfin_facets[facet] = bc_pts_map[locator]
                break
    if split_avj_extend is not None:
        store_avj_points=np.unique(store_avj_points,axis=0)
        init_store_avj_points_num=store_avj_points.shape[0]
        store_avj_points_prev_cnt=0
        store_avj_points_cnt=store_avj_points.shape[0]
        print('extending from',store_avj_points_cnt)
        print('split_avj_extend',split_avj_extend)
        for n in range(50):
            store_avj_points_prev_cnt=store_avj_points_cnt
            for facet in dolfin.SubsetIterator(dolfin_facets, bc_pts_map[EndoAtrId]):
                cnt=0
                cnt_exceed=0
                new_tri=np.zeros((3,3))
                new_tri_bool=[True,True,True]
                for p in range(0,3):
                    new_tri[p,0] =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(0)
                    new_tri[p,1] =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(1)
                    new_tri[p,2] =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(2)
                    dist = np.sum((store_avj_points-new_tri[p:p+1,:])**2.,axis=-1)
                    dist_0 = dist[:init_store_avj_points_num].min()
                    dist = dist.min()
                    if (dist < 1e-2*tol):
                        cnt = cnt + 1
                        new_tri_bool[p]=False
                    if (dist_0>split_avj_extend):# or (np.abs(new_tri[p,2])>split_avj_extend):
                        cnt_exceed=cnt_exceed+1
                if (cnt > 1) and cnt_exceed==0:
                    dolfin_facets[facet] = 9
                    if np.any(new_tri_bool):
                        store_avj_points=np.concatenate((store_avj_points,new_tri[new_tri_bool]),axis=0)
            for facet in dolfin.SubsetIterator(dolfin_facets, bc_pts_map[EndoVenId]):
                cnt=0
                cnt_exceed=0
                new_tri=np.zeros((3,3))
                new_tri_bool=[True,True,True]
                for p in range(0,3):
                    new_tri[p,0] =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(0)
                    new_tri[p,1] =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(1)
                    new_tri[p,2] =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(2)
                    dist = np.sum((store_avj_points-new_tri[p:p+1,:])**2.,axis=-1)
                    dist_0 = dist[:init_store_avj_points_num].min()
                    dist = dist.min()
                    if (dist < 1e-2*tol):
                        cnt = cnt + 1
                        new_tri_bool[p]=False
                    if (dist_0>split_avj_extend):
                        cnt_exceed=cnt_exceed+1
                if (cnt > 1) and cnt_exceed==0:
                    dolfin_facets[facet] = 10
                    if np.any(new_tri_bool):
                        store_avj_points=np.concatenate((store_avj_points,new_tri[new_tri_bool]),axis=0)
            store_avj_points_cnt=store_avj_points.shape[0]
            if store_avj_points_cnt==store_avj_points_prev_cnt:
                break
        print('extending to',store_avj_points_cnt,'with',n,'loop')
	#dolfin_edges = dolfin.EdgeFunction('size_t', dolfin_mesh)
    dolfin_edges = dolfin.MeshFunction('size_t', dolfin_mesh,1, dolfin_mesh.domains())
    dolfin_edges.set_all(0)

    #For AVsplit, avj:1,2 inlet:3,4 ; outlet:5,6 , lower number is epi, 
    epilocator = EpiAtrId
    endolocator = EndoAtrId
    epiVenlocator = EpiVenId
    endoVenlocator = EndoVenId
    avjlocator = AVJid
    inlocator = Inletid
    outlocator = Outletid
    for edge in dolfin.SubsetIterator(dolfin_edges, 0):
        cnt_in_epi = 0; cnt_in_endo = 0; cnt_avj_epi = 0; cnt_avj_endo = 0; cnt_out_epi = 0; cnt_out_endo = 0;
        for p in range(0,2):
            v0 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(0)
            v1 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(1)
            v2 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(2)

            epiptid = bc_pts_locator[epilocator].FindClosestPoint(v0, v1, v2)
            epix0 =  bc_pts[epilocator].GetPoints().GetPoint(epiptid)
            epidist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], epix0)
            
            laendoptid = bc_pts_locator[endolocator].FindClosestPoint(v0, v1, v2)
            laendox0 =  bc_pts[endolocator].GetPoints().GetPoint(laendoptid)
            laendodist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], laendox0)
            
            epiVenptid = bc_pts_locator[epiVenlocator].FindClosestPoint(v0, v1, v2)
            epiVenx0 =  bc_pts[epiVenlocator].GetPoints().GetPoint(epiVenptid)
            epiVendist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], epiVenx0)
            
            lvendoptid = bc_pts_locator[endoVenlocator].FindClosestPoint(v0, v1, v2)
            lvendox0 =  bc_pts[endoVenlocator].GetPoints().GetPoint(lvendoptid)
            lvendodist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], lvendox0)

            avjptid = bc_pts_locator[avjlocator].FindClosestPoint(v0, v1, v2)
            avjx0 =  bc_pts[avjlocator].GetPoints().GetPoint(avjptid)
            avjdist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], avjx0)
            
            inptid = bc_pts_locator[inlocator].FindClosestPoint(v0, v1, v2)
            inx0 =  bc_pts[inlocator].GetPoints().GetPoint(inptid)
            indist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], inx0)
            
            outptid = bc_pts_locator[outlocator].FindClosestPoint(v0, v1, v2)
            outx0 =  bc_pts[outlocator].GetPoints().GetPoint(outptid)
            outdist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], outx0)
            if(indist < 1e-5*tol and epidist < 1e-5*tol):
                cnt_in_epi = cnt_in_epi + 1

            if(indist < 1e-5*tol and laendodist < 1e-5*tol):
                cnt_in_endo = cnt_in_endo + 1
            
            if(avjdist < 1e-5*tol and epidist < 1e-5*tol):
                cnt_avj_epi = cnt_avj_epi + 1

            if(avjdist < 1e-5*tol and laendodist < 1e-5*tol):
                cnt_avj_endo = cnt_avj_endo + 1
            
            if(outdist < 1e-5*tol and epiVendist < 1e-5*tol):
                cnt_out_epi = cnt_out_epi + 1

            if(outdist < 1e-5*tol and lvendodist < 1e-5*tol):
                cnt_out_endo = cnt_out_endo + 1

        if(cnt_in_epi == 2):
            dolfin_edges[edge] = 1
        if(cnt_in_endo == 2):
            dolfin_edges[edge] = 2
        if(cnt_avj_epi == 2):
            dolfin_edges[edge] = 3
        if(cnt_avj_endo == 2):
            dolfin_edges[edge] = 4
        if(cnt_out_epi == 2):
            dolfin_edges[edge] = 5
        if(cnt_out_endo == 2):
            dolfin_edges[edge] = 6
    #For AVsplit, avj:1,2 inlet:3,4 ; outlet:5,6 , lower number is epi, 
    
    dolfin.File(savePath+"temp.pvd") << dolfin_facets

    return dolfin_mesh, dolfin_facets, dolfin_edges	, thresd[EndoAtrId], thresd[EndoVenId], thresd[EpiAtrId], thresd[EpiVenId], thresd[Inletid], thresd[AVJid], thresd[Outletid]