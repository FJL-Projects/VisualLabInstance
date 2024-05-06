#include "Spline/ClosedMeshSpline.h"
#include <set>
#include <algorithm>


ClosedMeshSpline::ClosedMeshSpline()
{
	SplineActor = vtkSmartPointer<vtkActor>::New();
	SplinePolydata = vtkSmartPointer<vtkPolyData>::New();
}

ClosedMeshSpline::~ClosedMeshSpline()
{

}

void ClosedMeshSpline::initial(SurfaceMesh* sm, std::map<unsigned int, face_descriptor>& fmap, std::map<unsigned int, vertex_descriptor>& vmap, std::map<unsigned int, edge_descriptor>& emap, std::map<unsigned int, halfedge_descriptor>& hemap,  SurfaceMesh::Property_map<vertex_descriptor, Point_2>& uvmap)
{
	pMesh = sm;
	FaceMap = fmap;
	VertexMap = vmap;
	EdgeMap = emap;
	HalfedgeMap= hemap;
	UVmap = uvmap;
	bDrawOnBorderFlag = 0;
	bClosed = 0;
	nIndexType = 0;
	vtSplinePoints.reserve(3000);
	vtCtrlPoints.reserve(100);
}

void ClosedMeshSpline::UpdateSpline(std::vector<MeshPoint>& SplinePoints)
{
	bDrawOnBorderFlag = 0;
	if (vtCtrlPoints.size() < 2) return;
	//点数量等于2
	if (vtCtrlPoints.size() == 2)
	{
		//裁剪直线
		vtSplinePoints.clear();
		vtSplinePoints.emplace_back(vtCtrlPoints[0]);
		FindAllCrosses(pMesh->halfedge(FaceMap[vtCtrlPoints[0].nTriId]), vtCtrlPoints[0].uv, vtCtrlPoints[1].uv, true);
		vtEquiSubscript.resize(2);
		vtEquiSubscript[0] = 0;
		vtCtrlSubscript.resize(2);
		vtCtrlSubscript[0] = 0;
		vtCtrlSubscript[1] = vtSplinePoints.size() - 1;
	}
	//点数量大于2
	if (vtCtrlPoints.size() > 2)
	{
		//计算样条线
		vtSplinePoints.clear();
		vtkSmartPointer<vtkPoints> Ctrl_Points = vtkSmartPointer<vtkPoints>::New();
		for (int i = 0; i < vtCtrlPoints.size(); i++)
		{
			double pt[3];
			pt[0] = vtCtrlPoints[i].uv[0];
			pt[1] = vtCtrlPoints[i].uv[1];
			pt[2] = static_cast<double>(0);
			Ctrl_Points->InsertNextPoint(pt);
		}
		vtkSmartPointer<vtkParametricSpline> ParametricSpline = vtkSmartPointer<vtkParametricSpline>::New();
		ParametricSpline->SetPoints(Ctrl_Points);
		if (bClosed)
		{
			ParametricSpline->ClosedOn();
		}
		vtkNew<vtkParametricFunctionSource> ParametricFunctionSource;
		ParametricFunctionSource->SetParametricFunction(ParametricSpline);
		ParametricFunctionSource->SetWResolution(100);
		ParametricFunctionSource->SetUResolution(100);
		ParametricFunctionSource->SetVResolution(100);
		//ParametricFunctionSource->SetScalarModeToDistance();
		ParametricFunctionSource->Update();
		uvSpline.clear();
		uvSpline.resize(ParametricFunctionSource->GetOutput()->GetNumberOfPoints(), Point_2(0, 0));
		for (int i = 0; i < ParametricFunctionSource->GetOutput()->GetNumberOfPoints(); i++)
		{
			double pfspt[3];
			ParametricFunctionSource->GetOutput()->GetPoint(i, pfspt);
			uvSpline[i] = Point_2(pfspt[0], pfspt[1]);
		}
		halfedge_descriptor he = pMesh->halfedge(FaceMap[vtCtrlPoints[0].nTriId]);
		vtSplinePoints.push_back(vtCtrlPoints[0]);
		//裁剪样条线
		vtCtrlSubscript.resize(vtCtrlPoints.size());
		vtEquiSubscript.resize(vtCtrlPoints.size());
		vtEquiSubscript[0] = 0;
		for (int i = 0; i < uvSpline.size() - 1; i++)
		{
			/* if (i % 20 == 0)
			{
				vtCtrlSubscript[static_cast<unsigned int>(i / 20)] = vtSplinePoints.size() - 1;
			}*/
			he = FindAllCrosses(he, Point_2(uvSpline[i][0], uvSpline[i][1]), Point_2(uvSpline[i + 1][0], uvSpline[i + 1][1]), true);
		}
		//闭合判断 塞入最后一个下标判断
		/*if (bClosed)
			vtCtrlSubscript.push_back(vtSplinePoints.size() - 1);
		else
			vtCtrlSubscript.back() = vtSplinePoints.size() - 1;*/
	}
	if (bDrawOnBorderFlag) return;
	//等距采样
	vtEquidistantSpline.clear();
	vtEquidistantSpline.emplace_back(vtSplinePoints.front());
	std::vector<MeshPoint> vtEquidistantLine;
	std::vector<MeshPoint> vtCrossPoints;
	GetEquidistantSamplingPoints(vtSplinePoints[0], vtSplinePoints.back(), vtSplinePoints, vtEquidistantSpline, 0.2);
	//for (auto& pt0 : vtCtrlPoints)
	for (int i = 0; i < vtCtrlPoints.size(); i++)
	{
		double dis = 1000.0;
		for (int j = 0; j < vtEquidistantSpline.size(); j++)
		{
			if (dis > vtkMath::Distance2BetweenPoints(vtCtrlPoints[i].xyz, vtEquidistantSpline[j].xyz))
			{
				dis = vtkMath::Distance2BetweenPoints(vtCtrlPoints[i].xyz, vtEquidistantSpline[j].xyz);
				vtEquiSubscript[i] = j;
			}
		}
	}
	if (bClosed)
		vtEquiSubscript.push_back(vtEquidistantSpline.size() - 1);
	vtEquiSubscript[0] = 0;
	//for (int i = 0; i < vtCtrlPoints.size() - 1; i++)
	//{
	//	std::vector<MeshPoint> vtCrossPoints;
	//	vtCrossPoints.assign(vtSplinePoints.begin() + vtCtrlSubscript[i], vtSplinePoints.begin() + vtCtrlSubscript[i + 1] + 1);
	//	GetEquidistantSamplingPoints(vtSplinePoints[vtCtrlSubscript[i]], vtSplinePoints[vtCtrlSubscript[i + 1]], vtCrossPoints, vtEquidistantLine, 0.2);
	//	vtEquidistantSpline.insert(vtEquidistantSpline.end(), vtEquidistantLine.begin(), vtEquidistantLine.end());
	//	vtEquiSubscript[i + 1] = vtEquidistantSpline.size() - 1;
	//}
	////闭合判断 塞入最后一个下标判断
	//if (bClosed)
	//{
	//	std::vector<MeshPoint> vtCrossPoints;
	//	vtCrossPoints.assign(vtSplinePoints.begin() + vtCtrlSubscript[vtCtrlPoints.size() - 1], vtSplinePoints.begin() + vtCtrlSubscript[vtCtrlPoints.size()] + 1);
	//	GetEquidistantSamplingPoints(vtSplinePoints[vtCtrlSubscript[vtCtrlPoints.size() - 1]], vtSplinePoints[vtCtrlSubscript[vtCtrlPoints.size()]], vtCrossPoints, vtEquidistantLine, 0.2);
	//	vtEquidistantSpline.insert(vtEquidistantSpline.end(), vtEquidistantLine.begin(), vtEquidistantLine.end());
	//	vtEquiSubscript.push_back(vtEquidistantSpline.size() - 1);
	//}


	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> PolyLine = vtkSmartPointer<vtkPolyLine>::New();
	Points->SetNumberOfPoints(vtEquidistantSpline.size());
	PolyLine->GetPointIds()->SetNumberOfIds(vtEquidistantSpline.size());
//#pragma omp parallel for
	for (int i = 0; i < vtEquidistantSpline.size(); i++)
	{
		Points->InsertPoint(i, vtEquidistantSpline[i].xyz);
		PolyLine->GetPointIds()->SetId(i, i);
	}

	vtkNew<vtkCellArray> cells;
	cells->InsertNextCell(PolyLine);
	SplinePolydata->SetPoints(Points);
	SplinePolydata->SetLines(cells);
	vtSplinePoints.pop_back();//避免首尾重复
}

/// <summary>
/// 添加样条线控制点
/// 点数量小于2:添加控制点,不做其他操作
/// 点数量等于2:添加控制点,裁剪直线,等距采样,垂直采样
/// 点数量大于2:添加控制点,计算样条线,裁剪样条线,等距采样,垂直采样
/// </summary>
/// <param name="TriangleID"> 三角面ID</param>
/// <param name="PickPosition">三角面上的点坐标</param>
/// <returns></returns>
int ClosedMeshSpline::add(unsigned int TriangleID, double PickPosition[3], double cameraDir[3])
{
	unsigned int pointSize = vtCtrlPoints.size();
	UpdateCtrlPoint(pointSize, TriangleID, PickPosition, State::ADD);
	UpdateSpline(vtEquidistantSpline);
	return 0;
}


int ClosedMeshSpline::move(unsigned int CtrlPtIndex, unsigned int TriangleID, double PickPosition[3], double cameraDir[3])
{
	UpdateCtrlPoint(CtrlPtIndex, TriangleID, PickPosition, State::MOVE);
	UpdateSpline(vtEquidistantSpline);
	return 0;
}

int ClosedMeshSpline::remove(unsigned int CtrlPtIndex)
{
	if (vtCtrlPoints.size() <= 2) return -1;
	double fTmp[3];
	UpdateCtrlPoint(CtrlPtIndex, 0, fTmp, State::REMOVE);
	UpdateSpline(vtEquidistantSpline);
	return 0;
}


/**************************************************************************************************
 *函数名：
 *时间：   2023年1月14日19:22:50
 *用户：
 *参数： unsigned int &CtrlPtIndex  鼠标按下的后一个点 序号
		 unsigned int TriangleID	鼠标按下点三角面的 ID
		 double PickPosition[3]		鼠标按下点的位置
 *返回值：
 *描述： 添加新的点
 *		注意： 这里的点位置 并不是模型上的点位置， 是管线道的位置，空间上的点。
**************************************************************************************************/
int ClosedMeshSpline::addF(unsigned int& CtrlPtIndex, unsigned int TriangleID, double PickPosition[3], double cameraDir[3])
{
	UpdateCtrlPoint(CtrlPtIndex, TriangleID, PickPosition, State::ADDF);
	UpdateSpline(vtEquidistantSpline);
	return 0;
}

/// <summary>
/// 根据裁剪出的样条线等距采样
/// </summary>
/// <param name="startCutP">输入:起始控制点用于控制宽度</param>
/// <param name="endP">输入:终止控制点用于控制宽度</param>
/// <param name="crossPt">输入:的裁剪样条线</param>
/// <param name="EquidistantSamplingPt">输出:等距采样点</param>
/// <param name="fSegmentLength">输入:等距采样点的距离</param>
void ClosedMeshSpline::GetEquidistantSamplingPoints(MeshPoint startCutP, MeshPoint endP, std::vector<MeshPoint>& crossPt, std::vector<MeshPoint>& EquidistantSamplingPt, double fSegmentLength)
{
	EquidistantSamplingPt.clear();

	//计算总长度
	double fSplineLen = 0.0;
	for (int i = 0; i < crossPt.size() - 1; i++)
		fSplineLen += (double3(crossPt[i].xyz) - double3(crossPt[i + 1].xyz)).getLength();

	//重新计算点的数量和分段的数量以及分段的长度
	int nSegmentPointNum = static_cast<int>(fSplineLen / fSegmentLength) + static_cast<int>(1);
	int nSegmentNum = nSegmentPointNum - 1;
	fSegmentLength = fSplineLen / nSegmentPointNum;
	//预留vector空间
	EquidistantSamplingPt.reserve(nSegmentNum + 10);
	//手动添加第一个点
	//EquidistantSamplingPt.emplace_back(crossPt.front());

	//fProcessLength:当前进行的长度
	double fProcessLength = 0.0;
	//nCurNum:当前的点数量
	int nCurNum = 1;
	double3 vmCur, vmNext;
	double fCurrentLength = 0.0f;
	int nNumPts = 0;
	double fLastLength = 0.0f;
	double fPreserveLength = 0.0f;
	double fMoveLen = 0.0f;
	double flag = 0.0f;
	for (int i = 0; i < crossPt.size() - 1; i++)
	{
		//vmCur为当前点 vmNext为下一个点
		//先计算两点之间距离,判断当前距离应有的点数量和已经存储的点数量是否一致
		vmCur.data[0] = crossPt[i].xyz[0];
		vmCur.data[1] = crossPt[i].xyz[1];
		vmCur.data[2] = crossPt[i].xyz[2];
		vmNext.data[0] = crossPt[i + 1].xyz[0];
		vmNext.data[1] = crossPt[i + 1].xyz[1];
		vmNext.data[2] = crossPt[i + 1].xyz[2];
		//fCurrentLength:cur和next之间的距离
		fCurrentLength = (vmNext - vmCur).getLength();
		//nNumPts:应有点的数量
		nNumPts = (fProcessLength + fCurrentLength) / fSegmentLength + 1;
		//fLastLength:上一个采样点所处的距离
		fLastLength = (nCurNum - 1) * fSegmentLength;
		//上一个采样点距离cur的距离
		fPreserveLength = fProcessLength - fLastLength;

		for (int j = 0; j < nNumPts - nCurNum; j++)
		{
			//fMoveLen:移动距离
			fMoveLen = fSegmentLength * (j + 1) - fPreserveLength;
			//vmDir:移动方向
			double3 vmDir = vmNext - vmCur;
			vmDir.normalize();
			vmDir = vmDir * fMoveLen;
			//vmEquidistantPt:采样点三维坐标
			double3 vmEquidistantPt = vmCur + vmDir;
			//计算uv和宽度并添加等距采样点
			MeshPoint SegmentPoint(pMesh->face(crossPt[i].he).idx(), GetPointUV(crossPt[i].he, Point_3(vmEquidistantPt.data[0], vmEquidistantPt.data[1], vmEquidistantPt.data[2])), vmEquidistantPt.data, crossPt[i].he);
			flag += 1.0f;
			EquidistantSamplingPt.emplace_back(SegmentPoint);
		}
		//更新点数量和当前进行的长度
		nCurNum = nNumPts;
		fProcessLength += fCurrentLength;
	}
}

void ClosedMeshSpline::BFSEdge(edge_descriptor& e, std::vector<bool>& vt_fvisit)
{
	if (pMesh->is_border(pMesh->halfedge(e)) || pMesh->is_border(pMesh->opposite(pMesh->halfedge(e))))return;
	if (vt_fvisit[pMesh->face(pMesh->halfedge(e))] == 1 || vt_fvisit[pMesh->face(pMesh->opposite(pMesh->halfedge(e)))] == 1) return;
	if (eVisit[e.idx()] == -1)
		eVisit[e.idx()] = eIndex++;
	else return;
	if (vVisit[pMesh->source(e.halfedge()).idx()] == -1)
		vVisit[pMesh->source(e.halfedge()).idx()] = vIndex++;
	if (vVisit[pMesh->target(e.halfedge()).idx()] == -1)
		vVisit[pMesh->target(e.halfedge()).idx()] = vIndex++;
	for (auto& he : pMesh->halfedges_around_target(e.halfedge()))
	{
		BFSEdge(pMesh->edge(he), vt_fvisit);
	}

	for (auto& he : pMesh->halfedges_around_target(pMesh->next(pMesh->next(e.halfedge()))))
	{
		BFSEdge(pMesh->edge(he), vt_fvisit);
	}
	return;
}

void ClosedMeshSpline::cutMesh()
{
	std::vector<bool> vt_fvisit(pMesh->number_of_faces(), 0);
	std::vector<bool> vt_vvisit(pMesh->number_of_vertices(), 0);
	for (auto v : vtSplinePoints)
	{
		vt_fvisit[v.nTriId] = 1;
	}
	vVisit.resize(pMesh->number_of_vertices(), -1);
	eVisit.resize(pMesh->number_of_edges(), -1);
	BFSEdge(pMesh->edge(HalfedgeMap[0]), vt_fvisit);

	std::vector<Point_3> ptList;
	for (int i = 0; i < vtSplinePoints.size(); i++)
	{
		ptList.push_back(Point_3(vtSplinePoints[i].xyz[0], vtSplinePoints[i].xyz[1], vtSplinePoints[i].xyz[2]));
	}
	std::vector<vertex_descriptor>insideVer;
	std::map<vertex_descriptor, int> insideVertexMap;
	for (auto v : pMesh->vertices())
	{
		if (vVisit[v.idx()] == -1 && !pMesh->is_border(v))
		{
			insideVer.push_back(v);
			insideVertexMap.insert(std::make_pair(v, insideVer.size() - 1));
			ptList.push_back(get(CGAL::get(CGAL::vertex_point, *pMesh), v));
		}
	}
	std::vector<std::vector<int>> insideSeg;
	for (auto e : pMesh->edges())
	{
		if (vVisit[pMesh->source(e.halfedge())] == -1 && vVisit[pMesh->target(e.halfedge())] == -1
			&& !pMesh->is_border(pMesh->source(e.halfedge())) && !pMesh->is_border(pMesh->target(e.halfedge())))
		{
			auto p0 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->source(e.halfedge()));
			auto p1 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->source(e.halfedge()));
			bool toonearflag = 0;
			for (int i = 0; i < vtSplinePoints.size(); i++)
			{
				if (CGAL::squared_distance(p0, Point_3(vtSplinePoints[i].xyz[0], vtSplinePoints[i].xyz[1], vtSplinePoints[i].xyz[2])) < 0.25)
					toonearflag = 1;
				if (CGAL::squared_distance(p1, Point_3(vtSplinePoints[i].xyz[0], vtSplinePoints[i].xyz[1], vtSplinePoints[i].xyz[2])) < 0.25)
					toonearflag = 1;
			}
			if (toonearflag) continue;
			insideSeg.push_back(std::vector<int>{insideVertexMap[pMesh->source(e.halfedge())], insideVertexMap[pMesh->target(e.halfedge())]});
		}
	}

	Triangulator2D tri2D;
	for (int i = 0; i < vtSplinePoints.size(); i++)
	{
		tri2D.AddPoint(vtSplinePoints[i].uv[0], vtSplinePoints[i].uv[1], i);
	}
	for (int i = 0; i < insideVer.size(); i++)
	{
		tri2D.AddPoint(UVmap[insideVer[i]][0], UVmap[insideVer[i]][1], i + vtSplinePoints.size());
	}
	for (int i = 0; i < vtSplinePoints.size() - 1; i++)
	{
		tri2D.AddSegment(i, i + 1);
	}
	//cout << insideVer.size()<<endl;
	for (int i = 0; i < insideSeg.size(); i++)
	{
		tri2D.AddSegment(insideSeg[i][0] + vtSplinePoints.size(), insideSeg[i][1] + vtSplinePoints.size());
		//cout<<insideSeg[i][0]<<" "<<insideSeg[i][1] <<endl;
	}
	tri2D.AddSegment(vtSplinePoints.size() - 1, 0);
	tri2D.SetEnclosingSegmentsProvided(01);
	tri2D.SetSubdivideAnySegments(0);
	tri2D.SetSubdivideOuterSegments(0);
	tri2D.Compute();
	std::vector<vertex_descriptor> vlist(tri2D.GetOutputVertexCount());
	double point[3];
	int nmarker;
	for (int i = 0; i < tri2D.GetOutputVertexCount(); i++)
	{
		tri2D.GetOutputVertex(i, point, nmarker);

		vlist[i] = cut.add_vertex(ptList[i]);
	}
	CGAL::IO::write_OFF("cut.off", cut);
	//cout << tri2D.GetOutputTriCount();
	int tri[3];
	for (int i = 0; i < tri2D.GetOutputTriCount(); i++)
	{
		tri2D.GetOutputTri(i, tri);
		cut.add_face(vlist[tri[0]], vlist[tri[1]], vlist[tri[2]]);
	}
	CGAL::IO::write_polygon_mesh("cut.stl", cut);
	//auto cutpolydata = CGALSurfaceMesh2VTKPolyData(cut);
	//Renderer->AddActor(MakeActor(cutpolydata, 1, 0, 0));
	//RenderWindow->Render();
}



int ClosedMeshSpline::ChangeType(int IndexType)
{
	if (IndexType == nIndexType) return 1;
	nIndexType = IndexType;
	return 0;
}

/**************************************************************************************************
 *函数名： UpdateCtrlPoint
 *时间：   2023年1月29日15:13:58
 *用户：   张旭坤
 *参数：	unsigned int &nIndex,    添加新的点时，传入的索引是样条线上的索引。其他均为控制点索引。
		int TriangleID,			 三角面 索引所在 ID
		double PickPosition[3],   点击点的位置
		State state				 状态标志
 *返回值：无
 *描述：  更新并显示控制点的位置。
**************************************************************************************************/
void ClosedMeshSpline::UpdateCtrlPoint(unsigned int& nIndex, int TriangleID, double PickPosition[3], State state)
{
	// add新增控制点
	if (state == State::ADD)
	{
		vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
		sphere->SetCenter(PickPosition[0], PickPosition[1], PickPosition[2]);
		sphere->SetRadius(0.3);
		sphere->SetPhiResolution(32);
		sphere->SetThetaResolution(32);
		sphere->Update();
		CtrlPointSphere.push_back(sphere);
		CtrlPointActor.push_back(vtkSmartPointer<vtkActor>::New());
		/*std::cout << "TriangleID: " << TriangleID << std::endl;
		std::cout << "FaceMap[TriangleID]: " << FaceMap[TriangleID] << std::endl;
		std::cout << "Halfedge: " << pMesh->halfedge(FaceMap[TriangleID]) << std::endl;
		std::cout << "GetPointUV: " << GetPointUV(pMesh->halfedge(FaceMap[TriangleID]), Point_3(PickPosition[0], PickPosition[1], PickPosition[2])) << std::endl;*/
		vtCtrlPoints.emplace_back(TriangleID, GetPointUV(pMesh->halfedge(FaceMap[TriangleID]), Point_3(PickPosition[0], PickPosition[1], PickPosition[2])), PickPosition, pMesh->halfedge(FaceMap[TriangleID]));
	}
	if (state == State::ADDF)
	{
		unsigned int nIndex1 = nIndex;
		nIndex = std::distance(vtEquiSubscript.begin(), std::lower_bound(vtEquiSubscript.begin(), vtEquiSubscript.end(), nIndex));
		vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
		sphere->SetCenter(vtEquidistantSpline[nIndex1].xyz[0], vtEquidistantSpline[nIndex1].xyz[1], vtEquidistantSpline[nIndex1].xyz[2]);
		sphere->SetRadius(0.3);
		sphere->SetPhiResolution(32);
		sphere->SetThetaResolution(32);
		sphere->Update();
		std::vector< vtkSmartPointer<vtkSphereSource>>::iterator itSphere = CtrlPointSphere.begin();
		std::vector<vtkSmartPointer<vtkActor>>::iterator itActor = CtrlPointActor.begin();
		std::vector<MeshPoint>::iterator itMeshPoint = vtCtrlPoints.begin();

		std::advance(itSphere, nIndex);  //advance  迭代器前移多少
		std::advance(itActor, nIndex);
		std::advance(itMeshPoint, nIndex);
		CtrlPointSphere.emplace(itSphere, sphere);  //emplace 在 哪一个位置上，直接生成一个元素
		CtrlPointActor.emplace(itActor, vtkSmartPointer<vtkActor>::New());
		vtCtrlPoints.emplace(itMeshPoint, vtEquidistantSpline[nIndex1]);


	}
	if (state == State::MOVE)
	{
		CtrlPointSphere[nIndex]->SetCenter(PickPosition[0], PickPosition[1], PickPosition[2]);
		vtCtrlPoints[nIndex].he = pMesh->halfedge(FaceMap[TriangleID]);
		vtCtrlPoints[nIndex].nTriId = TriangleID;
		vtCtrlPoints[nIndex].uv = GetPointUV(vtCtrlPoints[nIndex].he, Point_3(PickPosition[0], PickPosition[1], PickPosition[2]));
		vtCtrlPoints[nIndex].xyz[0] = PickPosition[0]; vtCtrlPoints[nIndex].xyz[1] = PickPosition[1]; vtCtrlPoints[nIndex].xyz[2] = PickPosition[2];

	}
	if (state == State::REMOVE)
	{
		std::vector< vtkSmartPointer<vtkSphereSource>>::iterator itSphere = CtrlPointSphere.begin();
		std::vector<vtkSmartPointer<vtkActor>>::iterator itActor = CtrlPointActor.begin();
		std::vector<MeshPoint>::iterator itMeshPoint = vtCtrlPoints.begin();
		std::advance(itSphere, nIndex);
		std::advance(itActor, nIndex);
		std::advance(itMeshPoint, nIndex);
		CtrlPointSphere.erase(itSphere);
		CtrlPointActor.erase(itActor);
		vtCtrlPoints.erase(itMeshPoint);
	}
}

//void ClosedMeshSpline::cutMesh(std::vector<rms::rms::VFTriangleMesh*>newTwoMeshes)
//{
//	//深拷贝mesh
//	SurfaceMesh* pCutMesh = new SurfaceMesh;
//	*pCutMesh = *pMesh;
//	//样条线点数量
//	volatile int numPoints = vtSplinePoints.size();
//	//新增样条线顶点的索引表
//	std::vector<vertex_descriptor> vtAddVertices;
//	vtAddVertices.reserve(numPoints);
//	//样条线经过原有mesh的半边描述符集合
//	std::vector<halfedge_descriptor> vtHalfedges;
//	vtHalfedges.reserve(numPoints);
//	// 遍历样条线上的点，添加到新mesh上，并记录经过的半边
//	for (MeshPoint& v : vtSplinePoints)
//	{
//		vtAddVertices.emplace_back(pMesh->add_vertex(Point_3(v.xyz[0], v.xyz[1], v.xyz[2])));
//		if (v.he == vtHalfedges.back())
//			continue;
//		else
//			vtHalfedges.emplace_back(v.he);
//	}
//
//	// 对于每个经过的半边，找到它所在的三角形的顶点索引
//	std::vector<std::vector<vertex_descriptor>> vt2TriVertices;
//	vt2TriVertices.resize(vtHalfedges.size(), std::vector<vertex_descriptor>(3));
//	volatile int nTri = 0;
//	for (halfedge_descriptor& he : vtHalfedges)
//	{
//		vt2TriVertices[nTri][0] = pMesh->target(he);
//		vt2TriVertices[nTri][1] = pMesh->target(pMesh->next(he));
//		vt2TriVertices[nTri][2] = pMesh->source(he);
//		nTri++;
//	}
//	// 删除原mesh中经过的半边所在的三角形
//	//for (halfedge_descriptor& he : vtHalfedges)
//	//{
//	//	pCutMesh->remove_face(pCutMesh->face(he));
//	//}
//
//	std::vector<bool> eFlag;
//	eFlag.resize(vtSplinePoints.size(), 0);
//	for (int i = 0; i < vtSplinePoints.size(); i++)
//	{
//		if (point_on_edge(vtSplinePoints[i].uv, UVmap[pMesh->source(vtSplinePoints[i].he)], UVmap[pMesh->target(vtSplinePoints[i].he)]))
//		{
//			vtSplinePoints[i].uv = Point_2(UVmap[pMesh->source(vtSplinePoints[i].he)].x()+ UVmap[pMesh->target(vtSplinePoints[i].he)].x(), UVmap[pMesh->source(vtSplinePoints[i].he)].y()*0.5 + UVmap[pMesh->target(vtSplinePoints[i].he)].y());
//			eFlag[i] = 1;
//		}
//		else
//		{
//			eFlag[i] = 0;
//		}
//		//cout << eFlag[i] << std::endl;
//		std::cout << vtSplinePoints[i].he;
//	}
//
//	nTri = 0;
//	bool first = 1;
//	//对于每一个经过的三角形，找到与样条线相交的点并使用Triangulator2D进行三角剖分
//	for (halfedge_descriptor& he : vtHalfedges)
//	{
//		// 获取三角形顶点的UV坐标
//		Point_2 uv[3];
//		uv[0] = UVmap[vt2TriVertices[nTri][0]];
//		uv[1] = UVmap[vt2TriVertices[nTri][1]];
//		uv[2] = UVmap[vt2TriVertices[nTri][2]];
//
//		// 使用rms的Triangulator2D工具类，进行三角剖分
//		rms::Triangulator2D tri2D;
//		tri2D.AddPoint(CGAL::to_double(uv[0].x()), CGAL::to_double(uv[0].y()), 0);
//		tri2D.AddPoint(CGAL::to_double(uv[1].x()), CGAL::to_double(uv[1].y()), 1);
//		tri2D.AddPoint(CGAL::to_double(uv[2].x()), CGAL::to_double(uv[2].y()), 2);
//
//		std::cout << CGAL::to_double(uv[0].x()) << " " << CGAL::to_double(uv[0].y()) << " " << 0 << std::endl;
//		std::cout << CGAL::to_double(uv[1].x()) << " " << CGAL::to_double(uv[1].y()) << " " << 0 << std::endl;
//		std::cout << CGAL::to_double(uv[2].x()) << " " << CGAL::to_double(uv[2].y()) << " " << 0 << std::endl;
//		tri2D.AddSegment(0, 1);
//		tri2D.AddSegment(1, 2);
//		tri2D.AddSegment(2, 0);
//
//		// 循环遍历样条线上的所有点，并在三角形平面上进行插值计算，得到插值点的UV坐标
//		volatile int start;
//
//		start= -1;
//
//		for (int nflag=-10; nflag < numPoints; nflag++)
//		{
//			if ((vtSplinePoints[(nflag+ numPoints) % numPoints].he == he || pMesh->next(vtSplinePoints[(nflag + numPoints) % numPoints].he) == he || pMesh->next(pMesh->next(vtSplinePoints[(nflag + numPoints) % numPoints].he)) == he) && eFlag[(nflag + numPoints) % numPoints])
//			{
//				start = (nflag + numPoints) % numPoints;
//				break;
//			}
//		}
//	
//		volatile int end = -1;
//		for (end = start + 1;; end = (end + 1) % numPoints)
//		{
//			if (eFlag[end % numPoints] == 1)
//				break;
//		}
//		//std::cout << start << " " << end << std::endl;
//		volatile int index = 3;
//			if (start>end)
//			{
//				for (int i = start; i < numPoints; i++)
//				{
//					tri2D.AddPoint(CGAL::to_double(vtSplinePoints[i].uv.x()), CGAL::to_double(vtSplinePoints[i].uv.y()), index);
//					std::cout << CGAL::to_double(vtSplinePoints[i].uv.x()) << " " << CGAL::to_double(vtSplinePoints[i].uv.y()) << " " << 0 << std::endl;
//					if (index > 3)
//					{
//						tri2D.AddSegment(index - 1, index);
//					}
//					index++;
//				}
//				for (int i = 0; i < end; i++)
//				{
//					tri2D.AddPoint(CGAL::to_double(vtSplinePoints[i].uv.x()), CGAL::to_double(vtSplinePoints[i].uv.y()), index);
//					std::cout << CGAL::to_double(vtSplinePoints[i].uv.x()) << " " << CGAL::to_double(vtSplinePoints[i].uv.y()) << " " << 0 << std::endl;
//					if (index > 3)
//					{
//						tri2D.AddSegment(index - 1, index);
//					}
//					index++;
//				}
//			}
//			else 
//			{
//				for (int i = start; i < end + 1; i = (i + 1) % numPoints)
//				{
//					tri2D.AddPoint(CGAL::to_double(vtSplinePoints[i].uv.x()), CGAL::to_double(vtSplinePoints[i].uv.y()), index);
//					std::cout << CGAL::to_double(vtSplinePoints[i].uv.x()) << " " << CGAL::to_double(vtSplinePoints[i].uv.y()) << " " << 0 << std::endl;
//					if (index > 3)
//					{
//						tri2D.AddSegment(index - 1, index);
//					}
//					index++;
//				}
//			}
//
//		tri2D.SetEnclosingSegmentsProvided(true);
//		tri2D.SetSubdivideAnySegments(false);
//		volatile bool bOK = tri2D.Compute();
//		if (!bOK)
//		{
//			continue;
//		}
//
//		for (volatile int i = 0; i < tri2D.GetOutputTriCount(); i++)
//		{
//			int p[3];
//			volatile int p1, p2, p3;
//			tri2D.GetOutputTri(i, p);
//			p1=p[0], p2=p[1], p3=p[2];
//			vertex_descriptor v[3];
//			for (int j = 0; j < 3; j++)
//			{
//				if (p[j] == 0)
//					v[j] = vt2TriVertices[nTri][0];
//				else if (p[j] == 1)
//					v[j] = vt2TriVertices[nTri][1];
//				else if (p[j] == 2)
//					v[j] = vt2TriVertices[nTri][2];
//				else
//				{
//					v[j] = vtAddVertices[(start + p[j] - 3) % numPoints];
//				}
//
//			}
//			//auto success=pCutMesh->add_face(v[0], v[1], v[2]);
//			std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
//			std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
//		}
//		//CGAL::IO::write_polygon_mesh("pCutMesh.ply", *pCutMesh);
//		nTri++;
//	}
//	
//	CGAL::IO::write_polygon_mesh("pCutMesh.ply", *pCutMesh);
//}
