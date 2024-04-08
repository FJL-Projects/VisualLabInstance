#include"stdafx.h"

vtkSmartPointer<vtkPolyData> CGAL_Surface_Mesh2VTK_PolyData(SurfaceMesh& pmesh)
{
	vtkSmartPointer<vtkPoints>  Points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> TrianglePolys = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyData> Polydata = vtkSmartPointer<vtkPolyData>::New();
	std::map<vertex_descriptor, int> VertexDescriptorIndexMap;
	int VertexIndex = 0;
	int TriIndex = 0;

	for (vertex_descriptor v : CGAL::vertices(pmesh))
	{
		const boost::property_map_value<SurfaceMesh, CGAL::vertex_point_t>::type& p = get(CGAL::get(CGAL::vertex_point, pmesh), v);
		Points->InsertNextPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
		VertexDescriptorIndexMap[v] = VertexIndex++;
	}

	for (face_descriptor f : CGAL::faces(pmesh))
	{
		TriIndex = 0;
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		for (halfedge_descriptor h : CGAL::halfedges_around_face(CGAL::halfedge(f, pmesh), pmesh))
			triangle->GetPointIds()->SetId(TriIndex++, VertexDescriptorIndexMap[CGAL::target(h, pmesh)]);
		TrianglePolys->InsertNextCell(triangle);
	}
	Polydata->SetPoints(Points);
	Polydata->SetPolys(TrianglePolys);
	return Polydata;
}

int SurfaceMeshToPolyData(SurfaceMesh& mesh, vtkSmartPointer<vtkPolyData>& polyData)
{
	//polyData->Initialize();
	polyData = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints>  VTKPoints = vtkSmartPointer<vtkPoints>::New();
	VTKPoints->SetNumberOfPoints(mesh.vertices().size());

#pragma omp parallel for
	for (int i = 0; i < mesh.vertices().size(); i++)
	{
		const Kernel::Point_3& point = mesh.point(vertex_descriptor(i));
		VTKPoints->SetPoint(i, point.x(), point.y(), point.z());
	}

	vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
	connectivity->SetNumberOfComponents(4);
	connectivity->SetNumberOfTuples(mesh.faces().size());

	//#pragma omp parallel for
	for (int j = 0; j < mesh.faces().size(); j++)
	{
		int ids[3];
		int index = 0;
		for (const halfedge_descriptor& h : mesh.halfedges_around_face(mesh.halfedge(face_descriptor(j))))
		{
			ids[index] = mesh.target(h).idx();
			index++;
		}
		connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
	}

	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	cells->SetCells(mesh.faces().size(), connectivity);

	polyData->SetPoints(VTKPoints);
	polyData->SetPolys(cells);
	//polyData->Modified();
	return 0;

}

/**
	 Converts a VTK m_polydata object to a CGAL surface mesh.
	 @param P The VTK m_polydata object to convert.
	 @return The converted CGAL surface mesh.
 */
SurfaceMesh VTK_PolyData2CGAL_Surface_Mesh(const vtkSmartPointer<vtkPolyData>& P)
{
	// Check input validity
	if (!P || P->GetNumberOfPoints() == 0 || P->GetNumberOfCells() == 0) return SurfaceMesh();

	SurfaceMesh M;
	// List of vertex indices for CGAL mesh
	std::vector< SurfaceMesh::Vertex_index> vlist;

	// Add vertices
	for (int i = 0; i < P->GetNumberOfPoints(); i++)
		vlist.push_back(M.add_vertex(Kernel::Point_3(P->GetPoints()->GetPoint(i)[0], P->GetPoints()->GetPoint(i)[1], P->GetPoints()->GetPoint(i)[2])));

	// Add faces
	for (int i = 0; i < P->GetNumberOfCells(); i++)
		M.add_face(vlist[P->GetCell(i)->GetPointId(0)], vlist[P->GetCell(i)->GetPointId(1)], vlist[P->GetCell(i)->GetPointId(2)]);

	return M;
}
