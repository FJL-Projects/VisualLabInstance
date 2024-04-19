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
	 Converts a VTK PolyData object to a CGAL surface mesh.
	 @param P The VTK PolyData object to convert.
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

/**
 * @brief Converts a CGAL SurfaceMesh to Eigen matrix representations.
 *
 * This function takes a CGAL SurfaceMesh object and converts it into two Eigen matrices, one
 * for the vertices (V) and one for the faces (F). The vertex matrix (V) contains the x, y, and z
 * coordinates of each vertex, and the face matrix (F) contains indices into V that define each face.
 *
 * @param sm The input CGAL SurfaceMesh to convert.
 * @param[out] V An Eigen::MatrixXd that will contain the vertex coordinates after conversion.
 *               It will have as many rows as there are vertices in the mesh, and each row will
 *               have 3 columns corresponding to the x, y, and z coordinates of the vertex.
 * @param[out] F An Eigen::MatrixXi that will contain the indices of the vertices defining each face
 *               of the mesh after conversion. It will have as many rows as there are faces in the
 *               mesh, and each row will have 3 columns corresponding to the indices of the vertices
 *               that form the triangular face.
 */
void CGALSurfaceMeshToEigen(const SurfaceMesh& sm, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	// Get the number of vertices and faces
	size_t num_vertices = sm.number_of_vertices();
	size_t num_faces = sm.number_of_faces();

	// Initialize Eigen matrices
	V.resize(num_vertices, 3);
	F.resize(num_faces, 3);

	// Map for storing vertex indices (CGAL vertex descriptor -> continuous index)
	std::unordered_map<vertex_descriptor, Eigen::Index> vertex_indices;
	// Extract vertices
	Eigen::Index v_index = 0;
	for (auto vd : sm.vertices())
	{
		const auto& point = sm.point(vd);
		V.row(v_index) << point.x(), point.y(), point.z();
		vertex_indices[vd] = v_index++; // Save the index mapping
	}

	// Extract faces
	int f_index = 0;
	for (auto f : sm.faces())
	{
		int inner_index = 0;
		for (auto v : CGAL::vertices_around_face(sm.halfedge(f), sm))
		{
			F(f_index, inner_index++) = vertex_indices[v];
		}
		++f_index;
	}
}