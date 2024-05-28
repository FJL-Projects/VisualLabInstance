#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"
#include <MRMesh/MRMeshLoad.h>
#include <MRMesh/MRMeshSave.h>
#include <MRMesh/MRId.h>
#include <MRMesh/MRMesh.h>
#include <MRMesh/MRBitSetParallelFor.h>
#include <MRMesh/MRMeshTopology.h>
#include <MRMesh/MRExpected.h>
#include <MRMesh/MRMeshFillHole.h>
#include <MRMesh/MRpositionVertsSmoothly.h>
#include <MRMesh/MRMeshDecimate.h>
#include <MRMesh/MRRingIterator.h>
#include <MRMesh/MRMeshFillHole.h>
#include <MRMesh/MRContoursCut.h>
#include <MRMesh/MRFillContour.h>
#include <MRMesh/MR2DContoursTriangulation.h>
#include <MRMesh/MRPointCloud.h>
#include <MRMesh/MRPointCloudTriangulation.h>
#include <MRMesh/MRMeshFwd.h>

vtkRenderPipeline* pipeline;

void RightPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Press" << endl;
}
void RightRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Released" << endl;
}

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Left Press" << endl;
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	
}

void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}

//bool FillHolesInSurfaceMesh(SurfaceMesh& mesh) {
//	// 找到一个洞的半边描述符
//	halfedge_descriptor h = *CGAL::halfedges(mesh).first;
//
//	// 查找洞
//	if (!CGAL::is_border(h, mesh)) {
//		for (; !CGAL::is_border(h, mesh); h = next(h, mesh));
//	}
//
//	// 存储用于填补洞的顶点
//	std::vector<face_descriptor> patched_f;
//	std::vector<vertex_descriptor> patched_v;
//	CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(
//		mesh,
//		h,
//		CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh))
//		.geom_traits(Kernel()));
//
//	return true;
//}

int main()
{
	using namespace MR;
	pipeline = new vtkRenderPipeline();

	Mesh mesh_mr = *MR::MeshLoad::fromAnySupportedFormat("data/16.stl");
	std::vector<Vector3f> boundary_points;
	PointCloud point_cloud;
	TriangulationParameters triangulation_parameters;
	triangulation_parameters.numNeighbours = 20;
	//triangulation_parameters.critAngle = PI;
	triangulation_parameters.critHoleLength = 20.0;
	// Outer params is old_mesh_face_size and ref cylinder_faces;
	// Each element in getValidFaces() is passed to the lambda function's FaceId f
	for (auto v : mesh_mr.topology.getValidVerts())
	{
		if (mesh_mr.topology.isBdVertex(v))
		{
			boundary_points.push_back(mesh_mr.points[v]);
		}
		point_cloud.addPoint(mesh_mr.points[v]);
	}

	std::cout << "Boundary points size: " << boundary_points.size() << std::endl;
	Vector3f points_sum;
	for (auto p : boundary_points)
	{
		points_sum += p;
	}
	Vector3f center = Vector3f(points_sum.x / boundary_points.size(), points_sum.y / boundary_points.size(), points_sum.z / boundary_points.size());
	point_cloud.addPoint(center);

	for (auto p : boundary_points)
	{
		auto new_point0 = Vector3f(p.x * 0.9 + center.x * 0.1, p.y * 0.9 + center.y * 0.1, p.z * 0.9 + center.z * 0.1);
		auto new_point1 = Vector3f(p.x * 0.8 + center.x * 0.2, p.y * 0.8 + center.y * 0.2, p.z * 0.8 + center.z * 0.2);
		auto new_point2 = Vector3f(p.x * 0.6 + center.x * 0.4, p.y * 0.6 + center.y * 0.4, p.z * 0.6 + center.z * 0.4);
		point_cloud.addPoint(new_point0);
		point_cloud.addPoint(new_point1);
		point_cloud.addPoint(new_point2);
	}

	std:ofstream file("data/point_cloud.xyz");
	for (auto p : point_cloud.points)
	{
		file << p.x << " " << p.y << " " << p.z << std::endl;
	}
	file.close();

	auto triangulated_mesh_res = triangulatePointCloud(point_cloud, triangulation_parameters);
	auto triangulated_mesh = triangulated_mesh_res.value();

	MeshSave::toAnySupportedFormat(triangulated_mesh, "data/triangulated_mesh.stl");

	vtkSmartPointer<vtkPolyData> triangulated_pd;
	SurfaceMeshToPolyData(MRMeshToSurfaceMesh(triangulated_mesh), triangulated_pd);
	RenderPolydata(triangulated_pd, pipeline->Renderer);

	//for (auto f : mesh.topology.getValidFaces())
	//{
	//	if (f.get() > old_mesh_face_size)
	//	{
	//		cylinder_faces.set(f);
	//	}
	//}

	//RemeshSettings remesh_settings;
	//remesh_settings.targetEdgeLen = 0.5f;
	//remesh_settings.region = &cylinder_faces;
	//remesh(mesh, remesh_settings);

	//MeshSave::toAnySupportedFormat(mesh, "data/stitched_holes.stl");
	//size_t old_mesh_vert_size = old_mesh.topology.vertSize();
	//std::cout << "Old mesh vert size: " << old_mesh_vert_size << std::endl;
 //   std::cout << "New mesh vert size: " << mesh.topology.vertSize() << std::endl;
	//VertBitSet cylinder_points(mesh.topology.vertSize());
	//BitSetParallelFor(mesh.topology.getValidVerts(), [old_mesh_vert_size, &cylinder_points](VertId v)
	//	{
	//		if (v.get() > old_mesh_vert_size)
	//		{
	//			cylinder_points.set(v);
	//		}
	//	});

	/*MR::positionVertsSmoothly(mesh, cylinder_points);
	MeshSave::toAnySupportedFormat(mesh, "data/smoothed_stitched_holes.stl");*/

	// MRMeshToPolyData may cause crashes. Still need to investigate why.
	//vtkSmartPointer<vtkPolyData> polydata = CGAL_Surface_Mesh2VTK_PolyData(MRMeshToSurfaceMesh(mesh));
	//RenderPolydata(polydata, pipeline->Renderer);

	// Set up the camera and interactor.
	pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
	pipeline->Renderer->ResetCamera();
	// Set up the callback functions for mouse events.
	pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
	pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
	pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);
	pipeline->addObserver(vtkCommand::RightButtonPressEvent, RightPress);
	pipeline->addObserver(vtkCommand::RightButtonReleaseEvent, RightRelease);

	pipeline->RenderWindowInteractor->Start();

	// Clean up the pipeline after each run.
	delete pipeline;
}