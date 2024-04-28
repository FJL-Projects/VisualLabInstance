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
	

int main()
{
	using namespace MR;
	pipeline = new vtkRenderPipeline();

	Mesh mesh = *MR::MeshLoad::fromAnySupportedFormat("data/stitch_holes.stl");

	Mesh old_mesh(mesh);
	std::vector<EdgeId> edges = mesh.topology.findHoleRepresentiveEdges();

	if (edges.size() != 2)
	{
		std::cerr << "Two holes is required in the mesh." << std::endl;
		return 1;
	}

	buildCylinderBetweenTwoHoles(mesh);

	MeshSave::toAnySupportedFormat(mesh, "data/not_remeshed_stitched_holes.stl");

	FaceBitSet cylinder_faces(mesh.topology.faceSize());

	size_t old_mesh_face_size = old_mesh.topology.faceSize();
	std::cout << "Old mesh face size: " << old_mesh_face_size << std::endl;
	std::cout << "New mesh face size: " << mesh.topology.faceSize() << std::endl;

	// Outer params is old_mesh_face_size and ref cylinder_faces;
	// Each element in getValidFaces() is passed to the lambda function's FaceId f
	BitSetParallelFor(mesh.topology.getValidFaces(), [old_mesh_face_size, &cylinder_faces](FaceId f)
		{
			if (f.get() > old_mesh_face_size)
			{
				cylinder_faces.set(f);
			}
		});

	//for (auto f : mesh.topology.getValidFaces())
	//{
	//	if (f.get() > old_mesh_face_size)
	//	{
	//		cylinder_faces.set(f);
	//	}
	//}

	RemeshSettings remesh_settings;
	remesh_settings.targetEdgeLen = 0.5f;
	remesh_settings.region = &cylinder_faces;
	remesh(mesh, remesh_settings);

	MeshSave::toAnySupportedFormat(mesh, "data/stitched_holes.stl");
	size_t old_mesh_vert_size = old_mesh.topology.vertSize();
	std::cout << "Old mesh vert size: " << old_mesh_vert_size << std::endl;
    std::cout << "New mesh vert size: " << mesh.topology.vertSize() << std::endl;
	VertBitSet cylinder_points(mesh.topology.vertSize());
	BitSetParallelFor(mesh.topology.getValidVerts(), [old_mesh_vert_size, &cylinder_points](VertId v)
		{
			if (v.get() > old_mesh_vert_size)
			{
				cylinder_points.set(v);
			}
		});

	MR::positionVertsSmoothly(mesh, cylinder_points);
	MeshSave::toAnySupportedFormat(mesh, "data/smoothed_stitched_holes.stl");

	// MRMeshToPolyData may cause crashes. Still need to investigate why.
	vtkSmartPointer<vtkPolyData> polydata = CGAL_Surface_Mesh2VTK_PolyData(MRMeshToSurfaceMesh(mesh));
	RenderPolydata(polydata, pipeline->Renderer);

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