#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"

vtkRenderPipeline* pipeline;

/**
 * @brief Generate a 4-digit number string with leading zeros.
 *
 * This function takes an integer and converts it into a string representation
 * with a fixed width of 4 characters. If the number has less than 4 digits,
 * leading zeros are added to pad the string to the desired width.
 *
 * @param number The input integer to be converted.
 * @return A string representation of the input number with leading zeros.
 *
 * @note The function uses std::ostringstream, std::setw(), and std::setfill()
 *       to format the output string.
 *
 * @example
 *   int num = 42;
 *   std::string num_str = generate_leading_zero_number_str(num);
 *   // num_str will be "0042"
 */
std::string generate_leading_zero_number_str(int number)
{
	std::ostringstream stream;
	stream << std::setw(4) << std::setfill('0') << number;
	return stream.str();
}

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
	pipeline = new vtkRenderPipeline();

	SurfaceMesh abutment;
	CGAL::IO::read_PLY("data/extracted_abutment.ply", abutment);
	SurfaceMesh inflated_abutment = abutment;

	auto& normal_map = inflated_abutment.add_property_map<vertex_descriptor, Vector_3>("v:normal").first;
	PMP::compute_vertex_normals(inflated_abutment, normal_map);

	for (auto& v : inflated_abutment.vertices())
	{
		auto& n = normal_map[v];
		inflated_abutment.point(v) = inflated_abutment.point(v) + n * 0.05;
	}



	vtkSmartPointer<vtkPolyData> abutment_pd = CGAL_Surface_Mesh2VTK_PolyData(abutment);
	RenderPolydata(abutment_pd, pipeline->Renderer, 0, 0, 1, 1);

	SurfaceMesh convex_abutment;

	CGAL::convex_hull_3(inflated_abutment, convex_abutment);

	PMP::isotropic_remeshing(convex_abutment.faces(), 0.25, convex_abutment);

	vtkSmartPointer<vtkPolyData> convex_abutment_pd = CGAL_Surface_Mesh2VTK_PolyData(convex_abutment);
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(convex_abutment_pd);
	mapper->Update();

	vtkSmartPointer<vtkActor> convex_abutment_actor = vtkSmartPointer<vtkActor>::New();
	convex_abutment_actor->SetMapper(mapper);
	convex_abutment_actor->GetProperty()->SetColor(0.51, 0.51, 0.51);
	convex_abutment_actor->GetProperty()->SetOpacity(0.7);
	convex_abutment_actor->GetProperty()->SetEdgeVisibility(true);
	pipeline->Renderer->AddActor(convex_abutment_actor);

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