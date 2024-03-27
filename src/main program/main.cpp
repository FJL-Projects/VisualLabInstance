#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"

#include "CervicalMarginLineWrapper.h"
#include "ClosedSplineDesignInteractorStyle.h"


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

	vtkSmartPointer<vtkPolyData> arch_pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName("data/1.vtp");
	reader->Update();
	arch_pd = reader->GetOutput();

	SurfaceMesh arch_sm = VTK_PolyData2CGAL_Surface_Mesh(arch_pd);
	
	vtkSmartPointer<vtkActor> colored_actor = vtkSmartPointer<vtkActor>::New();


	ClosedSplineDesignInteractorStyle closed_spline_design_interactor_style;
	closed_spline_design_interactor_style.setSurfaceMesh(&arch_sm);
	closed_spline_design_interactor_style.setPolyData(arch_pd);
	closed_spline_design_interactor_style.setRenderer(pipeline->Renderer);
	closed_spline_design_interactor_style.setRenderWindow(pipeline->RenderWindow);
	closed_spline_design_interactor_style.setPolydataActor(colored_actor);
	closed_spline_design_interactor_style.Init();

	CervicalMarginLineWrapper cervical_margin_line_wrapper;
	cervical_margin_line_wrapper.SetRenderer(pipeline->Renderer);
	cervical_margin_line_wrapper.SetRenderWindow(pipeline->RenderWindow);
	cervical_margin_line_wrapper.SetArchSurfaceMesh(&arch_sm);
	cervical_margin_line_wrapper.SetArchPolyData(arch_pd);
	cervical_margin_line_wrapper.SetMinCurvatureThreshold(0.2);
	cervical_margin_line_wrapper.SetMeanCurvatureThreshold(0.2);
	cervical_margin_line_wrapper.SetCervicalMarginLineInteractorStyle(&closed_spline_design_interactor_style);

	cervical_margin_line_wrapper.SetAbutmentPolyData(teeth_polydata);
	cervical_margin_line_wrapper.SetSelectedId(selected_id);
	cervical_margin_line_wrapper.SetProjectionDirection(projection_direction_v3);

	cervical_margin_line_wrapper.Init();

	cervical_margin_line_wrapper.ExtractAbutmentSurfaceMesh();
	cervical_margin_line_wrapper.GenerateAbutmentEdgeSpline();
	cervical_margin_line_wrapper.GenerateCervicalMarginLine();

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