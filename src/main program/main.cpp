#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"
#include "PolygonMovementInteractor.h"

vtkRenderPipeline* pipeline;
PolygonMovementInteractorStyle generated_crown_scale_interactor;
Vector_3 projection_direction(0.194717, 0.911011, -0.363515);
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

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	generated_crown_scale_interactor.caller = caller;
	generated_crown_scale_interactor.OnMouseMove();
}

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	generated_crown_scale_interactor.caller = caller;
	generated_crown_scale_interactor.OnLeftButtonDown();
}
	
void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	generated_crown_scale_interactor.caller = caller;
	generated_crown_scale_interactor.OnLeftButtonUp();
}
int main()
{
	pipeline = new vtkRenderPipeline();

	SurfaceMesh teeth_sm;
	CGAL::IO::read_PLY("data/teeth.ply", teeth_sm);

	vtkSmartPointer<vtkPolyData> teeth_pd = CGAL_Surface_Mesh2VTK_PolyData(teeth_sm);

	double3 direction = { projection_direction.x(), projection_direction.y(), projection_direction.z() };
	
	double3 center(teeth_pd->GetCenter());
	Point_3 source(center[0], center[1], center[2]);
	Vector_3 direction_vector(direction[0], direction[1], direction[2]);
	direction_vector *= 10;
	Point_3 destination = source + direction_vector;

	// Create a line
	vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
	lineSource->SetPoint1(source[0], source[1], source[2]);
	lineSource->SetPoint2(destination[0], destination[1], destination[2]);
	lineSource->Update();

	RenderPolydata(lineSource->GetOutput(), pipeline->m_renderer);

	//// Create an arrow
	//vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();

	//// Create a transform that aligns the Z-axis with the desired direction
	//vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	//transform->Translate(teeth_pd->GetCenter());
	//transform->Scale(20, 20, 20);
	//double zAxis[3] = { 0.0, 0.0, 1.0 };
	//// Compute the angle and axis for the rotation
	//double angle = vtkMath::AngleBetweenVectors(zAxis, direction.data);
	//double axis[3];
	//vtkMath::Cross(zAxis, direction.data, axis);

	//// Apply the rotation if the angle is not 0
	//if (angle >= 1e-6)
	//{
	//	transform->RotateWXYZ(vtkMath::DegreesFromRadians(angle), axis);
	//}
	//transform->Update();
	//// Mapper
	//vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	//mapper->SetInputConnection(arrowSource->GetOutputPort());
	//mapper->Update();

	//// Actor
	//vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	//actor->SetMapper(mapper);
	//actor->SetUserTransform(transform);

	//pipeline->m_renderer->AddActor(actor);

	//RenderPolydata(teeth_pd, pipeline->m_renderer);
	generated_crown_scale_interactor.SetRenderer(pipeline->m_renderer);
	generated_crown_scale_interactor.SetRenderWindow(pipeline->m_render_window);
	generated_crown_scale_interactor.SetSurfaceMesh(teeth_sm);
	generated_crown_scale_interactor.SetBlockMove(true);
	generated_crown_scale_interactor.SetBlockRotate(true);
	generated_crown_scale_interactor.SetConstrainBorder(true);
	generated_crown_scale_interactor.SetCorrectedOcclusalDirection(direction);
	// Set up the camera and interactor.
	pipeline->m_renderer->GetActiveCamera()->SetParallelProjection(1);
	pipeline->m_renderer->ResetCamera();
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