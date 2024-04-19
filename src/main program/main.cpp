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
	using namespace Eigen;

	pipeline = new vtkRenderPipeline();

	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName("data/sample_upper_left_upsampled.vtp");
	reader->Update();
	vtkSmartPointer<vtkPolyData> arch_pd = reader->GetOutput();

	SurfaceMesh arch_sm;
	CGAL::IO::read_VTP("data/sample_upper_left_upsampled.vtp", arch_sm);
	MatrixXd HN;
	SparseMatrix<double> L, M, Minv;
	Eigen::MatrixXd	m_V; ///< Eigen matrix for vertices.
	Eigen::MatrixXi	m_F; ///< Eigen matrix for faces.

	CGALSurfaceMeshToEigen(arch_sm, m_V, m_F);
	igl::cotmatrix(m_V, m_F, L);
	igl::massmatrix(m_V, m_F, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);



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