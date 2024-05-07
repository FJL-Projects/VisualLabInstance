#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"
#include <MRMesh/MRMeshLoad.h>
#include <MRMesh/MRId.h>
#include <MRMesh/MRMesh.h>
#include <MRMesh/MRBitSetParallelFor.h>
#include <MRMesh/MRMeshTopology.h>
#include <MRMesh/MRExpected.h>

#include "CervicalMarginLineWrapper.h"
#include "ClosedSplineDesignInteractorStyle.h"
#include "TeethDataInitialization.h"

#include <CGAL/IO/Color.h>

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
	std::map<std::string, int> filename_teeth_id_map;
	//filename_teeth_id_map["5410.vtp"] = 5;
	//filename_teeth_id_map["4746.vtp"] = 6;
	//filename_teeth_id_map["4422.vtp"] = 6;
	//filename_teeth_id_map["1900.vtp"] = 2;
	//filename_teeth_id_map["1826.vtp"] = 3;
	//filename_teeth_id_map["11737.vtp"] = 2;
	//filename_teeth_id_map["2310.vtp"] = 3;
	//filename_teeth_id_map["3383.vtp"] = 7;
	//filename_teeth_id_map["1773.vtp"] = 3;
	//filename_teeth_id_map["1781.vtp"] = 6;

	//filename_teeth_id_map["13067.vtp"] = 6;
	//filename_teeth_id_map["19486.vtp"] = 4;
	//filename_teeth_id_map["34554.vtp"] = 5;
	//filename_teeth_id_map["38381.vtp"] = 1;
	filename_teeth_id_map["38278.vtp"] = 7;
	//filename_teeth_id_map["39221.vtp"] = 6;


	for (auto& pair : filename_teeth_id_map)
	{
		pipeline = new vtkRenderPipeline();

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
}