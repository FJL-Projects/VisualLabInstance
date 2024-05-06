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
#include <MRMesh/MRRingIterator.h>

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
	using namespace MR;
	pipeline = new vtkRenderPipeline();

	constexpr double MEAN_CURVATURE_THRESHOLD = -1.0;

	Mesh mrmesh = *MR::MeshLoad::fromAnySupportedFormat("data/38461 UpperJawScan.stl");
	VertBitSet verts_under_threshold(mrmesh.topology.vertSize());
	VertBitSet verts_visited(mrmesh.topology.vertSize());
	VertBitSet verts_extracted(mrmesh.topology.vertSize());

	std::ofstream curvature_ofs("data/curvature.xyz");

	for (auto& v : mrmesh.topology.getValidVerts())
	{
		if (mrmesh.discreteMeanCurvature(v) < MEAN_CURVATURE_THRESHOLD)
		{
			curvature_ofs << mrmesh.points[v].x << " " << mrmesh.points[v].y << " " << mrmesh.points[v].z << "\n";
			verts_under_threshold.set(v);
		}
	}
	curvature_ofs.close();

	for (auto& v : mrmesh.topology.getValidVerts())
	{
		if (verts_under_threshold.test(v))
		{
			if (verts_visited.test(v))
			{
				continue;
			}
			else
			{
				std::vector<VertId> visited_verts;
				visited_verts.push_back(v);
				verts_visited.set(v);
				std::queue<VertId> vert_queue;
				vert_queue.push(v);
				while (!vert_queue.empty())
				{
					VertId current_vert = vert_queue.front();
					vert_queue.pop();
					for (auto& next_edge : orgRing(mrmesh.topology, current_vert))
					{
						VertId& next_vert = mrmesh.topology.dest(next_edge);
						if (verts_visited.test(next_vert))
						{
							continue;
						}
						else if (verts_under_threshold.test(next_vert))
						{
							verts_visited.set(next_vert);
							vert_queue.push(next_vert);
							visited_verts.push_back(next_vert);
						}
					}
				}
				if (visited_verts.size() > 20)
				{
					std::cout << "Visited Verts: " << visited_verts.size() << std::endl;
					for (auto& vert : visited_verts)
					{
						verts_extracted.set(vert);
					}
				}
			}
		}
	}

	std::ofstream filtered_verts_ofs("data/filtered_verts.xyz");
	for (auto& v : mrmesh.topology.getValidVerts())
	{
		if (verts_extracted.test(v))
		{
			filtered_verts_ofs << mrmesh.points[v].x << " " << mrmesh.points[v].y << " " << mrmesh.points[v].z << "\n";
		}
	}
	filtered_verts_ofs.close();

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