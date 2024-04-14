#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"
#include <sstream>
#include <vtkKochanekSpline.h>
#include <vtkGlyph3DMapper.h>
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
	
	std::vector<std::vector<SurfaceMesh> > meshes_vec;
	std::vector<vtkSmartPointer<vtkPoints> > teeth_center_vec;
	std::vector<std::vector<vtkSmartPointer<vtkPolyData> > > meshes_pd_vec;
	for (int i = 1; i < 5; ++i)
	{
		std::vector<SurfaceMesh> meshes;
		std::vector<vtkSmartPointer<vtkPolyData> > meshes_pd;
		for (int j = 1; j < 8; ++j)
		{
			std::stringstream ss;
			ss << "data/" << i * 10 + j;
			std::string file_name = ss.str() + ".stl";

			SurfaceMesh mesh;
			CGAL::IO::read_STL(file_name, mesh);
			meshes.push_back(mesh);
			meshes_pd.push_back(CGAL_Surface_Mesh2VTK_PolyData(mesh));
		}
		meshes_vec.push_back(meshes);
		meshes_pd_vec.push_back(meshes_pd);
	}

	double r = 0.5;
	double g = 0.25;
	double b = 0.25;
	double opacity = 0.7;
	for (auto& meshes_pd : meshes_pd_vec)
	{
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkIdType id = 0;
		for (auto& mesh_pd : meshes_pd)
		{
			RenderPolydata(mesh_pd, pipeline->Renderer, r, g, b, opacity);
			double* bound = mesh_pd->GetBounds();
			double center[3] = {
				(bound[0] + bound[1]) / 2,
				(bound[2] + bound[3]) / 2,
				(bound[4] + bound[5]) / 2 
			};
			points->InsertPoint(id++, center);
		}
		teeth_center_vec.push_back(points);
		g += 0.25;
		b += 0.25;

		for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
		{
			double p[3];
			points->GetPoint(i, p);
			RenderSphere(p, 0.5, pipeline->Renderer, 1, 0, 0, 1);
		}
	}

	vtkIdType spline_id = 0;
	for (auto& meshes_pd : meshes_pd_vec)
	{
		vtkSmartPointer<vtkPoints> teeth_center = teeth_center_vec[spline_id++];

		std::vector<std::vector<double>> all_spline_points;
		vtkNew<vtkNamedColors> color;

		// Create splines for x, y, and z coordinates
		vtkNew<vtkKochanekSpline> xSpline;
		vtkNew<vtkKochanekSpline> ySpline;
		vtkNew<vtkKochanekSpline> zSpline;

		vtkNew<vtkParametricSpline> spline;
		spline->SetXSpline(xSpline);
		spline->SetYSpline(ySpline);
		spline->SetZSpline(zSpline);
		spline->SetPoints(teeth_center);

		vtkNew<vtkParametricFunctionSource> functionSource;
		functionSource->SetParametricFunction(spline);
		functionSource->SetUResolution(200);
		functionSource->SetVResolution(200);
		functionSource->SetWResolution(200);
		functionSource->Update();

		auto splinepolydata = functionSource->GetOutput();
		vtkIdList* ptId = vtkIdList::New();

		// Get the points from the spline polydata
		for (vtkIdType i = 0; i < splinepolydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			ptId->InsertNextId(i);
		}

		vtkPoints* fp = vtkPoints::New();
		splinepolydata->GetPoints()->GetPoints(ptId, fp);

		// Store the spline points in a vector
		for (vtkIdType i = 0; i < splinepolydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			std::vector<double> tmp;
			tmp.push_back(fp->GetPoint(i)[0]);
			tmp.push_back(fp->GetPoint(i)[1]);
			tmp.push_back(fp->GetPoint(i)[2]);
			all_spline_points.push_back(tmp);
		}

		// Create spheres for each spline point and add them as actors
		for (int i = 0; i < all_spline_points.size(); i++)
		{
			vtkNew<vtkSphereSource> vtk_sphere;
			vtk_sphere->SetCenter(all_spline_points[i][0], all_spline_points[i][1], all_spline_points[i][2]);
			vtk_sphere->SetRadius(static_cast<double>(0.2));
			vtk_sphere->SetPhiResolution(21);
			vtk_sphere->SetThetaResolution(21);
			vtk_sphere->Update();
			vtkNew<vtkPolyDataMapper> vtk_dot_mapper;
			vtk_dot_mapper->SetInputConnection(vtk_sphere->GetOutputPort());
			vtkSmartPointer<vtkActor> dot_actor = vtkSmartPointer<vtkActor>::New();
			dot_actor->GetProperty()->SetColor(color->GetColor3d("Green").GetData());
			dot_actor->SetMapper(vtk_dot_mapper);
			pipeline->Renderer->AddActor(dot_actor);
			//m_spline_actor_vec.push_back(dot_actor);
		}

		// Create an actor for the spline line
		vtkNew<vtkPolyDataMapper> LineMapper;
		LineMapper->SetInputConnection(functionSource->GetOutputPort());

		vtkNew<vtkActor> LineActor;
		LineActor->SetMapper(LineMapper);
		LineActor->GetProperty()->SetColor(color->GetColor3d("Red").GetData());
		LineActor->GetProperty()->SetLineWidth(1.0);

		// Glyph the points
		vtkNew<vtkSphereSource> PointsSphere;
		PointsSphere->SetPhiResolution(21);
		PointsSphere->SetThetaResolution(21);
		PointsSphere->SetRadius(1);

		// Create a polydata to store the points
		vtkNew<vtkPolyData> PointPolyData;
		PointPolyData->SetPoints(teeth_center);

		vtkNew<vtkGlyph3DMapper> pointMapper;
		pointMapper->SetInputData(PointPolyData);
		pointMapper->SetSourceConnection(PointsSphere->GetOutputPort());

		// Create an actor for the points
		vtkNew<vtkActor> PointActor;
		PointActor->SetMapper(pointMapper);
		PointActor->GetProperty()->SetColor(color->GetColor3d("Peacock").GetData());
		PointActor->GetProperty()->SetOpacity(1);

		//m_spline_actor_vec.push_back(PointActor);
		//m_spline_actor_vec.push_back(LineActor);
		pipeline->Renderer->AddActor(PointActor);
		pipeline->Renderer->AddActor(LineActor);
	}
	

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