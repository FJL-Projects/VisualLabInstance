#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"
#include <sstream>
#include <vtkKochanekSpline.h>
#include <vtkGlyph3DMapper.h>
#include <vtkRenderWindowInteractor.h>
#include "PolygonMovementInteractor.h"

vtkRenderPipeline* pipeline;

PolygonMovementInteractorStyle plane_movement_interactor;

vtkSmartPointer<vtkPolyData> plane_pd;
vtkSmartPointer<vtkActor> plane_actor;
vtkSmartPointer<vtkActor> line_actor;
Vector_3 normal(0, 1, 0);

void compute_rotation_axis(int start[3], int end[3], double axis[3])
{
	double p1[3] = { double(start[0]), double(start[1]), 0.0 };
	double p2[3] = { double(end[0]), double(end[1]), 0.0 };

	// 这里假设Z坐标为0，因为通常在2D平面内操作
	vtkMath::Subtract(p2, p1, axis);
	vtkMath::Normalize(axis);
	axis[2] = 0.0;  // 因为是在2D平面内，所以Z分量为0
}

void RightPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Press" << endl;
}
void RightRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Released" << endl;
	std::cout << "Normal: " << normal << std::endl;
}

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	plane_movement_interactor.caller = caller;
	plane_movement_interactor.OnLeftButtonDown();
	pipeline->Renderer->RemoveActor(line_actor);
	pipeline->RenderWindow->Render();
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	plane_movement_interactor.caller = caller;
	plane_movement_interactor.OnMouseMove();
}
void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	plane_movement_interactor.caller = caller;
	plane_movement_interactor.OnLeftButtonUp();
	double3 center;
	vtkSmartPointer<vtkPolyData> pd = plane_movement_interactor.GetPolyData();
	pd->GetCenter(center.data);
	SurfaceMesh sm = VTK_PolyData2CGAL_Surface_Mesh(pd);
	face_descriptor face = *(sm.faces_begin());
	normal = PMP::compute_face_normal(face, sm);
	Point_3 source_point(center[0], center[1], center[2]);
	Point_3 target_point = source_point + normal * 20;

	line_actor = RenderLine(source_point, target_point, pipeline->Renderer, 1, 1, 0, 1);
	pipeline->Renderer->RemoveActor(line_actor);
	line_actor->GetProperty()->SetLineWidth(2);
	pipeline->Renderer->AddActor(line_actor);
	pipeline->RenderWindow->Render();
}
	

int main()
{
	pipeline = new vtkRenderPipeline();

	vtkNew<vtkPlaneSource> plane;
	plane->SetNormal(0.0, 1.0, 0.0);
	plane->Update();

	double width = 100.0; // 宽度
	double height = 100.0; // 高度

	// 设置平面的角点
	plane->SetOrigin(-width / 2, -height / 2, 0); // 原点在中心下方和左边
	plane->SetPoint1(width / 2, -height / 2, 0);  // X方向的终点
	plane->SetPoint2(-width / 2, height / 2, 0);  // Y方向的终点
	plane->Update();

	plane_pd = plane->GetOutput();
	
	double3 center;
	plane_pd->GetCenter(center.data);
	SurfaceMesh sm = VTK_PolyData2CGAL_Surface_Mesh(plane_pd);
	face_descriptor face = *(sm.faces_begin());
	Vector_3 normal = PMP::compute_face_normal(face, sm);
	Point_3 source_point(center[0], center[1], center[2]);
	Point_3 target_point = source_point + normal * 20;

	line_actor = RenderLine(source_point, target_point, pipeline->Renderer, 1, 1, 0, 1);
	pipeline->Renderer->RemoveActor(line_actor);
	line_actor->GetProperty()->SetLineWidth(2);
	pipeline->Renderer->AddActor(line_actor);

	//vtkNew<vtkPolyDataMapper> mapper;
	//mapper->SetInputData(plane_pd);

	//plane_actor = vtkSmartPointer<vtkActor>::New();
	//plane_actor->SetMapper(mapper);
	//plane_actor->GetProperty()->SetColor(0, 0.5, 0.5);
	//plane_actor->GetProperty()->SetOpacity(1);
	//pipeline->Renderer->AddActor(plane_actor);
	
	plane_movement_interactor.SetRenderer(pipeline->Renderer);
	plane_movement_interactor.SetRenderWindow(pipeline->RenderWindow);
	plane_movement_interactor.SetCorrectedOcclusalDirection(double3(0, 1, 0));
	plane_movement_interactor.SetColor(CGAL::Color(0, 255, 255));
	plane_movement_interactor.SetOpacity(1);
	plane_movement_interactor.SetPolyDataAndRender(plane_pd);


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

		bool render_flag = false;
		if (meshes_pd == *(meshes_pd_vec.begin()))
		{
			render_flag = true;
		}
		else
		{
			render_flag = false;
		}
		for (auto& mesh_pd : meshes_pd)
		{
			if (render_flag)
			{
				RenderPolydata(mesh_pd, pipeline->Renderer, r, g, b, opacity);
			}
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

		//for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
		//{
		//	double p[3];
		//	points->GetPoint(i, p);
		//	RenderSphere(p, 0.5, pipeline->Renderer, 1, 0, 0, 1);
		//}
	}

	vtkIdType spline_id = 0;
	vtkIdType arch_id = 1;
	vtkIdType teeth_id = 1;
	std::ofstream ofs("data/teeth_data.ini");
	for (auto& meshes_pd : meshes_pd_vec)
	{
		vtkSmartPointer<vtkPoints> teeth_center = teeth_center_vec[spline_id++];
		vtkSmartPointer<vtkPoints> spline_teeth_center_nearest_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkIntArray> teeth_center_id = vtkSmartPointer<vtkIntArray>::New();
		std::vector<vtkIdType> teeth_center_id_in_spline;
		teeth_center_id->SetName("teeth_center");

		std::vector<std::vector<double>> all_spline_points;
		vtkNew<vtkNamedColors> color;

		// Create splines for x, y, and z coordinates
		vtkNew<vtkKochanekSpline> x_spline;
		vtkNew<vtkKochanekSpline> y_spline;
		vtkNew<vtkKochanekSpline> z_spline;

		vtkNew<vtkParametricSpline> spline;
		spline->SetXSpline(x_spline);
		spline->SetYSpline(y_spline);
		spline->SetZSpline(z_spline);
		spline->SetPoints(teeth_center);

		vtkNew<vtkParametricFunctionSource> function_source;
		function_source->SetParametricFunction(spline);
		function_source->SetUResolution(200);
		function_source->SetVResolution(200);
		function_source->SetWResolution(200);
		function_source->Update();

		auto spline_polydata = function_source->GetOutput();
		vtkIdList* pt_id = vtkIdList::New();

		// Get the points from the spline polydata
		for (vtkIdType i = 0; i < spline_polydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			pt_id->InsertNextId(i);
		}

		vtkPoints* fp = vtkPoints::New();
		spline_polydata->GetPoints()->GetPoints(pt_id, fp);

		// Store the spline points in a vector
		for (vtkIdType i = 0; i < spline_polydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			std::vector<double> tmp;
			tmp.push_back(fp->GetPoint(i)[0]);
			tmp.push_back(fp->GetPoint(i)[1]);
			tmp.push_back(fp->GetPoint(i)[2]);
			all_spline_points.push_back(tmp);
		}

		
		for (vtkIdType teeth_center_id = 0; teeth_center_id < teeth_center->GetNumberOfPoints(); ++teeth_center_id)
		{
			double min_distance = std::numeric_limits<double>::max();
			vtkIdType min_distance_id = 0;
			for (vtkIdType i = 0; i < spline_polydata->GetPoints()->GetNumberOfPoints(); i++)
			{
				double p[3];
				fp->GetPoint(i, p);
				double q[3];
				teeth_center->GetPoint(teeth_center_id, q);
				double distance = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
				if (distance < min_distance)
				{
					min_distance = distance;
					min_distance_id = i;
				}
				if (distance > min_distance)
				{
					break;
				}
			}
			teeth_center_id_in_spline.push_back(min_distance_id);
		}
		teeth_id = 1;
		for (auto it = teeth_center_id_in_spline.begin(); it != teeth_center_id_in_spline.end(); ++it)
		{	
			auto id = *it;
			spline_teeth_center_nearest_points->InsertNextPoint(all_spline_points[id][0], all_spline_points[id][1], all_spline_points[id][2]);
			int teeth_fdi = 10 * arch_id + teeth_id;
			std::cout << "Teeth " << teeth_fdi << " ";
			ofs << teeth_fdi << "="; 
			Vector_3 direction;
			if (it == teeth_center_id_in_spline.end() - 1)  // Reached the last element
			{
				vtkIdType prev_id = id - 1;
				Point_3 prev_point(all_spline_points[prev_id][0], all_spline_points[prev_id][1], all_spline_points[prev_id][2]);
				Point_3 curr_point(all_spline_points[id][0], all_spline_points[id][1], all_spline_points[id][2]);
				direction = Vector_3(prev_point, curr_point);
				direction = direction / std::sqrt(direction.squared_length());
				std::cout << direction << std::endl;
				ofs << direction << std::endl;

				Point_3 source_point(prev_point);
				Point_3 target_point(prev_point + direction * 10);
				RenderLine(source_point, target_point, pipeline->Renderer, 0, 1, 0, 1);
			}
			else
			{
				vtkIdType next_id = id + 1;
				Point_3 curr_point(all_spline_points[id][0], all_spline_points[id][1], all_spline_points[id][2]);
				Point_3 next_point(all_spline_points[next_id][0], all_spline_points[next_id][1], all_spline_points[next_id][2]);
				direction = Vector_3(curr_point, next_point);
				direction = direction / std::sqrt(direction.squared_length());
				std::cout << direction << std::endl;
				ofs << direction << std::endl;

				Point_3 source_point(curr_point);
				Point_3 target_point(curr_point + direction * 10);
				RenderLine(source_point, target_point, pipeline->Renderer, 0, 1, 0, 1);
			}
			ofs.flush();
			teeth_id++;
		}
		std::cout << std::endl;
		arch_id++;
		// Glyph the points
		vtkNew<vtkSphereSource> points_actor;
		points_actor->SetPhiResolution(21);
		points_actor->SetThetaResolution(21);
		points_actor->SetRadius(0.3);

		// Create a polydata to store the points
		vtkNew<vtkPolyData> point_polydata;
		point_polydata->SetPoints(spline_teeth_center_nearest_points);

		vtkNew<vtkGlyph3DMapper> point_mapper;
		point_mapper->SetInputData(point_polydata);
		point_mapper->SetSourceConnection(points_actor->GetOutputPort());

		 //Create an plane_actor for the points
		vtkNew<vtkActor> point_actor;
		point_actor->SetMapper(point_mapper);
		point_actor->GetProperty()->SetColor(color->GetColor3d("Peacock").GetData());
		point_actor->GetProperty()->SetOpacity(1);
		pipeline->Renderer->AddActor(point_actor);

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

		// Create an plane_actor for the spline line
		vtkNew<vtkPolyDataMapper> line_mapper;
		line_mapper->SetInputConnection(function_source->GetOutputPort());

		vtkNew<vtkActor> line_actor;
		line_actor->SetMapper(line_mapper);
		line_actor->GetProperty()->SetColor(color->GetColor3d("Red").GetData());
		line_actor->GetProperty()->SetLineWidth(1.0);

		//// Glyph the points
		//vtkNew<vtkSphereSource> points_actor;
		//points_actor->SetPhiResolution(21);
		//points_actor->SetThetaResolution(21);
		//points_actor->SetRadius(1);

		//// Create a polydata to store the points
		//vtkNew<vtkPolyData> point_polydata;
		//point_polydata->SetPoints(teeth_center);

		//vtkNew<vtkGlyph3DMapper> point_mapper;
		//point_mapper->SetInputData(point_polydata);
		//point_mapper->SetSourceConnection(points_actor->GetOutputPort());

		// //Create an plane_actor for the points
		//vtkNew<vtkActor> point_actor;
		//point_actor->SetMapper(point_mapper);
		//point_actor->GetProperty()->SetColor(color->GetColor3d("Peacock").GetData());
		//point_actor->GetProperty()->SetOpacity(1);

		//m_spline_actor_vec.push_back(point_actor);
		//m_spline_actor_vec.push_back(line_actor);
		//pipeline->Renderer->AddActor(point_actor);
		pipeline->Renderer->AddActor(line_actor);

	}
	
	ofs.close();

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