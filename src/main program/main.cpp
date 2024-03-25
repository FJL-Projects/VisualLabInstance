#include"stdafx.h"
#include"vtkRenderPipeline.h"
#include"meshTransform.h"
#include"simpleRender.h"

vtkRenderPipeline* pipeline;
SurfaceMesh toothmesh0;
SurfaceMesh toothmesh1;
SurfaceMesh crownmesh;
SurfaceMesh bitemesh;

int x_dim;
int z_dim;
double3 centroid;

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//cout<<"Left Press" << endl;
	vtkSmartPointer<vtkCamera> camera = pipeline->Renderer->GetActiveCamera();
	double3 camera_orientation = camera->GetOrientation();
	std::cout << "Camera orientation angles: " << camera_orientation[0] << ", " << camera_orientation[1] << ", " << camera_orientation[2] << std::endl;
	double angle_x = -camera_orientation[0] * CGAL_PI / 180;
	double angle_y = -camera_orientation[1] * CGAL_PI / 180;
	double angle_z = -camera_orientation[2] * CGAL_PI / 180;

	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->RotateX(-camera_orientation[0]);
	transform->RotateY(-camera_orientation[1]);
	transform->RotateZ(-camera_orientation[2]);

	vtkSmartPointer<vtkTransformPolyDataFilter> transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	vtkSmartPointer<vtkPolyData> tooth_polydata = CGAL_Surface_Mesh2VTK_PolyData(toothmesh0);
	transform_filter->SetInputData(tooth_polydata);
	transform_filter->SetTransform(transform);
	transform_filter->Update();
	vtkSmartPointer<vtkPolyData> transformed_tooth_polydata = transform_filter->GetOutput();
	vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
	writer->SetFileName("rotated_toothmesh.ply");
	writer->SetInputData(transformed_tooth_polydata);
	writer->Write();
	//Eigen::Matrix3d rotation_matrix;
	//rotation_matrix = Eigen::AngleAxisd(-angle_z, Eigen::Vector3d::UnitX())
	//	* Eigen::AngleAxisd(-angle_y, Eigen::Vector3d::UnitY())
	//	* Eigen::AngleAxisd(-angle_x, Eigen::Vector3d::UnitZ());

	//auto rotate_mesh = [&rotation_matrix](SurfaceMesh sm) -> SurfaceMesh
	//	{
	//		for (auto& v : sm.vertices())
	//		{
	//			Point_3& p = sm.point(v);
	//			Eigen::Vector3d vec(p.x(), p.y(), p.z());
	//			vec = rotation_matrix * vec;
	//			p = Point_3(vec.x(), vec.y(), vec.z());
	//		}

	//		return sm;
	//	};
	//
	//CGAL::IO::write_PLY("rotated_toothmesh.ply", rotate_mesh(toothmesh0));
	//CGAL::IO::write_PLY("rotated_toothmesh1.ply", rotate_mesh(toothmesh1));
	//CGAL::IO::write_PLY("rotated_crownmesh.ply", rotate_mesh(crownmesh));
	//CGAL::IO::write_PLY("rotated_bitemesh.ply", rotate_mesh(bitemesh));
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	
}

std::string select_folder()
{
	BROWSEINFO  bi;
	bi.hwndOwner = NULL;
	bi.pidlRoot = CSIDL_DESKTOP; 
	bi.pszDisplayName = NULL;
	bi.lpszTitle = NULL; 
	bi.ulFlags = BIF_DONTGOBELOWDOMAIN | BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE; 
	bi.lpfn = NULL;
	bi.iImage = 0;
	LPITEMIDLIST pidl = SHBrowseForFolder(&bi); 
	if (pidl == NULL)
	{
		std::cout << "" << std::endl;
		return std::string();
	}
	TCHAR folder_tchar[MAX_PATH];
	SHGetPathFromIDList(pidl, folder_tchar);
	std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
	std::string selected_folder = converter.to_bytes(folder_tchar);

	return selected_folder;
}

void writePNG(SurfaceMesh sm,double3 dir,std::string path)
{
	dir = -dir;
	dir.normalize();
	double3 y(0, 1, 0);
	double3 axis = double3::crossProduct(dir,y);
	axis.normalize();
	axis = axis * acos(double3::dotProduct(dir, y));
	double3 center(0, 0, 0);
	for (auto v : sm.vertices())
		center = center + double3(sm.point(v).x(), sm.point(v).y(), sm.point(v).z());
	center = center / sm.number_of_vertices();
	for (auto v : sm.vertices())
	{
		double3 pt(sm.point(v).x(), sm.point(v).y(), sm.point(v).z());
		pt = pt - center;
		AngleAxisRotatePoint(axis.data, pt.data, pt.data);
		pt = pt + center;
		sm.point(v) = Point_3(pt[0] , pt[1] , pt[2] );
	}
	
	Tree tree(faces(sm).first, faces(sm).second, sm);
	double x_min = std::numeric_limits<double>::max();
	double x_max = std::numeric_limits<double>::min();
	double z_min = std::numeric_limits<double>::max();
	double z_max = std::numeric_limits<double>::min();
	double y_max = std::numeric_limits<double>::min();
	for (auto v : sm.vertices())
	{
		Point_3 p = sm.point(v);
		x_min = std::min(x_min, p.x());
		x_max = std::max(x_max, p.x());
		z_min = std::min(z_min, p.z());
		z_max = std::max(z_max, p.z());
		y_max = std::max(y_max, p.y());
	}
	double max;
	if ((x_max - x_min) > (z_max - z_min))
		max = x_max - x_min;
	else
		max = z_max - z_min;

	vtkSmartPointer< vtkImageData> image = vtkSmartPointer< vtkImageData>::New();
	image->SetDimensions(1000, 1000, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	int dim[3];
	image->GetDimensions(dim);
	cout << dim[0] << " " << dim[1] << " " << dim[2] << endl;
	std::vector<std::vector<double>> depth(dim[0], std::vector<double>(dim[1], 0));
	for (int x = 0; x < dim[0]; x++)
	{
		for (int z = 0; z < dim[1]; z++)
		{
			double x_ = x_min + (x_max - x_min) * x / dim[0];
			double z_ = z_min + (z_max - z_min) * z / dim[1];
			Ray_3 ray_query(Point_3(x_, y_max, z_), Vector_3(0, -1, 0));
			auto intersection = tree.first_intersection(ray_query);
			const Point_3* p;
			if (intersection && boost::get<Point_3>(&(intersection->first)))
			{
				p = boost::get<Point_3>(&(intersection->first));
				depth[x][z] = y_max - p->y();
			}
			else
			{
				depth[x][z] = 0;
			}
		}
	}
	double depth_max = 0;
	for (int x = 0; x < 1000; x++)
	{
		for (int z = 0; z < 1000; z++)
		{
			depth_max = std::max(depth_max, depth[x][z]);
		}
	}
	for (int x = 0; x < 1000; x++)
	{
		for (int z = 0; z < 1000; z++)
		{
			if (depth[x][z] == 0)
			{
				unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, z, 0));
				pixel[0] = static_cast<int>(depth[x][z] / depth_max * 255); 
				pixel[1] = static_cast<int>(depth[x][z] / depth_max * 255);  
				pixel[2] = static_cast<int>(depth[x][z] / depth_max * 255); 
			}
			else
			{

					unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, z, 0));
					pixel[0] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);
					pixel[1] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);
					pixel[2] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);

			}
		}
	}
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(path.c_str());
	writer->SetInputData(image);
	writer->Write();
}

SurfaceMesh writePNG(SurfaceMesh sm, double3 dir, std::string path, int x_dim, int z_dim, double3 centroid)
{
	dir = -dir;
	dir.normalize();
	std::cout << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;

	double3 x(1, 0, 0);
	double3 y(0, 1, 0);
	double3 z(0, 0, 1);
	double3 axis_y = double3::crossProduct(dir, y);
	axis_y.normalize();
	axis_y = axis_y * acos(double3::dotProduct(dir, y));

	//double3 axis_z = double3::crossProduct(dir, z);
	//axis_z.normalize();
	//axis_z = axis_z * acos(double3::dotProduct(dir, z));

	//double3 axis_x = double3::crossProduct(dir, x);
	//axis_x.normalize();
	//axis_x = axis_x * acos(double3::dotProduct(dir, x));

	double3 center(0, 0, 0);
	//for (auto v : sm.vertices())
	//	center = center + double3(sm.point(v).x(), sm.point(v).y(), sm.point(v).z());
	//center = center / sm.number_of_vertices();
	for (auto v : sm.vertices())
	{
		double3 pt(sm.point(v).x(), sm.point(v).y(), sm.point(v).z());
		pt -= centroid;
		AngleAxisRotatePoint(axis_y.data, pt.data, pt.data);
		//AngleAxisRotatePoint(axis_x.data, pt.data, pt.data);
		//AngleAxisRotatePoint(axis_z.data, pt.data, pt.data);

		pt += centroid;
		sm.point(v) = Point_3(pt[0], pt[1], pt[2]);
	}

	Tree tree(faces(sm).first, faces(sm).second, sm);
	double x_min = std::numeric_limits<double>::max();
	double x_max = std::numeric_limits<double>::min();
	double z_min = std::numeric_limits<double>::max();
	double z_max = std::numeric_limits<double>::min();
	double y_max = std::numeric_limits<double>::min();
	for (auto v : sm.vertices())
	{
		Point_3 p = sm.point(v);
		x_min = std::min(x_min, p.x());
		x_max = std::max(x_max, p.x());
		z_min = std::min(z_min, p.z());
		z_max = std::max(z_max, p.z());
		y_max = std::max(y_max, p.y());
	}
	double max;
	if ((x_max - x_min) > (z_max - z_min))
		max = x_max - x_min;
	else
		max = z_max - z_min;

	vtkSmartPointer< vtkImageData> image = vtkSmartPointer< vtkImageData>::New();
	image->SetDimensions(x_dim, z_dim, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	int dim[3];
	image->GetDimensions(dim);
	cout << dim[0] << " " << dim[1] << " " << dim[2] << endl;
	std::vector<std::vector<double>> depth(dim[0], std::vector<double>(dim[1], 0));
	for (int x = 0; x < dim[0]; x++)
	{
		for (int z = 0; z < dim[1]; z++)
		{
			double x_ = x_min + (x_max - x_min) * x / dim[0];
			double z_ = z_min + (z_max - z_min) * z / dim[1];
			Ray_3 ray_query(Point_3(x_, y_max, z_), Vector_3(0, -1, 0));
			//Ray_3 ray_query(Point_3(x, y_max, z), Vector_3(0, -1, 0));
			auto intersection = tree.first_intersection(ray_query);
			const Point_3* p;
			if (intersection && boost::get<Point_3>(&(intersection->first)))
			{
				p = boost::get<Point_3>(&(intersection->first));
				depth[x][z] = y_max - p->y();
			}
			else
			{
				depth[x][z] = 0;
			}
		}
	}
	double depth_max = 0;
	for (int x = 0; x < x_dim; x++)
	{
		for (int z = 0; z < z_dim; z++)
		{
			depth_max = std::max(depth_max, depth[x][z]);
		}
	}
	for (int x = 0; x < x_dim; x++)
	{
		for (int z = 0; z < z_dim; z++)
		{
			if (depth[x][z] == 0)
			{
				unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, z, 0));
				pixel[0] = static_cast<int>(depth[x][z] / depth_max * 255);
				pixel[1] = static_cast<int>(depth[x][z] / depth_max * 255);
				pixel[2] = static_cast<int>(depth[x][z] / depth_max * 255);
			}
			else
			{

				unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, z, 0));
				pixel[0] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);
				pixel[1] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);
				pixel[2] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);

			}
		}
	}
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(path.c_str());
	writer->SetInputData(image);
	writer->Write();

	return sm;
}

SurfaceMesh rotate_sm(SurfaceMesh sm, double3 dir, std::string path, int x_dim, int z_dim, double3 centroid)
{
	dir = -dir;
	dir.normalize();
	double3 x(1, 0, 0);
	double3 y(0, 1, 0);
	double3 z(0, 0, 1);
	double3 axis_y = double3::crossProduct(dir, y);
	axis_y.normalize();
	axis_y = axis_y * acos(double3::dotProduct(dir, y));

	double3 axis_z = double3::crossProduct(dir, z);
	axis_z.normalize();
	axis_z = axis_z * acos(double3::dotProduct(dir, z));

	double3 axis_x = double3::crossProduct(dir, x);
	axis_x.normalize();
	axis_x = axis_x * acos(double3::dotProduct(dir, x));

	double3 center(0, 0, 0);
	//for (auto v : sm.vertices())
	//	center = center + double3(sm.point(v).x(), sm.point(v).y(), sm.point(v).z());
	//center = center / sm.number_of_vertices();
	for (auto v : sm.vertices())
	{
		double3 pt(sm.point(v).x(), sm.point(v).y(), sm.point(v).z());
		pt -= centroid;
		AngleAxisRotatePoint(axis_y.data, pt.data, pt.data);
		AngleAxisRotatePoint(axis_x.data, pt.data, pt.data);
		AngleAxisRotatePoint(axis_z.data, pt.data, pt.data);

		pt += centroid;
		sm.point(v) = Point_3(pt[0], pt[1], pt[2]);
	}

	Tree tree(faces(sm).first, faces(sm).second, sm);
	double x_min = std::numeric_limits<double>::max();
	double x_max = std::numeric_limits<double>::min();
	double z_min = std::numeric_limits<double>::max();
	double z_max = std::numeric_limits<double>::min();
	double y_max = std::numeric_limits<double>::min();
	for (auto v : sm.vertices())
	{
		Point_3 p = sm.point(v);
		x_min = std::min(x_min, p.x());
		x_max = std::max(x_max, p.x());
		z_min = std::min(z_min, p.z());
		z_max = std::max(z_max, p.z());
		y_max = std::max(y_max, p.y());
	}
	double max;
	if ((x_max - x_min) > (z_max - z_min))
		max = x_max - x_min;
	else
		max = z_max - z_min;

	vtkSmartPointer< vtkImageData> image = vtkSmartPointer< vtkImageData>::New();
	image->SetDimensions(x_dim, z_dim, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	int dim[3];
	image->GetDimensions(dim);
	cout << dim[0] << " " << dim[1] << " " << dim[2] << endl;
	std::vector<std::vector<double>> depth(dim[0], std::vector<double>(dim[1], 0));
	for (int x = 0; x < dim[0]; x++)
	{
		for (int z = 0; z < dim[1]; z++)
		{
			double x_ = x_min + (x_max - x_min) * x / dim[0];
			double z_ = z_min + (z_max - z_min) * z / dim[1];
			Ray_3 ray_query(Point_3(x_, y_max, z_), Vector_3(0, -1, 0));
			//Ray_3 ray_query(Point_3(x, y_max, z), Vector_3(0, -1, 0));
			auto intersection = tree.first_intersection(ray_query);
			const Point_3* p;
			if (intersection && boost::get<Point_3>(&(intersection->first)))
			{
				p = boost::get<Point_3>(&(intersection->first));
				depth[x][z] = y_max - p->y();
			}
			else
			{
				depth[x][z] = 0;
			}
		}
	}
	double depth_max = 0;
	for (int x = 0; x < x_dim; x++)
	{
		for (int z = 0; z < z_dim; z++)
		{
			depth_max = std::max(depth_max, depth[x][z]);
		}
	}
	for (int x = 0; x < x_dim; x++)
	{
		for (int z = 0; z < z_dim; z++)
		{
			if (depth[x][z] == 0)
			{
				unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, z, 0));
				pixel[0] = static_cast<int>(depth[x][z] / depth_max * 255);
				pixel[1] = static_cast<int>(depth[x][z] / depth_max * 255);
				pixel[2] = static_cast<int>(depth[x][z] / depth_max * 255);
			}
			else
			{

				unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, z, 0));
				pixel[0] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);
				pixel[1] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);
				pixel[2] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);

			}
		}
	}
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(path.c_str());
	writer->SetInputData(image);
	writer->Write();

	return sm;
}


void checkfloder(std::string path)
{
	// Attempt to create the directory.
	if (CreateDirectoryA(path.c_str(), NULL)) {
		printf("Directory created successfully.\n");
	}
	else {
		// If the directory could not be created, print an error message.
		if (GetLastError() == ERROR_ALREADY_EXISTS) {
			printf("Directory already exists.\n");
		}
		else {
			printf("Failed to create directory. Error code: %ld\n", GetLastError());
		}
	}
}

std::string output_folder_path = "F:/.tmp/output/";
int cur_folder = 0;

void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSmartPointer<vtkCoordinate> coordinate = vtkSmartPointer<vtkCoordinate>::New();
	coordinate->SetCoordinateSystemToDisplay();
	coordinate->SetValue(pipeline->RenderWindow->GetSize()[0] / 2, pipeline->RenderWindow->GetSize()[1] / 2, 0);
	double3 dir(coordinate->GetComputedWorldValue(pipeline->Renderer));
	dir = dir - double3(pipeline->Renderer->GetActiveCamera()->GetPosition());

	checkfloder(output_folder_path + std::to_string(cur_folder));

	/*CGAL::IO::write_PLY("rotated_toothmesh.ply", writePNG(toothmesh0, dir, output_folder_path + std::to_string(cur_folder) + "/toothmesh.png", x_dim, z_dim, centroid));
	CGAL::IO::write_PLY("rotated_toothmesh1.ply", writePNG(toothmesh1, dir, output_folder_path + std::to_string(cur_folder) + "/toothmesh1.png", x_dim, z_dim, centroid));
	CGAL::IO::write_PLY("rotated_crownmesh.ply", writePNG(crownmesh, dir, output_folder_path + std::to_string(cur_folder) + "/crownmesh.png", x_dim, z_dim, centroid));
	CGAL::IO::write_PLY("rotated_bitemesh.ply", writePNG(bitemesh, dir, output_folder_path + std::to_string(cur_folder) + "/bitemesh.png", x_dim, z_dim, centroid));*/
}


int main()
{
	std::string selected_folder_path = select_folder();

	std::cout << "selected_folder_path: " << selected_folder_path << std::endl;
	int num_folders = 35;
	x_dim = 1024;
	z_dim = 1024;
	centroid = double3(0, 0, 0);

	for (int i = 1; i <= num_folders; i++)
	{
		cur_folder = i;
		std::string folder_path;
		if(i<10)
			folder_path = selected_folder_path + "/000" + std::to_string(i);
		else
			folder_path = selected_folder_path + "/00" + std::to_string(i);


		std::string m1 = folder_path + "/m2.stl";
		std::string crown = folder_path + "/c.stl";
		std::string bite = folder_path + "/b.stl";
		std::string m0 = folder_path + "/m1.ply";

		pipeline = new vtkRenderPipeline();
		toothmesh0.clear();
		toothmesh1.clear();
		crownmesh.clear();
		bitemesh.clear();

		bool success=CGAL::IO::read_polygon_mesh(m0, toothmesh0);
		std::cout<< success <<toothmesh0.number_of_vertices()<<std::endl;
		CGAL::IO::read_polygon_mesh(m1, toothmesh1);
		CGAL::IO::read_polygon_mesh(crown, crownmesh);
		CGAL::IO::read_polygon_mesh(bite, bitemesh);
		RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(toothmesh0), pipeline->Renderer, 1, 1, 1, 1);
		RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(toothmesh1), pipeline->Renderer, 1, 1, 1, 1);
		
		vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
		axes->SetTotalLength(10.0, 10.0, 10.0);

		vtkSmartPointer<vtkTextProperty> text_prop = vtkSmartPointer<vtkTextProperty>::New();
		text_prop->SetFontSize(1);

		axes->GetXAxisCaptionActor2D()->SetCaptionTextProperty(text_prop);
		axes->GetYAxisCaptionActor2D()->SetCaptionTextProperty(text_prop);
		axes->GetZAxisCaptionActor2D()->SetCaptionTextProperty(text_prop);

		pipeline->Renderer->AddActor(axes);

		pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
		pipeline->Renderer->ResetCamera();
		pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
		pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
		pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);

		pipeline->RenderWindowInteractor->Start();
		delete pipeline;
	}
}