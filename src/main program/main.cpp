#include"stdafx.h"
#include"vtkRenderPipeline.h"
#include"meshTransform.h"
#include"simpleRender.h"

vtkRenderPipeline* pipeline;
SurfaceMesh toothmesh0;
SurfaceMesh toothmesh1;
SurfaceMesh crownmesh;
SurfaceMesh bitemesh;
SurfaceMesh rotated_toothmesh0;
SurfaceMesh rotated_toothmesh1;
SurfaceMesh rotated_crownmesh;
SurfaceMesh rotated_bitemesh;

std::string output_folder_path = "F:\\.tmp\\output\\";  // The path to save the output files.
int current_folder_num = 0; 
int resolution = 4096;  // The resolution of the depth image.

std::string generate_leading_zero_number_str(int number)
{
	std::ostringstream stream;
	stream << std::setw(4) << std::setfill('0') << number;
	return stream.str();
}

SurfaceMesh rotate_mesh_copy(SurfaceMesh& sm, Eigen::Matrix3d& rotation_matrix)
{
	SurfaceMesh sm_copy = sm;
	for (auto& v : sm_copy.vertices())
	{
		Point_3& p = sm_copy.point(v);
		Eigen::Vector3d vec(p.x(), p.y(), p.z());
		vec = rotation_matrix * vec;
		p = Point_3(vec.x(), vec.y(), vec.z());
	}
	return sm_copy;
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

int CALLBACK BrowseCallbackProc(HWND hwnd, UINT uMsg, LPARAM lParam, LPARAM lpData)
{
	if (uMsg == BFFM_INITIALIZED)
	{
		// lpData is the lParam value passed to SHBrowseForFolder
		SendMessage(hwnd, BFFM_SETSELECTION, TRUE, lpData);
	}
	return 0;
}

std::string select_folder()
{
	std::string last_path;
	read_ini_file("last_path.ini", last_path);

	BROWSEINFO bi = { 0 };
	bi.lpszTitle = L"Browse for folder...";
	bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
	bi.lpfn = BrowseCallbackProc;

	std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
	std::wstring wide_last_path = converter.from_bytes(last_path);

	bi.lParam = reinterpret_cast<LPARAM>(wide_last_path.c_str());

	LPITEMIDLIST pidl = SHBrowseForFolder(&bi);

	if (pidl != 0)
	{
		wchar_t path[MAX_PATH];
		if (SHGetPathFromIDList(pidl, path))
		{
			IMalloc* imalloc = 0;
			if (SUCCEEDED(SHGetMalloc(&imalloc)))
			{
				imalloc->Free(pidl);
				imalloc->Release();
			}
			std::string selected_folder = converter.to_bytes(path);
			write_ini_file("last_path.ini", selected_folder);
			return selected_folder;
		}
	}

	return std::string();
}

BOOL create_directories_recursively(const std::string& path) 
{
	DWORD dwAttrib = GetFileAttributesA(path.c_str());

	// Check if the path exists and is not a file
	if (dwAttrib != INVALID_FILE_ATTRIBUTES &&
		!(dwAttrib & FILE_ATTRIBUTE_DIRECTORY)) 
	{
		return FALSE;
	}

	// Try to create the directory
	if (dwAttrib == INVALID_FILE_ATTRIBUTES) 
	{
		// Recursively create the parent directory
		size_t slashIndex = path.find_last_of("/\\");
		if (slashIndex != std::string::npos) {
			if (!create_directories_recursively(path.substr(0, slashIndex)))
			{
				return FALSE;
			}
		}

		// Create the last directory
		if (!CreateDirectoryA(path.c_str(), NULL)) 
		{
			return FALSE;
		}
	}

	return TRUE;
}

void generate_depth_image(
	const std::string& path,
	const SurfaceMesh& sm,
	const double& x_min,
	const double& y_min,
	const double& z_max, 
	const double& step
)
{
	Tree tree(faces(sm).first, faces(sm).second, sm);
	Vector_3 negative_z_axis(0, 0, -1); 
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
	image->SetDimensions(resolution, resolution, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	int dim[3];
	image->GetDimensions(dim);

	std::vector<std::vector<double>> depth(dim[0], std::vector<double>(dim[1], 0));
	double depth_max = std::numeric_limits<double>::min();
	double x_pos = x_min;
	for (int x = 0; x < dim[0]; x++, x_pos += step)
	{
		double y_pos = y_min; 
		for (int y = 0; y < dim[1]; y++, y_pos += step) 
		{
			Ray_3 ray_query(Point_3(x_pos, y_pos, z_max), negative_z_axis); 
			auto intersection = tree.first_intersection(ray_query);
			const Point_3* p;
			if (intersection)
			{
				p = boost::get<Point_3>(&(intersection->first));
				depth_max = std::max(depth_max, depth[x][y] = z_max - p->z()); 
			}
			else
			{
				depth[x][y] = 0;
			}
		}
	}

	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			if (depth[x][y] == 0) 
			{
				pixel[0] = static_cast<int>(depth[x][y] / depth_max * 255); 
				pixel[1] = static_cast<int>(depth[x][y] / depth_max * 255);
				pixel[2] = static_cast<int>(depth[x][y] / depth_max * 255);
			}
			else
			{
				pixel[0] = 255 - static_cast<int>(depth[x][y] / depth_max * 255); 
				pixel[1] = 255 - static_cast<int>(depth[x][y] / depth_max * 255);
				pixel[2] = 255 - static_cast<int>(depth[x][y] / depth_max * 255);
			}
		}
	}
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(path.c_str());
	writer->SetInputData(image);
	writer->Write();
	std::cout << "Depth image saved to " << path << std::endl;
}

void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSmartPointer<vtkCamera> camera = pipeline->Renderer->GetActiveCamera();
	auto transform_matrix = camera->GetModelViewTransformMatrix();

	// Do not use transform_matrix directly, it involves translation. 
	// Use Eigen::Matrix3d transform_matrix_eigen instead.
	Eigen::Matrix3d transform_matrix_eigen;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			transform_matrix_eigen(i, j) = transform_matrix->GetElement(i, j);
		}
	}

	rotated_toothmesh0 = rotate_mesh_copy(toothmesh0, transform_matrix_eigen);
	rotated_toothmesh1 = rotate_mesh_copy(toothmesh1, transform_matrix_eigen);
	rotated_crownmesh = rotate_mesh_copy(crownmesh, transform_matrix_eigen);
	rotated_bitemesh = rotate_mesh_copy(bitemesh, transform_matrix_eigen);

	std::string output_folder_path_prefix = output_folder_path + generate_leading_zero_number_str(current_folder_num);
	create_directories_recursively(output_folder_path_prefix);

	CGAL::IO::write_PLY(output_folder_path_prefix + "\\toothmesh.ply", rotated_toothmesh0);
	CGAL::IO::write_PLY(output_folder_path_prefix + "\\toothmesh1.ply", rotated_toothmesh1);
	CGAL::IO::write_PLY(output_folder_path_prefix + "\\crownmesh.ply", rotated_crownmesh);
	CGAL::IO::write_PLY(output_folder_path_prefix + "\\bitemesh.ply", rotated_bitemesh);


	double x_min = std::numeric_limits<double>::max();
	double x_max = std::numeric_limits<double>::min();
	double y_min = std::numeric_limits<double>::max();
	double y_max = std::numeric_limits<double>::min();
	double rotated_toothmesh0_z_max = std::numeric_limits<double>::min();
	double rotated_toothmesh1_z_max = std::numeric_limits<double>::min();
	double rotated_crownmesh_z_max = std::numeric_limits<double>::min();
	double rotated_bitemesh_z_max = std::numeric_limits<double>::min();

	auto get_dimension = [&x_min, &x_max, &y_min, &y_max](SurfaceMesh& sm, double& z_max)
		{
			for (auto v : sm.vertices())
			{
				Point_3 p = sm.point(v);
				x_min = std::min(x_min, p.x());
				x_max = std::max(x_max, p.x());
				y_min = std::min(y_min, p.y());
				y_max = std::max(y_max, p.y());
				z_max = std::max(z_max, p.z());
			}
			double max;
			if ((x_max - x_min) > (y_max - y_min))
				max = x_max - x_min;
			else
				max = y_max - y_min;
			return max;
		};

	double max = get_dimension(rotated_toothmesh0, rotated_toothmesh0_z_max);
	max = std::max(max, get_dimension(rotated_toothmesh1, rotated_toothmesh1_z_max));
	max = std::max(max, get_dimension(rotated_crownmesh, rotated_crownmesh_z_max));
	max = std::max(max, get_dimension(rotated_bitemesh, rotated_bitemesh_z_max));

	// Make sure the step is great enough to avoid missing points.
	double step = max / (resolution - 1);

	generate_depth_image(output_folder_path_prefix + "\\toothmesh.png", rotated_toothmesh0, x_min, y_min, rotated_toothmesh0_z_max, step);
	generate_depth_image(output_folder_path_prefix + "\\toothmesh1.png", rotated_toothmesh1, x_min, y_min, rotated_toothmesh1_z_max, step);
	generate_depth_image(output_folder_path_prefix + "\\crownmesh.png", rotated_crownmesh, x_min, y_min, rotated_crownmesh_z_max, step);
	generate_depth_image(output_folder_path_prefix + "\\bitemesh.png", rotated_bitemesh, x_min, y_min, rotated_bitemesh_z_max, step);
}

int main()
{
	std::string selected_folder_path = select_folder();

	std::cout << "selected_folder_path: " << selected_folder_path << std::endl;
	int num_folders = 35;

	for (int i = 1; i <= num_folders; i++)
	{
		current_folder_num = i;
		std::string folder_path = selected_folder_path + "\\" + generate_leading_zero_number_str(i);

		// CGAL's read_STL may fail to read some shattered files.
		std::string toothmesh_path = folder_path + "\\m1.ply";
		std::string toothmesh1_path = folder_path + "\\m2.ply";
		std::string crown_path = folder_path + "\\c.ply";
		std::string bite_path = folder_path + "\\b.ply";

		pipeline = new vtkRenderPipeline();
		toothmesh0.clear();
		toothmesh1.clear();
		crownmesh.clear();
		bitemesh.clear();

		CGAL::IO::read_polygon_mesh(toothmesh_path, toothmesh0);
		CGAL::IO::read_polygon_mesh(toothmesh1_path, toothmesh1);
		CGAL::IO::read_polygon_mesh(crown_path, crownmesh);
		CGAL::IO::read_polygon_mesh(bite_path, bitemesh);
		
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
		pipeline->addObserver(vtkCommand::RightButtonPressEvent, RightPress);
		pipeline->addObserver(vtkCommand::RightButtonReleaseEvent, RightRelease);

		pipeline->RenderWindowInteractor->Start();
		delete pipeline;
	}
}