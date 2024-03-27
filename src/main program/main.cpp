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

std::string output_folder_path = "F:\\.tmp\\output\\";  // Be advised: ATTACH an ending '\\'. The path to save the output files.
int current_folder_num = 0; 
int resolution = 4096;  // The resolution of the depth image.

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

/**
 * @brief Rotate a copy of the input mesh using the given rotation matrix.
 *
 * This function creates a copy of the input SurfaceMesh and applies a rotation
 * to all its vertices using the provided rotation matrix. The rotated mesh is
 * then returned as a new SurfaceMesh object.
 *
 * @param sm The input SurfaceMesh to be rotated.
 * @param rotation_matrix The 3x3 rotation matrix to be applied to the mesh.
 * @return A new SurfaceMesh object containing the rotated mesh.
 *
 * @note The function creates a copy of the input mesh to avoid modifying the
 *       original mesh. The rotation is applied to each vertex of the mesh
 *       using Eigen library's matrix-vector multiplication.
 */
SurfaceMesh rotate_mesh_copy(const SurfaceMesh& sm, const Eigen::Matrix3d& rotation_matrix)
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

/**
 * @brief Open a folder selection dialog and return the selected folder path.
 *
 * This function displays a folder selection dialog using the Windows Shell API.
 * It allows the user to browse and select a folder. The function returns the
 * path of the selected folder as a string. If no folder is selected, an empty
 * string is returned.
 *
 * The function also reads and writes the last selected folder path to an INI file
 * named "last_path.ini". If the file exists, the last selected path is used as the
 * initial directory for the folder selection dialog. After the user selects a folder,
 * the selected path is written back to the INI file.
 *
 * @return The path of the selected folder, or an empty string if no folder is selected.
 *
 * @note This function uses the Windows Shell API (ShBrowseForFolder) to display the
 *       folder selection dialog. It requires the Windows-specific headers and libraries.
 *
 * @note The function uses the C++11 `<codecvt>` library for string conversion between
 *       UTF-8 and wide strings.
 *
 * @note The `BrowseCallbackProc` function is a callback function used by the folder
 *       selection dialog to set the initial directory. It is not shown in this code snippet.
 */
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

/**
 * @brief Create directories recursively.
 *
 * @param path The full path of the directory to be created.
 *
 * @return TRUE if the directory is successfully created or already exists, FALSE otherwise.
 *
 * @note This function uses Windows API (GetFileAttributesA, CreateDirectoryA) and requires Windows-specific headers.
 *
 * @note The function recursively creates parent directories if they don't exist.
 *
 * @note If the specified path already exists and is a file, the function returns FALSE.
 *
 * @example
 *   std::string path = "C:\\parent\\child1\\child2";
 *   BOOL result = create_directories_recursively(path);
 *   // If result is TRUE, the directory structure is created successfully or already exists.
 */
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

/**
 * @brief Generate a depth image from a surface mesh.
 *
 * This function generates a depth image by casting rays from a virtual camera position (x, y, z_max)
 * towards the negative z-direction and finding the intersection points with the given surface mesh.
 * The depth values are computed based on the distance between the camera position and the intersection points.
 * The generated depth image is then saved as a PNG file using the VTK library.
 *
 * @param path The file path to save the generated depth image.
 * @param sm The input surface mesh.
 * @param x_min The minimum x-coordinate of the bounding box.
 * @param y_min The minimum y-coordinate of the bounding box.
 * @param z_max The maximum z-coordinate of the bounding box (camera position).
 * @param step The step size for ray casting (determines the resolution of the depth image).
 *
 * @note This function uses the CGAL library for ray-mesh intersection tests.
 *       - It constructs an AABB tree (Tree) from the faces of the input surface mesh (sm).
 *       - The tree is used to efficiently find the intersection points between the rays and the mesh.
 *
 * @note The depth image is generated by iterating over each pixel in the image plane.
 *       - For each pixel, a ray is cast from the camera position (x, y, z_max) towards the negative z-direction.
 *       - The first_intersection function of the AABB tree is used to find the intersection point between the ray and the mesh.
 *       - If an intersection point is found, the depth value is computed as the difference between z_max and the z-coordinate of the intersection point.
 *       - If no intersection is found, the depth value is set to 0.
 *
 * @note The depth values are normalized and mapped to grayscale pixel intensities.
 *       - The maximum depth value (depth_max) is computed during the ray casting process.
 *       - Pixels with depth value 0 (no intersection) are set to black.
 *       - Pixels with non-zero depth values are assigned grayscale intensities based on the normalized depth value.
 *       - The pixel intensities are inverted (255 - intensity) for better visualization, so that closer objects appear brighter.
 *
 * @note The generated depth image is saved as a PNG file using the VTK library.
 *       - The VTK image data (vtkImageData) is created with the same dimensions as the depth image.
 *       - The depth values are mapped to grayscale pixel intensities and stored in the VTK image data.
 *       - The VTK PNG writer (vtkPNGWriter) is used to save the image data as a PNG file.
 */
void generate_depth_image(
	const std::string& path,
	const SurfaceMesh& sm,
	const double& x_min,
	const double& y_min,
	const double& z_max,
	const double& step
)
{
	// Construct an AABB tree from the faces of the input surface mesh
	Tree tree(faces(sm).first, faces(sm).second, sm);
	Vector_3 negative_z_axis(0, 0, -1);

	// Create a VTK image data to store the depth image
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
	image->SetDimensions(resolution, resolution, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	int dim[3];
	image->GetDimensions(dim);

	// Create a 2D vector to store the depth values for each pixel
	std::vector<std::vector<double>> depth(dim[0], std::vector<double>(dim[1], 0));
	double depth_max = std::numeric_limits<double>::min();
	double x_pos = x_min;

	// Iterate over each pixel in the image plane
	for (int x = 0; x < dim[0]; x++, x_pos += step)
	{
		double y_pos = y_min;
		for (int y = 0; y < dim[1]; y++, y_pos += step)
		{
			// Cast a ray from the camera position (x_pos, y_pos, z_max) towards the negative z-direction
			Ray_3 ray_query(Point_3(x_pos, y_pos, z_max), negative_z_axis);

			// Find the first intersection point between the ray and the mesh using the AABB tree
			auto intersection = tree.first_intersection(ray_query);

			const Point_3* p;
			if (intersection)
			{
				// If an intersection point is found, get the point coordinates
				p = boost::get<Point_3>(&(intersection->first));

				// Compute the depth value as the difference between z_max and the z-coordinate of the intersection point
				depth_max = std::max(depth_max, depth[x][y] = z_max - p->z());
			}
			else
			{
				// If no intersection is found, set the depth value to 0
				depth[x][y] = 0;
			}
		}
	}

	// Map the depth values to grayscale pixel intensities
	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			// Get the pointer to the pixel data
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));

			if (depth[x][y] == 0)
			{
				// If the depth value is 0 (no intersection), set the pixel to black
				pixel[0] = static_cast<int>(depth[x][y] / depth_max * 255);
				pixel[1] = static_cast<int>(depth[x][y] / depth_max * 255);
				pixel[2] = static_cast<int>(depth[x][y] / depth_max * 255);
			}
			else
			{
				// If the depth value is non-zero, compute the grayscale intensity based on the normalized depth value
				// Invert the intensity (255 - intensity) for better visualization
				pixel[0] = 255 - static_cast<int>(depth[x][y] / depth_max * 255);
				pixel[1] = 255 - static_cast<int>(depth[x][y] / depth_max * 255);
				pixel[2] = 255 - static_cast<int>(depth[x][y] / depth_max * 255);
			}
		}
	}

	// Create a VTK PNG writer and save the depth image as a PNG file
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
	// Select the input folder containing the mesh files.
	/* The selected file hierarchy should be as follows:
	* 	 * - selected_folder_path
	* 	 *   - 0001
	* 	 *     - m1.ply
	* 	 *     - m2.ply
	* 	 *     - c.ply
	* 	 *     - b.ply
	* 	 *   - 0002
	* 	 *     - m1.ply
	* 	 *     - m2.ply
	* 	 *     - c.ply
	* 	 *     - b.ply
	* 	 *   ...
	* 	 *   - 0035
	* 	 *     - m1.ply
	* 	 *     - m2.ply
	* 	 *     - c.ply
	* 	 *     - b.ply
	*/
	std::string selected_folder_path = select_folder();

	std::cout << "selected_folder_path: " << selected_folder_path << std::endl;
	int num_folders = 35;

	// Iterate over each folder in the selected directory.
	for (int i = 1; i <= num_folders; i++)
	{
		current_folder_num = i;
		std::string folder_path = selected_folder_path + "\\" + generate_leading_zero_number_str(i);

		// CGAL's read_STL may fail to read some shattered files.
		std::string toothmesh_path = folder_path + "\\m1.ply";
		std::string toothmesh1_path = folder_path + "\\m2.ply";
		std::string crown_path = folder_path + "\\c.ply";
		std::string bite_path = folder_path + "\\b.ply";

		// Create a new render pipeline and load the mesh files.
		pipeline = new vtkRenderPipeline();
		// Clear the mesh data at each run.
		toothmesh0.clear();
		toothmesh1.clear();
		crownmesh.clear();
		bitemesh.clear();

		CGAL::IO::read_polygon_mesh(toothmesh_path, toothmesh0);
		CGAL::IO::read_polygon_mesh(toothmesh1_path, toothmesh1);
		CGAL::IO::read_polygon_mesh(crown_path, crownmesh);
		CGAL::IO::read_polygon_mesh(bite_path, bitemesh);
		
		// Render the tooth mesh's arch and the abutment for the user to align.
		RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(toothmesh0), pipeline->Renderer, 1, 1, 1, 1);
		RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(toothmesh1), pipeline->Renderer, 1, 1, 1, 1);

		// Render the axes.
		vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
		axes->SetTotalLength(10.0, 10.0, 10.0);

		vtkSmartPointer<vtkTextProperty> text_prop = vtkSmartPointer<vtkTextProperty>::New();
		text_prop->SetFontSize(1);

		axes->GetXAxisCaptionActor2D()->SetCaptionTextProperty(text_prop);
		axes->GetYAxisCaptionActor2D()->SetCaptionTextProperty(text_prop);
		axes->GetZAxisCaptionActor2D()->SetCaptionTextProperty(text_prop);

		pipeline->Renderer->AddActor(axes);

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