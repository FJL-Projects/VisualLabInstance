#include"stdafx.h"
#include"vtkRenderPipeline.h"
#include"meshTransform.h"
#include"simpleRender.h"

#include <opencv2/opencv.hpp>
#include <limits>

#include "TeethDataInitialization.h"

vtkRenderPipeline* pipeline;
SurfaceMesh arch_sm;
SurfaceMesh rotated_toothmesh0;

namespace fs = std::filesystem;

bool disable_left_key = false;

std::string output_folder_path = "D:\\Code\\VisualLabExperiment\\data\\output\\";  // Be advised: ATTACH an ending '\\'. The path to save the output files.
fs::path output_image_with_boxes_path = fs::path(output_folder_path) / "images" / "with_boxes";
fs::path output_image_path = fs::path(output_folder_path) / "images" / "test";
fs::path output_depth_image_path = fs::path(output_folder_path) / "images" / "depth";

fs::path output_label_path = fs::path(output_folder_path) / "labels" / "test";
int current_folder_num = 0; 
constexpr int RESOLUTION = 2048;  // The resolution of the depth image.
std::string file_name_stem;

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
 *   std::string num_str = GenerateLeadingZeroNumberStr(num);
 *   // num_str will be "0042"
 */
std::string GenerateLeadingZeroNumberStr(int number)
{
	std::ostringstream stream;
	stream << std::setw(1) << std::setfill('0') << number;
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
SurfaceMesh RotateMeshCopy(const SurfaceMesh& sm, const Eigen::Matrix3d& rotation_matrix)
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
std::string SelectFolder()
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
 *   BOOL result = CreateDirectoriesRecursively(path);
 *   // If result is TRUE, the directory structure is created successfully or already exists.
 */
BOOL CreateDirectoriesRecursively(const std::string& path) 
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
			if (!CreateDirectoriesRecursively(path.substr(0, slashIndex)))
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

BOOL CreateDirectoriesRecursively(const fs::path& fs_path)
{
	std::string path = fs_path.string();
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
			if (!CreateDirectoriesRecursively(path.substr(0, slashIndex)))
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


cv::Mat GenerateDepthImage(
	const SurfaceMesh& sm,
	const double& x_min,
	const double& y_min,
	const double& z_max,
	double max,
	const int resolution
) 
{
	using namespace cv;

	// Construct an AABB tree from the faces of the input surface mesh
	Tree tree(faces(sm).first, faces(sm).second, sm);
	Vector_3 negative_z_axis(0, 0, -1);

	// Create an OpenCV Mat to store the depth image
	Mat depth_image(resolution, resolution, CV_8UC3, Scalar(0, 0, 0));
	double depth_max = std::numeric_limits<double>::min();
	std::vector<std::vector<double>> depth(resolution, std::vector<double>(resolution, 0));
	double step = max / (resolution - 1);

	double x_pos = x_min;
	for (int x = 0; x < resolution; x++, x_pos += step) 
	{
		double y_pos = y_min;
		for (int y = 0; y < resolution; y++, y_pos += step) 
		{
			// Cast a ray from the camera position (x_pos, y_pos, z_max) towards the negative z-direction
			Ray_3 ray_query(Point_3(x_pos, y_pos, z_max), negative_z_axis);
			auto intersection = tree.first_intersection(ray_query);

			if (intersection) 
			{
				const Point_3* p = boost::get<Point_3>(&(intersection->first));
				if (p) 
				{
					depth[x][y] = z_max - p->z();
					depth_max = std::max(depth_max, depth[x][y]);
				}
			}
		}
	}

	uchar intensity;
	// Map the depth values to grayscale pixel intensities
	for (int x = 0; x < resolution; x++) 
	{
		for (int y = 0; y < resolution; y++) 
		{
			intensity = static_cast<uchar>((depth[x][y] / depth_max) * 255);
			if (!(depth[x][y] == 0))
			{
				intensity = 255 - intensity;
			}
			depth_image.at<Vec3b>(resolution - y - 1, x) = Vec3b(intensity, intensity, intensity);
		}
	}

	cv::cvtColor(depth_image, depth_image, COLOR_RGB2GRAY);
	return depth_image;
}

cv::Mat ConvertContourImage(const cv::Mat& depth_image)
{
	cv::Mat contour_image = depth_image;
	if (depth_image.type() != CV_8UC1) {
		cv::Mat temp;
		std::cout << "Converting to CV_8UC1\n";
		depth_image.convertTo(temp, CV_8UC1);  // 调整为合适的scale因子
		contour_image = temp;
	}

	adaptiveThreshold(contour_image, contour_image, 255, cv::ADAPTIVE_THRESH_MEAN_C,
		cv::THRESH_BINARY_INV, 31, 1.5);

	return contour_image;
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
void GenerateDepthImage(
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
	image->SetDimensions(RESOLUTION, RESOLUTION, 1);
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
	for (int x = 0; x < RESOLUTION; x++)
	{
		for (int y = 0; y < RESOLUTION; y++)
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
	if (!disable_left_key)
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

		rotated_toothmesh0 = RotateMeshCopy(arch_sm, transform_matrix_eigen);

		std::string output_folder_path_prefix = output_folder_path;
		fs::path output_folder_fspath(output_folder_path_prefix);
		CreateDirectoriesRecursively(output_folder_fspath);

		CGAL::IO::write_PLY(fs::path(output_folder_fspath / (file_name_stem + "_rotated.ply")).string(), rotated_toothmesh0);

		double x_min = std::numeric_limits<double>::max();
		double x_max = std::numeric_limits<double>::min();
		double y_min = std::numeric_limits<double>::max();
		double y_max = std::numeric_limits<double>::min();
		double rotated_toothmesh0_z_max = std::numeric_limits<double>::min();

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

		// Make sure the step is great enough to avoid missing points.
		double step = max / (RESOLUTION - 1);

		fs::path output_depth_path = fs::path(output_folder_fspath / (file_name_stem + ".png"));
		fs::path output_contour_path = fs::path(output_folder_fspath / (file_name_stem + "_contour.png"));

		cv::Mat depth_image = GenerateDepthImage(rotated_toothmesh0, x_min, y_min, rotated_toothmesh0_z_max, max, RESOLUTION);
		cv::imwrite(output_depth_path.string(), depth_image);
		std::cout << "Depth image saved to " << output_depth_path << std::endl;

		cv::Mat contour_image = ConvertContourImage(depth_image);
		cv::imwrite(output_contour_path.string(), contour_image);
	}
}

int main()
{
	using namespace cv;

	// Select the input folder containing the mesh files.
	/* The selected file hierarchy should be as follows:
	* 	 * - selected_folder_path
	* 	 *   - 1
	* 	 *     - xxx.stl
	* 	 *   - 2
	* 	 *     - yyy.stl
	* 	 *   ...
	* 	 *   - 35
	* 	 *     - m1.stl
	*/
	std::string selected_folder_path = SelectFolder();

	std::cout << "selected_folder_path: " << selected_folder_path << std::endl;
	int num_folders = 1;
	fs::path toothmesh_path;

	fs::path directory_path(selected_folder_path);

	CreateDirectoriesRecursively(output_image_with_boxes_path);
	CreateDirectoriesRecursively(output_image_path);
	CreateDirectoriesRecursively(output_label_path);
	CreateDirectoriesRecursively(output_depth_image_path);
	vtkSmartPointer<vtkPolyData> arch_pd;
	for (const auto& entry : fs::directory_iterator(directory_path))
	{
		// Create a new render pipeline and load the mesh files.
		pipeline = new vtkRenderPipeline();

		const auto& path = entry.path();
		file_name_stem = path.stem().string();
		if (path.extension() == ".vtp")
		{
			disable_left_key = true;

			toothmesh_path /= path;
			std::cout << toothmesh_path << std::endl;

			vtkSmartPointer<vtkXMLPolyDataReader> vtp_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
			vtp_reader->SetFileName(toothmesh_path.string().c_str());
			vtp_reader->Update();
			arch_pd = vtp_reader->GetOutput();

			// Clear the mesh data at each run.
			arch_sm.clear();
			arch_sm = PolyDataToSurfaceMesh(arch_pd);

			TeethDataInitialization teeth_data_initialization(arch_pd, arch_sm);
			teeth_data_initialization.SetRenderer(pipeline->Renderer);
			teeth_data_initialization.SetRenderWindow(pipeline->RenderWindow);
			teeth_data_initialization.Execute();

			auto& teeth_poly_data = teeth_data_initialization.m_labeledPolyData;

			std::vector<std::vector<double> > teeth_dimension_vec;
			for (size_t tooth_num = 1; tooth_num < 8; ++tooth_num)
			{
				auto tooth_pd = teeth_poly_data[tooth_num];
				auto tooth_vertices_number = tooth_pd->GetNumberOfPoints();
				//if (tooth_vertices_number < 300)
				//{
				//	continue;
				//}


				auto [x_min, x_max, y_min, y_max] = [](vtkSmartPointer<vtkPolyData> pd)
					{
						auto* bounding_box = pd->GetBounds();
						return std::tuple<double, double, double, double>(bounding_box[0], bounding_box[1], bounding_box[2], bounding_box[3]);
					}(tooth_pd);

					/*std::cout << "Tooth " << tooth_num << " has vertices: " << tooth_vertices_number << std::endl;
					std::cout << "Bounding box: x_min: " << x_min << " x_max: " << x_max << " y_min: " << y_min << " y_max: " << y_max << std::endl;*/
					teeth_dimension_vec.push_back({ x_min, x_max, y_min, y_max });
					RenderPolydata(tooth_pd, pipeline->Renderer, 1, 0, 0, 1);
			}

			double x_min = std::numeric_limits<double>::max();
			double x_max = std::numeric_limits<double>::min();
			double y_min = std::numeric_limits<double>::max();
			double y_max = std::numeric_limits<double>::min();
			double z_max = std::numeric_limits<double>::min();

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

			double max = get_dimension(arch_sm, z_max);

			auto contour_image_with_boxes_path = output_image_with_boxes_path / std::string(file_name_stem + "with_boxes.png");
			auto contour_image_path = output_image_path / std::string(file_name_stem + ".png");
			auto label_txt_path = output_label_path / std::string(file_name_stem + ".txt");
			auto depth_image_path = output_depth_image_path / std::string(file_name_stem + ".png");

			auto depth_image = GenerateDepthImage(arch_sm, x_min, y_min, z_max, max, RESOLUTION);
			cv::imwrite(depth_image_path.string(), depth_image);
			Mat contour_image = ConvertContourImage(depth_image);

			auto contour_image_with_boxes = contour_image.clone();
			std::ofstream label_txt(label_txt_path.string());
			double step = max / (RESOLUTION - 1);
			for (auto& box : teeth_dimension_vec)
			{
				if (box.size() == 4)
				{
					auto x_min_bound_pixel = (box[0] - x_min) / step;
					auto x_max_bound_pixel = (box[1] - x_min) / step;
					auto y_min_bound_pixel = RESOLUTION - (box[2] - y_min) / step;
					auto y_max_bound_pixel = RESOLUTION - (box[3] - y_min) / step;

					auto x_center = (x_min_bound_pixel + x_max_bound_pixel) / 2;
					auto y_center = (y_min_bound_pixel + y_max_bound_pixel) / 2;

					std::cout << "x_min: " << x_min_bound_pixel << " x_max: " << x_max_bound_pixel << " y_min: " << y_min_bound_pixel << " y_max: " << y_max_bound_pixel << " x_center: " << x_center << " y_center: " << y_center << std::endl;

					label_txt << "0 " << x_center / RESOLUTION << " " << y_center / RESOLUTION << " " << abs(x_max_bound_pixel - x_min_bound_pixel) / RESOLUTION << " " << abs(y_max_bound_pixel - y_min_bound_pixel) / RESOLUTION << std::endl;

					// 创建一个矩形（方框）
					cv::rectangle(
						contour_image_with_boxes,
						cv::Point(x_min_bound_pixel, y_min_bound_pixel),
						cv::Point(x_max_bound_pixel, y_max_bound_pixel),
						cv::Scalar(255, 0, 0), 2
					);
				}
			}
			label_txt.close();
			//// Render the axes.
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
			delete pipeline;
			cv::imwrite(contour_image_with_boxes_path.string(), contour_image_with_boxes);
			cv::imwrite(contour_image_path.string(), contour_image);
		}
		else if (path.extension() == ".stl")
		{
			toothmesh_path /= path;
			std::cout << toothmesh_path << std::endl;
			file_name_stem = toothmesh_path.stem().string();

			arch_sm.clear();
			CGAL::IO::read_STL(toothmesh_path.string(), arch_sm);
			arch_pd = CGAL_Surface_Mesh2VTK_PolyData(arch_sm);

			std::cout << "Arch has vertices: " << arch_sm.number_of_vertices() << std::endl;
			// Render the tooth mesh's arch and the abutment for the user to align.
			RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(arch_sm), pipeline->Renderer, 1, 1, 1, 1);

			//// Render the axes.
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

			//// Clean up the pipeline after each run.
			delete pipeline;
		}
		else
		{
			continue;
		}
	}
}