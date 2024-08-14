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
#include <MRMesh/MRMeshBuilder.h>

#include <limits>
#include <cstdio>

vtkRenderPipeline* pipeline;
SurfaceMesh arch_sm;
SurfaceMesh rotated_toothmesh0;

namespace fs = std::filesystem;

bool disable_left_key = false;
std::string accessing_data_path;

std::string input_folder = "F:\\0x0000\\Model\\classified\\Nature Workshop 2\\exocad\\Crown&Coping";

std::string output_folder = "F:\\0x0000\\pointr_dataset\\output2\\";  // Be advised: ATTACH an ending '\\'. The path to save the output files.
fs::path output_folder_path = fs::path(output_folder);
fs::path current_input_folder;
fs::path current_output_folder;
std::string fdi_str;
int file_count = 0;
fs::path last_processed_path;

constexpr int RESOLUTION = 2048;  // The resolution of the depth image.
std::string file_name_stem;

BOOL CreateDirectoriesRecursively(const std::string& path);
BOOL CreateDirectoriesRecursively(const fs::path& fs_path);
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
auto GenerateLeadingZeroNumberStr(int number, size_t len = 4) -> std::string
{
	std::ostringstream stream;
	stream << std::setw(len) << std::setfill('0') << number;
	return stream.str();
}

auto ExtractNumberFromFileName(const std::string& filename) -> int
{
	std::stringstream ss(filename);
	std::string temp;
	int number = -1;

	while (ss >> temp)
	{
		for (char c : temp) 
		{
			if (std::isdigit(c))
			{
				ss.seekg(-static_cast<int>(temp.size()), std::ios::cur); 
				ss >> number;
				return number;
			}
		}
	}

	throw std::runtime_error("No number found in the filename");
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
	//std::cout << "deleting the file: " << accessing_data_path << std::endl;
	//std::remove(accessing_data_path.c_str());
}
void RightRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Released" << endl;
}

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	if (current_input_folder == last_processed_path)
	{
		std::cout << "Already processed. Skipped overwriting." << std::endl;
		return;
	}

	int output_num = 0;
	++file_count;
	current_output_folder = fs::path(output_folder) / fdi_str / GenerateLeadingZeroNumberStr(file_count, 4);  // "output_folder/11/0001/"
	std::cout << "Copied to " << current_output_folder << std::endl;

	CreateDirectoriesRecursively(current_output_folder);
	for (const auto& entry : fs::directory_iterator(current_input_folder))
	{
		const auto& path = entry.path();
		file_name_stem = path.stem().string();
		if (file_name_stem == "Antagonist")
		{
			fs::copy(entry, current_output_folder / std::string("bite.stl"), fs::copy_options::overwrite_existing);
			++output_num;
		}
		else if (file_name_stem.find("Jaw_scan") != std::string::npos)
		{
			fs::copy(entry, current_output_folder / std::string("jaw.stl"), fs::copy_options::overwrite_existing);
			++output_num;
		}
		else if (file_name_stem.find("Preparation") != std::string::npos)
		{
			fs::copy(entry, current_output_folder / std::string("abutment.stl"), fs::copy_options::overwrite_existing);
			++output_num;
		}
		else if (file_name_stem.find("Tooth_scan") != std::string::npos)
		{
			fs::copy(entry, current_output_folder / std::string("full_abutment.stl"), fs::copy_options::overwrite_existing);
			++output_num;
		}
		else if (file_name_stem.find("Restoration") != std::string::npos)
		{
			fs::copy(entry, current_output_folder / std::string("crown.stl"), fs::copy_options::overwrite_existing);
			++output_num;
		}
	}
	if (output_num != 5)
	{
		std::cout << "Files not copied successfully" << std::endl;
	}
	last_processed_path = current_input_folder;
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{

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


auto ExtractPolyDataByLabelSurfaceMesh(vtkSmartPointer<vtkPolyData> polydata) -> std::array<SurfaceMesh, 8>
{
	std::array<SurfaceMesh, 8> surface_meshes;

	vtkSmartPointer<vtkCellData> celldata = polydata->GetCellData();
	vtkSmartPointer<vtkDataArray> labels = celldata->GetArray("Label");

	std::array<std::map<vtkIdType, vertex_descriptor>, 8> point_maps;

	for (vtkIdType i = 0; i < polydata->GetNumberOfCells(); i++)
	{
		int label = static_cast<int>(labels->GetTuple1(i)) % 8;
		vtkCell* cell = polydata->GetCell(i);

		if (label >= 0 && label <= 7)
		{
			std::vector<vertex_descriptor> vertex_indices;

			for (int j = 0; j < cell->GetNumberOfPoints(); j++)
			{
				vtkIdType point_id = cell->GetPointId(j);
				if (point_maps[label].find(point_id) == point_maps[label].end())
				{
					double p[3];
					polydata->GetPoint(point_id, p);
					Point_3 point(p[0], p[1], p[2]);
					vertex_descriptor vertex = surface_meshes[label].add_vertex(point);
					point_maps[label][point_id] = vertex;
				}
				vertex_indices.push_back(point_maps[label][point_id]);
			}

			if (vertex_indices.size() == 3)
			{
				surface_meshes[label].add_face(vertex_indices[0], vertex_indices[1], vertex_indices[2]);
			}
		}
	}

	return surface_meshes;
}

void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}

int main(int argc, char* argv[])
{
	using namespace MR;

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

	input_folder = select_folder();
	output_folder = select_folder();

	file_count = 0;
	fdi_str = fs::path(input_folder).filename().string();
	current_input_folder = fs::path(input_folder);
	for (const auto& batches : fs::directory_iterator(current_input_folder))
	{
		//if (file_count > 30)
		//{
		//	break;
		//}

		// Create a new render pipeline and load the mesh files.
		pipeline = new vtkRenderPipeline();

		std::cout << batches << std::endl;
		current_input_folder = batches.path();
		for (const auto& each_batch : fs::directory_iterator(batches))
		{
			const auto& path = each_batch.path();
			file_name_stem = path.stem().string();
			if (file_name_stem.find("Restoration") != std::string::npos)
			{
				auto stl_reader = vtkSmartPointer<vtkSTLReader>::New();
				stl_reader->SetFileName(path.string().c_str());
				stl_reader->Update();

				auto polydata = stl_reader->GetOutput();
				RenderPolydata(polydata, pipeline->Renderer);
			}
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
	}
	delete pipeline;
}
