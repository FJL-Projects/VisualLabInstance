#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"

#include <CGAL/mesh_segmentation.h>
#define ENABLE_TIMER_H

vtkRenderPipeline* pipeline;
using Facet_double_map = SurfaceMesh::Property_map<face_descriptor, double>;
// create a property-map for segment-ids
using Facet_int_map = SurfaceMesh::Property_map<face_descriptor, std::size_t>;
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

int SurfaceMeshToPolyDataWithColor(
    const SurfaceMesh& mesh,
    vtkSmartPointer<vtkPolyData>& polydata
)
{
    //polyData->Initialize();
    typedef SurfaceMesh::Property_map<vertex_descriptor, CGAL::Color> Color_map;
    polydata = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkFloatArray> curvature_color_array = vtkSmartPointer<vtkFloatArray>::New();
    curvature_color_array->SetName("curvature");
    Color_map color_map;
    bool found;
    boost::tie(color_map, found) = mesh.property_map<vertex_descriptor, CGAL::Color>("v:color");

    if (!found)
    {
        std::cout << "Property map v:color not found!" << std::endl;

        vtkSmartPointer<vtkPoints> VTKPoints = vtkSmartPointer<vtkPoints>::New();
        VTKPoints->SetNumberOfPoints(mesh.vertices().size());

        for (int i = 0; i < mesh.vertices().size(); i++)
        {
            vertex_descriptor vd = vertex_descriptor(i);
            const Kernel::Point_3& point = mesh.point(vd);
            VTKPoints->SetPoint(i, point.x(), point.y(), point.z());
        }

        vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
        connectivity->SetNumberOfComponents(4);
        connectivity->SetNumberOfTuples(mesh.faces().size());

        //#pragma omp parallel for
        for (int j = 0; j < mesh.faces().size(); j++)
        {
            int ids[3] = { 0 };
            int index = 0;
            for (const halfedge_descriptor& h : mesh.halfedges_around_face(mesh.halfedge(face_descriptor(j))))
            {
                ids[index] = mesh.target(h).idx();
                index++;
            }
            if (index > 3)
                std::cerr << "SurfaceMeshToPolyData Error";

            connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
        }

        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        cells->SetCells(mesh.faces().size(), connectivity);

        polydata->SetPoints(VTKPoints);
        polydata->SetPolys(cells);

    }
    else
    {
        // Property map v:color found!
        vtkSmartPointer<vtkPoints> VTKPoints = vtkSmartPointer<vtkPoints>::New();
        VTKPoints->SetNumberOfPoints(mesh.vertices().size());

        for (int i = 0; i < mesh.vertices().size(); i++)
        {
            vertex_descriptor vd = vertex_descriptor(i);
            const Kernel::Point_3& point = mesh.point(vd);
            VTKPoints->SetPoint(i, point.x(), point.y(), point.z());
            if (color_map[vd] == CGAL::Color(255, 0, 0))
            {
                curvature_color_array->InsertNextValue(0.0);
            }
            else if (color_map[vd] == CGAL::Color(0, 255, 0))
            {
                curvature_color_array->InsertNextValue(1.0);
            }
            else
            {
                curvature_color_array->InsertNextValue(2.0);
            }
        }

        vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
        connectivity->SetNumberOfComponents(4);
        connectivity->SetNumberOfTuples(mesh.faces().size());

        //#pragma omp parallel for
        for (int j = 0; j < mesh.faces().size(); j++)
        {
            int ids[3] = { 0 };
            int index = 0;
            for (const halfedge_descriptor& h : mesh.halfedges_around_face(mesh.halfedge(face_descriptor(j))))
            {
                ids[index] = mesh.target(h).idx();
                index++;
            }
            if (index > 3)
                std::cerr << "SurfaceMeshToPolyData Error";

            connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
        }

        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        cells->SetCells(mesh.faces().size(), connectivity);

        polydata->SetPoints(VTKPoints);
        polydata->SetPolys(cells);
        polydata->GetPointData()->SetScalars(curvature_color_array);
    }

    return 0;
}

std::pair<vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkFloatArray> > SurfaceMeshToPolyDataWithColor(const SurfaceMesh& mesh)
{
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> VTKPoints = vtkSmartPointer<vtkPoints>::New();
    VTKPoints->SetNumberOfPoints(mesh.vertices().size());
    const Facet_int_map segment_property_map = mesh.property_map<face_descriptor, std::size_t>("f:sid").first;

    vtkSmartPointer<vtkFloatArray> color_array = vtkSmartPointer<vtkFloatArray>::New();
    color_array->SetName("color");

    for (size_t i = 0; i < mesh.vertices().size(); i++)
    {
		vertex_descriptor vd = vertex_descriptor(i);
		const Kernel::Point_3& point = mesh.point(vd);
		VTKPoints->SetPoint(i, point.x(), point.y(), point.z());

        const halfedge_descriptor& hd = *(CGAL::halfedges_around_source(vd, mesh).begin());
        const face_descriptor& fd = mesh.face(hd);

        const std::size_t segment_id = segment_property_map[fd];
        color_array->InsertNextValue(static_cast<double>(segment_id));
	}

    vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
    connectivity->SetNumberOfComponents(4);
    connectivity->SetNumberOfTuples(mesh.faces().size());

    for (size_t j = 0; j < mesh.faces().size(); j++)
    {
        int ids[3] = { 0 };
        int index = 0;
        for (const halfedge_descriptor& h : mesh.halfedges_around_face(mesh.halfedge(face_descriptor(j))))
        {
            ids[index] = mesh.target(h).idx();
            index++;
        }
        if (index > 3)
            std::cerr << "SurfaceMeshToPolyData Error";

        connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
    }

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(mesh.faces().size(), connectivity);

    polydata->SetPoints(VTKPoints);
    polydata->SetPolys(cells);
    polydata->GetPointData()->SetScalars(color_array);

    return { polydata, color_array };
}
int main()
{
    using namespace Eigen;

    pipeline = new vtkRenderPipeline();

    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName("data/sample_lower_left_upsampled.vtp");
    reader->Update();
    vtkSmartPointer<vtkPolyData> arch_pd = reader->GetOutput();

    SurfaceMesh arch_sm;
    CGAL::IO::read_VTP("data/sample_lower_left_upsampled.vtp", arch_sm);

    std::vector<halfedge_descriptor> border_halfedges;

    for (auto& h : arch_sm.halfedges())
    {
        if (arch_sm.is_border(h))
        {
			border_halfedges.push_back(h);
		}
    }

    std::vector<face_descriptor> new_faces;
    CGAL::Polygon_mesh_processing::triangulate_hole(
        arch_sm,
        *(border_halfedges.begin()),
        std::back_inserter(new_faces)
    );

    RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(arch_sm), pipeline->Renderer);
    //Facet_double_map sdf_property_map;
    //sdf_property_map = arch_sm.add_property_map<face_descriptor, double>("f:sdf").first;
    //CGAL::sdf_values(arch_sm, sdf_property_map);

    //// create a property-map for segment-ids
    //Facet_int_map segment_property_map = arch_sm.add_property_map<face_descriptor, std::size_t>("f:sid").first;;
    //// segment the mesh using default parameters for number of levels, and smoothing lambda
    //// Any other scalar values can be used instead of using SDF values computed using the CGAL function
    //std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(arch_sm, sdf_property_map, segment_property_map);
    //typedef CGAL::Face_filtered_graph<SurfaceMesh> Filtered_graph;
    ////print area of each segment and then put it in a Mesh and print it in an OFF file
    //Filtered_graph segment_mesh(arch_sm);
    //for (std::size_t id = 0; id < number_of_segments; ++id)
    //{
    //    segment_mesh.set_selected_faces(id, segment_property_map);
    //    std::cout << "Segment " << id << "'s area is : " << CGAL::Polygon_mesh_processing::area(segment_mesh) << std::endl;
    //    SurfaceMesh out;
    //    CGAL::copy_face_graph(segment_mesh, out);
    //    std::ostringstream oss;
    //    oss << "data/Segment_" << id << ".off";
    //    std::ofstream os(oss.str().data());
    //    os << out;
    //}





 //   //RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(arch_sm), pipeline->Renderer);
 //   Facet_double_map sdf_property_map;
 //   sdf_property_map = arch_sm.add_property_map<face_descriptor, double>("f:sdf").first;
 //   // compute SDF values
 //   // We can't use default parameters for number of rays, and cone angle
 //   // and the postprocessing
 //   CGAL::sdf_values(arch_sm, sdf_property_map, 2.0 / 3.0 * CGAL_PI, 40, true);
 //
 //   Facet_int_map segment_property_map = arch_sm.add_property_map<face_descriptor, std::size_t>("f:sid").first;
 //   // segment the mesh using default parameters for number of levels, and smoothing lambda
 //   // Any other scalar values can be used instead of using SDF values computed using the CGAL function
 //   //std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(arch_sm, sdf_property_map, segment_property_map);

 //   //std::cout << "Number of segments: " << number_of_segments << std::endl;

 //   const std::size_t number_of_clusters = 8;       // use 4 clusters in soft clustering
 //   const double smoothing_lambda = 0.3;  // importance of surface features, suggested to be in-between [0,1]
 //   // Note that we can use the same SDF values (sdf_property_map) over and over again for segmentation.
 //   // This feature is relevant for segmenting the mesh several times with different parameters.
 //   CGAL::segmentation_from_sdf_values(arch_sm, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);

 //   std::set<std::size_t> segment_ids;

 //   for (face_descriptor fd : arch_sm.faces())
 //   {
	//	segment_ids.insert(segment_property_map[fd]);
	//}

 //   std::cout << "Segment count: " << segment_ids.size() << std::endl;

 //   vtkSmartPointer<vtkPolyData> colored_arch_pd;
 //   vtkSmartPointer<vtkFloatArray> color_array;
 //   std::tie(colored_arch_pd, color_array) = SurfaceMeshToPolyDataWithColor(arch_sm);

 //   vtkSmartPointer<vtkLookupTable> color_lut = vtkSmartPointer<vtkLookupTable>::New();
 //   color_lut->SetNumberOfColors(segment_ids.size()); // Set the number of colors
 //   color_lut->SetHueRange(1.0, 0.0); // Green to Red
 //   color_lut->Build(); // Generate lookup table

 //   //color_lut->SetTableValue(color_lut->GetIndex(0), 1.0, 0.0, 1.0, 1.0);  
 //   color_lut->SetTableValue(color_lut->GetIndex(-1), 0 / 255.0, 238 / 255.0, 238 / 255.0, 1.0);  


 //   int num = color_array->GetNumberOfValues();
 //   std::cout << "Number of values: " << num << std::endl;

 //   vtkSmartPointer<vtkPolyDataMapper> colored_arch_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
 //   colored_arch_mapper->SetInputData(colored_arch_pd);
 //   colored_arch_mapper->SetScalarRange(-1, segment_ids.size());
 //   colored_arch_mapper->SetLookupTable(color_lut);
 //   colored_arch_mapper->SetInterpolateScalarsBeforeMapping(1);

 //   vtkSmartPointer<vtkActor> colored_arch_actor = vtkSmartPointer<vtkActor>::New();
 //   colored_arch_actor->SetMapper(colored_arch_mapper);
 //   colored_arch_actor->GetProperty()->SetOpacity(1);

 //   pipeline->Renderer->AddActor(colored_arch_actor);

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