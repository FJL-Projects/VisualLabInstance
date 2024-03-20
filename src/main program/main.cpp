#include"stdafx.h"
#include"vtkRenderPipeline.h"
#include"meshTransform.h"
#include"simpleRender.h"

#include <sstream>
#include <vtkSelection.h>
#include <vtkInformation.h>
#include <vtkSelectionNode.h>
#include <vtkCellArray.h>
#include <vtkLineSource.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include "vectorAlgorithm.h"
#include <vtkLookupTable.h>
#include <vtkFloatArray.h>

#include <memory>
#include "TeethWrapper.h"
#include "Timer.hpp"
vtkSmartPointer<vtkActor> PolyDataActor;
vtkSmartPointer<vtkActor> CutPolyDataActor;
vtkSmartPointer<vtkActor> TeethPolyDataActor;
vtkSmartPointer<vtkActor> RemeshedTeethPolyDataActor;
vtkSmartPointer<vtkActor> BitePolyDataActor;
vtkSmartPointer<vtkActor> OccludedPolyDataActor;


typedef Kernel::Segment_3 Segment_3;

typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Kernel::Ray_3>::Type> Ray_intersection;
typedef boost::optional<Tree::Intersection_and_primitive_id<Segment_3>::Type> Segment_intersection;
typedef Kernel::Ray_3 Ray;

std::unique_ptr<Tree> bite_tree_ptr;
Vector_3 projection_direction(-0.037686, 0.379984, -0.203274);

std::set<uint32_t> intersected_v_idx_set;
std::set<vertex_descriptor> intersected_v_set;

double squared_distance_limit = 64.0;
std::map<uint32_t, double> vertex_idx_closest_distance_map;
double max_distance = 2.0;

TeethWrapper teeth_wrapper;
double shrink_distance = 0.0;
bool rotate_flag = true;


SurfaceMesh mesh;
vtkRenderPipeline* pipeline;

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
    if (rotate_flag)
    {
        pipeline->Renderer->RemoveActor(TeethPolyDataActor);
        pipeline->Renderer->AddActor(RemeshedTeethPolyDataActor);
        pipeline->RenderWindow->Render();
        rotate_flag = false;
    }
    else
    {
        pipeline->Renderer->RemoveActor(RemeshedTeethPolyDataActor);
        pipeline->Renderer->AddActor(TeethPolyDataActor);
        pipeline->RenderWindow->Render();
        rotate_flag = true;
    }
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}

void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
    std::cout << "Switched Actor" << std::endl;
}

int count_adjacent_vertices(SurfaceMesh& mesh, vertex_descriptor v)
{
    int count = 0;
    auto range = mesh.vertices_around_target(mesh.halfedge(v));
    for (auto vcirc = range.begin(); vcirc != range.end(); ++vcirc)
    {
        ++count;
    }
    return count;
}

void smooth_mesh(SurfaceMesh& teeth_sm) {
    int vNum = teeth_sm.number_of_vertices();
    for (int i = 0; i < 15; ++i) {
        std::vector<Vector_3> vDeltas(vNum, Vector_3(0, 0, 0));
        int j = 0;
        for (auto& v : teeth_sm.vertices())
        {
            Vector_3& delta = vDeltas[j];
            Point_3 vf = teeth_sm.point(v);
            int vNerN = count_adjacent_vertices(teeth_sm, v);
            auto vcirc = teeth_sm.vertices_around_target(teeth_sm.halfedge(v)).begin(), done(teeth_sm.vertices_around_target(teeth_sm.halfedge(v)).end());
            do { delta = delta + ((teeth_sm.point(*vcirc) - vf) / vNerN); } while (++vcirc != done);
            ++j;
        }
        j = 0;
        for (auto& v : teeth_sm.vertices())
        {

            Vector_3& delta = vDeltas[j];
            Point_3 ans = teeth_sm.point(v);
            ans = ans + delta * 0.05;
            if (!teeth_sm.is_border(v))
                teeth_sm.point(v) = ans;
            ++j;
        }
    }
}
Vector_3 calculate_smooth_normal(SurfaceMesh& mesh, SurfaceMesh::Vertex_index v) {
    Vector_3 normal(0, 0, 0);
    auto circ = mesh.halfedge(v);
    auto done = circ;

    do {
        Vector_3 vec = CGAL::normal(
            mesh.point(mesh.source(circ)),
            mesh.point(mesh.target(circ)),
            mesh.point(mesh.target(mesh.next(circ)))
        );

        normal = normal + vec;
        circ = mesh.opposite(mesh.next(circ));

    } while (circ != done);

    normal = normal / std::sqrt(normal.squared_length());

    return normal;
}

void smooth_all_normal(SurfaceMesh& sm, VNMap& map)
{
    for (auto& v : sm.vertices())
    {
        map[v] = calculate_smooth_normal(sm, v);
    }
}

void GetIntersectedVertices(
    SurfaceMesh& teeth_sm,
    std::unique_ptr<Tree>& bite_tree_ptr,
    VLMap& distance_property_map,
    std::set<vertex_descriptor>& intersected_vertex_set,
    double max_distance
)
{
    Timer t("GetIntersectedVertices");
    intersected_vertex_set.clear();
    double distance = 0.0;

    for (auto vd : teeth_sm.vertices())
    {
        Point_3 source_point = teeth_sm.point(vd);

        Ray ray(source_point, projection_direction);
        std::vector<Primitive_id> intersections;
        bite_tree_ptr->all_intersected_primitives(ray, std::back_inserter(intersections));
        if (intersections.size() % 2 == 0)
        {
            continue;
        }
        auto intersection = bite_tree_ptr->first_intersection(ray);

        if (intersection)
        {
            if (const Point_3* intersected_point = boost::get<Point_3>(&(intersection->first)))
            {
                intersected_vertex_set.insert(vd);
            }
        }

        // This vertex is not intersected with the bite or newly added
        // Unchanged vertices shall have unchanged values other than 0.0
        if (distance_property_map[vd] == 0.0)
        {
            std::pair<Point_3, Primitive_id> pp = bite_tree_ptr->closest_point_and_primitive(source_point);
            distance = sqrt(CGAL::squared_distance(source_point, pp.first));
            max_distance = max_distance > distance ? max_distance : distance;
            //vertex_idx_closest_distance_map[vd.idx()] = distance;  // Old values of removed vertices will be overwritten
            distance_property_map[vd] = distance;
        }

    }
    std::cout << "intersected_vertex_set.size(): " << intersected_vertex_set.size() << std::endl;
}

void writePNG(SurfaceMesh& sm)
{
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
				pixel[0] = static_cast<int>(depth[x][z] / depth_max * 255); // 修改红色通道的值
				pixel[1] = static_cast<int>(depth[x][z] / depth_max * 255);   // 修改绿色通道的值
				pixel[2] = static_cast<int>(depth[x][z] / depth_max * 255);   // 修改蓝色通道的值
			}
			else
			{
				unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, z, 0));
				pixel[0] = 255 - static_cast<int>(depth[x][z] / depth_max * 255); // 修改红色通道的值
				pixel[1] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);   // 修改绿色通道的值
				pixel[2] = 255 - static_cast<int>(depth[x][z] / depth_max * 255);   // 修改蓝色通道的值
			}
		}
	}
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName("output.png");
	writer->SetInputData(image);
	writer->Write();
}

int main(int argc, char* argv[])
{
	pipeline = new vtkRenderPipeline();

	//CGAL::IO::read_polygon_mesh("data/test.stl",mesh);
	//RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(mesh), pipeline->Renderer,1,1,1,1);
	
        //vtkNew<vtkSTLReader> reader;
    //reader->SetFileName("arch.stl");
    //reader->Update();
    //PolyData = reader->GetOutput();
    for (int i = 0; i < argc; ++i)
    {
        std::cout << argv[i] << std::endl;
    }
    std::cout << "main" << std::endl;

    SurfaceMesh teeth_sm;

    vtkNew<vtkXMLPolyDataReader> vtp_reader;
    vtp_reader->SetFileName("data/arch_38057.vtp");
    vtp_reader->Update();
    vtkSmartPointer<vtkPolyData> PolyData = vtp_reader->GetOutput();
    std::cout << "arch.vtp" << std::endl;
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(PolyData);
    mapper->Update();
    PolyDataActor = vtkSmartPointer<vtkActor>::New();
    PolyDataActor->SetMapper(mapper);
    PolyDataActor->GetProperty()->SetColor(1, 1, 1);
    PolyDataActor->GetProperty()->SetAmbient(0.5);
    PolyDataActor->GetProperty()->SetSpecularPower(100);
    PolyDataActor->GetProperty()->SetSpecular(0.5);
    PolyDataActor->GetProperty()->SetDiffuse(0.5);
    PolyDataActor->GetProperty()->EdgeVisibilityOff();
    PolyDataActor->PickableOn();
    //Renderer->AddActor(PolyDataActor);

    SurfaceMesh arch_sm;
    PolyDataToSurfaceMesh(PolyData, arch_sm);

    std::ifstream visited_vd_reader("data/visited_vd_38057.txt");
    std::set<vertex_descriptor> arch_without_abutment_set;
    {
        std::string line;
        while (std::getline(visited_vd_reader, line))
        {
            std::istringstream iss(line);
            int x;
            if (!(iss >> x)) { break; } // error
            arch_without_abutment_set.insert(vertex_descriptor(x));
        }
    }

    std::cout << arch_without_abutment_set.size() << std::endl;
    // PolyDataToSurfaceMesh(PolyData, teeth);
    double thickness = 0.4;
    double offset = 0.2;
    double transform_radius = 0.3;



    vtkSmartPointer<vtkPLYReader> teeth_reader = vtkSmartPointer<vtkPLYReader>::New();
    teeth_reader->SetFileName("data/teeth_38057.ply");
    teeth_reader->Update();
    vtkSmartPointer<vtkPolyData> TeethPolydata = teeth_reader->GetOutput();

    PolyDataToSurfaceMesh(TeethPolydata, teeth_sm);

    teeth_wrapper.SetTeethSurfaceMesh(teeth_sm);
    teeth_wrapper.SetTeethActor(TeethPolyDataActor);
    teeth_wrapper.SetRenderer(pipeline->Renderer);
    teeth_wrapper.SetRenderWindow(pipeline->RenderWindow);
    teeth_wrapper.SetProjectionDirection(projection_direction);
    //teeth_wrapper.SetMaxDistance(max_distance);
    teeth_wrapper.SetOpacity(1);

    teeth_wrapper.SetSelectedId(6);
    teeth_wrapper.SetArchPolyData(PolyData);
    teeth_wrapper.SetArchSurfaceMesh(arch_sm);
    teeth_wrapper.SetArchWithoutAbutmentSet(arch_without_abutment_set);
    teeth_wrapper.ExtractAdjacentTeethWithoutAbutment();
    teeth_wrapper.ProximalShaving(0.08);
    SurfaceMesh remeshed_sm = teeth_wrapper.GetTeethSurfaceMesh();
    vtkSmartPointer<vtkPolyData> remeshed_polydata = CGAL_Surface_Mesh2VTK_PolyData(remeshed_sm);
    vtkSmartPointer<vtkPolyDataMapper> remeshed_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    remeshed_mapper->SetInputData(remeshed_polydata);
    remeshed_mapper->Update();
    RemeshedTeethPolyDataActor = vtkSmartPointer<vtkActor>::New();
    RemeshedTeethPolyDataActor->SetMapper(remeshed_mapper);
    RemeshedTeethPolyDataActor->GetProperty()->SetColor(0, 1, 1);
    RemeshedTeethPolyDataActor->GetProperty()->SetAmbient(0.5);
    RemeshedTeethPolyDataActor->GetProperty()->SetSpecularPower(100);
    RemeshedTeethPolyDataActor->GetProperty()->SetSpecular(0.5);
    RemeshedTeethPolyDataActor->GetProperty()->SetDiffuse(0.5);
    RemeshedTeethPolyDataActor->GetProperty()->SetOpacity(1.0);
    RemeshedTeethPolyDataActor->GetProperty()->EdgeVisibilityOff();
    //RemeshedTeethPolyDataActor->PickableOn();

    std::shared_ptr<SurfaceMesh> adjacent_teeth_sm = teeth_wrapper.GetAdjacentTeethSurfaceMesh();
    vtkSmartPointer<vtkPolyData> adjacent_teeth_polydata = CGAL_Surface_Mesh2VTK_PolyData(*adjacent_teeth_sm);

    vtkSmartPointer<vtkPolyDataMapper> adjacent_teeth_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    adjacent_teeth_mapper->SetInputData(adjacent_teeth_polydata);
    adjacent_teeth_mapper->Update();

    vtkSmartPointer<vtkActor> adjacent_teeth_actor = vtkSmartPointer<vtkActor>::New();
    adjacent_teeth_actor->SetMapper(adjacent_teeth_mapper);
    adjacent_teeth_actor->GetProperty()->SetColor(1, 1, 1);
    adjacent_teeth_actor->GetProperty()->SetAmbient(0.5);
    adjacent_teeth_actor->GetProperty()->SetSpecularPower(100);
    adjacent_teeth_actor->GetProperty()->SetSpecular(0.5);
    adjacent_teeth_actor->GetProperty()->SetDiffuse(0.5);
    adjacent_teeth_actor->GetProperty()->SetOpacity(1.0);
    adjacent_teeth_actor->GetProperty()->EdgeVisibilityOff();
    //Renderer->AddActor(adjacent_teeth_actor);

    vtkSmartPointer<vtkPolyDataMapper> teeth_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    teeth_mapper->SetInputData(TeethPolydata);
    teeth_mapper->Update();
    TeethPolyDataActor = vtkSmartPointer<vtkActor>::New();
    TeethPolyDataActor->SetMapper(teeth_mapper);
    TeethPolyDataActor->GetProperty()->SetColor(0, 1, 1);
    TeethPolyDataActor->GetProperty()->SetAmbient(0.5);
    TeethPolyDataActor->GetProperty()->SetSpecularPower(100);
    TeethPolyDataActor->GetProperty()->SetSpecular(0.5);
    TeethPolyDataActor->GetProperty()->SetDiffuse(0.5);
    TeethPolyDataActor->GetProperty()->SetOpacity(1.0);
    //TeethPolyDataActor->GetProperty()->EdgeVisibilityOn();
    //TeethPolyDataActor->PickableOn();
    pipeline->Renderer->AddActor(TeethPolyDataActor);

    //PolyData = TeethPolydata;
    //vtkSmartPointer<vtkPolyDataMapper> teeth_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    //teeth_mapper->SetInputData(TeethPolydata);
    //teeth_mapper->Update();
    //TeethPolyDataActor = vtkSmartPointer<vtkActor>::New();
    //TeethPolyDataActor->SetMapper(teeth_mapper);
    //TeethPolyDataActor->GetProperty()->SetColor(0, 1, 1);
    //TeethPolyDataActor->GetProperty()->SetAmbient(0.5);
    //TeethPolyDataActor->GetProperty()->SetSpecularPower(100);
    //TeethPolyDataActor->GetProperty()->SetSpecular(0.5);
    //TeethPolyDataActor->GetProperty()->SetDiffuse(0.5);
    //TeethPolyDataActor->GetProperty()->SetOpacity(0.8);
    ////TeethPolyDataActor->GetProperty()->EdgeVisibilityOn();
    //TeethPolyDataActor->PickableOn();
    //Renderer->AddActor(TeethPolyDataActor);

    pipeline->Renderer->SetBackground(0.41, 0.41, 0.41);
	pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
	pipeline->Renderer->ResetCamera();
	pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
	pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
	pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);

	pipeline->RenderWindowInteractor->Start();
}