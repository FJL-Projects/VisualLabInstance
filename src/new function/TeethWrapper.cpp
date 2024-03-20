#include"stdafx.h"
#include "TeethWrapper.h"

/**
 * @brief Default constructor for TeethWrapper
 * Initializes all member pointers to nullptr
 */
TeethWrapper::TeethWrapper() :
    m_teeth_sm(nullptr),
    m_bite_sm(nullptr),
    m_cut_sm(nullptr),
    m_bite_tree_ptr(nullptr),
    m_cut_tree_ptr(nullptr),
    m_teeth_halfedge_tree_ptr(nullptr),
    m_teeth_slicer(nullptr)
{}

/**
 * @brief Sets the teeth surface mesh by making a deep copy of the input surface mesh
 * @param sm The surface mesh to be set
 */
void TeethWrapper::SetTeethSurfaceMesh(SurfaceMesh& sm)
{
    //deep copy
    m_teeth_sm = new SurfaceMesh(sm);
}

/**
 * @brief Gets the teeth surface mesh
 * @return A reference to the teeth surface mesh
 * @throws An assertion error if the teeth surface mesh is not set
 */
SurfaceMesh& TeethWrapper::GetTeethSurfaceMesh()
{
    assert(m_teeth_sm && "Teeth surface mesh not set.");
    return *m_teeth_sm;
}

/**
 * @brief Sets the bite surface mesh by making a deep copy of the input surface mesh
 * @param sm The surface mesh to be set
 */
void TeethWrapper::SetBiteSurfaceMesh(SurfaceMesh& sm)
{
    //deep copy
    m_bite_sm = new SurfaceMesh(sm);
}

/**
 * @brief Converts the internal SurfaceMesh of teeth to vtkPolyData format with color.
 *
 * This method converts the CGAL SurfaceMesh representing the teeth into a VTK PolyData object,
 * which is useful for visualization. Each vertex of the mesh is added to the VTK PolyData as a point.
 * For each vertex that is part of the intersected vertex set, its distance from bite is used as a color value.
 * If a vertex is not part of the intersected vertex set, a color value of 0.0 is assigned.
 * The connectivity between vertices is also preserved in the conversion process.
 *
 * @return A smart pointer to the created vtkPolyData object.
 */
vtkSmartPointer<vtkPolyData>& TeethWrapper::SurfaceMeshToPolyDataWithColor()
{
    m_teeth_pd = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> VTKPoints = vtkSmartPointer<vtkPoints>::New();
    VTKPoints->SetNumberOfPoints((*m_teeth_sm).vertices().size());
    VLMap distance_property_map = GetDistanceFromBitePropertyMap();

    // Create a color array to store the colors
    m_distance_color_array = vtkSmartPointer<vtkFloatArray>::New();
    m_distance_color_array->SetName("Distances");

    for (uint32_t i = 0; i < (*m_teeth_sm).number_of_vertices(); i++)
    {
        vertex_descriptor vd = vertex_descriptor(i);
        const Point_3& point = (*m_teeth_sm).point(vd);
        VTKPoints->SetPoint(i, point.x(), point.y(), point.z());

        // If the point is intersected, set the color to the distance from bite
        if (m_intersected_v_set.find(vd) != m_intersected_v_set.end())
        {
            m_distance_color_array->InsertNextValue(distance_property_map[vd]);
        }
        else
        {
            m_distance_color_array->InsertNextValue(0.0);
        }
    }

    vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
    connectivity->SetNumberOfComponents(4);
    connectivity->SetNumberOfTuples((*m_teeth_sm).faces().size());

    for (int j = 0; j < (*m_teeth_sm).faces().size(); j++)
    {
        int ids[3] = { 0 };
        int index = 0;
        for (const halfedge_descriptor& h : (*m_teeth_sm).halfedges_around_face((*m_teeth_sm).halfedge(face_descriptor(j))))
        {
            ids[index] = (*m_teeth_sm).target(h).idx();
            index++;
        }
        if (index > 3)
            std::cerr << "SurfaceMeshToPolyData Conversion Error";

        connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
    }

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells((*m_teeth_sm).faces().size(), connectivity);

    m_teeth_pd->SetPoints(VTKPoints);
    m_teeth_pd->SetPolys(cells);
    m_teeth_pd->GetPointData()->SetScalars(m_distance_color_array);

    return m_teeth_pd;
}

/**
 * @brief Set the set of intersected vertices.
 *
 * @param intersected_v_set A set of vertex descriptors representing the intersected vertices.
 */
void TeethWrapper::SetIntersectedVertexSet(const std::set<vertex_descriptor>& intersected_v_set)
{
    m_intersected_v_set = intersected_v_set;
}

/**
 * @brief Get the set of intersected vertices.
 *
 * @return A reference to a set of vertex descriptors representing the intersected vertices.
 */
std::set<vertex_descriptor>& TeethWrapper::GetIntersectedVertexSet()
{
    return m_intersected_v_set;
}

/**
 * @brief Set the VTK actor for teeth visualization.
 *
 * @param actor A smart pointer to a VTK actor.
 */
void TeethWrapper::SetTeethActor(vtkSmartPointer<vtkActor>& actor)
{
    m_teeth_actor = actor;
}

/**
 * @brief Get the VTK actor for teeth visualization.
 *
 * @return A reference to a smart pointer of the VTK actor.
 */
vtkSmartPointer<vtkActor>& TeethWrapper::GetTeethActor()
{
    return m_teeth_actor;
}

/**
 * @brief Get all polylines for each layer.
 *
 * @return A reference to a vector of polylines.
 */
std::vector<Polylines>& TeethWrapper::GetAllLayerPolylines()
{
    return m_all_layer_polylines;
}

/**
 * @brief Set the renderer for the visualization.
 *
 * @param renderer A smart pointer to a VTK renderer.
 */
void TeethWrapper::SetRenderer(vtkSmartPointer<vtkRenderer> renderer)
{
    m_renderer = renderer;
}

/**
 * @brief Set the render window.
 *
 * @param render_win A smart pointer to a VTK render window.
 */
void TeethWrapper::SetRenderWindow(vtkSmartPointer<vtkRenderWindow> render_win)
{
    m_render_win = render_win;
}

/**
 * @brief Set the gap between slicing layers.
 *
 * @param gap The gap between layers, typically a double.
 */
void TeethWrapper::SetSlicingLayerGap(double gap)
{
    m_slicing_layer_gap = gap;
}

/**
 * @brief Get the vertices at the cusps.
 *
 * @return A reference to a vector of vertex descriptors representing the cusp vertices.
 */
std::vector<vertex_descriptor>& TeethWrapper::GetCuspVertices()
{
    return m_cusp_vertices;
}

/**
 * @brief Calculate vertex weights based on distance from a cusp vertex within a specified radius.
 *
 * This function performs a Breadth-First Search (BFS) starting from a cusp vertex and assigns
 * weights to vertices within a certain radius based on their distance from the cusp vertex.
 * The weight is determined using a sigmoid function of the distance, such that vertices closer
 * to the cusp vertex have higher weight. The modified faces that are adjacent to visited vertices
 * are also recorded in the process.
 *
 * @param cusp_vertex The vertex descriptor of the cusp vertex.
 * @param radius The radius within which to calculate vertex weights.
 * @param modified_faces A reference to an unordered_set that will be populated with the descriptors
 *                       of faces adjacent to the vertices visited during the BFS.
 * @return A map of vertex descriptors to their corresponding weights.
 */
std::unordered_map<vertex_descriptor, double> TeethWrapper::GetVerticesDistanceWeightsMap(
    vertex_descriptor cusp_vertex,
    double radius,
    std::unordered_set<face_descriptor>& modified_faces
)
{
    std::unordered_map<vertex_descriptor, double> weight_map;

    // Queue for Breadth-First Search (BFS)
    std::queue<vertex_descriptor> queue;
    double squared_radius = radius * radius; // Squared radius for distance comparison

    // Initialize the cusp vertex with the maximum weight
    weight_map[cusp_vertex] = 1.0;
    queue.push(cusp_vertex);

    while (!queue.empty())
    {
        vertex_descriptor u = queue.front();
        queue.pop();

        // Iterate over the halfedges around the current vertex
        for (halfedge_descriptor h : CGAL::halfedges_around_target(u, (*m_teeth_sm)))
        {
            // Skip border halfedges
            if (CGAL::is_border(h, (*m_teeth_sm)))
            {
                continue;
            }

            // Record the face of the halfedge as modified
            modified_faces.insert((*m_teeth_sm).face(h));

            // Get the opposite vertex of the halfedge
            vertex_descriptor v = source(h, (*m_teeth_sm));

            // Calculate the squared distance from the cusp vertex to the current vertex
            double distance = CGAL::squared_distance((*m_teeth_sm).point(cusp_vertex), (*m_teeth_sm).point(v));

            // If the distance is within the specified radius and the vertex hasn't been visited yet
            if (distance < squared_radius && weight_map.count(v) == 0)
            {
                queue.push(v);
                // Use a sigmoid function to calculate the weight based on the distance
                //weight_map[v] = 2.0 / (1.0 + exp(distance / radius)); // Sigmoid function for weight

                // Use a Gaussian function to calculate the weight based on the distance
                double sigma = 0.3 * radius;
                weight_map[v] = std::exp(-0.5 * std::pow(std::sqrt(distance) / sigma, 2)) / (sigma * std::sqrt(2 * M_PI));
            }
        }
    }

    return weight_map;
}

/**
 * @brief Get the maximum overlapping depth.
 *
 * This function retrieves the maximum depth of overlap among the teeth components.
 * Overlapping depth is a pre-calculated value that represents the maximum extent to which
 * one part of the teeth geometry might intersect with another part.
 *
 * If there are no overlapping parts, it returns std::numeric_limits<double>::lowest().
 *
 * @return The maximum overlapping depth as a double. If there is no overlap, it returns std::numeric_limits<double>::lowest().
 */
double TeethWrapper::GetMaxOverlappingDepth()
{
    return m_max_overlapping_depth;
}

std::shared_ptr<SurfaceMesh> TeethWrapper::GetAdjacentTeethSurfaceMesh()
{
    return m_adjacent_teeth_sm;
}

/**
 * @brief Adds a property map to the internal SurfaceMesh for storing vertex normals.
 *
 * This method creates a property map that associates each vertex of the internal SurfaceMesh
 * with a 3D vector representing its normal. The initial normal for each vertex is set to
 * CGAL::NULL_VECTOR (a zero vector). The property map is stored in the member variable
 * `m_teeth_vertices_normal`.
 *
 * @return A reference to the created property map.
 */
VNMap& TeethWrapper::AddVerticesNormalsPropertyMap()
{
    auto pair_map = (*m_teeth_sm).template add_property_map<vertex_descriptor, Vector_3>("v:teeth_vertices_normal", CGAL::NULL_VECTOR);
    assert(pair_map.second && "Property map 'v:teeth_vertices_normal' doesn't exists.");

    return m_vertices_normals_map = pair_map.first;
}

/**
 * @brief Retrieves the property map for the vertex normals from the internal SurfaceMesh.
 *
 * This method fetches the property map that associates each vertex of the internal SurfaceMesh
 * with a 3D vector representing its normal. The property map is stored in the member variable
 * `m_teeth_vertices_normal`.
 *
 * @return A reference to the fetched property map.
 *
 * @note If the property map does not exist, the SurfaceMesh's `property_map` method will
 * throw a `std::runtime_error`.
 */
VNMap& TeethWrapper::GetVerticesNormalsPropertyMap()
{
    auto pair_map = (*m_teeth_sm).property_map<vertex_descriptor, Vector_3>("v:teeth_vertices_normal");
    if (!pair_map.second)
    {
        std::cerr << "Property map 'v:teeth_vertices_normal' doesn't exists." << std::endl;

        return VNMap();
    }

    return m_vertices_normals_map = pair_map.first;
}

/**
 * @brief Adds a property map to the internal SurfaceMesh for storing distances from bite.
 *
 * This method creates a property map that associates each vertex of the internal SurfaceMesh
 * with a double value representing its distance from the bite. The initial distance for each vertex is set to
 * 0.0. The property map is stored in the member variable `m_distance_from_bite`.
 *
 * @return A reference to the created property map.
 */
VLMap& TeethWrapper::AddDistanceFromBitePropertyMap()
{
    auto pair_map = (*m_teeth_sm).template add_property_map<vertex_descriptor, double>("v:distance_from_bite", 0.0);
    assert(pair_map.second && "Property map 'v:distance_from_bite' doesn't exists.");

    return m_distance_from_bite_map = pair_map.first;
}

/**
 * @brief Retrieves the property map for the distances from bite from the internal SurfaceMesh.
 *
 * This method fetches the property map that associates each vertex of the internal SurfaceMesh
 * with a double value representing its distance from the bite. The property map is stored in
 * the member variable `m_distance_from_bite`.
 *
 * @return A reference to the fetched property map.
 *
 * @note If the property map does not exist, the SurfaceMesh's `property_map` method will
 * throw a `std::runtime_error`.
 */
VLMap& TeethWrapper::GetDistanceFromBitePropertyMap()
{
    auto pair_map = (*m_teeth_sm).property_map<vertex_descriptor, double>("v:distance_from_bite");
    assert(pair_map.second && "Property map 'v:distance_from_bite' doesn't exists.");

    return m_distance_from_bite_map = pair_map.first;
}

double TeethWrapper::GetMaxDimension()
{
    return m_max_dimension;
}

// If the vertex v overlaps with the bite, return the positive shift distance according to the projection direction
// Else return the negative shift distance
// If none of intersections is found, return the minimum double value
double TeethWrapper::GetProjectionDistanceFromBite(vertex_descriptor v)
{
    if (!m_bite_tree_ptr)
    {
        std::cerr << "m_bite_tree_ptr is not set" << std::endl;

        return std::numeric_limits<double>::lowest();
    }

    Point_3 origin_point = (*m_teeth_sm).point(v);
    Ray projection_ray(origin_point, m_projection_direction);
    Ray reversed_projection_ray(origin_point, -m_projection_direction);
    std::vector<Primitive_id> intersections;
    m_bite_tree_ptr->all_intersected_primitives(projection_ray, std::back_inserter(intersections));
    if (intersections.size() > 0)
    {
		auto intersection = m_bite_tree_ptr->first_intersection(projection_ray);
        if (intersection)
        {
             if (const Point_3* intersected_point = boost::get<Point_3>(&(intersection->first)))
             {
				return sqrt(CGAL::squared_distance(origin_point, *intersected_point));
			 }
		}
	}
	else
	{
		m_bite_tree_ptr->all_intersected_primitives(reversed_projection_ray, std::back_inserter(intersections));
        if (intersections.size() > 0)
        {
			auto intersection = m_bite_tree_ptr->first_intersection(reversed_projection_ray);
            if (intersection)
            {
                if (const Point_3* intersected_point = boost::get<Point_3>(&(intersection->first)))
                {
					return -sqrt(CGAL::squared_distance(origin_point, *intersected_point));
				}
			}
		}   
	}

    return std::numeric_limits<double>::lowest();
}

/**
 * @brief Renders the teeth model as a colored polydata.
 *
 * This method first checks if the max distance has been set and if a previous actor for the teeth exists,
 * removing the previous actor if necessary. Then it creates a lookup table for mapping distances to colors
 * and applies this color mapping to the polydata generated from the teeth's SurfaceMesh. Lastly, it creates
 * an actor from the colored polydata, adds it to the renderer and renders the scene.
 *
 * @note If the max distance is not set, an error message is printed and the function returns without doing anything.
 */
void TeethWrapper::RenderTeethPolydata()
{
    if (m_teeth_actor)
    {
        m_renderer->RemoveActor(m_teeth_actor);
    }

    // Create a color array to store the colors
    vtkSmartPointer<vtkFloatArray> distances_color_array = vtkSmartPointer<vtkFloatArray>::New();
    distances_color_array->SetName("Distances");

    vtkSmartPointer<vtkLookupTable> distance_lut = vtkSmartPointer<vtkLookupTable>::New();
    distance_lut->SetNumberOfColors(5); // Set the number of colors
    distance_lut->SetHueRange(1.0 / 3, 0.0); // Green to Red
    distance_lut->Build(); // Generate lookup table

    //distance_lut->SetTableValue(distance_lut->GetIndex(0), 1.0, 0.0, 1.0, 1.0);  
    distance_lut->SetTableValue(distance_lut->GetIndex(-1), 0 / 255.0, 238 / 255.0, 238 / 255.0, 1.0);  // Unintersected part

    m_teeth_pd = SurfaceMeshToPolyDataWithColor();

    // Map the colors using a vtkPolyDataMapper
    vtkSmartPointer<vtkPolyDataMapper> teeth_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    teeth_mapper->SetInputData(m_teeth_pd);
    teeth_mapper->SetScalarRange(-1, 5);
    teeth_mapper->SetLookupTable(distance_lut);
    teeth_mapper->SetInterpolateScalarsBeforeMapping(true);

    // Apply the mapper to the actor
    m_teeth_actor = vtkSmartPointer<vtkActor>::New();
    m_teeth_actor->SetMapper(teeth_mapper);
    m_teeth_actor->GetProperty()->SetOpacity(m_opacity);
    //m_teeth_actor->GetProperty()->EdgeVisibilityOn();
    m_renderer->AddActor(m_teeth_actor);
    m_render_win->Render();
}

/**
 * @brief Calculates the bite tree pointer by deep copying the input surface mesh and generating a tree from its faces
 * @param sm The surface mesh to be used for calculation
 */
void TeethWrapper::CalculateBiteTreePtr(SurfaceMesh sm)
{
    m_bite_sm = new SurfaceMesh(sm);
    m_bite_tree_ptr = std::make_shared<FaceTree>(CGAL::faces(*m_bite_sm).first, CGAL::faces(*m_bite_sm).second, *m_bite_sm);
}

/**
 * @brief Calculates the cut tree pointer by deep copying the input surface mesh and generating a tree from its faces
 * @param sm The surface mesh to be used for calculation
 */
void TeethWrapper::CalculateCutTreePtr(SurfaceMesh sm)
{
    m_cut_sm = new SurfaceMesh(sm);
    m_cut_tree_ptr = std::make_shared<FaceTree>(CGAL::faces(*m_cut_sm).first, CGAL::faces(*m_cut_sm).second, *m_cut_sm);
}

/**
 * @brief Calculates the slicer for the given SurfaceMesh.
 *
 * This function initializes a new SurfaceMesh based on the input, and then creates a HalfedgeTree
 * and a Polygon_mesh_slicer based on this new SurfaceMesh. These objects are stored for later use.
 *
 * @param sm The SurfaceMesh to calculate the slicer for.
 */
void TeethWrapper::CalculateTeethSlicer(SurfaceMesh sm)
{
    Timer timer("CalculateTeethSlicer");
    m_teeth_sm = new SurfaceMesh(sm);
	m_teeth_halfedge_tree_ptr = std::make_shared<HalfedgeTree>(edges(*m_teeth_sm).first, edges(*m_teeth_sm).second, *m_teeth_sm);
    m_teeth_slicer = std::make_shared<CGAL::Polygon_mesh_slicer<SurfaceMesh, Kernel>>(*m_teeth_sm, *m_teeth_halfedge_tree_ptr);
}

/**
 * @brief Calculates the polylines for all layers of the teeth.
 *
 * This function calculates a number of properties for the teeth SurfaceMesh, such as the maximum
 * and minimum x, y, and z coordinates and the maximum dimension. It then calculates the number of
 * layers and resizes the m_all_layer_polylines vector to accommodate these layers. Finally, it
 * calculates the polylines for each layer and adds them to m_all_layer_polylines.
 */
void TeethWrapper::CalculateAllLayerPolylines()
{
    Timer timer("CalculateAllLayerPolylines");
    double maxx = std::numeric_limits<double>::lowest();
    double minx = std::numeric_limits<double>::max();
    double maxy = std::numeric_limits<double>::lowest();
    double miny = std::numeric_limits<double>::max();
    double maxz = std::numeric_limits<double>::lowest();
    double minz = std::numeric_limits<double>::max();
    double symbol = 1.0;
    if (m_projection_direction * Vector_3(0, 1, 0) < 0)
    {
        symbol = -1.0;
    }
    for (auto v : (*m_teeth_sm).vertices())
    {
        Point_3 point = (*m_teeth_sm).point(v);
        maxx = std::max(maxx, point.x());
        minx = std::min(minx, point.x());
        maxz = std::max(maxz, point.z());
        minz = std::min(minz, point.z());
        double projection = m_projection_direction * (point - CGAL::ORIGIN) * symbol;
        maxy = std::max(maxy, projection);
        miny = std::min(miny, projection);
    }
    m_max_dimension = std::max(maxx - minx, maxz - minz);  // The maximum dimension of the teeth in x and z directions
    double starting_y = (m_projection_direction * Vector_3(0, 1, 0) > 0) ? miny : maxy;  // The starting y coordinate for slicing

    size_t number_of_layers = (std::abs(maxy - miny) / m_slicing_layer_gap) / 3;
    std::cout << "Maxy: " << maxy << " Miny: " << miny << std::endl;
    std::cout << "Slicing layer gap: " << m_slicing_layer_gap << std::endl;
    std::cout << "Number of layers: " << number_of_layers << std::endl;
    m_all_layer_polylines.clear();
    m_all_layer_polylines.resize(number_of_layers + 1);  // The iteration starts with (i - 1)

    Point_3 point(0, 0, 0);
    Vector_3 y_vector = Vector_3(0, starting_y, 0);
    for (int i = 0; i < number_of_layers; i++)
    {
        // The iteration starts with (i - 1) in order to avoid the failure to slice the first layer
        Point_3 current_point = point + y_vector + ((i - 1) * m_slicing_layer_gap) * m_projection_direction;
        Plane_3 plane(current_point, m_projection_direction);
        //std::cout << "Current point: " << current_point << " m_projection_direction: " << m_projection_direction << std::endl;
        (*m_teeth_slicer)(plane, std::back_inserter(m_all_layer_polylines[i]));

        //// Create a sphere
        //vtkSmartPointer<vtkSphereSource> sphere_source =
        //    vtkSmartPointer<vtkSphereSource>::New();
        //sphere_source->SetCenter(current_point[0], current_point[1], current_point[2]);
        //sphere_source->SetPhiResolution(20);
        //sphere_source->SetThetaResolution(20);
        //sphere_source->SetRadius(0.05);

        //// Create a sphere_mapper and actor
        //vtkSmartPointer<vtkPolyDataMapper> sphere_mapper =
        //    vtkSmartPointer<vtkPolyDataMapper>::New();
        //sphere_mapper->SetInputConnection(sphere_source->GetOutputPort());

        //vtkSmartPointer<vtkActor> sphere_actor =
        //    vtkSmartPointer<vtkActor>::New();
        //sphere_actor->SetMapper(sphere_mapper);
        //sphere_actor->GetProperty()->SetColor(0.0, 0.0, 1.0); // Set color to blue

        //// Add the actor to the scene
        //m_renderer->AddActor(sphere_actor);
    }
}

/**
 * @brief Calculates the cusp vertices for the teeth.
 *
 * This function iterates over the layers of polylines calculated in CalculateAllLayerPolylines and
 * identifies the cusp vertices. These vertices are then added to m_cusp_vertices. The function
 * stops searching once a quota of cusp points have been found.
 */
void TeethWrapper::CalculateCusp()
{
    Timer timer("CalculateCusp");
    m_cusp_vertices.clear(); 
    
    for (size_t k = 1; k < m_all_layer_polylines.size(); ++k)
    {  
        // If the enough quota of cusp points have been found, stop searching
        if (m_cusp_vertices.size() >= 5)
        {
            break;
        }

        auto polyline_prev = m_all_layer_polylines[k - 1];
        auto polyline_curr = m_all_layer_polylines[k];

        size_t size_prev = polyline_prev.size();
        size_t size_curr = polyline_curr.size();
        //std::cout << "Size of polyline_prev: " << size_prev << " " << "Size of polyline_curr: " << size_curr << std::endl;

        std::vector<std::vector<double>> polyline_curr_range(size_curr, std::vector<double>(4));

        std::vector<vtkNew<vtkPoints>> points_prev(size_prev), points_curr(size_curr);
        for (int i = 0; i < size_prev; i++)
        {
            for (int j = 0; j < polyline_prev[i].size(); j++)
            {
                points_prev[i]->InsertNextPoint(polyline_prev[i][j][0], polyline_prev[i][j][1], polyline_prev[i][j][2]);
            }
        }
        for (int i = 0; i < size_curr; i++)
        {
            double max_curr_x = -std::numeric_limits<double>::max();
            double max_curr_y = -std::numeric_limits<double>::max();
            double min_curr_x = std::numeric_limits<double>::max();
            double min_curr_y = std::numeric_limits<double>::max();

            for (int j = 0; j < polyline_curr[i].size(); j++)
            {
                double x = polyline_curr[i][j][0];
                double y = polyline_curr[i][j][1];
                points_curr[i]->InsertNextPoint(x, y, polyline_curr[i][j][2]);

                max_curr_x = std::max(max_curr_x, x);
                max_curr_y = std::max(max_curr_y, y);
                min_curr_x = std::min(min_curr_x, x);
                min_curr_y = std::min(min_curr_y, y);
            }
            std::vector<double> range = { min_curr_x, max_curr_x, min_curr_y, max_curr_y };
            polyline_curr_range[i] = range;
        }

        // New circle found
        if (size_curr - size_prev > 0)
        {
            std::vector<char> included(size_curr, 0);
            for (size_t i = 0; i < size_prev; ++i)
            {
                // New point found
                double point[3] = { 0.0 };
                points_prev[i]->GetPoint(0, point);

                double x = point[0];
                double y = point[1];

                // Check if the point is encircled by any polyline_curr
                for (size_t j = 0; j < size_curr; ++j)
                {
                    double min_curr_x = polyline_curr_range[j][0];
                    double max_curr_x = polyline_curr_range[j][1];
                    double min_curr_y = polyline_curr_range[j][2];
                    double max_curr_y = polyline_curr_range[j][3];

                    if (min_curr_x < x && x < max_curr_x && min_curr_y < y && y < max_curr_y)
                    {
                        included[j] = true;
                        break;
                    }
                }
            }

            for (size_t i = 0; i < size_curr; ++i)
            {
                if (!included[i])
                {
                    double point[3] = { 0.0 };
                    points_curr[i]->GetPoint(0, point);

                    // Check if the new point is too close to previous cusp points
                    Point_3 cusp_point = { point[0], point[1], point[2] };
                    //for (size_t j = 0; j < m_cusp_vertices.size(); ++j)
                    //{
                    //    Point_3 cusp_point_prev = m_cusp_vertices[j];
                    //    double3 cusp_point_prev = cusp_points[j];
                    //    double distance = (cusp_point - cusp_point_prev).getLength();
                    //    if (distance < 0.3)
                    //    {
                    //        return;
                    //    }
                    //}
                    for (auto prev_cusp_vertex : m_cusp_vertices)
                    {
                        Point_3 prev_cusp_point = m_teeth_sm->point(prev_cusp_vertex);
						double distance = std::sqrt(CGAL::squared_distance(cusp_point, prev_cusp_point));
                        if (distance < 1.5)
                        {
							return;
						}
					}
                    // The new point is not too close to previous cusp points
                    // , continue to locate the nearest vertex on the m_teeth_sm
                    std::pair<Point_3, HalfedgeTree::Primitive_id> closest = m_teeth_halfedge_tree_ptr->closest_point_and_primitive(cusp_point);
                    edge_descriptor closest_halfedge_id = closest.second;
                    vertex_descriptor cusp_vertex = target(closest_halfedge_id, *m_teeth_sm);
                    m_cusp_vertices.push_back(cusp_vertex);
                }
            }
        }
	}
}

/**
 * @brief Gets the cut tree
 * @return A reference to the cut tree
 * @throws An assertion error if the cut tree is not set
 */
FaceTree& TeethWrapper::GetCutTree()
{
    assert(m_cut_tree_ptr && "Cut tree not set.");
    return *m_cut_tree_ptr;
}

/**
 * @brief Calculates the intersected vertices and their distances from the bite, and the maximum overlapping depth.
 *
 * Please perform this method each time after adjusting teeth wall thickness and occlusion and before rendering.
 * 
 * This method first clears the set of intersected vertices. Then for each vertex in the teeth's SurfaceMesh,
 * it constructs a ray from the vertex in the direction of projection and tests for intersection with the bite.
 * If the vertex is found to intersect with the bite, it is added to the set of intersected vertices.
 *
 * The distance from each vertex to the bite is also calculated and stored in the `m_distance_from_bite` property map.
 * This distance is calculated as the minimum distance from the vertex to any point on the bite,
 * regardless of the direction of projection.
 *
 * In addition to this, the function now also calculates the maximum overlapping depth (`m_max_overlapping_depth`),
 * which is the maximum distance any vertex is found to be intersecting with the bite in the direction of `m_projection_direction`.
 *
 * @throws If `m_bite_tree_ptr` is not set, the function will print an error message to standard error and return immediately.
 */
void TeethWrapper::CalculateIntersectedVertexAndDistanceFromBite()
{
    if (!m_bite_tree_ptr)
    {
        std::cerr << "m_bite_tree_ptr is not set" << std::endl;
        
        return;
    }

    m_intersected_v_set.clear();
    m_distance_from_bite_map = GetDistanceFromBitePropertyMap();
    m_max_overlapping_depth = std::numeric_limits<double>::lowest();

    for (auto vd : (*m_teeth_sm).vertices())
    {
        Point_3 source_point = (*m_teeth_sm).point(vd);

        Ray ray(source_point, m_projection_direction);
        std::vector<Primitive_id> intersections;
        m_bite_tree_ptr->all_intersected_primitives(ray, std::back_inserter(intersections));
        if (intersections.size() % 2 == 0)
        {
            continue;
        }
        auto intersection = m_bite_tree_ptr->first_intersection(ray);

        if (intersection)
        {
            if (const Point_3* intersected_point = boost::get<Point_3>(&(intersection->first)))
            {
                m_intersected_v_set.insert(vd);
                double depth = sqrt(CGAL::squared_distance(source_point, *intersected_point));
                m_max_overlapping_depth = std::max(m_max_overlapping_depth, depth);
            }
        }

        // This vertex is not intersected with the bite or newly added
        // Unchanged vertices shall have unchanged values other than 0.0
        if (m_distance_from_bite_map[vd] == 0.0)
        {
            std::pair<Point_3, Primitive_id> pp = m_bite_tree_ptr->closest_point_and_primitive(source_point);
            m_distance_from_bite_map[vd] = sqrt(CGAL::squared_distance(source_point, pp.first));
        }
    }
}

void TeethWrapper::ExtractAdjacentTeethWithoutAbutment()
{
    //std::cout << "ExtractAdjacentTeethWithoutAbutment" << std::endl;
    if (!m_selected_id)
    {
        std::cerr << "m_selected_id is not set" << std::endl;
        return;
    }
    else if (!m_arch_pd)
    {
		std::cerr << "m_arch_pd is not set" << std::endl;
        return;
	}
    else if (!m_arch_sm)
    {
        std::cerr << "m_arch_sm is not set" << std::endl;
        return;
    }
    else if (m_arch_without_abutment_set.size() == 0)
    {
        std::cerr << "m_arch_without_abutment_set is not set" << std::endl;
        return;
    }

    ofstream arch_without_abutment_file("arch_without_abutment.xyz");
    for (auto& vd : m_arch_without_abutment_set)
	{
		Point_3 point = (*m_arch_sm).point(vd);
		arch_without_abutment_file << point.x() << " " << point.y() << " " << point.z() << std::endl;
	}
    arch_without_abutment_file.close();
    std::unordered_map<face_descriptor, face_descriptor> new_to_origin_face_map;
    unsigned expansion_level = 200;
    SurfaceMesh* expanded_abutment_sm = AreaExpander(*m_arch_sm, m_arch_pd, m_selected_id, new_to_origin_face_map, expansion_level);
    CGAL::IO::write_PLY("expanded_abutment_" + std::to_string(expansion_level) + ".ply", *expanded_abutment_sm);
    std::set<face_descriptor> expanded_abutment_face_set;
    for (auto& pair : new_to_origin_face_map)
	{
		expanded_abutment_face_set.insert(pair.second);
	}
    //std::cout << "expanded_abutment_face_set" << std::endl;
    std::set<vertex_descriptor> expanded_abutment_vertex_set;
    for (auto& fd : expanded_abutment_face_set)
	{
		for (auto he : halfedges_around_face(halfedge(fd, *m_arch_sm), *m_arch_sm))
		{
			expanded_abutment_vertex_set.insert(target(he, *m_arch_sm));
		}
	}
    std::ofstream expanded_abutment_file("expanded_abutment.xyz");
    std::ofstream expanded_abutment_vd_file("expanded_abutment_vd.txt");
    for (auto& vd : expanded_abutment_vertex_set)
    {
		Point_3 point = (*m_arch_sm).point(vd);
		expanded_abutment_file << point.x() << " " << point.y() << " " << point.z() << std::endl;
        expanded_abutment_vd_file << vd.idx() << std::endl;
    }
    expanded_abutment_file.close();
    expanded_abutment_vd_file.close();
    //std::cout << "expanded_abutment_vertex_set" << std::endl;

    //std::cout << "The number of vertices of m_arch_without_abutment_set: " << m_arch_without_abutment_set.size() << std::endl;
    //std::cout << "The number of vertices of expanded abutment: " << expanded_abutment_vertex_set.size() << std::endl;
    std::set<vertex_descriptor> adjacent_teeth_without_abutment_set; // The set of vertices of adjacent teeth without abutment
	std::set_intersection(m_arch_without_abutment_set.begin(), m_arch_without_abutment_set.end(), expanded_abutment_vertex_set.begin(), expanded_abutment_vertex_set.end(), std::inserter(adjacent_teeth_without_abutment_set, adjacent_teeth_without_abutment_set.begin()));
    //std::cout << "The number of vertices of adjacent_teeth_without_abutment_set: " << adjacent_teeth_without_abutment_set.size() << std::endl;
    std::ofstream adjacent_teeth_without_abutment_file("adjacent_teeth_without_abutment.xyz");
    for (auto& vd : adjacent_teeth_without_abutment_set)
	{
		Point_3 point = (*m_arch_sm).point(vd);
		adjacent_teeth_without_abutment_file << point.x() << " " << point.y() << " " << point.z() << std::endl;
	}
	adjacent_teeth_without_abutment_file.close();

    std::vector<vertex_descriptor> adjacent_teeth_without_abutment_vector(adjacent_teeth_without_abutment_set.begin(), adjacent_teeth_without_abutment_set.end());
    std::pair<SurfaceMesh, std::unordered_map<vertex_descriptor, vertex_descriptor>> 
        arch_without_abutment_pair = ExtractSubmeshAndVertexMapping(*m_arch_sm, adjacent_teeth_without_abutment_vector);
    //std::cout << arch_without_abutment_pair.first.number_of_vertices();
    m_adjacent_teeth_sm = std::make_shared<SurfaceMesh>(arch_without_abutment_pair.first);

    m_adjacent_teeth_face_tree_ptr = std::make_shared<FaceTree>(CGAL::faces(*m_adjacent_teeth_sm).first, CGAL::faces(*m_adjacent_teeth_sm).second, *m_adjacent_teeth_sm);
}

/**
 * @brief Sets the projection direction for the TeethWrapper.
 *
 * This method accepts a vector representing the projection direction, normalizes it, and then sets `m_projection_direction` with it.
 *
 * @param projection_direction The vector representing the desired projection direction. It does not need to be normalized.
 */
void TeethWrapper::SetProjectionDirection(const Vector_3 projection_direction)
{
    m_projection_direction = projection_direction;
    m_projection_direction = m_projection_direction / std::sqrt(m_projection_direction.squared_length());
    std::cout << "Projection_direciton: " << m_projection_direction << std::endl;
}

/**
 * @brief Sets the opacity for the TeethWrapper.
 *
 * This method accepts a double representing the desired opacity and sets `m_opacity` with it.
 *
 * @param opacity The desired opacity, typically a value between 0.0 (completely transparent) and 1.0 (completely opaque).
 */
void TeethWrapper::SetOpacity(double opacity)
{
    m_opacity = opacity;
}

/**
 * @brief Sets the maximum cusp offset threshold for the TeethWrapper.
 *
 * This method accepts a double representing the desired maximum cusp offset threshold, takes its absolute value, and then sets `m_max_cusp_offset_threshold` with it.
 *
 * @param threshold The desired maximum cusp offset threshold. The absolute value will be taken if a negative value is passed.
 */
void TeethWrapper::SetMaxCuspOffsetThreshold(double threshold)
{
    m_max_cusp_offset_threshold = threshold >= 0 ? threshold : -threshold;
}

void TeethWrapper::SetArchSurfaceMesh(SurfaceMesh& sm)
{
	m_arch_sm = &sm;
}

void TeethWrapper::SetArchPolyData(vtkSmartPointer<vtkPolyData> pd)
{
	m_arch_pd = pd;
}

void TeethWrapper::SetArchWithoutAbutmentSet(const std::set<vertex_descriptor>& arch_without_abutment_set)
{
    m_arch_without_abutment_set = arch_without_abutment_set;
}

void TeethWrapper::SetSelectedId(int id)
{
    m_selected_id = id;
}

/**
 * @brief Counts the number of vertices directly connected to a given vertex in a 3D surface mesh.
 *
 * @param v The vertex whose adjacent vertices are to be counted.
 * @return The number of vertices directly connected to the given vertex.
 */
int TeethWrapper::CountAdjacentVertices(vertex_descriptor v)
{
    // Initialize a counter for the adjacent vertices.
    int count = 0;

    // Get the range of vertices around the target vertex.
    auto range = (*m_teeth_sm).vertices_around_target((*m_teeth_sm).halfedge(v));

    // Iterate over the range of vertices.
    for (auto vcirc = range.begin(); vcirc != range.end(); ++vcirc)
    {
        // Increment the counter for each adjacent vertex.
        ++count;
    }

    // Return the count of adjacent vertices.
    return count;
}

/**
 * @brief Counts the number of vertices directly connected to a given vertex in a 3D surface mesh.
 *
 * @param sm The SurfaceMesh refers to.
 * @param v The vertex whose adjacent vertices are to be counted.
 * @return The number of vertices directly connected to the given vertex.
 */
int TeethWrapper::CountAdjacentVertices(SurfaceMesh& sm, vertex_descriptor v)
{
    // Initialize a counter for the adjacent vertices.
    int count = 0;

    // Get the range of vertices around the target vertex.
    auto range = sm.vertices_around_target(sm.halfedge(v));

    // Iterate over the range of vertices.
    for (auto vcirc = range.begin(); vcirc != range.end(); ++vcirc)
    {
        // Increment the counter for each adjacent vertex.
        ++count;
    }

    // Return the count of adjacent vertices.
    return count;
}

/**
 * @brief Smooths a 3D surface mesh using an iterative smoothing scheme.
 *
 * This function performs 15 iterations of a smoothing operation on the vertices of the given surface mesh.
 * In each iteration, it computes an offset for each vertex based on the average of the differences between
 * the vertex's position and the positions of its adjacent vertices. It then applies a fraction of this offset
 * to the vertex's position, unless the vertex is on the border of the mesh.
 *
 */
void TeethWrapper::SmoothMesh()
{
    // Number of vertices in the mesh.
    int vertex_number = (*m_teeth_sm).number_of_vertices();

    // Perform 15 iterations of the smoothing operation.
    for (int i = 0; i < 15; ++i)
    {
        std::vector<Vector_3> vDeltas(vertex_number, Vector_3(0, 0, 0));

        // Compute the offsets.
        int j = 0;
        for (auto& v : (*m_teeth_sm).vertices())
        {
            // Reference to the offset for the current vertex.
            Vector_3& delta = vDeltas[j];

            // Current vertex's position.
            Point_3 vf = (*m_teeth_sm).point(v);

            // Number of vertices adjacent to the current vertex.
            int vNerN = CountAdjacentVertices(v);

            // Iterate over the vertices adjacent to the current vertex.
            auto vcirc = (*m_teeth_sm).vertices_around_target((*m_teeth_sm).halfedge(v)).begin(), done((*m_teeth_sm).vertices_around_target((*m_teeth_sm).halfedge(v)).end());
            do
            {
                // Add the difference between the adjacent vertex's position and the current vertex's position 
                // to the offset, divided by the number of adjacent vertices.
                delta = delta + (((*m_teeth_sm).point(*vcirc) - vf) / vNerN);
            } while (++vcirc != done);

            ++j;
        }

        // Apply the offsets.
        j = 0;
        for (auto& v : (*m_teeth_sm).vertices())
        {
            // Reference to the offset for the current vertex.
            Vector_3& delta = vDeltas[j];

            // New position for the current vertex.
            Point_3 ans = (*m_teeth_sm).point(v) + delta * 0.05;

            // If the vertex is not on the border of the mesh, update its position.
            if (!(*m_teeth_sm).is_border(v))
            {
                (*m_teeth_sm).point(v) = ans;
            }

            ++j;
        }
    }
}


/**
 * @brief Smooths the mesh surface in the specified area.
 *
 * This function iterates over a given set of vertices, provided by a container (std::vector, std::set, std::unordered_set, ...), and applies a smoothing algorithm.
 * It adjusts the position of each vertex in the area based on the average of its adjacent
 * vertices. This process is repeated for the specified number of iteration times.
 * The smoothing intensity is lower for border vertices to prevent shrinkage.
 *
 * @tparam VertexContainer A container type that holds vertex descriptors. It must support iteration.
 * @param sm Reference to the SurfaceMesh object which contains the mesh to be smoothed.
 * @param area_vertices A container of vertex descriptors indicating the vertices in the mesh
 *                      that should be considered for smoothing.
 * @param radius The radius within which to expand the area of smoothing. Defaults to 1.0 if not specified.
 * @param iteration_times The number of times the smoothing operation is to be performed. Defaults to 15 if not specified.
 *
 * @note Border vertices are detected and handled differently to prevent the mesh from
 *       collapsing at the borders. The container type for area_vertices can be any iterable container
 *       such as std::vector, std::set, or std::unordered_set.
 *
 */
template <typename VertexContainer>
void TeethWrapper::SmoothMesh(
    SurfaceMesh& sm, 
    const VertexContainer& area_vertices,
    const double radius,
    const int iteration_times
)
{
    double squared_radius = radius * radius;

    std::set<vertex_descriptor> expanded_area_vertices;
    for (auto& vd : area_vertices)
    {
        Point_3 p0 = (*m_teeth_sm).point(vd);
        for (auto vd : (*m_teeth_sm).vertices())
        {
            Point_3 p1 = (*m_teeth_sm).point(vd);

            if (CGAL::squared_distance(p0, p1) < squared_radius)
            {
                expanded_area_vertices.insert(vd);
            }
        }
    }

    for (int k = 0; k < iteration_times; ++k)
    {
        std::vector<Vector_3> all_deltas(expanded_area_vertices.size(), Vector_3(0, 0, 0));
        int j = 0;
        for (auto& v : expanded_area_vertices)
        {
            Vector_3& delta = all_deltas[j];
            Point_3 p = sm.point(v);
            int num = CountAdjacentVertices(sm, v);
            auto circ = sm.halfedge(v);
            if (!sm.is_valid(circ)) 
            {
                continue;
            }
            if (sm.is_border(circ))
            {
                delta = delta + ((sm.point(sm.source(circ)) - p) / num);
                circ = sm.prev(sm.opposite(circ));
            }
            auto done = circ;
            do 
            {
                delta = delta + ((sm.point(sm.source(circ)) - p) / num);
                circ = sm.prev(sm.opposite(circ));
            } while (sm.is_border(circ) || circ != done);
            j++;
        }
        j = 0;
        for (auto& vd : expanded_area_vertices)
        {
            if (sm.is_border(vd))
            {
                Vector_3& delta = all_deltas[j];
                ++j;
                Point_3 ans = sm.point(vd);
                ans = ans + delta * 0.0005;
                sm.point(vd) = ans;
            }
            else
            {
                Vector_3& delta = all_deltas[j];
                ++j;
                Point_3 ans = sm.point(vd);
                ans = ans + delta * 0.05;
                sm.point(vd) = ans;
            }
        }
    }
}

/**
 * @brief Calculates the smooth normal of a vertex in a 3D surface mesh.
 *
 * The function iterates over all halfedges around the vertex, calculates the normal of the face they belong to,
 * and then averages all these normals to get the smooth normal. The normal is then normalized before being returned.
 *
 * @param v The vertex index in the mesh for which the smooth normal is being calculated.
 *
 * @return The normalized smooth normal as a 3D vector.
 */
Vector_3 TeethWrapper::CalculateSmoothNormal(vertex_descriptor v)
{
    Vector_3 normal(0, 0, 0);
    auto circ = (*m_teeth_sm).halfedge(v);

    // Keep a reference to the first halfedge to know when to stop iterating.
    auto done = circ;

    do {
        // Calculate the normal of the face that the current halfedge belongs to.
        Vector_3 vec = CGAL::normal(
            (*m_teeth_sm).point((*m_teeth_sm).source(circ)), // Vertex 1 of the face.
            (*m_teeth_sm).point((*m_teeth_sm).target(circ)), // Vertex 2 of the face.
            (*m_teeth_sm).point((*m_teeth_sm).target((*m_teeth_sm).next(circ))) // Vertex 3 of the face.
        );

        // Add the face's normal to the total.
        normal = normal + vec;

        // Move to the next halfedge around the vertex.
        circ = (*m_teeth_sm).opposite((*m_teeth_sm).next(circ));

        // Continue until all halfedges around the vertex have been handled.
    } while (circ != done);

    // Normalize the normal vector.
    normal = normal / std::sqrt(normal.squared_length());

    return normal;
}

/**
 * @brief Calculates the smooth normal for all vertices in a 3D surface mesh and stores them in a map.
 *
 * This function iterates over all vertices in the given surface mesh. For each vertex,
 * it calculates the smooth normal using the `calculate_smooth_normal` function and stores
 * the result in the provided map. The key in the map is the vertex index, and the value is the calculated smooth normal.
 *
 * @param map A map to store the calculated smooth normals for each vertex in the mesh.
 */
void TeethWrapper::SmoothAllNormal(VNMap& map)
{
    // Iterate over all vertices in the mesh.
    for (auto& v : (*m_teeth_sm).vertices())
    {
        // Calculate the smooth normal for the current vertex and store it in the map.
        map[v] = CalculateSmoothNormal(v);
    }
}

/**
 * @brief Expands a selected area within a mesh based on a label and creates a new mesh from the expanded area.
 *
 * This function works by first marking faces of the input mesh that correspond to a specific label
 * from a vtkPolyData object. It then iteratively expands the selection by adding adjacent faces.
 * After expansion, a new mesh is created containing only the selected faces. A mapping between
 * the faces of the new mesh and the faces of the original mesh is also provided.
 *
 * @param mesh The input mesh from which the area will be expanded.
 * @param pd A vtkSmartPointer to a vtkPolyData object that contains cell labels.
 * @param n The label indicating the initial area to expand from.
 * @param[out] face_map An unordered_map mapping faces in the new mesh to the corresponding faces in the input mesh.
 * @param[in] expansion_level The number of iterations for the expansion process.
 *
 * @return Returns a pointer to a new SurfaceMesh object that contains only the expanded selection of faces.
 */
SurfaceMesh* TeethWrapper::AreaExpander(SurfaceMesh& mesh, const vtkSmartPointer<vtkPolyData> pd, int n, std::unordered_map<face_descriptor, face_descriptor>& face_map, unsigned expansion_level = 50)
{
    auto selection = mesh.add_property_map<face_descriptor, bool>("f:select", false).first;
    vtkSmartPointer<vtkDataArray> labels = pd->GetCellData()->GetArray("Label");
    for (vtkIdType i = 0; i < labels->GetSize(); i++)
    {
        if (static_cast<int>(labels->GetTuple(i)[0]) == n)
        {
            selection[face_descriptor(i)] = true;
        }
    }
    std::vector<face_descriptor> expand_area;
    for (int i = 0; i < expansion_level; i++)
    {
        expand_area.clear();
        for (auto he : mesh.halfedges())
        {
            if (mesh.is_border(he) || mesh.is_border(mesh.opposite(he)))
            {
                continue;
            }
            if (selection[mesh.face(he)] == true && selection[mesh.face(mesh.opposite(he))] == false)
            {
                expand_area.push_back(mesh.face(mesh.opposite(he)));
            }
        }
        for (auto f : expand_area)
        {
            selection[f] = true;
        }
    }
    SurfaceMesh* sm_new = new SurfaceMesh;
    std::vector<int> visit(mesh.number_of_vertices(), -1);

    for (auto f : mesh.faces())
    {
        if (selection[f] == true)
        {
            auto v0 = mesh.source(mesh.halfedge(f));
            auto v1 = mesh.target(mesh.halfedge(f));
            auto v2 = mesh.target(mesh.next(mesh.halfedge(f)));
            vertex_descriptor v0_new, v1_new, v2_new;
            if (visit[v0.idx()] == -1)
            {
                v0_new = sm_new->add_vertex(mesh.point(v0));
                visit[v0.idx()] = v0_new.idx();
            }
            else
            {
                v0_new = vertex_descriptor(visit[v0.idx()]);
            }
            if (visit[v1.idx()] == -1)
            {
                v1_new = sm_new->add_vertex(mesh.point(v1));
                visit[v1.idx()] = v1_new.idx();
            }
            else
            {
                v1_new = vertex_descriptor(visit[v1.idx()]);
            }
            if (visit[v2.idx()] == -1)
            {
                v2_new = sm_new->add_vertex(mesh.point(v2));
                visit[v2.idx()] = v2_new.idx();
            }
            else
            {
                v2_new = vertex_descriptor(visit[v2.idx()]);
            }
            auto fd = sm_new->add_face(v0_new, v1_new, v2_new);
            face_map[fd] = f;
        }
    }

    return sm_new;
}

/**
 * Extracts a submesh from the original mesh using a list of vertices and creates a mapping
 * between vertex descriptors of the original mesh and the new submesh.
 *
 * @param original_mesh The original surface mesh from which to extract the submesh.
 * @param extracted_vertices A vector containing the descriptors of the vertices to include in the submesh.
 * @return A pair consisting of the new submesh and a map from the new mesh's vertex descriptors to the original mesh's vertex descriptors.
 */
std::pair<SurfaceMesh, std::unordered_map<vertex_descriptor, vertex_descriptor>> 
TeethWrapper::ExtractSubmeshAndVertexMapping(
    const SurfaceMesh& original_mesh,
    const std::vector<vertex_descriptor>& extracted_vertices
)
{
    SurfaceMesh new_mesh;
    std::unordered_map<vertex_descriptor, vertex_descriptor> original_to_new_vertex_map;
    std::unordered_map<vertex_descriptor, vertex_descriptor> new_to_original_vertex_map;

    // Add vertices to the new mesh and create the vertex mapping
    for (vertex_descriptor v : extracted_vertices)
    {
        vertex_descriptor new_v = new_mesh.add_vertex(original_mesh.point(v));
        original_to_new_vertex_map[v] = new_v;
        new_to_original_vertex_map[new_v] = v;
    }

    // Iterate over faces in the original mesh
    for (face_descriptor f : original_mesh.faces())
    {
        std::vector<vertex_descriptor> face_vertex_descriptors;
        for (vertex_descriptor v : vertices_around_face(original_mesh.halfedge(f), original_mesh))
        {
            // Check if the current vertex is in the extracted vertices map
            if (original_to_new_vertex_map.find(v) != original_to_new_vertex_map.end())
            {
                face_vertex_descriptors.push_back(original_to_new_vertex_map[v]);
            }
            else
            {
                // If any vertex of the face is not in the extracted list, discard the face
                face_vertex_descriptors.clear();
                break;
            }
        }
        // Add the face to the new mesh if all its vertices are in the extracted vertices
        if (!face_vertex_descriptors.empty())
        {
            new_mesh.add_face(face_vertex_descriptors);
        }
    }

    // Return the new mesh and the mapping from new to original vertex descriptors
    return { new_mesh, new_to_original_vertex_map };
}

/**
 * @brief Finds an adjacent face to a given vertex in a surface mesh.
 *
 * Given a vertex descriptor `vd`, this function searches for a non-border face adjacent to `vd` in the `mesh`.
 * If such a face is found, the descriptor of the first encountered adjacent face is returned. If the vertex is
 * associated with a border halfedge or there are no adjacent faces, a null face descriptor is returned.
 *
 * @param mesh The surface mesh to search within.
 * @param vd The vertex descriptor of the vertex whose adjacent face is to be found.
 * @return The descriptor of an adjacent face, or SurfaceMesh::null_face() if no adjacent face is found or if
 *         the vertex is on the border of the mesh.
 *
 * @note This function assumes that `mesh` is a valid SurfaceMesh object and `vd` is a valid vertex descriptor
 *       for this mesh. The function returns the first adjacent face found and does not guarantee that it is the
 *       "next" face in any particular traversal order. If the vertex is at the border of the mesh, no adjacent
 *       face will be found and the null face descriptor will be returned. The SurfaceMesh is assumed to adhere to
 *       the CGAL Surface Mesh data structure or similar, where `null_face()` and `null_halfedge()` are valid members
 *       indicating null or border conditions.
 */
face_descriptor TeethWrapper::FindAdjacentFace(const SurfaceMesh& mesh, const vertex_descriptor& vd)

{
    auto h = mesh.halfedge(vd);

    if (h == SurfaceMesh::null_halfedge() || mesh.is_border(h))
    {
        return SurfaceMesh::null_face();
    }

    auto start = h;
    auto face = mesh.face(h);
    if (face != SurfaceMesh::null_face())
    {
        return face;
    }
    return SurfaceMesh::null_face();
}

/**
 * @brief Adjusts the wall thickness of a tooth crown in a 3D surface mesh.
 *
 * @param thickness The desired thickness of the tooth crown wall.
 */
void TeethWrapper::AdjustCrownWallThickness(double thickness)
{
    if (!m_cut_tree_ptr)
	{
		std::cerr << "m_cut_tree_ptr is not set" << std::endl;

		return;
	}

    // Check if m_teeth_border_points is initialized.
    if (m_teeth_border_points.size() == 0)
    {
        for (auto vd : (*m_teeth_sm).vertices())
        {
            if (m_teeth_sm->is_border(vd))
            {
                Point_3 point = m_teeth_sm->point(vd);
                m_teeth_border_points.push_back(point);
            }
        }
    }

    VNMap vnormals_teeth = AddVerticesNormalsPropertyMap();
    PMP::compute_vertex_normals((*m_teeth_sm), vnormals_teeth);

    // Set to store the intersected face descriptors
    // std::set<face_descriptor> intersected_faces;
    std::map<vertex_descriptor, double> offset_map;
    std::vector<char> chosen_point((*m_teeth_sm).number_of_vertices(), 0);
    std::vector<vertex_descriptor> intersected_vertices;
    std::map<vertex_descriptor, bool> constraind_vertex_map;
    {
        // Timer intersection_timer("Intersection");
        for (auto v : (*m_teeth_sm).vertices())
        {
            Vector_3 n = get(vnormals_teeth, v);
            Point_3 origin_point = (*m_teeth_sm).point(v);
            Point_3 a = origin_point - thickness * n;
            Point_3 b = origin_point + 3.0 * n;

            Segment_3 segment_query(a, b);

            // computes first encountered intersection with segment query
            Segment_intersection intersection = m_cut_tree_ptr->any_intersection(segment_query);
            if (intersection)
            {
                //for (auto h : CGAL::halfedges_around_target(v, teeth))
                //{
                //	auto f = CGAL::face(h, teeth);

                //	if (f != boost::graph_traits<SurfaceMesh>::null_face() && !teeth.is_border(h))
                //	{
                //		intersected_faces.insert(f);
                //	}
                //}
                intersected_vertices.push_back(v);
                chosen_point[v.idx()] = true;
                constraind_vertex_map[v] = false;
                if (const Point_3* p = boost::get<Point_3>(&(intersection->first)))
                {
                    auto absolute_distance = sqrt(CGAL::squared_distance(origin_point, *p));
                    Vector_3 op = *p - origin_point;
                    //dotproduct
                    if (op * n > 0)
                    {
                        offset_map.insert(std::make_pair(v, absolute_distance));  // Teeth point locates inside cut
                    }
                    else
                    {
                        offset_map.insert(std::make_pair(v, -absolute_distance)); // Teeth point locates outside cut, need to check if the point is near the border
                        // Get distance from border
                        double min_distance = std::numeric_limits<double>::max();
                        for (auto border_point : m_teeth_border_points)
                        {
							double distance = sqrt(CGAL::squared_distance(origin_point, border_point));
                            if (distance < min_distance)
                            {
								min_distance = distance;
							}
						}
                        if (min_distance < 1.0)
                        {
							offset_map[v] = 0.0; // Do not offset
						}
                        else if (min_distance < 2.0)
                        {
                            offset_map[v] *= (min_distance - 1.0) * 2; // Gradually offset
                        }
                    }
                }
            }
            else
            {
                constraind_vertex_map[v] = true;
            }
        }
    }

    // To check if the offset_map is empty
    if (offset_map.size())
    {
        int iteration_times = 100;
        int smooth_interval = 20;

        for (int i = 0; i < iteration_times; i++)
        {
            for (auto v : (*m_teeth_sm).vertices())
            {
                if (chosen_point[v.idx()] && !(*m_teeth_sm).is_border(v))
                {
                    auto offset_vector = get(vnormals_teeth, v) * (offset_map[v] + thickness) / iteration_times;
                    (*m_teeth_sm).point(v) += offset_vector;
                }
            }
            //smooth_mesh(teeth);
            if (PMP::does_self_intersect((*m_teeth_sm)))
            {
                continue;
            }
            if (i % smooth_interval == 0)
            {
                SmoothMesh(*m_teeth_sm, intersected_vertices, 1.0, 15);
                //PMP::smooth_shape(intersected_faces, teeth, 1e-2/*, CGAL::parameters::vertex_is_constrained_map(constraint_vertex_map)*/);
            }
            //PMP::compute_vertex_normals(teeth, vnormals_teeth);
            SmoothAllNormal(vnormals_teeth);
        }

        SmoothMesh(*m_teeth_sm, intersected_vertices, 1.0, 15);
        //PMP::smooth_shape(intersected_faces, teeth, 1e-2/*, CGAL::parameters::vertex_is_constrained_map(constraint_vertex_map)*/);
    }

    CGAL::IO::write_polygon_mesh("teeth_after_thickness.ply", *m_teeth_sm);
}

//void TeethWrapper::AdjustCuspHeight()
//{
//    if (m_cusp_vertices.size() == 0)
//    {
//		std::cerr << "No cusp point is captured." << std::endl;
//		return;
//	}
//
//    if (!m_max_cusp_offset_threshold)
//    {
//        std::cerr << "m_max_cusp_offset_threshold is not set" << std::endl;
//        return;
//    }
//
//    double cusp_min_distance = std::numeric_limits<double>::max();
//    Surface_mesh_deformation deform_mesh(*m_teeth_sm);
//    std::unordered_set<face_descriptor> modified_faces;
//
//    // Find the minimum distance between cusp points, in order to avoid overlapping
//    for (size_t i = 0; i < m_cusp_vertices.size(); ++i)
//    {
//        for (size_t j = i + 1; j < m_cusp_vertices.size(); ++j)
//        {
//            Point_3 p1 = m_teeth_sm->point(m_cusp_vertices[i]);
//            Point_3 p2 = m_teeth_sm->point(m_cusp_vertices[j]);
//            double distance = std::sqrt(CGAL::squared_distance(p1, p2));
//            if (distance < cusp_min_distance)
//            {
//                cusp_min_distance = distance;
//            }
//        }
//    }
//
//    //for (vertex_descriptor v : m_cusp_vertices)
//    for (size_t i = 0; i < m_cusp_vertices.size(); ++i) // Only adjust the highest two cusp points
//    {
//        vertex_descriptor v = m_cusp_vertices[i];
//        // Set up the deformation for m_teeth_sm and add all points as ROI
//        deform_mesh.insert_roi_vertices(m_teeth_sm->vertices_begin(), m_teeth_sm->vertices_end());
//        double distance = GetProjectionDistanceFromBite(v);
//        std::cout << "Vertex " << i << " distance from bite: " << distance << std::endl;
//        // No intersection found
//        if (distance == std::numeric_limits<double>::lowest())
//        {
//            std::cout << "Cusp vertex " << i << " has no intersection with the bite." << std::endl;
//        }
//        else if (distance > m_max_cusp_offset_threshold || distance < -m_max_cusp_offset_threshold)
//        {
//			// If the offset distance exceeds the m_max_cusp_offset_threshold, do not offset
//            std::cout << "Cusp vertex " << i << " exceeds the threshold." << std::endl;
//		}
//		else
//        {
//            // If overlaps, the distance is negative, needs to offset cusp vertex along the projection direction
//            // Else the distance shall be positive, needs to offset cusp vertex along the reversed projection direction, the symbol is included in the distance
//            std::unordered_map<vertex_descriptor, double> total_vertices_with_weights;
//            total_vertices_with_weights = GetVerticesDistanceWeightsMap(v, cusp_min_distance / 2, modified_faces);
//            for (const auto& item : total_vertices_with_weights)
//            {
//                vertex_descriptor offset_v = item.first;
//                deform_mesh.insert_control_vertex(offset_v);
//                deform_mesh.set_target_position(offset_v, m_teeth_sm->point(offset_v) + m_projection_direction * distance * item.second);
//            }
//        }
//    }
//    deform_mesh.preprocess();
//    deform_mesh.deform(10000, 1e-4);
//    PMP::isotropic_remeshing(modified_faces, 0.05, *m_teeth_sm);
//}

/**
 * @brief Adjusts the height of the cusp points on the teeth mesh based on their distance from the bite.
 *
 * This method adjusts the height of the cusp vertices in the teeth's surface mesh. If no cusp vertices were found or
 * the maximum cusp offset threshold is not set, the function will print an error message and return immediately.
 * Otherwise, it calculates the minimum distance between cusp points to avoid overlapping, as well as the average distance
 * of all pairs of cusp points.
 *
 * For each cusp vertex, the function gets the projection distance from the bite. If the vertex does not intersect with
 * the bite or if the distance exceeds the maximum cusp offset threshold, the vertex is not offset. Otherwise, the vertex
 * is marked as modifiable and its offset is stored.
 *
 * Then, for each vertex in the teeth's surface mesh that is not a border vertex, the function calculates a target offset
 * based on the modifiable cusp vertices within a certain radius. If the target offset is not zero, the vertex is marked
 * as a control vertex, and its target position is set to its current position plus the target offset.
 *
 * Finally, the function preprocesses the deformation and applies it to the teeth's surface mesh. The mesh is then
 * isotropically remeshed for a better quality mesh.
 */
void TeethWrapper::AdjustCuspHeight()
{
    Timer timer("AdjustCuspHeight");
    if (m_cusp_vertices.size() == 0)
    {
        std::cerr << "No cusp point is captured." << std::endl;
        return;
    }

    if (!m_max_cusp_offset_threshold)
    {
        std::cerr << "m_max_cusp_offset_threshold is not set" << std::endl;
        return;
    }

    double cusp_min_distance = std::numeric_limits<double>::max();
    double cusp_distance_sum = 0.0;
    //Surface_mesh_deformation deform_mesh(*m_teeth_sm);
    //deform_mesh.insert_roi_vertices(m_teeth_sm->vertices_begin(), m_teeth_sm->vertices_end());

    std::unordered_set<face_descriptor> modified_faces;
    std::unordered_set<vertex_descriptor> modified_vertices;

    std::unordered_map<vertex_descriptor, double> modifiable_cusp_vertices_offset;

    // Find the minimum distance between cusp points, in order to avoid overlapping
    for (size_t i = 0; i < m_cusp_vertices.size(); ++i)
    {
        for (size_t j = i + 1; j < m_cusp_vertices.size(); ++j)
        {
            Point_3 p1 = m_teeth_sm->point(m_cusp_vertices[i]);
            Point_3 p2 = m_teeth_sm->point(m_cusp_vertices[j]);
            double distance = std::sqrt(CGAL::squared_distance(p1, p2));
            if (distance < cusp_min_distance)
            {
                cusp_min_distance = distance;
                cusp_distance_sum += distance;
            }
        }
    }
    //double radius = cusp_min_distance / 2;
    //double radius = cusp_distance_sum / (m_cusp_vertices.size() * (m_cusp_vertices.size() - 1) / 2);
    //double radius = cusp_distance_sum / m_cusp_vertices.size();
    double radius = m_max_dimension / 2;
    double squared_radius = radius * radius;
    //std::cout << radius << std::endl;
    for (vertex_descriptor v : m_cusp_vertices)
    {
        double distance = GetProjectionDistanceFromBite(v);
        std::cout << "Vertex " << v.idx() << " distance from bite: " << distance << std::endl;
        // No intersection found
        if (distance == std::numeric_limits<double>::lowest())
        {
            std::cout << "Cusp vertex " << v.idx() << " has no intersection with the bite." << std::endl;
        }
        else if (distance > m_max_cusp_offset_threshold || distance < -m_max_cusp_offset_threshold)
        {
            // If the offset distance exceeds the m_max_cusp_offset_threshold, do not offset
            std::cout << "Cusp vertex " << v.idx() << " exceeds the threshold." << std::endl;
        }
        else
        {
            modifiable_cusp_vertices_offset[v] = distance;
        }
    }

    if (modifiable_cusp_vertices_offset.size() == 0)
	{
		std::cout << "No cusp vertex is modifiable." << std::endl;
		return;
	}

    for (vertex_descriptor offset_v : m_teeth_sm->vertices())
    {
        if (CGAL::is_border(offset_v, *m_teeth_sm))
        {
            continue;
        }
        //std::cout << "Offset vertex " << offset_v.idx() << '\n';
        Vector_3 target_offset = Vector_3(0, 0, 0);
        double total_weight = 0.0;
        for (auto pair : modifiable_cusp_vertices_offset)
        {
            vertex_descriptor cusp_v = pair.first;
            double cusp_offset = pair.second;
            double distance = CGAL::squared_distance(m_teeth_sm->point(offset_v), m_teeth_sm->point(cusp_v));
            if (distance < squared_radius)
            {
                //auto h = m_teeth_sm->halfedge(offset_v);  // get a halfedge associated with the vertex
                //auto h_end = h;

                //do 
                //{
                //    // Get the face associated with the halfedge
                //    face_descriptor fd = m_teeth_sm->face(h);

                //    // Check if the face is valid (not a border halfedge)
                //    if (fd != m_teeth_sm->null_face()) 
                //    {
                //        modified_faces.insert(fd);
                //    }

                //    // Move to the next halfedge around the vertex
                //    h = m_teeth_sm->next(h);
                //} while (h != h_end);

                double sigma = 0.3 * radius;
                double weight = std::exp(-0.5 * std::pow(std::sqrt(distance) / sigma, 2)) / (sigma * std::sqrt(2 * M_PI));
                //double weight = 2.0 / (1.0 + std::exp(std::sqrt(distance) / radius));
                //double weight = 1.0 - sqrt(distance) / radius;


                target_offset += weight * cusp_offset * m_projection_direction;
                total_weight += weight;
            }
        }
        // If target_offset or total_weight equals to zero, do not offset
        if (target_offset == Vector_3(0, 0, 0))
        {
            continue;
        }
        //target_offset /= total_weight;
        //std::cout << "Target offset: " << target_offset << '\n';
        //deform_mesh.insert_control_vertex(offset_v);
        //deform_mesh.set_target_position(offset_v, m_teeth_sm->point(offset_v) + target_offset);
    }

    //deform_mesh.preprocess();
    //deform_mesh.deform(1000, 1e-3);
    (*m_teeth_sm).collect_garbage();
    //PMP::isotropic_remeshing(modified_faces, 0.05, *m_teeth_sm);
}

void TeethWrapper::OcclusionShaving(double offset)
{
    if (!m_bite_tree_ptr)
    {
        std::cerr << "m_bite_tree_ptr is not set" << std::endl;

        return;
    }

    m_intersected_v_set.clear();
    std::set<vertex_descriptor> m_changed_v_set;
    std::set<face_descriptor> extracted_faces;

    for (auto vd : m_teeth_sm->vertices())
    {
        Point_3 source_point = m_teeth_sm->point(vd);

        uint32_t index = vd.idx();

        Ray ray(source_point, m_projection_direction);
        std::vector<Primitive_id> intersections;
        m_bite_tree_ptr->all_intersected_primitives(ray, std::back_inserter(intersections));
        if (intersections.size() % 2 == 0)
        {
            continue;
        }
        auto intersection = m_bite_tree_ptr->first_intersection(ray);

        if (intersection)
        {
            if (const Point_3* intersected_point = boost::get<Point_3>(&(intersection->first)))
            {
                m_intersected_v_set.insert(vd);

                auto absolute_distance = sqrt(CGAL::squared_distance(source_point, *intersected_point));
                //m_max_intersection_distance = std::max(m_max_intersection_distance, absolute_distance);
                vertex_descriptor occluded_teeth_vd = vertex_descriptor(index);

                m_changed_v_set.insert(occluded_teeth_vd);
                m_teeth_sm->point(occluded_teeth_vd) = *intersected_point + m_projection_direction * offset;
            }
        }
    }

    for (auto f : m_teeth_sm->faces()) 
    {
        for (auto v : vertices_around_face(halfedge(f, *m_teeth_sm), *m_teeth_sm)) 
        {
            if (m_changed_v_set.find(v) != m_changed_v_set.end()) 
            {
                extracted_faces.insert(f);
                break; // Once we find a matching vertex, no need to check the rest for this face
            }
        }
    }

    PMP::isotropic_remeshing(
        extracted_faces,
        0.2,
        *m_teeth_sm,
        PMP::parameters::number_of_iterations(1)
    );
    m_teeth_sm->collect_garbage();
}

void TeethWrapper::ProximalShaving(double offset)
{
    if (!m_adjacent_teeth_face_tree_ptr)
    {
        std::cerr << "m_adjacent_teeth_sm is not set" << std::endl;

        return;
    }

    // Check if m_teeth_border_points is initialized.
    if (m_teeth_border_points.size() == 0)
    {
        for (auto vd : (*m_teeth_sm).vertices())
        {
            if (m_teeth_sm->is_border(vd))
            {
                Point_3 point = m_teeth_sm->point(vd);
                m_teeth_border_points.push_back(point);
            }
        }
    }


    int iteration_times = 10;
    int smooth_iteration = 10;
    int smooth_interval = 20;
    double protected_threshold = 1.0;
    double gradual_threshold = 2.0;
    std::set<vertex_descriptor> modified_vertices;


    VNMap teeth_vertices_normal = AddVerticesNormalsPropertyMap();
    PMP::compute_vertex_normals(*m_teeth_sm, teeth_vertices_normal);
        
    auto vertical_intersection_map = m_teeth_sm->add_property_map<vertex_descriptor, bool>("v:vertical_intersection", false).first;
    auto chosen_point_map = m_teeth_sm->add_property_map<vertex_descriptor, bool>("v:chosen_point", false).first;
    auto offset_map = m_teeth_sm->add_property_map<vertex_descriptor, double>("v:offset", 0.0).first;

    for (auto& v : (*m_teeth_sm).vertices())
    {
        Point_3 source_point = (*m_teeth_sm).point(v);
        Ray ray(source_point, -m_projection_direction);
        std::vector<Primitive_id> intersections;
        m_adjacent_teeth_face_tree_ptr->all_intersected_primitives(ray, std::back_inserter(intersections));
        if (intersections.size() % 2 == 0)
        {
            continue;
        }
        auto intersection = m_adjacent_teeth_face_tree_ptr->first_intersection(ray);
        if (intersection)
        {
            vertical_intersection_map[v] = true;
        }
    }

    for (int i = 0; i < iteration_times; ++i)
    {
        offset_map = m_teeth_sm->add_property_map<vertex_descriptor, double>("v:offset", 0.0).first;

        {
            //for (auto v : (*m_teeth_sm).vertices())
            for (auto v : m_teeth_sm->vertices())
            {
                if (vertical_intersection_map[v])
                {
                    Vector_3 n = get(teeth_vertices_normal, v);
                    Point_3 origin_point = (*m_teeth_sm).point(v);
                    Point_3 a = origin_point;
                    Point_3 b = origin_point - 5.0 * n;

                    /*vtkSmartPointer<vtkLineSource> line_source = vtkSmartPointer<vtkLineSource>::New();
                    line_source->SetPoint1(a.x(), a.y(), a.z());
                    line_source->SetPoint2(b.x(), b.y(), b.z());
                    line_source->Update();

                    vtkSmartPointer<vtkPolyDataMapper> line_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                    line_mapper->SetInputData(line_source->GetOutput());
                    line_mapper->Update();

                    vtkSmartPointer<vtkActor> line_actor = vtkSmartPointer<vtkActor>::New();
                    line_actor->SetMapper(line_mapper);
                    line_actor->GetProperty()->SetColor(1, 1, 0);
                    line_actor->GetProperty()->SetLineWidth(1);
                    m_renderer->AddActor(line_actor);*/

                    Segment_3 segment_query(a, b);

                    // computes first encountered intersection with segment query
                    Segment_intersection intersection = m_adjacent_teeth_face_tree_ptr->any_intersection(segment_query);
                    if (intersection)
                    {
                        chosen_point_map[v] = true;
                    }
                }

            }
        }

        PMP::compute_vertex_normals(*m_teeth_sm, teeth_vertices_normal);
        SmoothAllNormal(teeth_vertices_normal);
        for (auto v : (*m_teeth_sm).vertices())
        {
            if (chosen_point_map[v] && !(*m_teeth_sm).is_border(v))
            {
                Vector_3 n = get(teeth_vertices_normal, v);
                Point_3 origin_point = (*m_teeth_sm).point(v);
                Point_3 a = origin_point;
                Point_3 b = origin_point - 5.0 * n;
                Segment_3 segment_query(a, b);

                // computes first encountered intersection with segment query
                Segment_intersection intersection = m_adjacent_teeth_face_tree_ptr->any_intersection(segment_query);
                if (intersection)
                {
                    if (const Point_3* p = boost::get<Point_3>(&(intersection->first)))
                    {
                        auto absolute_distance = sqrt(CGAL::squared_distance(origin_point, *p));
                        Vector_3 op = *p - origin_point;
                        //dotproduct
                        if (op * n > 0)
                        {
                            offset_map[v] = absolute_distance;  // Teeth point locates inside arch
                        }
                        else
                        {
                            offset_map[v] = -absolute_distance; // Teeth point locates outside arch, need to check if the point is near the border
                            // Get distance from border
                            double min_distance = std::numeric_limits<double>::max();
                            for (auto& border_point : m_teeth_border_points)
                            {
                                double distance = sqrt(CGAL::squared_distance(origin_point, border_point));
                                if (distance < min_distance)
                                {
                                    min_distance = distance;
                                }
                            }
                            if (min_distance < protected_threshold)
                            {
                                //offset_map[v] = 0.0; // Do not offset
                                std::cout << min_distance << std::endl;
                                chosen_point_map[v] = false;
                                std::cout << "Skipped vertex: " << v.idx() << "\n";
                                continue;

                            }
                            else if (min_distance < gradual_threshold)
                            {
                                offset_map[v] *= (min_distance - protected_threshold) * protected_threshold / gradual_threshold; // Gradually offset
                            }
                            else
                            {
                            }
                        }
                    }
                }

                auto offset_vector = get(teeth_vertices_normal, v) * (offset_map[v]) / iteration_times;
                modified_vertices.insert(v);
                (*m_teeth_sm).point(v) += offset_vector;
            }
        }
    }
    {
        offset_map = m_teeth_sm->add_property_map<vertex_descriptor, double>("v:offset", 0.0).first;
        {
            //for (auto v : (*m_teeth_sm).vertices())
            for (auto v : m_teeth_sm->vertices())
            {
                if (vertical_intersection_map[v])
                {
                    Vector_3 n = get(teeth_vertices_normal, v);
                    Point_3 origin_point = (*m_teeth_sm).point(v);
                    Point_3 a = origin_point;
                    Point_3 b = origin_point - 5.0 * n;

                    /*vtkSmartPointer<vtkLineSource> line_source = vtkSmartPointer<vtkLineSource>::New();
                    line_source->SetPoint1(a.x(), a.y(), a.z());
                    line_source->SetPoint2(b.x(), b.y(), b.z());
                    line_source->Update();

                    vtkSmartPointer<vtkPolyDataMapper> line_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                    line_mapper->SetInputData(line_source->GetOutput());
                    line_mapper->Update();

                    vtkSmartPointer<vtkActor> line_actor = vtkSmartPointer<vtkActor>::New();
                    line_actor->SetMapper(line_mapper);
                    line_actor->GetProperty()->SetColor(1, 1, 0);
                    line_actor->GetProperty()->SetLineWidth(1);
                    m_renderer->AddActor(line_actor);*/

                    Segment_3 segment_query(a, b);

                    // computes first encountered intersection with segment query
                    Segment_intersection intersection = m_adjacent_teeth_face_tree_ptr->any_intersection(segment_query);
                    if (intersection)
                    {
                        chosen_point_map[v] = true;
                    }
                }

            }
        }

        PMP::compute_vertex_normals(*m_teeth_sm, teeth_vertices_normal);
        SmoothAllNormal(teeth_vertices_normal);
        for (auto v : (*m_teeth_sm).vertices())
        {
            if (chosen_point_map[v] && !(*m_teeth_sm).is_border(v))
            {
                Vector_3 n = get(teeth_vertices_normal, v);
                Point_3 origin_point = (*m_teeth_sm).point(v);
                Point_3 a = origin_point;
                Point_3 b = origin_point - 5.0 * n;
                Segment_3 segment_query(a, b);

                // computes first encountered intersection with segment query
                Segment_intersection intersection = m_adjacent_teeth_face_tree_ptr->any_intersection(segment_query);
                if (intersection)
                {

                    if (const Point_3* p = boost::get<Point_3>(&(intersection->first)))
                    {
                        auto absolute_distance = sqrt(CGAL::squared_distance(origin_point, *p));
                        Vector_3 op = *p - origin_point;
                        //dotproduct
                        if (op * n > 0)
                        {
                            offset_map[v] = absolute_distance;  // Teeth point locates inside arch
                        }
                        else
                        {
                            offset_map[v] = -absolute_distance; // Teeth point locates outside arch, need to check if the point is near the border
                            // Get distance from border
                            double min_distance = std::numeric_limits<double>::max();
                            for (auto& border_point : m_teeth_border_points)
                            {
                                double distance = sqrt(CGAL::squared_distance(origin_point, border_point));
                                if (distance < min_distance)
                                {
                                    min_distance = distance;
                                }
                            }
                            if (min_distance < protected_threshold)
                            {
                                //offset_map[v] = 0.0; // Do not offset
                                chosen_point_map[v] = false;
                                continue;
                            }
                            else if (min_distance < gradual_threshold)
                            {
                                offset_map[v] *= (min_distance - protected_threshold) * gradual_threshold / protected_threshold; // Gradually offset
                            }
                        }
                    }
                }

                auto offset_vector = get(teeth_vertices_normal, v) * (offset_map[v] - offset);
                modified_vertices.insert(v);
                (*m_teeth_sm).point(v) += offset_vector;
            }
        }
    }
    //PMP::compute_vertex_normals(*m_teeth_sm, teeth_vertices_normal);
    //for (auto v : (*m_teeth_sm).vertices())
    //{
    //    if (chosen_point_map[v] && !(*m_teeth_sm).is_border(v))
    //    {
    //        auto offset_vector = get(teeth_vertices_normal, v) * offset;
    //        (*m_teeth_sm).point(v) -= offset_vector;
    //    }
    //}
    std::set<face_descriptor> intersected_faces;
    for (auto& vd : modified_vertices) 
    {
        for (auto h : halfedges_around_target(vd, *m_teeth_sm))
        {
            face_descriptor fd = face(h, *m_teeth_sm);
            intersected_faces.insert(fd);
        }
    }
    std::cout << "Intersected faces: " << intersected_faces.size() << std::endl;
    // Check if the face contains border edge
    auto has_border_edge_or_not_valid = [](SurfaceMesh& sm, face_descriptor fd) -> bool 
    {
        if (!sm.is_valid(fd))
        {
			return true;
		}
        halfedge_descriptor hd = sm.halfedge(fd);
        for (int i = 0; i < 3; ++i, hd = sm.next(hd)) 
        {
            if (sm.is_border(sm.target(hd))) 
            {
                return true;
            }
        }
        return false;
    };

    // Find border faces around the modified faces.
    std::set<face_descriptor> modified_faces_border;
    std::set<vertex_descriptor> modified_vertices_border;
    {
        for (auto& fd : intersected_faces)
        {
            if (has_border_edge_or_not_valid(*m_teeth_sm, fd))
            {
                continue;
			}
            halfedge_descriptor hd = m_teeth_sm->halfedge(fd);
            for (int i = 0; i < 3; ++i, hd = m_teeth_sm->next(hd))
            {
                halfedge_descriptor hd_opposite = m_teeth_sm->opposite(hd);
                face_descriptor fd_opposite = m_teeth_sm->face(hd_opposite);

                if (
                    !has_border_edge_or_not_valid(*m_teeth_sm, fd_opposite)
                    && intersected_faces.find(fd_opposite) == intersected_faces.end() 
                )
                {
					modified_faces_border.insert(fd);
                    modified_faces_border.insert(fd_opposite);
				}
            }
        }
    }
    //std::cout << "Faces: " << m_teeth_sm->number_of_faces() << std::endl;
    //std::cout << "Vertices: " << m_teeth_sm->number_of_vertices() << std::endl;
    //std::cout << modified_faces_border.size() << std::endl;
    //for (auto& f : modified_faces_border)
    //{
    //    //std::cout << "Face: " << f.idx() << " ";
    //    halfedge_descriptor hd = m_teeth_sm->halfedge(f);
    //    vertex_descriptor v1 = m_teeth_sm->target(hd);
    //    vertex_descriptor v2 = m_teeth_sm->target(m_teeth_sm->next(hd));
    //    vertex_descriptor v3 = m_teeth_sm->target(m_teeth_sm->next(m_teeth_sm->next(hd)));
    //    Point_3 p1 = m_teeth_sm->point(v1);
    //    Point_3 p2 = m_teeth_sm->point(v2);
    //    Point_3 p3 = m_teeth_sm->point(v3);
    //    
    //    modified_vertices_border.insert(v1);
    //    modified_vertices_border.insert(v2);
    //    modified_vertices_border.insert(v3);

    //    //std::cout << v1.idx() << p1 << " " << v2.idx() << p2 << " " << v3.idx() << p3 << std::endl;

    //    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    //    points->InsertNextPoint(p1.x(), p1.y(), p1.z());
    //    points->InsertNextPoint(p2.x(), p2.y(), p2.z());
    //    points->InsertNextPoint(p3.x(), p3.y(), p3.z());

    //    vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
    //    line1->GetPointIds()->SetId(0, 0);
    //    line1->GetPointIds()->SetId(1, 1);

    //    vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
    //    line2->GetPointIds()->SetId(0, 1);
    //    line2->GetPointIds()->SetId(1, 2);

    //    vtkSmartPointer<vtkLine> line3 = vtkSmartPointer<vtkLine>::New();
    //    line3->GetPointIds()->SetId(0, 2);
    //    line3->GetPointIds()->SetId(1, 0);

    //    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    //    lines->InsertNextCell(line1);
    //    lines->InsertNextCell(line2);
    //    lines->InsertNextCell(line3);

    //    vtkSmartPointer<vtkPolyData> face_polydata = vtkSmartPointer<vtkPolyData>::New();
    //    face_polydata->SetPoints(points);
    //    face_polydata->SetLines(lines);
    //    vtkSmartPointer<vtkPolyDataMapper> face_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    //    face_mapper->SetInputData(face_polydata);
    //    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    //    actor->SetMapper(face_mapper);
    //    m_renderer->AddActor(actor);
    //}

    SmoothMesh(*m_teeth_sm, modified_vertices_border, 1.0, 20);
    PMP::isotropic_remeshing(intersected_faces, 0.1, *m_teeth_sm);

    m_teeth_sm->collect_garbage();

    CGAL::IO::write_polygon_mesh("teeth_after_proximal_shaving.ply", *m_teeth_sm);
}

/**
 * @brief Scales the teeth mesh towards the projection direction by a given offset.
 *
 * This method scales the teeth's surface mesh by a given offset. The offset is applied in several steps. For each step,
 * the function calculates a distance for each vertex to the closest border point, and then normalizes the distances so
 * that the maximum distance is 1. The vertices with distances less than 1/8 of the maximum distance are scaled by
 * their normalized distance, while the vertices with larger distances are not scaled.
 *
 * The vertices are then displaced along the projection direction by the offset multiplied by their scaling factor. The
 * faces that have at least one vertex with a scaling factor less than 1 are isotropically remeshed for a better quality
 * mesh. After each step, the surface mesh is cleaned up by removing unused vertices and faces.
 *
 * @param offset The distance to scale the teeth mesh towards the projection direction.
 */
void TeethWrapper::ScaleTeeth(double offset)
{
    Timer timer("ScaleTeeth");
    std::vector<Point_3> border;
    for (auto v : (*m_teeth_sm).vertices())
    {
        if ((*m_teeth_sm).is_border(v))
            border.push_back((*m_teeth_sm).point(v));
    }

    size_t steps = 5;
    double offset_per_step = offset / steps;

    for (size_t i = 0; i < steps; ++i)
    {
        std::vector<double> distance;
        for (auto v : (*m_teeth_sm).vertices())
        {
            double min = 100000;
            for (auto& b : border)
            {
                double d = CGAL::squared_distance((*m_teeth_sm).point(v), b);
                if (d < min) min = d;
            }
            distance.push_back(min);
        }
        double max = *std::max_element(distance.begin(), distance.end());

        for (auto& v : distance)
        {
            if (v * 8.0 < max)
            {
                v = v / max * 8.0;
            }
            else
            {
                v = 1.0;
            }
        }

        for (auto& v : (*m_teeth_sm).vertices())
        {
            auto& point = (*m_teeth_sm).point(v);
            auto displacement = distance[v.idx()] * offset_per_step;

            // Create a displacement vector along corrected_occlusal_direction
            Vector_3 displacement_vector = m_projection_direction * displacement;

            // Add displacement vector to the point
            point += displacement_vector;
            (*m_teeth_sm).point(v) = point;
        }

        std::vector<face_descriptor> area;
        for (auto f : (*m_teeth_sm).faces())
        {
            // Extract common halfedges and vertices
            auto halfedge0 = (*m_teeth_sm).halfedge(f);
            auto halfedge1 = (*m_teeth_sm).next(halfedge0);
            auto vertex0 = (*m_teeth_sm).source(halfedge0);
            auto vertex1 = (*m_teeth_sm).target(halfedge0);
            auto vertex2 = (*m_teeth_sm).target(halfedge1);

            // Check if any vertex is on the border
            if ((*m_teeth_sm).is_border(vertex0) || (*m_teeth_sm).is_border(vertex1) || (*m_teeth_sm).is_border(vertex2)) continue;

            // Compute the squared distance between the pick point and each point
            // If the distance is less than squared_radiusfor all points, add the face to the area
            if (distance[vertex0.idx()] < 1 || distance[vertex1.idx()] < 1 || distance[vertex2.idx()] < 1)
            {
                area.push_back(f);
            }
        }
        PMP::isotropic_remeshing(area, 0.05, *m_teeth_sm);
        (*m_teeth_sm).collect_garbage();
    }
}

/**
 * @brief Checks if the teeth mesh has self-intersections.
 *
 * This method checks if the teeth's surface mesh has any self-intersections. The check is performed by the
 * `does_self_intersect` function from the Polygon Mesh Processing (PMP) package of the Computational Geometry
 * Algorithms Library (CGAL).
 *
 * @return True if the teeth mesh has self-intersections, false otherwise.
 */
bool TeethWrapper::IsSelfIntersected() const
{
    return PMP::does_self_intersect((*m_teeth_sm));
}

void TeethWrapper::RemoveSelfIntersections()
{
    Face_intersections intersections;
    PMP::self_intersections(*m_teeth_sm, std::back_inserter(intersections));
    std::set<vertex_descriptor> intersected_vertices_set;
    for (auto& intersection : intersections)
    {
        face_descriptor fd = intersection.first;
        halfedge_descriptor hd = m_teeth_sm->halfedge(fd);
        vertex_descriptor vd0 = m_teeth_sm->target(hd);
        vertex_descriptor vd1 = m_teeth_sm->target(m_teeth_sm->next(hd));
        vertex_descriptor vd2 = m_teeth_sm->target(m_teeth_sm->next(m_teeth_sm->next(hd)));

        intersected_vertices_set.insert(vd0);
        intersected_vertices_set.insert(vd1);
        intersected_vertices_set.insert(vd2);
    }
    std::vector<vertex_descriptor> intersected_vertices_vec(intersected_vertices_set.begin(), intersected_vertices_set.end());
    SmoothMesh(*m_teeth_sm, intersected_vertices_vec, 2.0, 5);

    //SmoothMesh(*m_teeth_sm, intersected_vertices_vec);

    std::cout << "Attempted to remove self-intersections on number of " << intersected_vertices_vec.size() << " vertices." << std::endl;
}
