#include "CervicalMarginLineWrapper.h"

/**
 * @brief Default constructor for CervicalMarginLineWrapper class.
 *
 * Initializes a new instance of the CervicalMarginLineWrapper class with default values.
 * It sets the dental arch surface mesh pointer, the cervical margin line interactor style,
 * and the abutment edge spline pointer to nullptr. It also initializes the is_initialed flag
 * to false, indicating that the wrapper is not yet initialized, and sets the default projection
 * direction vector to (0, 0, 0).
 */
CervicalMarginLineWrapper::CervicalMarginLineWrapper() :
    m_arch_sm(nullptr),
    m_cervical_margin_line_interactor_style(nullptr),
    m_abutment_edge_spline(nullptr),
    m_is_initialed(false),
    m_projection_direction(Vector_3(0, 0, 0))
{}


/**
 * @brief Sets the projection direction and normalizes the vector.
 *
 * This function sets the projection direction of the CervicalMarginLineWrapper to the given vector.
 * It then normalizes the direction vector so that it has a unit length, which is essential for
 * many calculations that assume the direction vector is a unit vector.
 *
 * @param projection_direction The desired direction vector, which will be normalized.
 */
void CervicalMarginLineWrapper::SetProjectionDirection(const Vector_3 projection_direction)
{
    // Set the member variable to the input vector
    m_projection_direction = projection_direction;
    // Normalize the member variable by dividing it by its length
    m_projection_direction = m_projection_direction / std::sqrt(m_projection_direction.squared_length());
}

void CervicalMarginLineWrapper::SetArchPolyData(vtkSmartPointer<vtkPolyData> arch_pd)
{
	m_arch_pd = arch_pd;
}

void CervicalMarginLineWrapper::SetCervicalMarginLineInteractorStyle(ClosedSplineDesignInteractorStyle* cervical_margin_line_interactor_style)
{
   m_cervical_margin_line_interactor_style = cervical_margin_line_interactor_style;
}

void CervicalMarginLineWrapper::SetExpansionDistance(const double distance)
{
    m_expansion_distance = distance;
}

void CervicalMarginLineWrapper::SetSelectedId(const int selected_id)
{
    m_selected_id = selected_id;
}

void CervicalMarginLineWrapper::SetCtrlPtDensityCoefficient(const double coefficient)
{
	m_ctrl_pt_density_coefficient = coefficient;
}

void CervicalMarginLineWrapper::SetMinCurvatureThreshold(const double threshold)
{
	m_min_curvature_threshold = threshold;
}

void CervicalMarginLineWrapper::SetMeanCurvatureThreshold(const double threshold)
{
    m_mean_curvature_threshold = threshold;
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
CervicalMarginLineWrapper::ExtractSubmeshAndVertexMapping(
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
face_descriptor CervicalMarginLineWrapper::FindAdjacentFace(const SurfaceMesh& mesh, const vertex_descriptor& vd)
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
 * @brief Calculates and stores the border vertices of each disjoint border in a surface mesh.
 *
 * This function identifies all the disjoint borders of a given mesh and stores the vertices in
 * `border_vertices`. Each disjoint border is represented as a vector of vertex descriptors, and
 * all such vectors are stored in `border_vertices`. The function clears any existing data in
 * `border_vertices` before populating it.
 *
 * @param mesh The surface mesh from which to calculate borders.
 * @param border_vertices An output parameter where each element is a vector containing the vertex
 *        descriptors of one border in the mesh.
 *
 * @note This function will clear the contents of `border_vertices` before populating it. It uses
 *       an unordered set to keep track of visited halfedges to ensure that each border is only
 *       processed once. The function assumes that `mesh` is a valid SurfaceMesh object with
 *       functions such as `halfedges()`, `is_border()`, `next()`, and `target()` that comply
 *       with the expected interface for a CGAL Surface Mesh or similar mesh data structure.
 *       The `SurfaceMesh` is expected to have `null_halfedge()` as a member that indicates a
 *       null or non-existent halfedge. The vertices in each border will be ordered according to
 *       their appearance along the border's halfedges.
 */
void CervicalMarginLineWrapper::CalculateBorders(
    const SurfaceMesh& mesh,
    std::vector<std::vector<vertex_descriptor>>& border_vertices
)
{
    std::unordered_set<halfedge_descriptor> visited;
    border_vertices.clear();
    halfedge_descriptor result_halfedge = SurfaceMesh::null_halfedge();
    int num_borders = 0;
    for (halfedge_descriptor h : mesh.halfedges())
    {
        if (visited.find(h) == visited.end())
        {
            if (is_border(h, mesh))
            {
                std::vector<vertex_descriptor> border;
                auto endflag = h;
                auto he = mesh.next(h);
                while (endflag != he)
                {
                    visited.insert(he);
                    he = next(he, mesh);
                    border.push_back(target(he, mesh));
                }
                border_vertices.push_back(border);
                num_borders++;
            }
        }
    }
}

/**
 * Identifies the longest border of a mesh and maps its vertices to their original descriptors.
 *
 * Clears any existing data in `border_vertices` before populating it. It uses an internal function
 * to identify all borders within the `m_abutment_sm` mesh, selects the longest one, and then maps
 * the vertex descriptors from this border to their original descriptors as defined in the
 * `m_abutment_to_arch_vd_map`.
 *
 * The function will first clear the contents of `border_vertices`. It operates under the assumption
 * that `m_abutment_sm` is a member representing an abutment surface mesh and that
 * `m_abutment_to_arch_vd_map` is a member providing a mapping from abutment mesh vertex descriptors
 * to the original mesh vertex descriptors. The longest border is found and considered the edge border
 * of the mesh, and its vertex descriptors are translated to their original counterparts using the
 * mapping provided by `m_abutment_to_arch_vd_map`.
 *
 * The vertex descriptors in the longest border are stored in the order they appear along the border's
 * halfedges. It is assumed that the mapping in `m_abutment_to_arch_vd_map` is complete and can provide
 * a corresponding original vertex descriptor for every vertex descriptor in the abutment mesh.
 */
void CervicalMarginLineWrapper::CalculateEdgeBorder(const SurfaceMesh& mesh, std::vector<vertex_descriptor>& border_vertices)
{
    border_vertices.clear();
    // Extract the border vertices of all the holes
    std::vector<std::vector<vertex_descriptor>> borders;
    CalculateBorders(m_abutment_sm, borders);

    // Find the longest border i.e. the edge of the extracted mesh
    std::vector<vertex_descriptor> longest_border; // Storing vertices from the extracted mesh

    size_t max = 0;
    for (const auto& border : borders)
    {
        if (border.size() > max)
        {
            max = border.size();
            longest_border = border;
        }
    }
    for (const auto& vd : longest_border)
    {
        vertex_descriptor origin_vd = m_abutment_to_arch_vd_map.at(vd);
        border_vertices.push_back(origin_vd);
    }
}

bool CervicalMarginLineWrapper::GetBary(
    double v1[3], double v2[3], double v3[3],
    double vp[3],
    double& Bary1, double& Bary2, double& Bary3
)
{
    // Compute vectors
    double v02[3];
    v02[0] = v1[0] - v3[0];
    v02[1] = v1[1] - v3[1];
    v02[2] = v1[2] - v3[2];

    double v12[3];
    v12[0] = v2[0] - v3[0];
    v12[1] = v2[1] - v3[1];
    v12[2] = v2[2] - v3[2];

    double vp2[3];
    vp2[0] = vp[0] - v3[0];
    vp2[1] = vp[1] - v3[1];
    vp2[2] = vp[2] - v3[2];

    // Compute dot products
    double d0202 = v02[0] * v02[0] + v02[1] * v02[1] + v02[2] * v02[2];
    double d0212 = v02[0] * v12[0] + v02[1] * v12[1] + v02[2] * v12[2];
    double d1212 = v12[0] * v12[0] + v12[1] * v12[1] + v12[2] * v12[2];
    double d02p2 = v02[0] * vp2[0] + v02[1] * vp2[1] + v02[2] * vp2[2];
    double d12p2 = v12[0] * vp2[0] + v12[1] * vp2[1] + v12[2] * vp2[2];

    // Compute determinant
    double Det = d0202 * d1212 - d0212 * d0212;

    // Check if the determinant is non-zero
    if (Det == 0) {
        return false;
    }

    double dInvDet = static_cast<double>(1.0) / Det;

    // Compute barycentric coordinates
    Bary1 = (d1212 * d02p2 - d0212 * d12p2) * dInvDet;
    Bary2 = (d0202 * d12p2 - d0212 * d02p2) * dInvDet;
    Bary3 = static_cast<double>(1.0) - Bary1 - Bary2;

    // Check if the barycentric coordinates are valid
    return ((Bary1 >= 0) && (Bary2 >= 0) && (Bary1 + Bary2 < static_cast<double>(1.0)));
}

/**
 * Projects a set of 3D vertices from a surface mesh onto a 2D plane and creates a mapping
 * between the projected 2D points and their original vertex descriptors.
 *
 * @param mesh The surface mesh containing the 3D vertices to be projected.
 * @param vertices A vector of vertex descriptors whose associated points will be projected.
 * @param projected_points A reference to a vector that will be populated with the projected 2D points.
 * @return A map where each key is a Point_2 representing the projected point on the plane y = 0,
 *         and the corresponding value is the original vertex_descriptor from the mesh.
 *
 * This function iterates over each vertex descriptor provided in the `vertices` vector, retrieves
 * the corresponding 3D point from the `mesh`, and projects it onto the 2D plane by discarding
 * the y-coordinate. The resulting 2D point is then stored in `projected_points` and also used
 * as a key in the returned map, with its value being the original vertex descriptor. The function
 * clears any existing data in `projected_points` before adding new points.
 *
 * Note: If multiple 3D points project to the same 2D point, the map will only contain the last
 * vertex descriptor associated with that 2D point. Clients using this function should be aware of
 * this potential one-to-many mapping and handle it accordingly if necessary.
 *
 * Assumptions:
 * - The input `mesh` is a valid SurfaceMesh object and `vertices` contains valid vertex descriptors.
 * - There is a direct correspondence between each `vertex_descriptor` and a unique `Point_3` in the mesh.
 * - The projection plane is fixed at y = 0, parallel to the XZ-plane.
 * - The precision of the point projection is subject to the limits of floating-point arithmetic.
 */
std::map<Point_2, vertex_descriptor> CervicalMarginLineWrapper::ProjectToPlane(const SurfaceMesh& mesh, const std::vector<vertex_descriptor>& vertices, std::vector<Point_2>& projected_points)
{
    std::map<Point_2, vertex_descriptor> projection_map;
    projected_points.clear();

    for (const vertex_descriptor& vd : vertices) {
        const Point_3& p = mesh.point(vd);
        Point_2 proj(p.x(), p.z());  // Project onto the plane y = 0
        projection_map[proj] = vd;
        projected_points.push_back(proj);
    }

    return projection_map;
}

/**
 * Redistributes vertices along the surface mesh to achieve a more uniform spacing.
 *
 * @param mesh The surface mesh containing the vertices to be redistributed.
 * @param input_vertices A vector of vertex descriptors to be redistributed.
 * @return A vector of vertex descriptors representing the redistributed vertices.
 *
 * This function calculates the total distance around the path defined by `input_vertices`
 * and then uses a control point density coefficient (`m_ctrl_pt_density_coefficient`) to
 * determine a minimum interval distance. It iterates over each segment formed by consecutive
 * vertices and includes the starting vertex of the segment in the output only if the segment's
 * length is greater than or equal to the minimum interval.
 *
 * The goal is to create a set of output vertices that are more evenly spaced than the input,
 * in accordance with the total path length and the density coefficient. Vertices from the input
 * that are closer together than the minimum interval are skipped in the output.
 *
 * Assumptions:
 * - The input `mesh` is a valid `SurfaceMesh` object and `input_vertices` contains valid vertex descriptors.
 * - The member `m_ctrl_pt_density_coefficient` is a pre-set density coefficient that determines the desired
 *   spacing relative to the total path length.
 * - The path is implicitly closed, meaning the vertex following the last one in the input is the first one.
 * - The function does not insert new vertices; it only filters the input based on the spacing criterion.
 *
 * Note: This function does not necessarily result in vertices being exactly `min_interval` apart,
 * but no two consecutive vertices will be closer than `min_interval`.
 */
std::vector<vertex_descriptor> CervicalMarginLineWrapper::EquallyDistributeVertices(const SurfaceMesh& mesh, const std::vector<vertex_descriptor>& input_vertices)
{
    std::vector<vertex_descriptor> output_vertices;
    double total_distance = 0.0;

    for (int i = 0; i < input_vertices.size(); ++i)
    {
        Point_3 point_a = mesh.point(input_vertices[i]);
        Point_3 point_b = mesh.point(input_vertices[(i + 1) % input_vertices.size()]);
        double distance = std::sqrt(CGAL::squared_distance(point_a, point_b));

        total_distance += distance;
    }

    double min_interval = m_ctrl_pt_density_coefficient * total_distance / input_vertices.size();

    for (int i = 0; i < input_vertices.size(); ++i)
    {
        Point_3 point_a = mesh.point(input_vertices[i]);
        Point_3 point_b = mesh.point(input_vertices[(i + 1) % input_vertices.size()]);
        double distance = std::sqrt(CGAL::squared_distance(point_a, point_b));

        if (distance < min_interval)
        {
            //std::cout << "Skipping vertice " << i << " with distance " << distance << std::endl;
            continue;
        }
        output_vertices.push_back(input_vertices[i]);
    }

    return output_vertices;
}

void CervicalMarginLineWrapper::SetAbutmentPolyData(vtkSmartPointer<vtkPolyData> abutment_pd)
{
	m_abutment_pd = abutment_pd;
}

/**
 * Initializes the CervicalMarginLineWrapper class by setting up necessary data structures
 * and ensuring that all required components are properly set.
 *
 * Preconditions:
 * The following member variables must be set before calling Init:
 * - m_arch_sm: A pointer to the arch surface mesh.
 * - m_abutment_pd: A pointer to the abutment polydata.
 * - m_arch_pd: A pointer to the arch polydata.
 * - m_renderer: A pointer to the renderer.
 * - m_render_win: A pointer to the render window.
 * - m_projection_direction: A non-zero vector representing the projection direction.
 * - m_min_curvature_threshold: A double value set to the minimum curvature threshold.
 * - m_mean_curvature_threshold: A double value set to the mean curvature threshold.
 * - m_cervical_margin_line_interactor_style: A pointer to the interactor style.
 * - m_selected_id: A pointer to the selected ID.
 *
 * Overview of Operations:
 * - Verifies that all required member variables are set, otherwise prints an error message and returns.
 * - Converts the abutment polydata to a CGAL surface mesh.
 * - Maps vertices, faces, and halfedges to their corresponding indices.
 * - Checks if the halfedge property map "h:uv" exists, and if not, it is created.
 * - Parameterizes the arch surface mesh using CGAL's Surface Mesh Parameterization package.
 * - Initializes a closed mesh spline for the abutment edge.
 * - Sets the flag indicating that the initialization process is complete.
 *
 * Side Effects:
 * - If any of the checks fail, an error message will be printed to stderr and the function will return early.
 * - The member variable `m_is_initialed` is set to true if all operations complete successfully.
 *
 * Note:
 * - The function uses conditional compilation to include timing code if `ENABLE_TIMER_H` is defined.
 * - The function uses CGAL's surface mesh and surface mesh parameterization capabilities.
 */
void CervicalMarginLineWrapper::Init()
{
#ifdef ENABLE_TIMER_H
    Timer timer("CervicalMarginLineWrapper::Init");
#endif
    if (!m_arch_sm)
    {
		std::cerr << "Arch surface mesh is not set." << std::endl;
		return;
	}
    if (!m_abutment_pd)
    {
        std::cerr << "Abutment polydata is not set." << std::endl;
        return;
    }
    if (!m_arch_pd)
    {
		std::cerr << "Arch polydata is not set." << std::endl;
		return;
	}
    if (!m_renderer)
    {
		std::cerr << "Renderer is not set." << std::endl;
        return;
	}
    if (!m_render_win)
    {
		std::cerr << "Render window is not set." << std::endl;
        return;
	}
    if (m_projection_direction == Vector_3(0, 0, 0))
    {
        std::cerr << "Projection direction is not set." << std::endl;
        return;
    }
    if (m_min_curvature_threshold == std::numeric_limits<double>::max() || m_mean_curvature_threshold == std::numeric_limits<double>::max())
    {
        std::cerr << "Curvature thresholds are not set." << std::endl;
        return;
    }
    if (!m_cervical_margin_line_interactor_style)
    {
		std::cerr << "Cervical margin line interactor style is not set." << std::endl;
		return;
	}
    if (!m_selected_id)
    {
        std::cerr << "Selected id is not set." << std::endl;
        return;
    }

	// Convert the abutment polydata to a CGAL surface mesh
    PolyDataToSurfaceMesh(m_abutment_pd, m_abutment_sm);

    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(*m_arch_sm).first;
    int ver_index = 0;
    for (auto v : m_arch_sm->vertices()) m_vmap[ver_index++] = v;
    int face_index = 0;
    for (auto f : m_arch_sm->faces()) m_fmap[face_index++] = f;
    int halfedge_index = 0;
    for (auto he : m_arch_sm->halfedges()) m_hemap[halfedge_index++] = he;
    // Check if the "h:uv" property map exists
    bool found;
    boost::tie(m_uv_map, found) = m_arch_sm->property_map<vertex_descriptor, Point_2>("h:uv");
    // If it doesn't exist, add the property map
    if (!found) 
    {
        m_uv_map = m_arch_sm->add_property_map<vertex_descriptor, Point_2>("h:uv").first;
    }
    CGAL::Surface_mesh_parameterization::parameterize(*m_arch_sm, bhd, m_uv_map);

    m_abutment_edge_spline = new ClosedMeshSpline;
    m_abutment_edge_spline->initial(m_arch_sm, m_fmap, m_vmap, m_emap, m_hemap, m_uv_map);

    m_is_initialed = true;
}

/**
 * @brief Extracts and processes the surface mesh of the abutment.
 *
 * This function performs several operations to extract the surface mesh of the abutment
 * based on curvature analysis. It first checks if the class is properly initialized.
 * If not, it exits with an error message. The function expands the arch surface mesh,
 * converts it to Eigen matrices, and finds a starting vertex for a breadth-first search (BFS).
 * It then calculates mean and minimum curvatures, and based on these values, extracts vertices
 * that meet the curvature threshold criteria. It concludes by extracting a submesh and vertex mapping.
 *
 * @pre m_is_initialed must be true, indicating that the class is properly initialized.
 * @post On success, m_abutment_sm is populated with the extracted submesh, and m_abutment_to_arch_vd_map
 *       contains the mapping of vertices from the extracted to the original mesh.
 *
 * @note The function uses igl library functions to compute curvature and extract the mesh based on it.
 *       It also uses a BFS for the extraction process and employs conditional compilation for timing code.
 *
 * @warning If the class is not initialized, the function will output an error message to std::cerr
 *          and no further processing will take place.
 */
void CervicalMarginLineWrapper::ExtractAbutmentSurfaceMesh()
{
    if (!m_is_initialed)
    {
		std::cerr << "The CervicalMarginLineWrapper is not initialed." << std::endl;
		return;
	}

    using namespace Eigen;

    m_expanded_abutment_sm = *AreaExpander(*m_arch_sm, m_arch_pd, m_selected_id, m_expanded_face_map, 100);

    CGAL::IO::write_PLY("expanded_abutment.ply", m_expanded_abutment_sm);
    // Convert the arch surface mesh to Eigen matrices
    CGALSurfaceMeshToEigen(m_expanded_abutment_sm, m_V, m_F);

    // Find the start vertex for BFS on abutment surface mesh
    double abutment_center[3];
    m_abutment_pd->GetCenter(abutment_center);
    Point_3 abutment_center_point(abutment_center[0], abutment_center[1], abutment_center[2]);
    Vector_3 ray_direction = m_projection_direction;
    Ray_3 ray(abutment_center_point, ray_direction);
    Tree tree(faces(m_expanded_abutment_sm).first, faces(m_expanded_abutment_sm).second, m_expanded_abutment_sm);
    auto intersection = tree.first_intersection(ray);
    if (intersection)
    {
        //std::cout << "Intersection with mesh at " << intersection->first << std::endl;
        face_descriptor fd = intersection->second;
        CGAL::Halfedge_around_face_circulator<SurfaceMesh> circ(m_expanded_abutment_sm.halfedge(fd), m_expanded_abutment_sm), done(circ);
        m_bfs_start_vd = source(*circ, m_expanded_abutment_sm);
    }
    else
    {
        ray_direction = -ray_direction;
        //std::cout << "Changing ray direction to " << ray_direction << std::endl;
        Ray_3 reverse_ray(abutment_center_point, ray_direction);
        auto reverse_intersection = tree.first_intersection(reverse_ray);
        if (reverse_intersection)
        {
            //std::cout << "Intersection with mesh at " << reverse_intersection->first << std::endl;
            face_descriptor fd = reverse_intersection->second;
            CGAL::Halfedge_around_face_circulator<SurfaceMesh> circ(m_expanded_abutment_sm.halfedge(fd), m_expanded_abutment_sm), done(circ);
            m_bfs_start_vd = source(*circ, m_expanded_abutment_sm);
        }
        else
        {
            std::cout << "No intersection with mesh" << std::endl;
            return;
        }
    }
    Point_3 bfs_start_point = m_expanded_abutment_sm.point(m_bfs_start_vd);

    auto mean_curvature = m_expanded_abutment_sm.add_property_map<vertex_descriptor, double>("v:mean_curvature", 0.0).first;
    auto min_curvature = m_expanded_abutment_sm.add_property_map<vertex_descriptor, double>("v:min_curvature", 0.0).first;

    MatrixXd HN;
    SparseMatrix<double> L, M, Minv;
    igl::cotmatrix(m_V, m_F, L);
    igl::massmatrix(m_V, m_F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, Minv);

    // Laplace-Beltrami of position
    HN = -Minv * (L * m_V);
    // Extract magnitude as mean curvature
    VectorXd H = HN.rowwise().norm();

    VectorXd H5 = H;
    VectorXd H10 = H;
    VectorXd H15 = H;
    VectorXd H20 = H;

    // Compute curvature directions via quadric fitting
    MatrixXd PD1, PD2;
    MatrixXd PD1_5, PD2_5;
    MatrixXd PD1_10, PD2_10;
    MatrixXd PD1_15, PD2_15;
    MatrixXd PD1_20, PD2_20;
    VectorXd PV1, PV2;
    VectorXd PV1_5, PV2_5;
    VectorXd PV1_10, PV2_10;
    VectorXd PV1_15, PV2_15;
    VectorXd PV1_20, PV2_20;
    {
#ifdef ENABLE_TIMER_H
        Timer timer("Calculate principal_curvature");
        {
            Timer timer("Calculate principal_curvature 5");
#endif
            igl::principal_curvature(m_V, m_F, PD1_5, PD2_5, PV1_5, PV2_5, 5U, true);
#ifdef ENABLE_TIMER_H
        }
        {
            Timer timer("Calculate principal_curvature 10");
#endif
            igl::principal_curvature(m_V, m_F, PD1_10, PD2_10, PV1_10, PV2_10, 10U, true);
#ifdef ENABLE_TIMER_H
        }
        {
            Timer timer("Calculate principal_curvature 15");
#endif
            igl::principal_curvature(m_V, m_F, PD1_15, PD2_15, PV1_15, PV2_15, 15U, true);
#ifdef ENABLE_TIMER_H
        }
#endif
    }

    // mean curvature
    H5 = 0.5 * (PV1_5 + PV2_5);
    H10 = 0.5 * (PV1_10 + PV2_10);
    H15 = 0.5 * (PV1_15 + PV2_15);

    // Smooth mean curvature
    for (vertex_descriptor& vd : m_expanded_abutment_sm.vertices())
    {
        double h5_sum = 0.0;
        double h10_sum = 0.0;
        double h15_sum = 0.0;
        int h5_count = 0;
        int h10_count = 0;
        int h15_count = 0;

        // Traverse adjacent vertices
        auto circ = m_expanded_abutment_sm.vertices_around_target(m_expanded_abutment_sm.halfedge(vd));
        for (auto vic = circ.begin(); vic != circ.end(); ++vic) 
        {
            h5_sum += H5[vic->idx()];
            h10_sum += H10[vic->idx()];
            h15_sum += H15[vic->idx()];

            h5_count++;
            h10_count++;
            h15_count++;
        }
        H5[vd.idx()] = h5_sum / h5_count;
        H10[vd.idx()] = h10_sum / h10_count;
        H15[vd.idx()] = h15_sum / h15_count;
    }

    for (vertex_descriptor& vd : m_expanded_abutment_sm.vertices())
    {
        min_curvature[vd] = std::min(std::min(H5[vd.idx()], H10[vd.idx()]), H15[vd.idx()]);
        mean_curvature[vd] = H15[vd.idx()];
    }

    std::vector<vertex_descriptor> extracted_vertices;
    std::set<vertex_descriptor> difference_vertices;

    // Extract vertices via either mean or min curvature
    bool is_successful = BFSMeanCurvatureExtraction(m_expanded_abutment_sm, m_bfs_start_vd, extracted_vertices, difference_vertices, m_mean_curvature_threshold);
    std::cout << "Extracted vertices: " << extracted_vertices.size() << std::endl;
    std::cout << "Teeth polydata vertices: " << m_abutment_pd->GetNumberOfPoints() << std::endl;
    if (!is_successful)
    {
        std::cout << "Min curvature extraction method is used!" << std::endl;
        BFSMinCurvatureExtraction(m_expanded_abutment_sm, m_bfs_start_vd, extracted_vertices, difference_vertices, m_min_curvature_threshold);
    }
    else
    {
        std::cout << "Mean curvature extraction method is used!" << std::endl;
    }

    std::ofstream extracted_vertices_ofs("extracted_vertices.xyz");
    for (const auto& vd : extracted_vertices)
	{
		extracted_vertices_ofs << m_expanded_abutment_sm.point(vd) << "\n";
	}
    extracted_vertices_ofs.close();

    // Extract the submesh and the vertex mapping of the extracted to the original mesh
    std::pair<SurfaceMesh, std::unordered_map<vertex_descriptor, vertex_descriptor>> extracted_mesh_pair = ExtractSubmeshAndVertexMapping(m_expanded_abutment_sm, extracted_vertices);
    m_abutment_sm = extracted_mesh_pair.first;
    m_abutment_to_arch_vd_map = extracted_mesh_pair.second;
    CGAL::IO::write_PLY("extracted_abutment.ply", m_abutment_sm);
}

/**
 * @brief Generates a spline along the edge of an abutment.
 *
 * This method computes the edge of the abutment's surface mesh and creates a spline
 * that outlines the margin line. It is primarily used in dental CAD/CAM applications
 * for modeling dental restorations.
 *
 * The process begins by calculating the convex hull of the projected abutment border.
 * Vertices are then equally distributed along the convex hull to prevent a zigzag pattern
 * in the resulting spline. A spline is generated through these points by adding control points
 * that are slightly inset from the actual vertex positions to ensure the spline runs along the edge.
 *
 * @note If ENABLE_TIMER_H is defined, this method will also time and report the duration
 *       of the spline generation process.
 * @note This method contains commented code for VTK-based visualization, which can be
 *       used for debugging purposes.
 *
 * Usage:
 * @code
 * CervicalMarginLineWrapper wrapper;
 * wrapper.GenerateAbutmentEdgeSpline();
 * @endcode
 *
 * @pre m_abutment_sm must be a valid surface mesh.
 * @post m_abutment_edge_spline will contain the generated spline.
 */
void CervicalMarginLineWrapper::GenerateAbutmentEdgeSpline()
{
#ifdef ENABLE_TIMER_H
    Timer timer("Generate abutment edge spline");
#endif

    CalculateEdgeBorder(m_abutment_sm, m_abutment_border);
    std::cout << "Abutment border vertices: " << m_abutment_border.size() << std::endl;
    std::ofstream border_vertices_ofs("border_vertices.xyz");
    for (const auto& vd : m_abutment_border)
    {
        // m_abutment_to_arch_vd_map may not contain the vertex descriptor, so we need to check and drop it
        if (m_abutment_to_arch_vd_map.find(vd) == m_abutment_to_arch_vd_map.end())
        {
            m_abutment_to_arch_vd_map.erase(vd);
            continue;
        }
        vertex_descriptor origin_vd = m_abutment_to_arch_vd_map.at(vd);
        border_vertices_ofs << m_expanded_abutment_sm.point(origin_vd) << "\n";
    }
    border_vertices_ofs.close();
    std::cout << "Abutment border vertices(cleansed): " << m_abutment_border.size() << std::endl;

    std::vector<Point_2> projected_points;
    std::vector<Point_2> convex_hull_points;
    std::vector<vertex_descriptor> convex_hull_vertices;

    // Project the border vertices to the plane and find the convex hull
    std::map<Point_2, vertex_descriptor> point2_vd_map = ProjectToPlane(m_expanded_abutment_sm, m_abutment_border, projected_points);
    CGAL::convex_hull_2(projected_points.begin(), projected_points.end(), std::back_inserter(convex_hull_points));

    // Map the convex hull points to the original mesh
    for (const auto& point : convex_hull_points)
    {
        if (point2_vd_map.find(point) == point2_vd_map.end())
        {
            continue;
        }
        vertex_descriptor vd = point2_vd_map.at(point);
        convex_hull_vertices.push_back(vd);
    }
    std::cout << "Convex hull vertices: " << convex_hull_vertices.size() << std::endl;
    // Filter the convex_hull_vertices to avoid zigzagging spline cause by proximal vertices
    std::vector<vertex_descriptor> equally_distributed_vertices = EquallyDistributeVertices(m_expanded_abutment_sm, convex_hull_vertices);
    std::cout << "Equally distributed vertices: " << equally_distributed_vertices.size() << std::endl;

    auto weighted_point = [](Point_3& p1, Point_3& p2, Point_3& p3, double w1, double w2, double w3) -> Point_3
        {
            assert((1.0 - w1 - w2 - w3) < 1E-3);
            return Point_3(w1 * p1.x() + w2 * p2.x() + w3 * p3.x(),
                w1 * p1.y() + w2 * p2.y() + w3 * p3.y(),
                w1 * p1.z() + w2 * p2.z() + w3 * p3.z());
        };

    // Add the points from equally_distributed_vertices as control points to m_abutment_edge_spline 
    for (const auto& vd : equally_distributed_vertices)
    {
        //std::cout << "vd: " << vd.idx() << std::endl;
        face_descriptor fd = FindAdjacentFace(m_expanded_abutment_sm, vd);
        //std::cout << "fd: " << fd.idx() << std::endl;

        halfedge_descriptor hd = halfedge(fd, m_expanded_abutment_sm);
        //std::cout << "hd: " << hd.idx() << std::endl;
        vertex_descriptor v1 = m_expanded_abutment_sm.target(hd);
        //std::cout << "v1: " << v1.idx() << std::endl;
        vertex_descriptor v2 = m_expanded_abutment_sm.target(m_expanded_abutment_sm.next(hd));
        //std::cout << "v2: " << v2.idx() << std::endl;
        vertex_descriptor v3 = m_expanded_abutment_sm.target(m_expanded_abutment_sm.next(m_expanded_abutment_sm.next(hd)));
        //std::cout << "v3: " << v3.idx() << std::endl;
        vertex_descriptor vd1, vd2;
        MeshPoint ctrl_pt;
        ctrl_pt.nTriId = fd.idx();

        if (vd == v1)
        {
            vd1 = v2;
            vd2 = v3;
        }
        else if (vd == v2)
        {
            vd1 = v3;
            vd2 = v1;
        }
        else
        {
            vd1 = v1;
            vd2 = v2;
        }

        if (fd != SurfaceMesh::null_face())
        {
            Point_3 p1 = m_expanded_abutment_sm.point(vd);
            Point_3 p2 = m_expanded_abutment_sm.point(vd1);
            Point_3 p3 = m_expanded_abutment_sm.point(vd2);

            // Spline cannot take into the points at the vertices, so we need to add a control point slightly shifted from the vertex


            Point_3 weighted_p = weighted_point(p1, p2, p3, 0.98, 0.01, 0.01);
            ctrl_pt.xyz[0] = weighted_p.x();
            ctrl_pt.xyz[1] = weighted_p.y();
            ctrl_pt.xyz[2] = weighted_p.z();

            // Render the adjacent faces of each control point for debugging
            /*vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            points->InsertNextPoint(p1.x(), p1.y(), p1.z());
            points->InsertNextPoint(p2.x(), p2.y(), p2.z());
            points->InsertNextPoint(p3.x(), p3.y(), p3.z());

            vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
            line1->GetPointIds()->SetId(0, 0);
            line1->GetPointIds()->SetId(1, 1);

            vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
            line2->GetPointIds()->SetId(0, 1);
            line2->GetPointIds()->SetId(1, 2);

            vtkSmartPointer<vtkLine> line3 = vtkSmartPointer<vtkLine>::New();
            line3->GetPointIds()->SetId(0, 2);
            line3->GetPointIds()->SetId(1, 0);

            vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
            lines->InsertNextCell(line1);
            lines->InsertNextCell(line2);
            lines->InsertNextCell(line3);

            vtkSmartPointer<vtkPolyData> face_polydata = vtkSmartPointer<vtkPolyData>::New();
            face_polydata->SetPoints(points);
            face_polydata->SetLines(lines);
            vtkSmartPointer<vtkPolyDataMapper> face_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            face_mapper->SetInputData(face_polydata);
            vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(face_mapper);
            Renderer->AddActor(actor);*/

            //std::cout << "Add ctrl pt: face_descriptor: " << fd.idx() << std::endl;

            m_abutment_edge_spline->add(m_expanded_face_map.at(face_descriptor(ctrl_pt.nTriId)).idx(), ctrl_pt.xyz);
            /*vtkSmartPointer<vtkPolyDataNormals> vtkNormal = vtkSmartPointer<vtkPolyDataNormals>::New();
            vtkNormal->SetInputConnection(abutment_edge_spline.CtrlPointSphere.back()->GetOutputPort());
            vtkNormal->SetComputePointNormals(1);
            vtkNormal->SetComputeCellNormals(1);
            vtkNormal->SetAutoOrientNormals(1);
            vtkNormal->SetSplitting(0);
            vtkNormal->FlipNormalsOff();
            vtkNormal->Update();
            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(vtkNormal->GetOutputPort());
            mapper->Update();
            abutment_edge_spline.CtrlPointActor.back()->SetMapper(mapper);
            abutment_edge_spline.CtrlPointActor.back()->GetProperty()->SetColor(0, 1, 0);
            abutment_edge_spline.CtrlPointActor.back()->GetProperty()->SetAmbient(0.5);
            abutment_edge_spline.CtrlPointActor.back()->GetProperty()->SetSpecularPower(100);
            abutment_edge_spline.CtrlPointActor.back()->GetProperty()->SetSpecular(0.5);
            abutment_edge_spline.CtrlPointActor.back()->GetProperty()->SetDiffuse(0.5);*/
            //abutment_edge_spline.CtrlPointActor.back()->PickableOn();
            //Renderer->AddActor(abutment_edge_spline.CtrlPointActor.back());

        }
    }
    m_abutment_edge_spline->bClosed = true;
    m_abutment_edge_spline->UpdateSpline(m_abutment_edge_spline->vtEquidistantSpline);
}

/**
 * @brief Generates a cervical margin line for a dental abutment.
 *
 * This function computes the cervical margin line around a dental abutment's edge spline
 * and visualizes it using VTK (Visualization Toolkit). It creates a 3D tube representation
 * around the spline that represents the margin line.
 *
 * @pre m_abutment_edge_spline must be initialized with control points.
 * @pre m_arch_sm, m_fmap, m_vmap, m_emap, and m_hemap must be properly set.
 * @pre m_renderer and m_render_win must be properly configured.
 *
 * @post The cervical margin line is visualized in the render window.
 * @post If ENABLE_TIMER_H is defined, the function's execution time is measured.
 *
 * @par Algorithm:
 *      - Calculates the orientation of the polygon formed by the spline points.
 *      - Expands the spline using MeshSplineExpander to fit the contour of the dental model.
 *      - Visualizes the expanded spline points as spheres and adds them to the renderer.
 *      - Creates a tube around the expanded spline to represent the margin line.
 *      - Sets visual properties of the tube and adds it to the renderer.
 *      - Finalizes the interaction with OnLeftButtonUp call.
 *
 * @note This function prints the size of the edge spline to the standard output.
 * @note The function is part of the CervicalMarginLineWrapper class.
 * @note The expansion distance and other visualization parameters are currently hardcoded.
 *
 * @todo Implement expansion distance adjustment.
 *
 */
void CervicalMarginLineWrapper::GenerateCervicalMarginLine()
{
#ifdef ENABLE_TIMER_H
    Timer timer("Generate cervical margin line");
#endif
    std::cout << "Abutment Edge Spline size: " << m_abutment_edge_spline->uvSpline.size() << std::endl;
    Polygon_2 polygon(m_abutment_edge_spline->uvSpline.begin(), m_abutment_edge_spline->uvSpline.end());

    bool is_clockwise = polygon.is_clockwise_oriented();
    //uv_map = m_cervical_margin_line_interactor_style->uv_map;
    // Todo: Add the expansion distance adjustment
    MeshSplineExpander mesh_spline_expander(
        m_abutment_edge_spline->vtCtrlPoints,
        *m_arch_sm,
        m_expansion_distance,
        1,
        is_clockwise,
        m_abutment_edge_spline->vtEquidistantSpline,
        m_fmap,
        m_vmap,
        m_emap,
        m_hemap,
        m_uv_map
    );
    //MeshSplineExpander mesh_spline_expander(
    //    m_abutment_edge_spline->vtCtrlPoints,
    //    *m_arch_sm,
    //    2,
    //    is_clockwise,
    //    m_abutment_edge_spline->vtEquidistantSpline,
    //    m_fmap,
    //    m_vmap,
    //    m_emap,
    //    m_hemap,
    //    m_uv_map
    //);
    mesh_spline_expander.SetRenderer(m_renderer);
    mesh_spline_expander.SetRenderWin(m_render_win);
    bool success = mesh_spline_expander.ExpandSpline();
    //bool success = mesh_spline_expander.ExpandToLowestCurvature();
    //std::vector<ClosedMeshSpline> expanded_splines = mesh_spline_expander.GetExpandedSplines();
    //for (auto expanded_spline : expanded_splines)
    //{
    //    Renderer->AddActor(expanded_spline.SplineActor);
    //}
    std::vector<MeshPoint> expanded_spline_mp = mesh_spline_expander.GetExpandedSplinesMP()[0];
    for (MeshPoint& mp : expanded_spline_mp)
    {
        m_cervical_margin_line_interactor_style->spline->add(mp.nTriId, mp.xyz);
        vtkSmartPointer<vtkPolyDataNormals> mp_normal = vtkSmartPointer<vtkPolyDataNormals>::New();
        mp_normal->SetInputConnection(m_cervical_margin_line_interactor_style->spline->CtrlPointSphere.back()->GetOutputPort());
        mp_normal->SetComputePointNormals(1);
        mp_normal->SetComputeCellNormals(1);
        mp_normal->SetAutoOrientNormals(1);
        mp_normal->SetSplitting(0);
        mp_normal->FlipNormalsOff();
        mp_normal->Update();
        vtkSmartPointer<vtkPolyDataMapper> mp_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mp_mapper->SetInputConnection(mp_normal->GetOutputPort());
        mp_mapper->Update();
        m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->SetMapper(mp_mapper);
        m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetColor(0, 0, 1);
        m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetAmbient(0.5);
        m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetSpecularPower(100);
        m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetSpecular(0.5);
        m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetDiffuse(0.5);
        m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetOpacity(1.0);
        m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->PickableOn();
        m_renderer->AddActor(m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back());
    }
    m_cervical_margin_line_interactor_style->spline->bClosed = true;
    m_cervical_margin_line_interactor_style->spline->UpdateSpline(m_cervical_margin_line_interactor_style->spline->vtEquidistantSpline);

    vtkSmartPointer<vtkTubeFilter> margin_line_tube_filter = vtkSmartPointer<vtkTubeFilter>::New();
    margin_line_tube_filter->SetInputData(m_cervical_margin_line_interactor_style->spline->SplinePolydata);
    margin_line_tube_filter->SetRadius(0.025);
    margin_line_tube_filter->SetNumberOfSides(16);
    margin_line_tube_filter->Update();
    vtkSmartPointer<vtkPolyDataNormals> margin_line_tube_normal = vtkSmartPointer<vtkPolyDataNormals>::New();
    margin_line_tube_normal->SetInputConnection(margin_line_tube_filter->GetOutputPort());
    margin_line_tube_normal->SetComputePointNormals(1);
    margin_line_tube_normal->SetComputeCellNormals(1);
    margin_line_tube_normal->SetAutoOrientNormals(1);
    margin_line_tube_normal->SetSplitting(0);
    margin_line_tube_normal->FlipNormalsOff();
    margin_line_tube_normal->Update();
    vtkSmartPointer<vtkPolyDataMapper> abutment_line_tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    abutment_line_tube_mapper->SetInputConnection(margin_line_tube_normal->GetOutputPort());
    abutment_line_tube_mapper->Update();

    vtkSmartPointer<vtkPolyData> margin_line_tube = abutment_line_tube_mapper->GetInput();
    vtkSmartPointer<vtkPLYWriter> margin_line_tube_writer = vtkSmartPointer<vtkPLYWriter>::New();
    margin_line_tube_writer->SetFileName("margin_line.ply");
    margin_line_tube_writer->SetInputData(margin_line_tube);
    margin_line_tube_writer->Write();

    m_cervical_margin_line_interactor_style->spline->SplineActor->SetMapper(abutment_line_tube_mapper);
    m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetColor(1, 1, 0);
    m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetAmbient(0.5);
    m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetSpecularPower(100);
    m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetSpecular(0.5);
    m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetDiffuse(0.5);
    m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetOpacity(1.0);
    m_cervical_margin_line_interactor_style->spline->SplineActor->PickableOff();

    m_renderer->AddActor(m_cervical_margin_line_interactor_style->spline->SplineActor);
    m_cervical_margin_line_interactor_style->OnLeftButtonUp();
}

void CervicalMarginLineWrapper::GenerateImprovedMarginLine()
{
#ifdef ENABLE_TIMER_H
    Timer timer("Generate cervical margin line");
#endif
    std::cout << "Abutment Edge Spline size: " << m_abutment_edge_spline->uvSpline.size() << std::endl;
    Polygon_2 polygon(m_abutment_edge_spline->uvSpline.begin(), m_abutment_edge_spline->uvSpline.end());

    bool is_clockwise = polygon.is_clockwise_oriented();
    //uv_map = m_cervical_margin_line_interactor_style->uv_map;
    // Todo: Add the expansion distance adjustment
    
    //SurfaceMesh::Property_map<vertex_descriptor, double> min_curvature = m_arch_sm->add_property_map<vertex_descriptor, double>("v:min_curvature", 0.0).first;
    MeshSplineExpander mesh_spline_expander(
        m_abutment_edge_spline->vtCtrlPoints,
        *m_arch_sm,
        2.5,
        is_clockwise,
        m_abutment_edge_spline->vtEquidistantSpline,
        m_fmap,
        m_vmap,
        m_emap,
        m_hemap,
        m_uv_map
    );

    SurfaceMesh::Property_map<vertex_descriptor, double> min_curvature;
    bool success_get;
    boost::tie(min_curvature, success_get) = m_arch_sm->property_map<vertex_descriptor, double>("v:min_curvature");
    //auto mean_curvature = m_arch_sm->property_map<vertex_descriptor, double>("v:mean_curvature").first;
    std::cout << "Fetched min_curvature property map: " << success_get << std::endl;
    for (auto& vd : m_arch_sm->vertices())
    {
        std::cout << min_curvature[vd] << std::endl;

    }
    mesh_spline_expander.SetRenderer(m_renderer);
    mesh_spline_expander.SetRenderWin(m_render_win);
    bool success = mesh_spline_expander.ExpandToLowestCurvature();
    //bool success = mesh_spline_expander.ExpandToLowestCurvature();
    //std::vector<ClosedMeshSpline> expanded_splines = mesh_spline_expander.GetExpandedSplines();
    //for (auto expanded_spline : expanded_splines)
    //{
    //    Renderer->AddActor(expanded_spline.SplineActor);
    //}
    //std::vector<MeshPoint> expanded_spline_mp = mesh_spline_expander.GetExpandedSplinesMP()[0];
    //for (MeshPoint& mp : expanded_spline_mp)
    //{
    //    m_cervical_margin_line_interactor_style->spline->add(mp.nTriId, mp.xyz);
    //    vtkSmartPointer<vtkPolyDataNormals> mp_normal = vtkSmartPointer<vtkPolyDataNormals>::New();
    //    mp_normal->SetInputConnection(m_cervical_margin_line_interactor_style->spline->CtrlPointSphere.back()->GetOutputPort());
    //    mp_normal->SetComputePointNormals(1);
    //    mp_normal->SetComputeCellNormals(1);
    //    mp_normal->SetAutoOrientNormals(1);
    //    mp_normal->SetSplitting(0);
    //    mp_normal->FlipNormalsOff();
    //    mp_normal->Update();
    //    vtkSmartPointer<vtkPolyDataMapper> mp_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    //    mp_mapper->SetInputConnection(mp_normal->GetOutputPort());
    //    mp_mapper->Update();
    //    m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->SetMapper(mp_mapper);
    //    m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetColor(0, 0, 1);
    //    m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetAmbient(0.5);
    //    m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetSpecularPower(100);
    //    m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetSpecular(0.5);
    //    m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetDiffuse(0.5);
    //    m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->GetProperty()->SetOpacity(1.0);
    //    m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back()->PickableOn();
    //    m_renderer->AddActor(m_cervical_margin_line_interactor_style->spline->CtrlPointActor.back());
    //}
    //m_cervical_margin_line_interactor_style->spline->bClosed = true;
    //m_cervical_margin_line_interactor_style->spline->UpdateSpline(m_cervical_margin_line_interactor_style->spline->vtEquidistantSpline);

    //vtkSmartPointer<vtkTubeFilter> margin_line_tube_filter = vtkSmartPointer<vtkTubeFilter>::New();
    //margin_line_tube_filter->SetInputData(m_cervical_margin_line_interactor_style->spline->SplinePolydata);
    //margin_line_tube_filter->SetRadius(0.025);
    //margin_line_tube_filter->SetNumberOfSides(16);
    //margin_line_tube_filter->Update();
    //vtkSmartPointer<vtkPolyDataNormals> margin_line_tube_normal = vtkSmartPointer<vtkPolyDataNormals>::New();
    //margin_line_tube_normal->SetInputConnection(margin_line_tube_filter->GetOutputPort());
    //margin_line_tube_normal->SetComputePointNormals(1);
    //margin_line_tube_normal->SetComputeCellNormals(1);
    //margin_line_tube_normal->SetAutoOrientNormals(1);
    //margin_line_tube_normal->SetSplitting(0);
    //margin_line_tube_normal->FlipNormalsOff();
    //margin_line_tube_normal->Update();
    //vtkSmartPointer<vtkPolyDataMapper> abutment_line_tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    //abutment_line_tube_mapper->SetInputConnection(margin_line_tube_normal->GetOutputPort());
    //abutment_line_tube_mapper->Update();

    //vtkSmartPointer<vtkPolyData> margin_line_tube = abutment_line_tube_mapper->GetInput();
    //vtkSmartPointer<vtkPLYWriter> margin_line_tube_writer = vtkSmartPointer<vtkPLYWriter>::New();
    //margin_line_tube_writer->SetFileName("margin_line.ply");
    //margin_line_tube_writer->SetInputData(margin_line_tube);
    //margin_line_tube_writer->Write();

    //m_cervical_margin_line_interactor_style->spline->SplineActor->SetMapper(abutment_line_tube_mapper);
    //m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetColor(1, 1, 0);
    //m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetAmbient(0.5);
    //m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetSpecularPower(100);
    //m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetSpecular(0.5);
    //m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetDiffuse(0.5);
    //m_cervical_margin_line_interactor_style->spline->SplineActor->GetProperty()->SetOpacity(1.0);
    //m_cervical_margin_line_interactor_style->spline->SplineActor->PickableOff();

    //m_renderer->AddActor(m_cervical_margin_line_interactor_style->spline->SplineActor);
    //m_cervical_margin_line_interactor_style->OnLeftButtonUp();
}

void CervicalMarginLineWrapper::SetArchSurfaceMesh(SurfaceMesh* arch_sm)
{
	m_arch_sm = arch_sm;
}

/**
 * @brief Set the renderer for the visualization.
 *
 * @param renderer A smart pointer to a VTK renderer.
 */
void CervicalMarginLineWrapper::SetRenderer(vtkSmartPointer<vtkRenderer> renderer)
{
    m_renderer = renderer;
}

/**
 * @brief Set the render window.
 *
 * @param render_win A smart pointer to a VTK render window.
 */
void CervicalMarginLineWrapper::SetRenderWindow(vtkSmartPointer<vtkRenderWindow> render_win)
{
    m_render_win = render_win;
}

/**
 * @brief Converts a CGAL SurfaceMesh to Eigen matrix representations.
 *
 * This function takes a CGAL SurfaceMesh object and converts it into two Eigen matrices, one
 * for the vertices (V) and one for the faces (F). The vertex matrix (V) contains the x, y, and z
 * coordinates of each vertex, and the face matrix (F) contains indices into V that define each face.
 *
 * @param sm The input CGAL SurfaceMesh to convert.
 * @param[out] V An Eigen::MatrixXd that will contain the vertex coordinates after conversion.
 *               It will have as many rows as there are vertices in the mesh, and each row will
 *               have 3 columns corresponding to the x, y, and z coordinates of the vertex.
 * @param[out] F An Eigen::MatrixXi that will contain the indices of the vertices defining each face
 *               of the mesh after conversion. It will have as many rows as there are faces in the
 *               mesh, and each row will have 3 columns corresponding to the indices of the vertices
 *               that form the triangular face.
 */
void CervicalMarginLineWrapper::CGALSurfaceMeshToEigen(const SurfaceMesh& sm, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    // Get the number of vertices and faces
    size_t num_vertices = sm.number_of_vertices();
    size_t num_faces = sm.number_of_faces();

    // Initialize Eigen matrices
    V.resize(num_vertices, 3);
    F.resize(num_faces, 3);

    // Map for storing vertex indices (CGAL vertex descriptor -> continuous index)
    std::unordered_map<vertex_descriptor, Eigen::Index> vertex_indices;
    // Extract vertices
    Eigen::Index v_index = 0;
    for (auto vd : sm.vertices())
    {
        const auto& point = sm.point(vd);
        V.row(v_index) << point.x(), point.y(), point.z();
        vertex_indices[vd] = v_index++; // Save the index mapping
    }

    // Extract faces
    int f_index = 0;
    for (auto f : sm.faces())
    {
        int inner_index = 0;
        for (auto v : CGAL::vertices_around_face(sm.halfedge(f), sm))
        {
            F(f_index, inner_index++) = vertex_indices[v];
        }
        ++f_index;
    }
}

/**
 * @brief Performs a breadth-first search (BFS) to extract vertices based on mean curvature.
 *
 * This function uses a BFS approach starting from a given vertex to explore and collect vertices
 * that satisfy a curvature threshold. The mean curvature values are read from a property map
 * associated with the vertices of the mesh. The function also calculates the difference between
 * all vertices in the mesh and the extracted vertices based on curvature.
 *
 * @note If the ENABLE_TIMER_H macro is defined, the function will time its execution duration.
 *
 * @param mesh The input mesh as a SurfaceMesh object, which contains the vertices and their associated mean curvature values.
 * @param start_vd The starting vertex descriptor from which the BFS begins.
 * @param[out] extracted_vertices A vector to store the extracted vertices that meet the curvature threshold.
 * @param[out] difference_vertices A set to store the vertices that do not meet the curvature threshold.
 * @param threshold The curvature threshold used to determine whether a vertex should be extracted or not. A vertex is considered to meet the threshold if its mean curvature is greater than or equal to the negative of this threshold value.
 *
 * @return Returns true if the BFS completes without reaching a border vertex; otherwise, returns false.
 *
 * @warning If a border vertex is reached during the search, the function outputs a message to std::cout and returns false.
 */
bool CervicalMarginLineWrapper::BFSMeanCurvatureExtraction(
    SurfaceMesh& mesh,
    vertex_descriptor start_vd,
    std::vector<vertex_descriptor>& extracted_vertices,
    std::set<vertex_descriptor>& difference_vertices,
    const double& threshold
)
{
#ifdef ENABLE_TIMER_H
    Timer timer("BFSMeanCurvatureExtraction");
#endif
    extracted_vertices.clear();
    difference_vertices.clear();

    // Get the property map containing mean curvature values
    auto curvature_map = mesh.property_map<vertex_descriptor, double>("v:mean_curvature").first;

    // Define the curvature comparator
    auto CurvatureComparator = [&curvature_map, &threshold](vertex_descriptor vd)
        {
            return curvature_map[vd] >= -threshold;
        };

    // If the start vertex does not satisfy the curvature threshold, find the closest vertex that does
    if (!CurvatureComparator(start_vd))
    {
        //std::cout << "Start vertex does not satisfy curvature threshold!" << std::endl;

        std::stack<vertex_descriptor> stack;
        std::set<vertex_descriptor> visited;
        stack.push(start_vd);
        visited.insert(start_vd);

        while (!stack.empty()) {
            vertex_descriptor vd = stack.top();
            stack.pop();

            if (CurvatureComparator(vd)) {
                // Found a vertex that satisfies the curvature threshold
                start_vd = vd;
                break;
            }

            // Traverse adjacent vertices
            auto circ = mesh.vertices_around_target(mesh.halfedge(vd));
            for (auto vic = circ.begin(); vic != circ.end(); ++vic) {
                if (!visited.count(*vic)) {
                    stack.push(*vic);
                    visited.insert(*vic);
                }
            }
        }
    }
    // Vector to keep track of visited vertices - make sure vertex_descriptor can be used as an index
    std::vector<bool> visited_vertices(num_vertices(mesh), false);

    // Queue for BFS
    std::queue<vertex_descriptor> queue;

    // Start BFS from the start vertex
    queue.push(start_vd);
    visited_vertices[start_vd] = true; // Mark the start vertex as visited

    while (!queue.empty())
    {
        vertex_descriptor current_vd = queue.front();
        queue.pop();

        if (mesh.is_border(current_vd))
        {
            std::cout << "Border vertex reached!" << std::endl;
            return false;
        }

        // Check curvature and if not already added to extracted_vertices, add it
        if (CurvatureComparator(current_vd))
        {
            extracted_vertices.push_back(current_vd);

            // Iterate through out-edges to visit adjacent vertices
            for (const halfedge_descriptor& out_edge : CGAL::halfedges_around_target(current_vd, mesh))
            {
                vertex_descriptor next_vd = mesh.source(out_edge);
                if (!visited_vertices[next_vd])
                {
                    visited_vertices[next_vd] = true; // Mark as visited
                    if (CurvatureComparator(next_vd))
                    {
                        queue.push(next_vd);
                    }
                }
            }
        }
    }
    std::set<vertex_descriptor> mesh_vertex_set(mesh.vertices().begin(), mesh.vertices().end());
    std::set<vertex_descriptor> extracted_vertex_set(extracted_vertices.begin(), extracted_vertices.end());
    std::set_difference(mesh_vertex_set.begin(), mesh_vertex_set.end(), extracted_vertex_set.begin(), extracted_vertex_set.end(), std::inserter(difference_vertices, difference_vertices.begin()));
	
    return true;
}

    /**
 * @brief Extracts vertices with minimum curvature above a threshold using BFS.
 *
 * This function performs a breadth-first search (BFS) starting from a given vertex, extracting
 * vertices where the minimum curvature is above a given threshold. It uses a property map from
 * the mesh to read minimum curvature values. Additionally, it computes the set difference between
 * all vertices and the extracted ones to identify vertices not meeting the curvature criteria.
 *
 * @note If the ENABLE_TIMER_H macro is defined, the function will track the execution time.
 *
 * @param mesh The input mesh as a SurfaceMesh object, which contains the vertices and their associated minimum curvature values.
 * @param start_vd The starting vertex descriptor from which BFS will start.
 * @param[out] extracted_vertices A vector that will be populated with the vertices that meet or exceed the curvature threshold.
 * @param[out] difference_vertices A set that will be filled with vertices that do not meet the curvature threshold.
 * @param threshold The minimum curvature threshold. Vertices with minimum curvature greater than or equal to -threshold are considered.
 *
 * @warning The function outputs a message to std::cout if the start vertex does not meet the curvature threshold.
 */
void CervicalMarginLineWrapper::BFSMinCurvatureExtraction(
    SurfaceMesh& mesh,
    vertex_descriptor start_vd,
    std::vector<vertex_descriptor>& extracted_vertices,
    std::set<vertex_descriptor>& difference_vertices,
    const double& threshold
)
{
#ifdef ENABLE_TIMER_H
    Timer timer("BFSMinCurvatureExtraction");
#endif
    extracted_vertices.clear();
    difference_vertices.clear();

    // Get the property map containing mean curvature values
    auto curvature_map = mesh.property_map<vertex_descriptor, double>("v:min_curvature").first;

    // Define the curvature comparator
    auto CurvatureComparator = [&curvature_map, &threshold](vertex_descriptor vd)
        {
            return curvature_map[vd] >= -threshold;
        };

    // If the start vertex does not satisfy the curvature threshold, find the closest vertex that does
    if (!CurvatureComparator(start_vd))
    {
        std::cout << "Start vertex does not satisfy curvature threshold!" << std::endl;

        std::stack<vertex_descriptor> stack;
        std::set<vertex_descriptor> visited;
        stack.push(start_vd);
        visited.insert(start_vd);

        while (!stack.empty()) {
            vertex_descriptor vd = stack.top();
            stack.pop();

            if (CurvatureComparator(vd)) {
                // Found a vertex that satisfies the curvature threshold
                start_vd = vd;
                break;
            }

            // Traverse adjacent vertices
            auto circ = mesh.vertices_around_target(mesh.halfedge(vd));
            for (auto vic = circ.begin(); vic != circ.end(); ++vic) {
                if (!visited.count(*vic)) {
                    stack.push(*vic);
                    visited.insert(*vic);
                }
            }
        }
    }
    // Vector to keep track of visited vertices - make sure vertex_descriptor can be used as an index
    std::vector<bool> visited_vertices(num_vertices(mesh), false);

    // Queue for BFS
    std::queue<vertex_descriptor> queue;

    // Start BFS from the start vertex
    queue.push(start_vd);
    visited_vertices[start_vd] = true; // Mark the start vertex as visited

    while (!queue.empty())
    {
        vertex_descriptor current_vd = queue.front();
        queue.pop();

        // Check curvature and if not already added to extracted_vertices, add it
        if (CurvatureComparator(current_vd))
        {
            extracted_vertices.push_back(current_vd);

            // Iterate through out-edges to visit adjacent vertices
            for (const halfedge_descriptor& out_edge : CGAL::halfedges_around_target(current_vd, mesh))
            {
                vertex_descriptor next_vd = mesh.source(out_edge);
                if (!visited_vertices[next_vd])
                {
                    visited_vertices[next_vd] = true; // Mark as visited
                    if (CurvatureComparator(next_vd))
                    {
                        queue.push(next_vd);
                    }
                }
            }
        }
    }
    std::set<vertex_descriptor> mesh_vertex_set(mesh.vertices().begin(), mesh.vertices().end());
    std::set<vertex_descriptor> extracted_vertex_set(extracted_vertices.begin(), extracted_vertices.end());
    std::set_difference(mesh_vertex_set.begin(), mesh_vertex_set.end(), extracted_vertex_set.begin(), extracted_vertex_set.end(), std::inserter(difference_vertices, difference_vertices.begin()));
}

int CervicalMarginLineWrapper::PolyDataToSurfaceMesh(vtkPolyData* polyData, SurfaceMesh& surfaceMesh)
{
    vtkPoints* points = polyData->GetPoints();
    int numberOfPoints = polyData->GetNumberOfPoints();

    if (0 == numberOfPoints || points == NULL)
        return -1;

    for (int i = 0; i < numberOfPoints; i++)
    {
        double point[3];
        points->GetPoint(i, point);
        surfaceMesh.add_vertex(Kernel::Point_3(point[0], point[1], point[2]));
    }

    vtkCellArray* triangles = polyData->GetPolys();
    triangles->InitTraversal();

    //vtkIdType npts = 3, * cell;
    vtkIdType npts;               // 存储单元的点数
    vtkIdType* cell;         // 存储单元的点索引数组
    while (triangles->GetNextCell(npts, cell))
    {
        surfaceMesh.add_face(CGAL::SM_Vertex_index(cell[0]), CGAL::SM_Vertex_index(cell[1]), CGAL::SM_Vertex_index(cell[2]));
    }

    return 0;
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
SurfaceMesh* CervicalMarginLineWrapper::AreaExpander(SurfaceMesh& mesh, const vtkSmartPointer<vtkPolyData> pd, int n, std::unordered_map<face_descriptor, face_descriptor>& face_map, unsigned expansion_level = 50)
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
    for (unsigned i = 0; i < expansion_level; i++)
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
