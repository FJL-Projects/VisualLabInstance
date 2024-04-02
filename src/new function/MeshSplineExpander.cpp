#include <vector>
#include "Spline/ClosedMeshSpline.h" 
#include "MeshSplineExpander.h"

// Constructor
MeshSplineExpander::MeshSplineExpander(
	const ClosedMeshSpline& closed_mesh_spline,
	const SurfaceMesh& sm,
	double interval,
	const size_t splines_quota,
	const bool is_clockwise,
	const std::map<unsigned int, face_descriptor>& fmap,
	const std::map<unsigned int, vertex_descriptor>& vmap,
	const std::map<unsigned int, edge_descriptor>& emap,
	const std::map<unsigned int, halfedge_descriptor>& hemap
) : m_closed_mesh_spline(closed_mesh_spline), m_sm(sm), m_interval(interval), m_splines_quota(splines_quota), m_is_clockwise(is_clockwise),
	m_fmap(fmap), m_vmap(vmap), m_emap(emap), m_hemap(hemap), m_base_spline(closed_mesh_spline.vtCtrlPoints), m_equal_distance_spline(closed_mesh_spline.vtEquidistantSpline)
{}

MeshSplineExpander::MeshSplineExpander(
	const ClosedMeshSpline& closed_mesh_spline,
	const SurfaceMesh& sm,
	const double max_distance,
	const bool is_clockwise,
	const std::map<unsigned int, face_descriptor>& fmap,
	const std::map<unsigned int, vertex_descriptor>& vmap,
	const std::map<unsigned int, edge_descriptor>& emap,
	const std::map<unsigned int, halfedge_descriptor>& hemap
) : m_closed_mesh_spline(closed_mesh_spline), m_sm(sm), m_max_distance(max_distance), m_splines_quota(1), m_is_clockwise(is_clockwise),
m_fmap(fmap), m_vmap(vmap), m_emap(emap), m_hemap(hemap), m_base_spline(closed_mesh_spline.vtCtrlPoints), m_equal_distance_spline(closed_mesh_spline.vtEquidistantSpline)
{}

// SetSm member function
void MeshSplineExpander::SetSm(const SurfaceMesh& sm)
{
    m_sm = sm;
}

// SetInterval member function
void MeshSplineExpander::SetInterval(double interval) 
{
    m_interval = interval;
}

// SetSplinesQuota member function
void MeshSplineExpander::SetSplinesQuota(size_t splines_quota) 
{
    m_splines_quota = splines_quota;
}

void MeshSplineExpander::SetIsClockwise(bool is_clockwise)
{
    m_is_clockwise = is_clockwise;
}

void MeshSplineExpander::SetRenderWin(vtkSmartPointer<vtkRenderWindow> render_win)
{
	m_render_win = render_win;
}

void MeshSplineExpander::SetRenderer(vtkSmartPointer<vtkRenderer> renderer)
{
	m_renderer = renderer;
}

// GetBaseSpline member function
const std::vector<MeshPoint>& MeshSplineExpander::GetBaseSpline() const 
{
    return m_base_spline;
}

// GetSm member function
const SurfaceMesh& MeshSplineExpander::GetSm() const
{
    return m_sm;
}

// GetInterval member function
double MeshSplineExpander::GetInterval() const 
{
    return m_interval;
}

// GetSplinesQuota member function
size_t MeshSplineExpander::GetSplinesQuota() const 
{
    return m_splines_quota;
}

const std::vector<ClosedMeshSpline>& MeshSplineExpander::GetExpandedSplines() const
{
	return m_expanded_splines;
}

const std::vector<std::vector<MeshPoint>>& MeshSplineExpander::GetExpandedSplinesMP() const
{
	return m_expanded_splines_mp;
}

const std::vector<MeshPoint>& MeshSplineExpander::GetEqualDistanceSpline() const
{
	return m_equal_distance_spline;
}

const std::vector<std::vector<MeshPoint>>& MeshSplineExpander::GetLinkedPointsVec() const
{
	return m_linked_points_vec;
}

Point_2 MeshSplineExpander::GetPointUV(
	const SurfaceMesh& sm,
	const SurfaceMesh::Property_map<vertex_descriptor, Point_2>& uv_map,
	const halfedge_descriptor& he,
	const Point_3& pt
)
{
	vertex_descriptor vd0 = sm.target(he);
	vertex_descriptor vd1 = sm.target(sm.next(he));
	vertex_descriptor vd2 = sm.source(he);

	Point_3& p0 = get(CGAL::get(CGAL::vertex_point, sm), vd0);
	Point_3& p1 = get(CGAL::get(CGAL::vertex_point, sm), vd1);
	Point_3& p2 = get(CGAL::get(CGAL::vertex_point, sm), vd2);

	Point_2 uv[3] = { uv_map[vd0], uv_map[vd1], uv_map[vd2] };
	double bary0, bary1, bary2;
	GetBary(p0, p1, p2, pt, bary0, bary1, bary2);
	return Point_2(uv[0].x() * bary0 + uv[1].x() * bary1 + uv[2].x() * bary2, uv[0].y() * bary0 + uv[1].y() * bary1 + uv[2].y() * bary2);
}

bool MeshSplineExpander::GetBary(Point_3& v0, Point_3& v1, Point_3& v2, const Point_3& vp, double& bary0, double& bary1, double& bary2)
{
	Point_3 v02(v0.x() - v2.x(), v0.y() - v2.y(), v0.z() - v2.z());
	Point_3 v12(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
	Point_3 vp2(vp.x() - v2.x(), vp.y() - v2.y(), vp.z() - v2.z());

	double d0202 = v02.x() * v02.x() + v02.y() * v02.y() + v02.z() * v02.z();
	double d0212 = v02.x() * v12.x() + v02.y() * v12.y() + v02.z() * v12.z();
	double d1212 = v12.x() * v12.x() + v12.y() * v12.y() + v12.z() * v12.z();
	double d02p2 = v02.x() * vp2.x() + v02.y() * vp2.y() + v02.z() * vp2.z();
	double d12p2 = v12.x() * vp2.x() + v12.y() * vp2.y() + v12.z() * vp2.z();

	bary0 = (d1212 * d02p2 - d0212 * d12p2) / (d0202 * d1212 - d0212 * d0212);
	bary1 = (d0202 * d12p2 - d0212 * d02p2) / (d0202 * d1212 - d0212 * d0212);
	bary2 = static_cast<double>(1.0) - bary0 - bary1;

	return ((bary0 >= 0) && (bary1 >= 0) && (bary0 + bary1 < static_cast<double>(1.0)));
}

/**
 * @brief Computes and updates the expansion tangential vector directions for each control point on the base spline.
 *
 * This function is responsible for calculating the direction in which each control point on the base spline should be expanded.
 * It starts by clearing and resizing the container that stores the tangential vectors' directions. Then, it iterates over each control point
 * on the base spline. For each control point, it computes the indices of the previous and next control points in order to calculate the current
 * control point's expansion direction. Once the calculation is complete, the direction is saved in the member variable `m_spline_expand_directions`.
 *
 * @date 2024-02-26
 * @param None
 * @return void
 */
void MeshSplineExpander::CalculateSplineExpansionDirections()
{
    m_spline_expand_directions.clear();
    m_spline_expand_directions.resize(m_base_spline.size());
    for (size_t i = 0; i < m_base_spline.size(); ++i)
    {
        const MeshPoint& source_mesh_point = m_base_spline[i];
        int prev_index = (i - 1 + m_base_spline.size()) % m_base_spline.size();
        int next_index = (i + 1) % m_base_spline.size();

        // 计算方向
        Vector_3 direction = CalculateDirection(source_mesh_point, i, prev_index, next_index);

        // 设置方向
        m_spline_expand_directions[i] = direction;

		//Point_3 start_point(source_mesh_point.xyz[0], source_mesh_point.xyz[1], source_mesh_point.xyz[2]);
		//Point_3 end_point = start_point + direction;

		//vtkSmartPointer<vtkLineSource> line_source = vtkSmartPointer<vtkLineSource>::New();

		//// 设置射线的起点和终点
		//line_source->SetPoint1(start_point.x(), start_point.y(), start_point.z());
		//line_source->SetPoint2(end_point.x(), end_point.y(), end_point.z());
		//line_source->Update();

		//// 创建mapper
		//vtkSmartPointer<vtkPolyDataMapper> line_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		//line_mapper->SetInputConnection(line_source->GetOutputPort());

		//// 创建actor
		//vtkSmartPointer<vtkActor> line_actor = vtkSmartPointer<vtkActor>::New();
		//line_actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		//line_actor->GetProperty()->SetLineWidth(2.0);
		//line_actor->SetMapper(line_mapper);

		//// 将actor添加到渲染器中
		//m_renderer->AddActor(line_actor);
    }
}

/**
 * @brief Calculate expansion directions for each point on the base spline.
 *
 * This method computes and sets the direction vectors for expanding the spline into a mesh. It iterates
 * through each point on the base spline, calculates the direction of expansion using adjacent points, and
 * stores the result in a member variable.
 *
 * @date 2024-02-26
 */
Vector_3 MeshSplineExpander::CalculateDirection(const MeshPoint& source_mesh_point, size_t curr_index, size_t prev_index, size_t next_index)
{
    const Point_3 source_point(source_mesh_point.xyz[0], source_mesh_point.xyz[1], source_mesh_point.xyz[2]);
    const Point_3 prev_point(m_base_spline[prev_index].xyz[0], m_base_spline[prev_index].xyz[1], m_base_spline[prev_index].xyz[2]);
    const Point_3 next_point(m_base_spline[next_index].xyz[0], m_base_spline[next_index].xyz[1], m_base_spline[next_index].xyz[2]);

    Vector_3 source_to_prev = Vector_3(source_point, prev_point);
    Vector_3 prev_to_source = -source_to_prev;
    Vector_3 source_to_next = Vector_3(source_point, next_point);
    Vector_3 prev_source_next_cross = CGAL::cross_product(prev_to_source, source_to_next);

    source_to_prev /= source_to_prev.squared_length();
    source_to_next /= source_to_next.squared_length();

    Vector_3 direction = Vector_3(0, 0, 0);
    // 计算面的法线
    halfedge_descriptor hd = source_mesh_point.he;
    face_descriptor fd = m_sm.face(hd);
    Vector_3 face_normal = PMP::compute_face_normal(fd, m_sm);

    if ((!m_is_clockwise && (prev_source_next_cross * face_normal < 0))
        || (m_is_clockwise && (prev_source_next_cross * face_normal > 0)))
    {
        // 凹点处理逻辑
		//std::cout << "凸点\n";
        direction = CalculateConcavePointDirection(source_mesh_point, curr_index, source_point, source_to_prev, source_to_next, face_normal);
    }
    else
    {
        // 凸点处理逻辑
		//std::cout << "凹点\n";
        direction = CalculateConvexPointDirection(source_mesh_point, curr_index, source_point, source_to_prev, source_to_next, face_normal);
    }

    return direction;
}

/**
 * @brief Calculates the outward expansion direction vector at a concave point on a spline.
 *
 * This function uses different algorithmic strategies based on the angle formed by the concave point and its adjacent points.
 * If the angle is close to straight (between 135 and 180 degrees), it uses a neighborhood normal algorithm; 
 * if the angle is less than 135 degrees, it employs a bisector algorithm.
 * The neighborhood normal algorithm considers the order of two points on the spline equidistant from the concave point 
 * and calculates the expansion direction based on these points.
 * The bisector algorithm simply adds the vectors from the concave point to its adjacent points and projects this 
 * onto the face normal to determine the outward direction.
 *
 * @param source_mesh_point The mesh point considered as the source concave point.
 * @param curr_index The current index of the source mesh point in the spline.
 * @param source_point The 3D point representing the source mesh point.
 * @param source_to_prev Vector from the source point to the previous point in the spline.
 * @param source_to_next Vector from the source point to the next point in the spline.
 * @param face_normal The normal vector of the face that the spline lies on.
 * @return Vector_3 The calculated direction vector for the outward expansion from the concave point.
 * @date 2024-02-26
 */
Vector_3 MeshSplineExpander::CalculateConcavePointDirection(
	const MeshPoint& source_mesh_point,
	const size_t& curr_index,
	const Point_3& source_point,
	const Vector_3& source_to_prev,
	const Vector_3& source_to_next,
	const Vector_3& face_normal
)
{
	Vector_3 direction;
	double angle_rad = CalculateAngleRad(source_to_prev, source_to_next);
	if (angle_rad > 3 * M_PI / 4) // 接近平角,(135°, 180°]
	{
		//std::cout << "邻域法向算法\n";
		std::vector<std::pair<size_t, double>> equal_to_ctrl_distance;

		for (size_t k = 0; k < GetEqualDistanceSpline().size(); ++k)
		{
			MeshPoint mp2 = GetEqualDistanceSpline()[k];
			Point_3 equal_spline_point(mp2.xyz[0], mp2.xyz[1], mp2.xyz[2]);
			equal_to_ctrl_distance.push_back(std::make_pair(k, std::sqrt(CGAL::squared_distance(source_point, equal_spline_point))));
		}
		// Use lambda expression to sort by distance, which is the second element of the pair.
		std::sort(equal_to_ctrl_distance.begin(), equal_to_ctrl_distance.end(), [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) { return a.second < b.second; });
		// equal_to_ctrl_distance[0] is the control point itself
		size_t nearest_spline_point_index1 = equal_to_ctrl_distance[1].first;
		size_t nearest_spline_point_index2 = equal_to_ctrl_distance[2].first;

		// 对于每个控制点最近的两个等距样条线点，不满足顺序要求的要交换下标
		if (curr_index > 0)
		{
			// Swap
			if ((!m_is_clockwise && (nearest_spline_point_index1 > nearest_spline_point_index2))
				|| (m_is_clockwise && (nearest_spline_point_index1 < nearest_spline_point_index2)))
			{
				std::swap(nearest_spline_point_index1, nearest_spline_point_index2);
			}
		}
		else if (curr_index == 0)
		{
			MeshPoint mp1 = GetEqualDistanceSpline()[nearest_spline_point_index1];
			MeshPoint mp2 = GetEqualDistanceSpline()[nearest_spline_point_index2];
			Point_3 nearest_spline_point1(mp1.xyz[0], mp1.xyz[1], mp1.xyz[2]);
			Point_3 nearest_spline_point2(mp2.xyz[0], mp2.xyz[1], mp2.xyz[2]);
			Vector_3 source_to_point1 = Vector_3(source_point, nearest_spline_point1);
			Vector_3 source_to_point2 = Vector_3(source_point, nearest_spline_point2);

			if (source_to_point1 * source_to_point2 < 0) // 两个最近点分居0号ctrl_point点两侧, 若为逆时针，则大的指向小的
			{
				if ((!m_is_clockwise && nearest_spline_point_index1 < nearest_spline_point_index2)
					|| (m_is_clockwise && nearest_spline_point_index1 > nearest_spline_point_index2))
				{
					std::swap(nearest_spline_point_index1, nearest_spline_point_index2);
				}
			}
			else // 同侧， 逆时针则小的指向大的
			{
				if ((!m_is_clockwise && nearest_spline_point_index1 > nearest_spline_point_index2)
					|| (m_is_clockwise && nearest_spline_point_index1 < nearest_spline_point_index2))
				{
					std::swap(nearest_spline_point_index1, nearest_spline_point_index2);
				}
			}
		}
		// 因为可能交换过index所以需要重新提取MeshPoint
		MeshPoint mp1 = GetEqualDistanceSpline()[nearest_spline_point_index1];
		MeshPoint mp2 = GetEqualDistanceSpline()[nearest_spline_point_index2];
		Point_3 nearest_spline_point1(mp1.xyz[0], mp1.xyz[1], mp1.xyz[2]);
		Point_3 nearest_spline_point2(mp2.xyz[0], mp2.xyz[1], mp2.xyz[2]);
		Vector_3 point1_to_point2 = Vector_3(nearest_spline_point1, nearest_spline_point2);
		point1_to_point2 /= point1_to_point2.squared_length();
		direction = CGAL::cross_product(point1_to_point2, face_normal);
		direction /= std::sqrt(direction.squared_length());
	}
	else
	{
		// [0°, 135°]
		//std::cout << "角平分线算法\n";
		direction = source_to_prev + source_to_next;
		direction = direction - (direction * face_normal) * face_normal;
		direction /= std::sqrt(direction.squared_length());
	}
	
	return direction;
} 

/**
 * @fn Vector_3 MeshSplineExpander::CalculateConvexPointDirection(const MeshPoint& source_mesh_point, const size_t& curr_index, const Point_3& source_point, const Vector_3& source_to_prev, const Vector_3& source_to_next, const Vector_3& face_normal)
 * @brief Calculates the outward expansion direction vector at a convex point on a spline.
 * @date 2024-02-26
 *
 * @details This function computes the direction in which to move a convex point on a spline outward during an expansion process.
 * It uses a weighted bisector method for angles between the point and its neighbors that are less than or equal to 45 degrees,
 * and a neighborhood normal method for angles greater than 45 degrees.
 *
 * The weighted bisector method considers the lengths of the vectors from the convex point to its previous and next neighbors on the spline.
 * The result is projected onto the plane perpendicular to the face normal.
 *
 * The neighborhood normal method finds the two points on the spline closest to the convex point that are equidistant from it.
 * It ensures these points are in the correct order relative to the spline's direction and calculates the normal of the line segment joining these two points.
 * This normal is then projected onto the face normal to determine the outward direction, ensuring it points away from the surface.
 *
 * The function returns a normalized vector representing the direction from the convex point outward from the spline.
 *
 * @param source_mesh_point Current convex control point as a MeshPoint object.
 * @param curr_index Index of the current control point on the spline.
 * @param source_point Position of the current convex control point.
 * @param source_to_prev Vector from the current control point to the previous control point.
 * @param source_to_next Vector from the current control point to the next control point.
 * @param face_normal The normal vector of the face on which the spline resides.
 * @return Vector_3 The calculated direction vector for the outward expansion at the convex point.
 */
Vector_3 MeshSplineExpander::CalculateConvexPointDirection(
	const MeshPoint& source_mesh_point,
	const size_t& curr_index,
	const Point_3& source_point,
	const Vector_3& source_to_prev,
	const Vector_3& source_to_next,
	const Vector_3& face_normal
)
{
	Vector_3 direction;
	double angle_rad = CalculateAngleRad(source_to_prev, source_to_next);
	if (angle_rad <= M_PI / 4) // [0°, 45°]
	{
		//std::cout << "距离权重平分线算法\n";

		direction = -(source_to_prev * std::sqrt(source_to_prev.squared_length()) + source_to_next * std::sqrt(source_to_next.squared_length()));
		direction = direction - (direction * face_normal) * face_normal;
		direction /= std::sqrt(direction.squared_length());
	}
	else // (45°, 180°]
	{
		//std::cout << "邻域法向算法\n";
		std::vector<std::pair<size_t, double>> equal_to_ctrl_distance;

		for (size_t k = 0; k < GetEqualDistanceSpline().size(); ++k)
		{
			MeshPoint mp2 = GetEqualDistanceSpline()[k];
			Point_3 equal_spline_point(mp2.xyz[0], mp2.xyz[1], mp2.xyz[2]);
			equal_to_ctrl_distance.push_back(std::make_pair(k, std::sqrt(CGAL::squared_distance(source_point, equal_spline_point))));
		}
		// Use lambda expression to sort by distance, which is the second element of the pair.
		std::sort(equal_to_ctrl_distance.begin(), equal_to_ctrl_distance.end(), [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) { return a.second < b.second; });
		// equal_to_ctrl_distance[0] is the control point itself
		size_t nearest_spline_point_index1 = equal_to_ctrl_distance[1].first;
		size_t nearest_spline_point_index2 = equal_to_ctrl_distance[2].first;

		// 对于每个控制点最近的两个等距样条线点，不满足顺序要求的要交换下标
		if (curr_index > 0)
		{
			// Swap
			if ((!m_is_clockwise && (nearest_spline_point_index1 > nearest_spline_point_index2))
				|| (m_is_clockwise && (nearest_spline_point_index1 < nearest_spline_point_index2)))
			{
				std::swap(nearest_spline_point_index1, nearest_spline_point_index2);
			}
		}
		else if (curr_index == 0)
		{
			MeshPoint mp1 = GetEqualDistanceSpline()[nearest_spline_point_index1];
			MeshPoint mp2 = GetEqualDistanceSpline()[nearest_spline_point_index2];
			Point_3 nearest_spline_point1(mp1.xyz[0], mp1.xyz[1], mp1.xyz[2]);
			Point_3 nearest_spline_point2(mp2.xyz[0], mp2.xyz[1], mp2.xyz[2]);
			Vector_3 source_to_point1 = Vector_3(source_point, nearest_spline_point1);
			Vector_3 source_to_point2 = Vector_3(source_point, nearest_spline_point2);
			if (source_to_point1 * source_to_point2 < 0) // 两个最近点分居0号ctrl_point点两侧, 若为逆时针，则大的指向小的
			{
				if ((!m_is_clockwise && nearest_spline_point_index1 < nearest_spline_point_index2)
					|| (m_is_clockwise && nearest_spline_point_index1 > nearest_spline_point_index2))
				{
					std::swap(nearest_spline_point_index1, nearest_spline_point_index2);
				}
			}
			else // 同侧， 逆时针则小的指向大的
			{
				if ((!m_is_clockwise && nearest_spline_point_index1 > nearest_spline_point_index2)
					|| (m_is_clockwise && nearest_spline_point_index1 < nearest_spline_point_index2))
				{
					std::swap(nearest_spline_point_index1, nearest_spline_point_index2);
				}
			}
		}

		MeshPoint mp1 = GetEqualDistanceSpline()[nearest_spline_point_index1];
		MeshPoint mp2 = GetEqualDistanceSpline()[nearest_spline_point_index2];
		Point_3 nearest_spline_point1(mp1.xyz[0], mp1.xyz[1], mp1.xyz[2]);
		Point_3 nearest_spline_point2(mp2.xyz[0], mp2.xyz[1], mp2.xyz[2]);
		Vector_3 point1_to_point2 = Vector_3(nearest_spline_point1, nearest_spline_point2);
		point1_to_point2 /= point1_to_point2.squared_length();
		direction = CGAL::cross_product(point1_to_point2, face_normal);  // If counter-clockwise
		direction /= std::sqrt(direction.squared_length());
	}

	return direction;
}

/**
 * @brief Expand a spline by cutting a SurfaceMesh to generate equidistant points.
 *
 * This function cuts outward through each control point on the base spline against a SurfaceMesh.
 * It aims to produce equidistant points and, ultimately, connect these points to form an expanded spline.
 * The process starts by calculating the direction of the spline expansion. For each point on m_base_spline,
 * a cutting plane is defined and intersected with the triangles of the SurfaceMesh.
 * During the intersection, the sum of polyline segment lengths is calculated.
 * If the sum reaches the preset m_interval, the point is recorded as an equidistant point.
 * Also, the intersection points on the triangle edges are recorded.
 * If equidistant points cannot be successfully obtained at any control point, the function returns false.
 * If successful, a visualization pipeline is built using the VTK library, and the spline display is updated.
 * The function returns a boolean value indicating whether the expansion operation was entirely successful.
 *
 * @date 2024-02-26
 * @return bool True if the spline expansion was successful, false otherwise.
 */
bool MeshSplineExpander::ExpandSpline()
{
	bool success = false;

	SurfaceMesh::Property_map<vertex_descriptor, Point_2> uv_map;
	bool is_initialized = false;
	boost::tie(uv_map, is_initialized) = m_sm.property_map<vertex_descriptor, Point_2>("h:uv");
	if (!is_initialized)
	{
		std::cerr << "uv_map is not initialized!" << std::endl;
		return false;
	}

	CalculateSplineExpansionDirections();

	for (size_t i = 0; i < m_base_spline.size(); ++i)
	{
		//std::cout << "i: " << i << std::endl;
		const MeshPoint& source_mesh_point = m_base_spline[i];
		const Point_3 source_point(source_mesh_point.xyz[0], source_mesh_point.xyz[1], source_mesh_point.xyz[2]);
		int prev_index = (i - 1 + m_base_spline.size()) % m_base_spline.size();
		int next_index = (i + 1) % m_base_spline.size();
		const Point_3 prev_point(m_base_spline[prev_index].xyz[0], m_base_spline[prev_index].xyz[1], m_base_spline[prev_index].xyz[2]);
		const Point_3 next_point(m_base_spline[next_index].xyz[0], m_base_spline[next_index].xyz[1], m_base_spline[next_index].xyz[2]);
		std::vector<MeshPoint> equal_distance_points_vec;
		std::vector<MeshPoint> edge_points_vec;

		halfedge_descriptor hd = source_mesh_point.he;
		face_descriptor fd = m_sm.face(hd);

		//if (fd == m_last_fd)
		//{
		//	continue;
		//}
		//m_last_fd = fd;

		//std::cout << "fd: " << fd.idx() << std::endl;
		Vector_3 face_normal = PMP::compute_face_normal(fd, m_sm);

		Vector_3 spline_direction = m_spline_expand_directions[i];
		Vector_3 ray_direction = spline_direction - (spline_direction * face_normal) * face_normal;
		Vector_3 cutting_plane_normal = CGAL::cross_product(face_normal, ray_direction);
		cutting_plane_normal /= std::sqrt(cutting_plane_normal.squared_length());
		Plane_3 cutting_plane(source_point, cutting_plane_normal);
		
		//Point_3 start_point(source_mesh_point.xyz[0], source_mesh_point.xyz[1], source_mesh_point.xyz[2]);
		//Point_3 end_point = start_point + face_normal;

		//vtkSmartPointer<vtkLineSource> line_source = vtkSmartPointer<vtkLineSource>::New();

		//// 设置射线的起点和终点
		//line_source->SetPoint1(start_point.x(), start_point.y(), start_point.z());
		//line_source->SetPoint2(end_point.x(), end_point.y(), end_point.z());
		//line_source->Update();

		//// 创建mapper
		//vtkSmartPointer<vtkPolyDataMapper> line_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		//line_mapper->SetInputConnection(line_source->GetOutputPort());

		//// 创建actor
		//vtkSmartPointer<vtkActor> line_actor = vtkSmartPointer<vtkActor>::New();
		//line_actor->GetProperty()->SetColor(0.0, 0.0, 1.0);
		//line_actor->GetProperty()->SetLineWidth(2.0);
		//line_actor->SetMapper(line_mapper);

		//m_renderer->AddActor(line_actor);

		// 试相交运算，检测交点与ray_direction是否同向，如果不同向则调整起始hd为下一条边。
		halfedge_descriptor next_hd = next(hd, m_sm);
		halfedge_descriptor next_next_hd = next(next_hd, m_sm);
		Segment_3 next_segment = Segment_3(m_sm.point(source(next_hd, m_sm)), m_sm.point(target(next_hd, m_sm)));
		Segment_3 next_next_segment = Segment_3(m_sm.point(source(next_next_hd, m_sm)), m_sm.point(target(next_next_hd, m_sm)));

		if (CGAL::do_intersect(cutting_plane, next_segment))
		{
			auto intersection = CGAL::intersection(cutting_plane, next_segment);
			if (Point_3* edge_point = boost::get<Point_3>(&*intersection))
			{
				Vector_3 source_intersection_vector = Vector_3(source_point, *edge_point);
				if (source_intersection_vector * ray_direction < 0)
				{
					hd = next_hd;
				}
			}
		}
		else if (CGAL::do_intersect(cutting_plane, next_next_segment))
		{
			auto intersection = CGAL::intersection(cutting_plane, next_next_segment);
			if (Point_3* edge_point = boost::get<Point_3>(&*intersection))
			{
				Vector_3 source_intersection_vector = Vector_3(source_point, *edge_point);
				if (source_intersection_vector * ray_direction < 0)
				{
					hd = next_next_hd;
				}
			}
		}

		ExtractPointByInterval(
			hd,
			source_point,
			cutting_plane,
			m_interval,
			m_interval,
			equal_distance_points_vec,
			edge_points_vec
		);

		//if (equal_distance_points_vec.size() == 0) // 取点失败
		//{
		//	return false;
		//}
		success = true;

		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkPolyLine> poly_line = vtkSmartPointer<vtkPolyLine>::New();

		//edge_points_vec.clear();
		for (auto mp : equal_distance_points_vec)
		{
			if (mp.nTriId > m_sm.number_of_faces())
			{
				std::cout << "Exceeded" << std::endl;
			}

		}
		m_linked_points_vec.push_back(equal_distance_points_vec);
		m_edge_points_vec.push_back(edge_points_vec);
	}

	m_expanded_splines.resize(m_splines_quota);
	m_expanded_splines_mp.resize(m_splines_quota);
	for (size_t i = 0; i < m_splines_quota; ++i)
	{
		m_expanded_splines[i].initial(&m_sm, m_fmap, m_vmap, m_emap, m_hemap, uv_map);
		for (size_t j = 0; j < m_base_spline.size(); ++j)
		{
			if (m_linked_points_vec[j].size() < m_splines_quota)
			{
				continue;
			}
			MeshPoint mp = m_linked_points_vec[j][i];
			m_expanded_splines_mp[i].push_back(mp);
			m_expanded_splines[i].add(mp.nTriId, mp.xyz);
		}
		m_expanded_splines[i].bClosed = true;
		m_expanded_splines[i].UpdateSpline(m_expanded_splines[i].vtEquidistantSpline);

		vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
		tubeFilter->SetInputData(m_expanded_splines[i].SplinePolydata);
		tubeFilter->SetRadius(0.05);
		tubeFilter->SetNumberOfSides(16);
		tubeFilter->Update();
		vtkSmartPointer<vtkPolyDataNormals> vtkNormal1 = vtkSmartPointer<vtkPolyDataNormals>::New();
		vtkNormal1->SetInputConnection(tubeFilter->GetOutputPort());
		vtkNormal1->SetComputePointNormals(1);
		vtkNormal1->SetComputeCellNormals(1);
		vtkNormal1->SetAutoOrientNormals(1);
		vtkNormal1->SetSplitting(0);
		vtkNormal1->FlipNormalsOn();
		vtkNormal1->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper1->SetInputConnection(vtkNormal1->GetOutputPort());
		mapper1->Update();
		m_expanded_splines[i].SplineActor->SetMapper(mapper1);
		m_expanded_splines[i].SplineActor->GetProperty()->SetColor(1, 1, 0);
		m_expanded_splines[i].SplineActor->GetProperty()->SetAmbient(0.5);
		m_expanded_splines[i].SplineActor->GetProperty()->SetSpecularPower(100);
		m_expanded_splines[i].SplineActor->GetProperty()->SetSpecular(0.5);
		m_expanded_splines[i].SplineActor->GetProperty()->SetDiffuse(0.5);
		m_expanded_splines[i].SplineActor->PickableOff();
	}

	return true;
}

bool MeshSplineExpander::ExpandToLowestCurvature()
{
	bool success = false;

	bool is_initialized = false;

	SurfaceMesh::Property_map<vertex_descriptor, Point_2> uv_map;
	boost::tie(uv_map, is_initialized) = m_sm.property_map<vertex_descriptor, Point_2>("h:uv");
	if (!is_initialized)
	{
		std::cerr << "uv_map is not initialized!" << std::endl;
		return false;
	}

	bool min_curvature_is_initialized = false;
	SurfaceMesh::Property_map<vertex_descriptor, double> min_curvature;
	boost::tie(min_curvature, min_curvature_is_initialized) = m_sm.property_map<vertex_descriptor, double>("v:min_curvature");
	if (!min_curvature_is_initialized)
	{
		std::cerr << "min_curvature is not initialized!" << std::endl;
		return false;
	}

	CalculateSplineExpansionDirections();

	for (size_t i = 0; i < m_base_spline.size(); ++i)
	{
		//std::cout << "i: " << i << std::endl;
		const MeshPoint& source_mesh_point = m_base_spline[i];
		const Point_3 source_point(source_mesh_point.xyz[0], source_mesh_point.xyz[1], source_mesh_point.xyz[2]);
		int prev_index = (i - 1 + m_base_spline.size()) % m_base_spline.size();
		int next_index = (i + 1) % m_base_spline.size();
		const Point_3 prev_point(m_base_spline[prev_index].xyz[0], m_base_spline[prev_index].xyz[1], m_base_spline[prev_index].xyz[2]);
		const Point_3 next_point(m_base_spline[next_index].xyz[0], m_base_spline[next_index].xyz[1], m_base_spline[next_index].xyz[2]);
		std::vector<MeshPoint> equal_distance_points_vec;
		std::vector<MeshPoint> edge_points_vec;

		halfedge_descriptor hd = source_mesh_point.he;
		face_descriptor fd = m_sm.face(hd);

		//if (fd == m_last_fd)
		//{
		//	continue;
		//}
		//m_last_fd = fd;

		//std::cout << "fd: " << fd.idx() << std::endl;
		Vector_3 face_normal = PMP::compute_face_normal(fd, m_sm);

		Vector_3 spline_direction = m_spline_expand_directions[i];
		Vector_3 ray_direction = spline_direction - (spline_direction * face_normal) * face_normal;
		Vector_3 cutting_plane_normal = CGAL::cross_product(face_normal, ray_direction);
		cutting_plane_normal /= std::sqrt(cutting_plane_normal.squared_length());
		Plane_3 cutting_plane(source_point, cutting_plane_normal);

		//Point_3 start_point(source_mesh_point.xyz[0], source_mesh_point.xyz[1], source_mesh_point.xyz[2]);
		//Point_3 end_point = start_point + face_normal;

		//vtkSmartPointer<vtkLineSource> line_source = vtkSmartPointer<vtkLineSource>::New();

		//// 设置射线的起点和终点
		//line_source->SetPoint1(start_point.x(), start_point.y(), start_point.z());
		//line_source->SetPoint2(end_point.x(), end_point.y(), end_point.z());
		//line_source->Update();

		//// 创建mapper
		//vtkSmartPointer<vtkPolyDataMapper> line_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		//line_mapper->SetInputConnection(line_source->GetOutputPort());

		//// 创建actor
		//vtkSmartPointer<vtkActor> line_actor = vtkSmartPointer<vtkActor>::New();
		//line_actor->GetProperty()->SetColor(0.0, 0.0, 1.0);
		//line_actor->GetProperty()->SetLineWidth(2.0);
		//line_actor->SetMapper(line_mapper);

		//m_renderer->AddActor(line_actor);

		// 试相交运算，检测交点与ray_direction是否同向，如果不同向则调整起始hd为下一条边。
		halfedge_descriptor next_hd = next(hd, m_sm);
		halfedge_descriptor next_next_hd = next(next_hd, m_sm);
		Segment_3 next_segment = Segment_3(m_sm.point(source(next_hd, m_sm)), m_sm.point(target(next_hd, m_sm)));
		Segment_3 next_next_segment = Segment_3(m_sm.point(source(next_next_hd, m_sm)), m_sm.point(target(next_next_hd, m_sm)));

		if (CGAL::do_intersect(cutting_plane, next_segment))
		{
			auto intersection = CGAL::intersection(cutting_plane, next_segment);
			if (Point_3* edge_point = boost::get<Point_3>(&*intersection))
			{
				Vector_3 source_intersection_vector = Vector_3(source_point, *edge_point);
				if (source_intersection_vector * ray_direction < 0)
				{
					hd = next_hd;
				}
			}
		}
		else if (CGAL::do_intersect(cutting_plane, next_next_segment))
		{
			auto intersection = CGAL::intersection(cutting_plane, next_next_segment);
			if (Point_3* edge_point = boost::get<Point_3>(&*intersection))
			{
				Vector_3 source_intersection_vector = Vector_3(source_point, *edge_point);
				if (source_intersection_vector * ray_direction < 0)
				{
					hd = next_next_hd;
				}
			}
		}

		ExtractPointByInterval(
			hd,
			source_point,
			cutting_plane,
			m_max_distance,
			m_max_distance,
			equal_distance_points_vec,
			edge_points_vec
		);

		//for (auto mp : edge_points_vec)
		//{
		//	std::cout << mp.nTriId << " ";
		//}
		//std::cout << std::endl;
		//if (equal_distance_points_vec.size() == 0) // 取点失败
		//{
		//	return false;
		//}
		success = true;

		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkPolyLine> poly_line = vtkSmartPointer<vtkPolyLine>::New();

		//edge_points_vec.clear();
		for (auto mp : equal_distance_points_vec)
		{
			if (mp.nTriId > m_sm.number_of_faces())
			{
				std::cout << "Exceeded" << std::endl;
			}

		}
		m_linked_points_vec.push_back(equal_distance_points_vec);
		m_edge_points_vec.push_back(edge_points_vec);
		//std::cout << i << " edge_points_vec.size(): " << edge_points_vec.size() << std::endl;
	}



	//m_expanded_splines[0].initial(&m_sm, m_fmap, m_vmap, m_emap, m_hemap, uv_map);

	std::vector<std::pair<face_descriptor, double> > face_curvature_sum_vec;
	for (size_t i = 0; i < m_base_spline.size(); ++i)
	{
		auto edge_points_vec = m_edge_points_vec[i];
		//std::cout << i << " edge_points_vec.size(): " << edge_points_vec.size() << std::endl;
		std::map<face_descriptor, double> face_curvature_sum;
		double lowest_curvature = std::numeric_limits<double>::max();
		face_descriptor lowest_curvature_fd;
		for (size_t j = 0; j < edge_points_vec.size(); ++j)
		{
			halfedge_descriptor traversed_hd = edge_points_vec[j].he;
			face_descriptor traversed_fd = m_sm.face(traversed_hd);
			//std::cout << traversed_hd.idx() << " " << traversed_fd.idx() << " " << edge_points_vec[j].nTriId << std::endl;
			double curvature_sum = 0;
			if (face_curvature_sum.find(traversed_fd) == face_curvature_sum.end()) // Face is not included
			{
				halfedge_descriptor hd = traversed_hd;
				curvature_sum += min_curvature[m_sm.source(hd)];
				curvature_sum += min_curvature[m_sm.target(hd)];
				hd = next(hd, m_sm);
				curvature_sum += min_curvature[m_sm.target(hd)];

				face_curvature_sum[traversed_fd] = curvature_sum;
				if (curvature_sum < lowest_curvature)
				{
					lowest_curvature = curvature_sum;
					lowest_curvature_fd = traversed_fd;
				}
			}

			//std::cout << edge_points_vec[j].nTriId << " ";
		}
		face_curvature_sum_vec.push_back(std::make_pair(lowest_curvature_fd, lowest_curvature));
	}

	std::vector<Point_3> lowest_curvature_points;

	auto weighted_point = [](Point_3& p1, Point_3& p2, Point_3& p3, double w1, double w2, double w3) -> Point_3
		{
			assert(std::fabs(1.0 - w1 - w2 - w3) < 1E-3);
			return Point_3(w1 * p1.x() + w2 * p2.x() + w3 * p3.x(),
				w1 * p1.y() + w2 * p2.y() + w3 * p3.y(),
				w1 * p1.z() + w2 * p2.z() + w3 * p3.z());
		};

	for (size_t i = 0; i < face_curvature_sum_vec.size(); ++i)
	{
		auto pair = face_curvature_sum_vec[i];

		const face_descriptor& fd = pair.first;
		const double& curvature_sum = pair.second;
		if (std::fabs(curvature_sum) < std::numeric_limits<double>::epsilon())
		{
			continue;
		}
		
		halfedge_descriptor hd = m_sm.halfedge(fd);
		vertex_descriptor v1 = m_sm.source(hd);
		vertex_descriptor v2 = m_sm.target(hd);
		vertex_descriptor v3 = m_sm.target(next(hd, m_sm));

		Point_3 p1 = m_sm.point(v1);
		Point_3 p2 = m_sm.point(v2);
		Point_3 p3 = m_sm.point(v3);

		double w1 = min_curvature[v1] / curvature_sum;
		double w2 = min_curvature[v2] / curvature_sum;
		double w3 = min_curvature[v3] / curvature_sum;

		auto lowest_curvature_point = weighted_point(p1, p2, p3, w1, w2, w3);
		lowest_curvature_points.push_back(lowest_curvature_point);
		
	}

	for (size_t i = 0; i < lowest_curvature_points.size(); ++i)
	{
		size_t prev_index = (i + face_curvature_sum_vec.size() - 1) % face_curvature_sum_vec.size();
		size_t next_index = (i + 1) % face_curvature_sum_vec.size();

		Point_3 prev_point = lowest_curvature_points[prev_index];
		Point_3 curr_point = lowest_curvature_points[i];
		Point_3 next_point = lowest_curvature_points[next_index];

		Vector_3 curr_to_prev = Vector_3(curr_point, prev_point);
		Vector_3 curr_to_next = Vector_3(curr_point, next_point);

		double angle = CalculateAngleRad(curr_to_prev, curr_to_next) / CGAL_PI * 180;
		std::cout << "Angle: " << angle << std::endl;

		// Render a blue sphere at hte lowest curvature point
		vtkSmartPointer<vtkSphereSource> sphere_source = vtkSmartPointer<vtkSphereSource>::New();
		sphere_source->SetCenter(curr_point.x(), curr_point.y(), curr_point.z());
		sphere_source->SetRadius(0.2);
		sphere_source->SetPhiResolution(16);
		sphere_source->SetThetaResolution(16);
		sphere_source->Update();

		vtkSmartPointer<vtkPolyDataMapper> sphere_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		sphere_mapper->SetInputConnection(sphere_source->GetOutputPort());

		double r, g, b;
		r = 0.0;
		g = angle / 180;
		b = g;

		vtkSmartPointer<vtkActor> sphere_actor = vtkSmartPointer<vtkActor>::New();
		sphere_actor->SetMapper(sphere_mapper);
		if (angle < 125.0)
		{
			sphere_actor->GetProperty()->SetColor(1, 0, 0);
		}
		else
		{
			sphere_actor->GetProperty()->SetColor(r, g, b);
		}
		sphere_actor->GetProperty()->SetAmbient(0.5);
		sphere_actor->GetProperty()->SetSpecularPower(100);
		sphere_actor->GetProperty()->SetSpecular(0.5);
		sphere_actor->GetProperty()->SetDiffuse(0.5);
		sphere_actor->GetProperty()->SetOpacity(1.0);
		sphere_actor->PickableOff();

		m_renderer->AddActor(sphere_actor);
	}

	//m_expanded_splines[0].bClosed = true;
	//m_expanded_splines[0].UpdateSpline(m_expanded_splines[0].vtEquidistantSpline);
	

	/*m_expanded_splines.resize(m_splines_quota);
	m_expanded_splines_mp.resize(m_splines_quota);
	for (size_t i = 0; i < m_splines_quota; ++i)
	{
		m_expanded_splines[i].initial(&m_sm, m_fmap, m_vmap, m_emap, m_hemap, uv_map);
		for (size_t j = 0; j < m_base_spline.size(); ++j)
		{
			MeshPoint mp = m_linked_points_vec[j][i];
			m_expanded_splines_mp[i].push_back(mp);
			m_expanded_splines[i].add(mp.nTriId, mp.xyz);
		}
		m_expanded_splines[i].bClosed = true;
		m_expanded_splines[i].UpdateSpline(m_expanded_splines[i].vtEquidistantSpline);

		vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
		tubeFilter->SetInputData(m_expanded_splines[i].SplinePolydata);
		tubeFilter->SetRadius(0.05);
		tubeFilter->SetNumberOfSides(16);
		tubeFilter->Update();
		vtkSmartPointer<vtkPolyDataNormals> vtkNormal1 = vtkSmartPointer<vtkPolyDataNormals>::New();
		vtkNormal1->SetInputConnection(tubeFilter->GetOutputPort());
		vtkNormal1->SetComputePointNormals(1);
		vtkNormal1->SetComputeCellNormals(1);
		vtkNormal1->SetAutoOrientNormals(1);
		vtkNormal1->SetSplitting(0);
		vtkNormal1->FlipNormalsOn();
		vtkNormal1->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper1->SetInputConnection(vtkNormal1->GetOutputPort());
		mapper1->Update();
		m_expanded_splines[i].SplineActor->SetMapper(mapper1);
		m_expanded_splines[i].SplineActor->GetProperty()->SetColor(1, 1, 0);
		m_expanded_splines[i].SplineActor->GetProperty()->SetAmbient(0.5);
		m_expanded_splines[i].SplineActor->GetProperty()->SetSpecularPower(100);
		m_expanded_splines[i].SplineActor->GetProperty()->SetSpecular(0.5);
		m_expanded_splines[i].SplineActor->GetProperty()->SetDiffuse(0.5);
		m_expanded_splines[i].SplineActor->PickableOff();
	}*/

	return true;
}

/**
 * @brief Extract equidistant points by intersecting with a cutting plane starting from a specified half-edge.
 *
 * This function starts from a given half-edge and detects intersections with a provided cutting plane,
 * computing points equidistant from the original spline. The workflow includes checking if the starting
 * half-edge is a boundary half-edge, in which case the process terminates. It then checks if the number
 * of equidistant points found satisfies the requirement, terminating if it does. The function continues
 * to test the intersection of the cutting plane with the next two edges of the starting half-edge. If
 * an intersection is found, the ExtractPoint function is used for further processing. If neither edge
 * intersects, the function ends. Note that the starting half-edge is not considered for intersection
 * by default. The function modifies equal_distance_points_vec and edge_points_vec to return results
 * instead of returning a value directly.
 *
 * @date 2024-02-26
 * @param start_hd The starting half-edge descriptor from which to begin detection.
 * @param source_point The starting point, used as a reference for calculating the position of equidistant points.
 * @param cutting_plane The cutting plane used for intersection tests with mesh edges.
 * @param remaining_interval The remaining distance to the next equidistant point.
 * @param required_interval The predetermined distance between two equidistant points.
 * @param equal_distance_points_vec Vector of equidistant points to store the computed points.
 * @param edge_points_vec Vector of edge points to store the computed intersection points on the edges.
 */
void MeshSplineExpander::ExtractPointByInterval(
	const halfedge_descriptor   start_hd,
	const Point_3               source_point,
	const Plane_3&			cutting_plane,
	double                remaining_interval,
	const double          required_interval,
	std::vector<MeshPoint>& equal_distance_points_vec,
	std::vector<MeshPoint>& edge_points_vec
)
{
	if (m_sm.is_border(start_hd)) // 到达边界
	{
		std::cout << "Border reached\n";
		return;
	}
	if (m_splines_quota == equal_distance_points_vec.size())  // 找到足够的点了
	{
		return;
	}

	face_descriptor fd = m_sm.face(start_hd);
	Vector_3 face_normal = PMP::compute_face_normal(fd, m_sm);

	// 寻找交点和距离
	Point_3 destination_point;

	// 对start_hd的下一条边next_hd，计算与射线的交点
	halfedge_descriptor next_hd = next(start_hd, m_sm);
	auto v1 = m_sm.point(source(next_hd, m_sm));
	auto v2 = m_sm.point(target(next_hd, m_sm));
	auto segment = Segment_3(v1, v2);

	if (CGAL::do_intersect(cutting_plane, segment))
	{
		//std::cout << "Intersected with edge 1\n";
		ExtractPoint(
			next_hd,
			start_hd,
			segment,
			source_point,
			cutting_plane,
			remaining_interval,
			required_interval,
			equal_distance_points_vec,
			edge_points_vec
		);
	}
	else // 不相交，检查下一条边
	{
		next_hd = next(next_hd, m_sm);
		v1 = m_sm.point(source(next_hd, m_sm));
		v2 = m_sm.point(target(next_hd, m_sm));
		auto segment = Segment_3(v1, v2);

		if (CGAL::do_intersect(cutting_plane, segment))
		{
			//std::cout << "Intersected with edge 2\n";

			ExtractPoint(
				next_hd,
				start_hd,
				segment,
				source_point,
				cutting_plane,
				remaining_interval,
				required_interval,
				equal_distance_points_vec,
				edge_points_vec
			);
		}
	}
}

/**
 * @brief Intersects with a cutting plane starting from a specified half-edge and calculates equidistant points along the original spline.
 *
 * The function `ExtractPoint` is designed to initiate from a specific half-edge and detect intersections
 * with a given cutting plane to compute points that are equidistant from the original spline curve. The function
 * does not return a value but modifies the `equal_distance_points_vec` and `edge_points_vec` vectors to
 * provide the results. Firstly, it checks if the starting half-edge is on a boundary, in which case the
 * function returns immediately. If the size of the equidistant points vector has already reached the quota,
 * it will also return directly. The function then calculates the normal vector of the face associated with
 * the half-edge and proceeds to test for intersections with the next edge. If an intersection is found,
 * the `ExtractPoint` function is invoked to handle the intersection point. If the next edge does not intersect,
 * the function will check the following edge. In any case of intersection, `ExtractPoint` is called to process
 * those intersection points. The function does not return any values directly but returns the computation
 * results through reference parameters.
 *
 * @date 2024-02-26
 * @param start_hd The starting half-edge descriptor from which detection begins.
 * @param source_point The reference point for calculating the position of equidistant points.
 * @param cutting_plane The cutting plane used for intersection tests with mesh edges.
 * @param remaining_interval The remaining distance to the next equidistant point.
 * @param required_interval The predetermined distance between two equidistant points.
 * @param equal_distance_points_vec A vector to store the computed equidistant points.
 * @param edge_points_vec A vector to store the computed intersection points on the edges.
 */
void MeshSplineExpander::ExtractPoint(const halfedge_descriptor& hd,
	const halfedge_descriptor& start_hd,
	const Segment_3& segment,
	const Point_3& source_point,
	const Plane_3& cutting_plane,
	double remaining_interval,
	const double& required_interval,
	std::vector<MeshPoint>& equal_distance_points_vec,
	std::vector<MeshPoint>& edge_points_vec
)
{
	auto intersection = CGAL::intersection(cutting_plane, segment);
	face_descriptor fd = m_sm.face(hd);
	Vector_3 face_normal = PMP::compute_face_normal(fd, m_sm);
	auto& uv_map = m_sm.property_map<vertex_descriptor, Point_2>("h:uv").first;

	if (Point_3* edge_point = boost::get<Point_3>(&*intersection))
	{
		double distance = std::sqrt(CGAL::squared_distance(source_point, *edge_point));

		if (distance < remaining_interval)
		{
			// 距离不够取得下一个等距点
			remaining_interval -= distance;
			double point_data[3] = { edge_point->x(), edge_point->y(), edge_point->z() };
			MeshPoint mesh_point(static_cast<unsigned>(fd.idx()), GetPointUV(m_sm, uv_map, hd, *edge_point), point_data, hd);
			edge_points_vec.push_back(mesh_point);
			//std::cout << "Edge point: " << *edge_point << std::endl;

			halfedge_descriptor opposite_hd = opposite(hd, m_sm);

			// 需要对相交边的对边所在三角面进行相交运算。所以opposite_hd是递归起始边。
			ExtractPointByInterval(
				opposite_hd,
				*edge_point,
				cutting_plane,
				remaining_interval,
				required_interval,
				equal_distance_points_vec,
				edge_points_vec
			);
		}
		else
		{
			// 距离够了，取得等距点
			Vector_3 ray_direction = Vector_3(source_point, *edge_point);
			Vector_3 cross_product = -CGAL::cross_product(face_normal, cutting_plane.orthogonal_vector());
			ray_direction = ray_direction / std::sqrt(ray_direction.squared_length());
			Point_3 destination_point = source_point + ray_direction * remaining_interval;

			double point_data[3] = { destination_point.x(), destination_point.y(), destination_point.z() };
			MeshPoint mesh_point(static_cast<unsigned>(fd.idx()), GetPointUV(m_sm, uv_map, hd, destination_point), point_data, hd);
			equal_distance_points_vec.push_back(mesh_point);
			//std::cout << "Point added: " << destination_point << std::endl;

			Vector_3 next_v(destination_point, *edge_point);
			auto next_v_len = std::sqrt(next_v.squared_length());

			// 重置remaining_interval，同时也要传入start_hd，因为继续在这个三角面进行相交运算，仍然默认忽略这条边的相交情况。
			ExtractPointByInterval(
				start_hd,
				destination_point,
				cutting_plane,
				required_interval,
				required_interval,
				equal_distance_points_vec,
				edge_points_vec
			);
		}
	}
}

/**
 * @brief Calculates the angle between two vectors in radians.
 *
 * This function computes the angle between two vectors in radians. It begins by calculating the cross product and dot product of the two vectors.
 * Then, the cosine of the angle is computed using the dot product and the magnitudes of both vectors. To prevent numerical errors that could arise from
 * the cosine value falling outside its defined domain of [-1, 1], the cosine value is clamped within this range. Finally, the arccosine function `acos()`
 * is used to obtain the angle in radians. The function has no side effects and returns the calculated angle value.
 *
 * @param source_to_prev The vector pointing from the source to the previous point.
 * @param source_to_next The vector pointing from the source to the next point.
 * @return double The angle between the two vectors, in radians.
 *
 * @date 2024-02-26
 */double MeshSplineExpander::CalculateAngleRad(const Vector_3& source_to_prev, const Vector_3& source_to_next)
{
	Vector_3 cross_product = CGAL::cross_product(source_to_prev, source_to_next);
	double dot_product = source_to_prev * source_to_next;
	double cos_angle = dot_product / (CGAL::sqrt(source_to_prev.squared_length()) * CGAL::sqrt(source_to_next.squared_length()));

	// Clamp the cosine value to the range [-1, 1] to avoid numerical issues
	cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
	double angle_rad = std::acos(cos_angle);

	return angle_rad;
}



