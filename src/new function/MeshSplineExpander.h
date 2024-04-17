#pragma once
#ifndef MESH_SPLINE_EXPANDER_H
#define MESH_SPLINE_EXPANDER_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include "Spline/ClosedMeshSpline.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>


using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using SurfaceMesh = CGAL::Surface_mesh<Kernel::Point_3>;
using vertex_descriptor = boost::graph_traits<SurfaceMesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<SurfaceMesh>::halfedge_descriptor;
using face_descriptor = boost::graph_traits<SurfaceMesh>::face_descriptor;
using edge_descriptor = boost::graph_traits<SurfaceMesh>::edge_descriptor;
using Ray_3 = Kernel::Ray_3;
using Segment_3 = Kernel::Segment_3;
using Plane_3 = Kernel::Plane_3;

/**
 * @class MeshSplineExpander
 * @brief Expands a base ClosedMeshSpline by creating a series of splines at a specified interval.
 *
 * This class takes a ClosedMeshSpline object as input and generates a number of additional
 * ClosedMeshSpline objects spaced at a given interval distance, effectively creating an
 * expanded series of splines.
 */
class MeshSplineExpander
{
private:
    ClosedMeshSpline m_closed_mesh_spline; ///< Base closed mesh spline, including control points and spline
    std::vector<MeshPoint> m_base_spline; ///< Base closed mesh spline
    std::vector<MeshPoint> m_equal_distance_spline; ///< Base closed mesh spline with equal distance of points
    SurfaceMesh m_sm; ///< Surface mesh data
    double m_interval = 0; ///< Interval between expanded splines
    size_t m_splines_quota = 0; ///< Number of splines to generate
    std::map<unsigned int, face_descriptor> m_fmap;
    std::map<unsigned int, vertex_descriptor> m_vmap;
    std::map<unsigned int, edge_descriptor> m_emap;
    std::map<unsigned int, halfedge_descriptor> m_hemap;

    std::vector<ClosedMeshSpline> m_expanded_splines; ///< Collection of expanded splines
    std::vector<std::vector<MeshPoint>> m_expanded_splines_mp; ///< Collection of MeshPoints of each expanded splines
    bool m_is_clockwise; ///< Whether the base spline is clockwise
    std::vector<std::vector<MeshPoint> > m_linked_points_vec;
    std::vector<std::vector<MeshPoint> > m_edge_points_vec;
    vtkSmartPointer<vtkRenderWindow> m_render_win = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderer> m_renderer = vtkSmartPointer<vtkRenderer>::New();
    std::vector<Vector_3> m_spline_expand_directions;  ///< Directions for spline expansion
    const double M_PI = 3.14159265358979323846;
    // std::vector<MeshPoint> m_equal_distance_spline;
    face_descriptor m_last_fd = SurfaceMesh::null_face();
    double m_max_distance = 0.0;  ///< Maximum distance for lowest curvature expansion detection
    bool m_ctrl_pts_neighbor_direction = true; ///< Whether to use the direction of the neighboring control points or use the centroid of the points
    Point_3 m_expansion_source_center; ///< Center of the expansion model

public:
    // Constructor for multiple splines with equal distance
    MeshSplineExpander(
        const ClosedMeshSpline& closed_mesh_spline,
        const SurfaceMesh& sm,
        double interval,
        const size_t splines_quota,
        const bool is_clockwise,
        const std::map<unsigned int, face_descriptor>& fmap,
        const std::map<unsigned int, vertex_descriptor>& vmap,
        const std::map<unsigned int, edge_descriptor>& emap,
        const std::map<unsigned int, halfedge_descriptor>& hemap,
        const bool ctrl_pts_neighbor_direction = true
    );

    // Constructor for single cervical margin spline
    MeshSplineExpander(
        const ClosedMeshSpline& closed_mesh_spline,
        const SurfaceMesh& sm,
        const double max_distance,
        const bool is_clockwise,
        const std::map<unsigned int, face_descriptor>& fmap,
        const std::map<unsigned int, vertex_descriptor>& vmap,
        const std::map<unsigned int, edge_descriptor>& emap,
        const std::map<unsigned int, halfedge_descriptor>& hemap,
        const bool ctrl_pts_neighbor_direction = true
    );
    // Setters
    void SetSm(const SurfaceMesh& sm);
    void SetInterval(double interval);
    void SetSplinesQuota(size_t splines_quota);
    void SetIsClockwise(bool is_clockwise);
    void SetRenderWin(vtkSmartPointer<vtkRenderWindow> render_win);
    void SetRenderer(vtkSmartPointer<vtkRenderer> renderer);
    void SetExpansionSourceCenter(const Point_3& expansion_source_center);

    // Getters
    const std::vector<MeshPoint>& GetBaseSpline() const;
    const SurfaceMesh& GetSm() const;
    double GetInterval() const;
    size_t GetSplinesQuota() const;
    const std::vector<ClosedMeshSpline>& GetExpandedSplines() const;
    const std::vector<std::vector<MeshPoint>>& GetExpandedSplinesMP() const;
    const std::vector<MeshPoint>& GetEqualDistanceSpline() const;
    const std::vector<std::vector<MeshPoint>>& GetLinkedPointsVec() const;
    Point_2 GetPointUV(const SurfaceMesh& sm, const SurfaceMesh::Property_map<vertex_descriptor, Point_2>& uv_map, const halfedge_descriptor& he, const Point_3& pt);
    bool GetBary(Point_3& v0, Point_3& v1, Point_3& v2, const Point_3& vp, double& bary0, double& bary1, double& bary2);

    void CalculateSplineExpansionDirections();
    Vector_3 CalculateDirection(const MeshPoint& source_mesh_point, size_t curr_index, size_t prev_index, size_t next_index);
    Vector_3 CalculateConcavePointDirection(
        const MeshPoint& source_mesh_point,
        const size_t& curr_index,
        const Point_3& source_point,
        const Vector_3& source_to_prev,
        const Vector_3& source_to_next,
        const Vector_3& face_normal
    );
    Vector_3 CalculateConvexPointDirection(
        const MeshPoint& source_mesh_point,
        const size_t& curr_index,
        const Point_3& source_point,
        const Vector_3& source_to_prev,
        const Vector_3& source_to_next,
        const Vector_3& face_normal
    );
    double CalculateAngleRad(const Vector_3& source_to_prev, const Vector_3& source_to_next);

    bool ExpandSpline();

    bool ExpandToLowestCurvature();

    void ExtractPointByInterval(
        const halfedge_descriptor   start_hd,
        const Point_3               source_point,
        const Plane_3& cutting_plane,
        double                remaining_interval,
        const double          required_interval,
        std::vector<MeshPoint>& equal_distance_points_vec,
        std::vector<MeshPoint>& edge_points_vec
    );
    void ExtractPoint(
        const halfedge_descriptor& hd,
        const halfedge_descriptor& start_hd,
        const Segment_3& segment,
        const Point_3& source_point,
        const Plane_3& cutting_plane,
        double                      remaining_interval,
        const double& required_interval,
        std::vector<MeshPoint>& equal_distance_points_vec,
        std::vector<MeshPoint>& edge_points_vec
    );
};
#endif // MESH_SPLINE_EXPANDER_H