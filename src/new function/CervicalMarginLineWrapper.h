#pragma once
#ifndef CERVICALMARGINLINEWRAPPER_H
#define CERVICALMARGINLINEWRAPPER_H
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/IO/Color.h>

#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkFloatArray.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>

#include <algorithm>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <limits>
#include <exception>

#include "ClosedSplineDesignInteractorStyle.h"
#include "MeshSplineExpander.h"
#include "Timer.hpp"
//#define ENABLE_TIMER_H

typedef CGAL::Simple_cartesian<double>											Kernel;
typedef Kernel::Point_2															Point_2;
typedef Kernel::Point_3															Point_3;
typedef Kernel::Vector_3														Vector_3;
typedef Kernel::Ray_3															Ray_3;
typedef Kernel::Segment_3														Segment_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>										SurfaceMesh;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor						vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor					halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_iterator						vertex_iterator;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor						face_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor						edge_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh>					Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive>									Traits;
typedef CGAL::AABB_tree<Traits>													Tree;
typedef Tree::Primitive_id														Primitive_id;
typedef Kernel::Ray_3 Ray_3;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray_3>::Type>		Ray_intersection;
typedef boost::optional<Tree::Intersection_and_primitive_id<Segment_3>::Type>	Segment_intersection;

/**
 * @class CervicalMarginLineWrapper
 * @brief This class is designed to handle the extraction and processing of the cervical margin line from dental abutment models.
 *
 * The CervicalMarginLineWrapper class provides functionality to extract the abutment surface mesh from a dental model,
 * generate a spline along the abutment edge, and produce a cervical margin line by expanding the boundary of the abutment.
 * It integrates with VTK for rendering and allows for various customization options such as expansion.
 *
 */
class CervicalMarginLineWrapper
{
private:
	SurfaceMesh*												m_arch_sm; ///< Pointer to the dental arch surface mesh.
	SurfaceMesh													m_abutment_sm; ///< Surface mesh of the abutment.
	vtkSmartPointer<vtkPolyData>								m_abutment_pd; ///< PolyData of the abutment.
	vtkSmartPointer<vtkPolyData>								m_arch_pd; ///< PolyData of the dental arch.
	std::unordered_map<vertex_descriptor, vertex_descriptor>	m_abutment_to_arch_vd_map; ///< Map from abutment vertices to arch vertices.
	SurfaceMesh													m_expanded_abutment_sm; ///< Expanded abutment surface mesh.
	vtkSmartPointer<vtkPolyData>								m_expanded_abutment_pd; ///< PolyData of the expanded abutment.
	std::unordered_map<face_descriptor, face_descriptor>		m_expanded_to_arch_fd_map; ///< Map of expanded abutment faces.
	std::unordered_map<vertex_descriptor, vertex_descriptor>	m_expanded_to_arch_vd_map; ///< Map of expanded abutment vertices.
	std::vector<vertex_descriptor>								m_abutment_border; ///< Vertices of the original mesh at the abutment border.
	std::vector<vertex_descriptor>								m_abutment_convex_hull_border; ///< Convex hull vertices on the abutment border.
	double														m_expansion_distance = 0.6; ///< Distance to expand the abutment border to create the margin line.
	double														m_max_detection_distance = 2.0; ///< Maximum detection distance for margin line expansion.
	Vector_3													m_projection_direction; ///< Projection direction for the margin line.
	int															m_selected_id = 0; ///< Selected identifier for internal use.
	
	vtkSmartPointer<vtkRenderer>								m_renderer; ///< VTK renderer for visualization.
	vtkSmartPointer<vtkRenderWindow>							m_render_win; ///< VTK render window for visualization.

	double														m_min_curvature_threshold = std::numeric_limits<double>::max(); ///< Minimum curvature threshold for mesh extraction.
	double														m_mean_curvature_threshold = std::numeric_limits<double>::max(); ///< Mean curvature threshold for mesh extraction.
	std::map<unsigned int, face_descriptor>						m_fmap; ///< Face map for internal use.
	std::map<unsigned int, vertex_descriptor>					m_vmap; ///< Vertex map for internal use.
	std::map<unsigned int, edge_descriptor>						m_emap; ///< Edge map for internal use.
	std::map<unsigned int, halfedge_descriptor>					m_hemap; ///< Half-edge map for internal use.

	ClosedSplineDesignInteractorStyle*							m_cervical_margin_line_interactor_style; ///< Interactor style for margin line design.
	ClosedMeshSpline*											m_abutment_edge_spline; ///< Spline along the abutment edge.
	vertex_descriptor											m_bfs_start_vd; ///< Starting vertex descriptor for BFS operations.
	double														m_ctrl_pt_density_coefficient = 0.5; ///< Control point density coefficient for spline generation [0.2, 2.0].
	Eigen::MatrixXd												m_V; ///< Eigen matrix for vertices.
	Eigen::MatrixXi												m_F; ///< Eigen matrix for faces.
	bool														m_is_initialed; ///< Flag indicating if the wrapper has been initialized.
	CGAL::Color													m_margin_line_color = CGAL::yellow(); ///< Color of the margin line.
	CGAL::Color 												m_ctrl_point_color = CGAL::blue(); ///< Color of the control points.
	double                                                      m_margin_line_opacity = 1.0; ///< Opacity of the margin line.
	double 													    m_ctrl_point_opacity = 1.0; ///< Opacity of the control points.
	
	void CGALSurfaceMeshToEigen(const SurfaceMesh& sm, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	bool BFSMeanCurvatureExtraction(
		SurfaceMesh& mesh,
		vertex_descriptor start_vd,
		std::vector<vertex_descriptor>& extracted_vertices,
		std::set<vertex_descriptor>& difference_vertices,
		const double& threshold
	);
	void BFSMinCurvatureExtraction(
		SurfaceMesh& mesh,
		vertex_descriptor start_vd,
		std::vector<vertex_descriptor>& extracted_vertices,
		std::set<vertex_descriptor>& difference_vertices,
		const double& threshold
	);
	std::pair<SurfaceMesh, std::unordered_map<vertex_descriptor, vertex_descriptor>> ExtractSubmeshAndVertexMapping(
		const SurfaceMesh& original_mesh,
		const std::vector<vertex_descriptor>& extracted_vertices
	);
	face_descriptor FindAdjacentFace(const SurfaceMesh& mesh, const vertex_descriptor& vd);
	void CalculateBorders(const SurfaceMesh& mesh, std::vector<std::vector<vertex_descriptor>>& border_vertice);
	void CalculateEdgeBorder(const SurfaceMesh& mesh, std::vector<vertex_descriptor>& border_vertices);
	bool GetBary(
		double v1[3],
		double v2[3],
		double v3[3],
		double vp[3],
		double& Bary1,
		double& Bary2,
		double& Bary3
	);
	std::map<Point_2, vertex_descriptor> ProjectToPlane(
		const SurfaceMesh& mesh,
		const std::vector<vertex_descriptor>& vertices,
		std::vector<Point_2>& projected_points
	);
	std::vector<vertex_descriptor> EquallyDistributeVertices(
		const SurfaceMesh& mesh,
		const std::vector<vertex_descriptor>& input_vertices
	);

	int PolyDataToSurfaceMesh(vtkPolyData* polyData, SurfaceMesh& surfaceMesh);
	SurfaceMesh* AreaExpander(
		SurfaceMesh& mesh,
		const vtkSmartPointer<vtkPolyData> pd,
		int n,
		std::unordered_map<face_descriptor,face_descriptor>& face_map,
		std::unordered_map<vertex_descriptor, vertex_descriptor>& vertex_map,
		unsigned expansion_level
	);

public:
	CervicalMarginLineWrapper();
	void SetProjectionDirection(const Vector_3 projection_direction);
	void SetArchSurfaceMesh(SurfaceMesh* arch_sm);
	void SetAbutmentPolyData(vtkSmartPointer<vtkPolyData> abutment_pd);
	void SetArchPolyData(vtkSmartPointer<vtkPolyData> arch_pd);
	void SetMinCurvatureThreshold(const double threshold);
	void SetMeanCurvatureThreshold(const double threshold);
	void SetCervicalMarginLineInteractorStyle(ClosedSplineDesignInteractorStyle* cervical_margin_line_interactor_style);
	void SetExpansionDistance(const double distance);
	void SetMaxDetectionDistance(const double distance);
	void SetSelectedId(const int selected_id);
	void SetCtrlPtDensityCoefficient(const double coefficient);
	void SetMarginLineColor(const CGAL::Color& color);
	void SetMarginLineColor(const double r, const double g, const double b);
	void SetCtrlPointColor(const CGAL::Color& color);
	void SetCtrlPointColor(const double r, const double g, const double b);
	void SetMarginLineOpacity(const double opacity);
	void SetCtrlPointOpacity(const double opacity);
	void Init();

	void ExtractAbutmentSurfaceMesh();
	void GenerateAbutmentEdgeSpline();
	void GenerateCervicalMarginLine();
	void GenerateImprovedMarginLine();

	void SetRenderer(vtkSmartPointer<vtkRenderer> renderer);
	void SetRenderWindow(vtkSmartPointer<vtkRenderWindow> render_win);

};

#endif // !CERVICALMARGINLINEWRAPPER_H