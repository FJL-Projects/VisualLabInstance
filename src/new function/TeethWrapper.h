#pragma once
#ifndef TEETH_WRAPPER_H
#define TEETH_WRAPPER_H
#include"stdafx.h"


/**
 * @class TeethWrapper
 *
 * @brief This class serves as a wrapper for handling operations related to teeth.
 *
 * The TeethWrapper class encapsulates SurfaceMesh objects representing the geometry of teeth, cut and bite,
 * a vtkPolyData object for visualization, and other related properties and methods. It provides methods
 * for converting between SurfaceMesh and vtkPolyData, setting and getting various properties (like
 * intersected vertices, actor, renderer, etc.), and for performing calculations related to intersection
 * of vertices and distances from bite points.
 *
 * @note It uses CGAL library for geometric operations and VTK library for visualization.
 *
 *
 */

class TeethWrapper
{
private:
	SurfaceMesh*                     m_teeth_sm;
	SurfaceMesh*					 m_bite_sm;
	SurfaceMesh*					 m_cut_sm;
	std::shared_ptr<SurfaceMesh>	 m_adjacent_teeth_sm; ///< The adjacent teeth in arch without abutment.
	SurfaceMesh*					 m_arch_sm; ///< The dental arch.
	vtkSmartPointer<vtkPolyData>     m_teeth_pd;
	vtkSmartPointer<vtkPolyData>     m_arch_pd; ///< PolyData of the dental arch.
	vtkSmartPointer<vtkActor>        m_teeth_actor;
	std::set<vertex_descriptor>      m_intersected_v_set;
	vtkSmartPointer<vtkFloatArray>   m_distance_color_array;
	vtkSmartPointer<vtkRenderer>     m_renderer;
	vtkSmartPointer<vtkRenderWindow> m_render_win;

	VNMap                            m_vertices_normals_map;
	VLMap                            m_distance_from_bite_map;
	std::shared_ptr<FaceTree>        m_bite_tree_ptr;
	std::shared_ptr<FaceTree>        m_cut_tree_ptr;
	std::shared_ptr<HalfedgeTree>    m_teeth_halfedge_tree_ptr;
	std::shared_ptr<FaceTree>		 m_adjacent_teeth_face_tree_ptr;

	Vector_3                         m_projection_direction;
	double                           m_opacity = 0.8;
	double							 m_slicing_layer_gap = 0.03;
	double							 m_max_cusp_offset_threshold = 0.0;
	std::shared_ptr<CGAL::Polygon_mesh_slicer<SurfaceMesh, Kernel>> m_teeth_slicer;
	std::vector<Polylines>			 m_all_layer_polylines;
	std::vector<vertex_descriptor>   m_cusp_vertices;
	double							 m_max_overlapping_depth = std::numeric_limits<double>::lowest();
	double							 m_max_dimension = 0.0;
	const double					 M_PI = 3.14159265358979323846;
	std::vector<Point_3>			 m_teeth_border_points;
	std::set<vertex_descriptor>      m_arch_without_abutment_set;
	int								 m_selected_id = 0;

	void SmoothMesh();
	template <typename VertexContainer>
	void SmoothMesh(
		SurfaceMesh& sm,
		const VertexContainer& area_vertices,
		const double radius = 2.0,
		const int iteration_times = 5
	);
	Vector_3 CalculateSmoothNormal(vertex_descriptor v);
	void SmoothAllNormal(VNMap& map);
	SurfaceMesh* AreaExpander(SurfaceMesh& mesh, const vtkSmartPointer<vtkPolyData> pd, int n, std::unordered_map<face_descriptor, face_descriptor>& face_map, unsigned expansion_level);
	std::pair<SurfaceMesh, std::unordered_map<vertex_descriptor, vertex_descriptor>> ExtractSubmeshAndVertexMapping(
		const SurfaceMesh& original_mesh,
		const std::vector<vertex_descriptor>& extracted_vertices
	);
	face_descriptor FindAdjacentFace(const SurfaceMesh& mesh, const vertex_descriptor& vd);

public:
	TeethWrapper();

	void SetTeethSurfaceMesh(SurfaceMesh& sm);
	SurfaceMesh& GetTeethSurfaceMesh();
	void SetBiteSurfaceMesh(SurfaceMesh& sm);
	void SetIntersectedVertexSet(const std::set<vertex_descriptor>& intersected_v_set);
	std::set<vertex_descriptor>& GetIntersectedVertexSet();
	void SetTeethActor(vtkSmartPointer<vtkActor>& actor);
	vtkSmartPointer<vtkActor>& GetTeethActor();
	std::vector<Polylines>& GetAllLayerPolylines();
	void SetRenderer(vtkSmartPointer<vtkRenderer> renderer);
	void SetRenderWindow(vtkSmartPointer<vtkRenderWindow> render_win);
	void SetSlicingLayerGap(double gap);
	void SetProjectionDirection(const Vector_3 projection_direction);
	void SetOpacity(double opacity);
	void SetMaxCuspOffsetThreshold(double threshold);
	void SetArchSurfaceMesh(SurfaceMesh& sm);
	void SetArchPolyData(vtkSmartPointer<vtkPolyData> pd);
	void SetArchWithoutAbutmentSet(const std::set<vertex_descriptor>& arch_without_abutment_set);
	void SetSelectedId(int id);

	vtkSmartPointer<vtkPolyData>& SurfaceMeshToPolyDataWithColor();
	VNMap& AddVerticesNormalsPropertyMap();
	VNMap& GetVerticesNormalsPropertyMap();
	VLMap& AddDistanceFromBitePropertyMap();
	VLMap& GetDistanceFromBitePropertyMap();
	double GetMaxDimension();
	double GetProjectionDistanceFromBite(vertex_descriptor v);
	FaceTree& GetCutTree();
	std::vector<vertex_descriptor>& GetCuspVertices();
	std::unordered_map<vertex_descriptor, double> GetVerticesDistanceWeightsMap(vertex_descriptor cusp_vertex, double radius, std::unordered_set<face_descriptor>& modified_faces);
	double GetMaxOverlappingDepth();
	std::shared_ptr<SurfaceMesh> GetAdjacentTeethSurfaceMesh();

	void RenderTeethPolydata();

	void CalculateBiteTreePtr(SurfaceMesh sm);
	void CalculateCutTreePtr(SurfaceMesh sm);
	void CalculateTeethSlicer(SurfaceMesh sm);
	void CalculateAllLayerPolylines();
	void CalculateCusp();
	void CalculateIntersectedVertexAndDistanceFromBite();
	void ExtractAdjacentTeethWithoutAbutment();

	int CountAdjacentVertices(vertex_descriptor v);
	int CountAdjacentVertices(SurfaceMesh& sm, vertex_descriptor v);
	void AdjustCrownWallThickness(double thickness);
	void AdjustCuspHeight();
	void OcclusionShaving(double offset);
	void ProximalShaving(double offset);
	void ScaleTeeth(double offset);
	bool IsSelfIntersected() const;
	void RemoveSelfIntersections();
};
#endif // !TEETH_WRAPPER_H