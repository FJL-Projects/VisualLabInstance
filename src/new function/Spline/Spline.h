#pragma once
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

#include"vectorAlgorithm.h"
#include <stdio.h>
#include <tchar.h>
#include <mutex>
#include <thread>
#include <ctime>
#include <set>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include<vtkLineSource.h>
#include <vtkNamedColors.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkNamedColors.h>
#include <vtkCubeAxesActor.h>
#include <vtkJPEGReader.h>
#include <vtkImageReader2Factory.h>
#include <vtkTexture.h>
#include <map>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkInformation.h>
#include <vtkPolyDataNormals.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkXMLDataParser.h>
#include <vtkXMLDataElement.h>
#include <vtkXMLReader.h>
#include <vtkCellData.h>
#include <vtkDataArraySelection.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkDICOMImageReader.h>
#include <algorithm>
#include <vtkCell.h>
#include <vtkTriangle.h>
#include <vtkLineSource.h>
#include <vtkOBBTree.h>
#include <vtkGlyph3DMapper.h>
#include <vtkKochanekSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricSpline.h>
#include <vtkSphereSource.h>
#include <vtkPlaneSource.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkImagePlaneWidget.h>
#include <vtkCellArray.h>
#include <vtkDataArray.h>
#include <vtkCallbackCommand.h>
#include <vtkCellPicker.h>
#include <vtkMapper.h>
#include <vtkPlaneSource.h>
#include <vtkSphereSource.h>
#include <vtkTextMapper.h>
#include <vtkTextProperty.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkSTLReader.h>
#include <vtkAutoInit.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <vtkFillHolesFilter.h>
#include <vtkImplicitBoolean.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkMaskPoints.h>
#include <vtkProperty.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkKdTreePointLocator.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricSpline.h>
#include <vtkCutter.h>
#include <vtkTubeFilter.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkLookupTable.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkCurvatures.h>
#include <vtkColorSeries.h>
//#include "torch/script.h"
//#include "torch/torch.h"
//#include <Eigen/Dense>
//#include <unsupported/Eigen/CXX11/Tensor>
#include<algorithm>
#include<vtkCenterOfMass.h>
#include<vtkCellCenters.h>
#include<vtkTextActor.h>
#include <vtkLine.h>
#include <vtkPolyDataSilhouette.h>
#include<vtkRegularPolygonSource.h>
#include<vtkPlane.h>
#include<vtkArrowSource.h>
#include <vtkCleanPolyData.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkScalarBarActor.h>
#include<vtkPolyLine.h>
#include <vtkBox.h>
#include <vtkExtractPolyDataGeometry.h>

#include <vtkPNGReader.h>

#include <vtkCubeSource.h>
#include <vtkSelectEnclosedPoints.h>

#include <mutex>
#include <condition_variable>
#include<cstring>
#include<tchar.h> 

typedef CGAL::Simple_cartesian<double>			Kernel;
typedef Kernel::Point_2                         Point_2;
typedef Kernel::Point_3                         Point_3;
typedef Kernel::Plane_3											Plane_3;
typedef Kernel::Segment_2										Segment_2;
typedef Kernel::Segment_3										Segment_3;
typedef Kernel::Vector_3										Vector_3;
typedef CGAL::Surface_mesh<Point_3>     SurfaceMesh;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor       face_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor       edge_descriptor;
typedef CGAL::Polygon_2<Kernel>             Polygon_2;

namespace SMP = CGAL::Surface_mesh_parameterization;
namespace PMP = CGAL::Polygon_mesh_processing;

struct MeshPoint
{
	// modified
	MeshPoint() : nTriId(0), uv(), xyz{ 0.0, 0.0, 0.0 }, he() {}
	MeshPoint(unsigned int a, Point_2 b, double c[3], halfedge_descriptor d
		, double dir[3] = NULL)
		:nTriId(a), uv(b), he(d)
	{
		xyz[0] = c[0]; xyz[1] = c[1]; xyz[2] = c[2];
	}
	unsigned int nTriId;
	Point_2 uv;
	double xyz[3];
	halfedge_descriptor he;
	//unsigned int TriId;
};

struct MeshParameterization
{
	SurfaceMesh* pMesh;
	std::map<unsigned int, face_descriptor> FaceMap;
	std::map<unsigned int, vertex_descriptor> VertexMap;
	std::map<unsigned int, edge_descriptor> EdgeMap;
	SurfaceMesh::Property_map<vertex_descriptor, Point_2> UVmap;
};

enum State { ADD = 0, MOVE, REMOVE, ADDF };

class SplineBase
{
public:
	void BaseSpline()
	{
		bIsSelected = false;
		nMoveIndex = -1;
		bRightPressSelectedSpline = false;
		nRightPressSelectedCtrlPoint = -1;
		m_bRebuildDeleteFlag = false;
	}
	//CGAL Mesh
	SurfaceMesh* pMesh;                                                                               
	std::map<unsigned int, face_descriptor> FaceMap;
	std::map<unsigned int, vertex_descriptor> VertexMap;
	std::map<unsigned int, edge_descriptor> EdgeMap;
	std::map<unsigned int, halfedge_descriptor> HalfedgeMap;
	SurfaceMesh::Property_map<vertex_descriptor, Point_2> UVmap;

	vtkSmartPointer<vtkPolyData> SplinePolydata; /**< Polydata for the spline. */

	vtkSmartPointer<vtkActor> SplineActor; /**< Actor for the spline. */

	std::vector<vtkSmartPointer<vtkSphereSource>> CtrlPointSphere; /**< Sphere sources representing control points of the spline. */

	std::vector<vtkSmartPointer<vtkActor>> CtrlPointActor; /**< Actors for the control point spheres. */

	std::vector<MeshPoint> vtSplinePoints; /**< All points along the spline. */

	std::vector<MeshPoint> vtCtrlPoints; /**< Control points of the spline. */

	std::vector<int> vtCtrlSubscript, vtEquiSubscript; /**< Indexes of control points and equally spaced points in the global coordinate system. */

	State nState; /**< Editing state of the spline. */

	bool bIsSelected; /**< Indicates whether the spline is selected. */

	int nMoveIndex; /**< Index of the control point being moved. */

	bool bRightPressSelectedSpline; /**< Indicates if the spline is selected with right mouse button (used when right-clicking to add a point). */

	bool m_bRebuildDeleteFlag; /**< Flag to determine if reconstruction is needed upon deletion of points. */

	int nRightPressSelectedCtrlPoint; /**< Index of the control point selected with a right-click. */

	int nIndexType;

	bool bDrawOnBorderFlag;


	static bool GetBary(Point_3& v1, Point_3& v2, Point_3& v3, Point_3& vp, double& Bary1, double& Bary2, double& Bary3);

	halfedge_descriptor FindAllCrosses(halfedge_descriptor he, Point_2 uv1, Point_2 uv2, bool first);

	bool inTriangle(halfedge_descriptor& TriId, Point_2& vp);

	double* UV2XYZ(halfedge_descriptor& he, Point_2& uv);

	Point_2 GetPointUV(halfedge_descriptor& he, Point_3& pt);

	Point_2 EdgeCross(Point_2& p1, Point_2& p2, Point_2& p3, Point_2& p4);

	Point_3 PlaneCrossSegment(Plane_3& P, Segment_3& S);

	double3 GetHalfEdgeNormal(halfedge_descriptor& he);

	
	virtual void initial(SurfaceMesh* sm, std::map<unsigned int, face_descriptor>& fmap, std::map<unsigned int, vertex_descriptor>& vmap, std::map<unsigned int, edge_descriptor>& emap, std::map<unsigned int, halfedge_descriptor>&hemap, SurfaceMesh::Property_map<vertex_descriptor, Point_2>& uvmap) = 0;
	
	virtual int add(unsigned int TriangleID, double PickPosition[3], double cameraDir[3] = NULL) = 0;
		
	virtual int move(unsigned int CtrlPtIndex, unsigned int TriangleID, double PickPosition[3], double cameraDir[3] = NULL) = 0;
	
	virtual int ChangeType(int IndexType) = 0;
	
	virtual void UpdateSpline(std::vector<MeshPoint>& SplinePoints) = 0;
	
	virtual void UpdateCtrlPoint(unsigned int& nIndex, int TriangleID, double PickPosition[3], State state) = 0;
	
	virtual int addF(unsigned int& CtrlPtIndex, unsigned int TriangleID, double PickPosition[3], double cameraDir[3] = NULL) = 0;
	
	virtual int remove(unsigned int CtrlPtIndex) = 0;
};