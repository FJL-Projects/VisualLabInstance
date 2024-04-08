#ifndef POLYGONMOVEMENTINTERACTOR_H
#define POLYGONMOVEMENTINTERACTOR_H
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <vtkTubeFilter.h>
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
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCallbackCommand.h>
#include <vtkCellPicker.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>
#include <vtkPolyLine.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkMath.h>
#include <vtkCellData.h>
#include <vtkAutoInit.h>
#include <vtkSTLWriter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/basic.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/intersections.h>
#include "Rotation.h"
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <vtkArrowSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include "vectorAlgorithm.h"

using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Point_2 = Kernel::Point_2;
using Plane_3 = Kernel::Plane_3;
using Vector_3 = Kernel::Vector_3;
using SurfaceMesh = CGAL::Surface_mesh<Point_3>;
using vertex_descriptor = boost::graph_traits<SurfaceMesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<SurfaceMesh>::halfedge_descriptor;
using face_descriptor = boost::graph_traits<SurfaceMesh>::face_descriptor;
enum class MoveState
{
	none, move, rotate, all_zoom, zoom_left, zoom_right, zoom_up, zoom_down, move_x, move_y, move_z
};

class PolygonMovementInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static PolygonMovementInteractorStyle* New()
	{
		return new PolygonMovementInteractorStyle;
	}
	vtkTypeMacro(PolygonMovementInteractorStyle, vtkInteractorStyleTrackballCamera);

	PolygonMovementInteractorStyle(){}

	void SetPolyData(vtkSmartPointer<vtkPolyData>);
	void SetRenderer(vtkSmartPointer<vtkRenderer>);
	void SetRenderWindow(vtkSmartPointer<vtkRenderWindow>);

	virtual void OnLeftButtonDown()  override;
	virtual void OnLeftButtonUp() override;
	virtual void OnMouseMove() override;
	void SetCorrectedOcclusalDirection(double3);
	double3 GetCorrectedOcclusalDirection();

	vtkSmartPointer<vtkRenderWindowInteractor> GetRayIntersection(vtkObject* caller, int& nTriId, double* pWorld, vtkSmartPointer<vtkActor>& actor)
	{
		vtkSmartPointer<vtkRenderWindowInteractor> vtkInter = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		vtkInter = vtkRenderWindowInteractor::SafeDownCast(caller);
		if (!vtkInter) return vtkSmartPointer<vtkRenderWindowInteractor>();
		int* pEvtPos = vtkInter->GetEventPosition();
		vtkInter->FindPokedRenderer(pEvtPos[0], pEvtPos[1]);
		vtkInter->SetPicker(m_cell_picker);
		vtkInter->GetPicker()->Pick(pEvtPos[0], pEvtPos[1], 0, m_renderer.GetPointer());
		nTriId = m_cell_picker->GetCellId();
		m_cell_picker->GetPickPosition(pWorld);
		actor = m_cell_picker->GetActor();
		return vtkInter;
	}
	vtkSmartPointer<vtkPolyData> Mesh2PolyData(SurfaceMesh& pmesh)
	{
		vtkSmartPointer<vtkPoints>  Points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> TrianglePolys = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkPolyData> Polydata = vtkSmartPointer<vtkPolyData>::New();
		std::map<vertex_descriptor, int> VertexDescriptorIndexMap;
		vtkIdType VertexIndex = 0;
		vtkIdType TriIndex = 0;

		for (vertex_descriptor v : CGAL::vertices(pmesh))
		{
			const boost::property_map_value<SurfaceMesh, CGAL::vertex_point_t>::type& p = get(CGAL::get(CGAL::vertex_point, pmesh), v);
			Points->InsertNextPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
			VertexDescriptorIndexMap[v] = VertexIndex++;
		}

		for (face_descriptor f : CGAL::faces(pmesh))
		{
			TriIndex = 0;
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
			for (halfedge_descriptor h : CGAL::halfedges_around_face(CGAL::halfedge(f, pmesh), pmesh))
				triangle->GetPointIds()->SetId(TriIndex++, VertexDescriptorIndexMap[CGAL::target(h, pmesh)]);
			TrianglePolys->InsertNextCell(triangle);
		}
		Polydata->SetPoints(Points);
		Polydata->SetPolys(TrianglePolys);
		return Polydata;
	}
	vtkObject* caller;

public:
	vtkSmartPointer<vtkCellPicker> m_cell_picker = vtkSmartPointer<vtkCellPicker>::New();
	vtkSmartPointer<vtkPolyData> m_polydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkRenderer> m_renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkActor> m_polydata_actor = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkRenderWindow> m_render_window = vtkSmartPointer<vtkRenderWindow>::New();
	MoveState m_state;
	int m_last_x = -1000;
	int m_last_y = -1000;
	std::array<double, 3> m_center_position;
	std::array<double, 3> m_last_pick_position;
	const std::array<double, 3> m_axis_x = { 1, 0, 0 };
	const std::array<double, 3> m_axis_y = { 0, 1, 0 };
	const std::array<double, 3> m_axis_z = { 0, 0, 1 };
	vtkSmartPointer<vtkTransform> m_transform_x;
	vtkSmartPointer<vtkTransform> m_transform_y;
	vtkSmartPointer<vtkTransform> m_transform_z;

	vtkSmartPointer<vtkActor> m_arrow_actor_x = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkActor> m_arrow_actor_y = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkActor> m_arrow_actor_z = vtkSmartPointer<vtkActor>::New();

	std::array<double, 3> m_corrected_occlusal_direction = { 0, 0, 0 };
};

#endif // POLYGONMOVEMENTINTERACTOR_H