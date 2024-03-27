#ifndef CLOSEDSPLINEDESIGNINTERACTORSTYLE_H
#define CLOSEDSPLINEDESIGNINTERACTORSTYLE_H
#include "Spline/ClosedMeshSpline.h"
#include "Spline/Spline.h"
#include "Triangulator2D.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
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

typedef Kernel::Point_2                         Point_2;

class ClosedSplineDesignInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static ClosedSplineDesignInteractorStyle* New()
	{
		return new ClosedSplineDesignInteractorStyle;
	}
	vtkTypeMacro(ClosedSplineDesignInteractorStyle, vtkInteractorStyleTrackballCamera);

	ClosedSplineDesignInteractorStyle() :
		CellPicker(vtkSmartPointer<vtkCellPicker>::New()),
		PolyData(vtkSmartPointer<vtkPolyData>::New()),
		Renderer(vtkSmartPointer<vtkRenderer>::New()),
		PolyDataActor(vtkSmartPointer<vtkActor>::New()),
		RenderWindow(vtkSmartPointer<vtkRenderWindow>::New()),
		nMoveIndex(0) {}

	virtual ~ClosedSplineDesignInteractorStyle() {}
	void Init();
	void RemoveActors();
	void setSurfaceMesh(SurfaceMesh *newSurfaceMesh);
	void setPolyData(vtkSmartPointer<vtkPolyData> newPolyData);
	void setRenderer(vtkSmartPointer<vtkRenderer> newRenderer);
	void setPolydataActor(vtkSmartPointer<vtkActor> newPolydataActor);
	void setRenderWindow(vtkSmartPointer<vtkRenderWindow> newRenderWindow);
	// vtkInteractorStyle interface
	virtual void OnLeftButtonDown()  override;
	virtual void OnLeftButtonUp() override;
	virtual void OnRightButtonDown() override;// 避免vtk的GrabFocus接口占用交互命令
	//virtual void OnRightButtonUp() override;
	virtual void OnMouseMove() override;
	//virtual void OnMouseWheelForward()override;
	//virtual void OnMouseWheelBackward()override;
	//重写Rotate 函数  将旋转中心放在actor的中心上
	//virtual void Rotate() override;
	void GetRayIntersection(vtkObject* caller, int& nTriId, double* pWorld, vtkSmartPointer<vtkActor>& actor)
	{
		vtkSmartPointer<vtkRenderWindowInteractor>vtkInter = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		vtkInter = vtkRenderWindowInteractor::SafeDownCast(caller);
		if (!vtkInter) return;
		int* pEvtPos = vtkInter->GetEventPosition();
		vtkInter->FindPokedRenderer(pEvtPos[0], pEvtPos[1]);
		vtkInter->SetPicker(CellPicker);
		vtkInter->GetPicker()->Pick(pEvtPos[0], pEvtPos[1], 0, Renderer.GetPointer());
		nTriId = CellPicker->GetCellId();
		CellPicker->GetPickPosition(pWorld);
		actor = CellPicker->GetActor();
	}
	vtkObject* caller;
	ClosedMeshSpline* spline;
	std::vector<MeshPoint> GetSpline();

public:
	std::map<unsigned int, face_descriptor> fmap;
	std::map<unsigned int, vertex_descriptor> vmap;
	std::map<unsigned int, edge_descriptor> emap;
	std::map<unsigned int, halfedge_descriptor> hemap;
	SurfaceMesh *sm;
	vtkSmartPointer<vtkCellPicker> CellPicker;

	vtkSmartPointer<vtkPolyData> PolyData;
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkActor> PolyDataActor;
	vtkSmartPointer<vtkRenderWindow> RenderWindow;
	bool control_point_selected_mouse_move_flag;
	int nMoveIndex;

};

#endif // CLOSEDSPLINEDESIGNINTERACTORSTYLE_H