#include "ClosedSplineDesignInteractorStyle.h"

void ClosedSplineDesignInteractorStyle::setSurfaceMesh(SurfaceMesh* newSurfaceMesh)
{
	sm = newSurfaceMesh;
}

void ClosedSplineDesignInteractorStyle::setPolyData(vtkSmartPointer<vtkPolyData> newPolyData)
{
	PolyData = newPolyData;
}

void ClosedSplineDesignInteractorStyle::Init()
{
	halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(*sm).first;
	unsigned verIndex = 0;
	for (auto v : sm->vertices()) vmap[verIndex++] = v;
	unsigned faceIndex = 0;
	for (auto f : sm->faces()) fmap[faceIndex++] = f;
	unsigned halfedgeIndex = 0;
	for (auto he : sm->halfedges()) hemap[halfedgeIndex++] = he;
	SurfaceMesh::Property_map<vertex_descriptor, Point_2> uv_map = sm->add_property_map<vertex_descriptor, Point_2>("h:uv").first;
	CGAL::Surface_mesh_parameterization::parameterize(*sm, bhd, uv_map);

	control_point_selected_mouse_move_flag = false;

	spline = new ClosedMeshSpline;
	spline->initial(sm, fmap, vmap, emap, hemap, uv_map);
	vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
	tubeFilter->SetInputData(spline->SplinePolydata);
	tubeFilter->SetRadius(0.075);
	tubeFilter->SetNumberOfSides(16);
	tubeFilter->Update();
	vtkSmartPointer<vtkPolyDataNormals> vtkNormal1 = vtkSmartPointer<vtkPolyDataNormals>::New();
	vtkNormal1->SetInputConnection(tubeFilter->GetOutputPort());
	vtkNormal1->SetComputePointNormals(1);
	vtkNormal1->SetComputeCellNormals(1);
	vtkNormal1->SetAutoOrientNormals(1);
	vtkNormal1->SetSplitting(0);
	vtkNormal1->FlipNormalsOff();
	vtkNormal1->Update();
	vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper1->SetInputConnection(vtkNormal1->GetOutputPort());
	mapper1->Update();
	spline->SplineActor->SetMapper(mapper1);
	spline->SplineActor->GetProperty()->SetColor(1, 1, 0);
	spline->SplineActor->GetProperty()->SetAmbient(0.5);
	spline->SplineActor->GetProperty()->SetSpecularPower(100);
	spline->SplineActor->GetProperty()->SetSpecular(0.5);
	spline->SplineActor->GetProperty()->SetDiffuse(0.5);
	spline->SplineActor->PickableOff();

	Renderer->AddActor(spline->SplineActor);
}

void ClosedSplineDesignInteractorStyle::RemoveActors()
{
	if (spline->CtrlPointActor.size() > 0)
	{
		//std::cout << "Removed ctrl point actors" << std::endl;
		for (auto actor = spline->CtrlPointActor.begin(); actor != spline->CtrlPointActor.end(); ++actor)
		{
			if (*actor)
			{
				Renderer->RemoveActor(*actor);
			}
		}
		spline->CtrlPointActor.clear();
	}

	if (spline->SplineActor)
	{
		//std::cout << "Removed spline actor" << std::endl;
		Renderer->RemoveActor(spline->SplineActor);
	}
}

void ClosedSplineDesignInteractorStyle::setRenderer(vtkSmartPointer<vtkRenderer> newRenderer)
{
	Renderer = newRenderer;
}

void ClosedSplineDesignInteractorStyle::setPolydataActor(vtkSmartPointer<vtkActor> newPolydataActor)
{
	PolyDataActor = newPolydataActor;
}

void ClosedSplineDesignInteractorStyle::setRenderWindow(vtkSmartPointer<vtkRenderWindow> newRenderWindow)
{
	RenderWindow = newRenderWindow;
}
void ClosedSplineDesignInteractorStyle::OnLeftButtonDown()
{
	int TriID = -1;
	double Pick_Pos[3] = { 0, 0, 0 };
	vtkSmartPointer<vtkActor> Pick_actor;
	GetRayIntersection(caller, TriID, Pick_Pos, Pick_actor);

	if (TriID == -1) return;
	if (Pick_actor == PolyDataActor && spline->bClosed == false)
	{
		spline->add(TriID, Pick_Pos);
		vtkSmartPointer<vtkPolyDataNormals> vtkNormal = vtkSmartPointer<vtkPolyDataNormals>::New();
		vtkNormal->SetInputConnection(spline->CtrlPointSphere.back()->GetOutputPort());
		vtkNormal->SetComputePointNormals(1);
		vtkNormal->SetComputeCellNormals(1);
		vtkNormal->SetAutoOrientNormals(1);
		vtkNormal->SetSplitting(0);
		vtkNormal->FlipNormalsOff();
		vtkNormal->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(vtkNormal->GetOutputPort());
		mapper->Update();
		spline->CtrlPointActor.back()->SetMapper(mapper);
		spline->CtrlPointActor.back()->GetProperty()->SetColor(0, 0, 1);
		spline->CtrlPointActor.back()->GetProperty()->SetAmbient(0.5);
		spline->CtrlPointActor.back()->GetProperty()->SetSpecularPower(100);
		spline->CtrlPointActor.back()->GetProperty()->SetSpecular(0.5);
		spline->CtrlPointActor.back()->GetProperty()->SetDiffuse(0.5);
		spline->CtrlPointActor.back()->PickableOn();
		Renderer->AddActor(spline->CtrlPointActor.back());
	}
	if (Pick_actor == spline->CtrlPointActor.front())
	{
		spline->bClosed = true;
		spline->UpdateSpline(spline->vtEquidistantSpline);
	}
	if (spline->bClosed)
	{
		for (auto actor : spline->CtrlPointActor)
			if (Pick_actor == actor)
			{
				control_point_selected_mouse_move_flag = true;
				nMoveIndex = std::distance(spline->CtrlPointActor.begin(), std::find(spline->CtrlPointActor.begin(), spline->CtrlPointActor.end(), actor));
				actor->PickableOff();
				break;
			}
	}
	RenderWindow->Render();
}

void ClosedSplineDesignInteractorStyle::OnMouseMove()
{
	// To avoid the mouse move event being called without control point being selected.
	if (!control_point_selected_mouse_move_flag)
	{
		return;
	}

	int TriID = -1;
	double Pick_Pos[3] = { 0,0,0 };
	vtkSmartPointer<vtkActor>Pick_actor;
	GetRayIntersection(caller, TriID, Pick_Pos, Pick_actor);
	if (TriID == -1) return;
	if (spline->bClosed)
	{
		if (nMoveIndex != -1)
		{
			spline->move(nMoveIndex, TriID, Pick_Pos);
			RenderWindow->Render();
		}
	}
}

std::vector<MeshPoint> ClosedSplineDesignInteractorStyle::GetSpline()
{
	return spline->vtEquidistantSpline;
}

void ClosedSplineDesignInteractorStyle::OnLeftButtonUp()
{
	control_point_selected_mouse_move_flag = false;

	int TriID = -1;
	double Pick_Pos[3] = { 0,0,0 };
	vtkSmartPointer<vtkActor>Pick_actor;
	GetRayIntersection(caller, TriID, Pick_Pos, Pick_actor);
	if (TriID == -1) return;
	if (nMoveIndex != -1)
	{
		spline->CtrlPointActor[nMoveIndex]->PickableOn();
	}
	nMoveIndex = -1;
}

void ClosedSplineDesignInteractorStyle::OnRightButtonDown()
{
	spline->SplineActor->PickableOn();
	int TriID = -1;
	double Pick_Pos[3] = { 0,0,0 };
	vtkSmartPointer<vtkActor>Pick_actor;
	GetRayIntersection(caller, TriID, Pick_Pos, Pick_actor);
	spline->SplineActor->PickableOff();
	if (TriID == -1) return;
	if (Pick_actor == spline->SplineActor)
	{
		unsigned int nID = (TriID / 2) % ((spline->vtEquidistantSpline.size() - 1));
		spline->addF(nID, TriID, Pick_Pos);
		vtkSmartPointer<vtkPolyDataNormals> vtkNormal = vtkSmartPointer<vtkPolyDataNormals>::New();
		vtkNormal->SetInputConnection(spline->CtrlPointSphere[nID]->GetOutputPort());
		vtkNormal->SetComputePointNormals(1);
		vtkNormal->SetComputeCellNormals(1);
		vtkNormal->SetAutoOrientNormals(1);
		vtkNormal->SetSplitting(0);
		vtkNormal->FlipNormalsOff();
		vtkNormal->Update();
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(vtkNormal->GetOutputPort());
		mapper->Update();
		spline->CtrlPointActor[nID]->SetMapper(mapper);
		spline->CtrlPointActor[nID]->GetProperty()->SetColor(0, 0, 1);
		spline->CtrlPointActor[nID]->GetProperty()->SetAmbient(0.5);
		spline->CtrlPointActor[nID]->GetProperty()->SetSpecularPower(100);
		spline->CtrlPointActor[nID]->GetProperty()->SetSpecular(0.5);
		spline->CtrlPointActor[nID]->GetProperty()->SetDiffuse(0.5);
		spline->CtrlPointActor[nID]->PickableOn();
		Renderer->AddActor(spline->CtrlPointActor[nID]);
	}
	else if (spline->bClosed)
	{
		for (auto actor : spline->CtrlPointActor)
			if (Pick_actor == actor)
			{
				Renderer->RemoveActor(actor);
				spline->remove(std::distance(spline->CtrlPointActor.begin(), std::find(spline->CtrlPointActor.begin(), spline->CtrlPointActor.end(), actor)));

				break;
			}
	}
	RenderWindow->Render();
}