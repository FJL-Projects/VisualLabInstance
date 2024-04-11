#include "PolygonMovementInteractor.h"

void PolygonMovementInteractorStyle::SetPolyDataAndRender(vtkSmartPointer<vtkPolyData> pd)
{
	m_polydata = pd;
    vtkSmartPointer<vtkPolyDataMapper> mesh_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mesh_mapper->SetInputData(m_polydata);
    mesh_mapper->Update();

    std::cout << m_color.g() << std::endl;

    m_polydata_actor = vtkSmartPointer<vtkActor>::New();
    m_polydata_actor->SetMapper(mesh_mapper);
    m_polydata_actor->GetProperty()->SetColor(m_color.r() / 255.0, m_color.g() / 255.0, m_color.b() / 255.0);
    m_polydata_actor->GetProperty()->SetAmbient(0.5);
    m_polydata_actor->GetProperty()->SetSpecularPower(100);
    m_polydata_actor->GetProperty()->SetSpecular(0.5);
    m_polydata_actor->GetProperty()->SetDiffuse(0.5);
    m_polydata_actor->GetProperty()->SetOpacity(m_opacity);
    m_renderer->AddActor(m_polydata_actor);
    m_renderer->GetActiveCamera()->SetFocalPoint(m_polydata->GetCenter());

    //m_transform_x = vtkSmartPointer<vtkTransform>::New();
    //m_transform_x->Translate(m_polydata->GetCenter());
    //m_transform_x->Scale(15.0, 15.0, 15.0);
    //m_transform_x->Update();

    //vtkSmartPointer<vtkArrowSource> arrow_source_x = vtkSmartPointer<vtkArrowSource>::New();
    //arrow_source_x->Update();

    //vtkSmartPointer<vtkTransformPolyDataFilter> transform_filter_x = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    //transform_filter_x->SetInputConnection(arrow_source_x->GetOutputPort());
    //transform_filter_x->SetTransform(m_transform_x);
    //transform_filter_x->Update();

    //vtkSmartPointer<vtkPolyDataMapper> arrow_mapper_x = vtkSmartPointer<vtkPolyDataMapper>::New();
    //arrow_mapper_x->SetInputConnection(transform_filter_x->GetOutputPort());
    //arrow_mapper_x->Update();

    //m_arrow_actor_x = vtkSmartPointer<vtkActor>::New();
    //m_arrow_actor_x->SetMapper(arrow_mapper_x);
    //m_arrow_actor_x->GetProperty()->SetColor(1, 0, 0);
    //m_arrow_actor_x->GetProperty()->SetAmbient(1.0);
    //m_arrow_actor_x->GetProperty()->SetSpecularPower(0);
    //m_arrow_actor_x->GetProperty()->SetSpecular(0);
    //m_arrow_actor_x->GetProperty()->SetDiffuse(0.5);
    //m_renderer->AddActor(m_arrow_actor_x);

    //m_transform_y = vtkSmartPointer<vtkTransform>::New();
    //m_transform_y->Translate(m_polydata->GetCenter());
    //m_transform_y->Scale(15.0, 15.0, 15.0);
    //m_transform_y->RotateZ(90.0);
    //m_transform_y->Update();

    //vtkSmartPointer<vtkArrowSource> arrow_source_y = vtkSmartPointer<vtkArrowSource>::New();
    //arrow_source_y->Update();

    //vtkSmartPointer<vtkTransformPolyDataFilter> transform_filter_y = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    //transform_filter_y->SetInputConnection(arrow_source_y->GetOutputPort());
    //transform_filter_y->SetTransform(m_transform_y);
    //transform_filter_y->Update();

    //vtkSmartPointer<vtkPolyDataMapper> arrow_mapper_y = vtkSmartPointer<vtkPolyDataMapper>::New();
    //arrow_mapper_y->SetInputConnection(transform_filter_y->GetOutputPort());
    //arrow_mapper_y->Update();

    //m_arrow_actor_y = vtkSmartPointer<vtkActor>::New();
    //m_arrow_actor_y->SetMapper(arrow_mapper_y);
    //m_arrow_actor_y->GetProperty()->SetColor(0, 1, 0);
    //m_arrow_actor_y->GetProperty()->SetAmbient(1.0);
    //m_arrow_actor_y->GetProperty()->SetSpecularPower(0);
    //m_arrow_actor_y->GetProperty()->SetSpecular(0);
    //m_arrow_actor_y->GetProperty()->SetDiffuse(0.5);
    //m_renderer->AddActor(m_arrow_actor_y);

    //m_transform_z = vtkSmartPointer<vtkTransform>::New();
    //m_transform_z->Translate(m_polydata->GetCenter());
    //m_transform_z->Scale(15.0, 15.0, 15.0);
    //m_transform_z->RotateY(90.0);
    //m_transform_z->Update();

    //vtkSmartPointer<vtkArrowSource> arrow_source_z = vtkSmartPointer<vtkArrowSource>::New();
    //arrow_source_z->Update();

    //vtkSmartPointer<vtkTransformPolyDataFilter> transform_filter_z = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    //transform_filter_z->SetInputConnection(arrow_source_z->GetOutputPort());
    //transform_filter_z->SetTransform(m_transform_z);
    //transform_filter_z->Update();

    //vtkSmartPointer<vtkPolyDataMapper> arrow_mapper_z = vtkSmartPointer<vtkPolyDataMapper>::New();
    //arrow_mapper_z->SetInputConnection(transform_filter_z->GetOutputPort());
    //arrow_mapper_z->Update();

    //m_arrow_actor_z = vtkSmartPointer<vtkActor>::New();
    //m_arrow_actor_z->SetMapper(arrow_mapper_z);
    //m_arrow_actor_z->GetProperty()->SetColor(0, 0, 1);
    //m_arrow_actor_z->GetProperty()->SetAmbient(1.0);
    //m_arrow_actor_z->GetProperty()->SetSpecularPower(0);
    //m_arrow_actor_z->GetProperty()->SetSpecular(0);
    //m_arrow_actor_z->GetProperty()->SetDiffuse(0.5);
    //m_renderer->AddActor(m_arrow_actor_z);

}

void PolygonMovementInteractorStyle::SetRenderer(vtkSmartPointer<vtkRenderer> r)
{
	m_renderer = r;
}

void PolygonMovementInteractorStyle::SetRenderWindow(vtkSmartPointer<vtkRenderWindow> rw)
{
	m_render_window = rw;
}

void PolygonMovementInteractorStyle::SetSurfaceMeshAndRender(SurfaceMesh& mesh)
{
    m_mesh = std::make_shared<SurfaceMesh>(mesh);
    SetPolyDataAndRender(Mesh2PolyData(*m_mesh));
    m_border_points.clear();
    std::vector<Point_3> border_points;

    for (auto& v : m_mesh->vertices())
	{
		if (m_mesh->is_border(v))
		{
			m_border_points.insert(static_cast<vtkIdType>(v.idx()));
            border_points.push_back(m_mesh->point(v));
		}
	}

    auto zoom_ratio = m_mesh->add_property_map<vertex_descriptor, double>("v:zoom_ratio", 0.0).first;
    for (auto& v : m_mesh->vertices())
    {
		double min_distance = std::numeric_limits<double>::max();
        for (auto& p : border_points)
        {
			double distance = CGAL::squared_distance(m_mesh->point(v), p);
            if (distance < min_distance)
            {
				min_distance = distance;
			}
		}

        double min_threshold = 0.3;
        double max_threshold = 0.6;

        if (min_distance < min_threshold)
        {
            zoom_ratio[v] = 0.0;
        }
        else if (min_distance < max_threshold)
        {
            zoom_ratio[v] = (min_distance - min_threshold) / (max_threshold - min_threshold);
        }
        else
        {
            zoom_ratio[v] = 1.0;
        }
    }
}

void PolygonMovementInteractorStyle::SetBlockMove(bool block)
{
	m_block_move = block;
}

void PolygonMovementInteractorStyle::SetBlockRotate(bool block)
{
    m_block_rotate = block;
}

void PolygonMovementInteractorStyle::SetColor(CGAL::Color color)
{
    m_color = color;
}

void PolygonMovementInteractorStyle::SetOpacity(double opacity)
{
    m_opacity = opacity;
}

SurfaceMesh PolygonMovementInteractorStyle::GetSurfaceMesh()
{
    for (auto& v : m_mesh->vertices())
    {
        double3 polydata_point(m_polydata->GetPoint(static_cast<vtkIdType>(v.idx())));
        m_mesh->point(v) = Point_3(polydata_point.x(), polydata_point.y(), polydata_point.z());
    }

    return *m_mesh;
}

void PolygonMovementInteractorStyle::OnLeftButtonDown()
{
    int tri_id = -1;
    double pick_position[3];
    vtkSmartPointer<vtkActor> pick_actor;
    auto vtk_inter = GetRayIntersection(caller, tri_id, pick_position, pick_actor);
    if (!vtk_inter->GetControlKey() && !vtk_inter->GetShiftKey() && pick_actor == m_arrow_actor_x)
    {
        m_state = MoveState::move_x;
    }
    else if (!vtk_inter->GetControlKey() && !vtk_inter->GetShiftKey() && pick_actor == m_arrow_actor_y)
    {
        m_state = MoveState::move_y;
    }
    else  if (!vtk_inter->GetControlKey() && !vtk_inter->GetShiftKey() && pick_actor == m_arrow_actor_z)
    {
        m_state = MoveState::move_z;
    }
    else
    {
        m_renderer->RemoveActor(m_arrow_actor_x);
        m_renderer->RemoveActor(m_arrow_actor_y);
        m_renderer->RemoveActor(m_arrow_actor_z);
    }
    if (vtk_inter->GetControlKey() && !vtk_inter->GetShiftKey())
    {
        m_state = MoveState::rotate;
    }
    if (!vtk_inter->GetControlKey() && vtk_inter->GetShiftKey())
    {
        int event_pos[2];
        vtk_inter->GetEventPosition(event_pos);
        int size[2];
        vtk_inter->GetSize(size);
        event_pos[0] -= size[0] / 2;
        event_pos[1] -= size[1] / 2;
        if (abs(event_pos[1]) > abs(event_pos[0]))
        {
            if (event_pos[1] > 0)
            {
                m_state = MoveState::zoom_up;
            }
            else
            {
                m_state = MoveState::zoom_down;
            }
        }
        else
        {
            if (event_pos[0] > 0)
            {
                m_state = MoveState::zoom_right;
            }
            else
            {
                m_state = MoveState::zoom_left;
            }
        }

    }
    if (vtk_inter->GetControlKey() && vtk_inter->GetShiftKey())
    {
        m_state = MoveState::all_zoom;
    }
    if (!vtk_inter->GetControlKey() && !vtk_inter->GetShiftKey() && pick_actor == m_polydata)
    {
        m_state = MoveState::move;
    }

    m_polydata->GetCenter(m_center_position.data());
}

void PolygonMovementInteractorStyle::OnLeftButtonUp()
{
    /*if (state > MoveState::all_zoom)
    {
        m_renderer->RemoveActor(m_polydata);
        auto sm = VTK_PolyData2CGAL_Surface_Mesh(m_polydata);
        CGAL::Polygon_mesh_processing::isotropic_remeshing(sm.faces(), 0.2, sm);
        m_polydata = Mesh2PolyData(sm);
        vtkSmartPointer<vtkPolyDataMapper> mesh_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mesh_mapper->SetInputData(m_polydata);
        mesh_mapper->Update();

        m_polydata->SetMapper(mesh_mapper);
        m_polydata->GetProperty()->SetColor(1, 0, 1);
        m_polydata->GetProperty()->SetAmbient(0.5);
        m_polydata->GetProperty()->SetSpecularPower(100);
        m_polydata->GetProperty()->SetSpecular(0.5);
        m_polydata->GetProperty()->SetDiffuse(0.5);
        m_polydata->GetProperty()->SetOpacity(1);
        m_renderer->AddActor(m_polydata);
        m_render_window->Render();
    }*/
    m_state = MoveState::none;
    m_last_x = -1000;
    m_last_y = -1000;
    m_last_pick_position[0] = 0;
    m_last_pick_position[1] = 0;
    m_last_pick_position[2] = 0;
    m_renderer->AddActor(m_arrow_actor_x);
    m_renderer->AddActor(m_arrow_actor_y);
    m_renderer->AddActor(m_arrow_actor_z);
    //m_renderer->GetActiveCamera()->SetFocalPoint(m_polydata->GetCenter());
   
}

void PolygonMovementInteractorStyle::OnMouseMove()
{
    int tri_id = -1;
    double pick_position[3];
    vtkSmartPointer<vtkActor> pick_actor;
    auto vtk_inter = GetRayIntersection(caller, tri_id, pick_position, pick_actor);
    //if (pick_actor != m_polydata) return;

    int event_pos[2] = { vtk_inter->GetEventPosition()[0], vtk_inter->GetEventPosition()[1] };
    vtk_inter->FindPokedRenderer(event_pos[0], event_pos[1]);
    double camera_position[3];
    m_renderer->GetActiveCamera()->GetPosition(camera_position);
    Plane_3 parallel_plane(Point_3(m_center_position[0], m_center_position[1], m_center_position[2]), Vector_3(Point_3(camera_position[0], camera_position[1], camera_position[2]) - Point_3(m_center_position[0], m_center_position[1], m_center_position[2])));

    auto intersection_point = boost::get<Point_3>(*CGAL::intersection(parallel_plane, Kernel::Ray_3(Point_3(camera_position[0], camera_position[1], camera_position[2]), Point_3(pick_position[0], pick_position[1], pick_position[2]))));
    pick_position[0] = intersection_point.x();
    pick_position[1] = intersection_point.y();
    pick_position[2] = intersection_point.z();

    //double m_center_position[3];
    //m_polydata->GetCenter(m_center_position);

    if (m_last_x != -1000 && m_last_y != -1000)
    {
        int delta_x = event_pos[0] - m_last_x;
        int delta_y = event_pos[1] - m_last_y;
        if (m_state == MoveState::move)
        {
            auto p1 = parallel_plane.projection(Point_3(pick_position[0], pick_position[1], pick_position[2]));
            auto p0 = parallel_plane.projection(Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]));
            auto v = p1 - p0;
            std::array<double, 3> p;
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p.data());
                p[0] += v[0];
                p[1] += v[1];
                p[2] += v[2];
                m_polydata->GetPoints()->SetPoint(i, p.data());
            }
            m_polydata->GetPoints()->Modified();

            m_transform_x->Identity();
            m_transform_x->Translate(m_polydata->GetCenter());
            m_transform_x->Scale(15.0, 15.0, 15.0);
            m_transform_x->Update();

            m_transform_y->Identity();
            m_transform_y->Translate(m_polydata->GetCenter());
            m_transform_y->Scale(15.0, 15.0, 15.0);
            m_transform_y->RotateZ(90.0);
            m_transform_y->Update();

            m_transform_z->Identity();
            m_transform_z->Translate(m_polydata->GetCenter());
            m_transform_z->Scale(15.0, 15.0, 15.0);
            m_transform_z->RotateY(90.0);
            m_transform_z->Update();

        }
        else if (m_state == MoveState::zoom_up)
        {
            double view_up[3];
            m_renderer->GetActiveCamera()->GetViewUp(view_up);
            Vector_3 up(view_up[0], view_up[1], view_up[2]);
            auto last_pick2pick = Point_3(pick_position[0], pick_position[1], pick_position[2]) - Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]);
            double p[3];
            double max_up = -1;
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p);
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                Vector_3 center2p(p[0], p[1], p[2]);
                double angle = acos(center2p * up / sqrt(center2p.squared_length()) / sqrt(up.squared_length()));
                double dis = sqrt(center2p.squared_length()) * cos(angle);
                if (dis > max_up)
                {
                    max_up = dis;
                }
            }

            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                vertex_descriptor v = static_cast<vertex_descriptor>(i);
                double ratio = m_mesh->property_map<vertex_descriptor, double>("v:zoom_ratio").first[v];
                m_polydata->GetPoint(i, p);
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                Vector_3 center2p(p[0], p[1], p[2]);
                double angle = acos(center2p * up / sqrt(center2p.squared_length()) / sqrt(up.squared_length()));
                double dis = sqrt(center2p.squared_length()) * cos(angle) / max_up;
                //double ratio=1.0 / (1.0 * sqrt(2 * vtkMath::Pi()))*exp(-(1-dis)*(1-dis)/2/1.0/1.0);
                if (center2p * up < 0)
                {
                    continue;
                }
                if (last_pick2pick * up > 0)
                {
                    p[0] += ratio * up[0] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[1] += ratio * up[1] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[2] += ratio * up[2] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                }
                else
                {
                    p[0] -= ratio * up[0] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[1] -= ratio * up[1] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[2] -= ratio * up[2] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                }
                p[0] += m_center_position[0];
                p[1] += m_center_position[1];
                p[2] += m_center_position[2];
                m_polydata->GetPoints()->SetPoint(i, p);
            }
            m_polydata->GetPoints()->Modified();
        }
        else if (m_state == MoveState::zoom_down)
        {
            double view_up[3];
            m_renderer->GetActiveCamera()->GetViewUp(view_up);
            Vector_3 up(view_up[0], view_up[1], view_up[2]);
            auto last_pick2pick = Point_3(pick_position[0], pick_position[1], pick_position[2]) - Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]);
            std::array<double, 3> p;
            double max_up = -1;
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p.data());
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                Vector_3 center2p(p[0], p[1], p[2]);
                double angle = acos(center2p * (-up) / sqrt(center2p.squared_length()) / sqrt(up.squared_length()));
                double dis = sqrt(center2p.squared_length()) * cos(angle);
                if (dis > max_up)
                {
                    max_up = dis;
                }
            }
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                vertex_descriptor v = static_cast<vertex_descriptor>(i);
                double ratio = m_mesh->property_map<vertex_descriptor, double>("v:zoom_ratio").first[v];
                m_polydata->GetPoint(i, p.data());
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                Vector_3 center2p(p[0], p[1], p[2]);
                double angle = acos(center2p * up / sqrt(center2p.squared_length()) / sqrt(up.squared_length()));
                double dis = sqrt(center2p.squared_length()) * cos(angle) / max_up;
                if (center2p * up > 0)
                {
                    continue;
                }
                if (last_pick2pick * up < 0)
                {
                    p[0] += ratio * up[0] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[1] += ratio * up[1] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[2] += ratio * up[2] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                }
                else
                {
                    p[0] -= ratio * up[0] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[1] -= ratio * up[1] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[2] -= ratio * up[2] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                }
                p[0] += m_center_position[0];
                p[1] += m_center_position[1];
                p[2] += m_center_position[2];
                m_polydata->GetPoints()->SetPoint(i, p.data());
            }
            m_polydata->GetPoints()->Modified();
        }
        else if (m_state == MoveState::zoom_right)
        {
            double view_up[3];
            m_renderer->GetActiveCamera()->GetViewUp(view_up);
            Vector_3 up(view_up[0], view_up[1], view_up[2]);
            auto last_pick2pick = Point_3(pick_position[0], pick_position[1], pick_position[2]) - Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]);
            Vector_3 center2camera = Point_3(m_center_position[0], m_center_position[1], m_center_position[2]) - Point_3(camera_position[0], camera_position[1], camera_position[2]);
            Vector_3 left = CGAL::cross_product(up, center2camera);
            left /= sqrt(left.squared_length());
           
            double p[3];
            double max_right = -1;
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p);
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                Vector_3 center2p(p[0], p[1], p[2]);
                double angle = acos(center2p * (left) / sqrt(center2p.squared_length()) / sqrt(left.squared_length()));
                double dis = sqrt(center2p.squared_length()) * cos(angle);
                if (dis > max_right)
                {
                    max_right = dis;
                }
            }
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                vertex_descriptor v = static_cast<vertex_descriptor>(i);
                double ratio = m_mesh->property_map<vertex_descriptor, double>("v:zoom_ratio").first[v];
                m_polydata->GetPoint(i, p);
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                Vector_3 center2p(p[0], p[1], p[2]);
                double angle = acos(center2p * left / sqrt(center2p.squared_length()) / sqrt(left.squared_length()));
                double dis = sqrt(center2p.squared_length()) * cos(angle) / max_right;
                if (center2p * left > 0)
                {
                    continue;
                }
                if (last_pick2pick * left < 0)
                {
                    p[0] += ratio * left[0] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[1] += ratio * left[1] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[2] += ratio * left[2] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                }
                else
                {
                    p[0] -= ratio * left[0] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[1] -= ratio * left[1] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[2] -= ratio * left[2] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                }        
                p[0] += m_center_position[0];
                p[1] += m_center_position[1];
                p[2] += m_center_position[2];
                m_polydata->GetPoints()->SetPoint(i, p);
            }
            m_polydata->GetPoints()->Modified();
        }
        else if (m_state == MoveState::zoom_left)
        {
            double view_up[3];
            m_renderer->GetActiveCamera()->GetViewUp(view_up);
            Vector_3 up(view_up[0], view_up[1], view_up[2]);
            auto last_pick2pick = Point_3(pick_position[0], pick_position[1], pick_position[2]) - Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]);
            Vector_3 center2camera = Point_3(m_center_position[0], m_center_position[1], m_center_position[2]) - Point_3(camera_position[0], camera_position[1], camera_position[2]);
            Vector_3 left = CGAL::cross_product(up, center2camera);
            left /= sqrt(left.squared_length());
            std::array<double, 3> p;
            double max_left = -1;
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p.data());
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                Vector_3 center2p(p[0], p[1], p[2]);
                double angle = acos(center2p * (left) / sqrt(center2p.squared_length()) / sqrt(left.squared_length()));
                double dis = sqrt(center2p.squared_length()) * cos(angle);
                if (dis > max_left)
                {
                    max_left = dis;
                }
            }
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                vertex_descriptor v = static_cast<vertex_descriptor>(i);
                double ratio = m_mesh->property_map<vertex_descriptor, double>("v:zoom_ratio").first[v];
                m_polydata->GetPoint(i, p.data());
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                Vector_3 center2p(p[0], p[1], p[2]);
                double angle = acos(center2p * left / sqrt(center2p.squared_length()) / sqrt(left.squared_length()));
                double dis = sqrt(center2p.squared_length()) * cos(angle) / max_left;
                if (center2p * left < 0)
                {
                    continue;
                }
                if (last_pick2pick * left > 0)
                {
                    p[0] += ratio * left[0] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[1] += ratio * left[1] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[2] += ratio * left[2] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                }
                else
                {
                    p[0] -= ratio * left[0] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[1] -= ratio * left[1] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                    p[2] -= ratio * left[2] * sqrt(last_pick2pick.squared_length()) * dis / 10;
                }
                p[0] += m_center_position[0];
                p[1] += m_center_position[1];
                p[2] += m_center_position[2];
                m_polydata->GetPoints()->SetPoint(i, p.data());
            }
            m_polydata->GetPoints()->Modified();
        }
        else if (m_state == MoveState::all_zoom)
        {
            auto last_pick2pick = Point_3(pick_position[0], pick_position[1], pick_position[2]) - Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]);
            int sign=event_pos[0] - m_last_x + event_pos[1] - m_last_y;
            std::array<double, 3> p;
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p.data());
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                if (sign > 0)
                {
                    p[0] += (sqrt(last_pick2pick.squared_length()) / 100);
                    p[1] += (sqrt(last_pick2pick.squared_length()) / 100);
                    p[2] += (sqrt(last_pick2pick.squared_length()) / 100);
                }
                if (sign < 0)
                {
                    p[0] -= (sqrt(last_pick2pick.squared_length()) / 100);
                    p[1] -= (sqrt(last_pick2pick.squared_length()) / 100);
                    p[2] -= (sqrt(last_pick2pick.squared_length()) / 100);
                }
                p[0] += m_center_position[0];
                p[1] += m_center_position[1];
                p[2] += m_center_position[2];
                m_polydata->GetPoints()->SetPoint(i, p.data());
            }
            m_polydata->GetPoints()->Modified();
        }
        else if (m_state == MoveState::rotate && !m_block_rotate)
        {
            auto center2camera = Point_3(m_center_position[0], m_center_position[1], m_center_position[2]) - Point_3(camera_position[0], camera_position[1], camera_position[2]);
            auto center2pick = Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]) - Point_3(pick_position[0], pick_position[1], pick_position[2]);
            auto axis = CGAL::cross_product(center2camera, center2pick);
            axis = axis / sqrt(axis.squared_length());
            auto last_pick2pick = Point_3(pick_position[0], pick_position[1], pick_position[2]) - Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]);
            if (last_pick2pick.squared_length() < 0.0001) 
            { 
                return; 
            }
            axis = axis * sqrt(last_pick2pick.squared_length()) / 10;
            std::array<double, 3> p;
            std::array<double, 3> angle_axis = { axis[0], axis[1] , axis[2] };
            AngleAxisRotatePoint(angle_axis, m_corrected_occlusal_direction, m_corrected_occlusal_direction);
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p.data());
                p[0] -= m_center_position[0];
                p[1] -= m_center_position[1];
                p[2] -= m_center_position[2];
                AngleAxisRotatePoint(angle_axis, p, p);
                p[0] += m_center_position[0];
                p[1] += m_center_position[1];
                p[2] += m_center_position[2];
                m_polydata->GetPoints()->SetPoint(i, p.data());
            }
            m_polydata->GetPoints()->Modified();
        }
        else if (m_state == MoveState::move_x && !m_block_move)
        {
            Kernel::Line_3 line(Point_3(m_center_position[0], m_center_position[1], m_center_position[2]), Kernel::Vector_3(m_axis_x[0], m_axis_x[1], m_axis_x[2]));
            auto p1 = line.projection(Point_3(pick_position[0], pick_position[1], pick_position[2]));
            auto p0 = line.projection(Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]));
            auto v = p1 - p0;
            double p[3];
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p);
                p[0] += v[0];
                p[1] += v[1];
                p[2] += v[2];
                m_polydata->GetPoints()->SetPoint(i, p);
            }
            m_polydata->GetPoints()->Modified();

            m_transform_x->Identity();
            m_transform_x->Translate(m_polydata->GetCenter());
            m_transform_x->Scale(15.0, 15.0, 15.0);
            m_transform_x->Update();

            m_transform_y->Identity();
            m_transform_y->Translate(m_polydata->GetCenter());
            m_transform_y->Scale(15.0, 15.0, 15.0);
            m_transform_y->RotateZ(90.0);
            m_transform_y->Update();

            m_transform_z->Identity();
            m_transform_z->Translate(m_polydata->GetCenter());
            m_transform_z->Scale(15.0, 15.0, 15.0);
            m_transform_z->RotateY(90.0);
            m_transform_z->Update();
        }
        else if (m_state == MoveState::move_y && !m_block_move)
        {
            Kernel::Line_3 line(Point_3(m_center_position[0], m_center_position[1], m_center_position[2]), Kernel::Vector_3(m_axis_y[0], m_axis_y[1], m_axis_y[2]));
            auto p1 = line.projection(Point_3(pick_position[0], pick_position[1], pick_position[2]));
            auto p0 = line.projection(Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]));
            auto v = p1 - p0;
            double p[3];
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p);
                p[0] += v[0];
                p[1] += v[1];
                p[2] += v[2];
                m_polydata->GetPoints()->SetPoint(i, p);
            }
            m_polydata->GetPoints()->Modified();

            m_transform_x->Identity();
            m_transform_x->Translate(m_polydata->GetCenter());
            m_transform_x->Scale(15.0, 15.0, 15.0);
            m_transform_x->Update();

            m_transform_y->Identity();
            m_transform_y->Translate(m_polydata->GetCenter());
            m_transform_y->Scale(15.0, 15.0, 15.0);
            m_transform_y->RotateZ(90.0);
            m_transform_y->Update();

            m_transform_z->Identity();
            m_transform_z->Translate(m_polydata->GetCenter());
            m_transform_z->Scale(15.0, 15.0, 15.0);
            m_transform_z->RotateY(90.0);
            m_transform_z->Update();
        }
        else if (m_state == MoveState::move_z && !m_block_move)
        {
            Kernel::Line_3 line(Point_3(m_center_position[0], m_center_position[1], m_center_position[2]), Kernel::Vector_3(m_axis_z[0], m_axis_z[1], m_axis_z[2]));
            auto p1 = line.projection(Point_3(pick_position[0], pick_position[1], pick_position[2]));
            auto p0 = line.projection(Point_3(m_last_pick_position[0], m_last_pick_position[1], m_last_pick_position[2]));
            auto v = p1 - p0;
            double p[3];
            for (int i = 0; i < m_polydata->GetNumberOfPoints(); i++)
            {
                m_polydata->GetPoint(i, p);
                p[0] += v[0];
                p[1] += v[1];
                p[2] += v[2];
                m_polydata->GetPoints()->SetPoint(i, p);
            }
            m_polydata->GetPoints()->Modified();

            m_transform_x->Identity();
            m_transform_x->Translate(m_polydata->GetCenter());
            m_transform_x->Scale(15.0, 15.0, 15.0);
            m_transform_x->Update();

            m_transform_y->Identity();
            m_transform_y->Translate(m_polydata->GetCenter());
            m_transform_y->Scale(15.0, 15.0, 15.0);
            m_transform_y->RotateZ(90.0);
            m_transform_y->Update();

            m_transform_z->Identity();
            m_transform_z->Translate(m_polydata->GetCenter());
            m_transform_z->Scale(15.0, 15.0, 15.0);
            m_transform_z->RotateY(90.0);
            m_transform_z->Update();
        }
    }
    m_last_x = event_pos[0];
    m_last_y = event_pos[1];
    m_last_pick_position[0] = pick_position[0];
    m_last_pick_position[1] = pick_position[1];
    m_last_pick_position[2] = pick_position[2];
    m_render_window->Render();
}

void PolygonMovementInteractorStyle::SetCorrectedOcclusalDirection(double3 val)
{
    for (size_t i = 0; i < 3; ++i)
    {
        m_corrected_occlusal_direction[i] = val[i];
    }
}

double3 PolygonMovementInteractorStyle::GetCorrectedOcclusalDirection()
{
    double3 direction(m_corrected_occlusal_direction.data());
    return direction;
}

vtkSmartPointer<vtkActor> PolygonMovementInteractorStyle::GetPolydataActor()
{
    return m_polydata_actor;
}

void PolygonMovementInteractorStyle::SetConstrainBorder(bool b)
{
    m_constrain_border = b;
}
