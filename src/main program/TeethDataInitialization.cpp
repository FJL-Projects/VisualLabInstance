#include "TeethDataInitialization.h"
#include <vtkPLYWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkKochanekSpline.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkSphereSource.h>
#include <vtkNamedColors.h>
#include <vtkGlyph3DMapper.h>

/**
 * This method performs fragment clearing for teeth data.
 * 
 * It identifies and clears fragments in the teeth surface mesh based on the provided labels.
 */
void TeethDataInitialization::FragmentClearing()
{
    vtkSmartPointer<vtkCellData> celldata = m_teethPolyData->GetCellData();
    vtkSmartPointer<vtkDataArray> labels = celldata->GetArray("Label");
    std::vector<int> fragment;
    std::vector<std::vector<int>> group(20, std::vector<int>());
    std::vector<char> visit(m_teethSurfaceMesh->number_of_faces(), false);

    for (int i = 0; i < 20; i++)
    {
        for (SurfaceMesh::Face_index f : m_teethSurfaceMesh->faces())
        {
            int label = static_cast<int>(labels->GetTuple(f.idx())[0]);
            if (label != i || visit[f.idx()])
            {
                continue;
            }

            std::queue<int> bfsQueue;
            bfsQueue.push(f.idx());
            std::vector<int> tmp;
            while (!bfsQueue.empty())
            {
                int currentFaceIdx = bfsQueue.front();
                bfsQueue.pop();
                if (visit[currentFaceIdx] == true)
                {
                    continue;
                }
                tmp.push_back(currentFaceIdx);
                visit[currentFaceIdx] = true;

                SurfaceMesh::Halfedge_index he0 = m_teethSurfaceMesh->halfedge(SurfaceMesh::Face_index(currentFaceIdx));
                SurfaceMesh::Halfedge_index he1 = m_teethSurfaceMesh->next(he0);
                SurfaceMesh::Halfedge_index he2 = m_teethSurfaceMesh->next(he1);
                if (!m_teethSurfaceMesh->is_border(m_teethSurfaceMesh->opposite(he0)))
                {
                    SurfaceMesh::Face_index oppositeFace0 = m_teethSurfaceMesh->face(m_teethSurfaceMesh->opposite(he0));
                    if (!visit[oppositeFace0.idx()] && static_cast<int>(labels->GetTuple(oppositeFace0.idx())[0]) == i)
                    {
                        bfsQueue.push(oppositeFace0.idx());
                    }
                }
                if (!m_teethSurfaceMesh->is_border(m_teethSurfaceMesh->opposite(he1)))
                {
                    SurfaceMesh::Face_index oppositeFace1 = m_teethSurfaceMesh->face(m_teethSurfaceMesh->opposite(he1));
                    if (!visit[oppositeFace1.idx()] && static_cast<int>(labels->GetTuple(oppositeFace1.idx())[0]) == i)
                    {
                        bfsQueue.push(oppositeFace1.idx());
                    }
                }
                if (!m_teethSurfaceMesh->is_border(m_teethSurfaceMesh->opposite(he2)))
                {
                    SurfaceMesh::Face_index oppositeFace2 = m_teethSurfaceMesh->face(m_teethSurfaceMesh->opposite(he2));
                    if (!visit[oppositeFace2.idx()] && static_cast<int>(labels->GetTuple(oppositeFace2.idx())[0]) == i)
                    {
                        bfsQueue.push(oppositeFace2.idx());
                    }
                }
            }


            if (tmp.size() > group[i].size())
            {
                for (int j = 0; j < group[i].size(); j++)
                {
                    fragment.push_back(group[i][j]);
                }
                group[i] = tmp;
                //std::cout << i << " " << tmp.size() << endl;
            }
            else
            {
                for (int j = 0; j < tmp.size(); j++)
                {
                    fragment.push_back(tmp[j]);
                }
            }
        }

    }
    //std::cout << fragment.size();
    for (auto f : fragment)
    {
        double* label = new double[1] {8};
        labels->SetTuple(f, label);
    }
    m_teethPolyData->GetCellData()->Modified();
    m_teethPolyData->Modified();
    //vtkNew<vtkXMLPolyDataWriter> writer;
    //writer->SetInputData(m_teethPolyData);
    //writer->SetFileName("vtpFileName.vtp");
    //writer->Write();
}

/**
 * This method retrieves the labels from the teeth polydata and performs the necessary operations.
 */
void TeethDataInitialization::GetLabel()
{
    vtkSmartPointer<vtkCellData> celldata = m_teethPolyData->GetCellData();
    vtkSmartPointer<vtkDataArray> labels = celldata->GetArray("Label");
    vtkNew<vtkTriangle> triangle;
    //std::cout << labels->GetSize() << endl;

    // Split the polydata based on the label types
    for (vtkIdType i = 0; i < labels->GetSize(); i++)
    {
        // tmp is a double array that contains an integer value (in this case, only one value)
        int label = static_cast<int>(labels->GetTuple(i)[0]);

        // polydata->GetCell(i) retrieves the triangle face, and polydata->GetCell(i)->GetPointId(0) retrieves the id of the 0th point
        for (int j = 0; j < 3; j++)
        {
            if (m_splitPolydataVisit[label][m_teethPolyData->GetCell(i)->GetPointId(j)] == -1)
            {
                auto face_i = m_teethPolyData->GetCell(i); // Ids of the points of the i-th triangle face, which are three int point ids
                auto point_id_j = face_i->GetPointId(j); // The id of the 0th point of the triangle face (output is the 0th value of the three int point ids)
                auto new_point_id = m_splitPolyDataFlag[label]; // The initial value of flag is 0, and it outputs the id value in the new container
                m_labeledPoints[label]->InsertPoint(new_point_id, m_teethPolyData->GetPoint(point_id_j));
                triangle->GetPointIds()->SetId(j, new_point_id); // Set the id values of the new triangle face
                m_splitPolydataVisit[label][point_id_j] = new_point_id;
                m_splitPolyDataFlag[label]++;
            }
            else
            {
                triangle->GetPointIds()->SetId(j, m_splitPolydataVisit[label][m_teethPolyData->GetCell(i)->GetPointId(j)]);
            }
        }

        m_labeledTriangles[label]->InsertNextCell(triangle);
    }

    // Set the labeled points and triangles for each label
    for (int i = 0; i < 20; i++)
    {
        //std::cout << m_splitPolyDataFlag[i] << endl;
        //std::cout << i << " " << m_labeledPoints[i]->GetNumberOfPoints() << " " << m_labeledTriangles[i]->GetNumberOfCells() << endl;
        m_labeledPolyData[i]->SetPoints(m_labeledPoints[i]);
        m_labeledPolyData[i]->SetPolys(m_labeledTriangles[i]);
    }

    //// Write the labeled polydata to PLY files
    //for (int i = 0; i < 20; i++)
    //{
    //    vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();

    //    std::string tmp = "data" + std::to_string(i) + ".ply";
    //    const char* filename = tmp.c_str();
    //    writer->SetFileTypeToASCII();
    //    writer->SetFileName(filename);
    //    writer->SetInputData(m_labeledPolyData[i]);
    //    writer->Update();
    //    writer->Write();
    //}
}

/**
 * This method creates labeled actors for each tooth based on the labeled polydata.
 * It computes the point normals for each labeled polydata and assigns colors and properties to the actors.111111ww
 */
void TeethDataInitialization::MakeLabeledActor()
{
    for (int i = 0; i < 20; i++)
    {
        //std::cout << "Making labeled actor: "  << i << std::endl;
        if (m_labeledPolyData[i]->GetNumberOfPolys() == 0) // If the number of polygons is zero, the tooth does not exist
        {
            continue;
        }
        //std::cout << i << " " << m_labeledPolyData[i]->GetNumberOfPoints() << endl;

        // Compute point normals for the labeled polydata
        vtkNew<vtkPolyDataNormals> vtkNormal;
        vtkNormal->SetInputData(m_labeledPolyData[i]);
        vtkNormal->SetComputePointNormals(1); // Enable point normal calculation
        vtkNormal->SetComputeCellNormals(0); // Disable cell normal calculation
        vtkNormal->SetAutoOrientNormals(1); // Enable automatic adjustment of normals
        vtkNormal->SetSplitting(0); // By default, vtkPolyDataNormals splits sharp edges. Turning it off to prevent changes in the data.
        vtkNormal->FlipNormalsOn(); // Flip cell normals
        vtkNormal->Update();

        // Create a mapper for the labeled polydata
        vtkNew<vtkPolyDataMapper> vtkMapper;
        vtkMapper->SetInputData(vtkNormal->GetOutput());

        // Create a labeled actor and set the mapper, color, and properties
        m_labeledActors[i]->SetMapper(vtkMapper);
        double x = 0.1 * i;
        if (i == 0)
        {
            m_labeledActors[i]->GetProperty()->SetColor(222 / 255.0, 184 / 255.0, 135 / 255.0);
        }
        else
        {
            m_labeledActors[i]->GetProperty()->SetColor(1 - x, x, x);
        }
        m_labeledActors[i]->GetProperty()->SetAmbient(0.3);
        m_labeledActors[i]->GetProperty()->SetDiffuse(0.5);
        m_labeledActors[i]->GetProperty()->SetSpecular(0.1);
        m_labeledActors[i]->GetProperty()->SetOpacity(1);

        // Add the labeled actor to the render
        //m_vtkRender->AddActor(LabeledActors[i]);
    }
}


/**
 * Computes the oriented bounding boxes (OBBs) for teeth data.
 * - Sets the number of points in the teeth center to 7.
 * - Iterates through teeth 1 to 7 and calculates the OBB for each tooth.
 * - Handles missing teeth by interpolating their center points based on adjacent teeth.
 *
 * The OBB computation involves the following steps:
 * 1. Creates an OBB tree for the labeled polydata of each tooth.
 * 2. Computes the OBB using the OBB tree.
 * 3. Retrieves the bounds of the tooth.
 * 4. Calculates the center point of the tooth using the bounds.
 *
 * For missing teeth:
 * - Interpolates the center point of the missing tooth based on the center points of adjacent teeth.
 * - Handles the first missing tooth (index 0) by extrapolating the center point based on the next two adjacent teeth.
 * - Handles the last missing tooth (index 6) by extrapolating the center point based on the previous two adjacent teeth.
 */
void TeethDataInitialization::ObbBoundingBox()
{
    m_teeth_center->SetNumberOfPoints(7); // Set the number of points in the teeth center to 7

    for (int i = 1; i < 8; i++) // Iterate through teeth 1 to 7
    {
        if (m_labeledPolyData[i]->GetNumberOfPolys() < 10) // If the number of triangles in the tooth is less than 10, it is considered missing
        {
            m_missingTeeth.push_back(i - 1); // Add the missing tooth index to the missingTeeth vector
            continue; // Continue to the next iteration
        }
        m_valid_teeth_indices.push_back(i);
        ++m_valid_teeth_num;
        int maxLevel = 100;

        // Create an OBB tree for the labeled polydata of the tooth
        vtkNew<vtkOBBTree> obbTree;
        obbTree->SetDataSet(m_labeledPolyData[i]);
        obbTree->SetMaxLevel(maxLevel);
        obbTree->BuildLocator();

        double corner[3] = { 0.0, 0.0, 0.0 };
        double max[3] = { 0.0, 0.0, 0.0 };
        double mid[3] = { 0.0, 0.0, 0.0 };
        double min[3] = { 0.0, 0.0, 0.0 };
        double size[3] = { 0.0, 0.0, 0.0 };

        // Compute the oriented bounding box (OBB) for the tooth
        obbTree->ComputeOBB(m_labeledPolyData[i], corner, max, mid, min, size);

        auto bounds = m_labeledPolyData[i]->GetBounds(); // Get the bounding box of the tooth

        // Compute the center point of the tooth using the bounds of the tooth
        m_teeth_center->SetPoint(i - 1, (bounds[1] - bounds[0]) / 2.0 + bounds[0], (bounds[3] - bounds[2]) / 2.0 + bounds[2], (bounds[5] - bounds[4]) / 2.0 + bounds[4]);
    }

    // Handle missing teeth
    for (int i = 0; i < m_missingTeeth.size(); i++)
    {
        double pre[3], next[3];

        // Get the center point of the previous tooth
        if (m_missingTeeth[i] != 0)
        {
            m_teeth_center->GetPoint(m_missingTeeth[i] - 1, pre);
        }

        // Get the center point of the next tooth
        if (m_missingTeeth[i] != 6)
        {
            m_teeth_center->GetPoint(m_missingTeeth[i] + 1, next);
        }

        // Compute the center point of the missing tooth by averaging the center points of the adjacent teeth
        if (m_missingTeeth[i] != 0 && m_missingTeeth[i] != 6)
        {
            m_teeth_center->SetPoint(m_missingTeeth[i], (pre[0] + next[0]) / 2.0, (pre[1] + next[1]) / 2.0, (pre[2] + next[2]) / 2.0);
        }

        // Handle the first missing tooth (index 0)
        if (m_missingTeeth[i] == 0)
        {
            double3 next1(m_teeth_center->GetPoint(m_missingTeeth[i] + 1));
            double3 next2(m_teeth_center->GetPoint(m_missingTeeth[i] + 2));
            double3 pt0 = next1 * 2 - next2;
            m_teeth_center->SetPoint(0, pt0.data); // Set the center point of the missing tooth using extrapolated values
        }

        // Handle the last missing tooth (index 6)
        if (m_missingTeeth[i] == 6)
        {
            double3 pre1(m_teeth_center->GetPoint(m_missingTeeth[i] - 1));
            double3 pre2(m_teeth_center->GetPoint(m_missingTeeth[i] - 2));
            double3 pt6 = pre1 * 2 - pre2;
            m_teeth_center->SetPoint(6, pt6.data); // Set the center point of the missing tooth using extrapolated values
        }
    }
}

/**
 * @brief This method initializes and draws a spline representing the teeth center.
 */
void TeethDataInitialization::DrawTeethCenterSpline()
{
    std::vector<std::vector<double>> all_spline_points;
    vtkNew<vtkNamedColors> color;

    // Create splines for x, y, and z coordinates
    vtkNew<vtkKochanekSpline> xSpline;
    vtkNew<vtkKochanekSpline> ySpline;
    vtkNew<vtkKochanekSpline> zSpline;

    vtkNew<vtkParametricSpline> spline;
    spline->SetXSpline(xSpline);
    spline->SetYSpline(ySpline);
    spline->SetZSpline(zSpline);
    spline->SetPoints(m_teeth_center);

    vtkNew<vtkParametricFunctionSource> functionSource;
    functionSource->SetParametricFunction(spline);
    functionSource->SetUResolution(200);
    functionSource->SetVResolution(200);
    functionSource->SetWResolution(200);
    functionSource->Update();

    auto splinepolydata = functionSource->GetOutput();
    vtkIdList* ptId = vtkIdList::New();

    // Get the points from the spline polydata
    for (vtkIdType i = 0; i < splinepolydata->GetPoints()->GetNumberOfPoints(); i++)
    {
        ptId->InsertNextId(i);
    }

    vtkPoints* fp = vtkPoints::New();
    splinepolydata->GetPoints()->GetPoints(ptId, fp);

    // Store the spline points in a vector
    for (vtkIdType i = 0; i < splinepolydata->GetPoints()->GetNumberOfPoints(); i++)
    {
        std::vector<double> tmp;
        tmp.push_back(fp->GetPoint(i)[0]);
        tmp.push_back(fp->GetPoint(i)[1]);
        tmp.push_back(fp->GetPoint(i)[2]);
        all_spline_points.push_back(tmp);
    }

    // Create spheres for each spline point and add them as actors
    for (int i = 0; i < all_spline_points.size(); i++)
    {
        vtkNew<vtkSphereSource> vtk_sphere;
        vtk_sphere->SetCenter(all_spline_points[i][0], all_spline_points[i][1], all_spline_points[i][2]);
        vtk_sphere->SetRadius(static_cast<double>(0.2));
        vtk_sphere->SetPhiResolution(21);
        vtk_sphere->SetThetaResolution(21);
        vtk_sphere->Update();
        vtkNew<vtkPolyDataMapper> vtk_dot_mapper;
        vtk_dot_mapper->SetInputConnection(vtk_sphere->GetOutputPort());
        vtkSmartPointer<vtkActor> dot_actor = vtkSmartPointer<vtkActor>::New();
        dot_actor->GetProperty()->SetColor(color->GetColor3d("Green").GetData());
        dot_actor->SetMapper(vtk_dot_mapper);
        m_renderer->AddActor(dot_actor);
        m_spline_actor_vec.push_back(dot_actor);
    }

    // Create an actor for the spline line
    vtkNew<vtkPolyDataMapper> LineMapper;
    LineMapper->SetInputConnection(functionSource->GetOutputPort());

    vtkNew<vtkActor> LineActor;
    LineActor->SetMapper(LineMapper);
    LineActor->GetProperty()->SetColor(color->GetColor3d("Red").GetData());
    LineActor->GetProperty()->SetLineWidth(1.0);

    // Glyph the points
    vtkNew<vtkSphereSource> PointsSphere;
    PointsSphere->SetPhiResolution(21);
    PointsSphere->SetThetaResolution(21);
    PointsSphere->SetRadius(1);

    // Create a polydata to store the points
    vtkNew<vtkPolyData> PointPolyData;
    PointPolyData->SetPoints(m_teeth_center);

    vtkNew<vtkGlyph3DMapper> pointMapper;
    pointMapper->SetInputData(PointPolyData);
    pointMapper->SetSourceConnection(PointsSphere->GetOutputPort());

    // Create an actor for the points
    vtkNew<vtkActor> PointActor;
    PointActor->SetMapper(pointMapper);
    PointActor->GetProperty()->SetColor(color->GetColor3d("Peacock").GetData());
    PointActor->GetProperty()->SetOpacity(1);

    m_spline_actor_vec.push_back(PointActor);
    m_spline_actor_vec.push_back(LineActor);
    m_renderer->AddActor(PointActor);
    m_renderer->AddActor(LineActor);
}


/**
 * Calculates the direction of teeth.
 *
 * This method calculates the direction of teeth based on the provided data. It initializes the m_occlusal_direction vector
 * and performs various calculations to determine the direction of teeth. The calculated direction is stored in the m_average_occlusal variable.
 */
void TeethDataInitialization::CalculateDirection()
{
    int n_teeth_num = 7;
    for (int i = 0; i < 17; i++)
    {
        double3 t;
        m_occlusal_direction.push_back(t);
    }
    for (int i = 0; i < m_teeth_center->GetNumberOfPoints(); i++)
    {
        double* pt, * pt1, * pt2;
        double point[3], point1[3], point2[3];
        //cout << m_teeth_center->GetNumberOfPoints() << endl;
        pt = m_teeth_center->GetPoint(i);
        point[0] = pt[0]; point[1] = pt[1]; point[2] = pt[2];
        //cout << pt[0] << " " << pt[1] << " " << pt[2] << endl;
        pt1 = m_teeth_center->GetPoint((i + static_cast<int>(3)) % n_teeth_num);
        point1[0] = pt1[0]; point1[1] = pt1[1]; point1[2] = pt1[2];
        //cout << pt1[0] << " " << pt1[1] << " " << pt1[2] << endl;
        pt2 = m_teeth_center->GetPoint((i + static_cast<int>(6)) % n_teeth_num);
        point2[0] = pt2[0]; point2[1] = pt2[1]; point2[2] = pt2[2];
        //cout << pt2[0] << " " << pt2[1] << " " << pt2[2] << endl;
        double vec1[3], vec2[3], ans[3];

        vec1[0] = static_cast<double>(point1[0]) - static_cast<double>(point[0]);
        vec1[1] = static_cast<double>(point1[1]) - static_cast<double>(point[1]);
        vec1[2] = static_cast<double>(point1[2]) - static_cast<double>(point[2]);
        vec2[0] = static_cast<double>(point2[0]) - static_cast<double>(point[0]);
        vec2[1] = static_cast<double>(point2[1]) - static_cast<double>(point[1]);
        vec2[2] = static_cast<double>(point2[2]) - static_cast<double>(point[2]);
        double Occlusal[3] = { 0, 0, 0 };
        Occlusal[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
        Occlusal[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
        Occlusal[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

        double len = std::sqrt(Occlusal[0] * Occlusal[0] + Occlusal[1] * Occlusal[1] + Occlusal[2] * Occlusal[2]);
        Occlusal[0] /= len / 20;
        Occlusal[1] /= len / 20;
        Occlusal[2] /= len / 20;

        Occlusal[0] = -Occlusal[0];
        Occlusal[1] = -Occlusal[1];
        Occlusal[2] = -Occlusal[2];
        m_occlusal_direction[i] = Occlusal;
        m_average_occlusal += m_occlusal_direction[i];
    }
    if (m_checked_id == 2 || m_checked_id == 3)  // If checked_id is 2 or 3, then the teeth model should be rotated.
    {
        m_average_occlusal = -m_average_occlusal;
    }

    m_average_occlusal.normalize();
    std::cout << "average_occlusal: " << m_average_occlusal.x() << " " << m_average_occlusal.y() << " " << m_average_occlusal.z() << endl;
    //for (int i = 0; i < 17; i++)
    //{
    //    double3 t;
    //    m_teeth_arch.push_back(t);
    //}
    m_teeth_arch = std::vector<double3>(17, double3(0, 0, 0));
    for (int i = 1; i < 17; i++)
    {
        if (m_labeledPolyData[i]->GetNumberOfPolys() == 0)
        {
            continue;
        }
        double* pre_box, * cur_box, * next_box;
        pre_box = new double[3];
        next_box = new double[3];
        bool find_pre = 0, find_next = 0;
        cur_box = m_labeledPolyData[i]->GetBounds();
        for (int j = i - 1; j > 0; j--)
        {
            if (m_labeledPolyData[j]->GetNumberOfPolys() == 0)
            {
                continue;
            }
            pre_box = m_labeledPolyData[j]->GetBounds();
            find_pre = 1;
            break;
        }
        for (int j = i + 1; j < 17; j++)
        {
            if (m_labeledPolyData[j]->GetNumberOfPolys() == 0)
            {
                continue;
            }
            next_box = m_labeledPolyData[j]->GetBounds();
            find_next = 1;
            break;
        }
        if (!find_pre)
        {
            double3 cur_center((cur_box[0] + cur_box[1]) / 2.0, (cur_box[2] + cur_box[3]) / 2.0, (cur_box[4] + cur_box[5]) / 2.0);
            double3 next_center((next_box[0] + next_box[1]) / 2.0, (next_box[2] + next_box[3]) / 2.0, (next_box[4] + next_box[5]) / 2.0);
            m_teeth_arch[i] = next_center - cur_center;
        }
        else if (!find_next)
        {
            double3 pre_center((pre_box[0] + pre_box[1]) / 2.0, (pre_box[2] + pre_box[3]) / 2.0, (pre_box[4] + pre_box[5]) / 2.0);
            double3 cur_center((cur_box[0] + cur_box[1]) / 2.0, (cur_box[2] + cur_box[3]) / 2.0, (cur_box[4] + cur_box[5]) / 2.0);
            m_teeth_arch[i] = cur_center - pre_center;
        }
        else
        {
            double3 pre_center((pre_box[0] + pre_box[1]) / 2.0, (pre_box[2] + pre_box[3]) / 2.0, (pre_box[4] + pre_box[5]) / 2.0);
            double3 next_center((next_box[0] + next_box[1]) / 2.0, (next_box[2] + next_box[3]) / 2.0, (next_box[4] + next_box[5]) / 2.0);
            m_teeth_arch[i] = next_center - pre_center;
        }
        //std::cout << "Teeth arch " << i << ": " << m_teeth_arch[i].x() << " " << m_teeth_arch[i].y() << " " << m_teeth_arch[i].z() << endl;
        std::cout << i << " Teeth faces: " << m_labeledPolyData[i]->GetNumberOfPolys() << " Teeth vertices: " << m_labeledPolyData[i]->GetNumberOfPoints() << endl;
    }
}

double* TeethDataInitialization::GetTongueCheekDirection(int n)
{
   return double3::crossProduct(m_average_occlusal, -m_teeth_arch[n]).data;
}

int TeethDataInitialization::GetPointOnNthTeeth(int n)
{
    vtkSmartPointer<vtkCellData> celldata = m_teethPolyData->GetCellData();
    vtkSmartPointer<vtkDataArray> labels = celldata->GetArray("Label");
    for (vtkIdType i = 0; i < labels->GetSize(); i++)
    {
        // tmp is a double array that contains an integer value (in this case, only one value)
        int label = static_cast<int>(labels->GetTuple(i)[0]);
        if (label == n)
            return m_teethPolyData->GetCell(i)->GetPointIds()->GetId(0);
    }
    return -1;
}

void TeethDataInitialization::Execute()
{
    FragmentClearing();
    GetLabel();
    MakeLabeledActor();
    ObbBoundingBox();
    DrawTeethCenterSpline();
    CalculateDirection();
}

void TeethDataInitialization::RemoveAllSplineActors()
{
    for (auto actor : m_spline_actor_vec)
	{
		m_renderer->RemoveActor(actor);
	}
}

