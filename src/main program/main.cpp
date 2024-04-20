#include "stdafx.h"
//#include "vtkRenderPipeline.h"
//#include "meshTransform.h"
//#include "simpleRender.h"
//#include "IOManip.hpp"
//#define ENABLE_TIMER_H
//
//vtkRenderPipeline* pipeline;
//
///**
// * @brief Generate a 4-digit number string with leading zeros.
// *
// * This function takes an integer and converts it into a string representation
// * with a fixed width of 4 characters. If the number has less than 4 digits,
// * leading zeros are added to pad the string to the desired width.
// *
// * @param number The input integer to be converted.
// * @return A string representation of the input number with leading zeros.
// *
// * @note The function uses std::ostringstream, std::setw(), and std::setfill()
// *       to format the output string.
// *
// * @example
// *   int num = 42;
// *   std::string num_str = generate_leading_zero_number_str(num);
// *   // num_str will be "0042"
// */
//std::string generate_leading_zero_number_str(int number)
//{
//    std::ostringstream stream;
//    stream << std::setw(4) << std::setfill('0') << number;
//    return stream.str();
//}
//
//void RightPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//    //std::cout<<"Right Press" << endl;
//}
//void RightRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//    //std::cout<<"Right Released" << endl;
//}
//
//void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//    //std::cout<<"Left Press" << endl;
//}
//
//void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//
//}
//
//void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//}
//
//int SurfaceMeshToPolyDataWithColor(
//    const SurfaceMesh& mesh,
//    vtkSmartPointer<vtkPolyData>& polydata
//)
//{
//    //polyData->Initialize();
//    typedef SurfaceMesh::Property_map<vertex_descriptor, CGAL::Color> Color_map;
//    polydata = vtkSmartPointer<vtkPolyData>::New();
//    vtkSmartPointer<vtkFloatArray> curvature_color_array = vtkSmartPointer<vtkFloatArray>::New();
//    curvature_color_array->SetName("curvature");
//    Color_map color_map;
//    bool found;
//    boost::tie(color_map, found) = mesh.property_map<vertex_descriptor, CGAL::Color>("v:color");
//
//    if (!found)
//    {
//        std::cout << "Property map v:color not found!" << std::endl;
//
//        vtkSmartPointer<vtkPoints> VTKPoints = vtkSmartPointer<vtkPoints>::New();
//        VTKPoints->SetNumberOfPoints(mesh.vertices().size());
//
//        for (int i = 0; i < mesh.vertices().size(); i++)
//        {
//            vertex_descriptor vd = vertex_descriptor(i);
//            const Kernel::Point_3& point = mesh.point(vd);
//            VTKPoints->SetPoint(i, point.x(), point.y(), point.z());
//        }
//
//        vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
//        connectivity->SetNumberOfComponents(4);
//        connectivity->SetNumberOfTuples(mesh.faces().size());
//
//        //#pragma omp parallel for
//        for (int j = 0; j < mesh.faces().size(); j++)
//        {
//            int ids[3] = { 0 };
//            int index = 0;
//            for (const halfedge_descriptor& h : mesh.halfedges_around_face(mesh.halfedge(face_descriptor(j))))
//            {
//                ids[index] = mesh.target(h).idx();
//                index++;
//            }
//            if (index > 3)
//                std::cerr << "SurfaceMeshToPolyData Error";
//
//            connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
//        }
//
//        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
//        cells->SetCells(mesh.faces().size(), connectivity);
//
//        polydata->SetPoints(VTKPoints);
//        polydata->SetPolys(cells);
//
//    }
//    else
//    {
//        // Property map v:color found!
//        vtkSmartPointer<vtkPoints> VTKPoints = vtkSmartPointer<vtkPoints>::New();
//        VTKPoints->SetNumberOfPoints(mesh.vertices().size());
//
//        for (int i = 0; i < mesh.vertices().size(); i++)
//        {
//            vertex_descriptor vd = vertex_descriptor(i);
//            const Kernel::Point_3& point = mesh.point(vd);
//            VTKPoints->SetPoint(i, point.x(), point.y(), point.z());
//            if (color_map[vd] == CGAL::Color(255, 0, 0))
//            {
//                curvature_color_array->InsertNextValue(0.0);
//            }
//            else if (color_map[vd] == CGAL::Color(0, 255, 0))
//            {
//                curvature_color_array->InsertNextValue(1.0);
//            }
//            else
//            {
//                curvature_color_array->InsertNextValue(2.0);
//            }
//        }
//
//        vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
//        connectivity->SetNumberOfComponents(4);
//        connectivity->SetNumberOfTuples(mesh.faces().size());
//
//        //#pragma omp parallel for
//        for (int j = 0; j < mesh.faces().size(); j++)
//        {
//            int ids[3] = { 0 };
//            int index = 0;
//            for (const halfedge_descriptor& h : mesh.halfedges_around_face(mesh.halfedge(face_descriptor(j))))
//            {
//                ids[index] = mesh.target(h).idx();
//                index++;
//            }
//            if (index > 3)
//                std::cerr << "SurfaceMeshToPolyData Error";
//
//            connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
//        }
//
//        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
//        cells->SetCells(mesh.faces().size(), connectivity);
//
//        polydata->SetPoints(VTKPoints);
//        polydata->SetPolys(cells);
//        polydata->GetPointData()->SetScalars(curvature_color_array);
//    }
//
//    return 0;
//}
//
//int main()
//{
//    using namespace Eigen;
//
//    pipeline = new vtkRenderPipeline();
//
//    constexpr double RED_THER = -0.;
//    constexpr double BLUE_THER = 0.2;
//    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
//    reader->SetFileName("data/sample_lower_left_upsampled.vtp");
//    reader->Update();
//    vtkSmartPointer<vtkPolyData> arch_pd = reader->GetOutput();
//
//    SurfaceMesh arch_sm;
//    CGAL::IO::read_VTP("data/sample_lower_left_upsampled.vtp", arch_sm);
//    auto& min_curvature = arch_sm.add_property_map<vertex_descriptor, double>("v:min_curvature", 0.0).first;
//    auto& color_map = arch_sm.add_property_map<vertex_descriptor, CGAL::Color>("v:color", CGAL::Color(255, 255, 255)).first;
//
//    //RenderPolydata(CGAL_Surface_Mesh2VTK_PolyData(arch_sm), pipeline->Renderer);
//
//    MatrixXd HN;
//    SparseMatrix<double> L, M, Minv;
//    Eigen::MatrixXd	m_V; ///< Eigen matrix for vertices.
//    Eigen::MatrixXi	m_F; ///< Eigen matrix for faces.
//
//    CGALSurfaceMeshToEigen(arch_sm, m_V, m_F);
//    igl::cotmatrix(m_V, m_F, L);
//    igl::massmatrix(m_V, m_F, igl::MASSMATRIX_TYPE_VORONOI, M);
//    igl::invert_diag(M, Minv);
//
//    // Laplace-Beltrami of position
//    HN = -Minv * (L * m_V);
//    // Extract magnitude as mean curvature
//    VectorXd H = HN.rowwise().norm();
//
//    VectorXd H5 = H;
//    VectorXd H10 = H;
//    VectorXd H15 = H;
//    VectorXd H20 = H;
//
//    // Compute curvature directions via quadric fitting
//    MatrixXd PD1, PD2;
//    MatrixXd PD1_5, PD2_5;
//    MatrixXd PD1_10, PD2_10;
//    MatrixXd PD1_15, PD2_15;
//    MatrixXd PD1_20, PD2_20;
//    VectorXd PV1, PV2;
//    VectorXd PV1_5, PV2_5;
//    VectorXd PV1_10, PV2_10;
//    VectorXd PV1_15, PV2_15;
//    VectorXd PV1_20, PV2_20;
//    {
//#ifdef ENABLE_TIMER_H
//        Timer timer("Calculate principal_curvature");
//        {
//            Timer timer("Calculate principal_curvature 5");
//#endif
//            igl::principal_curvature(m_V, m_F, PD1_5, PD2_5, PV1_5, PV2_5, 5U, true);
//#ifdef ENABLE_TIMER_H
//        }
//        {
//            Timer timer("Calculate principal_curvature 10");
//#endif
//            igl::principal_curvature(m_V, m_F, PD1_10, PD2_10, PV1_10, PV2_10, 10U, true);
//#ifdef ENABLE_TIMER_H
//        }
//        {
//            Timer timer("Calculate principal_curvature 15");
//#endif
//            igl::principal_curvature(m_V, m_F, PD1_15, PD2_15, PV1_15, PV2_15, 15U, true);
//#ifdef ENABLE_TIMER_H
//        }
//#endif
//    }
//
//    //// mean curvature
//    H5 = 0.5 * (PV1_5 + PV2_5);
//    H10 = 0.5 * (PV1_10 + PV2_10);
//    H15 = 0.5 * (PV1_15 + PV2_15);
//
//    //// Smooth mean curvature
//    //for (vertex_descriptor& vd : arch_sm.vertices())
//    //{
//    //    double h5_sum = 0.0;
//    //    double h10_sum = 0.0;
//    //    double h15_sum = 0.0;
//    //    int h5_count = 0;
//    //    int h10_count = 0;
//    //    int h15_count = 0;
//
//    //    // Traverse adjacent vertices
//    //    auto circ = arch_sm.vertices_around_target(arch_sm.halfedge(vd));
//    //    for (auto vic = circ.begin(); vic != circ.end(); ++vic)
//    //    {
//    //        h5_sum += H5[vic->idx()];
//    //        h10_sum += H10[vic->idx()];
//    //        h15_sum += H15[vic->idx()];
//
//    //        h5_count++;
//    //        h10_count++;
//    //        h15_count++;
//    //    }
//    //    H5[vd.idx()] = h5_sum / h5_count;
//    //    H10[vd.idx()] = h10_sum / h10_count;
//    //    H15[vd.idx()] = h15_sum / h15_count;
//    //}
//
//    for (vertex_descriptor& vd : arch_sm.vertices())
//    {
//        try
//        {
//            double min_value = min_curvature[vd] = std::min(std::min(H5[vd.idx()], H10[vd.idx()]), H15[vd.idx()]);
//            if (min_value < RED_THER)
//            {
//                color_map[vd] = CGAL::Color(255, 0, 0);
//            }
//            else if (min_value > BLUE_THER)
//            {
//                color_map[vd] = CGAL::Color(0, 255, 0);
//            }
//            else
//            {
//                color_map[vd] = CGAL::Color(0, 0, 255);
//            }
//        }
//        catch (const std::out_of_range& e)
//        {
//            std::cerr << "Error: Vertex descriptor: " << vd.idx() << " not found in the mapping!" << std::endl;
//        }
//    }
//
//
//    //for (vertex_descriptor& vd : arch_sm.vertices())
//    //{
//    //    try
//    //    {
//    //        double mean_value = min_curvature[vd] = H15[vd.idx()];
//    //        if (mean_value < RED_THER)
//    //        {
//    //            color_map[vd] = CGAL::Color(255, 0, 0);
//    //        }
//    //        else if (mean_value > BLUE_THER)
//    //        {
//    //            color_map[vd] = CGAL::Color(0, 255, 0);
//    //        }
//    //        else
//    //        {
//    //            color_map[vd] = CGAL::Color(0, 0, 255);
//    //        }
//    //    }
//    //    catch (const std::out_of_range& e)
//    //    {
//    //        std::cerr << "Error: Vertex descriptor: " << vd.idx() << " not found in the mapping!" << std::endl;
//    //    }
//    //}
//
//    vtkSmartPointer<vtkPolyData> colored_arch_pd;
//    SurfaceMeshToPolyDataWithColor(arch_sm, colored_arch_pd);
//
//    vtkSmartPointer<vtkPolyDataMapper> colored_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    colored_mapper->SetInputData(colored_arch_pd);
//    colored_mapper->SetScalarRange(0, 2);
//    //colored_mapper->SetLookupTable(color_lut);
//    colored_mapper->SetInterpolateScalarsBeforeMapping(true);
//
//    vtkSmartPointer<vtkActor> colored_actor = vtkSmartPointer<vtkActor>::New();
//    colored_actor->SetMapper(colored_mapper);
//    colored_actor->GetProperty()->SetOpacity(1);
//    colored_actor->GetProperty()->EdgeVisibilityOff();
//
//    pipeline->Renderer->AddActor(colored_actor);
//
//    // Set up the camera and interactor.
//    pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
//    pipeline->Renderer->ResetCamera();
//    // Set up the callback functions for mouse events.
//    pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
//    pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
//    pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);
//    pipeline->addObserver(vtkCommand::RightButtonPressEvent, RightPress);
//    pipeline->addObserver(vtkCommand::RightButtonReleaseEvent, RightRelease);
//
//    pipeline->RenderWindowInteractor->Start();
//
//    // Clean up the pipeline after each run.
//    delete pipeline;
//}
