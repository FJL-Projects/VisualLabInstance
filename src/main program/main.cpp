#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
//#include "simpleRender.h"
#include "IOManip.hpp"
#include <MRMesh/MRMeshLoad.h>
#include <MRMesh/MRId.h>
#include <MRMesh/MRMesh.h>
#include <MRMesh/MRBitSetParallelFor.h>
#include <MRMesh/MRMeshTopology.h>
#include <MRMesh/MRExpected.h>

#include <QApplication>
#include <QMainWindow>
#include <QSurfaceFormat>

#include <QVTKOpenGLWidget.h>
#include <vtkAutoInit.h>
#include <vtkConeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <QHBoxLayout>
#include <vtkCommand.h>
#include "MainWindow.h"

vtkRenderPipeline* pipeline;



void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Press" << endl;
}
void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Released" << endl;
}
void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Released" << endl;
}
void RightPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
    //std::cout<<"Right Press" << endl;
}
void RightRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
    //std::cout<<"Right Released" << endl;
}

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

    MainWindow main_window;
    main_window.show();

    return app.exec();
}


//
//int main(int argc, char* argv[])
//{
//    QApplication app(argc, argv);
//
//    QMainWindow mainWin;
//    mainWin.resize(1600, 600); // 设置窗口大小以容纳两个并排窗口
//
//    QWidget* centralWidget = new QWidget(&mainWin);
//    mainWin.setCentralWidget(centralWidget);
//
//    QHBoxLayout* layout = new QHBoxLayout(centralWidget);
//
//    // 创建第一个 QVTKOpenGLWidget
//    QVTKOpenGLWidget* leftVtkWidget = new QVTKOpenGLWidget(centralWidget);
//    layout->addWidget(leftVtkWidget);
//
//    vtkNew<vtkRenderer> leftRenderer;
//    leftVtkWidget->GetRenderWindow()->AddRenderer(leftRenderer);
//    leftRenderer->SetBackground(0.1, 0.1, 0.1); // RGB color
//
//    // 创建第二个 QVTKOpenGLWidget
//    QVTKOpenGLWidget* rightVtkWidget = new QVTKOpenGLWidget(centralWidget);
//    layout->addWidget(rightVtkWidget);
//
//    vtkNew<vtkRenderer> rightRenderer;
//    rightVtkWidget->GetRenderWindow()->AddRenderer(rightRenderer);
//    rightRenderer->SetBackground(0.1, 0.1, 0.1); // RGB color
//
//    using namespace MR;
//
//    // 加载并处理网格数据
//    MR::Mesh mrmesh = *MR::MeshLoad::fromAnySupportedFormat("data/11.stl");
//    std::vector<vertex_descriptor> vertices_list(mrmesh.topology.vertSize());
//
//    SurfaceMesh sm = MRMeshToSurfaceMesh(mrmesh);
//    Mesh converted = SurfaceMeshToMRMesh(sm);
//
//    MeshSave::toAnySupportedFormat(converted, "data/converted.stl");
//
//    RenderPolydata(MRMeshToPolyData(mrmesh), leftRenderer);
//    RenderPolydata(MRMeshToPolyData(mrmesh), rightRenderer);
//
//    vtkNew<vtkSTLReader> reader;
//    reader->SetFileName("data/converted.stl");
//    reader->Update();
//    vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
//
//    Mesh converted_from_pd = PolyDataToMRMesh(polydata);
//    MeshSave::toAnySupportedFormat(converted_from_pd, "data/converted_from_pd.ply");
//
//    vtkCamera* leftCamera = leftRenderer->GetActiveCamera();
//    vtkCamera* rightCamera = rightRenderer->GetActiveCamera();
//
//    leftCamera->Azimuth(-4.0); // 左摄像机向左旋转4度
//    rightCamera->Azimuth(4.0); // 右摄像机向右旋转4度
//
//    // 使两台摄像机保持同步
//    rightCamera->SetViewUp(leftCamera->GetViewUp());
//    rightCamera->SetPosition(leftCamera->GetPosition());
//    rightCamera->SetFocalPoint(leftCamera->GetFocalPoint());
//
//    // 重新调整右摄像机的位置和焦点
//    rightCamera->Azimuth(8.0); // 增加8度的差异
//
//    // 创建自定义交互器并设置它们
//    vtkNew<CameraSyncInteractorStyle> leftInteractorStyle;
//    vtkNew<CameraSyncInteractorStyle> rightInteractorStyle;
//
//    leftInteractorStyle->SetDefaultRenderer(leftRenderer);
//    rightInteractorStyle->SetDefaultRenderer(rightRenderer);
//
//    leftInteractorStyle->SetOtherRenderer(rightRenderer);
//    rightInteractorStyle->SetOtherRenderer(leftRenderer);
//
//    leftVtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(leftInteractorStyle);
//    rightVtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(rightInteractorStyle);
//
//    // 渲染窗口
//    leftRenderer->ResetCamera();
//    rightRenderer->ResetCamera();
//
//    leftVtkWidget->GetRenderWindow()->Render();
//    rightVtkWidget->GetRenderWindow()->Render();
//
//    mainWin.show();
//
//    return app.exec();
//}


//int main(int argc, char* argv[])
//{
//    QApplication app(argc, argv);
//
//    QMainWindow mainWin;
//    mainWin.resize(800, 600);
//
//    QWidget* centralWidget = new QWidget(&mainWin);
//    mainWin.setCentralWidget(centralWidget);
//
//    QVTKOpenGLWidget* mainVtkWidget = new QVTKOpenGLWidget(centralWidget);
//    mainVtkWidget->resize(800, 600);
//
//    vtkNew<vtkRenderer> mainRenderer;
//    mainVtkWidget->GetRenderWindow()->AddRenderer(mainRenderer);
//    mainRenderer->SetBackground(0.1, 0.1, 0.1); // RGB color
//
//	using namespace MR;
//	//pipeline = new vtkRenderPipeline();
//
//	MR::Mesh mrmesh = *MR::MeshLoad::fromAnySupportedFormat("data/11.stl");
//	std::vector<vertex_descriptor> vertices_list(mrmesh.topology.vertSize());
//	
//	SurfaceMesh sm = MRMeshToSurfaceMesh(mrmesh);
//	Mesh converted = SurfaceMeshToMRMesh(sm);
//	
//	MeshSave::toAnySupportedFormat(converted, "data/converted.stl");
//
//	RenderPolydata(MRMeshToPolyData(mrmesh), mainRenderer);
//
//	vtkNew<vtkSTLReader> reader;
//	reader->SetFileName("data/converted.stl");
//	reader->Update();
//	vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
//
//	Mesh converted_from_pd = PolyDataToMRMesh(polydata);
//	MeshSave::toAnySupportedFormat(converted_from_pd, "data/converted_from_pd.ply");
//
//	//// Set up the camera and interactor.
//	//pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
//	//pipeline->Renderer->ResetCamera();
//	//// Set up the callback functions for mouse events.
//	//pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
//	//pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
//	//pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);
//	//pipeline->addObserver(vtkCommand::RightButtonPressEvent, RightPress);
//	//pipeline->addObserver(vtkCommand::RightButtonReleaseEvent, RightRelease);
//
//	//pipeline->RenderWindowInteractor->Start();
//
//	//// Clean up the pipeline after each run.
//	//delete pipeline;
//
//
//    mainWin.show();
//
//    return app.exec();
//}