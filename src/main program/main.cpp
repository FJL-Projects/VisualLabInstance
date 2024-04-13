#include "stdafx.h"
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
 
int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

    QMainWindow mainWin;
    mainWin.resize(800, 600);

    QWidget* centralWidget = new QWidget(&mainWin);
    mainWin.setCentralWidget(centralWidget);

    QVTKOpenGLWidget* mainVtkWidget = new QVTKOpenGLWidget(centralWidget);
    mainVtkWidget->resize(800, 600);

    vtkNew<vtkRenderer> mainRenderer;
    mainVtkWidget->GetRenderWindow()->AddRenderer(mainRenderer);
    mainRenderer->SetBackground(0.1, 0.2, 0.4); // RGB color

    vtkNew<vtkConeSource> coneSource;
    coneSource->SetHeight(3.0);
    coneSource->SetRadius(1.0);
    coneSource->SetResolution(10);

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(coneSource->GetOutputPort());

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    mainRenderer->AddActor(actor);

    // Create a smaller VTK widget for the bottom-right corner
    QVTKOpenGLWidget* smallVtkWidget = new QVTKOpenGLWidget(centralWidget);
    smallVtkWidget->setFixedSize(200, 150); // Fixed size
    smallVtkWidget->move(600, 450); // Move to bottom-right corner

    vtkNew<vtkRenderer> smallRenderer;
    smallVtkWidget->GetRenderWindow()->AddRenderer(smallRenderer);
    smallRenderer->SetBackground(0.2, 0.3, 0.5); // Different background color

    vtkNew<vtkConeSource> smallConeSource;
    smallConeSource->SetHeight(2.0);
    smallConeSource->SetRadius(0.8);
    smallConeSource->SetResolution(10);

    vtkNew<vtkPolyDataMapper> smallMapper;
    smallMapper->SetInputConnection(smallConeSource->GetOutputPort());

    vtkNew<vtkActor> smallActor;
    smallActor->SetMapper(smallMapper);
    smallRenderer->AddActor(smallActor);

    // Setting a different interactor style (e.g., Trackball Camera)
    vtkNew<vtkInteractorStyleTrackballCamera> smallInteractorStyle;
    smallVtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(smallInteractorStyle);

    mainWin.show();

    return app.exec();
}