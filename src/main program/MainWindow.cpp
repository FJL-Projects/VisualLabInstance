#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2); // Necessary for VTK rendering
VTK_MODULE_INIT(vtkInteractionStyle);

#pragma warning(push)
#pragma warning(disable : 4996) // 忽略QVTKWidget的deprecated警告

#include <QVTKWidget.h>

#include "MainWindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <vtkRenderWindow.h>
#include <vtkCamera.h>
#include <vtkTransform.h>
#include <vtkInteractorStyleTrackballActor.h>
#include "CustomInteractorStyle.h"

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent),
    m_p_main_layout(new QVBoxLayout),
    m_p_load_button(new QPushButton("Load STL File")),
    m_p_renderer_left(vtkSmartPointer<vtkRenderer>::New()),
    m_p_renderer_right(vtkSmartPointer<vtkRenderer>::New()),
    m_p_actor(vtkSmartPointer<vtkActor>::New())
{
    QWidget* centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    centralWidget->setLayout(m_p_main_layout);

    m_p_main_layout->addWidget(m_p_load_button);

    QHBoxLayout* h_layout = new QHBoxLayout;
    m_p_main_layout->addLayout(h_layout);

    m_p_vtk_widget_left = new QVTKWidget(this);
    m_p_vtk_widget_right = new QVTKWidget(this);

    m_p_vtk_widget_left->setFixedSize(800, 800);
    m_p_vtk_widget_right->setFixedSize(800, 800);

    h_layout->addWidget(m_p_vtk_widget_left);
    h_layout->addWidget(m_p_vtk_widget_right);

    InitializeRenderers();

    connect(m_p_load_button, &QPushButton::clicked, this, &MainWindow::LoadSTLFile);
}

MainWindow::~MainWindow()
{
}

void MainWindow::InitializeRenderers()
{
    m_p_vtk_widget_left->GetRenderWindow()->AddRenderer(m_p_renderer_left);
    m_p_vtk_widget_right->GetRenderWindow()->AddRenderer(m_p_renderer_right);

    // Setup default camera positions with parallax
    vtkSmartPointer<vtkCamera> camera_left = vtkSmartPointer<vtkCamera>::New();
    vtkSmartPointer<vtkCamera> camera_right = vtkSmartPointer<vtkCamera>::New();

    camera_left->Azimuth(-4.0);  // 左视角 -4度
    camera_right->Azimuth(4.0);  // 右视角 4度

    m_p_renderer_left->SetActiveCamera(camera_left);
    m_p_renderer_right->SetActiveCamera(camera_right);

    SetupInteractorStyle();
}

void MainWindow::SetupInteractorStyle()
{
    m_p_interactor_style_left = vtkSmartPointer<CustomInteractorStyle>::New();
    m_p_interactor_style_right = vtkSmartPointer<CustomInteractorStyle>::New();

    m_p_vtk_widget_left->GetInteractor()->SetInteractorStyle(m_p_interactor_style_left);
    m_p_vtk_widget_right->GetInteractor()->SetInteractorStyle(m_p_interactor_style_right);

    ConnectInteractorStyles();
}
void MainWindow::ConnectInteractorStyles()
{
    // Connect the signal and slots between the two interactor styles for synchronization
    connect(m_p_interactor_style_left, &CustomInteractorStyle::SignalSyncRotate,
        m_p_interactor_style_right, &CustomInteractorStyle::SlotSyncRotate);

    connect(m_p_interactor_style_right, &CustomInteractorStyle::SignalSyncRotate,
        m_p_interactor_style_left, &CustomInteractorStyle::SlotSyncRotate);
}

void MainWindow::LoadSTLFile()
{
    QString file_path = QFileDialog::getOpenFileName(this, "Open STL File", "", "STL Files (*.stl)");
    if (file_path.isEmpty())
    {
        return;
    }

    auto reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(file_path.toStdString().c_str());
    reader->Update();

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    m_p_actor->SetMapper(mapper);

    // Add the actor to both renderers
    m_p_renderer_left->AddActor(m_p_actor);
    m_p_renderer_right->AddActor(m_p_actor);

    m_p_renderer_left->ResetCamera();
    m_p_renderer_right->ResetCamera();

    m_p_vtk_widget_left->GetRenderWindow()->Render();
    m_p_vtk_widget_right->GetRenderWindow()->Render();
}

void MainWindow::RotateActorForParallax(vtkSmartPointer<vtkActor> actor, double angle)
{
    actor->RotateY(angle);
}
#pragma warning(pop)
