#include "mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#pragma warning(push)
#pragma warning(disable: 4996)

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent),
    m_p_main_layout(new QVBoxLayout), // 主布局改为垂直布局
    m_p_load_button(new QPushButton("Load STL File")),
    m_p_renderer_left(vtkSmartPointer<vtkRenderer>::New()),
    m_p_renderer_right(vtkSmartPointer<vtkRenderer>::New()),
    m_p_actor(vtkSmartPointer<vtkActor>::New()),
    m_p_camera_callback(vtkSmartPointer<vtkCallbackCommand>::New())
{
    m_p_vtk_widget_left = new QVTKWidget(this, Qt::WindowFlags());
    m_p_vtk_widget_right = new QVTKWidget(this, Qt::WindowFlags());

    // 设置 QVTKWidget 的初始大小为 800x800
    m_p_vtk_widget_left->setFixedSize(800, 800);
    m_p_vtk_widget_right->setFixedSize(800, 800);

    QWidget* central_widget = new QWidget(this);
    setCentralWidget(central_widget);
    central_widget->setLayout(m_p_main_layout);

    // 创建一个水平布局用于放置按钮
    QHBoxLayout* button_layout = new QHBoxLayout;
    button_layout->addWidget(m_p_load_button);
    button_layout->addStretch(); // 添加伸展以将按钮居左

    // 创建一个水平布局用于放置两个 QVTKWidget
    QHBoxLayout* vtk_layout = new QHBoxLayout;
    vtk_layout->addWidget(m_p_vtk_widget_left);
    vtk_layout->addWidget(m_p_vtk_widget_right);

    // 将按钮布局和 VTK 布局添加到主垂直布局中
    m_p_main_layout->addLayout(button_layout);
    m_p_main_layout->addLayout(vtk_layout);

    m_p_vtk_widget_left->GetRenderWindow()->AddRenderer(m_p_renderer_left);
    m_p_vtk_widget_right->GetRenderWindow()->AddRenderer(m_p_renderer_right);

    connect(m_p_load_button, &QPushButton::clicked, this, &MainWindow::OnLoadButtonClicked);

    m_p_camera_callback->SetCallback(MainWindow::CameraModifiedCallback);
    m_p_camera_callback->SetClientData(this);
    m_p_renderer_left->GetActiveCamera()->AddObserver(vtkCommand::ModifiedEvent, m_p_camera_callback);
}

MainWindow::~MainWindow()
{
}

void MainWindow::OnLoadButtonClicked()
{
    QString file_path = QFileDialog::getOpenFileName(this, tr("Open STL File"), "", tr("STL Files (*.stl)"));
    if (!file_path.isEmpty()) {
        LoadSTLFile(file_path);
    }
}

void MainWindow::LoadSTLFile(const QString& file_path)
{
    auto reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(file_path.toStdString().c_str());
    reader->Update();

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    m_p_actor->SetMapper(mapper);
    m_p_actor->RotateY(8); // Apply the 8-degree triangular parallax

    m_p_renderer_left->AddActor(m_p_actor);
    m_p_renderer_left->ResetCamera();

    // Create a new actor for the right renderer to avoid sharing the same actor
    auto actor_right = vtkSmartPointer<vtkActor>::New();
    actor_right->ShallowCopy(m_p_actor);
    m_p_renderer_right->AddActor(actor_right);
    m_p_renderer_right->ResetCamera();

    SynchronizeCamerasWithParallax(-8.0);

    m_p_vtk_widget_left->GetRenderWindow()->Render();
    m_p_vtk_widget_right->GetRenderWindow()->Render();
}

void MainWindow::SynchronizeCameras()
{
    vtkCamera* left_camera = m_p_renderer_left->GetActiveCamera();
    vtkCamera* right_camera = m_p_renderer_right->GetActiveCamera();

    right_camera->SetPosition(left_camera->GetPosition());
    right_camera->SetFocalPoint(left_camera->GetFocalPoint());
    right_camera->SetViewUp(left_camera->GetViewUp());
    right_camera->SetParallelScale(left_camera->GetParallelScale());

    m_p_vtk_widget_right->GetRenderWindow()->Render();
}

void MainWindow::CameraModifiedCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
    MainWindow* self = static_cast<MainWindow*>(clientData);

    // Synchronize the right camera with a parallax angle of 8 degrees
    self->SynchronizeCamerasWithParallax(-8.0);

    // Only render the right window once after synchronization
    self->m_p_vtk_widget_right->GetRenderWindow()->Render();
}

void MainWindow::SynchronizeCamerasWithParallax(double angle)
{
    vtkCamera* left_camera = m_p_renderer_left->GetActiveCamera();
    vtkCamera* right_camera = m_p_renderer_right->GetActiveCamera();

    // Get the current camera position and focal point
    double left_position[3];
    double left_focal_point[3];
    double right_position[3];

    left_camera->GetPosition(left_position);
    left_camera->GetFocalPoint(left_focal_point);

    // Compute the right camera position by rotating the left camera position around the focal point
    double angle_rad = vtkMath::RadiansFromDegrees(angle);
    double cos_angle = cos(angle_rad);
    double sin_angle = sin(angle_rad);

    // Compute the right camera position based on the 8 degree parallax
    right_position[0] = cos_angle * (left_position[0] - left_focal_point[0]) - sin_angle * (left_position[2] - left_focal_point[2]) + left_focal_point[0];
    right_position[1] = left_position[1];
    right_position[2] = sin_angle * (left_position[0] - left_focal_point[0]) + cos_angle * (left_position[2] - left_focal_point[2]) + left_focal_point[2];

    right_camera->SetPosition(right_position);
    right_camera->SetFocalPoint(left_focal_point);
    right_camera->SetViewUp(left_camera->GetViewUp());
    right_camera->SetParallelScale(left_camera->GetParallelScale());
}
#pragma warning(pop)