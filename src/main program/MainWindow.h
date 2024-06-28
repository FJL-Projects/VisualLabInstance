#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <QVTKWidget.h>
#include <vtkCallbackCommand.h>

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

private slots:
    void OnLoadButtonClicked();

private:
    QVBoxLayout* m_p_main_layout;
    QPushButton* m_p_load_button;
    QVTKWidget* m_p_vtk_widget_left;
    QVTKWidget* m_p_vtk_widget_right;
    vtkSmartPointer<vtkRenderer> m_p_renderer_left;
    vtkSmartPointer<vtkRenderer> m_p_renderer_right;
    vtkSmartPointer<vtkActor> m_p_actor;
    vtkSmartPointer<vtkCallbackCommand> m_p_camera_callback;

    void LoadSTLFile(const QString& file_path);
    void SynchronizeCameras();

    static void CameraModifiedCallback(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
    void SynchronizeCamerasWithParallax(double angle);
};

#endif // MAINWINDOW_H