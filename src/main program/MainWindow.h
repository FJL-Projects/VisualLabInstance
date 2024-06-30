#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <QVTKWidget.h>
#include "CustomInteractorStyle.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

private slots:
    void LoadSTLFile();

private:
    QVBoxLayout* m_p_main_layout;
    QPushButton* m_p_load_button;
    QVTKWidget* m_p_vtk_widget_left;
    QVTKWidget* m_p_vtk_widget_right;
    vtkSmartPointer<vtkRenderer> m_p_renderer_left;
    vtkSmartPointer<vtkRenderer> m_p_renderer_right;
    vtkSmartPointer<vtkActor> m_p_actor;

    vtkSmartPointer<CustomInteractorStyle> m_p_interactor_style_left;
    vtkSmartPointer<CustomInteractorStyle> m_p_interactor_style_right;

    void InitializeRenderers();
    void SetupInteractorStyle();
    void RotateActorForParallax(vtkSmartPointer<vtkActor> actor, double angle);

    void ConnectInteractorStyles();
};

#endif // MAINWINDOW_H