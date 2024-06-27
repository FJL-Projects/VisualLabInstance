#pragma once
#include <QApplication>
#include <QMainWindow>
#include <QSurfaceFormat>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPushButton>
#include <QFileDialog>
#include <QVTKOpenGLWidget.h>
#include <QDebug>
#ifndef Q_MOC_RUN
#if defined(emit)
#undef emit
#include <tbb/tbb.h>
#define emit // restore the macro definition of "emit", as it was defined in gtmetamacros.h
#else
#include <tbb/tbb.h>
#endif // defined(emit)
#endif // Q_MOC_RUN
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkNew.h>
#include <vtkCamera.h>
#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAxesActor.h>
#include <vtkTextProperty.h>
#include <vtkAxesActor.h>
#include <vtkTextProperty.h>
#include <vtkTextActor.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCenterOfMass.h>
#include <vtkMath.h>

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();

private slots:
    void openFile();

private:
    QVTKOpenGLWidget* leftVtkWidget;
    vtkRenderer* leftRenderer;
    QVTKOpenGLWidget* rightVtkWidget;
    vtkRenderer* rightRenderer;

    void AddAxesToRenderer(vtkRenderer* renderer);
};