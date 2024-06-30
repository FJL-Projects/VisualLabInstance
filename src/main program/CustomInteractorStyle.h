#ifndef CUSTOMINTERACTORSTYLE_H
#define CUSTOMINTERACTORSTYLE_H

#include <vtkInteractorStyleTrackballActor.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkObjectFactory.h>
#include <QObject>

class CustomInteractorStyle : public QObject, public vtkInteractorStyleTrackballActor
{
    Q_OBJECT
public:
    static CustomInteractorStyle* New();
    vtkTypeMacro(CustomInteractorStyle, vtkInteractorStyleTrackballActor);

    virtual void Rotate() override;

    void OnMouseMove() override;

    void OnLeftButtonDown() override;

    void OnLeftButtonUp() override;

signals:
    void SignalSyncRotate(double angle);

public slots:
    void SlotSyncRotate(double angle);

private:
    double MotionFactor = 0.1;
};

#endif // CUSTOMINTERACTORSTYLE_H