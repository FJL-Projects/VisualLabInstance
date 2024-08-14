#ifndef CAMERASYNCINTERACTORSTYLE_H
#define CAMERASYNCINTERACTORSTYLE_H

#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkMath.h>

class CameraSyncInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static CameraSyncInteractorStyle* New();
    vtkTypeMacro(CameraSyncInteractorStyle, vtkInteractorStyleTrackballCamera);

    void SetOtherRenderer(vtkRenderer* renderer)
    {
        this->OtherRenderer = renderer;
    }

    void SetAngleDifference(double angle)
    {
        this->AngleDifference = angle;
    }

    virtual void OnMouseMove() override
    {
        vtkInteractorStyleTrackballCamera::OnMouseMove();
        if (this->OtherRenderer)
        {
            SyncCameras();
        }
    }

    virtual void OnLeftButtonDown() override
    {
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
        if (this->OtherRenderer)
        {
            SyncCameras();
        }
    }

    virtual void OnLeftButtonUp() override
    {
        vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
        if (this->OtherRenderer)
        {
            SyncCameras();
        }
    }

    virtual void OnRightButtonDown() override
    {
        vtkInteractorStyleTrackballCamera::OnRightButtonDown();
        if (this->OtherRenderer)
        {
            SyncCameras();
        }
    }

    virtual void OnRightButtonUp() override
    {
        vtkInteractorStyleTrackballCamera::OnRightButtonUp();
        if (this->OtherRenderer)
        {
            SyncCameras();
        }
    }

    virtual void OnMiddleButtonDown() override
    {
        vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
        if (this->OtherRenderer)
        {
            SyncCameras();
        }
    }

    virtual void OnMiddleButtonUp() override
    {
        vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
        if (this->OtherRenderer)
        {
            SyncCameras();
        }
    }

private:
    vtkRenderer* OtherRenderer;
    double AngleDifference;

    void SyncCameras()
    {
        vtkRenderer* currentRenderer = this->GetDefaultRenderer();
        if (!currentRenderer || !this->OtherRenderer)
        {
            return;
        }

        vtkCamera* camera = currentRenderer->GetActiveCamera();
        vtkCamera* otherCamera = this->OtherRenderer->GetActiveCamera();

        // 复制相机参数
        double position[3];
        double focalPoint[3];
        double viewUp[3];

        camera->GetPosition(position);
        camera->GetFocalPoint(focalPoint);
        camera->GetViewUp(viewUp);

        // 应用角度差异
        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform->PostMultiply();
        transform->Translate(focalPoint);
        transform->RotateY(this->AngleDifference);
        transform->Translate(-focalPoint[0], -focalPoint[1], -focalPoint[2]);

        double newPosition[3];
        transform->TransformPoint(position, newPosition);

        otherCamera->SetPosition(newPosition);
        otherCamera->SetFocalPoint(focalPoint);
        otherCamera->SetViewUp(viewUp);
        otherCamera->SetClippingRange(camera->GetClippingRange());
        otherCamera->SetParallelScale(camera->GetParallelScale());

        this->OtherRenderer->ResetCameraClippingRange();
    }
};

vtkStandardNewMacro(CameraSyncInteractorStyle);

#endif // CAMERASYNCINTERACTORSTYLE_H
