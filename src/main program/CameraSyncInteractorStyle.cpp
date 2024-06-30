#include "CameraSyncInteractorStyle.h"

vtkStandardNewMacro(CameraSyncInteractorStyle);

void CameraSyncInteractorStyle::SyncCameras()
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
