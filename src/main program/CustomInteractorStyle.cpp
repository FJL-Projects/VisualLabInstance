#include "CustomInteractorStyle.h"
#include <vtkTransform.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

vtkStandardNewMacro(CustomInteractorStyle);

void CustomInteractorStyle::Rotate()
{
    if (this->CurrentRenderer == nullptr) {
        return;
    }

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    double deltaX = this->Interactor->GetLastEventPosition()[0] - this->Interactor->GetEventPosition()[0];
    double angle = this->MotionFactor * deltaX;

    // Apply rotation to all actors in the renderer
    vtkActorCollection* actors = this->CurrentRenderer->GetActors();
    actors->InitTraversal();
    vtkActor* actor;
    while ((actor = actors->GetNextActor()) != nullptr)
    {
        actor->RotateY(angle);
    }

    this->CurrentRenderer->ResetCameraClippingRange();
    this->Interactor->GetRenderWindow()->Render();

    emit SignalSyncRotate(angle);
}

void CustomInteractorStyle::SlotSyncRotate(double angle)
{
    if (this->CurrentRenderer == nullptr) {
        return;
    }

    // Apply rotation to all actors in the renderer
    vtkActorCollection* actors = this->CurrentRenderer->GetActors();
    actors->InitTraversal();
    vtkActor* actor;
    while ((actor = actors->GetNextActor()) != nullptr)
    {
        actor->RotateY(angle);
    }

    this->CurrentRenderer->ResetCameraClippingRange();
    this->Interactor->GetRenderWindow()->Render();
}

void CustomInteractorStyle::OnMouseMove()
{
    this->CurrentRenderer = this->Interactor->FindPokedRenderer(
        this->Interactor->GetEventPosition()[0],
        this->Interactor->GetEventPosition()[1]
    );
    if (this->CurrentRenderer != nullptr) {
        vtkInteractorStyleTrackballActor::OnMouseMove();
        this->CurrentRenderer->ResetCameraClippingRange();
        this->Interactor->GetRenderWindow()->Render();
    }
}

void CustomInteractorStyle::OnLeftButtonDown()
{
    this->CurrentRenderer = this->Interactor->FindPokedRenderer(
        this->Interactor->GetEventPosition()[0],
        this->Interactor->GetEventPosition()[1]
    );
    if (this->CurrentRenderer != nullptr) {
        vtkInteractorStyleTrackballActor::OnLeftButtonDown();
    }
}

void CustomInteractorStyle::OnLeftButtonUp()
{
    this->CurrentRenderer = this->Interactor->FindPokedRenderer(
        this->Interactor->GetEventPosition()[0],
        this->Interactor->GetEventPosition()[1]
    );
    if (this->CurrentRenderer != nullptr) {
        vtkInteractorStyleTrackballActor::OnLeftButtonUp();
        this->CurrentRenderer->ResetCameraClippingRange();
        this->Interactor->GetRenderWindow()->Render();
    }
}
