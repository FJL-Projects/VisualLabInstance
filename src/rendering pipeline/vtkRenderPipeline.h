#include"stdafx.h"

class DesignInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static DesignInteractorStyle* New()
    {
        return new DesignInteractorStyle;
    }
    vtkTypeMacro(DesignInteractorStyle, vtkInteractorStyleTrackballCamera);

    DesignInteractorStyle() {}
    virtual ~DesignInteractorStyle() {}
    virtual void OnLeftButtonDown() {}
    virtual void OnLeftButtonUp() {}
    virtual void OnRightButtonDown() { this->StartRotate(); }
    virtual void OnRightButtonUp() { this->vtkInteractorStyleTrackballCamera::OnLeftButtonUp(); }
    virtual void OnMouseMove() { this->vtkInteractorStyleTrackballCamera::OnMouseMove(); }
    virtual void OnMouseWheelForward() { if (!this->Interactor->GetControlKey() && !this->Interactor->GetShiftKey()) this->vtkInteractorStyleTrackballCamera::OnMouseWheelForward(); }
    virtual void OnMouseWheelBackward() { if (!this->Interactor->GetControlKey() && !this->Interactor->GetShiftKey()) this->vtkInteractorStyleTrackballCamera::OnMouseWheelBackward(); }
};

class vtkRenderPipeline
{
public:
    vtkRenderPipeline()
    {
		this->m_render_window = vtkSmartPointer<vtkRenderWindow>::New();
		this->m_renderer = vtkSmartPointer<vtkRenderer>::New();
		this->m_cell_picker = vtkSmartPointer<vtkCellPicker>::New();
		this->m_renderer->SetBackground(0.41, 0.41, 0.41);
		this->m_render_window->AddRenderer(this->m_renderer);
		this->m_render_window->SetSize(1920, 1080);
		//this->m_render_window->Render();
		this->InteractorStyle = vtkSmartPointer<DesignInteractorStyle>::New();
		this->RenderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		this->RenderWindowInteractor->SetInteractorStyle(this->InteractorStyle);
		this->RenderWindowInteractor->SetRenderWindow(this->m_render_window);
	}
    void addObserver(unsigned long event,void (*f)(vtkObject* caller, unsigned long eid,void* clientdata, void* calldata))
    {
		vtkNew<vtkCallbackCommand> callback;
		callback->SetCallback(f);
		callback->InitializeObjectBase();
		this->RenderWindowInteractor->AddObserver(event, callback);
	}

    vtkSmartPointer<vtkRenderWindow> m_render_window;
    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkCellPicker> m_cell_picker;
    vtkSmartPointer<DesignInteractorStyle> InteractorStyle;
    vtkSmartPointer<vtkRenderWindowInteractor> RenderWindowInteractor;
};