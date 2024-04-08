#include"stdafx.h"
#include"vtkRenderPipeline.h"
#include"meshTransform.h"
#include"simpleRender.h"
#include "TeethWrapper.h"

vtkRenderPipeline* pipeline;
SurfaceMesh mesh;

TeethWrapper teeth_wrapper;

vtkSmartPointer<vtkPolyData> PolyData;
vtkSmartPointer<vtkActor> PolyDataActor;
vtkSmartPointer<vtkActor> TeethPolyDataActor;
vtkSmartPointer<vtkActor> RemeshedTeethPolyDataActor;
vtkSmartPointer<vtkPolyData> TeethPolydata;

Vector_3 projection_direction(-0.037686, 0.379984, -0.203274);

bool rotate_flag = true;
void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
    if (rotate_flag)
    {
        pipeline->Renderer->RemoveActor(TeethPolyDataActor);
        pipeline->Renderer->AddActor(RemeshedTeethPolyDataActor);
        pipeline->RenderWindow->Render();
        rotate_flag = false;
    }
    else
    {
        pipeline->Renderer->RemoveActor(RemeshedTeethPolyDataActor);
        pipeline->Renderer->AddActor(TeethPolyDataActor);
        pipeline->RenderWindow->Render();
        rotate_flag = true;
    }
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	
}


void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}


int main()
{
	pipeline = new vtkRenderPipeline();

    vtkNew<vtkXMLPolyDataReader> vtp_reader;
    vtp_reader->SetFileName("data/arch_38057.vtp");
    vtp_reader->Update();
    PolyData = vtp_reader->GetOutput();
    std::cout << "arch.vtp" << std::endl;
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(PolyData);
    mapper->Update();
    PolyDataActor = vtkSmartPointer<vtkActor>::New();
    PolyDataActor->SetMapper(mapper);
    PolyDataActor->GetProperty()->SetColor(1, 1, 1);
    PolyDataActor->GetProperty()->SetAmbient(0.5);
    PolyDataActor->GetProperty()->SetSpecularPower(100);
    PolyDataActor->GetProperty()->SetSpecular(0.5);
    PolyDataActor->GetProperty()->SetDiffuse(0.5);
    PolyDataActor->GetProperty()->EdgeVisibilityOff();
    PolyDataActor->PickableOn();

    SurfaceMesh arch_sm;
    //PolyDataToSurfaceMesh(PolyData, arch_sm);
    arch_sm = VTK_PolyData2CGAL_Surface_Mesh(PolyData);

    std::ifstream visited_vd_reader("data/visited_vd_38057.txt");
    std::set<vertex_descriptor> arch_without_abutment_set;
    {
        std::string line;
        while (std::getline(visited_vd_reader, line))
        {
            std::istringstream iss(line);
            int x;
            if (!(iss >> x)) { break; } // error
            arch_without_abutment_set.insert(vertex_descriptor(x));
        }
    }

    vtkSmartPointer<vtkPLYReader> teeth_reader = vtkSmartPointer<vtkPLYReader>::New();
    teeth_reader->SetFileName("data/teeth_38057.ply");
    teeth_reader->Update();
    TeethPolydata = teeth_reader->GetOutput();

    SurfaceMesh teeth_sm;
    teeth_sm = VTK_PolyData2CGAL_Surface_Mesh(TeethPolydata);

    teeth_wrapper.SetTeethSurfaceMesh(teeth_sm);
    teeth_wrapper.SetTeethActor(TeethPolyDataActor);
    teeth_wrapper.SetRenderer(pipeline->Renderer);
    teeth_wrapper.SetRenderWindow(pipeline->RenderWindow);
    teeth_wrapper.SetProjectionDirection(projection_direction);
    //teeth_wrapper.SetMaxDistance(max_distance);
    teeth_wrapper.SetOpacity(1);

    teeth_wrapper.SetSelectedId(6);
    teeth_wrapper.SetArchPolyData(PolyData);
    teeth_wrapper.SetArchSurfaceMesh(arch_sm);
    teeth_wrapper.SetArchWithoutAbutmentSet(arch_without_abutment_set);
    teeth_wrapper.ExtractAdjacentTeethWithoutAbutment();
    teeth_wrapper.ProximalShaving(0.08);
    SurfaceMesh remeshed_sm = teeth_wrapper.GetTeethSurfaceMesh();
    vtkSmartPointer<vtkPolyData> remeshed_polydata = CGAL_Surface_Mesh2VTK_PolyData(remeshed_sm);
    vtkSmartPointer<vtkPolyDataMapper> remeshed_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    remeshed_mapper->SetInputData(remeshed_polydata);
    remeshed_mapper->Update();
    RemeshedTeethPolyDataActor = vtkSmartPointer<vtkActor>::New();
    RemeshedTeethPolyDataActor->SetMapper(remeshed_mapper);
    RemeshedTeethPolyDataActor->GetProperty()->SetColor(0, 1, 1);
    RemeshedTeethPolyDataActor->GetProperty()->SetAmbient(0.5);
    RemeshedTeethPolyDataActor->GetProperty()->SetSpecularPower(100);
    RemeshedTeethPolyDataActor->GetProperty()->SetSpecular(0.5);
    RemeshedTeethPolyDataActor->GetProperty()->SetDiffuse(0.5);
    RemeshedTeethPolyDataActor->GetProperty()->SetOpacity(0.8);
    RemeshedTeethPolyDataActor->GetProperty()->EdgeVisibilityOff();
    //RemeshedTeethPolyDataActor->PickableOn();

    std::shared_ptr<SurfaceMesh> adjacent_teeth_sm = teeth_wrapper.GetAdjacentTeethSurfaceMesh();
    vtkSmartPointer<vtkPolyData> adjacent_teeth_polydata = CGAL_Surface_Mesh2VTK_PolyData(*adjacent_teeth_sm);

    vtkSmartPointer<vtkPolyDataMapper> adjacent_teeth_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    adjacent_teeth_mapper->SetInputData(adjacent_teeth_polydata);
    adjacent_teeth_mapper->Update();

    vtkSmartPointer<vtkActor> adjacent_teeth_actor = vtkSmartPointer<vtkActor>::New();
    adjacent_teeth_actor->SetMapper(adjacent_teeth_mapper);
    adjacent_teeth_actor->GetProperty()->SetColor(1, 1, 1);
    adjacent_teeth_actor->GetProperty()->SetAmbient(0.5);
    adjacent_teeth_actor->GetProperty()->SetSpecularPower(100);
    adjacent_teeth_actor->GetProperty()->SetSpecular(0.5);
    adjacent_teeth_actor->GetProperty()->SetDiffuse(0.5);
    adjacent_teeth_actor->GetProperty()->SetOpacity(0.5);
    adjacent_teeth_actor->GetProperty()->EdgeVisibilityOff();
    pipeline->Renderer->AddActor(adjacent_teeth_actor);

    vtkSmartPointer<vtkPolyDataMapper> teeth_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    teeth_mapper->SetInputData(TeethPolydata);
    teeth_mapper->Update();
    TeethPolyDataActor = vtkSmartPointer<vtkActor>::New();
    TeethPolyDataActor->SetMapper(teeth_mapper);
    TeethPolyDataActor->GetProperty()->SetColor(0, 1, 1);
    TeethPolyDataActor->GetProperty()->SetAmbient(0.5);
    TeethPolyDataActor->GetProperty()->SetSpecularPower(100);
    TeethPolyDataActor->GetProperty()->SetSpecular(0.5);
    TeethPolyDataActor->GetProperty()->SetDiffuse(0.5);
    TeethPolyDataActor->GetProperty()->SetOpacity(0.8);
    //TeethPolyDataActor->GetProperty()->EdgeVisibilityOn();
    //TeethPolyDataActor->PickableOn();
    pipeline->Renderer->AddActor(TeethPolyDataActor);

	pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
	pipeline->Renderer->ResetCamera();
	pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
	pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
	pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);

	pipeline->RenderWindowInteractor->Start();
}