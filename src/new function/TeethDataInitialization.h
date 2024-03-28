#ifndef TEETHDATAINITIALIZATION_H
#define TEETHDATAINITIALIZATION_H
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProperty.h>
#include <vtkOBBTree.h>
#include "vectorAlgorithm.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Point_3.h>
#include <queue>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

typedef CGAL::Simple_cartesian<double>          Kernel;
typedef Kernel::Point_3                         Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>     SurfaceMesh;
/**
 * Represents a TeethDataInitialization class.
 *
 * This class provides functionality for initializing teeth data, performing various operations on teeth, and rendering them.
 */

class TeethDataInitialization
{
public:
    TeethDataInitialization() = default;
    /**
      * Constructs a TeethDataInitialization object with the given parameters.
      *
      * @param pd A reference to a vtkPolyData object representing teeth data.
      * @param sm A reference to a SurfaceMesh object representing the teeth surface mesh.
      */
    TeethDataInitialization(vtkSmartPointer<vtkPolyData>& pd, SurfaceMesh& sm)
    {
        m_teethPolyData = pd;
        m_teethSurfaceMesh = &sm;
        m_labeledPolyData.resize(20);
        m_labeledPoints.resize(20);
        m_labeledTriangles.resize(20);
        m_labeledActors.resize(20);
        for (int i = 0; i < 20; ++i)
        {
            m_labeledPolyData[i] = vtkSmartPointer<vtkPolyData>::New();
            m_labeledPoints[i] = vtkSmartPointer<vtkPoints>::New();
            m_labeledTriangles[i] = vtkSmartPointer<vtkCellArray>::New();
            m_labeledActors[i]= vtkSmartPointer<vtkActor>::New();
        }
        m_splitPolyDataFlag.resize(20, 0);
        m_splitPolydataVisit.resize(20, std::vector<int>(m_teethPolyData->GetPoints()->GetNumberOfPoints(), -1));
        
        m_teeth_center = vtkSmartPointer<vtkPoints>::New();
    }

    void FragmentClearing();
    void GetLabel(); // initpolydata
    void MakeLabeledActor(); // RenderTeeth
    void ObbBoundingBox();
    void DrawTeethCenterSpline();
    void CalculateDirection();
    double* GetTongueCheekDirection(int n);
    int GetPointOnNthTeeth(int);
    void Execute();
    void RemoveAllSplineActors();


    void SetCheckId(int id)
    {
        m_checked_id = id;
    }

    void SetRenderer(vtkSmartPointer<vtkRenderer> renderer)
    {
        m_renderer = renderer;
    }

    vtkSmartPointer<vtkRenderer> GetRenderer() const
    {
        return m_renderer;
    }

    void SetRenderWindow(vtkSmartPointer<vtkRenderWindow> render_window)
    {
        m_renderWindow = render_window;
    }

    vtkSmartPointer<vtkRenderWindow> GetRenderWindow() const
    {
        return m_renderWindow;
    }

    void SetLabeledPolydata(const std::vector<vtkSmartPointer<vtkPolyData>>& labeledPolydata) 
    {
        m_labeledPolyData = labeledPolydata;
    }

    void SetLabeledActors(const std::vector<vtkSmartPointer<vtkActor>>& labeledActor) 
    {
        m_labeledActors = labeledActor;
    }

    const std::vector<vtkSmartPointer<vtkActor>>& GetLabeledActors() const 
    {
        return m_labeledActors;
    }

    void SetTeethCenter(vtkSmartPointer<vtkPoints> teeth_center)
    {
        m_teeth_center = teeth_center;
    }

    vtkSmartPointer<vtkPoints> GetTeethCenter() const 
    {
        return m_teeth_center;
    }

    const std::vector<vtkSmartPointer<vtkPolyData>>& GetLabeledPolyData() const
    {
        return m_labeledPolyData;
    }

    void SetLabeledPolyData(const std::vector<vtkSmartPointer<vtkPolyData>>& labeledPolyData) 
    {
        m_labeledPolyData = labeledPolyData;
    }

    double3 GetAverageOcclusal() const 
    {
        return m_average_occlusal;
    }

    void SetAverageOcclusal(const double3& value) 
    {
        m_average_occlusal = value;
    }

    // Getter
    const std::vector<double3> GetOcclusalDirection() const 
    {
        return m_occlusal_direction;
    }

    const std::vector<int>& GetValidTeethIndices() const 
	{
		return m_valid_teeth_indices;
	}
    
    const int& GetValidTeethNum() const
    {
        return m_valid_teeth_num;
    }

    vtkSmartPointer<vtkPolyData> GetTeethPolyData() const
    {
        return m_teethPolyData;
    }

public:
    vtkSmartPointer<vtkPolyData> m_teethPolyData;
    SurfaceMesh *m_teethSurfaceMesh;
    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkRenderWindow> m_renderWindow;
    std::vector<vtkSmartPointer<vtkPolyData>> m_labeledPolyData;
    std::vector<vtkSmartPointer<vtkPoints>> m_labeledPoints;
    std::vector<vtkSmartPointer<vtkCellArray>> m_labeledTriangles;
    std::vector<int> m_splitPolyDataFlag;
    std::vector<std::vector<int>> m_splitPolydataVisit;
    std::vector<vtkSmartPointer<vtkActor>> m_labeledActors;
    vtkSmartPointer<vtkPoints> m_teeth_center;
    int m_teethNum;
    std::vector<int> m_missingTeeth;
    std::vector<double3> m_occlusal_direction;
    double3 m_average_occlusal = { 0, 0, 0 };
    int m_checked_id = 1;
    std::vector<double3> m_teeth_arch;

    // All the actors in the center line, which is for later deletion.
    std::vector<vtkSmartPointer<vtkActor> > m_spline_actor_vec;
    std::vector<int> m_valid_teeth_indices;
    int m_valid_teeth_num = 0;
};

#endif
