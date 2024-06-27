#include "MainWindow.h"
#include "simpleRender.h"


class CameraSyncInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static CameraSyncInteractorStyle* New();
    vtkTypeMacro(CameraSyncInteractorStyle, vtkInteractorStyleTrackballCamera);

    void SetOtherRenderer(vtkRenderer* renderer)
    {
        this->OtherRenderer = renderer;
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

private:
    vtkRenderer* OtherRenderer;

    void SyncCameras()
    {
        vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
        vtkCamera* otherCamera = this->OtherRenderer->GetActiveCamera();
        otherCamera->SetPosition(camera->GetPosition());
        otherCamera->SetFocalPoint(camera->GetFocalPoint());
        otherCamera->SetViewUp(camera->GetViewUp());
        otherCamera->SetClippingRange(camera->GetClippingRange());
        otherCamera->SetParallelScale(camera->GetParallelScale());
        this->OtherRenderer->ResetCameraClippingRange();
    }
};

vtkStandardNewMacro(CameraSyncInteractorStyle);

MainWindow::MainWindow()
{
    QWidget* centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);

    QVBoxLayout* mainLayout = new QVBoxLayout(centralWidget);

    QPushButton* openButton = new QPushButton("Open File", centralWidget);
    mainLayout->addWidget(openButton);
    connect(openButton, &QPushButton::clicked, this, &MainWindow::openFile);

    QHBoxLayout* layout = new QHBoxLayout();
    mainLayout->addLayout(layout);

    // 创建第一个 QVTKOpenGLWidget
    leftVtkWidget = new QVTKOpenGLWidget(centralWidget);
    layout->addWidget(leftVtkWidget);

    leftRenderer = vtkRenderer::New();
    leftVtkWidget->GetRenderWindow()->AddRenderer(leftRenderer);
    leftRenderer->SetBackground(0.1, 0.1, 0.1); // RGB color

    // 创建第二个 QVTKOpenGLWidget
    rightVtkWidget = new QVTKOpenGLWidget(centralWidget);
    layout->addWidget(rightVtkWidget);

    rightRenderer = vtkRenderer::New();
    rightVtkWidget->GetRenderWindow()->AddRenderer(rightRenderer);
    rightRenderer->SetBackground(0.1, 0.1, 0.1); // RGB color

    AddAxesToRenderer(leftRenderer);
    AddAxesToRenderer(rightRenderer);

    // 设置自定义交互器
    vtkNew<CameraSyncInteractorStyle> leftInteractorStyle;
    vtkNew<CameraSyncInteractorStyle> rightInteractorStyle;

    leftInteractorStyle->SetDefaultRenderer(leftRenderer);
    rightInteractorStyle->SetDefaultRenderer(rightRenderer);

    leftInteractorStyle->SetOtherRenderer(rightRenderer);
    rightInteractorStyle->SetOtherRenderer(leftRenderer);

    leftVtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(leftInteractorStyle);
    rightVtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(rightInteractorStyle);
}

void MainWindow::AddAxesToRenderer(vtkRenderer* renderer)
{
    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
    axes->SetTotalLength(10.0, 10.0, 10.0); // 设置轴的长度

    vtkSmartPointer<vtkTextProperty> textProperty = vtkSmartPointer<vtkTextProperty>::New();
    textProperty->SetFontSize(3); // 设置字体大小
    
    axes->GetXAxisCaptionActor2D()->GetTextActor()->SetTextProperty(textProperty);
    axes->GetYAxisCaptionActor2D()->GetTextActor()->SetTextProperty(textProperty);
    axes->GetZAxisCaptionActor2D()->GetTextActor()->SetTextProperty(textProperty);

    renderer->AddActor(axes);
}

void MainWindow::openFile()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "", tr("STL Files (*.stl);;All Files (*)"));

    if (!fileName.isEmpty())
    {
        vtkNew<vtkSTLReader> reader;
        reader->SetFileName(fileName.toStdString().c_str());
        reader->Update();

        vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();

        // 计算polydata的质心
        vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter = vtkSmartPointer<vtkCenterOfMass>::New();
        centerOfMassFilter->SetInputData(polydata);
        centerOfMassFilter->Update();

        double center[3];
        centerOfMassFilter->GetCenter(center);

        // 创建一个平移变换，将质心移动到原点
        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform->Translate(-center[0], -center[1], -center[2]);

        // 应用变换到polydata
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter->SetTransform(transform);
        transformFilter->SetInputData(polydata);
        transformFilter->Update();

        // 获取移动后的polydata
        vtkSmartPointer<vtkPolyData> centeredPolydata = transformFilter->GetOutput();

        // 渲染左窗口
        RenderPolydata(centeredPolydata, leftRenderer);

        // 渲染右窗口
        RenderPolydata(centeredPolydata, rightRenderer);

        // 重置摄像机
        leftRenderer->ResetCamera();
        rightRenderer->ResetCamera();

        // 设置左摄像机的位置和方向
        vtkCamera* leftCamera = leftRenderer->GetActiveCamera();
        double leftPosition[3] = { 0.0, 0.0, 100.0 };
        double focalPoint[3] = { 0.0, 0.0, 0.0 };
        double viewUp[3] = { 0.0, 1.0, 0.0 };

        leftCamera->SetPosition(leftPosition);
        leftCamera->SetFocalPoint(focalPoint);
        leftCamera->SetViewUp(viewUp);

        // 设置右摄像机的位置和方向，使其有8度的视角差异
        vtkCamera* rightCamera = rightRenderer->GetActiveCamera();
        double angleRadians = vtkMath::RadiansFromDegrees(8.0);
        double rightPosition[3] = {
            leftPosition[0] * cos(angleRadians) - leftPosition[1] * sin(angleRadians),
            leftPosition[0] * sin(angleRadians) + leftPosition[1] * cos(angleRadians),
            leftPosition[2]
        };

        rightCamera->SetPosition(rightPosition);
        rightCamera->SetFocalPoint(focalPoint);
        rightCamera->SetViewUp(viewUp);

        // 打印相机位置和方向以检查差异
        double leftCameraPosition[3];
        double rightCameraPosition[3];
        leftCamera->GetPosition(leftCameraPosition);
        rightCamera->GetPosition(rightCameraPosition);
        qDebug() << "Left Camera Position:" << leftCameraPosition[0] << leftCameraPosition[1] << leftCameraPosition[2];
        qDebug() << "Right Camera Position:" << rightCameraPosition[0] << rightCameraPosition[1] << rightCameraPosition[2];

        leftVtkWidget->GetRenderWindow()->Render();
        rightVtkWidget->GetRenderWindow()->Render();
    }
}
