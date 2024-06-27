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

    // ������һ�� QVTKOpenGLWidget
    leftVtkWidget = new QVTKOpenGLWidget(centralWidget);
    layout->addWidget(leftVtkWidget);

    leftRenderer = vtkRenderer::New();
    leftVtkWidget->GetRenderWindow()->AddRenderer(leftRenderer);
    leftRenderer->SetBackground(0.1, 0.1, 0.1); // RGB color

    // �����ڶ��� QVTKOpenGLWidget
    rightVtkWidget = new QVTKOpenGLWidget(centralWidget);
    layout->addWidget(rightVtkWidget);

    rightRenderer = vtkRenderer::New();
    rightVtkWidget->GetRenderWindow()->AddRenderer(rightRenderer);
    rightRenderer->SetBackground(0.1, 0.1, 0.1); // RGB color

    AddAxesToRenderer(leftRenderer);
    AddAxesToRenderer(rightRenderer);

    // �����Զ��彻����
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
    axes->SetTotalLength(10.0, 10.0, 10.0); // ������ĳ���

    vtkSmartPointer<vtkTextProperty> textProperty = vtkSmartPointer<vtkTextProperty>::New();
    textProperty->SetFontSize(3); // ���������С
    
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

        // ����polydata������
        vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter = vtkSmartPointer<vtkCenterOfMass>::New();
        centerOfMassFilter->SetInputData(polydata);
        centerOfMassFilter->Update();

        double center[3];
        centerOfMassFilter->GetCenter(center);

        // ����һ��ƽ�Ʊ任���������ƶ���ԭ��
        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform->Translate(-center[0], -center[1], -center[2]);

        // Ӧ�ñ任��polydata
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter->SetTransform(transform);
        transformFilter->SetInputData(polydata);
        transformFilter->Update();

        // ��ȡ�ƶ����polydata
        vtkSmartPointer<vtkPolyData> centeredPolydata = transformFilter->GetOutput();

        // ��Ⱦ�󴰿�
        RenderPolydata(centeredPolydata, leftRenderer);

        // ��Ⱦ�Ҵ���
        RenderPolydata(centeredPolydata, rightRenderer);

        // ���������
        leftRenderer->ResetCamera();
        rightRenderer->ResetCamera();

        // �������������λ�úͷ���
        vtkCamera* leftCamera = leftRenderer->GetActiveCamera();
        double leftPosition[3] = { 0.0, 0.0, 100.0 };
        double focalPoint[3] = { 0.0, 0.0, 0.0 };
        double viewUp[3] = { 0.0, 1.0, 0.0 };

        leftCamera->SetPosition(leftPosition);
        leftCamera->SetFocalPoint(focalPoint);
        leftCamera->SetViewUp(viewUp);

        // �������������λ�úͷ���ʹ����8�ȵ��ӽǲ���
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

        // ��ӡ���λ�úͷ����Լ�����
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
