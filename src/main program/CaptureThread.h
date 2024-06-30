#pragma once
#ifndef CAPTURETHREAD_H
#define CAPTURETHREAD_H

#include <QThread>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <string>
#include <sstream>
#include <iomanip>

class CaptureThread : public QThread
{
    Q_OBJECT

public:
    CaptureThread(int numFrames, QObject* parent = nullptr)
        : QThread(parent), m_numFrames(numFrames) {}

signals:
    void requestCapture(int frameNumber);

protected:
    void run() override
    {
        for (int i = 0; i < m_numFrames; ++i) {
            emit requestCapture(i);
            msleep(100); // Optional: Adjust sleep time as needed
        }
    }

private:
    int m_numFrames;
};

#endif // CAPTURETHREAD_H