#include <QApplication>
#include <QVTKWidget.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkCubeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>

#include "DirectXOutput.h"

DirectXOutput::DirectXOutput(HWND hwnd, int width, int height)
{
    InitializeDeviceAndSwapChain(hwnd, width, height);
}

void DirectXOutput::InitializeDeviceAndSwapChain(HWND hwnd, int width, int height)
{
    DXGI_SWAP_CHAIN_DESC swapChainDesc = {};
    swapChainDesc.BufferCount = 1;
    swapChainDesc.BufferDesc.Width = width;
    swapChainDesc.BufferDesc.Height = height;
    swapChainDesc.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    swapChainDesc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    swapChainDesc.OutputWindow = hwnd;
    swapChainDesc.SampleDesc.Count = 1;
    swapChainDesc.Windowed = TRUE;

    D3D11CreateDeviceAndSwapChain(
        nullptr,
        D3D_DRIVER_TYPE_HARDWARE,
        nullptr,
        0,
        nullptr,
        0,
        D3D11_SDK_VERSION,
        &swapChainDesc,
        &m_swapChain,
        &m_device,
        nullptr,
        &m_deviceContext);

    ComPtr<ID3D11Texture2D> backBuffer;
    m_swapChain->GetBuffer(0, IID_PPV_ARGS(&backBuffer));
    m_device->CreateRenderTargetView(backBuffer.Get(), nullptr, &m_renderTargetView);

    // ´´½¨ÎÆÀí
    D3D11_TEXTURE2D_DESC textureDesc = {};
    textureDesc.Width = width;
    textureDesc.Height = height;
    textureDesc.MipLevels = 1;
    textureDesc.ArraySize = 1;
    textureDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    textureDesc.SampleDesc.Count = 1;
    textureDesc.Usage = D3D11_USAGE_DYNAMIC;
    textureDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
    textureDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

    m_device->CreateTexture2D(&textureDesc, nullptr, &m_texture);

    D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
    srvDesc.Format = textureDesc.Format;
    srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
    srvDesc.Texture2D.MostDetailedMip = 0;
    srvDesc.Texture2D.MipLevels = 1;

    m_device->CreateShaderResourceView(m_texture.Get(), &srvDesc, &m_textureView);
}

void DirectXOutput::UpdateTexture(unsigned char* buffer, int width, int height)
{
    D3D11_MAPPED_SUBRESOURCE mappedResource;
    m_deviceContext->Map(m_texture.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedResource);

    unsigned char* textureData = reinterpret_cast<unsigned char*>(mappedResource.pData);
    for (int y = 0; y < height; ++y)
    {
        memcpy(textureData + y * mappedResource.RowPitch, buffer + (height - 1 - y) * width * 4, width * 4);
    }

    m_deviceContext->Unmap(m_texture.Get(), 0);
}

void DirectXOutput::RenderFrame(unsigned char* buffer, int width, int height)
{
    UpdateTexture(buffer, width, height);

    // Bind the render target view
    m_deviceContext->OMSetRenderTargets(1, m_renderTargetView.GetAddressOf(), nullptr);

    // Set up the viewport
    D3D11_VIEWPORT viewport = {};
    viewport.TopLeftX = 0;
    viewport.TopLeftY = 0;
    viewport.Width = static_cast<float>(width);
    viewport.Height = static_cast<float>(height);
    m_deviceContext->RSSetViewports(1, &viewport);

    // Clear the render target
    const float clearColor[4] = { 0.0f, 0.2f, 0.4f, 1.0f };
    m_deviceContext->ClearRenderTargetView(m_renderTargetView.Get(), clearColor);

    // Render the texture
    ComPtr<ID3D11ShaderResourceView> srvs[] = { m_textureView };
    m_deviceContext->PSSetShaderResources(0, 1, srvs->GetAddressOf());

    // Here you would typically set up your shaders and draw calls
    // This is a simplified example assuming you have a shader that draws a full-screen quad

    // Present the frame
    m_swapChain->Present(1, 0);
}