#pragma once
#include <d3d11.h>
#include <dxgi1_2.h>
#include <wrl/client.h>

// Link with d3d11.lib and dxgi.lib
#pragma comment (lib, "d3d11.lib")
#pragma comment (lib, "dxgi.lib")

using Microsoft::WRL::ComPtr;

class DirectXOutput
{
public:
    DirectXOutput(HWND hwnd, int width, int height);
    void RenderFrame(unsigned char* buffer, int width, int height);

private:
    void InitializeDeviceAndSwapChain(HWND hwnd, int width, int height);
    void UpdateTexture(unsigned char* buffer, int width, int height);

    ComPtr<ID3D11Device> m_device;
    ComPtr<ID3D11DeviceContext> m_deviceContext;
    ComPtr<IDXGISwapChain> m_swapChain;
    ComPtr<ID3D11RenderTargetView> m_renderTargetView;
    ComPtr<ID3D11Texture2D> m_texture;
    ComPtr<ID3D11ShaderResourceView> m_textureView;
};
