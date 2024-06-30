// WARNING: Please don't edit this file. It was generated by C++/WinRT v2.0.240405.15

#pragma once
#ifndef WINRT_Windows_System_RemoteDesktop_Input_0_H
#define WINRT_Windows_System_RemoteDesktop_Input_0_H
WINRT_EXPORT namespace winrt::Windows::System::RemoteDesktop::Input
{
    enum class RemoteKeyEventAttributes : uint32_t
    {
        None = 0,
        KeyUp = 0x1,
        Repeat = 0x2,
        Extended = 0x4,
        Extended1 = 0x8,
    };
    enum class RemoteTextConnectionOptions : uint32_t
    {
        None = 0,
        EnablePredictedKeyReporting = 0x1,
    };
    struct IRemoteTextConnection;
    struct IRemoteTextConnection2;
    struct IRemoteTextConnectionFactory;
    struct IRemoteTextConnectionFactory2;
    struct RemoteTextConnection;
    struct RemoteTextConnectionDataHandler;
}
namespace winrt::impl
{
    template <> struct category<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection>{ using type = interface_category; };
    template <> struct category<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection2>{ using type = interface_category; };
    template <> struct category<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory>{ using type = interface_category; };
    template <> struct category<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory2>{ using type = interface_category; };
    template <> struct category<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnection>{ using type = class_category; };
    template <> struct category<winrt::Windows::System::RemoteDesktop::Input::RemoteKeyEventAttributes>{ using type = enum_category; };
    template <> struct category<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionOptions>{ using type = enum_category; };
    template <> struct category<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionDataHandler>{ using type = delegate_category; };
    template <> inline constexpr auto& name_v<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnection> = L"Windows.System.RemoteDesktop.Input.RemoteTextConnection";
    template <> inline constexpr auto& name_v<winrt::Windows::System::RemoteDesktop::Input::RemoteKeyEventAttributes> = L"Windows.System.RemoteDesktop.Input.RemoteKeyEventAttributes";
    template <> inline constexpr auto& name_v<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionOptions> = L"Windows.System.RemoteDesktop.Input.RemoteTextConnectionOptions";
    template <> inline constexpr auto& name_v<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection> = L"Windows.System.RemoteDesktop.Input.IRemoteTextConnection";
    template <> inline constexpr auto& name_v<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection2> = L"Windows.System.RemoteDesktop.Input.IRemoteTextConnection2";
    template <> inline constexpr auto& name_v<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory> = L"Windows.System.RemoteDesktop.Input.IRemoteTextConnectionFactory";
    template <> inline constexpr auto& name_v<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory2> = L"Windows.System.RemoteDesktop.Input.IRemoteTextConnectionFactory2";
    template <> inline constexpr auto& name_v<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionDataHandler> = L"Windows.System.RemoteDesktop.Input.RemoteTextConnectionDataHandler";
    template <> inline constexpr guid guid_v<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection>{ 0x4E7BB02A,0x183E,0x5E66,{ 0xB5,0xE4,0x3E,0x6E,0x5C,0x57,0x0C,0xF1 } }; // 4E7BB02A-183E-5E66-B5E4-3E6E5C570CF1
    template <> inline constexpr guid guid_v<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection2>{ 0x05F99345,0x84C8,0x56C5,{ 0x93,0x4F,0x73,0xEA,0x00,0xF8,0xC2,0xD5 } }; // 05F99345-84C8-56C5-934F-73EA00F8C2D5
    template <> inline constexpr guid guid_v<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory>{ 0x88E075C2,0x0CAE,0x596C,{ 0x85,0x0F,0x78,0xD3,0x45,0xCD,0x72,0x8B } }; // 88E075C2-0CAE-596C-850F-78D345CD728B
    template <> inline constexpr guid guid_v<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory2>{ 0x9425C7D9,0xED9B,0x5D00,{ 0x99,0xCC,0xB0,0xB8,0xDC,0x9E,0x4C,0x60 } }; // 9425C7D9-ED9B-5D00-99CC-B0B8DC9E4C60
    template <> inline constexpr guid guid_v<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionDataHandler>{ 0x099FFBC8,0x8BCB,0x41B5,{ 0xB0,0x56,0x57,0xE7,0x70,0x21,0xBF,0x1B } }; // 099FFBC8-8BCB-41B5-B056-57E77021BF1B
    template <> struct default_interface<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnection>{ using type = winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection; };
    template <> struct abi<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection>
    {
        struct WINRT_IMPL_NOVTABLE type : inspectable_abi
        {
            virtual int32_t __stdcall get_IsEnabled(bool*) noexcept = 0;
            virtual int32_t __stdcall put_IsEnabled(bool) noexcept = 0;
            virtual int32_t __stdcall RegisterThread(uint32_t) noexcept = 0;
            virtual int32_t __stdcall UnregisterThread(uint32_t) noexcept = 0;
            virtual int32_t __stdcall ReportDataReceived(uint32_t, uint8_t*) noexcept = 0;
        };
    };
    template <> struct abi<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection2>
    {
        struct WINRT_IMPL_NOVTABLE type : inspectable_abi
        {
            virtual int32_t __stdcall ReportPredictedKeyEvent(uint16_t, uint32_t) noexcept = 0;
        };
    };
    template <> struct abi<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory>
    {
        struct WINRT_IMPL_NOVTABLE type : inspectable_abi
        {
            virtual int32_t __stdcall CreateInstance(winrt::guid, void*, void**) noexcept = 0;
        };
    };
    template <> struct abi<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory2>
    {
        struct WINRT_IMPL_NOVTABLE type : inspectable_abi
        {
            virtual int32_t __stdcall CreateInstance(winrt::guid, void*, uint32_t, void**) noexcept = 0;
        };
    };
    template <> struct abi<winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionDataHandler>
    {
        struct WINRT_IMPL_NOVTABLE type : unknown_abi
        {
            virtual int32_t __stdcall Invoke(uint32_t, uint8_t*, bool*) noexcept = 0;
        };
    };
    template <typename D>
    struct consume_Windows_System_RemoteDesktop_Input_IRemoteTextConnection
    {
        [[nodiscard]] auto IsEnabled() const;
        auto IsEnabled(bool value) const;
        auto RegisterThread(uint32_t threadId) const;
        auto UnregisterThread(uint32_t threadId) const;
        auto ReportDataReceived(array_view<uint8_t const> pduData) const;
    };
    template <> struct consume<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection>
    {
        template <typename D> using type = consume_Windows_System_RemoteDesktop_Input_IRemoteTextConnection<D>;
    };
    template <typename D>
    struct consume_Windows_System_RemoteDesktop_Input_IRemoteTextConnection2
    {
        auto ReportPredictedKeyEvent(uint16_t scanCode, winrt::Windows::System::RemoteDesktop::Input::RemoteKeyEventAttributes const& attributes) const;
    };
    template <> struct consume<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnection2>
    {
        template <typename D> using type = consume_Windows_System_RemoteDesktop_Input_IRemoteTextConnection2<D>;
    };
    template <typename D>
    struct consume_Windows_System_RemoteDesktop_Input_IRemoteTextConnectionFactory
    {
        auto CreateInstance(winrt::guid const& connectionId, winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionDataHandler const& pduForwarder) const;
    };
    template <> struct consume<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory>
    {
        template <typename D> using type = consume_Windows_System_RemoteDesktop_Input_IRemoteTextConnectionFactory<D>;
    };
    template <typename D>
    struct consume_Windows_System_RemoteDesktop_Input_IRemoteTextConnectionFactory2
    {
        auto CreateInstance(winrt::guid const& connectionId, winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionDataHandler const& pduForwarder, winrt::Windows::System::RemoteDesktop::Input::RemoteTextConnectionOptions const& options) const;
    };
    template <> struct consume<winrt::Windows::System::RemoteDesktop::Input::IRemoteTextConnectionFactory2>
    {
        template <typename D> using type = consume_Windows_System_RemoteDesktop_Input_IRemoteTextConnectionFactory2<D>;
    };
}
#endif
