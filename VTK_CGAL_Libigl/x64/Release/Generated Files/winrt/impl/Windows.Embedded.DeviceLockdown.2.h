// WARNING: Please don't edit this file. It was generated by C++/WinRT v2.0.240405.15

#pragma once
#ifndef WINRT_Windows_Embedded_DeviceLockdown_2_H
#define WINRT_Windows_Embedded_DeviceLockdown_2_H
#include "winrt/impl/Windows.Embedded.DeviceLockdown.1.h"
WINRT_EXPORT namespace winrt::Windows::Embedded::DeviceLockdown
{
    struct DeviceLockdownProfile
    {
        DeviceLockdownProfile() = delete;
        static auto GetSupportedLockdownProfiles();
        static auto GetCurrentLockdownProfile();
        static auto ApplyLockdownProfileAsync(winrt::guid const& profileID);
        static auto GetLockdownProfileInformation(winrt::guid const& profileID);
    };
    struct WINRT_IMPL_EMPTY_BASES DeviceLockdownProfileInformation : winrt::Windows::Embedded::DeviceLockdown::IDeviceLockdownProfileInformation
    {
        DeviceLockdownProfileInformation(std::nullptr_t) noexcept {}
        DeviceLockdownProfileInformation(void* ptr, take_ownership_from_abi_t) noexcept : winrt::Windows::Embedded::DeviceLockdown::IDeviceLockdownProfileInformation(ptr, take_ownership_from_abi) {}
    };
}
#endif
