// WARNING: Please don't edit this file. It was generated by C++/WinRT v2.0.240405.15

#pragma once
#ifndef WINRT_Windows_UI_Notifications_Management_2_H
#define WINRT_Windows_UI_Notifications_Management_2_H
#include "winrt/impl/Windows.UI.Notifications.Management.1.h"
WINRT_EXPORT namespace winrt::Windows::UI::Notifications::Management
{
    struct WINRT_IMPL_EMPTY_BASES UserNotificationListener : winrt::Windows::UI::Notifications::Management::IUserNotificationListener
    {
        UserNotificationListener(std::nullptr_t) noexcept {}
        UserNotificationListener(void* ptr, take_ownership_from_abi_t) noexcept : winrt::Windows::UI::Notifications::Management::IUserNotificationListener(ptr, take_ownership_from_abi) {}
        [[nodiscard]] static auto Current();
    };
}
#endif
