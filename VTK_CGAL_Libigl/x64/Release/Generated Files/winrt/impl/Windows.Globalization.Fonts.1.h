// WARNING: Please don't edit this file. It was generated by C++/WinRT v2.0.240405.15

#pragma once
#ifndef WINRT_Windows_Globalization_Fonts_1_H
#define WINRT_Windows_Globalization_Fonts_1_H
#include "winrt/impl/Windows.Globalization.Fonts.0.h"
WINRT_EXPORT namespace winrt::Windows::Globalization::Fonts
{
    struct WINRT_IMPL_EMPTY_BASES ILanguageFont :
        winrt::Windows::Foundation::IInspectable,
        impl::consume_t<ILanguageFont>
    {
        ILanguageFont(std::nullptr_t = nullptr) noexcept {}
        ILanguageFont(void* ptr, take_ownership_from_abi_t) noexcept : winrt::Windows::Foundation::IInspectable(ptr, take_ownership_from_abi) {}
    };
    struct WINRT_IMPL_EMPTY_BASES ILanguageFontGroup :
        winrt::Windows::Foundation::IInspectable,
        impl::consume_t<ILanguageFontGroup>
    {
        ILanguageFontGroup(std::nullptr_t = nullptr) noexcept {}
        ILanguageFontGroup(void* ptr, take_ownership_from_abi_t) noexcept : winrt::Windows::Foundation::IInspectable(ptr, take_ownership_from_abi) {}
    };
    struct WINRT_IMPL_EMPTY_BASES ILanguageFontGroupFactory :
        winrt::Windows::Foundation::IInspectable,
        impl::consume_t<ILanguageFontGroupFactory>
    {
        ILanguageFontGroupFactory(std::nullptr_t = nullptr) noexcept {}
        ILanguageFontGroupFactory(void* ptr, take_ownership_from_abi_t) noexcept : winrt::Windows::Foundation::IInspectable(ptr, take_ownership_from_abi) {}
    };
}
#endif
