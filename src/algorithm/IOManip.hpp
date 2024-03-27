#pragma once
#include "stdafx.h"

int CALLBACK BrowseCallbackProc(HWND hwnd, UINT uMsg, LPARAM lParam, LPARAM lpData)
{
	if (uMsg == BFFM_INITIALIZED)
	{
		// lpData is the lParam value passed to SHBrowseForFolder
		SendMessage(hwnd, BFFM_SETSELECTION, TRUE, lpData);
	}
	return 0;
}

/**
 * @brief Create directories recursively.
 *
 * @param path The full path of the directory to be created.
 *
 * @return TRUE if the directory is successfully created or already exists, FALSE otherwise.
 *
 * @note This function uses Windows API (GetFileAttributesA, CreateDirectoryA) and requires Windows-specific headers.
 *
 * @note The function recursively creates parent directories if they don't exist.
 *
 * @note If the specified path already exists and is a file, the function returns FALSE.
 *
 * @example
 *   std::string path = "C:\\parent\\child1\\child2";
 *   BOOL result = create_directories_recursively(path);
 *   // If result is TRUE, the directory structure is created successfully or already exists.
 */
BOOL create_directories_recursively(const std::string& path)
{
	DWORD dwAttrib = GetFileAttributesA(path.c_str());

	// Check if the path exists and is not a file
	if (dwAttrib != INVALID_FILE_ATTRIBUTES &&
		!(dwAttrib & FILE_ATTRIBUTE_DIRECTORY))
	{
		return FALSE;
	}

	// Try to create the directory
	if (dwAttrib == INVALID_FILE_ATTRIBUTES)
	{
		// Recursively create the parent directory
		size_t slashIndex = path.find_last_of("/\\");
		if (slashIndex != std::string::npos) {
			if (!create_directories_recursively(path.substr(0, slashIndex)))
			{
				return FALSE;
			}
		}

		// Create the last directory
		if (!CreateDirectoryA(path.c_str(), NULL))
		{
			return FALSE;
		}
	}

	return TRUE;
}

bool write_ini_file(const std::string& ini_file_name, const std::string& last_path)
{
	std::ofstream ini_file(ini_file_name);
	if (!ini_file.is_open())
	{
		std::cerr << "Unable to open file for writing: " << ini_file_name << std::endl;
		return false;
	}

	ini_file << "last_path=" << last_path << std::endl;

	ini_file.close();
	return true;
}

bool read_ini_file(const std::string& ini_file_name, std::string& last_path)
{
	std::ifstream ini_file(ini_file_name);
	if (!ini_file.is_open())
	{
		std::cerr << "Unable to open file for reading: " << ini_file_name << std::endl;
		return false;
	}

	std::string line;
	while (std::getline(ini_file, line))
	{
		size_t pos = line.find("last_path=");
		if (pos != std::string::npos) {
			last_path = line.substr(pos + std::string("last_path=").length());
			break;
		}
	}

	ini_file.close();
	return !last_path.empty();
}

/**
 * @brief Open a folder selection dialog and return the selected folder path.
 *
 * This function displays a folder selection dialog using the Windows Shell API.
 * It allows the user to browse and select a folder. The function returns the
 * path of the selected folder as a string. If no folder is selected, an empty
 * string is returned.
 *
 * The function also reads and writes the last selected folder path to an INI file
 * named "last_path.ini". If the file exists, the last selected path is used as the
 * initial directory for the folder selection dialog. After the user selects a folder,
 * the selected path is written back to the INI file.
 *
 * @return The path of the selected folder, or an empty string if no folder is selected.
 *
 * @note This function uses the Windows Shell API (ShBrowseForFolder) to display the
 *       folder selection dialog. It requires the Windows-specific headers and libraries.
 *
 * @note The function uses the C++11 `<codecvt>` library for string conversion between
 *       UTF-8 and wide strings.
 *
 * @note The `BrowseCallbackProc` function is a callback function used by the folder
 *       selection dialog to set the initial directory. It is not shown in this code snippet.
 */
std::string select_folder()
{
	std::string last_path;
	read_ini_file("last_path.ini", last_path);

	BROWSEINFO bi = { 0 };
	bi.lpszTitle = L"Browse for folder...";
	bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
	bi.lpfn = BrowseCallbackProc;

	std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
	std::wstring wide_last_path = converter.from_bytes(last_path);

	bi.lParam = reinterpret_cast<LPARAM>(wide_last_path.c_str());

	LPITEMIDLIST pidl = SHBrowseForFolder(&bi);

	if (pidl != 0)
	{
		wchar_t path[MAX_PATH];
		if (SHGetPathFromIDList(pidl, path))
		{
			IMalloc* imalloc = 0;
			if (SUCCEEDED(SHGetMalloc(&imalloc)))
			{
				imalloc->Free(pidl);
				imalloc->Release();
			}
			std::string selected_folder = converter.to_bytes(path);
			write_ini_file("last_path.ini", selected_folder);
			return selected_folder;
		}
	}

	return std::string();
}
