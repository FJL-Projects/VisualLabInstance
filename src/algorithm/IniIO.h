#pragma once
#include <iostream>
#include <fstream>
#include <string>

/**
 * Interacts with an INI file to read and write the last path selected by the user.
 */
bool write_ini_file(const std::string& iniFileName, const std::string& lastPath);
bool read_ini_file(const std::string& iniFileName, std::string& lastPath);
