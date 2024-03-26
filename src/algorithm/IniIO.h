#pragma once
#include <iostream>
#include <fstream>
#include <string>

bool write_ini_file(const std::string& iniFileName, const std::string& lastPath);
bool read_ini_file(const std::string& iniFileName, std::string& lastPath);
