#include "stdafx.h"
#include "IniIO.h"


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
