#include"PolygonIO.h"

PolygonIO::PolygonIO()
{
	b_FileType = false;
	str_Filename = "";
	pointer_PLY_Reader = nullptr;
	pointer_STL_Reader = nullptr;
	pointer_VTP_Reader = nullptr;
	pointer_PLY_Writer = nullptr;
	pointer_STL_Writer = nullptr;
	pointer_VTP_Writer = nullptr;
	pointer_Polydata = nullptr;
}

PolygonIO::~PolygonIO()
{
	// TO DO
	// memory management
}

int PolygonIO::ReadPLY(char* filename)
{
	// Set Filename and initialize Smart Pointer
	str_Filename = filename; 
	if (pointer_PLY_Reader != nullptr)
		pointer_PLY_Reader->Delete();
	pointer_PLY_Reader = vtkSmartPointer<vtkPLYReader>::New();

/*0 - ERROR:please check path, filename and number of points
* 1 - I think I can read the file but I cannot prove it
* 2 - I definitely can read the file
* 3 - I can read the file and I have validated that I am the correct reader for this file*/
	int n_Read_Flag = pointer_PLY_Reader->CanReadFile(filename);
	if (!n_Read_Flag)
	{
		//TO DO:
		// Add log handler
		return 0;
	}

	pointer_PLY_Reader->SetFileName(filename);
	pointer_PLY_Reader->Update();
	if (pointer_Polydata != nullptr)
		pointer_Polydata->Delete();
	pointer_Polydata = vtkSmartPointer<vtkPolyData>::New();
	pointer_Polydata = pointer_PLY_Reader->GetOutput();
	if (DEBUG)
	{
		pointer_PLY_Reader->Print(std::cout);
		pointer_Polydata->Print(std::cout);
	}
	//TO DO:
	// Add log handler
	return n_Read_Flag;
}

int PolygonIO::ReadPLY(std::string str)
{
	str_Filename = str;
	char* filename = (char*)str.data();

	if (pointer_PLY_Reader != nullptr)
		pointer_PLY_Reader->Delete();
	pointer_PLY_Reader = vtkSmartPointer<vtkPLYReader>::New();

	/*1 - I think I can read the file but I cannot prove it
	* 2 - I definitely can read the file
	* 3 - I can read the file and I have validated that I am the correct reader for this file*/
	int Read_Flag = pointer_PLY_Reader->CanReadFile(filename);
	if (!Read_Flag)
	{
		return 0;
	}

	pointer_PLY_Reader->SetFileName(filename);
	pointer_PLY_Reader->Update();
	if (pointer_Polydata != nullptr)
		pointer_Polydata->Delete();
	pointer_Polydata = vtkSmartPointer<vtkPolyData>::New();
	pointer_Polydata = pointer_PLY_Reader->GetOutput();
	if (DEBUG)
	{
		pointer_PLY_Reader->Print(std::cout);
		pointer_Polydata->Print(std::cout);
	}
	return Read_Flag;
}

int PolygonIO::ReadSTL(char* filename)
{
	if (pointer_STL_Reader != nullptr)
		pointer_STL_Reader->Delete();
	pointer_STL_Reader = vtkSmartPointer<vtkSTLReader>::New();

	pointer_STL_Reader->SetFileName(filename);
	pointer_STL_Reader->Update();
	if (pointer_Polydata != nullptr)
		pointer_Polydata->Delete();
	pointer_Polydata = vtkSmartPointer<vtkPolyData>::New();
	pointer_Polydata = pointer_STL_Reader->GetOutput();
	if (DEBUG)
	{
		pointer_STL_Reader->Print(std::cout);
		pointer_Polydata->Print(std::cout);
	}
	if (pointer_Polydata->GetNumberOfPoints() > 0)
		return true;
	else
		return false;
}

int PolygonIO::ReadSTL(std::string str)
{
	str_Filename = str;
	char* filename = (char*)str.data();

	if (pointer_STL_Reader != nullptr)
		pointer_STL_Reader->Delete();
	pointer_STL_Reader = vtkSmartPointer<vtkSTLReader>::New();

	pointer_STL_Reader->SetFileName(filename);
	pointer_STL_Reader->Update();
	if (pointer_Polydata != nullptr)
		pointer_Polydata->Delete();
	pointer_Polydata = vtkSmartPointer<vtkPolyData>::New();
	pointer_Polydata = pointer_STL_Reader->GetOutput();
	if (DEBUG)
	{
		pointer_STL_Reader->Print(std::cout);
		pointer_Polydata->Print(std::cout);
	}
	if (pointer_Polydata->GetNumberOfPoints() > 0)
		return true;
	else 
		return false;
}

int PolygonIO::ReadVTP(char* filename)
{
	str_Filename = filename;
	if (pointer_VTP_Reader != nullptr)
		pointer_VTP_Reader->Delete();
	pointer_VTP_Reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

	/*1 - I think I can read the file but I cannot prove it
	* 2 - I definitely can read the file
	* 3 - I can read the file and I have validated that I am the correct reader for this file*/
	int Read_Flag = pointer_VTP_Reader->CanReadFile(filename);
	if (!Read_Flag)
	{
		return 0;
	}

	pointer_VTP_Reader->SetFileName(filename);
	pointer_VTP_Reader->Update();
	if (pointer_Polydata != nullptr)
		pointer_Polydata->Delete();
	pointer_Polydata = vtkSmartPointer<vtkPolyData>::New();
	pointer_Polydata = pointer_VTP_Reader->GetOutput();
	if (DEBUG)
	{
		pointer_VTP_Reader->Print(std::cout);
		pointer_Polydata->Print(std::cout);
	}
	return Read_Flag;
}

int PolygonIO::ReadVTP(std::string str)
{
	str_Filename = str;
	char* filename = (char*)str.data();
	if (pointer_VTP_Reader != nullptr)
		pointer_VTP_Reader->Delete();
	pointer_VTP_Reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

	/*1 - I think I can read the file but I cannot prove it
	* 2 - I definitely can read the file
	* 3 - I can read the file and I have validated that I am the correct reader for this file*/
	int Read_Flag = pointer_VTP_Reader->CanReadFile(filename);
	if (!Read_Flag)
	{
		return 0;
	}

	pointer_VTP_Reader->SetFileName(filename);
	pointer_VTP_Reader->Update();
	if (pointer_Polydata != nullptr)
		pointer_Polydata->Delete();
	pointer_Polydata = vtkSmartPointer<vtkPolyData>::New();
	pointer_Polydata = pointer_VTP_Reader->GetOutput();
	if (DEBUG)
	{
		pointer_VTP_Reader->Print(std::cout);
		pointer_Polydata->Print(std::cout);
	}
	return Read_Flag;
}

int PolygonIO::Read(char* filename)
{
	str_Filename = filename;
	auto size = str_Filename.size();
	if (str_Filename[size - 4] == '.' && str_Filename[size - 3] == 'p' && str_Filename[size - 2] == 'l' && str_Filename[size - 1] == 'y') // Read ply file.
	{
		return ReadPLY(str_Filename);
	}
	if (str_Filename[size - 4] == '.' && str_Filename[size - 3] == 's' && str_Filename[size - 2] == 't' && str_Filename[size - 1] == 'l') // Read stl file.
	{
		return ReadSTL(str_Filename);
	}
	if (str_Filename[size - 4] == '.' && str_Filename[size - 3] == 'v' && str_Filename[size - 2] == 't' && str_Filename[size - 1] == 'p') // Read vtp file. 
	{
		return ReadVTP(str_Filename);
	}
	return 0;
}

int PolygonIO::Read(std::string str)
{
	auto size = str.size();
	if (str[size - 4] == '.' && str[size - 3] == 'p' && str[size - 2] == 'l' && str[size - 1] == 'y') // Read ply file.
	{
		return ReadPLY(str);
	}
	if (str[size - 4] == '.' && str[size - 3] == 's' && str[size - 2] == 't' && str[size - 1] == 'l') // Read stl file.
	{
		return ReadSTL(str);
	}
	if (str[size - 4] == '.' && str[size - 3] == 'v' && str[size - 2] == 't' && str[size - 1] == 'p') // Read vtp file. 
	{
		return ReadVTP(str);
	}
	return 0;
}

vtkSmartPointer<vtkPolyData> PolygonIO::GetPolyData()
{
	return pointer_Polydata ? pointer_Polydata : nullptr;
}

int PolygonIO::SetPolyData(vtkSmartPointer<vtkPolyData> pd)
{
	if (!pd) return 0;
	pointer_Polydata = pd;
	return 1;
}

void PolygonIO::SetModeToASCII()
{
	b_FileType = false;
}

void PolygonIO::SetModeToBinary()
{
	b_FileType = false;
}

int PolygonIO::WritePLY(char* filename)
{
	str_Filename = filename;
	vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
	if (b_FileType)
		writer->SetFileTypeToASCII();
	else
		writer->SetFileTypeToBinary();
	writer->SetFileName(filename);
	//Return 2 on invalid Polydata
	if (!pointer_Polydata)
	{
		//TO DO:
		//Add log handler 
		return 2;
	}
	writer->SetInputData(pointer_Polydata);
	writer->Update();
	//Returns 1 on success and 0 on failure.  
	return writer->Write(); 
}

int PolygonIO::WritePLY(std::string str)
{
	str_Filename = str;
	char* filename = (char*)str.data();
	vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
	if (b_FileType)
		writer->SetFileTypeToASCII();
	else
		writer->SetFileTypeToBinary();
	writer->SetFileName(filename);
	//Return 2 on invalid Polydata
	if (!pointer_Polydata)
	{
		//TO DO:
		//Add log handler 
		return 2;
	}
	writer->SetInputData(pointer_Polydata);
	writer->Update();
	//Returns 1 on success and 0 on failure.  
	return writer->Write();
}

int PolygonIO::WriteSTL(char* filename)
{
	str_Filename = filename;
	vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
	if (b_FileType)
		writer->SetFileTypeToASCII();
	else
		writer->SetFileTypeToBinary();
	writer->SetFileName(filename);
	//Return 2 on invalid Polydata
	if (!pointer_Polydata)
	{
		//TO DO:
		//Add log handler 
		return 2;
	}
	writer->SetInputData(pointer_Polydata);
	writer->Update();
	//Returns 1 on success and 0 on failure.  
	return writer->Write();
}

int PolygonIO::WriteSTL(std::string str)
{
	str_Filename = str;
	char* filename = (char*)str.data();
	vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
	if (b_FileType)
		writer->SetFileTypeToASCII();
	else
		writer->SetFileTypeToBinary();
	writer->SetFileName(filename);
	//Return 2 on invalid Polydata
	if (!pointer_Polydata)
	{
		//TO DO:
		//Add log handler 
		return 2;
	}
	writer->SetInputData(pointer_Polydata);
	writer->Update();
	//Returns 1 on success and 0 on failure.  
	return writer->Write();
}

int PolygonIO::WriteVTP(char* filename)
{
	str_Filename = filename;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename);
	//Return 2 on invalid Polydata
	if (!pointer_Polydata)
	{
		//TO DO:
		//Add log handler 
		return 2;
	}
	writer->SetInputData(pointer_Polydata);
	writer->Update();
	//Returns 1 on success and 0 on failure.  
	return writer->Write();
}

int PolygonIO::WriteVTP(std::string str)
{

	str_Filename = str;
	char* filename = (char*)str.data();
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename);
	//Return 2 on invalid Polydata
	if (!pointer_Polydata)
	{
		//TO DO:
		//Add log handler 
		return 2;
	}
	writer->SetInputData(pointer_Polydata);
	writer->Update();
	//Returns 1 on success and 0 on failure.  
	return writer->Write();
}


int PolygonIO::Write(char * filename)
{
	str_Filename = filename;

	auto size = str_Filename.size();
	if (str_Filename[size - 4] == '.' && str_Filename[size - 3] == 'p' && str_Filename[size - 2] == 'l' && str_Filename[size - 1] == 'y') // Write ply file.
	{
		return WritePLY(str_Filename);
	}
	if (str_Filename[size - 4] == '.' && str_Filename[size - 3] == 's' && str_Filename[size - 2] == 't' && str_Filename[size - 1] == 'l') // Write stl file.
	{
		return WriteSTL(str_Filename);
	}
	if (str_Filename[size - 4] == '.' && str_Filename[size - 3] == 'v' && str_Filename[size - 2] == 't' && str_Filename[size - 1] == 'p') // Write vtp file. 
	{
		return WriteVTP(str_Filename);
	}
	return 0;
}

int PolygonIO::Write(std::string str)
{
	str_Filename = str;
	auto size = str.size();
	if (str[size - 4] == '.' && str[size - 3] == 'p' && str[size - 2] == 'l' && str[size - 1] == 'y') // Write ply file.
	{
		return WritePLY(str);
	}
	if (str[size - 4] == '.' && str[size - 3] == 's' && str[size - 2] == 't' && str[size - 1] == 'l') // Write stl file.
	{
		return WriteSTL(str);
	}
	if (str[size - 4] == '.' && str[size - 3] == 'v' && str[size - 2] == 't' && str[size - 1] == 'p') // Write vtp file. 
	{
		return WriteVTP(str);
	}
	return 0;
}