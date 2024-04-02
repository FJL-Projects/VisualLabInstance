#pragma once
#include<stdafx.h>
#include <cstring>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include<vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkXMLPolyDataReader.h>
#include<vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include<iostream>
#include<vtkIOStream.h>
#include <vtkDataSetMapper.h>
#define DEBUG 1
// 20220807 add stl/ply/vtp Reader and Writer by Zhang Xukun
// Current stage: Not tested
// TO DO:
// Add log handler
// Check Polygon on write function
// Modify Header on write
// Read and Write polygon with face normal, vertex normal, face color(rgb format, rgba format and hex), point color, face texture, customize feature points and customize labels(points, edges, faces) 
// memory management
// Mesh management(adjacency point, edge, face list) 
// Mesh simplify and densification
// Polygon and Mesh hole repair
// Polygon and Mesh bubble repair

class PolygonIO
{
public:
	PolygonIO();

	~PolygonIO();
	// Method: 
	// PolygonIO Reader;
	// reader.ReadPLY('Filename.ply');
	//
	// string Filename = "filename"; 
	// Filename += ".vtp";
	// Read(Filename);
	// 
	// return 0: Error : File can't find or polygon contains zero point. May not support Korean(need to fix).
	// return 1 - I think I can read the file but I cannot prove it
	// return 2 - I definitely can read the file
	// return 3 - I can read the file and I have validated that I am the correct reader for this file
	//////////////////////////////////////////////////////////////////////////////////// Read File
	int ReadPLY(char*);

	int ReadPLY(std::string);

	int ReadSTL(char*);

	int ReadSTL(std::string);

	int ReadVTP(char*);

	int ReadVTP(std::string);

	int Read(char*);

	int Read(std::string);
	//////////////////////////////////////////////////////////////////////////////////// Read File


	//////////////////////////////////////////////////////////////////////////////////// Write File 
	int WritePLY(char*);

	int WritePLY(std::string);

	int WriteSTL(char*);

	int WriteSTL(std::string);

	int WriteVTP(char*);

	int WriteVTP(std::string);

	int Write(char*);

	int Write(std::string);
	//////////////////////////////////////////////////////////////////////////////////// Write File

	//////////////////////////////////////////////////////////////////////////////////// Set Polydata and Get Polydata
	vtkSmartPointer<vtkPolyData> GetPolyData();

	int SetPolyData(vtkSmartPointer<vtkPolyData>);
	//////////////////////////////////////////////////////////////////////////////////// Set Polydata and Get Polydata

	//////////////////////////////////////////////////////////////////////////////////// Write Method
	void SetModeToASCII();
	void SetModeToBinary();
	//////////////////////////////////////////////////////////////////////////////////// Write Method

	bool b_FileType; //Binary(Defult):0  ASCII:1

	std::string str_Filename; // Need to be set before Write

	vtkSmartPointer<vtkPLYReader> pointer_PLY_Reader;
	vtkSmartPointer<vtkPLYWriter> pointer_PLY_Writer;
	vtkSmartPointer<vtkSTLReader> pointer_STL_Reader;
	vtkSmartPointer<vtkSTLWriter> pointer_STL_Writer;
	vtkSmartPointer<vtkXMLPolyDataReader> pointer_VTP_Reader;
	vtkSmartPointer<vtkXMLPolyDataWriter> pointer_VTP_Writer;
	vtkSmartPointer<vtkPolyData> pointer_Polydata; // Need to be set before Write


private:

};





