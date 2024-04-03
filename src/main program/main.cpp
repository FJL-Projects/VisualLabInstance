#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"

vtkRenderPipeline* pipeline;

/**
 * @brief Generate a 4-digit number string with leading zeros.
 *
 * This function takes an integer and converts it into a string representation
 * with a fixed width of 4 characters. If the number has less than 4 digits,
 * leading zeros are added to pad the string to the desired width.
 *
 * @param number The input integer to be converted.
 * @return A string representation of the input number with leading zeros.
 *
 * @note The function uses std::ostringstream, std::setw(), and std::setfill()
 *       to format the output string.
 *
 * @example
 *   int num = 42;
 *   std::string num_str = generate_leading_zero_number_str(num);
 *   // num_str will be "0042"
 */
std::string generate_leading_zero_number_str(int number)
{
	std::ostringstream stream;
	stream << std::setw(4) << std::setfill('0') << number;
	return stream.str();
}

void RightPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Press" << endl;
}
void RightRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Released" << endl;
}

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Left Press" << endl;
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	
}

void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}

#include "stdafx.h"
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkGlyph3DMapper.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkLine.h>
#include <vtkMatrix4x4.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include<vtkLandmarkTransform.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataReader.h>
#include<vtkAppendPolyData.h>
//#include "Render.h"
#include"vectorAlgorithm.h"
#include"TeethDataInitialization.h"
//#include"Rotation.h"

#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/bounding_box.h>
#include"PolygonIO.h"
#include <codecvt>
typedef std::vector<K::Point_3>                               Polyline_type;
typedef std::vector<Polyline_type>                              Polylines;
typedef CGAL::AABB_halfedge_graph_segment_primitive<SurfaceMesh>     HGSP;
typedef CGAL::AABB_traits<K, HGSP>                            AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>                          AABB_tree;
typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;


std::string wide_string2utf8_string(std::wstring origin)
{
	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;

	return converter.to_bytes(origin);
}

int selected_id = 6;
int checked_id = 4;
std::wstring w_absolute_working_directory;
vtkSmartPointer<vtkPolyData> Crownpolydata = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyData> PolyData = vtkSmartPointer<vtkPolyData>::New();
SurfaceMesh mesh;



int main()
{
	std::map<std::string, std::pair<int, int> > file_name_map;  // file_name, checked_id, selected_id
	//file_name_map["1781.vtp"] = std::make_pair(1, 6);
	file_name_map["11737.vtp"] = std::make_pair(2, 7);
	//file_name_map["13067.vtp"] = std::make_pair(4, 6);
	file_name_map["33923.vtp"] = std::make_pair(1, 7);


	for (auto& files : file_name_map)
	{
		std::string file_path = "data\\" + files.first;


		checked_id = files.second.first;
		selected_id = files.second.second;
		pipeline = new vtkRenderPipeline();


		vtkSmartPointer<vtkXMLPolyDataReader> vtp_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
		vtp_reader->SetFileName(file_path.c_str());
		vtp_reader->Update();

		PolyData = vtp_reader->GetOutput();
		mesh = VTK_PolyData2CGAL_Surface_Mesh(PolyData);

		TeethDataInitialization teeth_data_initializer;
		teeth_data_initializer = TeethDataInitialization(PolyData, mesh);
		teeth_data_initializer.SetRenderer(pipeline->Renderer);
		teeth_data_initializer.SetRenderWindow(pipeline->RenderWindow);
		teeth_data_initializer.Execute();

		//crown Polydata
		std::wstring crown_name = w_absolute_working_directory + L"Teeth\\10.stl";
		crown_name[crown_name.length() - 6] = checked_id + '0';
		crown_name[crown_name.length() - 5] = static_cast<int>(selected_id % 10) + '0';
		PolygonIO reader3;
		reader3.Read(wide_string2utf8_string(crown_name));
		Crownpolydata = reader3.GetPolyData();
		std::cout << "Crownpolydata->GetNumberOfPoints():" << Crownpolydata->GetNumberOfPoints() << std::endl;
		//crown Polydata

		//source Polydata
		double3 SourcePre2Next;
		std::wstring pre_name = w_absolute_working_directory + L"Teeth\\Clip\\10.stl";
		pre_name[pre_name.length() - 6] = checked_id + '0';
		if (selected_id == 1)
			pre_name[pre_name.length() - 5] = static_cast<int>(selected_id % 10 + 1) + '0';
		else if (selected_id == 7)
			pre_name[pre_name.length() - 5] = static_cast<int>(selected_id % 10 - 2) + '0';
		else
			pre_name[pre_name.length() - 5] = static_cast<int>(selected_id % 10 - 1) + '0';

		PolygonIO pre_reader;
		pre_reader.Read(wide_string2utf8_string(pre_name));
		vtkSmartPointer<vtkPolyData> pre_polydata = vtkSmartPointer<vtkPolyData>::New();
		pre_polydata = pre_reader.GetPolyData();
		std::cout << "pre_polydata->GetNumberOfPoints():" << pre_polydata->GetNumberOfPoints() << std::endl;

		std::wstring next_name = w_absolute_working_directory + L"Teeth\\Clip\\10.stl";
		next_name[next_name.length() - 6] = checked_id + '0';
		if (selected_id == 1)
			next_name[next_name.length() - 5] = static_cast<int>(selected_id % 10 + 2) + '0';
		else if (selected_id == 7)
			next_name[next_name.length() - 5] = static_cast<int>(selected_id % 10 - 1) + '0';
		else
			next_name[next_name.length() - 5] = static_cast<int>(selected_id % 10 + 1) + '0';
		PolygonIO next_reader;
		next_reader.Read(wide_string2utf8_string(next_name));
		vtkSmartPointer<vtkPolyData> next_polydata = vtkSmartPointer<vtkPolyData>::New();
		next_polydata = next_reader.GetPolyData();
		std::cout << "next_polydata->GetNumberOfPoints():" << next_polydata->GetNumberOfPoints() << std::endl;



		vtkSmartPointer<vtkAppendPolyData> SourceAppendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
		SourceAppendFilter->AddInputData(pre_polydata);
		SourceAppendFilter->AddInputData(next_polydata);
		SourceAppendFilter->Update();

		vtkSmartPointer<vtkPolyData> source_polydata = vtkSmartPointer<vtkPolyData>::New();
		source_polydata = SourceAppendFilter->GetOutput();
		SourcePre2Next.data[0] = next_polydata->GetCenter()[0] - pre_polydata->GetCenter()[0];
		SourcePre2Next.data[1] = next_polydata->GetCenter()[1] - pre_polydata->GetCenter()[1];
		SourcePre2Next.data[2] = next_polydata->GetCenter()[2] - pre_polydata->GetCenter()[2];

		//source Polydata

		//target Polydata
		double3 TargetPre2Next;
		vtkSmartPointer<vtkAppendPolyData> TargetAppendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
		if (selected_id == 1)
		{
			TargetAppendFilter->AddInputData(teeth_data_initializer.m_labeledPolyData[selected_id + 1]);
			TargetAppendFilter->AddInputData(teeth_data_initializer.m_labeledPolyData[selected_id + 2]);
			TargetPre2Next.data[0] = teeth_data_initializer.m_labeledPolyData[selected_id + 2]->GetCenter()[0] - teeth_data_initializer.m_labeledPolyData[selected_id + 1]->GetCenter()[0];
			TargetPre2Next.data[1] = teeth_data_initializer.m_labeledPolyData[selected_id + 2]->GetCenter()[1] - teeth_data_initializer.m_labeledPolyData[selected_id + 1]->GetCenter()[1];
			TargetPre2Next.data[2] = teeth_data_initializer.m_labeledPolyData[selected_id + 2]->GetCenter()[2] - teeth_data_initializer.m_labeledPolyData[selected_id + 1]->GetCenter()[2];
		}
		else if (selected_id == 7)
		{
			TargetAppendFilter->AddInputData(teeth_data_initializer.m_labeledPolyData[selected_id - 1]);
			TargetAppendFilter->AddInputData(teeth_data_initializer.m_labeledPolyData[selected_id - 2]);
			TargetPre2Next.data[0] = teeth_data_initializer.m_labeledPolyData[selected_id - 1]->GetCenter()[0] - teeth_data_initializer.m_labeledPolyData[selected_id - 2]->GetCenter()[0];
			TargetPre2Next.data[1] = teeth_data_initializer.m_labeledPolyData[selected_id - 1]->GetCenter()[1] - teeth_data_initializer.m_labeledPolyData[selected_id - 2]->GetCenter()[1];
			TargetPre2Next.data[2] = teeth_data_initializer.m_labeledPolyData[selected_id - 1]->GetCenter()[2] - teeth_data_initializer.m_labeledPolyData[selected_id - 2]->GetCenter()[2];
		}
		else
		{
			TargetAppendFilter->AddInputData(teeth_data_initializer.m_labeledPolyData[selected_id - 1]);
			TargetAppendFilter->AddInputData(teeth_data_initializer.m_labeledPolyData[selected_id + 1]);
			TargetPre2Next.data[0] = teeth_data_initializer.m_labeledPolyData[selected_id + 1]->GetCenter()[0] - teeth_data_initializer.m_labeledPolyData[selected_id - 1]->GetCenter()[0];
			TargetPre2Next.data[1] = teeth_data_initializer.m_labeledPolyData[selected_id + 1]->GetCenter()[1] - teeth_data_initializer.m_labeledPolyData[selected_id - 1]->GetCenter()[1];
			TargetPre2Next.data[2] = teeth_data_initializer.m_labeledPolyData[selected_id + 1]->GetCenter()[2] - teeth_data_initializer.m_labeledPolyData[selected_id - 1]->GetCenter()[2];
		}
		TargetAppendFilter->Update();

		vtkSmartPointer<vtkPolyData> target_polydata = vtkSmartPointer<vtkPolyData>::New();
		target_polydata = TargetAppendFilter->GetOutput();
		//target Polydata

		//中心对其+法向量对其测试
		double3 source_center(source_polydata->GetCenter());
		double3 target_center(target_polydata->GetCenter());
		for (int i = 0; i < source_polydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double point[3];
			source_polydata->GetPoints()->GetPoint(i, point);
			point[0] = point[0] - source_center[0];
			point[1] = point[1] - source_center[1];
			point[2] = point[2] - source_center[2];
			source_polydata->GetPoints()->SetPoint(i, point);
		}
		source_polydata->Modified();

		for (int i = 0; i < Crownpolydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double point[3];
			Crownpolydata->GetPoints()->GetPoint(i, point);
			point[0] = point[0] - source_center[0];
			point[1] = point[1] - source_center[1];
			point[2] = point[2] - source_center[2];
			Crownpolydata->GetPoints()->SetPoint(i, point);
		}
		Crownpolydata->Modified();

		for (int i = 0; i < target_polydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double point[3];
			target_polydata->GetPoints()->GetPoint(i, point);
			point[0] = point[0] - target_center[0];
			point[1] = point[1] - target_center[1];
			point[2] = point[2] - target_center[2];
			target_polydata->GetPoints()->SetPoint(i, point);
		}
		target_polydata->Modified();

		SourcePre2Next.normalize();
		TargetPre2Next.normalize();

		double3 source_normal;
		double3 target_normal;
		for (int i = 0; i < source_polydata->GetNumberOfCells(); i++)
		{
			double pt0[3];
			source_polydata->GetCell(i)->GetPoints()->GetPoint(0, pt0);
			double pt1[3];
			source_polydata->GetCell(i)->GetPoints()->GetPoint(1, pt1);
			double pt2[3];
			source_polydata->GetCell(i)->GetPoints()->GetPoint(2, pt2);
			double3 vec1(pt2[0] - pt0[0], pt2[1] - pt0[1], pt2[2] - pt0[2]);
			double3 vec2(pt1[0] - pt0[0], pt1[1] - pt0[1], pt1[2] - pt0[2]);
			double3 normal = double3::crossProduct(vec1, vec2);
			normal.normalize();
			source_normal += normal;
		}
		source_normal = source_normal / static_cast<double>(source_polydata->GetNumberOfCells());

		for (int i = 0; i < target_polydata->GetNumberOfCells(); i++)
		{
			double pt0[3];
			target_polydata->GetCell(i)->GetPoints()->GetPoint(0, pt0);
			double pt1[3];
			target_polydata->GetCell(i)->GetPoints()->GetPoint(1, pt1);
			double pt2[3];
			target_polydata->GetCell(i)->GetPoints()->GetPoint(2, pt2);
			double3 vec1(pt2[0] - pt0[0], pt2[1] - pt0[1], pt2[2] - pt0[2]);
			double3 vec2(pt1[0] - pt0[0], pt1[1] - pt0[1], pt1[2] - pt0[2]);
			double3 normal = double3::crossProduct(vec1, vec2);
			normal.normalize();
			target_normal += normal;
		}
		target_normal = target_normal / static_cast<double>(target_polydata->GetNumberOfCells());

		//if (checked_id == 1 || checked_id == 2)
		//	source_normal = double3(0, 0, -1);

		//if (checked_id == 3 || checked_id == 4)
		//	source_normal = double3(0, 0, 1);
		//target_normal = teeth_data_initializer.GetAverageOcclusal();
		target_normal.normalize();
		double3 axis = double3::crossProduct(source_normal, target_normal);
		axis.normalize();
		axis = axis * acos(double3::dotProduct(source_normal, target_normal));
		for (int i = 0; i < source_polydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double pt[3];
			source_polydata->GetPoints()->GetPoint(i, pt);
			AngleAxisRotatePoint(axis.data, pt, pt);
			source_polydata->GetPoints()->SetPoint(i, pt);
		}
		source_polydata->Modified();
		for (int i = 0; i < Crownpolydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double pt[3];
			Crownpolydata->GetPoints()->GetPoint(i, pt);
			AngleAxisRotatePoint(axis.data, pt, pt);
			Crownpolydata->GetPoints()->SetPoint(i, pt);
		}
		Crownpolydata->Modified();
		AngleAxisRotatePoint(axis.data, SourcePre2Next.data, SourcePre2Next.data);


		TargetPre2Next.normalize();
		axis = double3::crossProduct(SourcePre2Next, TargetPre2Next);
		axis.normalize();
		axis = axis * acos(double3::dotProduct(SourcePre2Next, TargetPre2Next));
		for (int i = 0; i < source_polydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double pt[3];
			source_polydata->GetPoints()->GetPoint(i, pt);
			AngleAxisRotatePoint(axis.data, pt, pt);
			source_polydata->GetPoints()->SetPoint(i, pt);
		}
		source_polydata->Modified();
		for (int i = 0; i < Crownpolydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double pt[3];
			Crownpolydata->GetPoints()->GetPoint(i, pt);
			AngleAxisRotatePoint(axis.data, pt, pt);
			Crownpolydata->GetPoints()->SetPoint(i, pt);
		}
		Crownpolydata->Modified();

		//中心对其+法向量对其测试

		for (int i = 0; i < 3; i++)
		{
			vtkNew<vtkIterativeClosestPointTransform> icp1;
			icp1->SetSource(target_polydata);
			icp1->SetTarget(source_polydata);
			icp1->SetMaximumNumberOfIterations(50);
			icp1->GetLandmarkTransform()->SetModeToRigidBody();
			icp1->Modified();
			icp1->Update();

			vtkSmartPointer<vtkMatrix4x4> m1 = icp1->GetMatrix();
			m1->Invert();
			std::cout << "The resulting matrix is: " << *m1 << std::endl;
			SourcePre2Next.data[0] = m1->Element[0][0] * SourcePre2Next[0] + m1->Element[0][1] * SourcePre2Next[1] + m1->Element[0][2] * SourcePre2Next[2];
			SourcePre2Next.data[1] = m1->Element[1][0] * SourcePre2Next[0] + m1->Element[1][1] * SourcePre2Next[1] + m1->Element[1][2] * SourcePre2Next[2];
			SourcePre2Next.data[2] = m1->Element[2][0] * SourcePre2Next[0] + m1->Element[2][1] * SourcePre2Next[1] + m1->Element[2][2] * SourcePre2Next[2];
			SourcePre2Next.normalize();
			vtkNew<vtkTransformPolyDataFilter> icpTransformFilter1;
			icpTransformFilter1->SetInputData(source_polydata);
			icpTransformFilter1->SetTransform(icp1);
			icpTransformFilter1->Update();

			source_polydata = icpTransformFilter1->GetOutput();

			vtkNew<vtkTransformPolyDataFilter> icpTransformFilter11;
			icpTransformFilter11->SetInputData(Crownpolydata);
			icpTransformFilter11->SetTransform(icp1);
			icpTransformFilter11->Update();

			Crownpolydata = icpTransformFilter11->GetOutput();

			vtkNew<vtkIterativeClosestPointTransform> icp2;
			icp2->SetSource(target_polydata);
			icp2->SetTarget(source_polydata);
			icp2->SetMaximumNumberOfIterations(50);
			icp2->GetLandmarkTransform()->SetModeToSimilarity();
			icp2->Modified();
			icp2->Update();

			vtkSmartPointer<vtkMatrix4x4> m2 = icp2->GetMatrix();
			m2->Invert();
			std::cout << "The resulting matrix is: " << *m2 << std::endl;
			SourcePre2Next.data[0] = m2->Element[0][0] * SourcePre2Next[0] + m2->Element[0][1] * SourcePre2Next[1] + m2->Element[0][2] * SourcePre2Next[2];
			SourcePre2Next.data[1] = m2->Element[1][0] * SourcePre2Next[0] + m2->Element[1][1] * SourcePre2Next[1] + m2->Element[1][2] * SourcePre2Next[2];
			SourcePre2Next.data[2] = m2->Element[2][0] * SourcePre2Next[0] + m2->Element[2][1] * SourcePre2Next[1] + m2->Element[2][2] * SourcePre2Next[2];
			SourcePre2Next.normalize();
			vtkNew<vtkTransformPolyDataFilter> icpTransformFilter2;
			icpTransformFilter2->SetInputData(source_polydata);
			icpTransformFilter2->SetTransform(icp2);
			icpTransformFilter2->Update();

			source_polydata = icpTransformFilter2->GetOutput();

			vtkNew<vtkTransformPolyDataFilter> icpTransformFilter22;
			icpTransformFilter22->SetInputData(Crownpolydata);
			icpTransformFilter22->SetTransform(icp2);
			icpTransformFilter22->Update();

			Crownpolydata = icpTransformFilter22->GetOutput();
		}

		for (int i = 0; i < Crownpolydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double point[3];
			Crownpolydata->GetPoints()->GetPoint(i, point);
			point[0] = point[0] + target_center[0];
			point[1] = point[1] + target_center[1];
			point[2] = point[2] + target_center[2];
			Crownpolydata->GetPoints()->SetPoint(i, point);
		}
		Crownpolydata->Modified();

		for (int i = 0; i < source_polydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double point[3];
			source_polydata->GetPoints()->GetPoint(i, point);
			point[0] = point[0] + target_center[0];
			point[1] = point[1] + target_center[1];
			point[2] = point[2] + target_center[2];
			source_polydata->GetPoints()->SetPoint(i, point);
		}
		source_polydata->Modified();

		for (int i = 0; i < target_polydata->GetPoints()->GetNumberOfPoints(); i++)
		{
			double point[3];
			target_polydata->GetPoints()->GetPoint(i, point);
			point[0] = point[0] + target_center[0];
			point[1] = point[1] + target_center[1];
			point[2] = point[2] + target_center[2];
			target_polydata->GetPoints()->SetPoint(i, point);
		}
		target_polydata->Modified();

		////////////////////////////////////////////////////////////////获取前后牙齿的中心点
		Point_3 pre_center, next_center;
		vtkSmartPointer<vtkPolyData> pre_, next_;
		pre_ = vtkSmartPointer<vtkPolyData>::New();
		next_ = vtkSmartPointer<vtkPolyData>::New();
		if (selected_id == 1)
		{
			pre_ = teeth_data_initializer.m_labeledPolyData[selected_id + 1];
			next_ = teeth_data_initializer.m_labeledPolyData[selected_id + 2];
		}
		else if (selected_id == 7)
		{
			pre_ = teeth_data_initializer.m_labeledPolyData[selected_id - 1];
			next_ = teeth_data_initializer.m_labeledPolyData[selected_id - 2];
		}
		else
		{
			pre_ = teeth_data_initializer.m_labeledPolyData[selected_id - 1];
			next_ = teeth_data_initializer.m_labeledPolyData[selected_id + 1];
		}
		pre_center = Point_3(pre_->GetCenter()[0], pre_->GetCenter()[1], pre_->GetCenter()[2]);
		next_center = Point_3(next_->GetCenter()[0], next_->GetCenter()[1], next_->GetCenter()[2]);
		////////////////////////////////////////////////////////////////获取前后牙齿的中心点

		for (int t = 0; t < 10; t++)
		{
			SurfaceMesh source_mesh;
			source_mesh = VTK_PolyData2CGAL_Surface_Mesh(source_polydata);

			Tree target_tree(faces(mesh).first, faces(mesh).second, mesh);
			Tree source_tree(faces(source_mesh).first, faces(source_mesh).second, source_mesh);

			//五向十交点

			std::vector<Point_3> start_pts;
			start_pts.push_back(pre_center);
			start_pts.push_back(next_center);

			std::vector<Vector_3> Directions;
			Directions.push_back(Vector_3(pre_center[0] - next_center[0], pre_center[1] - next_center[1], pre_center[2] - next_center[2]));
			Directions.push_back(-Directions[0]);
			Directions.push_back(Vector_3(target_normal[0], target_normal[1], target_normal[2]));
			Directions.push_back(CGAL::cross_product(Directions[0], Directions[2]));
			Directions.push_back(-Directions[3]);
			Directions.push_back(-Directions[2]);
			std::vector<Ray_3> rays;
			for (int i = 0; i < start_pts.size(); i++)
			{
				for (int j = 0; j < Directions.size(); j++)
				{
					rays.push_back(Ray_3(start_pts[i], Directions[j]));
				}
			}
			std::vector<Point_3> target_5pts, source_5pts;
			double offset = 0.4;

			for (int i = 0; i < rays.size(); i++)
			{
				auto ray = rays[i];
				auto intersection1 = target_tree.first_intersection(ray);
				auto intersection2 = source_tree.first_intersection(ray);
				Point_3* p1 = new Point_3;  Point_3* p2 = new Point_3;
				if (intersection1 && intersection2)
				{
					p1 = boost::get<Point_3>(&(intersection1->first));

					//RenderSphere(*p, .25, Renderer, 0,1, 0, 1);
					if (i < 6)
					{
						double3 d(Directions[0][0], Directions[0][1], Directions[0][2]);
						d.normalize();
						p1 = &Point_3(p1->x() + d.data[0] * offset, p1->y() + d.data[0] * offset, p1->z() + d.data[0] * offset);
					}
					if (i >= 6)
					{
						double3 d(Directions[1][0], Directions[1][1], Directions[1][2]);
						d.normalize();
						p1 = &Point_3(p1->x() + d.data[0] * offset, p1->y() + d.data[0] * offset, p1->z() + d.data[0] * offset);
					}
					//target_5pts.push_back(*p1);

				}

				if (intersection1 && intersection2)
				{
					p2 = boost::get<Point_3>(&(intersection2->first));

					//cout << (p1->x() - p2->x()) * (p1->x() - p2->x()) + (p1->y() - p2->y()) * (p1->y() - p2->y()) + (p1->z() - p2->z()) * (p1->z() - p2->z()) <<endl;
					if ((p1->x() - p2->x()) * (p1->x() - p2->x()) + (p1->y() - p2->y()) * (p1->y() - p2->y()) + (p1->z() - p2->z()) * (p1->z() - p2->z()) < 4.0)
					{
						target_5pts.push_back(*p1);
						source_5pts.push_back(*p2);
						RenderSphere(*p1, .25, pipeline->Renderer, 1, 0, 0, 1);
						RenderSphere(*p2, .25, pipeline->Renderer, 0, 1, 0, 1);
					}

				}

			}
			//Renderer->AddActor(MakeActor(Crownpolydata, 0, 1, 0, 0.5));
			vtkSmartPointer<vtkPoints> target_5pts_points = vtkSmartPointer<vtkPoints>::New();
			for (int i = 0; i < target_5pts.size(); i++)
			{
				target_5pts_points->InsertNextPoint(target_5pts[i][0], target_5pts[i][1], target_5pts[i][2]);
			}

			vtkSmartPointer<vtkPoints> source_5pts_points = vtkSmartPointer<vtkPoints>::New();
			for (int i = 0; i < source_5pts.size(); i++)
			{
				source_5pts_points->InsertNextPoint(source_5pts[i][0], source_5pts[i][1], source_5pts[i][2]);
			}

			//五向十交点
			for (int i = 0; i < 1; i++)
			{
				vtkNew<vtkLandmarkTransform> icp;
				icp->SetSourceLandmarks(source_5pts_points);
				icp->SetTargetLandmarks(target_5pts_points);
				icp->SetModeToRigidBody();
				icp->Modified();
				icp->Update();

				vtkSmartPointer<vtkMatrix4x4> m1 = icp->GetMatrix();
				std::cout << "The resulting matrix is: " << *m1 << std::endl;
				//for (int j = 0; j < source_5pts_points->GetNumberOfPoints(); j++)
				//{
				//	double point[3];
				//	source_5pts_points->GetPoint(j, point);
				//	point[0] = m1->Element[0][0] * point[0] + m1->Element[0][1] * point[1] + m1->Element[0][2] * point[2] + m1->Element[0][3];
				//	point[1] = m1->Element[1][0] * point[0] + m1->Element[1][1] * point[1] + m1->Element[1][2] * point[2] + m1->Element[1][3];
				//	point[2] = m1->Element[2][0] * point[0] + m1->Element[2][1] * point[1] + m1->Element[2][2] * point[2] + m1->Element[2][3];
				//	source_5pts_points->SetPoint(j, point);
				//}
				//source_5pts_points->Modified();
				SourcePre2Next.data[0] = m1->Element[0][0] * SourcePre2Next[0] + m1->Element[0][1] * SourcePre2Next[1] + m1->Element[0][2] * SourcePre2Next[2];
				SourcePre2Next.data[1] = m1->Element[1][0] * SourcePre2Next[0] + m1->Element[1][1] * SourcePre2Next[1] + m1->Element[1][2] * SourcePre2Next[2];
				SourcePre2Next.data[2] = m1->Element[2][0] * SourcePre2Next[0] + m1->Element[2][1] * SourcePre2Next[1] + m1->Element[2][2] * SourcePre2Next[2];
				SourcePre2Next.normalize();
				vtkNew<vtkTransformPolyDataFilter> icpTransformFilter1;
				icpTransformFilter1->SetInputData(source_polydata);
				icpTransformFilter1->SetTransform(icp);
				icpTransformFilter1->Update();

				source_polydata = icpTransformFilter1->GetOutput();

				vtkNew<vtkTransformPolyDataFilter> icpTransformFilter11;
				icpTransformFilter11->SetInputData(Crownpolydata);
				icpTransformFilter11->SetTransform(icp);
				icpTransformFilter11->Update();


				Crownpolydata = icpTransformFilter11->GetOutput();

			}
		}



		//////////////////////////////////////////////////以下为测试代码和渲染部分 无需整合

		//auto crown_center = CGAL::midpoint(target_5pts[1], target_5pts[6]);







		//std::cout << file_path << std::endl;
		//std::cout << "checked_id: " << checked_id << ", selected_id: " << selected_id << std::endl;

		//SurfaceMesh crown_mesh;
		//crown_mesh = VTK_PolyData2CGAL_Surface_Mesh(Crownpolydata);
		//CGAL::Polygon_mesh_slicer<SurfaceMesh, K> slicer1(crown_mesh);
		//Polylines polyline1;
		//slicer1(K::Plane_3(crown_center, K::Vector_3(TargetPre2Next.data[0], TargetPre2Next.data[1], TargetPre2Next.data[2])), std::back_inserter(polyline1));
		//vtkNew<vtkPoints> points1;
		//for (int j = 0; j < polyline1[0].size(); j++)
		//{
		//	points1->InsertNextPoint(polyline1[0][j][0], polyline1[0][j][1], polyline1[0][j][2]);
		//}
		//vtkNew<vtkCellArray> cells1;
		//vtkNew<vtkPolyLine>Polyline1;
		//Polyline1->GetPointIds()->SetNumberOfIds(points1->GetNumberOfPoints());
		//for (int j = 0; j < points1->GetNumberOfPoints(); j++)
		//{
		//	Polyline1->GetPointIds()->SetId(j, j);
		//}
		//cells1->InsertNextCell(Polyline1);
		//vtkNew<vtkPolyData> polydata1;
		//polydata1->SetPoints(points1);
		//polydata1->SetLines(cells1);
		//vtkNew<vtkPolyDataMapper> mapper1;
		//mapper1->SetInputData(polydata1);
		//vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
		//actor1->SetMapper(mapper1);
		//actor1->GetProperty()->SetColor(1, 1, 0);
		////Renderer->AddActor(actor1);

		//auto dir = double3::crossProduct(target_normal, TargetPre2Next);

		//CGAL::Polygon_mesh_slicer<SurfaceMesh, K> slicer2(crown_mesh);
		//Polylines polyline2;
		//slicer2(K::Plane_3(crown_center, K::Vector_3(dir[0], dir[1], dir[2])), std::back_inserter(polyline2));
		//vtkNew<vtkPoints> points2;
		//for (int j = 0; j < polyline2[0].size(); j++)
		//{
		//	points2->InsertNextPoint(polyline2[0][j][0], polyline2[0][j][1], polyline2[0][j][2]);
		//}
		//vtkNew<vtkCellArray> cells2;
		//vtkNew<vtkPolyLine>Polyline2;
		//Polyline2->GetPointIds()->SetNumberOfIds(points2->GetNumberOfPoints());
		//for (int j = 0; j < points2->GetNumberOfPoints(); j++)
		//{
		//	Polyline2->GetPointIds()->SetId(j, j);
		//}
		//cells2->InsertNextCell(Polyline2);
		//vtkNew<vtkPolyData> polydata2;
		//polydata2->SetPoints(points2);
		//polydata2->SetLines(cells2);
		//vtkNew<vtkPolyDataMapper> mapper2;
		//mapper2->SetInputData(polydata2);
		//vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
		//actor2->SetMapper(mapper2);
		//actor2->GetProperty()->SetColor(0, 0, 1);
		////Renderer->AddActor(actor2);

		//CGAL::Polygon_mesh_slicer<SurfaceMesh, K> slicer3(crown_mesh);
		//Polylines polyline3;
		//slicer3(K::Plane_3(crown_center, K::Vector_3(target_normal[0], target_normal[1], target_normal[2])), std::back_inserter(polyline3));
		//vtkNew<vtkPoints> points3;
		//for (int j = 0; j < polyline3[0].size(); j++)
		//{
		//	points3->InsertNextPoint(polyline3[0][j][0], polyline3[0][j][1], polyline3[0][j][2]);
		//}
		//vtkNew<vtkCellArray> cells3;
		//vtkNew<vtkPolyLine>Polyline3;
		//Polyline3->GetPointIds()->SetNumberOfIds(points3->GetNumberOfPoints());
		//for (int j = 0; j < points3->GetNumberOfPoints(); j++)
		//{
		//	Polyline3->GetPointIds()->SetId(j, j);
		//}
		//cells3->InsertNextCell(Polyline3);
		//vtkNew<vtkPolyData> polydata3;
		//polydata3->SetPoints(points3);
		//polydata3->SetLines(cells3);
		//vtkNew<vtkPolyDataMapper> mapper3;
		//mapper3->SetInputData(polydata3);
		//vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
		//actor3->SetMapper(mapper3);
		//actor3->GetProperty()->SetColor(1, 0, 1);
		//Renderer->AddActor(actor3);
		std::cout << "source_polydata->GetNumberOfPoints():" << source_polydata->GetNumberOfPoints() << std::endl;
		// Visualize

		RenderPolydata(Crownpolydata, pipeline->Renderer, 0, 1, 0, 0.5);
		//Renderer->AddActor(MakeActor(source_polydata, 0, 0.5, 1, 0.2));
		//Renderer->AddActor(MakeActor(cut, 1, 1, 1, 1));
		RenderPolydata(PolyData, pipeline->Renderer, 1, 0.5, 1, 0.5);
		//Renderer->AddActor(MakeActor(target_polydata, 1, 1, 1, 1));
		//RenderPolydata(source_polydata, pipeline->Renderer, 0, 1, 0, 0.5);
		RenderPolydata(target_polydata, pipeline->Renderer, 1, 1, 1, 1);

		// Set up the camera and interactor.
		pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
		pipeline->Renderer->ResetCamera();
		// Set up the callback functions for mouse events.
		pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
		pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
		pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);
		pipeline->addObserver(vtkCommand::RightButtonPressEvent, RightPress);
		pipeline->addObserver(vtkCommand::RightButtonReleaseEvent, RightRelease);

		pipeline->RenderWindowInteractor->Start();

		// Clean up the pipeline after each run.
		delete pipeline;
	}
	
}