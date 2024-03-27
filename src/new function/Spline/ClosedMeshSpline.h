#ifndef CLOSEDMESHSPLINE_H
#define CLOSEDMESHSPLINE_H
#include "Spline.h"
#include <algorithm>
#include <list>
#include <map>
#include "Triangulator2D.h"


class ClosedMeshSpline :public SplineBase
{
public:
	ClosedMeshSpline();
	~ClosedMeshSpline();
	//��ʼ��
	void initial(
		SurfaceMesh* sm,
		std::map<unsigned int, face_descriptor>& fmap,
		std::map<unsigned int, vertex_descriptor>& vmap,
		std::map<unsigned int, edge_descriptor>& emap,
		std::map<unsigned int, halfedge_descriptor>& hemap,
		SurfaceMesh::Property_map<vertex_descriptor, Point_2>& uvmap) override;
	//��ӵ�
	int add(unsigned int TriangleID, double PickPosition[3], double cameraDir[3] = NULL) override;
	//�ƶ����Ƶ�
	int move(unsigned int CtrlPtIndex, unsigned int TriangleID, double PickPosition[3], double cameraDir[3] = NULL) override;
	//������ɺ�������
	int addF(unsigned int& CtrlPtIndex, unsigned int TriangleID, double PickPosition[3], double cameraDir[3] = NULL) override;
	//�Ƴ���
	int remove(unsigned int CtrlPtIndex) override;

	//�޸�����
	int ChangeType(int IndexType)override;
	//����������
	void UpdateSpline(std::vector<MeshPoint>& SplinePoints)override;
	//���¿��Ƶ�
	void UpdateCtrlPoint(unsigned int& nIndex, int TriangleID, double PickPosition[3], State state)override;

	void GetEquidistantSamplingPoints(MeshPoint startCutP, MeshPoint endP, std::vector<MeshPoint>& crossPt, std::vector<MeshPoint>& EquidistantSamplingPt, double fSegmentLength);

	void BFSEdge(edge_descriptor& e, std::vector<bool>& vt_fvisit);

	void cutMesh();

public:
	//�Ⱦ������
	std::vector<MeshPoint> vtEquidistantSpline;
	bool bClosed;
	std::vector<int> m_newTri;
	std::vector<int> vVisit;
	std::vector<int> eVisit;
	int eIndex = 0;
	int vIndex = 0;
	SurfaceMesh cut;
	std::vector<Point_2> uvSpline;
};

#endif