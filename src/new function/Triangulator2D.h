#ifndef TRIANGULATOR_2D_H
#define TRIANGULATOR_2D_H
#include <vector>

class Triangulator2D
{
public:
	struct InputData {
		std::vector<double> vPoints;
		std::vector<int> vPointMarkers;
		std::vector<int> vSegments;

		std::vector<double> vHoles;
	};
	struct OutputData_CDT {
		std::vector<double> vPoints;
		std::vector<int> vPointMarkers;
		std::vector<int> vTriangles;
	};

	Triangulator2D(void);
	~Triangulator2D(void);

	void ClearInput(bool bFreeMem = false);
	void PreAllocateInput(int nPoints, int nSegments);

	inline unsigned int AddPoint(double x, double y, int nMarker = 0)
	{
		m_input.vPoints.push_back(x);  m_input.vPoints.push_back(y);
		m_input.vPointMarkers.push_back(nMarker); return ((unsigned int)m_input.vPoints.size() - 2) / 2;
	}
	inline void AddSegment(int v1, int v2)
	{
		m_input.vSegments.push_back(v1); m_input.vSegments.push_back(v2);
	}
	inline void AddHole(double x, double y)
	{
		m_input.vHoles.push_back(x); m_input.vHoles.push_back(y);
	}



	bool Compute();

	//! if bCompact is true, will rewrite mesh to skip vertices which are not used in triangles



	const std::vector<int>& GetOutputMarkers_MeshVtx() { return m_VtxMarkers; }

	void SetSubdivideOuterSegments(bool bAllowSubdiv) { m_bNoSubdivdeOuterSegments = !bAllowSubdiv; }
	void SetSubdivideAnySegments(bool bAllowSubdiv) { m_bNoSubdivideAnySegments = !bAllowSubdiv; }
	void SetEnclosingSegmentsProvided(bool bProvided) { m_bEnclosingSegmentsProvided = bProvided; }

protected:
	InputData m_input;
	OutputData_CDT m_output;


	std::vector<int> m_VtxMarkers;

	bool m_bEnclosingSegmentsProvided;
	bool m_bNoSubdivdeOuterSegments;
	bool m_bNoSubdivideAnySegments;


public:
	// output stuff - shouldn't use because i != mesh IDs...
	size_t GetOutputVertexCount();
	bool GetOutputVertex(unsigned int i, double* pPoint, int& nMarker);
	int GetOutputMarker(unsigned int i);
	size_t GetOutputTriCount();
	bool GetOutputTri(unsigned int i, int* pTri);

	const std::vector<int>& GetPointMarkers() { return m_output.vPointMarkers; }

};
#endif  // __RMS_TRIANGULATOR_2D_H__