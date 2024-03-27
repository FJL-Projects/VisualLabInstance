#include "Spline/Spline.h"

/**
 * Find all the intersection points between the edge formed by uv1 and uv2 and the surface mesh.
 *
 * This function recursively traverses the mesh to find all intersection points, starting from the halfedge he.
 * If uv2 is inside the triangle of the current halfedge, add the two UV points to the MeshSpline and return the current halfedge.
 * Otherwise, try to find the intersection between the edge and the line connecting uv1 and uv2.
 * If an intersection is found, add the starting point uv1 to MeshSpline and call the function recursively on the opposite halfedge, with the intersection point as the new uv1.
 * The process is repeated on the next halfedge until all intersections are found, or no intersection is found on all three halfedges.
 *
 * @param he The starting halfedge for the traversal.
 * @param uv1 The starting point of the edge.
 * @param uv2 The end point of the edge.
 * @param first A boolean indicating whether this is the first call of the function.
 * @return The halfedge where the traversal ends.
 */
halfedge_descriptor SplineBase::FindAllCrosses(halfedge_descriptor he, Point_2 uv1, Point_2 uv2, bool first)
{
	if (!pMesh->is_valid(he)||pMesh->is_border(he) || bDrawOnBorderFlag == 1)
	{
		bDrawOnBorderFlag = 1;
		return he;
	}
	if (inTriangle(he, uv2))
	{
		vtSplinePoints.emplace_back(pMesh->face(he).idx(), uv2, UV2XYZ(he, uv2), he);
		return he;
	}
	halfedge_descriptor he0 = he;
	Point_2 crossVer;
	if (first)
	{
		crossVer = EdgeCross(uv1, uv2, UVmap[pMesh->source(he0)], UVmap[pMesh->target(he0)]);
		if (crossVer != Point_2(0, 0))
		{
			vtSplinePoints.emplace_back(pMesh->face(pMesh->next(he0)).idx(), crossVer, UV2XYZ(pMesh->opposite(he0), crossVer), pMesh->opposite(he0));
			return FindAllCrosses(pMesh->opposite(he0), crossVer, uv2, 0);
		}
	}
	halfedge_descriptor he1 = pMesh->next(he0);
	crossVer = EdgeCross(uv1, uv2, UVmap[pMesh->source(he1)], UVmap[pMesh->target(he1)]);
	if (crossVer != Point_2(0, 0))
	{
		vtSplinePoints.emplace_back(pMesh->face(pMesh->opposite(he1)).idx(), crossVer, UV2XYZ(pMesh->opposite(he1), crossVer), pMesh->opposite(he1));
		return FindAllCrosses(pMesh->opposite(he1), crossVer, uv2, 0);
	}
	halfedge_descriptor he2 = pMesh->next(he1);
	crossVer = EdgeCross(uv1, uv2, UVmap[pMesh->source(he2)], UVmap[pMesh->target(he2)]);
	if (crossVer != Point_2(0, 0))
	{
		vtSplinePoints.emplace_back(pMesh->face(pMesh->opposite(he2)).idx(), crossVer, UV2XYZ(pMesh->opposite(he2), crossVer), pMesh->opposite(he2));
		return FindAllCrosses(pMesh->opposite(he2), crossVer, uv2, 0);
	}
	std::cout << "error" << std::endl;
}

Point_2 SplineBase::EdgeCross(Point_2& p1, Point_2& p2, Point_2& p3, Point_2& p4)
{
	auto ans = CGAL::intersection(Kernel::Segment_2(p1, p2), Kernel::Segment_2(p3, p4));
	if (ans)
	{
		if (const Kernel::Point_2* s = boost::get<Kernel::Point_2>(&*ans))
		{
			return Point_2(s->x(), s->y());
		}
	}
	else
	{
		return Point_2(0, 0);
	}
}

Point_3 SplineBase::PlaneCrossSegment(Plane_3& P, Segment_3& S)
{
	auto ans = CGAL::intersection(P, S);
	if (ans)
	{
		if (const Kernel::Point_3* s = boost::get<Kernel::Point_3>(&*ans))
		{
			return Point_3(s->x(), s->y(), s->z());
		}
	}
	else
	{
		return Point_3(0, 0, 0);
	}
}

double3 SplineBase::GetHalfEdgeNormal(halfedge_descriptor& he)
{
	Point_3& P0 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->target(he));
	Point_3& P1 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->target(pMesh->next(he)));
	Point_3& P2 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->source(he));
	double3 Normal((P1.y() - P0.y()) * (P2.z() - P0.z()) - (P1.z() - P0.z()) * (P2.y() - P0.y()),
		(P1.z() - P0.z()) * (P2.x() - P0.x()) - (P1.x() - P0.x()) * (P2.z() - P0.z()),
		(P1.x() - P0.x()) * (P2.y() - P0.y()) - (P1.y() - P0.y()) * (P2.x() - P0.x()));
	Normal.normalize();
	return Normal;
}


Point_2 SplineBase::GetPointUV(halfedge_descriptor& he, Point_3& pt)
{
	Point_3& p0 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->target(he));
	Point_3& p1 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->target(pMesh->next(he)));
	Point_3& p2 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->source(he));
	Point_2 uv[3];
	uv[0] = UVmap[pMesh->target(he)];
	uv[1] = UVmap[pMesh->target(pMesh->next(he))];
	uv[2] = UVmap[pMesh->source(he)];
	double bary1, bary2, bary3;
	GetBary(p0, p1, p2, pt, bary1, bary2, bary3);
	return Point_2(uv[0].x() * bary1 + uv[1].x() * bary2 + uv[2].x() * bary3, uv[0].y() * bary1 + uv[1].y() * bary2 + uv[2].y() * bary3);
}

bool SplineBase::GetBary(Point_3& v1, Point_3& v2, Point_3& v3, Point_3& vp, double& Bary1, double& Bary2, double& Bary3)
{
	Point_3 v02(v1.x() - v3.x(), v1.y() - v3.y(), v1.z() - v3.z());
	Point_3 v12(v2.x() - v3.x(), v2.y() - v3.y(), v2.z() - v3.z());
	Point_3 vp2(vp.x() - v3.x(), vp.y() - v3.y(), vp.z() - v3.z());

	double d0202 = v02.x() * v02.x() + v02.y() * v02.y() + v02.z() * v02.z();
	double d0212 = v02.x() * v12.x() + v02.y() * v12.y() + v02.z() * v12.z();
	double d1212 = v12.x() * v12.x() + v12.y() * v12.y() + v12.z() * v12.z();
	double d02p2 = v02.x() * vp2.x() + v02.y() * vp2.y() + v02.z() * vp2.z();
	double d12p2 = v12.x() * vp2.x() + v12.y() * vp2.y() + v12.z() * vp2.z();

	Bary1 = (d1212 * d02p2 - d0212 * d12p2) / (d0202 * d1212 - d0212 * d0212);
	Bary2 = (d0202 * d12p2 - d0212 * d02p2) / (d0202 * d1212 - d0212 * d0212);
	Bary3 = static_cast<double>(1.0) - Bary1 - Bary2;

	return ((Bary1 >= 0) && (Bary2 >= 0) && (Bary1 + Bary2 < static_cast<double>(1.0)));
}

double* SplineBase::UV2XYZ(halfedge_descriptor& he, Point_2& uv)
{
	Point_3& p0 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->target(he));
	Point_3& p1 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->target(pMesh->next(he)));
	Point_3& p2 = get(CGAL::get(CGAL::vertex_point, *pMesh), pMesh->source(he));

	Point_2& v1 = UVmap[pMesh->target(he)];
	Point_2& v2 = UVmap[pMesh->target(pMesh->next(he))];
	Point_2& v3 = UVmap[pMesh->source(he)];
	Point_2 v02(v1.x() - v3.x(), v1.y() - v3.y());
	Point_2 v12(v2.x() - v3.x(), v2.y() - v3.y());
	Point_2 vp2(uv.x() - v3.x(), uv.y() - v3.y());

	double d0202 = v02.x() * v02.x() + v02.y() * v02.y();
	double d0212 = v02.x() * v12.x() + v02.y() * v12.y();
	double d1212 = v12.x() * v12.x() + v12.y() * v12.y();
	double d02p2 = v02.x() * vp2.x() + v02.y() * vp2.y();
	double d12p2 = v12.x() * vp2.x() + v12.y() * vp2.y();

	double Bary1 = (d1212 * d02p2 - d0212 * d12p2) / (d0202 * d1212 - d0212 * d0212);
	double Bary2 = (d0202 * d12p2 - d0212 * d02p2) / (d0202 * d1212 - d0212 * d0212);
	static double ans[3];
	ans[0] = p0.x() * Bary1 + p1.x() * Bary2 + p2.x() * (1.0f - Bary1 - Bary2);
	ans[1] = p0.y() * Bary1 + p1.y() * Bary2 + p2.y() * (1.0f - Bary1 - Bary2);
	ans[2] = p0.z() * Bary1 + p1.z() * Bary2 + p2.z() * (1.0f - Bary1 - Bary2);
	return ans;
}

bool SplineBase::inTriangle(halfedge_descriptor& he, Point_2& vp)
{
	if (pMesh->is_border(he))
	{
		cout << "Halfedge is on border!";
		return 0;
	}
	Point_2& v1 = UVmap[pMesh->target(he)];
	Point_2& v2 = UVmap[pMesh->target(pMesh->next(he))];
	Point_2& v3 = UVmap[pMesh->source(he)];
	Point_2 v02(v1.x() - v3.x(), v1.y() - v3.y());
	Point_2 v12(v2.x() - v3.x(), v2.y() - v3.y());
	Point_2 vp2(vp.x() - v3.x(), vp.y() - v3.y());

	double d0202 = v02[0] * v02[0] + v02[1] * v02[1];
	double d0212 = v02[0] * v12[0] + v02[1] * v12[1];
	double d1212 = v12[0] * v12[0] + v12[1] * v12[1];
	double d02p2 = v02[0] * vp2[0] + v02[1] * vp2[1];
	double d12p2 = v12[0] * vp2[0] + v12[1] * vp2[1];

	return (((d1212 * d02p2 - d0212 * d12p2) / (d0202 * d1212 - d0212 * d0212) >= 0) &&
		((d0202 * d12p2 - d0212 * d02p2) / (d0202 * d1212 - d0212 * d0212) >= 0) &&
		((d1212 * d02p2 - d0212 * d12p2) / (d0202 * d1212 - d0212 * d0212) + (d0202 * d12p2 - d0212 * d02p2) / (d0202 * d1212 - d0212 * d0212) <= static_cast<double>(1.0)));
}

