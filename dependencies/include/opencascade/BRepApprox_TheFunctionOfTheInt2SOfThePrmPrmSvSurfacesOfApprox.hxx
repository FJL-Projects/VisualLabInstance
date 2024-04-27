// Created on: 1995-06-06
// Created by: Jean Yves LEBEY
// Copyright (c) 1995-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#ifndef _BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox_HeaderFile
#define _BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox_HeaderFile

#include "Standard.hxx"
#include "Standard_DefineAlloc.hxx"
#include "Standard_Handle.hxx"

#include "Standard_Address.hxx"
#include "gp_Pnt.hxx"
#include "Standard_Real.hxx"
#include "Standard_Boolean.hxx"
#include "gp_Vec.hxx"
#include "IntImp_ConstIsoparametric.hxx"
#include "math_FunctionSetWithDerivatives.hxx"
#include "Standard_Integer.hxx"
#include "math_Vector.hxx"
#include "TColStd_Array1OfReal.hxx"
#include "gp_Dir.hxx"
#include "gp_Dir2d.hxx"
class StdFail_UndefinedDerivative;
class Standard_ConstructionError;
class BRepAdaptor_Surface;
class BRepApprox_SurfaceTool;
class math_Matrix;



class BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox  : public math_FunctionSetWithDerivatives
{
public:

  DEFINE_STANDARD_ALLOC

  Standard_EXPORT BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox(const BRepAdaptor_Surface& S1, const BRepAdaptor_Surface& S2);
  
  Standard_EXPORT Standard_Integer NbVariables() const;
  
  Standard_EXPORT Standard_Integer NbEquations() const;
  
  Standard_EXPORT Standard_Boolean Value (const math_Vector& X, math_Vector& F);
  
  Standard_EXPORT Standard_Boolean Derivatives (const math_Vector& X, math_Matrix& D);
  
  Standard_EXPORT Standard_Boolean Values (const math_Vector& X, math_Vector& F, math_Matrix& D);
  
  Standard_EXPORT void ComputeParameters (const IntImp_ConstIsoparametric ChoixIso, const TColStd_Array1OfReal& Param, math_Vector& UVap, math_Vector& BornInf, math_Vector& BornSup, math_Vector& Tolerance);
  
  //! returns somme des fi*fi
    Standard_Real Root() const;
  
    gp_Pnt Point() const;
  
  Standard_EXPORT Standard_Boolean IsTangent (const math_Vector& UVap, TColStd_Array1OfReal& Param, IntImp_ConstIsoparametric& BestChoix);
  
    gp_Dir Direction() const;
  
    gp_Dir2d DirectionOnS1() const;
  
    gp_Dir2d DirectionOnS2() const;
  
    const BRepAdaptor_Surface& AuxillarSurface1() const;
  
    const BRepAdaptor_Surface& AuxillarSurface2() const;




protected:





private:



  Standard_Address surf1;
  Standard_Address surf2;
  gp_Pnt pntsol1;
  gp_Pnt pntsol2;
  Standard_Real f[3];
  Standard_Boolean compute;
  Standard_Boolean tangent;
  Standard_Real tgduv[4];
  gp_Vec dpuv[4];
  IntImp_ConstIsoparametric chxIso;
  Standard_Real paramConst;
  Standard_Real ua0;
  Standard_Real va0;
  Standard_Real ua1;
  Standard_Real va1;
  Standard_Real ub0;
  Standard_Real vb0;
  Standard_Real ub1;
  Standard_Real vb1;
  Standard_Real ures1;
  Standard_Real ures2;
  Standard_Real vres1;
  Standard_Real vres2;


};

#define ThePSurface BRepAdaptor_Surface
#define ThePSurface_hxx <BRepAdaptor_Surface.hxx>
#define ThePSurfaceTool BRepApprox_SurfaceTool
#define ThePSurfaceTool_hxx <BRepApprox_SurfaceTool.hxx>
#define IntImp_ZerParFunc BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox
#define IntImp_ZerParFunc_hxx <BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox.hxx>

#include "IntImp_ZerParFunc.lxx"

#undef ThePSurface
#undef ThePSurface_hxx
#undef ThePSurfaceTool
#undef ThePSurfaceTool_hxx
#undef IntImp_ZerParFunc
#undef IntImp_ZerParFunc_hxx




#endif // _BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox_HeaderFile
