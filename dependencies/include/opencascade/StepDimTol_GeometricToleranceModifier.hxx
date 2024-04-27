// Created on: 2015-07-10
// Created by: Irina KRYLOVA
// Copyright (c) 2015 OPEN CASCADE SAS
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

#ifndef _StepDimTol_GeometricToleranceModifier_HeaderFile
#define _StepDimTol_GeometricToleranceModifier_HeaderFile

#include "Standard_PrimitiveTypes.hxx"

enum StepDimTol_GeometricToleranceModifier {
  StepDimTol_GTMAnyCrossSection,
  StepDimTol_GTMCommonZone,
  StepDimTol_GTMEachRadialElement,
  StepDimTol_GTMFreeState,
  StepDimTol_GTMLeastMaterialRequirement,
  StepDimTol_GTMLineElement,
  StepDimTol_GTMMajorDiameter,
  StepDimTol_GTMMaximumMaterialRequirement,
  StepDimTol_GTMMinorDiameter,
  StepDimTol_GTMNotConvex,
  StepDimTol_GTMPitchDiameter,
  StepDimTol_GTMReciprocityRequirement,
  StepDimTol_GTMSeparateRequirement,
  StepDimTol_GTMStatisticalTolerance,
  StepDimTol_GTMTangentPlane
};

#endif
