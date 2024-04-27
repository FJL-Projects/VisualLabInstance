// Created on: 1994-03-03
// Created by: Joelle CHAUVET
// Copyright (c) 1994-1999 Matra Datavision
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

#ifndef BRepFill_DataMapOfNodeShape_HeaderFile
#define BRepFill_DataMapOfNodeShape_HeaderFile

#include "MAT_Node.hxx"
#include "TopoDS_Shape.hxx"
#include "TColStd_MapTransientHasher.hxx"
#include "NCollection_DataMap.hxx"

typedef NCollection_DataMap<Handle(MAT_Node),TopoDS_Shape,TColStd_MapTransientHasher> BRepFill_DataMapOfNodeShape;
typedef NCollection_DataMap<Handle(MAT_Node),TopoDS_Shape,TColStd_MapTransientHasher>::Iterator BRepFill_DataMapIteratorOfDataMapOfNodeShape;


#endif
