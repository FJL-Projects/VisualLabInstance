// Created on: 1995-09-18
// Created by: Bruno DUMORTIER
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

#ifndef _BRepOffsetAPI_MakeEvolved_HeaderFile
#define _BRepOffsetAPI_MakeEvolved_HeaderFile

#include "Standard.hxx"
#include "Standard_DefineAlloc.hxx"
#include "Standard_Handle.hxx"

#include "BRepFill_Evolved.hxx"
#include "BRepFill_AdvancedEvolved.hxx"
#include "BRepBuilderAPI_MakeShape.hxx"
#include "GeomAbs_JoinType.hxx"
#include "Standard_Boolean.hxx"
#include "Standard_Real.hxx"
#include "TopTools_ListOfShape.hxx"
class TopoDS_Wire;
class TopoDS_Shape;


//! Describes functions to build evolved shapes.
//! An evolved shape is built from a planar spine (face or
//! wire) and a profile (wire). The evolved shape is the
//! unlooped sweep (pipe) of the profile along the spine.
//! Self-intersections are removed.
//! A MakeEvolved object provides a framework for:
//! - defining the construction of an evolved shape,
//! - implementing the construction algorithm, and
//! - consulting the result.
//! Computes an Evolved by
//! 1 - sweeping a profile along a spine.
//! 2 - removing the self-intersections.
//!
//! The Profile is expected to be planar and can be a line
//! (which lies in infinite number of planes).
//!
//! The profile is defined in a Referential R. The position of
//! the profile at the current point of the  spine is given by
//! confusing R  and the local  referential given by (  D0, D1
//! and the normal of the Spine).
//!
//! The coordinate system is determined by theIsAxeProf argument:
//! - if theIsAxeProf is true, R is the global coordinate system,
//! - if theIsAxeProf is false, R is computed so that:
//!     * its origin is given by the point on the spine which is
//!         closest to the profile,
//!     * its "X Axis" is given by the tangent to the spine at this point, and
//!     * its "Z Axis" is the normal to the plane which contains the spine.
//!
//! theJoinType defines the type of pipe generated by the salient
//! vertices of the spine. The default type is GeomAbs_Arc
//! where the vertices generate revolved pipes about the
//! axis passing along the vertex and the normal to the
//! plane of the spine. At present, this is the only
//! construction type implemented.
//! 
//! if <theIsSolid> is TRUE the Shape result is completed to be a
//! solid or a compound of solids.
//! 
//! If theIsProfOnSpine == TRUE then the profile must connect with the spine.
//!
//! If theIsVolume option is switched on then self-intersections
//! in the result of Pipe-algorithm will be removed by 
//! BOPAlgo_MakerVolume algorithm. At that the arguments
//! "theJoinType", "theIsAxeProf", "theIsProfOnSpine" are not used.

class BRepOffsetAPI_MakeEvolved  : public BRepBuilderAPI_MakeShape
{
public:

  DEFINE_STANDARD_ALLOC

  
  Standard_EXPORT BRepOffsetAPI_MakeEvolved();
  
  //! Constructs an evolved shape by sweeping the profile
  //! (theProfile) along the spine (theSpine).
  //! theSpine can be shape only of type wire or face.
  //! See description to this class for detailed information.
  Standard_EXPORT BRepOffsetAPI_MakeEvolved(const TopoDS_Shape& theSpine,
                                            const TopoDS_Wire& theProfile,
                                            const GeomAbs_JoinType theJoinType = GeomAbs_Arc,
                                            const Standard_Boolean theIsAxeProf = Standard_True,
                                            const Standard_Boolean theIsSolid = Standard_False,
                                            const Standard_Boolean theIsProfOnSpine = Standard_False,
                                            const Standard_Real theTol = 0.0000001,
                                            const Standard_Boolean theIsVolume = Standard_False,
                                            const Standard_Boolean theRunInParallel = Standard_False);
  
  Standard_EXPORT const BRepFill_Evolved& Evolved() const;
  
  //! Builds the resulting shape (redefined from MakeShape).
  Standard_EXPORT virtual void Build(const Message_ProgressRange& theRange = Message_ProgressRange()) Standard_OVERRIDE;
  
  //! Returns   the  shapes  created  from   a  subshape
  //! <SpineShape>  of     the  spine   and   a subshape
  //! <ProfShape> on the profile.
  Standard_EXPORT const TopTools_ListOfShape& GeneratedShapes (const TopoDS_Shape& SpineShape, const TopoDS_Shape& ProfShape) const;
  
  //! Return the face Top if <Solid> is True in the constructor.
  Standard_EXPORT const TopoDS_Shape& Top() const;
  
  //! Return the face Bottom  if <Solid> is True in the constructor.
  Standard_EXPORT const TopoDS_Shape& Bottom() const;




protected:





private:


  BRepFill_Evolved myEvolved;
  BRepFill_AdvancedEvolved myVolume;
  Standard_Boolean myIsVolume;

};







#endif // _BRepOffsetAPI_MakeEvolved_HeaderFile
