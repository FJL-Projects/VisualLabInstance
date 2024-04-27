// Created on: 1997-01-09
// Created by: Jean Yves LEBEY
// Copyright (c) 1997-1999 Matra Datavision
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

#ifndef _TopOpeBRep_Bipoint_HeaderFile
#define _TopOpeBRep_Bipoint_HeaderFile

#include "Standard.hxx"
#include "Standard_DefineAlloc.hxx"
#include "Standard_Handle.hxx"

#include "Standard_Integer.hxx"



class TopOpeBRep_Bipoint 
{
public:

  DEFINE_STANDARD_ALLOC

  
  Standard_EXPORT TopOpeBRep_Bipoint();
  
  Standard_EXPORT TopOpeBRep_Bipoint(const Standard_Integer I1, const Standard_Integer I2);
  
  Standard_EXPORT Standard_Integer I1() const;
  
  Standard_EXPORT Standard_Integer I2() const;




protected:





private:



  Standard_Integer myI1;
  Standard_Integer myI2;


};







#endif // _TopOpeBRep_Bipoint_HeaderFile
