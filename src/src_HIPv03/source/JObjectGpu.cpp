//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JObjectGpu.cpp \brief Implements the class \ref JObjectGpu.

#include "JObjectGpu.h"
#include "JException.h"
#include "Functions.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>

//##############################################################################
//# JObjectGpu
//##############################################################################
//==============================================================================
/// Throws exception related to a HIP error.
//==============================================================================
void JObjectGpu::RunExceptionHip(const std::string &srcfile,int srcline
  ,const std::string &classname,const std::string &method
  ,hipError_t hiperr,std::string msg)const
{
  msg=msg+fun::PrintStr(" (HIP error %d (%s)).\n",hiperr,hipGetErrorString(hiperr));
  throw JException(srcfile,srcline,classname,method,msg,"");
}

//==============================================================================
/// Checks HIP error and throws exception.
//==============================================================================
void JObjectGpu::CheckHipError(const std::string &srcfile,int srcline
  ,const std::string &classname,const std::string &method
  ,std::string msg)const
{
  hipError_t hiperr=hipGetLastError();
  if(hiperr!=hipSuccess)RunExceptionHip(srcfile,srcline,classname,method,hiperr,msg);
}
