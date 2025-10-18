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

/// \file FunSphEos_iker.h \brief Implementa funciones de dispositivo HIP para la Ecuación de Estado para SPH.

#include "TypesDef.h"
#include "DualSphDef.h"
#include <hip/hip_runtime.h>

/// Implementa funciones de dispositivo HIP para la Ecuación de Estado para SPH.
namespace cufsph{

//##############################################################################
//# Ecuación de Estado - Monaghan, 1994
//##############################################################################
//#define CTE_AVAILABLE
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Devuelve la presión a partir de la densidad usando la ecuación de estado 
/// basada en [Monaghan, 1994].
//------------------------------------------------------------------------------
__device__ float ComputePressMonaghanCte(float rhop){ 
  return(CTE.cteb*(pow(rhop*CTE.ovrhopzero,CTE.gamma)-1.0f));
}
#endif
//------------------------------------------------------------------------------
/// Devuelve la presión a partir de la densidad usando la ecuación de estado 
/// basada en [Monaghan, 1994].
//------------------------------------------------------------------------------
__device__ float ComputePressMonaghan(float rhop,float rhop0,float b,float gamma){ 
  return(b*(pow(rhop/rhop0,gamma)-1.0f));
}


//##############################################################################
//# Ecuación de Estado por Defecto.
//##############################################################################
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Devuelve la presión a partir de la densidad usando la ecuación de estado por defecto.
//------------------------------------------------------------------------------
__device__ float ComputePressCte(float rhop){ 
  return(ComputePressMonaghanCte(rhop));
}
#endif
//------------------------------------------------------------------------------
/// Devuelve la presión a partir de la densidad usando la ecuación de estado por defecto.
//------------------------------------------------------------------------------
__device__ float ComputePress(float rhop,float rhop0,float b,float gamma,float cs0){ 
  return(ComputePressMonaghan(rhop,rhop0,b,gamma));
}


}
