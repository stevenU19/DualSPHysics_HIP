//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 Este archivo es parte de DualSPHysics.

 DualSPHysics es software libre: puedes redistribuirlo y/o modificarlo bajo los términos de la GNU Lesser General Public License 
 publicada por la Free Software Foundation, ya sea la versión 2.1 de la licencia, o (a tu elección) cualquier versión posterior.

 DualSPHysics se distribuye con la esperanza de que sea útil, pero SIN NINGUNA GARANTÍA; sin siquiera la garantía implícita de
 COMERCIABILIDAD o APTITUD PARA UN PROPÓSITO PARTICULAR. Consulta la GNU Lesser General Public License para más detalles.

 Debiste haber recibido una copia de la GNU Lesser General Public License junto con DualSPHysics. Si no es así, consulta <http://www.gnu.org/licenses/>.
*/

/// \file JCellDivDataGpu.h \brief Declara estructuras e implementa funciones en línea para la búsqueda de vecinos.

#ifndef _JCellDivDataGpu_
#define _JCellDivDataGpu_

#include "DualSphDef.h"
#include <hip/hip_runtime.h>

/// Estructura con datos de la división celular para la búsqueda de vecinos en la GPU.
typedef struct{
  TpMgDivMode axis;
  tuint3 ncells;
  int scelldiv;         ///< Valor para dividir el KernelSize (1 o 2)
  int4 nc;
  unsigned cellfluid;
  int3 cellzero;
  const int2* beginendcell; //- int2*
  float scell;
  unsigned domcellcode;
  double3 domposmin;
  float kernelsize2;   ///< Distancia máxima de interacción al cuadrado (KernelSize^2).
  float poscellsize;   ///< Tamaño de las celdas usado para codificar PosCell (generalmente es KernelSize).
}StDivDataGpu;

//==============================================================================
/// Devuelve una estructura StDivDataGpu vacía.
//==============================================================================
inline StDivDataGpu DivDataGpuNull(){
  StDivDataGpu c={MGDIV_None,TUint3(0),0,{0,0,0,0},0,{0,0,0},NULL,0,0,{0,0,0},0,0};
  return(c);
}

//==============================================================================
/// Devuelve la estructura con datos para la búsqueda de vecinos en Single-GPU.
//==============================================================================
inline StDivDataGpu MakeDivDataGpu(int scelldiv,const tuint3 &ncells
  ,const tuint3 &cellmin,const int2* beginendcell,float scell,unsigned domcellcode
  ,const tdouble3 &domposmin,float kernelsize2,float poscellsize)
{
  StDivDataGpu ret;
  ret.axis=MGDIV_Z;
  ret.scelldiv=scelldiv;
  ret.ncells=ncells;
  ret.nc.x=int(ncells.x);
  ret.nc.y=int(ncells.y);
  ret.nc.z=int(ncells.z);
  ret.nc.w=int(ncells.x*ncells.y);  //- Para Single-GPU.
  ret.cellfluid=ret.nc.w*ret.nc.z+1;
  ret.cellzero.x=int(cellmin.x);
  ret.cellzero.y=int(cellmin.y);
  ret.cellzero.z=int(cellmin.z);
  ret.beginendcell=beginendcell;
  ret.scell=scell;
  ret.domcellcode=domcellcode;
  ret.domposmin.x=domposmin.x;
  ret.domposmin.y=domposmin.y;
  ret.domposmin.z=domposmin.z;
  ret.kernelsize2=kernelsize2;
  ret.poscellsize=poscellsize;
  return(ret);
}

#endif
