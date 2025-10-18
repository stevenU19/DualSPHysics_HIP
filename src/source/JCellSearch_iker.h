//HEAD_DSCODES
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

/// \file JCellSearch_iker.h \brief Implementa funciones del dispositivo HIP para la búsqueda de vecindades.

#ifndef _JCellSearch_iker_
#define _JCellSearch_iker_

#include "DualSphDef.h"

/// Implementa funciones del dispositivo HIP para la búsqueda de vecindades.
namespace cunsearch {
//#define CTE_AVAILABLE
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Devuelve los límites de la celda para la interacción.
/// Returns cell limits for the interaction.
//------------------------------------------------------------------------------
__device__ void InitCte(unsigned rcell
  ,const int &scelldiv,const int4 &nc,const int3 &cellzero
  ,int &ini1,int &fin1,int &ini2,int &fin2,int &ini3,int &fin3)
{
  // Obtiene los límites de interacción.
  const int cx=PC__Cellx(CTE.cellcode,rcell)-cellzero.x;
  const int cy=PC__Celly(CTE.cellcode,rcell)-cellzero.y;
  const int cz=PC__Cellz(CTE.cellcode,rcell)-cellzero.z;
  // Ordena los componentes según el eje utilizado en el orden de la celda.
  const int c1=(CTE.axis==MGDIV_X? cy: cx);
  const int c2=(CTE.axis==MGDIV_Z? cy: cz);
  const int c3=(CTE.axis==MGDIV_Z? cz: (CTE.axis==MGDIV_X? cx: cy));
  // Código para scelldiv 1 o 2 pero no cero.
  ini1=c1-min(c1,scelldiv);
  fin1=c1+min(nc.x-c1-1,scelldiv)+1;
  ini2=c2-min(c2,scelldiv);
  fin2=c2+min(nc.y-c2-1,scelldiv)+1;
  ini3=c3-min(c3,scelldiv);
  fin3=c3+min(nc.z-c3-1,scelldiv)+1;
  // Computa los límites absolutos.
  ini3*=nc.w; fin3*=nc.w;
  ini2*=nc.x; fin2*=nc.x;
}

//------------------------------------------------------------------------------
/// Devuelve los límites de la celda para la interacción.
/// Returns cell limits for the interaction.
//------------------------------------------------------------------------------
__device__ void InitCte(const double &px,const double &py,const double &pz
  ,const int &scelldiv,const int4 &nc,const int3 &cellzero
  ,int &ini1,int &fin1,int &ini2,int &fin2,int &ini3,int &fin3)
{
  // Obtiene los límites de interacción.
  const int cx=int((px-CTE.domposminx)/CTE.scell)-cellzero.x;
  const int cy=int((py-CTE.domposminy)/CTE.scell)-cellzero.y;
  const int cz=int((pz-CTE.domposminz)/CTE.scell)-cellzero.z;
  // Ordena los componentes según el eje utilizado en el orden de la celda.
  const int c1=(CTE.axis==MGDIV_X? cy: cx);
  const int c2=(CTE.axis==MGDIV_Z? cy: cz);
  const int c3=(CTE.axis==MGDIV_Z? cz: (CTE.axis==MGDIV_X? cx: cy));
  // Código para scelldiv 1 o 2 pero no cero.
  ini1=c1-min(c1,scelldiv);
  fin1=c1+min(nc.x-c1-1,scelldiv)+1;
  ini2=c2-min(c2,scelldiv);
  fin2=c2+min(nc.y-c2-1,scelldiv)+1;
  ini3=c3-min(c3,scelldiv);
  fin3=c3+min(nc.z-c3-1,scelldiv)+1;
  // Computa los límites absolutos.
  ini3*=nc.w; fin3*=nc.w;
  ini2*=nc.x; fin2*=nc.x;
}

#endif

//------------------------------------------------------------------------------
/// Devuelve los límites de la celda para la interacción.
/// Returns cell limits for the interaction.
//------------------------------------------------------------------------------
__device__ void Initsp(unsigned rcell
  ,const unsigned &axis,const unsigned &cellcode
  ,const int &scelldiv,const int4 &nc,const int3 &cellzero
  ,int &ini1,int &fin1,int &ini2,int &fin2,int &ini3,int &fin3)
{
  // Obtiene los límites de interacción.
  const int cx=PC__Cellx(cellcode,rcell)-cellzero.x;
  const int cy=PC__Celly(cellcode,rcell)-cellzero.y;
  const int cz=PC__Cellz(cellcode,rcell)-cellzero.z;
  // Ordena los componentes según el eje utilizado en el orden de la celda.
  const int c1=(axis==MGDIV_X? cy: cx);
  const int c2=(axis==MGDIV_Z? cy: cz);
  const int c3=(axis==MGDIV_Z? cz: (axis==MGDIV_X? cx: cy));
  // Código para scelldiv 1 o 2 pero no cero.
  ini1=c1-min(c1,scelldiv);
  fin1=c1+min(nc.x-c1-1,scelldiv)+1;
  ini2=c2-min(c2,scelldiv);
  fin2=c2+min(nc.y-c2-1,scelldiv)+1;
  ini3=c3-min(c3,scelldiv);
  fin3=c3+min(nc.z-c3-1,scelldiv)+1;
  // Computa los límites absolutos.
  ini3*=nc.w; fin3*=nc.w;
  ini2*=nc.x; fin2*=nc.x;
}

//------------------------------------------------------------------------------
/// Devuelve los límites de la celda para la interacción.
/// Returns cell limits for the interaction.
//------------------------------------------------------------------------------
__device__ void Initsp(const double &px,const double &py,const double &pz
  ,const unsigned &axis,const double3 &domposmin,const float &scell
  ,const int &scelldiv,const int4 &nc,const int3 &cellzero
  ,int &ini1,int &fin1,int &ini2,int &fin2,int &ini3,int &fin3)
{
  // Obtiene los límites de interacción.
  const int cx=int((px-domposmin.x)/scell)-cellzero.x;
  const int cy=int((py-domposmin.y)/scell)-cellzero.y;
  const int cz=int((pz-domposmin.z)/scell)-cellzero.z;
  // Ordena los componentes según el eje utilizado en el orden de la celda.
  const int c1=(axis==MGDIV_X? cy: cx);
  const int c2=(axis==MGDIV_Z? cy: cz);
  const int c3=(axis==MGDIV_Z? cz: (axis==MGDIV_X? cx: cy));
  // Código para scelldiv 1 o 2 pero no cero.
  ini1=c1-min(c1,scelldiv);
  fin1=c1+min(nc.x-c1-1,scelldiv)+1;
  ini2=c2-min(c2,scelldiv);
  fin2=c2+min(nc.y-c2-1,scelldiv)+1;
  ini3=c3-min(c3,scelldiv);
  fin3=c3+min(nc.z-c3-1,scelldiv)+1;
  // Computa los límites absolutos.
  ini3*=nc.w; fin3*=nc.w;
  ini2*=nc.x; fin2*=nc.x;
}


//------------------------------------------------------------------------------
/// Devuelve el rango de partículas límite para la búsqueda de vecinos.
/// Returns range of boundary particles for neighborhood search.
//------------------------------------------------------------------------------
__device__ void ParticleRange(const int &c2,const int &c3
  ,const int &ini1,const int &fin1,const int2 *begincell
  ,unsigned &pini,unsigned &pfin)
{
  const int v=c2+c3;
  for(int c1=ini1;c1<fin1;c1++){
    const int2 cbeg=begincell[c1+v];
    if(cbeg.y){
      if(!pfin)pini=cbeg.x;
      pfin=cbeg.y;
    }
  }
}

//==============================================================================
/// Devuelve la distancia al cuadrado entre partículas 1 y 2 (rr2).
/// Returns distance squared between particles 1 and 2 (rr2).
//==============================================================================
__device__ float Distance2(const float4 &pscellp1,const float4 &pscellp2
  ,const float &poscellsize)
{
  const float drx=pscellp1.x-pscellp2.x + poscellsize*(CEL_GetX(__float_as_int(pscellp1.w))-CEL_GetX(__float_as_int(pscellp2.w)));
  const float dry=pscellp1.y-pscellp2.y + poscellsize*(CEL_GetY(__float_as_int(pscellp1.w))-CEL_GetY(__float_as_int(pscellp2.w)));
  const float drz=pscellp1.z-pscellp2.z + poscellsize*(CEL_GetZ(__float_as_int(pscellp1.w))-CEL_GetZ(__float_as_int(pscellp2.w)));
  return(drx*drx + dry*dry + drz*drz);
}

}

#endif
