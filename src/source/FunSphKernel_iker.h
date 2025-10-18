//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is parte de DualSPHysics.

 DualSPHysics es software libre: puedes redistribuirlo y/o modificarlo bajo los términos de la GNU Lesser General Public License 
 publicada por la Free Software Foundation, ya sea la versión 2.1 de la licencia, o (a tu elección) cualquier versión posterior.

 DualSPHysics se distribuye con la esperanza de que sea útil, pero SIN NINGUNA GARANTÍA; sin siquiera la garantía implícita de
 COMERCIABILIDAD o APTITUD PARA UN PROPÓSITO PARTICULAR. Consulta la GNU Lesser General Public License para más detalles.

 Debiste haber recibido una copia de la GNU Lesser General Public License junto con DualSPHysics. Si no es así, consulta <http://www.gnu.org/licenses/>.
*/

/// \file FunSphKernel_iker.h \brief Implementa funciones de dispositivo HIP para los kernels de SPH.

#include "TypesDef.h"
#include "DualSphDef.h"
#include <hip/hip_runtime.h>

/// Implementa funciones de dispositivo HIP para los kernels de SPH.
namespace cufsph{

//##############################################################################
//# Kernel Cubic Spline
//##############################################################################
//#define CTE_AVAILABLE
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Devuelve wab del kernel.
//------------------------------------------------------------------------------
__device__ float GetKernelCubic_Wab(float rr2){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  if(rad>CTE.kernelh){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    return(CTE.cubic_a24*(wqq2*wqq1));
  }
  else{
    const float wqq2=qq*qq;
    return(CTE.cubic_a2*(1.0f+(0.75f*qq-1.5f)*wqq2));
  }
}
//------------------------------------------------------------------------------
/// Devuelve fac del kernel.
//------------------------------------------------------------------------------
__device__ float GetKernelCubic_Fac(float rr2){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  if(rad>CTE.kernelh){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    return(CTE.cubic_c2*wqq2/rad);
  }
  else{
    const float wqq2=qq*qq;
    return((CTE.cubic_c1*qq+CTE.cubic_d1*wqq2)/rad);
  }
}
//------------------------------------------------------------------------------
/// Devuelve wab y fac del kernel.
//------------------------------------------------------------------------------
__device__ float GetKernelCubic_WabFac(float rr2,float &fac){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  if(rad>CTE.kernelh){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    fac=(CTE.cubic_c2*wqq2/rad);
    return(CTE.cubic_a24*(wqq2*wqq1));
  }
  else{
    float wqq2=qq*qq;
    fac=((CTE.cubic_c1*qq+CTE.cubic_d1*wqq2)/rad);
    return(CTE.cubic_a2*(1.0f+(0.75f*qq-1.5f)*wqq2));
  }
}
//------------------------------------------------------------------------------
/// Devuelve la corrección tensil del kernel Cubic.
//------------------------------------------------------------------------------
__device__ float GetKernelCubic_Tensil(float rr2
  ,float rhopp1,float pressp1,float rhopp2,float pressp2)
{
  const float wab=GetKernelCubic_Wab(rr2);
  float fab=wab*CTE.cubic_odwdeltap;
  fab*=fab; fab*=fab; //fab=fab^4
  const float tensilp1=(pressp1/(rhopp1*rhopp1))*(pressp1>0? 0.01f: -0.2f);
  const float tensilp2=(pressp2/(rhopp2*rhopp2))*(pressp2>0? 0.01f: -0.2f);
  return(fab*(tensilp1+tensilp2));
}
#endif

//##############################################################################
//# Kernel Wendland
//##############################################################################
//------------------------------------------------------------------------------
/// Devuelve wab del kernel.
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_Wab(float rr2,float h,float awen){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  const float wqq=qq+qq+1.f;
  const float wqq1=1.f-0.5f*qq;
  const float wqq2=wqq1*wqq1;
  return(awen*wqq*wqq2*wqq2);
}
//------------------------------------------------------------------------------
/// Devuelve fac del kernel.
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_Fac(float rr2,float h,float bwen){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  const float wqq1=1.f-0.5f*qq;
  return(bwen*qq*wqq1*wqq1*wqq1/rad);
}

#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Devuelve wab y fac del kernel.
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_WabFac(float rr2,float &fac){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  const float wqq1=1.f-0.5f*qq;
  const float wqq2=wqq1*wqq1;
  fac=CTE.bwen*qq*wqq2*wqq1/rad;
  const float wqq=qq+qq+1.f;
  return(CTE.awen*wqq*wqq2*wqq2);
}
//------------------------------------------------------------------------------
/// Devuelve wab del kernel.
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_Wab(float rr2){ 
  return(GetKernelWendland_Wab(rr2,CTE.kernelh,CTE.awen)); 
}
//------------------------------------------------------------------------------
/// Devuelve fac del kernel.
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_Fac(float rr2){
  return(GetKernelWendland_Fac(rr2,CTE.kernelh,CTE.bwen));
}
#endif

//##############################################################################
//# Computa los valores del kernel usando plantillas.
//##############################################################################
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Devuelve wab del kernel de acuerdo con la plantilla.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_Wab(float rr2){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Wab  (rr2));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_Wab     (rr2));
  else return(0);
}
//------------------------------------------------------------------------------
/// Devuelve fac del kernel de acuerdo con la plantilla.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_Fac(float rr2){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Fac  (rr2));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_Fac     (rr2));
  else return(0);
}
//------------------------------------------------------------------------------
/// Devuelve wab y fac del kernel de acuerdo con la plantilla.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_WabFac(float rr2,float &fac){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_WabFac  (rr2,fac));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_WabFac     (rr2,fac));
  else return(0);
}
#endif
//------------------------------------------------------------------------------
/// Devuelve wab del kernel de acuerdo con la plantilla.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_Wab(float rr2,float h,float aker){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Wab  (rr2,h,aker));
  else return(0);
}
//------------------------------------------------------------------------------
/// Devuelve wab del kernel de acuerdo con la plantilla.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_Fac(float rr2,float h,float bker){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Fac  (rr2,h,bker));
  else return(0);
}

}
