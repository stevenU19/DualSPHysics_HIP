//HEAD_DSPH
/*
 <DUALSPHYSICS> Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/).

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 Este archivo es parte de DualSPHysics.

 DualSPHysics es software libre: puedes redistribuirlo y/o modificarlo bajo los términos de la GNU Lesser General Public License 
 publicada por la Free Software Foundation, ya sea la versión 2.1 de la licencia, o (a tu elección) cualquier versión posterior.

 DualSPHysics se distribuye con la esperanza de que sea útil, pero SIN NINGUNA GARANTÍA; sin siquiera la garantía implícita de
 COMERCIABILIDAD o APTITUD PARA UN PROPÓSITO PARTICULAR. Consulta la GNU Lesser General Public License para más detalles.

 Debiste haber recibido una copia de la GNU Lesser General Public License junto con DualSPHysics. Si no es así, consulta <http://www.gnu.org/licenses/>.
*/

/// \file JDsGauge_ker.h \brief Declara funciones y kernels de HIP para las clases JGauge.

#ifndef _JDsGauge_ker_
#define _JDsGauge_ker_

#include "DualSphDef.h"
#include "JCellDivDataGpu.h"
#include <hip/hip_runtime_api.h>  // Cambiado de CUDA a HIP

/// Implementa un conjunto de funciones y kernels de HIP para las clases que gestionan gauges.
namespace cugauge{

//-Kernel para JGaugeVelocity.
void Interaction_GaugeVel(const StCteSph &CSP,const StDivDataGpu &dvd
  ,tdouble3 ptpos,const double2 *posxy,const double *posz
  ,const typecode *code,const float4 *velrhop,float3 *ptvel);

//-Kernel para JGaugeSwl.
void Interaction_GaugeSwl(const StCteSph &CSP,const StDivDataGpu &dvd
  ,tdouble3 point0,tdouble3 pointdir,unsigned pointnp,float masslimit
  ,const double2 *posxy,const double *posz,const typecode *code
  ,const float4 *velrhop,float3 *ptres);

//-Kernel para JGaugeMaxZ.
void Interaction_GaugeMaxz(tdouble3 point0,float maxdist2,const StDivDataGpu &dvd
  ,int cxini,int cxfin,int yini,int yfin,int zini,int zfin
  ,const double2 *posxy,const double *posz,const typecode *code
  ,float3 *ptres);

//-Kernel para JGaugeForce.
void Interaction_GaugeForce(const StCteSph &CSP,const StDivDataGpu &dvd
  ,unsigned n,unsigned idbegin,typecode codesel,const double2 *posxy,const double *posz
  ,const typecode *code,const unsigned *idp,const float4 *velrhop
  ,float3 *partace);

}

#endif
