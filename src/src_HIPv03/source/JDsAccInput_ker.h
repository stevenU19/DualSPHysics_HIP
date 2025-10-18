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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Implementación en fichero independiente. (22-12-2019)
//:# - Cambio de nombre de J.SphAccInput a J.DsAccInput. (28-06-2020)
//:#############################################################################

/// \file JDsAccInput_ker.h \brief Declara funciones y kernels de HIP para fuerzas externas (JDsAccInput) en la GPU.

#ifndef _JDsAccInput_ker_
#define _JDsAccInput_ker_

#include "DualSphDef.h"
#include <hip/hip_runtime_api.h>  // Cambiado de CUDA a HIP

/// Implementa un conjunto de funciones y kernels de HIP para fuerzas externas (JDsAccInput) en la GPU.
namespace cuaccin{

//-Kernels para fuerzas externas (JDsAccInput).
void AddAccInput(unsigned n,unsigned pini,typecode codesel1,typecode codesel2
  ,tdouble3 acclin,tdouble3 accang,tdouble3 centre,tdouble3 velang,tdouble3 vellin,bool setgravity
  ,tfloat3 gravity,const typecode *code,const double2 *posxy,const double *posz
  ,const float4 *velrhop,float3 *ace,hipStream_t stm);  // Cambiado cudaStream_t a hipStream_t

}

#endif
