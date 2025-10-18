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

/// \file JCellDivGpuSingle_ker.h \brief Declara funciones y kernels HIP para calcular operaciones de la lista de vecinos.

#ifndef _JCellDivGpuSingle_ker_
#define _JCellDivGpuSingle_ker_

#include "JCellDivGpu_ker.h"

namespace hipdiv{

void PreSortFull(unsigned np,unsigned cellcode,const unsigned *dcell,const typecode *code,tuint3 cellmin,tuint3 ncells,unsigned *cellpart,unsigned *sortpart);
void PreSortFluid(unsigned npf,unsigned pini,unsigned cellcode,const unsigned *dcell,const typecode *code,tuint3 cellmin,tuint3 ncells,unsigned *cellpart,unsigned *sortpart);

}

#endif
