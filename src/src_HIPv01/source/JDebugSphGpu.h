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

/// \file JDebugSphGpu.h \brief Declares the class \ref JDebugSphGpu.

#ifndef _JDebugSphGpu_
#define _JDebugSphGpu_

#include "TypesDef.h"
#include "DualSphDef.h"

#include <string>
#include <cstring>
#include <hip/hip_runtime_api.h> // Cambiado de CUDA a HIP

//-Defines for normal exceptions for static methods.
#ifndef Run_ExceptionSta
#define Run_ExceptionSta(msg) RunExceptionStatic(__FILE__,__LINE__,__func__,msg)
#endif
#ifndef Run_ExceptionFileSta
#define Run_ExceptionFileSta(msg,file) RunExceptionStatic(__FILE__,__LINE__,__func__,msg,file)
#endif
//-Defines for HIP exceptions for static methods.
#ifndef Run_ExceptionHipSta
#define Run_ExceptionHipSta(hiperr,msg) RunExceptionHipStatic(__FILE__,__LINE__,__func__,hiperr,msg)
#endif
#ifndef Check_HipErrorSta
#define Check_HipErrorSta(msg) CheckHipErrorStatic(__FILE__,__LINE__,__func__,msg)
#endif

class JDataArrays;
class JSphGpuSingle;
class JSphMgpuNodeUnit;
class JSphDomain;

//##############################################################################
//# JDebugSphGpu
//##############################################################################
/// \brief Proporciona funciones para almacenar datos de partículas en formatos VTK, CSV, ASCII.

class JDebugSphGpu
{
protected:
  static void RunExceptioonStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,const std::string &msg,const std::string &file="");
  static void RunExceptionHipStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,hipError_t hiperr,std::string msg); // Cambiado cudaError_t a hipError_t
  static void CheckHipErrorStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,std::string msg);

public:

  static byte*     GetCodeType     (unsigned n,const typecode *code);
  static typecode* GetCodeTypeValue(unsigned n,const typecode *code);
  static tuint3*   GetCell3(unsigned n,const unsigned *dcell,unsigned cellcode);
  static tfloat3*  GetPosf3(unsigned n,const tdouble3 *pos);
  static tfloat3*  GetPosf3(unsigned n,const tdouble2 *posxy,const double *posz);
  static tdouble3* GetPosd3(unsigned n,const tdouble2 *posxy,const double *posz);
  static tfloat3*  GetPosCell_Pos (unsigned n,const tfloat4 *poscell);
  static tuint3*   GetPosCell_Cell(unsigned n,const tfloat4 *poscell);

  static std::string PrepareVars(const std::string &vlist);
  static bool FindVar(const std::string &var,const std::string &vlist){ return(int(vlist.find(std::string(",")+var+","))>=0); }
  static std::string CheckVars(std::string vlist);

  static std::string GetFileName(std::string filename,int numfile,int gid=-1);

  static void LoadParticlesData(const JSphGpuSingle *gp,unsigned pini,unsigned pfin,std::string vars,JDataArrays *arrays,std::string file="");
  static void SaveVtk(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string vars,const JSphGpuSingle *gp);

};


#endif
