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
//:# - GetAccValues() guarda entrada y salida para evitar cálculos con llamadas 
//:#   consecutivas iguales. (13-09-2019)
//:# - Objeto JXml pasado como const para operaciones de lectura. (18-03-2020)  
//:# - Comprueba opción active en elementos de primer y segundo nivel. (18-03-2020)  
//:# - Permite definir un intervalo de tiempo para su activación. (26-04-2020)  
//:# - Cambio de nombre de J.SphAccInput a J.DsAccInput. (28-06-2020)
//:#############################################################################

/// \file JDsAccInput.h \brief Declares the class \ref JDsAccInput.

#ifndef _JDsAccInput_
#define _JDsAccInput_

#include "JObject.h"
#include "DualSphDef.h"
#ifdef _WITHGPU
  #include <hip/hip_runtime_api.h> // Cambiado de CUDA a HIP
#endif

#include <string>
#include <vector>

class JLog2;
class JXml;
class TiXmlElement;
class JLinearValue;
class JSphMk;

//##############################################################################
//# JDsAccInputMk
//##############################################################################
/// \brief Proporciona la fuerza que se aplicará a diferentes bloques de partículas cargadas desde archivos.

class JDsAccInputMk : protected JObject
{
public:
  const unsigned Idx;     ///<Índice de la configuración.
  const bool Bound;       ///<Tipo de partículas objetivo (frontera, flotante o fluido).
  const word MkType1;     ///<El tipo MK frontera o fluido para seleccionar las partículas objetivo.
  const word MkType2;     ///<Rango final de valores MK para seleccionar las partículas objetivo.
  const double TimeIni;   ///<Tiempo inicial de activación (por defecto = 0).
  const double TimeEnd;   ///<Tiempo final de activación (por defecto = DBL_MAX).
  const bool GravityEnabled;  ///<Determina si la gravedad global está habilitada o deshabilitada para este conjunto de partículas.
  const tfloat3 AccCoG;       ///<Centro de gravedad que se utilizará para los cálculos de aceleración angular.

protected:
  JLog2* Log;

  typecode CodeSel1; ///<Primer código de las partículas objetivo (Tipo y Valor).
  typecode CodeSel2; ///<Último código de las partículas objetivo (Tipo y Valor).

  JLinearValue *AceData; ///<Datos de aceleración de entrada.
  JLinearValue *VelData; ///<Datos de velocidad de entrada.

  double LastTimestepInput; ///<Guarda el último valor utilizado con GetAccValues().
  StAceInput LastOutput;    ///<Guarda el último valor devuelto por GetAccValues().

  void Reset();

public:
  JDsAccInputMk(JLog2* log,unsigned idx,bool bound,word mktype1,word mktype2
    ,double tini,double tend,bool genabled,tfloat3 acccentre
    ,const JLinearValue &acedata,const JLinearValue &veldata);
  ~JDsAccInputMk();
  long long GetAllocMemory()const;

  void ConfigCodeSel(typecode codesel1,typecode codesel2){ 
    CodeSel1=codesel1; CodeSel2=codesel2; 
  }
  typecode GetCodeSel1()const{ return(CodeSel1); };
  typecode GetCodeSel2()const{ return(CodeSel2); };

  void GetConfig(std::vector<std::string> &lines)const;

  const StAceInput& GetAccValues(double timestep);
};

//##############################################################################
//# JDsAccInput
//##############################################################################
/// \brief Administra la aplicación de fuerzas externas a diferentes bloques de partículas (con el mismo MK).

class JDsAccInput : protected JObject
{
protected:
  JLog2* Log;
  std::string DirData;
  std::vector<JDsAccInputMk*> Inputs;
  long long MemSize;

  void Reset();
  bool ExistMk(bool bound,word mktype)const;
  void LoadXml(const JXml *sxml,const std::string &place);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void ComputeVelocity(const JLinearValue &acedata,JLinearValue &veldata)const;

public:
  JDsAccInput(JLog2* log,const std::string &dirdata,const JXml *sxml,const std::string &place);
  ~JDsAccInput();
  long long GetAllocMemory()const{ return(MemSize); }

  void VisuConfig(std::string txhead,std::string txfoot)const;
  void Init(const JSphMk *mkinfo);

  unsigned GetCount()const{ return(unsigned(Inputs.size())); };
  const StAceInput& GetAccValues(unsigned cinput,double timestep);

  void RunCpu(double timestep,tfloat3 gravity,unsigned n,unsigned pini
    ,const typecode *code,const tdouble3 *pos,const tfloat4 *velrhop,tfloat3 *ace);

#ifdef _WITHGPU
  void RunGpu(double timestep,tfloat3 gravity,unsigned n,unsigned pini
    ,const typecode *code,const double2 *posxy,const double *posz,const float4 *velrhop,float3 *ace);
#endif
};

#endif
