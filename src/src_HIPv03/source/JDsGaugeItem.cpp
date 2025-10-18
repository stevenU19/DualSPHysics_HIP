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

#include "JDsGaugeItem.h"
#include "JCellSearch_inline.h"
#include "FunSphKernel.h"
#include "FunSphEos.h"
#include "JException.h"
#include "JLog2.h"
#include "JSaveCsv2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsGeo3d.h"
#include "JDataArrays.h"
#include "JVtkLib.h"
#ifdef _WITHGPU
  #include "FunctionsHip.h"  // Cambiado de FunctionsCuda.h a FunctionsHip.h
  #include "JDsGauge_ker.h"
  #include "JReduSum_ker.h"
#endif
#include <cmath>
#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;

//##############################################################################
//# JGaugeItem
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeItem::JGaugeItem(TpGauge type,unsigned idx,std::string name,bool cpu,JLog2* log)
  :Type(type),Idx(idx),Name(name),Cpu(cpu),Log(log)
{
  ClassName="JGaugeItem";
  Reset();
}

#ifdef _WITHGPU
//==============================================================================
/// Lanza una excepción relacionada con un error HIP desde un método estático.
//==============================================================================
void JGaugeItem::RunExceptioonHip(const std::string &srcfile,int srcline
  ,const std::string &classname,const std::string &method
  ,hipError_t hiperr,std::string msg)const
{
  msg=msg+fun::PrintStr(" (HIP error %d (%s)).\n",hiperr,hipGetErrorString(hiperr));
  throw JException(srcfile,srcline,classname,method,msg,"");
}

//==============================================================================
/// Verifica el error de HIP y lanza una excepción si la hay.
//==============================================================================
void JGaugeItem::CheckHipError(const std::string &srcfile,int srcline
  ,const std::string &classname,const std::string &method
  ,std::string msg)const
{
  hipError_t hiperr=hipGetLastError();
  if(hiperr!=hipSuccess)RunExceptioonHip(srcfile,srcline,classname,method,hiperr,msg);
}
#endif

//==============================================================================
/// Inicialización de variables.
//==============================================================================
void JGaugeItem::Reset(){
  Config(CteSphNull(),false,TDouble3(0),TDouble3(0),0,0);
  SaveVtkPart=false;
  ConfigComputeTiming(0,0,0);
  ConfigOutputTiming(false,0,0,0);
  TimeStep=0;
  OutCount=0;
  OutFile="";
}

//==============================================================================
/// Configuración del objeto.
//==============================================================================
void JGaugeItem::Config(const StCteSph & csp,bool symmetry,tdouble3 domposmin
    ,tdouble3 domposmax,float scell,int scelldiv)
{
  CSP=csp;
  Symmetry=symmetry;
  DomPosMin=domposmin;
  DomPosMax=domposmax;
  Scell=scell;
  ScellDiv=scelldiv;
}

//==============================================================================
/// Configuración de tiempos de cálculo.
//==============================================================================
void JGaugeItem::ConfigComputeTiming(double start,double end,double dt){
  ComputeDt=dt;
  ComputeStart=start;
  ComputeEnd=end;
  ComputeNext=0;
}

//==============================================================================
/// Configuración de tiempos de salida.
//==============================================================================
void JGaugeItem::ConfigOutputTiming(bool save,double start,double end,double dt){
  OutputSave=save;
  OutputDt=dt;
  OutputStart=start;
  OutputEnd=end;
  OutputNext=0;
}

//==============================================================================
/// Retorna el tipo en formato de texto.
//==============================================================================
std::string JGaugeItem::GetNameType(TpGauge type){
  switch(type){
    case GAUGE_Vel:   return("Vel");
    case GAUGE_Swl:   return("SWL");
    case GAUGE_MaxZ:  return("MaxZ");
    case GAUGE_Force: return("Force");
  }
  return("???");
}

//==============================================================================
/// Obtiene las líneas con la configuración de la sonda.
//==============================================================================
void JGaugeItem::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back(fun::PrintStr("Type.......: %s",JGaugeItem::GetNameType(Type).c_str()));
  const string cpend=fun::DoublexStr(ComputeEnd,"%g");
  lines.push_back(fun::PrintStr("Compute....: %g - %s   dt:%g",ComputeStart,cpend.c_str(),ComputeDt));
  const string ouend=fun::DoublexStr(OutputEnd,"%g");
  lines.push_back(fun::PrintStr("Output.....:%s%g - %s   dt:%g",(OutputSave? " ": " (disabled)  "),OutputStart ,ouend.c_str(),OutputDt));
  lines.push_back(fun::PrintStr("SaveVtkPart: %s",(SaveVtkPart? "True": "False")));
  if(Type==GAUGE_Vel){
    const JGaugeVelocity* gau=(JGaugeVelocity*)this;
    lines.push_back(fun::PrintStr("Point......: (%g,%g,%g)",gau->GetPoint().x,gau->GetPoint().y,gau->GetPoint().z));
  }
  else if(Type==GAUGE_Swl){
    const JGaugeSwl* gau=(JGaugeSwl*)this;
    lines.push_back(fun::PrintStr("MassLimit..: %g",gau->GetMassLimit()));
    lines.push_back(fun::PrintStr("GaugePoints: %s",fun::Double3gRangeStr(gau->GetPoint0(),gau->GetPoint2()).c_str()));
    lines.push_back(fun::PrintStr("PointDp....: %g",gau->GetPointDp()));
  }
  else if(Type==GAUGE_MaxZ){
    const JGaugeMaxZ* gau=(JGaugeMaxZ*)this;
    lines.push_back(fun::PrintStr("Point0.....: (%g,%g,%g)   Height:%g",gau->GetPoint0().x,gau->GetPoint0().y,gau->GetPoint0().z,gau->GetHeight()));
    lines.push_back(fun::PrintStr("DistLimit..: %g",gau->GetDistLimit()));
  }
  else if(Type==GAUGE_Force){
    const JGaugeForce* gau=(JGaugeForce*)this;
    lines.push_back(fun::PrintStr("MkBound.....: %u (%s particles)",gau->GetMkBound(),TpPartGetStrCode(gau->GetTypeParts())));
    lines.push_back(fun::PrintStr("Particles id: %u - %u",gau->GetIdBegin(),gau->GetIdBegin()+gau->GetCount()-1));
  }
  else Run_Exceptioon("Type unknown.");
}

//==============================================================================
/// Guarda los resultados en archivo CSV y/o VTK.
//==============================================================================
void JGaugeItem::SaveResults(unsigned cpart){
  SaveResults();
  if(SaveVtkPart)SaveVtkResult(cpart);
}

//==============================================================================
/// Calcula la velocidad en puntos indicados (en GPU con HIP).
//==============================================================================
#ifdef _WITHGPU
void JGaugeVelocity::CalculeGpu(double timestep,const StDivDataGpu &dvd
  ,unsigned npbok,unsigned npb,unsigned np,const double2 *posxy,const double *posz
  ,const typecode *code,const unsigned *idp,const float4 *velrhop,float3 *aux)
{
  SetTimeStep(timestep);
  //-Inicio de la medida.
  tfloat3 ptvel=TFloat3(0);
  const bool ptout=PointIsOut(Point.x,Point.y,Point.z);//-Comprueba que el punto esté dentro de los límites del dominio.
  if(!ptout){
    cugauge::Interaction_GaugeVel(CSP,dvd,Point,posxy,posz,code,velrhop,aux);
    hipMemcpy(&ptvel,aux,sizeof(float3),hipMemcpyDeviceToHost);  // Cambiado cudaMemcpy a hipMemcpy
    CheckHipError("Failed in velocity calculation.");
  }
  //-Guarda el resultado.
  Result.Set(timestep,ToTFloat3(Point),ptvel);
  if(Output(timestep))StoreResult();
}
#endif

//==============================================================================
/// Calcula el nivel de agua superficial en los puntos indicados (en GPU con HIP).
//==============================================================================
#ifdef _WITHGPU
void JGaugeSwl::CalculeGpu(double timestep,const StDivDataGpu &dvd
  ,unsigned npbok,unsigned npb,unsigned np,const double2 *posxy,const double *posz
  ,const typecode *code,const unsigned *idp,const float4 *velrhop,float3 *aux)
{
  SetTimeStep(timestep);
  //-Inicio de la medida.
  cugauge::Interaction_GaugeSwl(CSP,dvd,Point0,PointDir,PointNp,MassLimit
    ,posxy,posz,code,velrhop,aux);
  tfloat3 ptsurf=TFloat3(0);
  hipMemcpy(&ptsurf,aux,sizeof(float3),hipMemcpyDeviceToHost);  // Cambiado cudaMemcpy a hipMemcpy
  CheckHipError("Failed in Swl calculation.");
  //-Guarda el resultado.
  Result.Set(timestep,ToTFloat3(Point0),ToTFloat3(Point2),ptsurf);
  if(Output(timestep))StoreResult();
}
#endif

//==============================================================================
/// Calcula la fuerza sumatoria en partículas seleccionadas usando solo partículas fluidas (en GPU con HIP).
//==============================================================================
#ifdef _WITHGPU
void JGaugeForce::CalculeGpu(double timestep,const StDivDataGpu &dvd
  ,unsigned npbok,unsigned npb,unsigned np,const double2 *posxy,const double *posz
  ,const typecode *code,const unsigned *idp,const float4 *velrhop,float3 *aux)
{
  SetTimeStep(timestep);
  //-Inicializa el array de aceleración a cero.
  hipMemset(PartAceg,0,sizeof(float3)*Count);  // Cambiado cudaMemset a hipMemset
  const int n=int(TypeParts==TpPartFixed || TypeParts==TpPartMoving? npbok: np);
  //-Calcula la aceleración en partículas límite seleccionadas.
  cugauge::Interaction_GaugeForce(CSP,dvd,n,IdBegin,Code
    ,posxy,posz,code,idp,velrhop,PartAceg);

  //-Calcula la aceleración total.
  tfloat3 acesum=TFloat3(0);
  const float3 result=curedus::ReduSumFloat3(Count,0,PartAceg,Auxg);
  acesum=TFloat3(result.x,result.y,result.z);
  CheckHipError("Failed in Force calculation.");

  //-Guarda el resultado.
  Result.Set(timestep,acesum*CSP.massbound);
  if(Output(timestep))StoreResult();
}
#endif
