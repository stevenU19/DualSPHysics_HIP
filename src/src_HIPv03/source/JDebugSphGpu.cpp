//HEAD_DSCODES
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

#include "JDebugSphGpu.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JException.h"
#include "Functions.h"
#include "FunctionsHip.h"  // Cambiado de FunctionsCuda.h a FunctionsHip.h
#include "JVtkLib.h"
#include "JDataArrays.h"

#include "JSphGpuSingle.h"

#pragma warning(disable : 4996) // Cancela el aviso de deprecación de sprintf().

using namespace std;

//##############################################################################
//# JDebugSphGpu
//##############################################################################
//==============================================================================
/// Lanza una excepción relacionada con un archivo desde un método estático.
//==============================================================================
void JDebugSphGpu::RunExceptioonStatic(const std::string &srcfile, int srcline
  , const std::string &method
  , const std::string &msg, const std::string &file)
{
  throw JException(srcfile, srcline, "JDebugSphGpu", method, msg, file);
}
//==============================================================================
/// Lanza una excepción relacionada con un error HIP desde un método estático.
//==============================================================================
void JDebugSphGpu::RunExceptioonHipStatic(const std::string &srcfile, int srcline
  , const std::string &method
  , hipError_t hiperr, std::string msg)
{
  msg = msg + fun::PrintStr(" (HIP error %d (%s)).\n", hiperr, hipGetErrorString(hiperr));
  throw JException(srcfile, srcline, "JDebugSphGpu", method, msg, "");
}
//==============================================================================
/// Verifica el error de HIP y lanza una excepción desde un método estático.
//==============================================================================
void JDebugSphGpu::CheckHipErrorStatic(const std::string &srcfile, int srcline
  , const std::string &method, std::string msg)
{
  hipError_t hiperr = hipGetLastError();
  if (hiperr != hipSuccess) RunExceptioonHipStatic(srcfile, srcline, method, hiperr, msg);
}

//==============================================================================
/// Devuelve un puntero dinámico con el tipo de código de las partículas. (este puntero debe eliminarse)
//==============================================================================
byte* JDebugSphGpu::GetCodeType(unsigned n, const typecode *code) {
  byte *codetype = JDataArrays::NewArrayByte(n, false);
  for (unsigned p = 0; p < n; p++) {
    const typecode type = CODE_GetType(code[p]);
    codetype[p] = (type == CODE_TYPE_FIXED ? 0 : (type == CODE_TYPE_MOVING ? 1 : (type == CODE_TYPE_FLOATING ? 2 : (type == CODE_TYPE_FLUID ? 3 : 99))));
  }
  return(codetype);
}

//==============================================================================
/// Devuelve un puntero dinámico con el valor del tipo de código de las partículas. (este puntero debe eliminarse)
//==============================================================================
typecode* JDebugSphGpu::GetCodeTypeValue(unsigned n, const typecode *code) {
#ifdef CODE_SIZE4
  typecode *codetval = JDataArrays::NewArrayUint(n, false);
#else
  typecode *codetval = JDataArrays::NewArrayWord(n, false);
#endif
  for (unsigned c = 0; c < n; c++) codetval[c] = CODE_GetTypeValue(code[c]);
  return(codetval);
}

//==============================================================================
/// Devuelve un puntero dinámico con las coordenadas de la celda de las partículas. (este puntero debe eliminarse)
//==============================================================================
tuint3* JDebugSphGpu::GetCell3(unsigned n, const unsigned *dcell, unsigned cellcode) {
  tuint3 *cell3 = JDataArrays::NewArrayUint3(n, false);
  for (unsigned c = 0; c < n; c++) {
    const unsigned dcel = dcell[c];
    cell3[c] = TUint3(unsigned(PC__Cellx(cellcode, dcel)), unsigned(PC__Celly(cellcode, dcel)), unsigned(PC__Cellz(cellcode, dcel)));
  }
  return(cell3);
}

//==============================================================================
/// Devuelve un puntero dinámico con la posición como tfloat3. (este puntero debe eliminarse)
//==============================================================================
tfloat3* JDebugSphGpu::GetPosf3(unsigned n, const tdouble3 *pos) {
  tfloat3 *posf = JDataArrays::NewArrayFloat3(n, false);
  for (unsigned c = 0; c < n; c++) posf[c] = ToTFloat3(pos[c]);
  return(posf);
}

//==============================================================================
/// Devuelve un puntero dinámico con la posición como tfloat3. (este puntero debe eliminarse)
//==============================================================================
tfloat3* JDebugSphGpu::GetPosf3(unsigned n, const tdouble2 *posxy, const double *posz) {
  tfloat3 *posf = JDataArrays::NewArrayFloat3(n, false);
  for (unsigned c = 0; c < n; c++) posf[c] = TFloat3(float(posxy[c].x), float(posxy[c].y), float(posz[c]));
  return(posf);
}

//==============================================================================
/// Devuelve un puntero dinámico con la posición como tdouble3. (este puntero debe eliminarse)
//==============================================================================
tdouble3* JDebugSphGpu::GetPosd3(unsigned n, const tdouble2 *posxy, const double *posz) {
  tdouble3 *posd = JDataArrays::NewArrayDouble3(n, false);
  for (unsigned c = 0; c < n; c++) posd[c] = TDouble3(posxy[c].x, posxy[c].y, posz[c]);
  return(posd);
}

//==============================================================================
/// Devuelve un puntero dinámico con la posición relativa desde PosCell. (este puntero debe eliminarse)
//==============================================================================
tfloat3* JDebugSphGpu::GetPosCell_Pos(unsigned n, const tfloat4 *poscell) {
  tfloat3 *posf = JDataArrays::NewArrayFloat3(n, false);
  for (unsigned c = 0; c < n; c++) {
    const tfloat4 ps = poscell[c];
    posf[c] = TFloat3(ps.x, ps.y, ps.z);
  }
  return(posf);
}

//==============================================================================
/// Devuelve un puntero dinámico con las coordenadas de la celda desde PosCell. (este puntero debe eliminarse)
//==============================================================================
tuint3* JDebugSphGpu::GetPosCell_Cell(unsigned n, const tfloat4 *poscell) {
  const tuint4* poscellu = (const tuint4*)poscell;
  tuint3 *cell = JDataArrays::NewArrayUint3(n, false);
  for (unsigned c = 0; c < n; c++) {
    const unsigned cellu = poscellu[c].w;
    cell[c] = TUint3(CEL_GetX(cellu), CEL_GetY(cellu), CEL_GetZ(cellu));
  }
  return(cell);
}

//==============================================================================
/// Prepara la lista de variables para su verificación.
//==============================================================================
std::string JDebugSphGpu::PrepareVars(const std::string &vlist) {
  return(string(",") + fun::StrLower(vlist) + ",");
}

//==============================================================================
/// Verifica la lista de variables y devuelve la variable desconocida.
//==============================================================================
std::string JDebugSphGpu::CheckVars(std::string vlist) {
  const string allvars = ",all,idp,gid,seq,cell,code,vel,rhop,ace,ar,poscell,";
  vlist = fun::StrLower(vlist);
  while (!vlist.empty()) {
    string var = fun::StrSplit(",", vlist);
    if (!var.empty() && !FindVar(var, allvars)) return(var);
  }
  return("");
}

//==============================================================================
/// Devuelve el nombre del archivo formateado.
//==============================================================================
std::string JDebugSphGpu::GetFileName(std::string filename, int numfile, int gid) {
  const string fext = fun::GetExtension(filename);
  string file = AppInfo.GetDirOut() + fun::GetWithoutExtension(filename) + (gid >= 0 ? fun::PrintStr("_g%d", gid) : "") + "." + fext;
  if (numfile >= 0) file = fun::FileNameSec(file, numfile);
  return(file);
}

//==============================================================================
/// Carga los datos de las partículas en el objeto ffdata.
//==============================================================================
void JDebugSphGpu::LoadParticlesData(const JSphGpuSingle *gp, unsigned pini, unsigned pfin
  , std::string vars, JDataArrays *arrays, std::string file)
{
  if (pini >= pfin) Run_ExceptioonFileSta(fun::PrintStr("Datos inválidos (pini:%u >= pfin:%u).", pini, pfin), file);
  const unsigned n = pfin - pini;

  vars = PrepareVars(vars);
  string errvar = CheckVars(vars);
  if (!errvar.empty()) Run_ExceptioonFileSta(fun::PrintStr("La variable \'%s\' es desconocida.", errvar.c_str()), file);
  const bool all = FindVar("all", vars);

  arrays->Reset();
  if (all || FindVar("idp", vars)) {
    const unsigned *idp = fhip::ToHostUint(pini, n, gp->Idpg); // Cambiado de fcuda a fhip
    arrays->AddArray("Idp", n, idp, true);
  }
  if (all || FindVar("seq", vars)) {
    arrays->AddArray("Seq", n, JDataArrays::NewArraySeqUint(n, pini, 1), true);
  }
  if (all || FindVar("cell", vars)) {
    const unsigned *dcell = fhip::ToHostUint(pini, n, gp->Dcellg);  // Cambiado de fcuda a fhip
    arrays->AddArray("Cell", n, GetCell3(n, dcell, gp->DomCellCode), true);
    delete[] dcell; dcell = NULL;
  }
  if (all || FindVar("vel", vars) || FindVar("rhop", vars)) {
    const tfloat4 *velrhop = fhip::ToHostFloat4(pini, n, gp->Velrhopg);  // Cambiado de fcuda a fhip
    if (all || FindVar("vel", vars)) arrays->AddArray("Vel", n, JDataArrays::NewArrayFloat3xyz(n, velrhop), true);
    if (all || FindVar("rhop", vars)) arrays->AddArray("Rhop", n, JDataArrays::NewArrayFloat1w(n, velrhop), true);
    delete[] velrhop; velrhop = NULL;
  }
  if (all || FindVar("code", vars)) {
#ifdef CODE_SIZE4
    const typecode *code = fhip::ToHostUint(pini, n, gp->Codeg);  // Cambiado de fcuda a fhip
#else
    const typecode *code = fhip::ToHostWord(pini, n, gp->Codeg);  // Cambiado de fcuda a fhip
#endif
    arrays->AddArray("Code", n, code, true);
    arrays->AddArray("Type", n, GetCodeType(n, code), true);
    arrays->AddArray("TypeValue", n, GetCodeTypeValue(n, code), true);
  }
  {
    const tdouble2 *posxy = fhip::ToHostDouble2(pini, n, gp->Posxyg);  // Cambiado de fcuda a fhip
    const double   *posz = fhip::ToHostDouble(pini, n, gp->Poszg);  // Cambiado de fcuda a fhip
    if (FindVar("posd", vars)) arrays->AddArray("Pos", n, GetPosd3(n, posxy, posz), true);
    else arrays->AddArray("Pos", n, GetPosf3(n, posxy, posz), true);
    delete[] posxy; posxy = NULL;
    delete[] posz; posz = NULL;
  }
  if (all || FindVar("poscell", vars)) {
    const float4 *ptrg = gp->PosCellg;
    if (ptrg) {
      const tfloat4 *poscell = fhip::ToHostFloat4(pini, n, ptrg);  // Cambiado de fcuda a fhip
      arrays->AddArray("poscell_Pos", n, GetPosCell_Pos(n, poscell), true);
      arrays->AddArray("poscell_Cell", n, GetPosCell_Cell(n, poscell), true);
      delete[] poscell; poscell = NULL;
    }
    else if (!all) Run_ExceptioonFileSta("La variable PosCellg es NULL.", file);
  }
  if (all || FindVar("ace", vars)) {
    const float3 *ptrg = gp->Aceg;
    if (ptrg) arrays->AddArray("Ace", n, fhip::ToHostFloat3(pini, n, ptrg), true);  // Cambiado de fcuda a fhip
    else if (!all) Run_ExceptioonFileSta("La variable Aceg es NULL.", file);
  }
  if (all || FindVar("ar", vars)) {
    const float *ptrg = gp->Arg;
    if (ptrg) arrays->AddArray("Ar", n, fhip::ToHostFloat(pini, n, ptrg), true);  // Cambiado de fcuda a fhip
    else if (!all) Run_ExceptioonFileSta("La variable Arg es NULL.", file);
  }
}

//==============================================================================
/// Almacena los datos en formato VTK.
//============================================================================== 
void JDebugSphGpu::SaveVtk(std::string filename, int numfile, unsigned pini, unsigned pfin
  , std::string vars, const JSphGpuSingle *gp)
{
  const string file = GetFileName(filename, numfile);
  JDataArrays arrays;
  LoadParticlesData(gp, pini, pfin, vars, &arrays, file);
  JVtkLib::SaveVtkData(file, arrays, "Pos");
  arrays.Reset();
}
