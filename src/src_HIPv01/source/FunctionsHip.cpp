//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file FunctionsCuda.cpp \brief Implements basic/general GPU functions for the entire application.

#include "FunctionsHip.h"
#include "Functions.h"
#include <algorithm>
#include <hip/hip_runtime.h>  // Reemplazamos CUDA por HIP

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

namespace fhip{  // Puedes cambiar el namespace si lo necesitas

//==============================================================================
/// Checks error and throws exception.
//==============================================================================
void CheckHipErrorFun(const char *const file,int const line,const char *const fun
  ,std::string msg)
{
  const hipError_t hiperr=hipGetLastError();
  if(hiperr!=hipSuccess){
    msg=msg+fun::PrintStr(" (HIP error %d (%s)).\n",hiperr,hipGetErrorString(hiperr)); 
    fun::RunExceptioonFun(file,line,fun,msg);
  }
}

//==============================================================================
/// Returns information about selected GPU (code from deviceQuery example).
//==============================================================================
inline bool IsGPUCapableP2P(const hipDeviceProp_t *pProp){
    return(pProp->major >= 2);
}

//==============================================================================
/// Returns name about selected GPU.
//==============================================================================
std::string GetHipDeviceName(int gid){
  hipSetDevice(gid);
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "Failed selecting device.");
  hipDeviceProp_t deviceProp;
  hipGetDeviceProperties(&deviceProp,gid);
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "Failed getting selected device info.");
  return(deviceProp.name);
}

//==============================================================================
/// Returns information about selected GPU (code from deviceQuery example).
//==============================================================================
StGpuInfo GetHipDeviceInfo(int gid){
  hipSetDevice(gid);
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "Failed selecting device.");
  hipDeviceProp_t deviceProp;
  hipGetDeviceProperties(&deviceProp,gid);
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "Failed getting selected device info.");
  StGpuInfo g;
  g.id=gid;
  g.name=deviceProp.name;
  g.ccmajor=deviceProp.major;
  g.ccminor=deviceProp.minor;
  g.globalmem=deviceProp.totalGlobalMem;
  g.mp=deviceProp.multiProcessorCount;
  g.coresmp=_ConvertSMVer2Cores(deviceProp.major,deviceProp.minor);
  g.cores=g.coresmp*g.mp;
  g.clockrate=deviceProp.clockRate;
  g.constantmem=deviceProp.totalConstMem;
  g.sharedmem=deviceProp.sharedMemPerBlock;
  g.regsblock=deviceProp.regsPerBlock;
  g.maxthmp=deviceProp.maxThreadsPerMultiProcessor;
  g.maxthblock=deviceProp.maxThreadsPerBlock;
  g.overlap=deviceProp.deviceOverlap;
  g.overlapcount=deviceProp.asyncEngineCount;
  g.limitrun=deviceProp.kernelExecTimeoutEnabled;
  g.integrated=deviceProp.integrated;
  g.maphostmem=deviceProp.canMapHostMemory;
  g.eccmode=deviceProp.ECCEnabled;
  g.uva=deviceProp.unifiedAddressing;
  g.pcidomain=deviceProp.pciDomainID;
  g.pcibus=deviceProp.pciBusID;
  g.pcidevice=deviceProp.pciDeviceID;
  return(g);
}

//==============================================================================
/// Returns information about detected GPUs (code from deviceQuery example).
//==============================================================================
int GetHipDevicesInfo(std::vector<std::string> *gpuinfo,std::vector<StGpuInfo> *gpuprops){
  if(gpuinfo)gpuinfo->push_back("[HIP Capable device(s)]");
  int deviceCount=0;
  hipGetDeviceCount(&deviceCount);
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "Failed getting devices info.");
  if(gpuinfo){
    if(!deviceCount)gpuinfo->push_back("  There are no available device(s) that support HIP");
    else gpuinfo->push_back(fun::PrintStr("  Detected %d HIP Capable device(s)",deviceCount));
  }
  int gid0=-10; hipGetDevice(&gid0);
  //-Driver information.
  int driverVersion=0,runtimeVersion=0;
  hipDriverGetVersion(&driverVersion);
  hipRuntimeGetVersion(&runtimeVersion);
  if(gpuinfo)gpuinfo->push_back(fun::PrintStr("  HIP Driver Version / Runtime Version: %d.%d / %d.%d",driverVersion/1000,(driverVersion%100)/10,runtimeVersion/1000,(runtimeVersion%100)/10));
  //-Devices information.
  for(int dev=0;dev<deviceCount;++dev){
    const StGpuInfo g=GetHipDeviceInfo(dev);
    if(gpuinfo){
      gpuinfo->push_back(" ");
      gpuinfo->push_back(fun::PrintStr("Device %d: \"%s\"",dev,g.name.c_str()));
      gpuinfo->push_back(fun::PrintStr("  HIP Capability Major....: %d.%d",g.ccmajor,g.ccminor));
      gpuinfo->push_back(fun::PrintStr("  Global memory............: %.0f MBytes",(float)g.globalmem/1048576.0f));
      gpuinfo->push_back(fun::PrintStr("  GPU Max Clock rate.......: %.0f MHz (%0.2f GHz)",1e-3f*g.clockrate,1e-6f*g.clockrate));
      gpuinfo->push_back(fun::PrintStr("  Constant memory..........: %.0f KBytes",(float)g.constantmem/1024.f));
      gpuinfo->push_back(fun::PrintStr("  Shared memory per block..: %.0f KBytes",(float)g.sharedmem/1024.f));
      gpuinfo->push_back(fun::PrintStr("  Registers per block......: %d",g.regsblock));
      gpuinfo->push_back(fun::PrintStr("  Maximum threads per MP...: %d",g.maxthmp));
      gpuinfo->push_back(fun::PrintStr("  Maximum threads per block: %d",g.maxthblock));
      gpuinfo->push_back(fun::PrintStr("  Device supports Unified Addressing (UVA): %s",(g.uva? "Yes": "No")));
      gpuinfo->push_back(fun::PrintStr("  Device PCI (Domain / Bus / location)....: %d / %d / %d",g.pcidomain,g.pcibus,g.pcidevice));
    }
    if(gpuprops)gpuprops->push_back(g);
  }
  int gid1=-10; hipGetDevice(&gid1);
  if(gid0>=0 && gid0!=gid1)hipSetDevice(gid0);
  return(deviceCount);
}

//==============================================================================
/// Returns cores per multiprocessor (code from helper_hip.h).
//==============================================================================
int _ConvertSMVer2Cores(int major, int minor){
  typedef struct {
    int SM; // 0xMm (hexadecimal notation), M = SM Major version, and m = SM minor version
    int Cores;
  } sSMtoCores;

  sSMtoCores nGpuArchCoresPerSM[] = {
    { 0x20, 32 }, // Fermi Generation (SM 2.0) - Unsupported by HIP
    { 0x30, 192}, // Kepler Generation (SM 3.0)
    { 0x50, 128}, // Maxwell Generation (SM 5.0)
    { 0x60, 64 }, // Pascal Generation (SM 6.0)
    { 0x70, 64 }, // Volta Generation (SM 7.0)
    { -1, -1 }
  };

  int index = 0;
  while (nGpuArchCoresPerSM[index].SM != -1){
    if(nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) return(nGpuArchCoresPerSM[index].Cores);
    index++;
  }
  return(0);
}

//##############################################################################
//## Functions to allocate GPU memory.
//##############################################################################

//==============================================================================
/// Allocates memory for word on GPU.
//==============================================================================
size_t Malloc(byte **ptr,unsigned count){
  size_t size=sizeof(byte)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for word on GPU.
//==============================================================================
size_t Malloc(word **ptr,unsigned count){
  size_t size=sizeof(word)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for unsigned on GPU.
//==============================================================================
size_t Malloc(unsigned **ptr,unsigned count){
  size_t size=sizeof(unsigned)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for uint4 on GPU.
//==============================================================================
size_t Malloc(uint4 **ptr,unsigned count){
  size_t size=sizeof(uint4)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for int on GPU.
//==============================================================================
size_t Malloc(int **ptr,unsigned count){
  size_t size=sizeof(int)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for int2 on GPU.
//==============================================================================
size_t Malloc(int2 **ptr,unsigned count){
  size_t size=sizeof(int2)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for float on GPU.
//==============================================================================
size_t Malloc(float **ptr,unsigned count){
  size_t size=sizeof(float)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for float2 on GPU.
//==============================================================================
size_t Malloc(float2 **ptr,unsigned count){
  size_t size=sizeof(float2)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for float3 on GPU.
//==============================================================================
size_t Malloc(float3 **ptr,unsigned count){
  size_t size=sizeof(float3)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for float4 on GPU.
//==============================================================================
size_t Malloc(float4 **ptr,unsigned count){
  size_t size=sizeof(float4)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for double on GPU.
//==============================================================================
size_t Malloc(double **ptr,unsigned count){
  size_t size=sizeof(double)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for double2 on GPU.
//==============================================================================
size_t Malloc(double2 **ptr,unsigned count){
  size_t size=sizeof(double2)*count;  hipMalloc((void**)ptr,size);  return(size);
}

//==============================================================================
/// Allocates memory for double3 on GPU.
//==============================================================================
size_t Malloc(double3 **ptr,unsigned count){
  size_t size=sizeof(double3)*count;  hipMalloc((void**)ptr,size);  return(size);
}


//##############################################################################
//## Functions to allocate pinned CPU memory.
//##############################################################################

//==============================================================================
/// Allocates pinned memory for byte on CPU.
//==============================================================================
size_t HostAlloc(byte **ptr,unsigned count){
  size_t size=sizeof(byte)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//==============================================================================
/// Allocates pinned memory for word on CPU.
//==============================================================================
size_t HostAlloc(word **ptr,unsigned count){
  size_t size=sizeof(word)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//==============================================================================
/// Allocates pinned memory for unsigned on CPU.
//==============================================================================
size_t HostAlloc(unsigned **ptr,unsigned count){
  size_t size=sizeof(unsigned)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//==============================================================================
/// Allocates pinned memory for int on CPU.
//==============================================================================
size_t HostAlloc(int **ptr,unsigned count){
  size_t size=sizeof(int)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//==============================================================================
/// Allocates pinned memory for int2 on CPU.
//==============================================================================
size_t HostAlloc(int2 **ptr,unsigned count){
  size_t size=sizeof(int2)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//==============================================================================
/// Allocates pinned memory for float on CPU.
//==============================================================================
size_t HostAlloc(float **ptr,unsigned count){
  size_t size=sizeof(float)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//==============================================================================
/// Allocates pinned memory for tfloat4 on CPU.
//==============================================================================
size_t HostAlloc(tfloat4 **ptr,unsigned count){
  size_t size=sizeof(tfloat4)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//==============================================================================
/// Allocates pinned memory for double on CPU.
//==============================================================================
size_t HostAlloc(double **ptr,unsigned count){
  size_t size=sizeof(double)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//==============================================================================
/// Allocates pinned memory for tdouble2 on CPU.
//==============================================================================
size_t HostAlloc(tdouble2 **ptr,unsigned count){
  size_t size=sizeof(tdouble2)*count;  hipHostMalloc((void**)ptr,size,hipHostMallocDefault);  return(size);
}

//##############################################################################
//## Functions to copy data to Host (debug).
//##############################################################################

//==============================================================================
/// Returns dynamic pointer with byte data. (this pointer must be deleted)
//==============================================================================
byte* ToHostByte(unsigned pini,unsigned n,const byte *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    byte *v=new byte[n];
    hipMemcpy(v,vg+pini,sizeof(byte)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with word data. (this pointer must be deleted)
//==============================================================================
word* ToHostWord(unsigned pini,unsigned n,const word *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    word *v=new word[n];
    hipMemcpy(v,vg+pini,sizeof(word)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with ushort4 data. (this pointer must be deleted)
//==============================================================================
ushort4* ToHostWord4(unsigned pini,unsigned n,const ushort4 *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    ushort4 *v=new ushort4[n];
    hipMemcpy(v,vg+pini,sizeof(ushort4)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with int data. (this pointer must be deleted)
//==============================================================================
int* ToHostInt(unsigned pini,unsigned n,const int *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    int *v=new int[n];
    hipMemcpy(v,vg+pini,sizeof(int)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with unsigned data. (this pointer must be deleted)
//==============================================================================
unsigned* ToHostUint(unsigned pini,unsigned n,const unsigned *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    unsigned *v=new unsigned[n];
    hipMemcpy(v,vg+pini,sizeof(unsigned)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with tint2 data. (this pointer must be deleted)
//==============================================================================
tint2* ToHostInt2(unsigned pini,unsigned n,const int2 *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    tint2 *v=new tint2[n];
    hipMemcpy(v,vg+pini,sizeof(tint2)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with float data. (this pointer must be deleted)
//==============================================================================
float* ToHostFloat(unsigned pini,unsigned n,const float *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    float *v=new float[n];
    hipMemcpy(v,vg+pini,sizeof(float)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with tfloat3 data. (this pointer must be deleted)
//==============================================================================
tfloat3* ToHostFloat3(unsigned pini,unsigned n,const float3 *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    tfloat3 *v=new tfloat3[n];
    hipMemcpy(v,vg+pini,sizeof(tfloat3)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with tfloat4 data. (this pointer must be deleted)
//==============================================================================
tfloat4* ToHostFloat4(unsigned pini,unsigned n,const float4 *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    tfloat4 *v=new tfloat4[n];
    hipMemcpy(v,vg+pini,sizeof(tfloat4)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with double data. (this pointer must be deleted)
//==============================================================================
double* ToHostDouble(unsigned pini,unsigned n,const double *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    double *v=new double[n];
    hipMemcpy(v,vg+pini,sizeof(double)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with tdouble2 data. (this pointer must be deleted)
//==============================================================================
tdouble2* ToHostDouble2(unsigned pini,unsigned n,const double2 *vg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    tdouble2 *v=new tdouble2[n];
    hipMemcpy(v,vg+pini,sizeof(tdouble2)*n,hipMemcpyDeviceToHost);
    CheckHipErrorFun(__FILE__, __LINE__, __func__, "After hipMemcpy().");
    return(v);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with position in tfloat3. (this pointer must be deleted)
//==============================================================================
tfloat3* ToHostPosf3(unsigned pini,unsigned n,const double2 *posxyg,const double *poszg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    tdouble2 *posxy=ToHostDouble2(pini,n,posxyg);
    double   *posz= ToHostDouble(pini,n,poszg);
    tfloat3 *posf=new tfloat3[n];
    for(unsigned p=0;p<n;p++)posf[p]=ToTFloat3(TDouble3(posxy[p].x,posxy[p].y,posz[p]));
    delete[] posxy;  posxy=NULL;
    delete[] posz;   posz =NULL;
    return(posf);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointer with position in tfloat3. (this pointer must be deleted)
//==============================================================================
tdouble3* ToHostPosd3(unsigned pini,unsigned n,const double2 *posxyg,const double *poszg){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    tdouble2 *posxy=ToHostDouble2(pini,n,posxyg);
    double   *posz= ToHostDouble(pini,n,poszg);
    tdouble3 *posd=new tdouble3[n];
    for(unsigned p=0;p<n;p++)posd[p]=TDouble3(posxy[p].x,posxy[p].y,posz[p]);
    delete[] posxy;  posxy=NULL;
    delete[] posz;   posz =NULL;
    return(posd);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

//==============================================================================
/// Returns dynamic pointers with x,y,z data and w values are saved in ptr_w. (theses pointers must be deleted)
//==============================================================================
tfloat3* ToHostFloatXYZ_W(unsigned pini,unsigned n,const float4 *ptrg,float **ptr_w){
  CheckHipErrorFun(__FILE__, __LINE__, __func__, "At the beginning.");
  try{
    const tfloat4 *v4=ToHostFloat4(pini,n,ptrg);
    tfloat3 *vxyz=new tfloat3[n];
    float   *vw  =new float[n];
    for(unsigned p=0;p<n;p++){
      const tfloat4 v=v4[p];
      vxyz[p]=TFloat3(v.x,v.y,v.z);
      vw  [p]=v.w;
    }
    delete[] v4; v4=NULL;
    if(ptr_w)*ptr_w=vw;
    else{ delete[] vw; vw=NULL; }
    return(vxyz);
  }
  catch(const std::bad_alloc){
    fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory (np=%u).",n));
  }
  return(NULL);
}

}








