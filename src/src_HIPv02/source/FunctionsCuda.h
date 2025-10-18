﻿//HEAD_DSCODES
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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Implementa funciones basicas de HIP. (17-02-2017)
//:# - Reconoce arquitectura Pascal hasta 6.2. (14-03-2017)
//:# - Ahora GetHipDevicesInfo() no cambia la seleccion de GPU. (01-07-2017)
//:# - Documentacion del codigo en ingles. (08-08-2017)
//:# - Nueva funcion ToHostInt(). (30-01-2018)
//:# - Incluye cores para SM 75. (20-12-2018)
//:# - Nueva funcion ToHostFloatXYZ_W(). (10-09-2019)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:# - Nuevas funciones HostAlloc(byte), ToHostByte(). (18-01-2020)
//:#############################################################################

/// \file FunctionsHip.h \brief Declares basic/general GPU functions for the entire application.

#ifndef _FunctionsHip_
#define _FunctionsHip_

#include <string>
#include <vector>
#include <hip/hip_runtime.h>  // Reemplazamos CUDA por HIP
#include "TypesDef.h"

#ifndef Check_HipErrorFun
#define Check_HipErrorFun(msg) CheckHipErrorFun(__FILE__,__LINE__,__func__,msg)
#endif

/// Implements a set of basic/general GPU functions.
namespace fhip{

/// Structure to save basic GPU information.
typedef struct StrGpuInfo{
  int id;              ///<Device number.
  std::string name;    ///<GPU name.
  int ccmajor;         ///<Capability major.
  int ccminor;         ///<Capability minor.
  ullong globalmem;    ///<Global memory.
  int cores;           ///<GPU Cores.
  int mp;              ///<Multiprocessors.
  int coresmp;         ///<GPU Cores/MP.
  int clockrate;       ///<GPU Max Clock rate.
  int clockratemem;    ///<Memory Clock rate.
  int busmem;          ///<Memory Bus Width.
  int cachelv2;        ///<L2 Cache Size.
  ullong constantmem;  ///<Constant memory.
  ullong sharedmem;    ///<Shared memory per block.
  int regsblock;       ///<Registers per block.
  int maxthmp;         ///<Maximum threads per MP.
  int maxthblock;      ///<Maximum threads per block.
  int overlap;         ///<Concurrent copy and kernel execution.
  int overlapcount;    ///<Number of concurrent copy and kernel execution.
  int limitrun;        ///<Run time limit on kernels.
  int integrated;      ///<Integrated GPU sharing Host Memory.
  int maphostmem;      ///<Support host page-locked memory mapping.
  int eccmode;         ///<Device has ECC support.
  int tccdriver;       ///<GPU Device Driver Mode.
  int uva;             ///<Device supports Unified Addressing (UVA).
  int pcidomain;       ///<PCI Domain ID.
  int pcibus;          ///<PCI Bus ID.
  int pcidevice;       ///<PCI location ID.
  bool rdma;           ///<Device supports P2P and RDMA.
  static const int sizep2pto=8; ///<Maximum number of GPUs with peer access to/from.
  int countp2pto;               ///<Number of GPUs with peer access to/from.
  int p2pto[sizep2pto];         ///<GPU list with peer access to/from.

  StrGpuInfo(){ 
    id=-1;
    name=""; 
    ccmajor=ccminor=0;
    globalmem=0;
    cores=mp=coresmp=clockrate=0;
    clockratemem=busmem=cachelv2=0;
    constantmem=sharedmem=0;
    regsblock=maxthmp=maxthblock=0;
    overlap=overlapcount=0;
    limitrun=integrated=maphostmem=eccmode=0;
    tccdriver=uva=pcidomain=pcibus=pcidevice=0;
    rdma=false;
    countp2pto=0;
    for(int c=0;c<sizep2pto;c++)p2pto[c]=-1;
  }
}StGpuInfo;

void CheckHipErrorFun(const char *const file,int const line
  ,const char *const fun,std::string msg);

std::string GetHipDeviceName(int gid);
StGpuInfo GetHipDeviceInfo(int gid);
int GetHipDevicesInfo(std::vector<std::string> *gpuinfo,std::vector<StGpuInfo> *gpuprops);
int _ConvertSMVer2Cores(int major, int minor);

//-Functions to allocate GPU memory.
size_t Malloc(byte     **,unsigned count);
size_t Malloc(word     **,unsigned count);
size_t Malloc(unsigned **,unsigned count);
size_t Malloc(uint4    **,unsigned count);
size_t Malloc(int      **,unsigned count);
size_t Malloc(int2     **,unsigned count);
size_t Malloc(float    **,unsigned count);
size_t Malloc(float2   **,unsigned count);
size_t Malloc(float3   **,unsigned count);
size_t Malloc(float4   **,unsigned count);
size_t Malloc(double   **,unsigned count);
size_t Malloc(double2  **,unsigned count);
size_t Malloc(double3  **,unsigned count);
//:hipFree(GpuMem);

//-Functions to allocate pinned CPU memory.
size_t HostAlloc(byte     **,unsigned count);
size_t HostAlloc(word     **,unsigned count);
size_t HostAlloc(unsigned **,unsigned count);
size_t HostAlloc(int      **,unsigned count);
size_t HostAlloc(int2     **,unsigned count);
size_t HostAlloc(float    **,unsigned count);
size_t HostAlloc(tfloat4  **,unsigned count);
size_t HostAlloc(double   **,unsigned count);
size_t HostAlloc(tdouble2 **,unsigned count);
//:hipFreeHost(PinnedMem);

//-Functions to copy data to Host (debug).
byte*     ToHostByte   (unsigned pini,unsigned n,const byte     *ptrg);
word*     ToHostWord   (unsigned pini,unsigned n,const word     *ptrg);
ushort4*  ToHostWord4  (unsigned pini,unsigned n,const ushort4  *ptrg);
int*      ToHostInt    (unsigned pini,unsigned n,const int      *ptrg);
unsigned* ToHostUint   (unsigned pini,unsigned n,const unsigned *ptrg);
tint2*    ToHostInt2   (unsigned pini,unsigned n,const int2     *ptrg);
float*    ToHostFloat  (unsigned pini,unsigned n,const float    *ptrg);
tfloat3*  ToHostFloat3 (unsigned pini,unsigned n,const float3   *ptrg);
tfloat4*  ToHostFloat4 (unsigned pini,unsigned n,const float4   *ptrg);
double*   ToHostDouble (unsigned pini,unsigned n,const double   *ptrg);
tdouble2* ToHostDouble2(unsigned pini,unsigned n,const double2  *ptrg);

tfloat3*  ToHostPosf3(unsigned pini,unsigned n,const double2 *posxyg,const double *poszg);
tdouble3* ToHostPosd3(unsigned pini,unsigned n,const double2 *posxyg,const double *poszg);
tfloat3*  ToHostFloatXYZ_W(unsigned pini,unsigned n,const float4 *ptrg,float **ptr_w);

}

//:-SYNCHRONIZATION:
//: hipDeviceSynchronize();           //-Blocks host until all issued HIP calls are complete.
//: hipStreamSynchronize(stream);     //-Blocks host until all issued HIP calls in stream are complete.
//: hipEventSynchronize(event);       //-Blocks host until event has occurred. 
//: hipStreamWaitEvent(stream,event,0); //-Blocks stream until event occurs. Only blocks launches after this call. Does not block the host!

//:-MANAGING EVENTS:
//: hipEvent_t event=NULL;
//: hipEventCreate(&event);    //-Creates an event.
//: hipEventDestroy(event);   //-Destroys an event.
//: hipEventCreateWithFlags(&event,hipEventDisableTiming);  //-Disables timing to increase performance and avoid synchronization issues.
//: hipEventRecord(event,stream); //-Set the event state to not occurred. Enqueue the event into a stream. Event state is set to occurred when it reaches the front of the stream.
//: float t=0; hipEventElapsedTime(&t,evt,evt2); //-Devuelve tiempo entre eventos. Entre el hipEventRecord() y hipEventElapsedTime() tiene que haber alguna llamada de sincronizacion (hipDeviceSynchronize(), hipStreamSynchronize() o hipEventSynchronize()).

//:-MANAGING STREAMS:
//: hipStream_t stm=NULL;
//: hipStreamCreate(&stm);
//: if(stm)hipStreamDestroy(stm);  stm=NULL;

//:-OVERLAPPING KERNELS AND MEMORY OPERATIONS:
//: - El tiempo total se reduce porque hay cierto solapamiento, pero el tiempo de cada operacion aumenta.

//:-OTHER FUNCTIONS:
//: hipMemset(datad,0,sizeof(unsigned)*n);
//: hipMemsetAsync(datad,0,sizeof(unsigned)*n,stm);
//: hipMemcpy(datah,datad,sizeof(unsigned)*n,hipMemcpyDeviceToHost);
//: hipMemcpyAsync(datah,datad,sizeof(unsigned)*n,hipMemcpyDeviceToHost,stm);
#endif
