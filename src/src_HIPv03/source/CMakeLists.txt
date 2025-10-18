#==================================================================
# DualSPHysics GPU/CPU v5.0.140 22-07-2020
#==================================================================

cmake_minimum_required(VERSION 2.8)

PROJECT(DualSPHysics)

#------------------------------------------------------------------
# Coupling options
#------------------------------------------------------------------
option(ENABLE_MOORDYN "Enable the MoorDyn+ library" ON)
option(ENABLE_CHRONO "Enable the Chrono Engine library" ON)
option(ENABLE_CHRONO_OMP "Enable the Chrono Parallel module" ON)

#------------------------------------------------------------------
# Source files
#------------------------------------------------------------------

# CPU Objects
set(OBJXML JXml.cpp tinystr.cpp tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp)
set(OBJSPHMOTION JMotion.cpp JMotionList.cpp JMotionMov.cpp JMotionObj.cpp JMotionPos.cpp JDsMotion.cpp)
set(OBCOMMON Functions.cpp FunctionsGeo3d.cpp FunSphKernelsCfg.cpp JAppInfo.cpp JBinaryData.cpp JCfgRunBase.cpp JDataArrays.cpp JException.cpp JLinearValue.cpp JLog2.cpp JMeanValues.cpp JObject.cpp JOutputCsv.cpp JRadixSort.cpp JRangeFilter.cpp JReadDatafile.cpp JSaveCsv2.cpp JTimeControl.cpp randomc.cpp)
set(OBCOMMONDSPH JDsphConfig.cpp JDsPips.cpp JPartDataBi4.cpp JPartDataHead.cpp JPartFloatBi4.cpp JPartOutBi4Save.cpp JCaseCtes.cpp JCaseEParms.cpp JCaseParts.cpp JCaseProperties.cpp JCaseUserVars.cpp JCaseVtkOut.cpp)
set(OBSPH JArraysCpu.cpp JCellDivCpu.cpp JSphCfgRun.cpp JDsDamping.cpp JDsGaugeItem.cpp JDsGaugeSystem.cpp JDsPartsOut.cpp JDsSaveDt.cpp JSphShifting.cpp JSph.cpp JDsAccInput.cpp JSphCpu.cpp JDsInitialize.cpp JSphMk.cpp JDsPartsInit.cpp JDsFixedDt.cpp JDsViscoInput.cpp JDsOutputTime.cpp JWaveAwasZsurf.cpp JWaveSpectrumGpu.cpp main.cpp)
set(OBSPHSINGLE JCellDivCpuSingle.cpp JPartsLoad4.cpp JSphCpuSingle.cpp)

# GPU Objects for ROCm/HIP
set(OBCOMMONGPU FunctionsHip.cpp JObjectGpu.cpp)
set(OBSPHGPU JArraysGpu.cpp JDebugSphGpu.cpp JCellDivGpu.cpp JSphGpu.cpp)
set(OBSPHSINGLEGPU JCellDivGpuSingle.cpp JSphGpuSingle.cpp)
set(OBHIP JCellDivGpu_ker.hip JCellDivGpuSingle_ker.hip JDsPips_ker.hip JDsGauge_ker.hip JReduSum_ker.hip JSphShifting_ker.hip JDsAccInput_ker.hip JSphGpu_ker.hip JSphGpuSimple_ker.hip JWaveOrder2_ker.hip)

# Other Objects
set(OBWAVERZ JMLPistonsGpu.cpp JRelaxZonesGpu.cpp)
set(OBCHRONO JChronoObjects.cpp)
set(OBMOORDYN JDsMooredFloatings.cpp JDsFtForcePoints.cpp)
set(OBINOUT JSphCpu_InOut.cpp JSphCpuSingle_InOut.cpp JSphBoundCorr.cpp JSphInOut.cpp JSphInOutZone.cpp JSphInOutGridData.cpp JSphInOutPoints.cpp JSimpleNeigs.cpp)
set(OBINOUTGPU JSphGpuSingle_InOut.cpp)
set(OBMDBC JPartNormalData.cpp JNormalsMarrone.cpp)

#------------------------------------------------------------------
# HIP Flags
#------------------------------------------------------------------
find_package(HIP REQUIRED)

if(HIP_FOUND)
  message("HIP libraries were found.")
  list(APPEND HIP_HIPCC_FLAGS "--amdgpu-target=gfx90a -ffast-math -O3")
else()
  message(FATAL_ERROR "HIP Libraries were not found.")
endif()

#------------------------------------------------------------------
# Static libraries linker path
#------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(PLATFORM "linux")
  message("Using libraries for gcc")
  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/linux_gcc)
elseif(MSVC)
  message("Windows version")
  set(PLATFORM "windows")
  #if(MSVC_VERSION VERSION_EQUAL 1900)
    LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/vs2015)
  #endif()
endif()

#--------------------------------------------
# OpenMP flags
#---------------------------------------------
find_package(OpenMP)

if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#------------------------------------------------------------------
# Don't use C++11 ABI for compatibility between GCC v4.9 and >=5
#------------------------------------------------------------------
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC")

#------------------------------------------------------------------
# Binaries
#------------------------------------------------------------------
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../bin/${PLATFORM}) # Generic output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})         # Release output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})           # Debug output directory

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_executable(DualSPHysics5.0CPU_linux64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBWAVERZ} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBMDBC})
  install(TARGETS	DualSPHysics5.0CPU_linux64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  if (HIP_FOUND)
    hip_add_executable(DualSPHysics5.0_linux64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBCOMMONGPU} ${OBSPHGPU} ${OBSPHSINGLEGPU} ${OBHIP} ${OBWAVERZ} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBINOUTGPU} ${OBMDBC})
    install(TARGETS DualSPHysics5.0_linux64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endif(HIP_FOUND)
elseif(MSVC) 
  add_executable(DualSPHysics5.0CPU_win64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBWAVERZ} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBMDBC})
  install(TARGETS	DualSPHysics5.0CPU_win64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  if (HIP_FOUND)
    hip_add_executable(DualSPHysics5.0_win64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBCOMMONGPU} ${OBSPHGPU} ${OBSPHSINGLEGPU} ${OBHIP} ${OBWAVERZ} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBINOUTGPU} ${OBMDBC})
    install(TARGETS DualSPHysics5.0_win64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endif(HIP_FOUND)
endif()

foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
  string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
endforeach(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

#------------------------------------------------------------------
# General linker flags
#------------------------------------------------------------------
set(LINKER_FLAGS "")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(LINKER_FLAGS ${LINKER_FLAGS} jvtklib_64 jwavegen_64 jnumexlib_64)
elseif(MSVC) 
  set(LINKER_FLAGS ${LINKER_FLAGS} LibJVtkLib_x64_v140_Release LibJWaveGen_x64_v140_Release LibJNumexLib_x64_v140_Release)
endif()

#------------------------------------------------------------------
# Coupled libraries
#------------------------------------------------------------------
# MoorDyn+
if(ENABLE_MOORDYN)
  message(STATUS "Using MoorDyn+")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(LINKER_FLAGS ${LINKER_FLAGS} dsphmoordyn_64)
  elseif(MSVC) 
    set(LINKER_FLAGS ${LINKER_FLAGS} LibDSphMoorDyn_x64_v140_Release)
  endif()
else()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_MOORDYN")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_MOORDYN")
endif(ENABLE_MOORDYN)

# Chrono Engine
if(ENABLE_CHRONO)
  message(STATUS "Using ChronoEngine Module")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(LINKER_FLAGS ${LINKER_FLAGS} ChronoEngine dsphchrono)
  elseif(MSVC) 
    set(LINKER_FLAGS ${LINKER_FLAGS} dsphchrono)
  endif()
else()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_CHRONO")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_CHRONO")
endif(ENABLE_CHRONO)

# Chrono Parallel
if(ENABLE_CHRONO_OMP)
  if(ENABLE_CHRONO)
    message(STATUS "Using ChronoEngine Parallel Module")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(LINKER_FLAGS ${LINKER_FLAGS} ChronoEngine_parallel)
    endif()
  else()
    message(FATAL_ERROR "You should enable CHRONO option to use CHRONO_OMP")
  endif(ENABLE_CHRONO)
else()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_CHRONO_OMP")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_CHRONO_OMP")
endif(ENABLE_CHRONO_OMP)

#------------------------------------------------------------------
# Linker flags
#------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(DualSPHysics5.0CPU_linux64 ${LINKER_FLAGS})
  set_target_properties(DualSPHysics5.0CPU_linux64 PROPERTIES COMPILE_FLAGS "-ffast-math -O3 -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC -std=c++0x")
  
  if (HIP_FOUND)
    target_link_libraries(DualSPHysics5.0_linux64 ${LINKER_FLAGS})
    set_target_properties(DualSPHysics5.0_linux64 PROPERTIES COMPILE_FLAGS "-ffast-math -O3 -D_WITHGPU -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC -std=c++0x")
  endif(HIP_FOUND)
  
elseif(MSVC)
  set_target_properties(DualSPHysics5.0_win64 PROPERTIES COMPILE_FLAGS "/D _WITHGPU")

  #if(MSVC_VERSION VERSION_EQUAL 1900)
    # MSVC 2015
    target_link_libraries(DualSPHysics5.0CPU_win64  ${LINKER_FLAGS})
    target_link_libraries(DualSPHysics5.0_win64 ${LINKER_FLAGS})
  #endif()
  
  SET(HIP_PROPAGATE_HOST_FLAGS OFF CACHE BOOL "Propagate C/CXX Flags and friends to the host compiler in HIP via -Xompile  " FORCE)
  
  foreach(CPP_FLAGS CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(${CPP_FLAGS} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${CPP_FLAGS} "${${CPP_FLAGS}}")
    endif(${CPP_FLAGS} MATCHES "/MD")
  endforeach(CPP_FLAGS)
endif()
