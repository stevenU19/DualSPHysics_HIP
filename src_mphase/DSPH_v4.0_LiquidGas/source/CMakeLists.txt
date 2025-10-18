cmake_minimum_required(VERSION 3.10)

PROJECT(DualSPHysics)

# Archivos fuente
set(OBJ_BASIC main.cpp Functions.cpp FunctionsMath.cpp JArraysCpu.cpp JBinaryData.cpp JCellDivCpu.cpp JCfgRun.cpp JException.cpp JLog2.cpp JObject.cpp JPartDataBi4.cpp JPartFloatBi4.cpp JPartOutBi4Save.cpp JPartsOut.cpp JRadixSort.cpp JRangeFilter.cpp JReadDatafile.cpp JSaveDt.cpp JSpaceCtes.cpp JSpaceEParms.cpp JSpaceParts.cpp JSpaceProperties.cpp JSph.cpp JSphAccInput.cpp JSphCpu.cpp JSphDtFixed.cpp JSphVisco.cpp randomc.cpp JTimeOut.cpp)
set(OBJ_CPU_SINGLE JCellDivCpuSingle.cpp JSphCpuSingle.cpp JPartsLoad4.cpp)
set(OBJ_GPU JArraysGpu.cpp JCellDivGpu.cpp JObjectGpu.cpp JSphGpu.cpp JBlockSizeAuto.cpp JMeanValues.cpp)
set(OBJ_GPU_SINGLE JCellDivGpuSingle.cpp JSphGpuSingle.cpp)

set(OBJ_HIP JCellDivGpu_ker.hip JSphGpu_ker.hip)
set(OBJ_HIP_SINGLE JCellDivGpuSingle_ker.hip)

# Encontrar y configurar HIP
find_package(HIP REQUIRED)

if(HIP_FOUND)
  message("HIP libraries were found.")
  list(APPEND HIP_HIPCC_FLAGS "--amdgpu-target=gfx908 -ffast-math -O3")
else()
  message(FATAL_ERROR "HIP Libraries were not found.")
endif()

# Configuración para GNU o MSVC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("Using libraries for gcc")
    LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
elseif(MSVC)
  message("Windows version")

  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Libs)

  # Comprobación de versión de MSVC (descomentada si es necesario)
  #if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "your.required.msvc.version")
  #  message(FATAL_ERROR "Insufficient MSVC version")
  #endif()
endif()

# Encontrar OpenMP
find_package(OpenMP)

if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Creación de los ejecutables
add_executable(dualsphysics4cpu ${OBJ_BASIC} ${OBJ_CPU_SINGLE})
hip_add_executable(dualsphysics4gpu ${OBJ_BASIC} ${OBJ_CPU_SINGLE} ${OBJ_GPU} ${OBJ_HIP} ${OBJ_GPU_SINGLE} ${OBJ_HIP_SINGLE})

# Instalación de los ejecutables
install(TARGETS dualsphysics4cpu dualsphysics4gpu DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../EXECS)
	
# Configuración de flags y enlazado
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(dualsphysics4cpu jxml_64 jformatfiles2_64 jsphmotion_64 jwavegen_64)
  target_link_libraries(dualsphysics4gpu jxml_64 jformatfiles2_64 jsphmotion_64 jwavegen_64)
  set_target_properties(dualsphysics4cpu PROPERTIES COMPILE_FLAGS "-ffast-math -O3")	
  set_target_properties(dualsphysics4gpu PROPERTIES COMPILE_FLAGS "-ffast-math -O3 -D_WITHGPU")	
elseif(MSVC)
  set_target_properties(dualsphysics4gpu PROPERTIES COMPILE_FLAGS "/D _WITHGPU")

  # MSVC 10
  if(MSVC_VERSION VERSION_EQUAL 1600)
    target_link_libraries(dualsphysics4cpu JXml_x64_v100_Release JFormatFiles2_x64_v100_Release JSphMotion_x64_v100_Release JWaveGen_x64_v100_Release)
    target_link_libraries(dualsphysics4gpu JXml_x64_v100_Release JFormatFiles2_x64_v100_Release JSphMotion_x64_v100_Release JWaveGen_x64_v100_Release)
  elseif(MSVC_VERSION VERSION_EQUAL 1800) 
    target_link_libraries(dualsphysics4cpu JXml_x64_v120_Release JFormatFiles2_x64_v120_Release JSphMotion_x64_v120_Release JWaveGen_x64_v120_Release)
    target_link_libraries(dualsphysics4gpu JXml_x64_v120_Release JFormatFiles2_x64_v120_Release JSphMotion_x64_v120_Release JWaveGen_x64_v120_Release)
  endif()  

  # Desactivar la propagación de flags de host en NVCC (no necesario para HIP)
  # SET(CUDA_PROPAGATE_HOST_FLAGS OFF CACHE BOOL "Propagate C/CXX Flags and friends to the host compiler in NVCC via -Xompile  " FORCE)

  foreach(CPP_FLAGS CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(${CPP_FLAGS} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${CPP_FLAGS} "${${CPP_FLAGS}}")
    endif(${CPP_FLAGS} MATCHES "/MD")
  endforeach(CPP_FLAGS)
endif()
