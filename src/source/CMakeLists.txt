#==================================================================
# DualSPHysics GPU/CPU v5.0.140 22-07-2020 Adaptado para ROCm/HIP
#==================================================================

cmake_minimum_required(VERSION 3.10)

project(DualSPHysics)

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

# GPU Objects (Adaptados para HIP)
set(OBCOMMONGPU FunctionsHip.cpp JObjectGpu.cpp)
set(OBSPHGPU JArraysGpu.cpp JDebugSphGpu.cpp JCellDivGpu.cpp JSphGpu.cpp JSphGpuSingle_InOut.cpp)
set(OBSPHSINGLEGPU JCellDivGpuSingle.cpp JSphGpuSingle.cpp)
set(OBHIP JCellDivGpu_ker.hip JCellDivGpuSingle_ker.hip JDsPips_ker.hip JDsGauge_ker.hip JReduSum_ker.hip JSphShifting_ker.hip JDsAccInput_ker.hip JSphGpu_ker.hip JSphGpuSimple_ker.hip JWaveOrder2_ker.hip JSphGpu_InOut_iker.hip)

# Otros objetos relacionados con la relajación en GPU
set(OBRELAXZONEGPU JRelaxZonesGpu.cpp JRelaxZone_ker.hip)

# Otros objetos relacionados con la división de celdas en GPU
set(OBCELLDIVGPU JCellDivGpu.cpp JCellDivGpuSingle.cpp)

# Otros Objetos
set(OBWAVERZ JMLPistonsGpu.cpp JRelaxZonesGpu.cpp)
set(OBWAVERZHIP JRelaxZone_ker.hip)
set(OBCHRONO JChronoObjects.cpp)
set(OBMOORDYN JDsMooredFloatings.cpp JDsFtForcePoints.cpp)
set(OBINOUT JSphCpu_InOut.cpp JSphCpuSingle_InOut.cpp JSphBoundCorr.cpp JSphInOut.cpp JSphInOutZone.cpp JSphInOutGridData.cpp JSphInOutPoints.cpp JSimpleNeigs.cpp)
set(OBINOUTGPU JSphGpuSingle_InOut.cpp)
set(OBMDBC JPartNormalData.cpp JNormalsMarrone.cpp)

# Inclusión de namespaces específicos para HIP
set(OBNAMESPACEINOUT JSphGpu_InOut_iker.hip)  # Archivos de 'cusphinout'
set(OBNAMESPACERELAXZONE JRelaxZone_ker.hip)  # Archivos de 'curelaxzone'
set(OBNAMESPACEHIPDIV JCellDivGpu_ker.hip JCellDivGpuSingle_ker.hip)  # Archivos de 'hipdiv'

# Combinación de todos los objetos GPU en un solo conjunto
set(ALLGPU ${OBCOMMONGPU} ${OBSPHGPU} ${OBSPHSINGLEGPU} ${OBHIP} ${OBRELAXZONEGPU} ${OBCELLDIVGPU} ${OBWAVERZHIP} ${OBINOUTGPU} ${OBNAMESPACEINOUT} ${OBNAMESPACERELAXZONE} ${OBNAMESPACEHIPDIV})


#------------------------------------------------------------------
# Set Compilation Flags
#------------------------------------------------------------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14 -fPIC -D_GLIBCXX_USE_CXX11_ABI=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

#------------------------------------------------------------------
# Find HIP and ROCm Libraries
#------------------------------------------------------------------
find_package(HIP REQUIRED)

if(HIP_FOUND)
    message(STATUS "Using HIP for GPU acceleration")
    include_directories(${HIP_INCLUDE_DIRS})
    link_directories(${HIP_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "HIP Libraries were not found.")
endif()

#------------------------------------------------------------------
# OpenMP Configuration
#------------------------------------------------------------------
find_package(OpenMP)

if(OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  message(STATUS "Found OpenMP: ${OpenMP_CXX_FLAGS}")
else()
  message(WARNING "OpenMP not found or incomplete")
endif()

#------------------------------------------------------------------
# Define Output Directory
#------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../bin/linux)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

#------------------------------------------------------------------
# Compilar ejecutable CPU
#------------------------------------------------------------------
add_executable(DualSPHysics5.0CPU_linux64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBWAVERZ} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBMDBC})

# Enlace para CPU
target_link_libraries(DualSPHysics5.0CPU_linux64 
    /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libdsphmoordyn_64.a 
    /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libChronoEngine.so 
    /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libdsphchrono.so 
    /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libChronoEngine_parallel.so 
    /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libjwavegen_64.a
    /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libjvtklib_64.a
    /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libjnumexlib_64.a
    ${OpenMP_CXX_LIBRARIES} -pthread)

set_target_properties(DualSPHysics5.0CPU_linux64 PROPERTIES COMPILE_FLAGS "-O3 -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC -std=c++14")
install(TARGETS DualSPHysics5.0CPU_linux64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

#------------------------------------------------------------------
# Compilar ejecutable GPU (HIP)
#------------------------------------------------------------------
if (HIP_FOUND)
  add_executable(DualSPHysics5.0_linux64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBCOMMONGPU} ${OBSPHGPU} ${OBSPHSINGLEGPU} ${OBHIP} ${OBRELAXZONEGPU} ${OBCELLDIVGPU} ${OBWAVERZ} ${OBWAVERZHIP} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBINOUTGPU} ${OBMDBC})


  # Enlace para GPU
target_link_libraries(DualSPHysics5.0_linux64 
      /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libdsphmoordyn_64.a 
      /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libChronoEngine.so 
      /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libdsphchrono.so 
      /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libChronoEngine_parallel.so 
      /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libjwavegen_64.a
      /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libjvtklib_64.a
      /home/fasanchezc/dualsphysics_hip/dualsphysics_hip/src/lib/linux_gcc/libjnumexlib_64.a
      ${HIP_LIBRARIES} ${OpenMP_CXX_LIBRARIES} -pthread)

  set_target_properties(DualSPHysics5.0_linux64 PROPERTIES COMPILE_FLAGS "-O3 -D_WITHGPU -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC -std=c++14")
  install(TARGETS DualSPHysics5.0_linux64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

#------------------------------------------------------------------
# Coupled Libraries
#------------------------------------------------------------------
if(ENABLE_MOORDYN)
  message(STATUS "Using MoorDyn+")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_MOORDYN")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_MOORDYN")
endif()

if(ENABLE_CHRONO)
  message(STATUS "Using ChronoEngine Module")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_CHRONO")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_CHRONO")
endif()

if(ENABLE_CHRONO_OMP)
  if(ENABLE_CHRONO)
    message(STATUS "Using ChronoEngine Parallel Module")
  else()
    message(FATAL_ERROR "You should enable CHRONO option to use CHRONO_OMP")
  endif()
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_CHRONO_OMP")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_CHRONO_OMP")
endif()







