# DualSPHysics GPU/CPU v5.0.140 23-07-2020

#=============== Compilation Options (YES/NO) ===============
USE_GCC5=NO
USE_DEBUG=NO
USE_FAST_MATH=YES
USE_NATIVE_CPU_OPTIMIZATIONS=NO
COMPILE_VTKLIB=YES
COMPILE_NUMEXLIB=YES
COMPILE_CHRONO=YES
COMPILE_CHRONO_OMP=YES
COMPILE_WAVEGEN=YES
COMPILE_MOORDYN=YES

LIBS_DIRECTORIES=-L./
LIBS_DIRECTORIES:=$(LIBS_DIRECTORIES) -L../lib/linux_gcc

EXECNAME=DualSPHysics5.0_linux64
EXECS_DIRECTORY=../../bin/linux

# -std=c++0x ---> Used to avoid errors for calls to enums
ifeq ($(USE_DEBUG), YES)
  CCFLAGS=-c -O0 -g -Wall -fopenmp -D_WITHGPU -std=c++0x
else
  CCFLAGS=-c -O3 -fopenmp -D_WITHGPU -std=c++0x
  ifeq ($(USE_FAST_MATH), YES)
    CCFLAGS+= -ffast-math
  endif
  ifeq ($(USE_NATIVE_CPU_OPTIMIZATIONS), YES)
    CCFLAGS+= -march=native
  endif
endif
CC=g++
CCLINKFLAGS=-fopenmp -lgomp

# Required for GCC versions >=5.0
ifeq ($(USE_GCC5), YES)
  CCFLAGS+=-D_GLIBCXX_USE_CXX11_ABI=0
  CCLINKFLAGS+=-D_GLIBCXX_USE_CXX11_ABI=0
  NCCFLAGS+=-D_GLIBCXX_USE_CXX11_ABI=0
endif

ifeq ($(COMPILE_VTKLIB), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_VTKLIB
endif
ifeq ($(COMPILE_NUMEXLIB), NO)
  CCMOREFLAGS:=$(CCMOREFLAGS) -DDISABLE_NUMEXLIB
endif
ifeq ($(COMPILE_CHRONO), NO)
  COMPILE_CHRONO_OMP=NO
  CCFLAGS:=$(CCFLAGS) -DDISABLE_CHRONO
endif
ifeq ($(COMPILE_CHRONO_OMP), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_CHRONO_OMP
endif
ifeq ($(COMPILE_WAVEGEN), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_WAVEGEN
endif
ifeq ($(COMPILE_MOORDYN), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_MOORDYN
endif

#=============== ROCm selection ===============
ROCM_VERSION=5.0

#=============== ROCm toolkit directory ===============
ROCM_PATH=/opt/rocm-$(ROCM_VERSION)

#=============== Select GPU architectures ===============
GENCODE=-amdgpu-target=gfx908

#=============== Files to compile ===============
OBJXML=JXml.o tinystr.o tinyxml.o tinyxmlerror.o tinyxmlparser.o
OBJSPHMOTION=JMotion.o JMotionList.o JMotionMov.o JMotionObj.o JMotionPos.o JDsMotion.o
OBCOMMON=Functions.o FunctionsGeo3d.o FunSphKernelsCfg.o JAppInfo.o JBinaryData.o JCfgRunBase.o JDataArrays.o JException.o JLinearValue.o JLog2.o JMeanValues.o JObject.o JOutputCsv.o JRadixSort.o JRangeFilter.o JReadDatafile.o JSaveCsv2.o JTimeControl.o randomc.o
OBCOMMONDSPH=JDsphConfig.o JDsPips.o JPartDataBi4.o JPartDataHead.o JPartFloatBi4.o JPartOutBi4Save.o JCaseCtes.o JCaseEParms.o JCaseParts.o JCaseProperties.o JCaseUserVars.o JCaseVtkOut.o
OBSPH=JArraysCpu.o JCellDivCpu.o JSphCfgRun.o JDsDamping.o JDsGaugeItem.o JDsGaugeSystem.o JDsPartsOut.o JDsSaveDt.o JSphShifting.o JSph.o JDsAccInput.o JSphCpu.o JDsInitialize.o JSphMk.o JDsPartsInit.o JDsFixedDt.o JDsViscoInput.o JDsOutputTime.o JWaveAwasZsurf.o JWaveSpectrumGpu.o main.o
OBSPHSINGLE=JCellDivCpuSingle.o JPartsLoad4.o JSphCpuSingle.o
OBCOMMONGPU=FunctionsHip.o JObjectGpu.o 
OBSPHGPU=JArraysGpu.o JDebugSphGpu.o JCellDivGpu.o JSphGpu.o 
OBSPHSINGLEGPU=JCellDivGpuSingle.o JSphGpuSingle.o
OBHIP=JCellDivGpu_ker.o JCellDivGpuSingle_ker.o JDsPips_ker.o JDsGauge_ker.o JReduSum_ker.o JSphShifting_ker.o JDsAccInput_ker.o JSphGpu_ker.o JSphGpuSimple_ker.o JWaveOrder2_ker.o

OBWAVERZ=JMLPistonsGpu.o JRelaxZonesGpu.o
OBWAVERZHIP=JRelaxZone_ker.o
OBCHRONO=JChronoObjects.o
OBMOORDYN=JDsMooredFloatings.o JDsFtForcePoints.o
OBINOUT=JSphCpu_InOut.o JSphCpuSingle_InOut.o JSphBoundCorr.o JSphInOut.o JSphInOutZone.o JSphInOutGridData.o JSphInOutPoints.o JSimpleNeigs.o
OBINOUTGPU=JSphGpuSingle_InOut.o
OBMDBC=JPartNormalData.o JNormalsMarrone.o

OBJECTS=$(OBJXML) $(OBJSPHMOTION) $(OBCOMMON) $(OBCOMMONDSPH) $(OBSPH) $(OBSPHSINGLE)
OBJECTS:=$(OBJECTS) $(OBCOMMONGPU) $(OBSPHGPU) $(OBSPHSINGLEGPU) $(OBHIP)
OBJECTS:=$(OBJECTS) $(OBWAVERZ) $(OBWAVERZHIP) $(OBCHRONO) $(OBMOORDYN) $(OBINOUT) $(OBINOUTGPU) $(OBMDBC)

#=============== DualSPHysics libs to be included ===============
JLIBS=${LIBS_DIRECTORIES}
ifeq ($(COMPILE_VTKLIB), YES)
  JLIBS:=$(JLIBS) -ljvtklib_64
endif
ifeq ($(COMPILE_NUMEXLIB), YES)
  JLIBS:=$(JLIBS) -ljnumexlib_64
endif
ifeq ($(COMPILE_CHRONO), YES)
  JLIBS:=$(JLIBS) -ldsphchrono -lChronoEngine 
endif
ifeq ($(COMPILE_CHRONO_OMP), YES)
  JLIBS:=$(JLIBS) -lChronoEngine_parallel
 endif
ifeq ($(COMPILE_WAVEGEN), YES)
  JLIBS:=$(JLIBS) -ljwavegen_64
endif
ifeq ($(COMPILE_MOORDYN), YES)
  JLIBS:=$(JLIBS) -ldsphmoordyn_64
endif

#=============== GPU Code Compilation ===============
CCFLAGS := $(CCFLAGS) -I./ -I$(ROCM_PATH)/include
CCLINKFLAGS := $(CCLINKFLAGS) -L$(ROCM_PATH)/lib -lamdhip64 -ldl -lrt
NCC=hipcc
ifeq ($(USE_DEBUG), NO)
  NCCFLAGS+=-c $(GENCODE) -O3 -ccbin $(CC)
else
  NCCFLAGS+=-c $(GENCODE) -O0 -ccbin $(CC) -g
endif
ifeq ($(USE_FAST_MATH), YES)
  NCCFLAGS+= -ffast-math
endif

all:$(EXECS_DIRECTORY)/$(EXECNAME)
	rm -rf *.o
ifeq ($(USE_DEBUG), NO)
	@echo "  --- Compiled Release GPU version ---"
else
	@echo "  --- Compiled Debug GPU version ---"
	mv $(EXECS_DIRECTORY)/$(EXECNAME) $(EXECNAME)_debug
endif

$(EXECS_DIRECTORY)/$(EXECNAME):  $(OBJECTS)
	$(CC) $(OBJECTS) $(CCLINKFLAGS) -o $@ $(JLIBS)

.cpp.o: 
	$(CC) $(CCFLAGS) $< 

# The following rules compile the .hip files with hipcc
JSphGpu_ker.o: JSphGpu_ker.hip
	$(NCC) $(NCCFLAGS) JSphGpu_ker.hip

JSphGpuSimple_ker.o: JSphGpuSimple_ker.hip
	$(NCC) $(NCCFLAGS) JSphGpuSimple_ker.hip

JCellDivGpu_ker.o: JCellDivGpu_ker.hip
	$(NCC) $(NCCFLAGS) JCellDivGpu_ker.hip

JCellDivGpuSingle_ker.o: JCellDivGpuSingle_ker.hip
	$(NCC) $(NCCFLAGS) JCellDivGpuSingle_ker.hip

JSphShifting_ker.o: JSphShifting_ker.hip
	$(NCC) $(NCCFLAGS) JSphShifting_ker.hip

JDsAccInput_ker.o: JDsAccInput_ker.hip
	$(NCC) $(NCCFLAGS) JDsAccInput_ker.hip

JDsPips_ker.o: JDsPips_ker.hip
	$(NCC) $(NCCFLAGS) JDsPips_ker.hip

JDsGauge_ker.o: JDsGauge_ker.hip
	$(NCC) $(NCCFLAGS) JDsGauge_ker.hip

JWaveOrder2_ker.o: JWaveOrder2_ker.hip
	$(NCC) $(NCCFLAGS) JWaveOrder2_ker.hip

JReduSum_ker.o: JReduSum_ker.hip
	$(NCC) $(NCCFLAGS) JReduSum_ker.hip

JRelaxZone_ker.o: JRelaxZone_ker.hip
	$(NCC) $(NCCFLAGS) JRelaxZone_ker.hip

clean:
	rm -rf *.o $(EXECNAME) $(EXECNAME)_debug
