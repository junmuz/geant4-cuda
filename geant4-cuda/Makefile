
# Description of different targets:
#
#	cpumain (default):	CPU/OpenMP version (track-parallel)
#	opencl:			OpenCL version (track-parallel)
#	cuda:			CUDA version (track or step parallel)
#

# --- Header and library directories

# The SDK root is, by default, in someone's home directory (e.g. $(HOME))
WHICHHOME=$(HOME)
# Path to NVIDIA GPU Computing SDK directory
SDKROOT=$(WHICHHOME)/NVIDIA_GPU_Computing_SDK
# Path to OpenCL headers
OPENCLINC = /opt/AMDAPP/include/
# Path to CUDPP headers
CUDPPINC=-I$(SDKROOT)/C/common/inc/cudpp
# Path to Boost headers
BOOSTINC=-I/home/junaid/Workspace/GPU/Geant4/g4gpu-code
# Path to CUDA libraries (e.g. /usr/local/cuda/lib)
CUDALIBDIR=/usr/local/cuda/lib
# Path to CUDPP libraries
CUDPPLIBDIR=$(SDKROOT)/C/common/lib/linux
# CUDPP library name (e.g. -lcudpp_i386)
CUDPPLIB=-lcudpp_i386
# All OpenCL code related header directories
OPENCLINCS=-I./OpenCL/common -I$(OPENCLINC)

# --- Command-line overrideable configuration parameters

# Number of OpenMP threads in the CPU version
NTHREADS=12
# Number of steps to compute before compact (track-parallel gpu versions)
ITERATIONS_BEFORE_COMPACT=100
# Maximum depth of G4 navigation history
K_NAVIGATION_HISTORY_DEPTH=16
# Enable voxel navigation?
ENABLE_VOXEL_NAVIGATION=0
# Track or step parallel?
STEP_PARALLEL=0
# Physics example
PHYSICS=0
# Enable/disable debugging
DEBUGLEVEL=0

# --- Translating parameters to preprocessor directives...

# Common flags //-DDOUBLE_PRECISION \
COMMOND=
	-DNTHREADS=$(NTHREADS) \
	-DITERATIONS_BEFORE_COMPACT=$(ITERATIONS_BEFORE_COMPACT) \
	-DK_NAVIGATION_HISTORY_DEPTH=$(K_NAVIGATION_HISTORY_DEPTH) \
	-DVERBOSE \
	-Wall -Wextra

# Debug & profile flags
DEBUGFLAGS=

# CPU flags (fast C99 "device" code)
CPUD=-DCPU_CODE -DINLINE_EVERYTHING -Wno-unused-function

ifeq ($(ENABLE_VOXEL_NAVIGATION),1)
	COMMOND:=$(COMMOND) -DENABLE_VOXEL_NAVIGATION
endif

ifeq ($(PHYSICS),1)
	COMMOND:=$(COMMOND) -DPHYSICS
endif

ifeq ($(DEBUGLEVEL),1)
	DEBUGFLAGS:=$(DEBUGFLAGS) -g -pg
else
	CPUD:=$(CPUD) -DNDEBUG -O3
endif

ifneq ($(NTHREADS),1)
	CPUD:=$(CPUD) -DUSE_OPENMP -fopenmp
endif

COMMOND:=$(COMMOND) $(DEBUGFLAGS)

# CUDA flags
CUDAD=-DINLINE_EVERYTHING -DCUDA_CODE # Device code
CUDAL=-L$(CUDALIBDIR) -lcudart # Host code

ifeq ($(STEP_PARALLEL),1)
	CUDAD:=$(CUDAD) -DSTEP_PARALLEL $(CUDPPINC)
	CUDAL:=$(CUDAL) $(CUDPPLIB) -L$(CUDPPLIBDIR)
endif

# OpenCL flags (reduced number of different solids)
OPENCLD= -DINLINE_EVERYTHING -DONLY_BOX_AND_ORB #-DFAST_CL_MATH

# --- Targets

# CPU/OpenMP version
cpumain:
	$(CC) -c -o cpuexec.o cpuexec.c -std=c99 $(CPUD) $(COMMOND)
	$(CXX) -o cpumain cpumain.cpp cpuexec.o $(COMMOND) -fopenmp -DCPU_HOST

# OpenCL version
prog.cl: gpu.c
	$(CC) -E -P gpu.c -o prog.cl $(COMMOND) $(OPENCLD) -DOPENCL_CODE 
openclmain: prog.cl
	$(CXX) openclmain.cpp -o openclmain $(COMMOND) $(OPENCLINCS) $(OPENCLD) -lOpenCL -DOPENCL_HOST

# CUDA version	
cuda:
	$(CC) -E -P cuda.cpp -o prog.cu $(COMMOND) $(CUDAD)
	nvcc -c -arch=sm_20 -Xptxas -dlcm=ca prog.cu -DNDEBUG $(DEBUGFLAGS)
	g++ -o cudamain cudamain.cpp $(COMMOND) prog.o $(CUDAL) -DCUDA_HOST $(BOOSTINC)

# Clean all
clean:
	rm -rf *.o cpumain openclmain cudamain prog.* prog2.* imggpu.txt imgcpu.txt imgcuda.txt *.log
