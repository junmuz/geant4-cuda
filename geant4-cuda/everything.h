
/** Device keyword & type definitions, some constant definitions */

#ifndef EVERYTHING_H__
#define EVERYTHING_H__

// ---------- Generic target platform defines
#ifdef CUDA_CODE
	#define GPU		// GPU related code (host or device)
	#define GPU_CODE	// GPU device code
	#define CUDA		// CUDA code (host or device)
#elif defined(OPENCL_CODE)
	#define GPU
	#define GPU_CODE
	#define OPENCL		// OpenCL code (host or device)
#elif defined(CUDA_HOST)
	#define GPU
	#define HOST_CODE	// host code
	#define CUDA
#elif defined(OPENCL_HOST)
	#define GPU
	#define HOST_CODE
	#define OPENCL
#elif defined(CPU_HOST)
	#define HOST_CODE
#endif

// ---------- address space and function qualifiers (GPU)
#ifdef OPENCL_CODE

	#define GLOBALFUNC __kernel		
				// host callable kernel function
	#define GLOBALTYPE __global		
				// global (slow/big) device memory
#if ( GLOBAL_MODE !=1)
	#define SHAREDTYPE __local
#else
	#define SHAREDTYPE __global	
		// shared (small/fast) device memory
#endif	
#define PRIVATEMEM __private
#define SHAREDMEM __local
#define SHAREDSHADOW __local
	#define CONSTTYPE   __constant		
               // constant device memory
	//EDIT
	//#define INLINEFUNC static inline
	#define INLINEFUNC inline
               // inline functions
    // Barrier functions
    #define BARRIER_LOCAL barrier( CLK_LOCAL_MEM_FENCE) 
    #define BARRIER_ALL barrier( CLK_GLOBAL_MEM_FENCE) 
#if (GLOBAL_MODE ==1)
    #define BARRIER_FLEXIBLE barrier( CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE) 
#else
    #define BARRIER_FLEXIBLE barrier( CLK_LOCAL_MEM_FENCE );
#endif
#elif defined(CUDA_CODE)

	#define GLOBALFUNC __global__
	#define GLOBALTYPE
	#define SHAREDTYPE __shared__
	#define PRIVATEMEM
	#define SHAREDMEM __shared__
	#define SHAREDSHADOW 
	#define CONSTTYPE const
	#define NULL 0
	#define GNULL 0
	#define INLINEFUNC __device__
    #define BARRIER_LOCAL __syncthreads()
    #define BARRIER_ALL __syncthreads()
    #define BARRIER_FLEXIBLE __syncthreads()

#else

	#define GLOBALFUNC
	#define GLOBALTYPE
	#define SHAREDTYPE
	#define CONSTTYPE
	//#define INLINEFUNC static inline
	//EDIT
	#define INLINEFUNC inline
	
#endif

// ---------- current thread ID (OpenCL / TODO: CUDA)
#ifdef CUDA_CODE
	#define get_global_id(X) (blockIdx.x * blockDim.x + threadIdx.x)
	#define get_local_id(X) threadIdx.x
	#define get_local_size(X) blockDim.x
#endif

// ---------- inline functions


#ifdef INLINE_EVERYTHING
#define MAYINLINE INLINEFUNC
#else
#define MAYINLINE
#endif

// ---------- standard C library components
#ifdef GPU_CODE
	#define myAssert(x) (void)0
	#define myAbort() (void)0
#else
	#define myAssert(x) assert(x)
	#define myAbort() abort()
	#include <assert.h>
	#include <math.h>
	#include <stdlib.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

// ---------- NULL pointer
#ifndef NULL
#define NULL ((void*)0)
#endif
#ifndef GNULL
#define GNULL ((GLOBALTYPE void*)0)
#endif
#ifndef CNULL
#define CNULL ((CONSTTYPE void*)0)
#endif

// ---------- C++ boolean type
#if !defined(__cplusplus) && !defined(GPU)
	#define true 1
	#define false 0
	typedef int bool;
#endif

// ---------- Numeric limits
#ifndef DBL_MAX
#ifdef OPENCL_CODE
#define DBL_MAX MAXFLOAT
#else
#ifdef DOUBLE_PRECISION
#define DBL_MAX 1.7976931348623157e+308 // TODO: check for CUDA!
#else
#define DBL_MAX 1e37
#endif
#endif
#endif

// ---------- basic Geant4 types
#ifdef DOUBLE_PRECISION
	#ifdef OPENCL_CODE
		#pragma OPENCL EXTENSION cl_khr_fp64: enable
	#endif
	typedef double G4double;
#else
	typedef float G4double;
#endif

typedef float G4float;
typedef int G4int;
typedef int G4bool;
typedef long G4long;	

// `Infinity' - Distance returned for no intersection etc.
#ifdef DOUBLE_PRECISION
static CONSTTYPE const G4double kInfinity = 9.0E99;
#else
//EDIT: Does not work on OpenCL
//static CONSTTYPE const G4double kInfinity = 1.0E37;
CONSTTYPE G4double kInfinity = 1.0E37;
#endif

CONSTTYPE int BlockSize = 32;
// The size of the block to be run on the GPU.
CONSTTYPE int Multiplier = 4;
// The multiplier factor for the size of Shared Memory to be used to store Solid information. See G4NewNavigation.c for more information.
CONSTTYPE G4double twopi = 2.0*M_PI;

// Minimum cosine of angle between surface normal & track direction
// for exiting normal optimisation
CONSTTYPE G4double kMinExitingNormalCosine = 1E-3;

// G4VSolid::Inside return codes
// kSurface => within tolerance of exact surface
typedef enum {kOutside,kSurface,kInside} EInside;

// kNormal = (G4PVPlacement) Conventional positioning
// kReplica = (G4PVReplica)  Consumed parameterised case
//                           => Distances & location computed with
//                              simple formulae & MOTHER volume(s)
//                              must also be checked
// kParameterised = (G4PVParameterised) General parameterised volume
//                           => Distance & location computed to volumes
//                              after setup/modification via user object
typedef enum {kNormal,kReplica,kParameterised} EVolume;

typedef enum {kXAxis,kYAxis,kZAxis,kRho,kRadial3D,kPhi,kUndefined} EAxis;

typedef enum { kBox = 0 , kOrb, kTubs, kCons, kPolyCone, Solidcount } ESolid;
		// Currently supported solids. When a new solid type is defined, add to this list. Don't forget to increment solid count.

// ---------- Geant4 constants
#define K_GEOMETRY_ANG_TOLERANCE 1E-9 //*rad
#ifdef DOUBLE_PRECISION
#define K_GEOMETRY_CAR_TOLERANCE 1E-9 //*mm
#define K_GEOMETRY_RAD_TOLERANCE 1E-9 //*mm
#else
#define K_GEOMETRY_CAR_TOLERANCE 1E-3 //*mm
#define K_GEOMETRY_RAD_TOLERANCE 1E-3 //*mm
#endif
#define K_MAX_VOXEL_NODES 1000
#define K_MIN_VOXEL_VOLUMES_LEVEL_1 2
#define K_MIN_VOXEL_VOLUMES_LEVEL_2 3
#define K_MIN_VOXEL_VOLUMES_LEVEL_3 4
#define K_NAVIGATOR_ACTION_THRESHOLD_NOZEROSTEPS 10
#define K_NAVIGATOR_ABANDON_THRESHOLD_NOZEROSTEPS 25

// ---------- placement of Geant4 components in memory
// Warning: only GLOBAL-something (currently) works 
#define GEOMETRYLOC GLOBALTYPE	// Location of geometry data
#define GEOMETRYNULL GNULL	// NULL pointer in geometry address space

// ---------- enable "combined navigation" on the GPU

// EDIT : Adding new qialifier for combined naigation
#if defined(ENABLE_COMBINED_NAVIGATION) && defined(GPU_CODE)
#define ENABLE_COMBINED_NAVIGATION
#endif

#endif
