

/** Main CUDA file (contains also host code) */

#include "everything.h"
#include "stubParticle.h"

#define _STDLIB_H // hack
//#define _MATH_H

//#define HUMAN_READABLE_STATS
#ifndef GPU_REFILL_THRESHOLD
#define GPU_REFILL_THRESHOLD 10000
#endif

#ifdef STEP_PARALLEL
#include <cudpp.h>
#include "G4Navigator2.h"
#include "gpu2.c"
#else
#include "gpu.c"
#endif

#include "geometry.hpp"

// Debug output / utility data structures
typedef struct { const char *err, *fn; int line, errcode; } my_cuda_err;
typedef struct { int secs; int usecs; } mytimet;

// Debug output / utility functions, defined in cudamain.cpp
extern "C"
{
	void myprint( const char *chr );
	void myprint1( const char *chr, int n );
	mytimet mytimer();
	void myprinttdiff(mytimet a, mytimet b);
	void mysleep(int n);
}

// Error handling (and disabling for no-debug modes)
#ifndef NDEBUG
#define CHECKERR( thing ) do {\
	cudaError_t errc = thing; if (errc != cudaSuccess) {\
	my_cuda_err r = { cudaGetErrorString(errc), __FILE__, __LINE__, errc }; \
	return r; } } while(0)
	
#define CUDPPCHECKERR( thing ) do {\
	CUDPPResult errc = thing; if (errc != CUDPP_SUCCESS) {\
	my_cuda_err r = { "CUDPP error", __FILE__, __LINE__, errc }; \
	return r; } } while(0)
	
#define CHECKLASTERR CHECKERR(cudaGetLastError())
#else
#define CHECKERR( thing ) thing
#define CUDPPCHECKERR( thing ) thing
#define CHECKLASTERR (void)0
#endif

#define RETURNOK my_cuda_err ok = { NULL, NULL, 0, cudaSuccess }; return ok

// Utility functions
static inline int ceilDiv( int a, int d )
{
	return a/d + ((a%d)?1:0);
}

/*static inline int powerOf2LowerBound( int n )
{
	int m = 1;
	while( 2*m <= n ) m *= 2;
	return m;
}

static inline int powerOf2UpperBound( int n )
{
	int m = 1;
	while( m < n ) m *= 2;
	return m;
}*/

// Global variables used by both step and track parallel versions
Particle *gpuInput;
G4double *gpuOutput;
Geometry::byte *gpuGeom;
int numInput, numOutput, numInputPerRound;

#ifdef STEP_PARALLEL

G4Navigator *navs;
int *navkeys, *navvals, *newnavvals, *particleIndices;
int *navvalsCpu, *particleIndicesCpu;
int *gpuNumInput;
Particle *particlesCpu;
G4double *outputCpu;
CUDPPHandle plan;

#endif

// Global variables / constants used by the step parallel version
const int WARP_SIZE = 32;

void createGrid( int numInput, dim3* grid, dim3* block )
{
	// max grid size (not used)
	const int MAXSIZE = 10000000; //64512; // = 63*1024
	
	// Hard-coded settings for NVIDIA GeForce 470 GTX
	// TODO: fetch at runtime
	//
	const int NUMCORES = 448;
	const int NUMMULTIPROC = 14;
	const int BLOCKS_PER_MULTIPROC = 8;
	//const int MAXBLOCKS = NUMMULTIPROC*BLOCKS_PER_MULTIPROC;
	const int MAX_WARPS_PER_MULTIPROC = 48;
	const int MAX_DATA_PER_MULTIPROC = MAX_WARPS_PER_MULTIPROC*WARP_SIZE;
	
	int size = numInput;
	if (size > MAXSIZE) size = MAXSIZE;
	
	int dataPerMultiproc = ceilDiv(size,NUMMULTIPROC);
	if ( dataPerMultiproc > MAX_DATA_PER_MULTIPROC )
		dataPerMultiproc = MAX_DATA_PER_MULTIPROC;
	
	int blockSize = ceilDiv(dataPerMultiproc,BLOCKS_PER_MULTIPROC);
	
	const int MAX_BLOCK_SIZE = 1024;
	if (blockSize > MAX_BLOCK_SIZE) blockSize = MAX_BLOCK_SIZE;
	
	int numBlocks = ceilDiv(size,blockSize);
	int numWarps = ceilDiv(blockSize,WARP_SIZE) * numBlocks;
	
	if (numWarps > NUMCORES)
	{
		blockSize = ceilDiv(blockSize,WARP_SIZE)*WARP_SIZE;
		dataPerMultiproc = blockSize * BLOCKS_PER_MULTIPROC;
		if ( dataPerMultiproc > MAX_DATA_PER_MULTIPROC )
			blockSize -= WARP_SIZE;
	}
		
	size = blockSize*ceilDiv(size,blockSize);
	if (size > MAXSIZE) size = MAXSIZE;

	block->x = blockSize;
	block->y = block->z = 1;
	grid->x = size/blockSize;
	grid->y = 1;
	grid->z = 1;
}

/** Initialization, setting input */
my_cuda_err cudainit( Geometry *geom, int N )
{
	const mytimet t0 = mytimer();

	numOutput = numInput = numInputPerRound = N;
	
	CHECKERR( cudaSetDeviceFlags(0) );
	
	CHECKERR( cudaMalloc( (void**)&gpuInput, sizeof(Particle)*numInput ) );
	CHECKERR( cudaMalloc( (void**)&gpuOutput, sizeof(G4double)*numOutput ) );
	CHECKERR( cudaMalloc( (void**)&gpuGeom, geom->size() ) );

	geom->relocate( gpuGeom );
	
#ifdef STEP_PARALLEL

	myprint1("Iterations before compact: %d\n", ITERATIONS_BEFORE_COMPACT);

	unsigned hostAllocFlags = cudaHostAllocDefault;

	CHECKERR( cudaMalloc( (void**)&navs, sizeof(G4Navigator)*numInput ) );
	CHECKERR( cudaMalloc( (void**)&navkeys, sizeof(int)*numInput ) );
	CHECKERR( cudaMalloc( (void**)&navvals, sizeof(int)*numInput ) );
	CHECKERR( cudaMalloc( (void**)&newnavvals, sizeof(int)*numInput ) );
	CHECKERR( cudaMalloc( (void**)&gpuNumInput, sizeof(int) ) );
	
	CHECKERR( cudaMalloc( (void**)&particleIndices, sizeof(int)*numInput ) );
	CHECKERR( cudaHostAlloc( (void**)&particlesCpu, sizeof(Particle)*numInput, hostAllocFlags ) );
	CHECKERR( cudaHostAlloc( (void**)&outputCpu, sizeof(G4double)*numInput, hostAllocFlags ) );
	CHECKERR( cudaHostAlloc( (void**)&particleIndicesCpu, sizeof(int)*numInput, hostAllocFlags ) );
	CHECKERR( cudaHostAlloc( (void**)&navvalsCpu, sizeof(int)*numInput, hostAllocFlags ) );
	
	for ( int i=0; i<numInput; ++i ) particleIndicesCpu[i] = -1;
	
	CUDPPConfiguration conf;
	conf.algorithm = CUDPP_COMPACT;
	conf.datatype = CUDPP_INT;
	conf.options = CUDPP_OPTION_FORWARD;
	CUDPPCHECKERR( cudppPlan( &plan, conf, numInput, 1, 0 ) );
		
	cudaFuncSetCacheConfig(doStep, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(init, cudaFuncCachePreferL1);
#else
	cudaFuncSetCacheConfig(trace, cudaFuncCachePreferL1);
#endif
	
	CHECKERR( cudaMemcpy( gpuGeom, geom->getBuffer(), geom->size(), cudaMemcpyHostToDevice ) );

	const mytimet t1 = mytimer();
	myprint("Initialization: ");
	myprinttdiff(t0, t1);

	RETURNOK;
}

#ifdef STEP_PARALLEL

/** Main computations, step parallel version */
my_cuda_err cudaexec( G4double phys_step, int totalInput, Particle *input, G4double *output )
{
   int inputIndex = 0;
   numInput = 0;
   while( inputIndex < totalInput )
   {
	// Init. round
	int firstFree0 = numInput;
	int firstFree = firstFree0;

	// Assign particles to free slots
	int addedP = 0, unused = 0;
	for ( int i=0; i < numInputPerRound; ++i )
	{
		if ( particleIndicesCpu[i] < 0 )
		{
			if (inputIndex < totalInput)
			{
				int idx = inputIndex++;
				particleIndicesCpu[i] = -idx-2; // -2 => 0, -3 => 1, -4 => 2, ...
				particlesCpu[i] = input[idx];
				navvalsCpu[firstFree] = i;
				firstFree++;
				addedP++;
			}
			else
			{
				particleIndicesCpu[i] = -1; // leave unused
				unused++;
			}
		}
	}
	myprint1("Added %d particles, ", addedP);
	myprint1("unused: %d, ", unused);
	myprint1("inputIndex = %d/", inputIndex);
	myprint1("%d\n", totalInput);
	
	numInput = firstFree;
	numOutput = numInputPerRound;
	CHECKERR( cudaMemcpy( gpuNumInput, &numInput, sizeof(int), cudaMemcpyHostToDevice ) );
	
	// Transfer
	CHECKERR( cudaMemcpyAsync( gpuInput, particlesCpu, sizeof(Particle)*numInputPerRound, cudaMemcpyHostToDevice ) );
	CHECKERR( cudaMemcpyAsync( particleIndices, particleIndicesCpu, sizeof(int)*numInputPerRound, cudaMemcpyHostToDevice ) );
	CHECKERR( cudaMemcpyAsync( navvals+firstFree0, navvalsCpu+firstFree0, sizeof(int)*(firstFree-firstFree0), cudaMemcpyHostToDevice ) );

	dim3 grid, block;

	createGrid( numInput, &grid, &block );
	const int gridRowSize = grid.x*block.x;

	CHECKERR( cudaThreadSynchronize() );
	
	for (int i=0; i<numInput; i+=gridRowSize)
	{
		init<<<grid, block>>>( gpuInput, gpuOutput, navs, particleIndices, (G4VPhysicalVolume*)gpuGeom, i, gpuNumInput );
		CHECKLASTERR;
	}

	CHECKERR( cudaThreadSynchronize() );
		
	const int MIN_THRESHOLD = GPU_REFILL_THRESHOLD;
	while( numInput > 0 && (numInput > MIN_THRESHOLD || inputIndex == totalInput) )
	{
		createGrid( numInput, &grid, &block );
		
		#ifdef HUMAN_READABLE_STATS
		myprint1("nInput: %d, ", numInput);
		myprint1("gridLen: %d, ", gridRowSize);
		myprint1("blockSize: %d,\t", block.x);
		myprint1("nGrids: %d, ", ceilDiv(numInput,gridRowSize));
		myprint1("nBlocks: %d, ", grid.x);
		myprint1("nWarps: %d\n", ceilDiv(block.x,WARP_SIZE)*grid.x);
		#else
		myprint1("%d\n", numInput);
		#endif
		
		for (int i=0; i<numInput; i += gridRowSize )
		{
			doStep<<< grid, block >>>( gpuOutput, navs, navkeys, navvals, particleIndices, phys_step, i, gpuNumInput );
			CHECKLASTERR;
		}
		
		CHECKERR( cudaThreadSynchronize() );
		CUDPPCHECKERR( cudppCompact( plan, newnavvals, (size_t*)gpuNumInput, navvals, (const unsigned int*)navkeys, numInput ) );
		
		CHECKERR( cudaThreadSynchronize() );
		CHECKERR( cudaMemcpy( &numInput, gpuNumInput, sizeof(int), cudaMemcpyDeviceToHost ) );
	
		CHECKERR( cudaThreadSynchronize() );
		int *tmp = navvals;
		navvals = newnavvals;
		newnavvals = tmp;
	}
	
	//myprint1("numInput = %d\n", numInput);
	//myprint1("numOutput = %d\n", numOutput);
	
	// Finish round
	CHECKERR( cudaMemcpyAsync( outputCpu, gpuOutput, sizeof(G4double)*numOutput, cudaMemcpyDeviceToHost ) );
	CHECKERR( cudaMemcpyAsync( particleIndicesCpu, particleIndices, sizeof(int)*numOutput, cudaMemcpyDeviceToHost ) );
	CHECKERR( cudaThreadSynchronize() );
		
	for ( int i=0; i<numOutput; ++i )
	{
		if ( particleIndicesCpu[i] < -1 )
		{
			output[ -particleIndicesCpu[i] - 2 ] = outputCpu[i];
		}
	}
   }

   RETURNOK;
}

#else

/** Main computations track parallel version*/
my_cuda_err cudaexec( G4double phys_step, int totalInput, Particle *input, G4double *output )
{
   for ( int i = 0; i < totalInput; i += numInput )
   {
	// Init
	if ( i + numInput > totalInput ) numInput = totalInput-i;
	CHECKERR( cudaMemcpy( gpuInput, input+i, sizeof(Particle)*numInput, cudaMemcpyHostToDevice ) );

	dim3 grid, block;
	createGrid( numInput, &grid, &block );

	// Compute
	trace <<< grid, block >>>( gpuInput, gpuOutput, (G4VPhysicalVolume*)gpuGeom, phys_step, numInput, NULL, NULL, NULL, NULL, NULL );
	CHECKLASTERR;
	
	// Finish round
	CHECKERR( cudaMemcpy( output+i, gpuOutput, sizeof(G4double)*numOutput, cudaMemcpyDeviceToHost ) );
   }

   RETURNOK;
}

#endif

/** Finalization, fetching output */
my_cuda_err cudafinish()
{
	const mytimet t0 = mytimer();
	
	CHECKERR( cudaFree( gpuInput ) );
	CHECKERR( cudaFree( gpuOutput ) );
	CHECKERR( cudaFree( gpuGeom ) );
	
#ifdef STEP_PARALLEL
	CUDPPCHECKERR( cudppDestroyPlan(plan) );
	
	CHECKERR( cudaFree( navs ) );
	CHECKERR( cudaFree( navkeys ) );
	CHECKERR( cudaFree( navvals ) );
	CHECKERR( cudaFree( newnavvals ) );
	CHECKERR( cudaFree( gpuNumInput ) );
	CHECKERR( cudaFree( particleIndices ) );
	CHECKERR( cudaFreeHost( particlesCpu ) );
	CHECKERR( cudaFreeHost( outputCpu ) );
	CHECKERR( cudaFreeHost( particleIndicesCpu ) );
	CHECKERR( cudaFreeHost( navvalsCpu ) );
#endif

	CHECKERR( cudaThreadExit() );

	const mytimet t1 = mytimer();
	myprint("Finalization: ");
	myprinttdiff(t0, t1);

	RETURNOK;
}
