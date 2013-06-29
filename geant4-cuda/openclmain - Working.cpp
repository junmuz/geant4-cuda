/*
Main file needed for execution of openCL port of Geant4 GPU code.

Used in close conjunction with cl.hpp and hostcommons.hpp
The main function is the handler class responsible for allocation and deallocation of memory on GPU ( and page locked mem on CPU), as well as setting kernel arguments and
execution of kernels. Memory handling uses the Buffer class from cl.hpp and kernels are defined in gpu.c which is then compiled into prog.cl (see Makefile).

Notes For Developers:
Some of the typical calls needed to allocate and transfer memory and other steps in a typical kernel call are defined below. These are only meant as a guideline so and do not ensure
proper functioning of the kernel. However, they may serve as a checklist to refer to before adding a new kernel. For proper definitions and notes of these functions, see cl.hpp
If a new context adn device queue are to be created then call SingleFileSingleGpuSetup again with new arguments.

1) Memory allocation and transfer(a-> Page Locked; b-> Normal Buffer type memory)
	a) CL::PinnedBufferPair gpuInput( CommandQueue& queue, size_t size, cl_mem_flags mem_flags, cl_map_flags map_flags);
	  Memory transfer can be done by using the transfertoDevice or transferFromDevice methods. The initialization of a PinnedBufferPair handles the creation of required memory
	  on both the device and the host. The host accessible pointer can be accessed with the getHostPtr function and the buffer in device memory can be accessed with the
	  getDeviceBuffer method. 
	  NOTE: The current implementation transfers the memory to the pointer returned from getHostPtr using a std::memcpy. 
	  //MODIFY
	  
	b) CL::Buffer gpuGeom( Context &context, cl_mem_flags flags, size_t size );
	  Data transfer fo Buffers is handled by the instance of SingleFileSingleGpuSetup. Default, use gpusetup.
	  NOTE: In most cases the default gpusetup can be passed in as the first argument.

2) Creating and executing the Kernel. 
	CL::Kernel kernelName( program &p, std::string funcname)
	NOTE: Here funcname is the name of the kernel passed in as a string. The kernel has to be defined in prog.cl. Example usage CL::KERNEL kernelTrace( gpusetup, "trace");
	Kernel argumets can be set by using the setArg method. Finally exection of kernel is done be using the  enqueueKernel( kernelname, numInput, blockSize ) method from an
	instance of SingleFileSingleGPUSetup class ( default uses gpusetup).
	In case kernel has to be executed as a task and not NDRange then use the function enqueueTask

Original Author -: Otto Seiskari


Changes by - Dhruva T.B


*/
/** Main OpenCL host file */
#define OPENCL_HOST

// Defining all parameters from make file here. Fix for Visual Studio Windows.
#define ENABLE_VOXEL_NAVIGATION
#define INLINE_EVERYTHING
#define NTHREADS 12
#define	ITERATIONS_BEFORE_COMPACT 100
#define	K_NAVIGATION_HISTORY_DEPTH 16
#define VERBOSE
//EDIT:2
#define PHYSICS 1

#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <fstream>
#include <cmath>

#include "utils.h"
#include "cl.hpp"

#include "hostcommons.hpp"

#include "G4Navigator.h"

// CL compiler flags

#ifdef VERBOSE
	#define CL_VERBOSE "-cl-single-precision-constant"
#else
	#define CL_VERBOSE ""
#endif

#define CL_WERROR //" -Werror"

#ifdef FAST_CL_MATH
	#define CL_OPTIMIZATIONS " -cl-fast-relaxed-math -cl-mad-enable"
#else
	#define CL_OPTIMIZATIONS
#endif


#define COMPILER_FLAGS \
	CL_OPTIMIZATIONS \
	CL_WERROR \
	CL_VERBOSE

//REMOVE
void print_stuff( void * pointer, int size)
{   int globalidx =0;
   for( globalidx=0;globalidx<size ;globalidx++)
	{int destoffs = *((int*)pointer + 2*globalidx);
    int targoffs = *((int*)pointer + 2*globalidx+1);
	// Addressability issue?
	//std::cout<<"i= "<< globalidx <<" and destination offset = "<<destoffs<<" and targoffs = "<<targoffs<<std::endl;
    }
}

int main( int argc, char *argv[] )
{
#ifndef NOCATCH
	try
	{
#endif
		// Initialize geometry and input/output buffers
		TestCase testCase( argc, argv );
#ifndef PHYSICS
		int numInput = testCase.xres * testCase.yres;
#else
		int numInput = testCase.getSize();
#endif
		
		const int numOutput = numInput;
		
		const int blockSize = 32;
		
		const char *BINARY_FILE_NAME = "prog.ptx";
		const char *SOURCE_FILE_NAME = "prog.cl";
		
		// Command line handling
		bool isBinary = false;
		
		for ( int i=1; i<argc; ++i )
		{
			const std::string arg( argv[i] );
			if ( arg == "-b" )
				isBinary = true;
			/*else
				throw std::runtime_error( "Invalid option "+arg );*/
		}
		
		// Load or compile program
		std::cerr << (isBinary ? "Loading" : "Compiling") << "...";
		//EDIT: Removing all of his clock functions, can use the VS Porfiler if necessary.
		//const my_clock_t tc0 = my_clock();
		
		CL::SingleFileSingleGPUSetup gpuSetup(
			isBinary ? BINARY_FILE_NAME : SOURCE_FILE_NAME,
			isBinary, COMPILER_FLAGS );
		
		//const my_clock_t tc1 = my_clock();
		std::cout<<"SETUP complete"<<std::endl;

		//std::cerr << "done in " << tdiff(tc0,tc1) << " seconds\n\n";
		
		#ifdef VERBOSE
		std::cerr << " ---- Build log\n" << gpuSetup.getBuildLog() << "\n";
		#endif
		
		// Save "binary" (PTX bytecode) for reuse
		if ( !isBinary )
		{
			std::ofstream f( BINARY_FILE_NAME );
			gpuSetup.writeBinary( f );
		}
		
		//EDIT: Changing code to get rid of kernel getptr
		//1) Define new kernel relocate 
		//2) Use function getPtrs to get the vector ptrs and use it to define size (ptrs.size()) and offset (sizeof(int)) on host.
		//3) Allocate some mem on GPU to store ptrs and use it in the relocate kernel.
		//4) Move the enqueue write here so that the WorldVolumePointer is defined.
		//5) As soon as the kernel is done, free the memory used for ptrs.
		// Import handles to OpenCL kernels (functions)

		CL::Kernel kernelTrace( gpuSetup, "trace" );
		CL::Kernel kernelRelocate( gpuSetup, "relocate" );
		//EDIT: New kernel check
		CL::Kernel kernelCheck ( gpuSetup, "check");

		// Reserve GPU & host buffers
		//EDIT Getting size of ptrs
		int size = testCase.geom->ptrs_size();
		int size_of_logical_checks = 1;
		
		const int auxBufSz = sizeof(cl_mem);
		
		// Page-locked buffers for fast DMA-IO
		CL::PinnedBufferPair gpuInput( gpuSetup, numInput*sizeof(StubParticle), CL_MEM_READ_WRITE, CL_MAP_WRITE );
		CL::PinnedBufferPair gpuOutput( gpuSetup, numOutput*sizeof(G4double), CL_MEM_WRITE_ONLY, CL_MAP_READ );
		CL::PinnedBufferPair gpuAux( gpuSetup, auxBufSz, CL_MEM_READ_WRITE, CL_MAP_READ | CL_MAP_WRITE );
		//EDIT
		CL::PinnedBufferPair ptrs( gpuSetup, size*2*sizeof(int), CL_MEM_READ_WRITE, CL_MAP_WRITE );
		CL::PinnedBufferPair result( gpuSetup,size_of_logical_checks*sizeof(cl_mem), CL_MEM_WRITE_ONLY, CL_MAP_READ );
			
		std::cout<<"Pinned Buffers allocation complete"<<std::endl;
		// GPU only buffers
		//EDIT
		
		CL::Buffer gpuGeom( gpuSetup, CL_MEM_READ_WRITE, testCase.geom->size());
		//EDIT2:
		//gpuSetup.enqueueWriteBuffer( gpuGeom, testCase.geom->getBuffer() );
		
		
		
		std::cout<<"Device Buffers allocation complete"<<std::endl;

		std::memcpy( ptrs.getHostPtr(), &(testCase.geom->ptrs[0]), size*2*sizeof(int) );
		//EDIT 2:
	   //print_stuff(ptrs.getHostPtr(), size);
		ptrs.transferToDevice();
		gpuSetup.finish();

		
		kernelCheck.setArg( 0 , gpuGeom);
		kernelCheck.setArg(	1, result.getDeviceBuffer());
		gpuSetup.enqueueKernel(kernelCheck, 1, 1);
		gpuSetup.finish();
		//EDIT
		result.transferFromDevice();
		gpuSetup.finish();
		int answer = *(int*)(result.getHostPtr());
		std::cout<<"Result before relocate is "<< answer<< std::endl;
	/*EDIT2
		kernelRelocate.setArg( 0, ptrs.getDeviceBuffer());
		kernelRelocate.setArg (1, gpuGeom);
		kernelRelocate.setArg(2, sizeof(int), &size);
    */
		//MODIFY
		//gpuSetup.enqueueKernel(kernelRelocate,size*2,2);
		//gpuSetup.finish();
		//OLD
		// Fetch address of gpuGeom in device memory space (kludge)
		//kernelGetPtr.setArg( 0, gpuGeom );
		//kernelGetPtr.setArg( 1, gpuAux.getDeviceBuffer() );
		//gpuSetup.enqueueTask( kernelGetPtr );
		//gpuAux.transferFromDevice( CL_TRUE, 0, sizeof(cl_mem) );
		//cl_mem gpuhandle = *(cl_mem*)gpuAux.getHostPtr();
		//EDIT2
		cl_mem gpuhandle = *(cl_mem*)result.getHostPtr();
		#ifdef VERBOSE
		std::cerr << "GPU geometry location: " << gpuhandle << "\n\n";
		#endif
		std::cout<< "About to run relocate, no problem so far\n";
		
		//OLD:EDIT2
		testCase.geom->relocate( gpuhandle );

		//EDIT2
		/*kernelCheck.setArg( 0 , gpuGeom);
		kernelCheck.setArg(	1, result.getDeviceBuffer());
		gpuSetup.enqueueKernel(kernelCheck, 1, 1);
		gpuSetup.finish();*/
		//EDIT
		result.transferFromDevice();
		gpuSetup.finish();
		int answer2 = *(int*)(result.getHostPtr());
		//std::cout<<"Result is "<< answer2<< std::endl;
		std::cout<<"size of pointer to void on GPU is "<< answer2<<" and on CPU it is "<<sizeof(void *)<< std::endl;

		// Generate particles (input) ... kind of stupid way to do this
		//EDIT2:

		// MODIFY: He also notes how this is not the optimal way of doing this. Memory can be saved here.
		std::memcpy( gpuInput.getHostPtr(), &(testCase.input[0]), gpuInput.size() );
		//EDIT2
		gpuSetup.enqueueWriteBuffer( gpuGeom, testCase.geom->getBuffer() );
		
		gpuInput.transferToDevice();
		gpuSetup.finish();
		// Set GPU kernel arguments
		kernelTrace.setArg( 0, gpuInput.getDeviceBuffer() );
		kernelTrace.setArg( 1, gpuOutput.getDeviceBuffer() );
		kernelTrace.setArg( 2, gpuGeom );
		kernelTrace.setArg( 3, sizeof(G4double), &(testCase.phys_step) );
		kernelTrace.setArg( 4, sizeof(cl_int), &numInput );
		
		
			std::cout<<"Arguments set"<<std::endl;
		// Write input to GPU memory
		// const my_clock_t t1 = my_clock();
		//OLD
		/*gpuSetup.enqueueWriteBuffer( gpuGeom, testCase.geom->getBuffer() );*/
		
		
		std::cout<< "Write complete, transfer done, finish\n";

		// Actual execution
		//const my_clock_t t2 = my_clock();
		gpuSetup.enqueueKernel( kernelTrace, numInput, blockSize );
		gpuSetup.finish();
		std::cout<<"Kernel trace done\n";
		
		// Transfer results back to host memory
		//const my_clock_t t3 = my_clock();
		gpuOutput.transferFromDevice();
		gpuSetup.finish();
		//const my_clock_t t4 = my_clock();
		
		// Print time summary
		//std::cerr << "Elapsed: " << tdiffms(t1,t4) << " ms"
			//<< "\n  Transfer: " << tdiffms( t1, t2 )+tdiffms(t3,t4)
			//<< "\n\tto GPU:\t" << tdiffms( t1, t2 )
			//<< "\n\tfrom GPU:\t" << tdiffms( t3, t4 )
			//<< "\n  Calculation: " << tdiffms( t2, t3 ) << "\n\n";

		//MODIFY
		// Output results (also a stupid copy)
		std::memcpy( &(testCase.output[0]), gpuOutput.getHostPtr(), gpuOutput.size() );
		testCase.outputData( "imggpu.txt" );

		return EXIT_SUCCESS;
		
#ifndef NOCATCH
	}
	catch ( const std::runtime_error &e )
	{
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
#endif
}
