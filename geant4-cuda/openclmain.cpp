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
	Kernel argumets can be set by using the setArg method. Finally exection of kernel is done be using the  enqueueKernel( kernelname, numInput, BlockSize ) method from an
	instance of SingleFileSingleGPUSetup class ( default uses gpusetup).
	In case kernel has to be executed as a task and not NDRange then use the function enqueueTask
	
Original Author -: Otto Seiskari


Changes by - Dhruva T.B


*/
/** Main OpenCL host file */
#define OPENCL_HOST


#include "gpuconf.h"

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

void check_navigation( void * pointer, std::vector<CheckPointer> check_pointer, int number_of_increments)
{   /*
		TYPE of check
		1-> Distance check.
		2-> Geometry check
		3-> Position check
		4-> Run all checks.
	*/
	float step;
	//Just one check to see if the distance travelled is consistent
	
	std::cout<<"--------------------------------------------------------------------------------"<<std::endl;

	if( (CHECK == 1 || CHECK == 4) && (number_of_increments == 0))
	{	
		float * ptr = (float * )pointer;
		std::cout<<"---Running Check for Distance---"<<std::endl;
		std::cout<<"\n\n";
		std::cout<<"Positions of particles before steps:"<<std::endl;
		for (int i=10;i<18;i=i+3)
		{
			std::cout<<"Particle "<<(i-10)/3<<" has X-: "<< ptr[i]<<" Y-: "<< ptr[i+1]<<" Z-: "<<ptr[i+2]<<std::endl;
		}
		std::cout<<"\n\n";
		std::cout<<"Positions of particles after steps:"<<std::endl;
		for (int i=20;i<28;i=i+3)
		{
			std::cout<<"Particle "<<(i-20)/3<<" has X-: "<< ptr[i]<<" Y-: "<< ptr[i+1]<<" Z-: "<<ptr[i+2]<<std::endl;
		}  
		std::cout<<"\n\n";
		std::cout<<"Directions of particles during the step:"<<std::endl;
		for (int i=30;i<38;i=i+3)
		{
			std::cout<<"Particle "<<(i-30)/3<<" has X-: "<< ptr[i]<<" Y-: "<< ptr[i+1]<<" Z-: "<<ptr[i+2]<<std::endl;
			
		}   
		std::cout<<"\n\n";

		for (int i=0;i<3;i++)
		{	
			float a = ptr[i*3 + 10]-ptr[i*3 + 20];
			float b = ptr[i*3 + 11]-ptr[i*3 + 21];
			float c = ptr[i*3 + 12]-ptr[i*3 + 22];
			float distance = sqrt(a*a + b*b + c*c);
			float step = ptr[i];
			float allowance = 0.01;
			//std::cout<<"At step number "<<i<<" the value of the step was "<<ptr[i]<<std::endl;
			assert(abs(distance-step)<allowance);
			std::cout<< "Difference between points = "<<distance<<" and step size returned was-> "<<step<<std::endl;
			
		}

		std::cout<<"The total number of steps that thread 0 took was "<< ptr[29] + 1<<std::endl;
		
		
		std::cout<<"\n---Check for Distance Passed---\n"<<std::endl;
	}

	else if ( CHECK == 2 || CHECK==4 )
	{
		unsigned long * ptr = (unsigned long*)pointer;
		std::cout<<"---Running Check for Geometry---"<<std::endl;
		int i;
		std::cout<<"\n\n";
		for( i=0; i<number_of_increments*3; i+=3)
		{
			std::cout<<"Count of Physical Volume = "<<ptr[i]<<" and corresponding offset for flogical = "<<ptr[i+1]<<" and property returned = "<<ptr[i+2]<<std::endl;
			int j;
			for(j=0; j < check_pointer[0].Phys_vector_pointer->counter ;j++)
			{
				if( check_pointer[j].Phys_vector_pointer->count == ptr[i])
				{	
					std::cout<<"Actual Value of Offset = "<< (unsigned long) check_pointer[j].offset <<" and that of Property = "<<check_pointer[j].property<<std::endl;
					assert( (check_pointer[j].offset) == ptr[i+1]); 
					assert( (check_pointer[j].property) == ptr[i+2]); 
					break;
				}
				
			}
		}
		std::cout<<"\n\n";
				
		std::cout<<"\n---Check for Geometry Passed---\n"<<std::endl;

	}
	
	else if(CHECK==3 || CHECK==4)
	{
		std::cout<<"---Running Check for Position---"<<std::endl;
	}

	std::cout<<"--------------------------------------------------------------------------------"<<std::endl;
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
		
		
		int number_of_increments = 5;
		

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
			// Kernel trace is the main kernel which does the navigation.
		CL::Kernel kernelRelocate( gpuSetup, "relocate" );
			// The kernel version of relocate. It is possible to run relocate without using the kernel straight on the host.
		CL::Kernel kernelTest( gpuSetup, "test");
			// Kernel to check for inconsistencies. 
		
		//EDIT: New kernel check
		CL::Kernel kernelCheck ( gpuSetup, "check");
			
		
		//EDIT2
		CL::Kernel kernelCheckGeometry(gpuSetup, "checkgeom");
			// Kernel to confirm gemoetry relocation happened as it should.	

		// Reserve GPU & host buffers
		//EDIT Getting size of ptrs
		int size = testCase.geom->ptrs_size();
		//REMOVE:
		std::cout<<"Size of ptrs is = "<<size<<std::endl;

		int size_of_logical_checks = 1000;
		
		const int auxBufSz = sizeof(cl_mem);
		
		// Page-locked buffers for fast DMA-IO
		CL::PinnedBufferPair gpuInput( gpuSetup, numInput*sizeof(StubParticle), CL_MEM_READ_WRITE, CL_MAP_WRITE );
		CL::PinnedBufferPair gpuOutput( gpuSetup, numOutput*sizeof(G4double), CL_MEM_WRITE_ONLY, CL_MAP_READ );
		CL::PinnedBufferPair gpuAux( gpuSetup, auxBufSz, CL_MEM_READ_WRITE, CL_MAP_READ | CL_MAP_WRITE );
		//EDIT
		CL::PinnedBufferPair ptrs( gpuSetup, size*2*sizeof(int), CL_MEM_READ_WRITE, CL_MAP_WRITE );
		CL::PinnedBufferPair result( gpuSetup,size_of_logical_checks*sizeof(cl_mem), CL_MEM_WRITE_ONLY, CL_MAP_READ );
			
		#if (GLOBAL_MODE ==1)
		// sizeof( size_t ) ?
		CL::Buffer Numbers_Of_Solid( gpuSetup, CL_MEM_READ_WRITE, numInput*sizeof(int));
		
		CL::Buffer Sum_Of_Solid( gpuSetup, CL_MEM_READ_WRITE, numInput*sizeof(int));
		
		CL::Buffer Solids( gpuSetup, CL_MEM_READ_WRITE, numInput*sizeof(SolidInfo));
		
		CL::Buffer Result_For_Current_Solid( gpuSetup, CL_MEM_READ_WRITE, numInput*sizeof(ResultInfo));
		
		CL::Buffer Compacter_Result( gpuSetup, CL_MEM_READ_WRITE, numInput*sizeof(FinalResult));

        #endif
		CL::Buffer nullVNode( gpuSetup, CL_MEM_READ_WRITE, 2 *sizeof(G4SmartVoxelNode ));
		CL::Buffer noStepArray ( gpuSetup, CL_MEM_READ_WRITE, numInput*sizeof ( bool ));
		CL::Buffer LocationArray( gpuSetup, CL_MEM_READ_WRITE, numInput * sizeof( PointInformation));
		std::cout<<"Pinned Buffers allocation complete"<<std::endl;
		// GPU only buffers
		//EDIT
		
		CL::Buffer gpuGeom( gpuSetup, CL_MEM_READ_WRITE, testCase.geom->size());
		//EDIT2:
		//gpuSetup.enqueueWriteBuffer( gpuGeom, testCase.geom->getBuffer() );
		
		


		std::cout<<"Device Buffers allocation complete"<<std::endl;

		std::memcpy( ptrs.getHostPtr(), &(testCase.geom->ptrs[0]), size*2*sizeof(GEOMTYPE) );
		//EDIT 2:
	   //check_navigation(ptrs.getHostPtr(), size);
		ptrs.transferToDevice();
		gpuSetup.finish();

		//EDIT
		kernelTest.setArg(0, result.getDeviceBuffer());
		if( GLOBAL_MODE  == 1)
		   kernelTest.setArg(1, noStepArray);
		

		gpuSetup.enqueueKernel( kernelTest, 8, 8);
		gpuSetup.finish();
		result.transferFromDevice();
		gpuSetup.finish();
		FinalResult * final;
		
		int * a = ( int *)(ptrs.getHostPtr());
		std::cout<<"Printing input: \n";
		for (int i=0; i<8; i++)
			std::cout<< a[i]  << " ";
		//EDIT : Printing the output array
		std::cout<<"Printing output: \n";
		for (int i=0; i<8; i++)
			std::cout<< (( int *)(result.getHostPtr()))[i] << " ";
		//EDIT: Changed kernel Test to fix Prefix Sum
		//std::cout<< " Values that were returned: ";
		//std::cout<<" For thread 1: ";
		//final = (FinalResult *)result.getHostPtr();
		//final += sizeof( FinalResult);
		//std::cout<<": Min. Step value = "<< final->step<< " and safety returned = "<< final->safety << std::endl;


		// MODIFY: have to loop through this and add as a check

		/*
		for( int i =0 ; i < 4; i++)
		{	
			std::cout<<" For thread "<< i;
			final = (FinalResult **)result.getHostPtr();
			std::cout<<": Min. Step value = "<< final[i]->step<< " and safety returned = "<< final[i]->safety << std::endl;
			
		}
		*/
		/*std::cout<<"\nOriginal values: ";
		for( int i =0 ; i<32; i++)
		{
		std::cout << ((int *)ptrs.getHostPtr())[i]<<" ";
		}
		*/

		//NOTE: kernelCheck is badly named. It was originally the replacement for kernel getPtr and was used to return the geometry start location on the GPU.
		// IT can also be used to check the sizes on CPU and GPu are consistent

		kernelCheck.setArg( 0 , gpuGeom);
		kernelCheck.setArg(	1, result.getDeviceBuffer());
		gpuSetup.enqueueKernel(kernelCheck, 1, 1);
		gpuSetup.finish();
		//EDIT
		result.transferFromDevice();
		gpuSetup.finish();

		std::cout<< "On the CPU, size of GEOMTYPE = " << sizeof( GEOMTYPE )<<"\n";
		//REMOVE:
		std::cout<<"Size of GEOMTYPE on GPU -> "<<*((int  *)result.getHostPtr())<<std::endl;
		//MODIFY:
			// Assert that these are equal here.

		//EDIT2:
		//int answer = *(int*)(result.getHostPtr());
		//std::cout<<"Result before is "<< answer<< std::endl;
		
		std::cout<< "About to run relocate, no problem so far\n";

		//EDIT2
		/*
		gpuSetup.enqueueWriteBuffer( gpuGeom, testCase.geom->getBuffer() );
		gpuSetup.finish();
		kernelRelocate.setArg( 0, ptrs.getDeviceBuffer());
		kernelRelocate.setArg (1, gpuGeom);
		kernelRelocate.setArg(2, sizeof(int), &size);
		


		//MODIFY
		gpuSetup.enqueueKernel(kernelRelocate,size*2,2);
		gpuSetup.finish();
		*/
		//OLD
		// Fetch address of gpuGeom in device memory space (kludge)
		//kernelGetPtr.setArg( 0, gpuGeom );
		//kernelGetPtr.setArg( 1, gpuAux.getDeviceBuffer() );
		//gpuSetup.enqueueTask( kernelGetPtr );
		//gpuAux.transferFromDevice( CL_TRUE, 0, sizeof(cl_mem) );
		//cl_mem gpuhandle = *(cl_mem*)gpuAux.getHostPtr();
		
		//EDIT2
		cl_mem gpuhandle = *(cl_mem*)result.getHostPtr();
	
		
		
		
		//OLD:EDIT2
		testCase.geom->relocate( gpuhandle );

		//EDIT2
		gpuSetup.enqueueWriteBuffer( gpuGeom, testCase.geom->getBuffer() );

		//EDIT2: Kernel which returns the checks
		if(CHECK == 2 || CHECK == 4)
		{
			kernelCheckGeometry.setArg( 0 , gpuGeom);
			kernelCheckGeometry.setArg( 1 , result.getDeviceBuffer());
			kernelCheckGeometry.setArg( 2, sizeof(int), &number_of_increments);
			gpuSetup.enqueueKernel(kernelCheckGeometry, 1, 1);
			gpuSetup.finish();
	
			result.transferFromDevice();
			gpuSetup.finish();

			check_navigation( result.getHostPtr(), testCase.geom->VolumeStore, number_of_increments);
						// print for Geometry test;
		}

		

		// MODIFY: He also notes how this is not the optimal way of doing this. Memory can be saved here.
		std::memcpy( gpuInput.getHostPtr(), &(testCase.input[0]), gpuInput.size() );
		//REMOVE
		std::cout<<"About to run trace\n";
		
		// THis part was written for test purpose to see if error was caused due to shared memory

	
		// This kernel is getting to have WAY too many arguments
		gpuInput.transferToDevice();
		gpuSetup.finish();
		// Set GPU kernel arguments
		kernelTrace.setArg( 0, gpuInput.getDeviceBuffer() );
		kernelTrace.setArg( 1, gpuOutput.getDeviceBuffer() );
		kernelTrace.setArg( 2, gpuGeom );
		kernelTrace.setArg( 3, sizeof(G4double), &(testCase.phys_step) );
		kernelTrace.setArg( 4, sizeof(cl_int), &numInput );	
#ifdef CHECK
		kernelTrace.setArg( 5, result.getDeviceBuffer());
#endif
		// Two uses -: One for debugging and one for checking. Remove at some point.
		if( CHECK == 1 || CHECK == 4)
			kernelTrace.setArg( 5, result.getDeviceBuffer());
#if( GLOBAL_MODE ==1)
		//kernelTrace.setArg( 6, Numbers_Of_Solid );
		//kernelTrace.setArg( 7, Sum_Of_Solid );
		kernelTrace.setArg( 6, Solids );
		kernelTrace.setArg( 7, Result_For_Current_Solid );
		kernelTrace.setArg( 8, Compacter_Result );		
		kernelTrace.setArg( 9, nullVNode);
#endif

/*
NOTE: The current way of setting kernel trace's arguments is bad and is bound to cause problems in future
Change the implementation to either removes one check or perhaps just replace the existing check with something
more useful.
*/
		
		std::cout<<"Arguments set and value of Physical step sent on CPU is = "<< (testCase.phys_step)<<std::endl;
		// Write input to GPU memory
		// const my_clock_t t1 = my_clock();
		//OLD
		/*gpuSetup.enqueueWriteBuffer( gpuGeom, testCase.geom->getBuffer() );*/
		
		
		std::cout<< "Write complete, transfer done, finish\n";
		

		// Actual execution
		//const my_clock_t t2 = my_clock();
		//EDIT
		//gpuSetup.enqueueKernel( kernelTrace, numInput, blockSize );
		gpuSetup.enqueueKernel( kernelTrace, BlockSize, BlockSize );

		
		gpuSetup.finish();
		std::cout<<"Kernel trace done\n";
		

		//EDIT2:
		result.transferFromDevice();
		gpuSetup.finish();
		
		if( CHECK==1 || CHECK==4)
			check_navigation( result.getHostPtr(), testCase.geom->VolumeStore, 0);
					// Run distance check

		// Transfer results back to host memory
		//const my_clock_t t3 = my_clock();
		gpuOutput.transferFromDevice();
		gpuSetup.finish();


		//for (int i=0; i<10; i++)
			//std::cout<< (( int *)(result.getHostPtr()))[i] << " ";
		
			std::cout<< "\n";
		
		for (int i=0; i<32; i++)
			std::cout<< (( G4double *)(result.getHostPtr()))[i] << " ";

		std::cout<<"From CPU -> the first particles position and direction are -" <<" X-: "<<testCase.input[0].pos.x<<" Y-: "<<testCase.input[0].pos.y<<" Z-: "<<testCase.input[0].pos.z;
		std::cout<<" \n and directions are -:  X-: "<<testCase.input[0].dir.x<<" Y-: "<<testCase.input[0].dir.y<<" Z-: "<<testCase.input[0].dir.z;
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
