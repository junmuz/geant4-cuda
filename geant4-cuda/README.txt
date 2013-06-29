
_______________________________________________________

	Geant4 Navigation for GPUs
	Otto Seiskari, August 2010 (updated Jul 2011)
	Dhruva T.B. , June 2012
  This product includes software developed by Members of
  the Geant4 Collaboration ( http://cern.ch/geant4 ).
________________________________________________________

----- General info -----

This software is an experimental port of Geant4's Navigation functionality to (NVIDIA) GPUs. It supports a slightly modified version of G4VoxelNavigation and four solids: G4Orb, G4Box, G4Tubs and G4Cons. The program can be used to
benchmark Geant4 navigation with dummy physics and data analysis.

The code has also been tested and confirmed to run on AMD 
(32-bit) GPU and Intel 64 CPU with OpenCL 1.1. 
The current implementation has two possible particle setups (controlled with the PHYSICS build variable). If PHYSICS=0, the program creates a bunch of particles with a common starting point, travelling towards different points in
a rectangular grid. The density of the materials pentrated by a particle are integrated along each particle track to create a volume raytraced "x-ray" image of the current geometry. If PHYSICS=1 the program creates N particles with a common starting point and "energy" (defined by a mass attenuation
coefficient), travelling to different random directions. Each line of the output contains the direction of a particle and the "distance to first interaction" 


----- Building -----

The program can be compiled as three different versions:

 * Track parallel CPU/OpenMP program, writes its output to "imgcpu.txt",
   built with (for example) the command:
	make cpumain NTHREADS=12 PHYSICS=1

 * Track parallel CUDA program:
	make cuda STEP_PARALLEL=0 ITERATIONS_BEFORE_COMPACT=50

 * Step parallel CUDA program:
	make cuda

In order to make for a specific platform ensure that the parameters specified in the Makefile for the install location of OpenCL headers and/or CUDA libraries correctly specify the install path on the machine to be tested on.
The PHYSICS flag can be set in all cases, the default is PHYSICS=0.

In the track parallel version, the full track of a particle is computed by a single thread. In the step parallel version, a thread computes a certain (small) number of steps, ITERATIONS_BEFORE_COMPACT, for a particle, after
which inactive particles are discarded by doing a stream compaction operation with CUDPP.

The CPU versions should be buildable on any (at least Linux) system with g++ and OpenMP. The CUDA version needs CUDA and CUDPP. Before building, check and adjust the configuration parameters in the Makefile.


----- Running -----

If PHYSICS=0 (example)
	./cpumain -test=cms -step=0.01 -vox1=10 -vox2=20 -vox3=30 -xres=800 -yres=600
If PHYSICS=1 (example)
	./cpumain -test=spheres -step=100 -vox1=10 -vox2=20 -vox3=30 -n=100000 -mean=0.0061

The first parameters are test type, physical step size, voxelization limits  and test size (xres*yres or n). If PHYSICS=1 there are also two additional (optional) parameters, 'mean' (the mass attenuation coefficient in m^2/kg)and 'rounds'. (see below for more info)

There are currently four available tests:

toy1:
	A very simple test geometry with a sphere in a box.

toy2:
	Test geometry with G4Cons/G4Orbs placed in a regular grid in a box.

spheres:
	Different solids scattered on spherical shells. Based on a voxelization benchmark by Ivaylo Boyadzhiev. Not available in the OpenCL version.

cms:
	A test based on the geometry of the CMS detector defined in benchmarks/cms/cms.gdml. The file is stripped down to supported types	of geometry with parse.py that outputs cms.txt. The voxN command line	parameters (see above) need to be adjusted to enable voxelization in this test. Default voxelization settings yield a bad_alloc on 32-bit systems. Not available in the OpenCL version.

The "physical step scale" is the size of "physical step" (maximum step size)as a propotion of world volume diameter, set to e.g. 100 for unlimited steps.

The parameters vox1, vox2 and vox3 are the "minimum volumes to voxelize at level N" limits. These correspond to constants kMinVoxelVolumesLevelN in the original Geant4 code.

There is also an optional parameter "-rounds=N", that is a multiplier for the number of particles, n. On CPU systems, the two numbers are redundant. On the GPU, however, n is the maximum number of particles in GPU memory at one time,
that is, there are n "slots" for particles on the device. The slots are periodically (at least N times) refilled with new particles as the results from old, completely tracked particles have been read back.

The executable programs output the x-ray images (PHYSICS=0) or directions & distances to first interactions (PHYSICS=1) as text files:

	cpumain:	imgcpu.txt
	cudamain:	imgcuda.txt

The x-ray image files (PHYSICS=0) contain a tab/linefeed separated table of values that are readily viewed in MATLAB/Octave, e.g.:
	load imgcpu.txt; imagesc(imgcpu), colormap "gray"

The step parallel version of the program also prints the amount of remaining
active particles after each compact operation.

The "Elapsed" value written to STDERR is the real time elapsed while running the simulation and moving data between host and device memory. This number does NOT include the time it took to compile or load the geometry or GPU
programs.

No space is allowed in parameters: "-test = cms" will not work. Note also that especially the track parallel GPU programs (if the GPU is attached to a display) will crash / terminate if the execution takes longer than 5 seconds.

Problems with nvcc crashing / overusing memory are avoided crudely by disabling functionality. This is why all tests are not available with all versions of the program.

---- Debugging ----
If the code does not run, it is advised to set a value for the macro CHECK defined in the configuration file gpuconf.h
Two tests are currently supported, with CHECK as 1 or 2.
CHECK 2 is the geometry check which checks whether the relocation ( transferring the entire geometry from CPU to GPU and updating pointers to point to GPU memory locations ) worked correctly. It is highly recommended to add this check while running for the first time to ensure correct operation.
CHECK 1 is the distance check which checks if the positions of particles after each step is being updated correctly. It is not essential to run this check and in the current version of the code this check might not run because of redefinition of Result.


----- Under the hood -----


-- Creating and storing the geometry --

The geometry is created by the host program and stored on a single continuous block of (virtual) memory. After this a memory block of same size is allocated from device memory. Then the device memory space address of this memory block
is fetched (requires very dirty hacks in OpenCL) and all the pointers in the geometry array in host memory are "relocated" to point to the corresponding places in device memory space, and the geometry buffer can finally be copied to device memory.


-- Transformed source: C++ to C --

The source code is C99 with special keywords (preprocessor macros) that allow the program to be compiled as C99, C++, CUDA or OpenCL. The source was created by transforming original Geant4 C++ source code to C using the following
renaming conventions:

	G4XXX::YYY(A a, B b)	==>	G4XXX_YYY(G4XXX *This, A a, B b)
	G4XXX::G4XXX()		==>	G4XXX_ctor(G4XXX *This)
	G4XXX::~G4XXX()		==>	G4XXX_dtor(G4XXX *This)

The files are renamed as follows:

	G4XXX.hh	==>	G4XXX.h
	G4XXX.icc	==>	G4XXX_inline.c
	G4XXX.cpp	==>	G4XXX.c

In addition, code is heavily moved between .icc and .cpp file and in the current version all code is actually inlined to the either the "host main file" or the "main device code" file.
The current convention will at some point change to specify GPU versions. This is to account for future usage of the code with single threaded Geant4 where funcitonality can be implemented for navigation of existing tests to run on the GPU.

Additionally:

 - Virtual function calls (in G4VSolid) are replaced by switch functions
 - G4VPhysicalVolume now equals G4PVPlacement
 - Replicated volumes are not supported
 - Regular navigation is not supported
 - Materials and particles are replaced with dummy classes
 - Blocked volume list in voxel navigation is removed
 - std::vectors have been replaced by flat arrays (allocated on host side)
 - Especially G4NavigationHistory::fNavHistory is always of constant
   length: K_NAVIGATION_HISTORY_DEPTH (set in Makefile)
 - New class: G4CombinedNavigation (see below)


-- Transformed source: C to CUDA or OpenCL --

C99 code can transformed to CUDA or OpenCL by adding special keywords such as address space qualifiers (e.g. OpenCL's __local, CUDA's __shared__).
These are implemented as preprocessor macros defined in everything.h. They are generally ignored in CPU builds.

Warning: there are some unnecessary / redundant macros. See comments.


-- Main files --

There are three different "main host files" written in C++: cpumain.cpp,cudamain.cpp and openclmain.cpp. They contain the main functions and initialization code for different versions (CPU, CUDA, OpenCL) of the program.
The OpenCL version uses my C++-wrapper found in OpenCL/common/cl.hpp.
Functionality that is common to all host programs is implemented/included in hostcommons.hpp.

The files cpuexec.c and gpu.c (and cuda.cpp) contain the track parallel main loops of the CPU and GPU versions of the program respectively. The step parallel version of the program is implemented in cuda.cpp and gpu2.c. It also uses a modified version of the G4Navigator (tracked particle added to the navigator data structure) defined in G4Navigator2.h. The file cuda.cpp also includes a function for calculating suitable CUDA block and grid sizes for different quantities of remaining active particles instead of using fixed
block and grid sizes.

I have also added a new navigator called G4CombinedNavigation which combines voxel and normal navigation. This aims to be more efficient on the SIMD architecture and reduces the amount of GPU code.

The allocation of geometry objects and relocation of pointers is done through the BasicGeometry class methods defined in geometry_common.hpp. All different benchmark geometries are created/imported in classes derived from BasicGeometry (and located in the benchmarks folder). Geometry (in
geometry.hpp) is a minimalistic interface to geometry creation and pointer relocation.


-- Known issues( August 2012) --
- New Navigation does not work on CPU ( Intel 64 bit). On AMD GPUs, code compiles and runs but output is not correct.

- GLOBAL_MODE has to be defined and used when running on GPU. For some reason the code crashes when using shared memory ( even though amount of memory used is within the maximum limit)

- CMS and Spheres examples only tested to work on Nvidia CUDA platform. Memory required by these examples is quite large and the AMD platform does not have that much global memory.



