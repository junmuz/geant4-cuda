// Defining all parameters from make file here. Fix for Visual Studio Windows.

//REMOVE:
// Defining a new macro only for test to check if the code crash is because of the use of shared memory.
// Im the final version used for profiling, set GLOBAL_MODE  to 0.

#define GLOBAL_MODE  1

#define INLINE_EVERYTHING
#define NTHREADS 12
#define	ITERATIONS_BEFORE_COMPACT 100
#define	K_NAVIGATION_HISTORY_DEPTH 16
#define VERBOSE
#define ONLY_BOX_AND_ORB

// Have to check for 32 or 64 bit GPU before this.
// 32-bit
#define GEOMTYPE unsigned int
// 64 bit
// #define GEOMTYPE unsigned long;

#define ENABLE_VOXEL_NAVIGATION
//#define NEW_NAVIGATION

// NOTE: To define NEW_NAVIGATION you MUST define ENABLE_VOXEL_NAVIGATION as well. NEW_NAVIGATION is a rewrite but at some level it is still a type 
// of Voxel navigation

//EDIT
#define CHECK 0
	// CHECK can help confirm that the kernel trace is executing correctly, and if the geometry relocation is correct.
	/*
		TYPE of check
		1-> Distance check.
		2-> Geometry check
		3-> Position check
		4-> Run all checks.
		0-> Do not run check.
	*/

//EDIT:2
// #define PHYSICS 0