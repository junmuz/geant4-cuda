
/** Track parallel GPU main loop */
#define OPENCL_CODE 1




#include "gpuconf.h"

#include "stubParticle.h"
/* NOTE: 
	Prefix_Sum is defined here so that it can be accessed by the functions in files from the #include below
*/

//------------------------------------
// Prefix sum function for the parallel calculation of sums 
// on an array in shared memory.
#ifdef CUDA
MAYINLINE void Prefix_Sum ( int * input, int * output, int length)
#else
MAYINLINE void Prefix_Sum ( SHAREDMEM int * input, SHAREDMEM int * output, int length)
#endif
{
	/* ------------------------------
		Prefix_Sum implementation for the parallel calculation of the sum of elements of an array. The input array is summed and the results are stored i output.
		Each element of the final output array consists of the sum upto and not including that element. Implementation is based on the one provided in the CUDA SDK,
		AMD APP SDK and Chapter 39 of GPU Computing Gems 3.

		This is NOT original code.
		----------------------------*/
#if (GLOBAL_MODE  == 1)
	int tid = get_global_id(0);
#else
	int tid = get_local_id(0);
#endif
	//int length = 32;
	int offset = 1;

	if ( tid< length)
		output[tid] = input[ tid ];

	
    /* build the sum in place up the tree */
	for(int d = length>>1; d > 0; d >>=1)
	{
BARRIER_FLEXIBLE;
		if(tid<d)
		{
			int ai = offset*(2*tid + 1) - 1;
			int bi = offset*(2*tid + 2) - 1;
			
			output[bi] += output[ai];
		}
		offset *= 2;
	}

    /* scan back down the tree */

    /* clear the last element */
	if(tid == 0)
	{	// Tweak to get the total sum stored as the last element plus 1	
		// EDIT
		//block[length] = block[length -1];
		output[length - 1] = 0;
		
	}

    /* traverse down the tree building the scan in the place */
	for(int d = 1; d < length ; d *= 2)
	{
		offset >>=1;
		BARRIER_FLEXIBLE;
		
		if(tid < d)
		{
			int ai = offset*(2*tid + 1) - 1;
			int bi = offset*(2*tid + 2) - 1;
			
			float t = output[ai];
			output[ai] = output[bi];
			output[bi] += t;
		}
	} 
BARRIER_FLEXIBLE;


}
#ifdef CUDA
MAYINLINE
G4bool NoStepReduction( G4bool * noStepArray, int length )
#else
MAYINLINE
G4bool NoStepReduction( SHAREDMEM G4bool * noStepArray, int length )
#endif
{
	// Reduction of No Step values from across threads in a block to return the boolean Or of all noStep values present in the noStepArray
	#if (GLOBAL_MODE  == 1)
	int tid = get_global_id(0);
#else
	int tid = get_local_id(0);
#endif
	//int length = 32;
	int offset = 1;


	
    /* build the sum in place up the tree */
	for(int d = length>>1; d > 0; d >>=1)
	{
     BARRIER_FLEXIBLE;
		if(tid<d)
		{
			int ai = offset*(2*tid + 1) - 1;
			int bi = offset*(2*tid + 2) - 1;
			
			noStepArray[bi] = (noStepArray[ai] || noStepArray[bi]);
		}
		offset *= 2;
	}
	G4bool result = noStepArray[ length - 1 ];
	BARRIER_FLEXIBLE;
	return result;
}
// NOTE: The function Find_Minimum can be found in G4Navigator.h.
// Its location may seem odd but it is there so that it can be accessed by the ComputeStep Method and 
// the structs it uses are defined.


#include "G4Navigator.h"




GLOBALFUNC void trace(
	GLOBALTYPE Particle *input,
	GLOBALTYPE G4double *output,
	GEOMETRYLOC G4VPhysicalVolume *worldVolumeAndGeomBuffer,
	G4double phys_step,
	int totalSize
	//REMOVE
#if (CHECK == 1) || (CHECK == 4)
	,GLOBALTYPE float * result
#endif
#ifdef CHECK
	,GLOBALTYPE G4double * Result
#endif

#if ( GLOBAL_MODE ==1)
	//,GLOBALTYPE int * Numbers_Of_Solid
	//,GLOBALTYPE int * Sum_Of_Solids
	,GLOBALTYPE SolidInfo * Solids
	,GLOBALTYPE ResultInfo * Result_For_Current_Solid
	,GLOBALTYPE FinalResult * Compacter_Result,
	GEOMETRYLOC G4SmartVoxelNode * nullVNode
	//,GLOBALTYPE bool * noStepArray
#endif
	)
{
	
	const unsigned globalIdx = get_global_id(0);
	const unsigned localIdx = get_local_id(0);
#if(GLOBAL_MODE  == 1)
	const unsigned locationId = globalIdx;
#else
	const unsigned locationId = localIdx;
#endif
	//const unsigned work_size = get_local_size(0);
	/* The next few lines of code have been added for testing purposes.
	The idea is that we only simulate 3 threads and store information regarding each of them
	in the result array. ( Not to be confused with Result, which is here only for testing the
	new navigation part. I do not plan to continue with this naming convention. They will be renamed
	when the new navigation is confirmed to be up and running.) The imformation stored is very specific
	and the way it is stored and printed is specific. The test itself is not essential and most likely will
	be removed in future versions.
	*/

#if (CHECK == 1) || (CHECK == 4) 
	if (globalIdx >= 3 ) return;
		int i=0;
#else	
	if (globalIdx >= totalSize ) return;
#endif

	// The shared memory arrays are defined here. Their pointers are then passed from here to the navigator ( see G4Navigator_ComputeStep)
	// At least with OpenCL 1.1 there seems to be an issue when I define these arrays directly into the voxelnavigation file.

	// MODIFY: Change the sizes of these arrays by computing how much you would need.
#ifdef NEW_NAVIGATION 
#if (GLOBAL_MODE !=1)
	SHAREDTYPE int Numbers_Of_Solid[ BlockSize ];
		// This array is required to store the number of a particular solid type in the voxel per thread.
		/* 
		For instance if thread 1 sees two boxes in its voxel and thread two has 5, the first element of the array would be 2
		and the second 5
		*/
	SHAREDTYPE int  Sum_Of_Solids[ BlockSize  ];

    SHAREDTYPE SolidInfo Solids[ BlockSize * Multiplier ]; 
		// This array stores information of type SolidInfo which contains a pointer to the Solid as well as the trackId.
		//
	SHAREDTYPE  ResultInfo Result_For_Current_Solid[ BlockSize * Multiplier ];
		// On each iteration the results of the distance to in method of that solid is stored in this array by all threads along with threadId

	SHAREDTYPE  FinalResult Compacter_Result[ BlockSize ];
		// This array will contain the minimums for every iteration for every track. Location 0 would have the minimum value for track 0, and so on.
	
	SHAREDTYPE bool  noStepArray [ BlockSize ];
	    // An array used to store the local values of noStep. See NewNavigation.c.
	// Compacter_Result is initialized to kInfinity. See the implementation of Find_Minimum.

	if ( localIdx < 32 ) // Hard-coded for now.
	{
		Compacter_Result[ localIdx ].step = kInfinity; 
		Compacter_Result [ localIdx ].safety = kInfinity;
	}
#else
	int number_of_threads = get_global_size(0);
	    // Declared here to be in scope for NewNavigation to access.
	Compacter_Result [ globalIdx ].step = kInfinity;
	Compacter_Result [ globalIdx ].safety = kInfinity;
#endif
#endif
	SHAREDMEM int Numbers_Of_Solid[ BlockSize ];
	SHAREDMEM int  Sum_Of_Solids[ BlockSize  ];
	SHAREDMEM bool  noStepArray [ BlockSize ];
	SHAREDMEM PointInformation LocationArray[ BlockSize ];
	SHAREDMEM G4VPhysicalVolume * info[ BlockSize ];
	G4VoxelNode_ctor( nullVNode ,1 );
	//SHAREDTYPE int abc [ 32 ] ;
	#ifdef CUDA
	G4bool Cur_Vol_Store [ BlockSize ];
	#else
	SHAREDMEM  G4bool Cur_Vol_Store [ BlockSize ];
	#endif
		// A store of the local values of cur_vol variable to decide if threads should exit the while ( cur_vol ) loop.
	G4Navigator navi;
	G4Navigator *nav = &navi;
	
	G4Navigator_ctor(nav);
	G4Navigator_SetWorldVolume( nav, worldVolumeAndGeomBuffer );

	Particle p = input[globalIdx];

	if( globalIdx == 0)
	{
	//	Result[ 0 ] = p.pos.x;
	//	Result[ 1 ] = p.pos.y;
	//	Result[ 2 ] = p.pos.z;
	//	Result[ 3 ] = p.dir.x;
	//	Result[ 4 ] = p.dir.y;
	//	Result[ 5 ] = p.dir.z;
	}

    GEOMETRYLOC const G4VPhysicalVolume * cur_vol =
		G4Navigator_LocateGlobalPointAndSetup(
			nav, p.pos, NULL, false, true,  Result );
	
	G4bool cur_vol_local = true, cur_vol_all = true;
	G4double step, safety = 0.1;
	G4double integratedDensity = 0;
	
	#ifdef PHYSICS
	G4double dist = 0;
	//unsigned numSteps = 0;
	#endif
	
	int temp = 0;
	
	while ( cur_vol_all )
	{
		//numSteps++;
		//if( cur_vol_local )
		{
		const G4double curDensity =
		  G4LogicalVolume_GetMaterial( G4VPhysicalVolume_GetLogicalVolume( cur_vol ))->property;
		
		PointInformation NewPoint = { p.pos, p.dir };
		LocationArray[ locationId ]  = NewPoint;
		
		if( temp == 1)	
		{
			Result[ locationId ] = step;//return;
			//Result[ temp*3 +1] = p.pos.y;
			//Result[ temp*3 +2] = p.pos.z;
		}
		
		step = G4Navigator_ComputeStep( nav, p.pos, p.dir, phys_step, &safety
			#ifdef NEW_NAVIGATION
					  , Numbers_Of_Solid, Sum_Of_Solids, Solids, Result_For_Current_Solid,
					  Compacter_Result, noStepArray, LocationArray, nullVNode
			#endif
					  , cur_vol_local
			#ifdef CHECK
					  , Result
			#endif
					  );	
	//	if( temp == 1) return;
/*
NOTE: The current implementation of passing in these shared memory pointers if NEW_NAVIGATION is defined looks clumsy.
There isn't much harm in just declaring shared memory arrays, and then using them only in the case of NEW_NAVIGATION.
In this way the code would be MUCH cleaner. 
The choice of whether to use these pointers or not could be left to the G4Navigator. 
*/		

		if ( step == kInfinity ) step = phys_step;
		
		const G4double nextStepIntegratedD = curDensity * step;
		int locationId = get_global_id(0);
		//if( temp == 0)
			//Result [ locationId ] = step;
		
#if ( CHECK == 1 ) || (CHECK == 4)
		if(i==0)
		{
			result[10 + globalIdx*3] = p.pos.x;
			result[11 + globalIdx*3] = p.pos.y;
			result[12 + globalIdx*3] = p.pos.z;
			result[30 + globalIdx*3] = p.dir.x;
			result[31 + globalIdx*3] = p.dir.y;
			result[32 + globalIdx*3] = p.dir.z;
		}
		// As stated earlier the check is very specific. It stores information for only three particles
		// and displays the 
#endif
		
		#if PHYSICS
		if ( nextStepIntegratedD + integratedDensity > p.t )
		{
			const G4double left = p.t - integratedDensity;
			const G4double lastStep = left / curDensity;
			dist += lastStep;
			break;
		}
		
		dist += step;
		#endif
	
		integratedDensity += nextStepIntegratedD;
		G4ThreeVector_sum_assign( &(p.pos), G4ThreeVector_mult( p.dir, step ) );
		
		G4Navigator_SetGeometricallyLimitedStep( nav );	
		//Result[ locationId ] = get_local_id(0);
	    
		if( globalIdx == 0 ){
			
			//Result[ locationId ]  = step;
			
		}
		
		cur_vol =
			G4Navigator_LocateGlobalPointAndSetup(
				nav, p.pos, &(p.dir), true, false, Result );
		//if( cur_vol  && globalIdx == 0)
			//Result[ 0 ] = temp;
		
		if ( !cur_vol )
			cur_vol_local = false;
		 
#if (CHECK == 1) || (CHECK == 4)
		if(globalIdx==0)
			result[29] = i;

		if(i==0)
		{		
			result[globalIdx] = step;
			result[20 + globalIdx*3] = p.pos.x;
			result[21 + globalIdx*3] = p.pos.y;
			result[22 + globalIdx*3] = p.pos.z;
			
		}
		
		i=i+1;
#endif
		
		}
		
		Cur_Vol_Store[ locationId ] = cur_vol_local;
		BARRIER_FLEXIBLE;
		cur_vol_all = NoStepReduction( Cur_Vol_Store, BlockSize );
			// Calling the OR reduction ( called NOStepreduction) on the cur_vol_store to decide whether all threads are ready to exit the while loop.
		BARRIER_FLEXIBLE;
		temp++;
		
	//	if( temp>1)
	//	return;
}
		
	#if PHYSICS
	output[globalIdx] = dist;
	#else
	output[globalIdx] = integratedDensity;
	#endif
	
}

/**
 * A kernel for fetching the device memory space address of a pointer
 * in OpenCL. An ugly hack.
 */
//OLD
/*
GLOBALFUNC void getPtr( GLOBALTYPE void *ptr, GLOBALTYPE void * GLOBALTYPE *out )
{
	*out = ptr;
}*/
//EDIT:
GLOBALFUNC void relocate ( GLOBALTYPE int * ptr, GLOBALTYPE void * buf, int size )
{ 
  typedef GLOBALTYPE unsigned char byte;
  const unsigned globalidx = get_global_id(0);
  if(globalidx>=size) return;
  //GLOBALTYPE void ** newbegin = &buf;
  int destoffs, targoffs;
  destoffs = *(ptr + 2*globalidx);
  targoffs = *(ptr + 2*globalidx+ 1);
  *((byte*)buf+destoffs) = (byte)((byte*)buf + targoffs);

}

//EDIT: New kernel

GLOBALFUNC void check( GLOBALTYPE  G4VPhysicalVolume *worldVolumeAndGeomBuffer, GLOBALTYPE unsigned long * result)
{ 
 GEOMTYPE hope = ( GEOMTYPE )worldVolumeAndGeomBuffer;

 *result = hope;
 // *result = sizeof( GEOMTYPE );
	
}

GLOBALFUNC void test ( GLOBALTYPE bool * output 
#if(GLOBAL_MODE ==1) 
	,GLOBALTYPE bool * input 
#endif
	)
{
	// TO TEST Find_Min use __global FinalResult * output and uncomment some of the stuff below
	
	// This kernel is meant for testing purposes only. 
	// Any new funcitionality or algorithm could, in part, be tested here to ensure that it works.
#if( GLOBAL_MODE   != 1)
	SHAREDTYPE bool  input[8];
#endif
	int tid = get_global_id(0);
	int offset = 1;
	

	G4bool result;
	if( tid == 0)
	{
	input[ 0] = true;
	input[ 1] = true;
	input[ 2] = true;
	input[ 3] = true;
	input[ 4] = true;
	input[ 5] = true;
    input[ 6] = false; 
    input[ 7] = true;
	}
	
	//result = NoStepReduction( input, 8);
	BARRIER_ALL;
	//output [ tid] = result;
	/*
	SHAREDTYPE ResultInfo block[64 ];
	SHAREDTYPE int other[24];
	SHAREDTYPE int size[24];
	SHAREDTYPE FinalResult out[ 64 ] ;
	out[ tid ].step = kInfinity;
	barrier( CLK_LOCAL_MEM_FENCE);


	block[0].safety = 1.0;
	block[0].step = 0.5;
	block[0].trackId = 0;
	
	block[1].safety = 1.5;
	block[1].step = 3.0;
	block[1].trackId = 0;

	
	block[2].safety = 2.6;
	block[2].step = 0.3;
	block[2].trackId = 0;

	block[3].safety = 1.0;
	block[3].step = 2.0;
	block[3].trackId = 1;
	
	block[4].safety = 1.6;
	block[4].step = 2.0;
	block[4].trackId = 1;
	
	block[5].safety = 1.8;
	block[5].step = 2.0;
	block[5].trackId = 1;
	
	block[6].safety = 1.4;
	block[6].step = 2.0;
	block[6].trackId = 1;

	block[7].safety = 1.2;
	block[7].step = 1.0;
	block[7].trackId = 2;
	
	block[8].safety = 1.0;
	block[8].step = 24.0;
	block[8].trackId = 2;
	
	block[9].safety = 1.2;
	block[9].step = 20.0;
	block[9].trackId = 2;
	
	block[10].safety = 1.21;
	block[10].step = 21.0;
	block[10].trackId = 3;
	
	block[11].safety = 1.4;
	block[11].step = 22.0;
	block[11].trackId = 3;
	
	block[12].safety = 1.354;
	block[12].step = 264.0;
	block[12].trackId = 3;

	block[13].safety = 1.45;
	block[13].step = 21.0;
	block[13].trackId = 3;

	block[14].safety = 1.075;
	block[14].step = 12.0;
	block[14].trackId = 3;

	block[15].safety = 1.088;
	block[15].step = 22.0;
	block[15].trackId = 3;
	

	
	block[1] = { 1.5, 3.0, 0, void* };
	block[2] = { 2.0, 1.0, 0, void* };
	block[3] = { 1.0, 2.0, 1, void* };
	block[4] = { 1.6, 2.0, 1, void* };
	block[5] = { 1.8, 2.0, 1, void* };
	block[6] = { 1.4, 2.0, 1, void* };
	block[7] = { 1.2, 1.0, 2, void* };
	block[8] = { 1.0, 24.0, 2, void* };
	block[9] = { 1.2, 20.0, 2, void* };
	block[10] = { 1.21, 21.0, 3, void* };
	block[11] = { 1.4, 22.0, 3, void* };
	block[12] = { 1.354, 264.0, 3, void* };
	block[13] = { 1.45, 21.0, 3, void* };
	block[14] = { 1.075, 12.0, 3, void* };
	block[15] = { 1.088, 22.0, 3, void* };
	

	other[0] = 0;
	other[1] = 3;
	other[2] = 7;
	other[3] = 10;
	
	size[0] = 3;
	size[1]	= 4;
	size[2] = 3;
	size[3] = 6;
	int PrevSum = other[tid];
	int siz = size[ tid];

	Find_minimum( block, out, PrevSum, siz );
	*/
	/*
	typedef struct{
		float safety;
		float step;
		int trackId;
		GLOBALTYPE const G4VSolid * Solid; // PhysicalVolume ?
	  }ResultInfo;

   typedef struct{
		float safety;
		float step;
		GLOBALTYPE const G4VSolid * Solid; // PhysicalVolume ?
	   }
   FinalResult;
   */
	//block[tid]   = tid;



	//__local G4bool input[64];
	
    /*write the results back to global memory */
	//output[tid]     = input[ 7];
	//BARRIER_ALL;
	//if(tid == 1)
		//output[ 0 ] = out[ 3 ] ;
    


}
//EDIT2
GLOBALFUNC void checkgeom( GLOBALTYPE  G4VPhysicalVolume *worldVolumeAndGeomBuffer, GLOBALTYPE int * result, int number_of_increments)
{	
	// Kernel that stores information about Physical Volumes encountered along a diagonal path
	// The offsets of the logical volume addresses encountered along this path are stored. These are then compared with offsets on GPU to confirm relocation was successful.
	const unsigned globalid = get_global_id(0);
	if(globalid>=1) return;

	PRIVATEMEM int i=0;
	G4Navigator navi;
	G4Navigator *nav = &navi;
	G4Navigator_ctor(nav);
		// Using the navigator to navigate through the geometry
	G4Navigator_SetWorldVolume( nav, worldVolumeAndGeomBuffer );
	G4ThreeVector pos = G4ThreeVector_create(  0.0, 0.0, 0.0);

	GLOBALTYPE const G4VPhysicalVolume * cur_vol;// = G4Navigator_LocateGlobalPointAndSetup( nav, pos , ((void*)0), false, true );
		// Set up
	
	GEOMTYPE geom_start = ( GEOMTYPE )worldVolumeAndGeomBuffer;
		// Starting address of the geometry
	pos =  G4ThreeVector_create(  0.7, 1.0, 0.7);
		// Starting position.
	float x_increment = 0.2, y_increment = 0.2, z_increment = 0.2;
		// Amounts to increment in x, y and z directions.
	
	
	
	for( i=0; i < number_of_increments*3 ; i+=3)
	{
		//cur_vol = G4Navigator_LocateGlobalPointAndSetup( nav, pos , ((void*)0), true, false );
			// Get the current Physical Volume
		result[i] = ( int ) cur_vol->count;
			// Identifier to find Physical Volume
		//result[i + 1] = ((int)(cur_vol->flogical) - geom_start);  
			// Offset from gpu start
		result[i + 1] = (( GEOMTYPE )(cur_vol->flogical) - geom_start);  

		result[i + 2] =  ( int ) G4LogicalVolume_GetMaterial( G4VPhysicalVolume_GetLogicalVolume( cur_vol ))->property;
			// The material property (density). Important test as it uses many pointers.

		pos.x+=x_increment;
		pos.y+=y_increment;
		pos.z+=z_increment;
			// Increment and continue
	}
	
}
