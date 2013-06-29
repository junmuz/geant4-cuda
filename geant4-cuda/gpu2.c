
/** GPU (practically CUDA) kernels for the step parallel GPU program */

#include "stubParticle.h"
#include "G4Navigator2.h"

// defined in Makefile
//#define ITERATIONS_BEFORE_COMPACT 20

/** Initialization of Navigators, input & output */
GLOBALFUNC void init(
	GLOBALTYPE Particle *particles,
	GLOBALTYPE G4double *output,
	GLOBALTYPE G4Navigator *navs,
	//GLOBALTYPE int *keys,
	//GLOBALTYPE int *vals,
	GLOBALTYPE int *particleIndices,
	GEOMETRYLOC G4VPhysicalVolume *worldVolumeAndGeomBuffer,
	int offset,
	GLOBALTYPE const int *totalSize )
{
	const unsigned globalIdx = get_global_id(0) + offset;
	if (globalIdx >= *totalSize) return;
	
	const int pIdx = particleIndices[globalIdx];
	
	// particleIndex < 0 means an uninitialized particle
	if ( pIdx < 0 )
	{
	  // -1 means: leave slot unused
	  if ( pIdx < -1 )
	  {
		// -2 => 0, -3 => 1, -4 = 2, ...
		particleIndices[globalIdx] = -pIdx - 2; 
		
		Particle p = particles[globalIdx];
		
		G4Navigator nav;// = navs+globalIdx;
		G4Navigator_ctor(&nav);
		nav.particle = p;
		G4Navigator_SetWorldVolume( &nav, worldVolumeAndGeomBuffer );

		GEOMETRYLOC const G4VPhysicalVolume * cur_vol =
			G4Navigator_LocateGlobalPointAndSetup(
				&nav, p.pos, NULL, false, true );
				
		output[globalIdx] = 0;
		navs[globalIdx] = nav;
	  }
	}
}

/** Step calculation. Does ITERATIONS_BEFORE_COMPACT steps. */
GLOBALFUNC void doStep(
	GLOBALTYPE G4double *output,
	GLOBALTYPE G4Navigator *navs,
	GLOBALTYPE int *keys,
	GLOBALTYPE const int *vals,
	GLOBALTYPE int *particleIndices,
	G4double phys_step,
	int offset,
	GLOBALTYPE const int *totalSize )
{
	const unsigned globalIdx = get_global_id(0) + offset;
	if (globalIdx >= *totalSize) return;
	
	const int origIdx = vals[globalIdx];
	if (particleIndices[origIdx] < 0)
	{
		keys[globalIdx] = 0;
		return; // < 0 means finished or unused slot
	}
	
	G4Navigator nav = navs[origIdx];
	Particle p = nav.particle;
	
	G4double localIntegratedDensity = 0.0;
	G4double step, safety = nav.fPreviousSafety;
	
	#ifdef PHYSICS
	G4double dist = output[origIdx];
	#endif
	
	GLOBALTYPE const G4VPhysicalVolume * cur_vol = 
		G4NavigationHistory_GetTopVolume(&(nav.fHistory));
		
	int newkey = 1;
		
	if ( cur_vol != GEOMETRYNULL )
	{
	  for (int i=0; i<ITERATIONS_BEFORE_COMPACT; ++i)
	  {
		const G4double curDensity = G4LogicalVolume_GetMaterial(
				G4VPhysicalVolume_GetLogicalVolume( cur_vol ))->property;
				
		step = G4Navigator_ComputeStep( &nav, p.pos, p.dir, phys_step, &safety );				
		if ( step == kInfinity ) step = phys_step;
				
		const G4double nextStepIntegratedD = curDensity * step;
	
		#ifdef PHYSICS
		if ( localIntegratedDensity + nextStepIntegratedD > p.t )
		{
			const G4double left = p.t - localIntegratedDensity;
			const G4double lastStep = left / curDensity;
			dist += lastStep;
			newkey = 0;
			break;
		}
		
		dist += step;
		#endif

		localIntegratedDensity += nextStepIntegratedD;
		
		G4ThreeVector_sum_assign( &(p.pos), G4ThreeVector_mult( p.dir, step ) );
		G4Navigator_SetGeometricallyLimitedStep( &nav );
		
		cur_vol = G4Navigator_LocateGlobalPointAndSetup(
				&nav, p.pos, &(p.dir), true, false );
		
		if ( cur_vol == GEOMETRYNULL )
		{
			newkey = 0;
			break;
		}
	  }
	}
	else
		newkey = 0;

	#ifdef PHYSICS
	p.t -= localIntegratedDensity;
	output[origIdx] = dist;
	#else
	output[origIdx] += localIntegratedDensity;
	#endif
		
	keys[globalIdx] = newkey;

	if (newkey != 0)
	{
		nav.particle = p;
		navs[origIdx] = nav;
	}
	else
	{
		if (particleIndices[origIdx] >= 0)
		{
			// Finished
			// 0 => -2, 1 => -3, 2 => -4, ...
			particleIndices[origIdx] = -particleIndices[origIdx] - 2;
		}
	}
}

#ifdef OPENCL_CODE
GLOBALFUNC void getPtr( GLOBALTYPE void *ptr, GLOBALTYPE void * GLOBALTYPE *out )
{
	*out = ptr;
}
#endif
