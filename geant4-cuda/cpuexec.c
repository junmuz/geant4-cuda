
/** CPU version main loop. Compile as C99 */

//#define PRINTDOTS

#ifdef PRINTDOTS
#include <stdio.h>
#endif

#include "everything.h"
#include "stubParticle.h"
#include "G4Navigator.h"
#include "G4VPhysicalVolume.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

// defined in Makefile
//#define NTHREADS 2

void cpuexec(
	int problemSz,
	const Particle *particles,
	G4VPhysicalVolume *geometryRoot,
	G4double *output,
	G4double phys_step)
{
	#ifdef USE_OPENMP
	#pragma omp parallel for num_threads(NTHREADS)
	#endif
	for ( unsigned tid=0; tid<NTHREADS; ++tid )
	{
		G4Navigator nav;
		G4Navigator_ctor(&nav);
		
		G4Navigator_SetWorldVolume( &nav, geometryRoot );
	
		for ( int i=tid; i<problemSz; i += NTHREADS )
		{
			#ifdef PRINTDOTS
			fprintf(stderr, ".");
			#endif
			
			Particle p = particles[i];
				
			const G4VPhysicalVolume * cur_vol =
				G4Navigator_LocateGlobalPointAndSetup(
					&nav, p.pos, NULL, false, true );
				
			G4double step, safety = 0;
			G4double integratedDensity = 0;
			
			unsigned numSteps = 0;
			G4double dist = 0;
			
			while ( cur_vol )
			{
				numSteps++;
				
				const double curDensity = G4LogicalVolume_GetMaterial(
					G4VPhysicalVolume_GetLogicalVolume( cur_vol ))->property;
				
				step = G4Navigator_ComputeStep( &nav, p.pos, p.dir, phys_step, &safety );				
				if ( step == kInfinity )
				{
					step = phys_step;
				}
				
				const G4double nextStepIntegratedD = curDensity * step;
				
				#if PHYSICS
				if ( nextStepIntegratedD + integratedDensity > p.t )
				{
					const G4double left = p.t - integratedDensity;
					const G4double lastStep = left / curDensity;
					dist += lastStep;
					break;
				}
				#endif
				
				G4Navigator_SetGeometricallyLimitedStep( &nav );
				
				integratedDensity += nextStepIntegratedD;
				dist += step;
				
				p.pos = G4ThreeVector_saxpy( step, p.dir, p.pos );
				
				cur_vol =
					G4Navigator_LocateGlobalPointAndSetup(
						&nav, p.pos, &(p.dir), true, false );
			}
			
			#if PHYSICS
			output[i] = dist;
			#else
			output[i] = integratedDensity;
			#endif
		}
	}
}
