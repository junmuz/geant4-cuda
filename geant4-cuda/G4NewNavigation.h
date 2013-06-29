
/**
 * G4NewNavigation header
 * More details can be found in G4NewNavigation.c
 * based on G4NewNavigation.hh of Geant 4.9.3
 */

#ifndef G4NEWNAVIGATION_H
#define G4NEWNAVIGATION_H

#include "G4Voxels.h"
#include "G4NavigationHistory.h"

#ifdef CPU_CODE
#define USE_BLIST
#include <string.h>
#endif

#define K_MAX_VOXEL_STACK_DEPTH 4 // mind alignment

typedef struct //__attribute__((__aligned__(16)))
{
	G4double fVoxelSliceWidthStack[K_MAX_VOXEL_STACK_DEPTH]; 
	// Width of voxels at each level 

	GEOMETRYLOC G4SmartVoxelHeader* fVoxelHeaderStack[K_MAX_VOXEL_STACK_DEPTH];
	// Voxel headers at each level

	G4int fVoxelNodeNoStack[K_MAX_VOXEL_STACK_DEPTH];    
	  // Node no point is inside at each level 

	G4int fVoxelNoSlicesStack[K_MAX_VOXEL_STACK_DEPTH];
	  // No slices per voxel at each level
	  
	EAxis fVoxelAxisStack[K_MAX_VOXEL_STACK_DEPTH];
	  // Voxel axes
	  
	G4int fVoxelDepth;
	// Note: fVoxelDepth==0+ => fVoxelAxisStack(0+) contains axes of voxel
	//       fVoxelDepth==-1 -> not in voxel

	GEOMETRYLOC G4SmartVoxelNode *fVoxelNode;
	
#ifdef USE_BLIST
	char *fBlist;
	int fBlistSz;
#endif
}
G4NewNavigation;

#ifndef HOST_CODE

MAYINLINE void G4NewNavigation_ctor( G4NewNavigation *This );

MAYINLINE G4bool G4NewNavigation_LevelLocate(
	G4NewNavigation *This,
	G4NavigationHistory *history,
	GEOMETRYLOC const G4VPhysicalVolume *blockedVol,
	G4ThreeVector globalPoint,
	const G4ThreeVector* globalDirection,
	const G4bool pLocatedOnEdge, 
	G4ThreeVector *localPoint );
	
MAYINLINE GEOMETRYLOC G4SmartVoxelNode* G4NewNavigation_VoxelLocate(
	G4NewNavigation *This,
	GEOMETRYLOC G4SmartVoxelHeader *voxelHeader,
	G4ThreeVector point);
	
MAYINLINE
G4double
G4NewNavigation_ComputeStep(
			G4NewNavigation *This,
			G4ThreeVector localPoint,
			G4ThreeVector localDirection,
			const G4double currentProposedStepLength,
			G4double *newSafety,
			G4NavigationHistory *history,
			G4bool *validExitNormal,
			G4ThreeVector *exitNormal,
			G4bool *exiting,
			G4bool *entering,
			GEOMETRYLOC G4VPhysicalVolume *(*pBlockedPhysical)
			,SHAREDMEM int * Numbers_Of_Solid,
		
			 SHAREDMEM int * Sum_Of_Solids,

		     SHAREDTYPE SolidInfo  * Solids,
		
			 SHAREDTYPE ResultInfo * Result_For_Current_Solid,
		
		     SHAREDTYPE FinalResult * Compacter_Result,

			 SHAREDMEM bool * noStepArray,

			 SHAREDMEM PointInformation * LocationArray,
			 GEOMETRYLOC G4SmartVoxelNode * nullVNode,
			 G4bool cur_vol_local
			
#ifdef CHECK 
			,GEOMETRYLOC float * Result
#endif		
			);
	
MAYINLINE G4double G4NewNavigation_ComputeSafety(
	G4NewNavigation *This,
	G4ThreeVector localPoint,
	const G4NavigationHistory *history);

#endif

#endif
