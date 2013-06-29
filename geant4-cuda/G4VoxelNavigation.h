
/**
 * G4VoxelNavigation header
 * based on G4VoxelNavigation.hh of Geant 4.9.3
 */

#ifndef G4VOXELNAVIGATION_H
#define G4VOXELNAVIGATION_H

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
G4VoxelNavigation;

#ifndef HOST_CODE

MAYINLINE void G4VoxelNavigation_ctor( G4VoxelNavigation *This );

MAYINLINE G4bool G4VoxelNavigation_LevelLocate(
	G4VoxelNavigation *This,
	G4NavigationHistory *history,
	GEOMETRYLOC const G4VPhysicalVolume *blockedVol,
	G4ThreeVector globalPoint,
	const G4ThreeVector* globalDirection,
	const G4bool pLocatedOnEdge, 
	G4ThreeVector *localPoint );
	
MAYINLINE GEOMETRYLOC G4SmartVoxelNode* G4VoxelNavigation_VoxelLocate(
	G4VoxelNavigation *This,
	GEOMETRYLOC G4SmartVoxelHeader *voxelHeader,
	G4ThreeVector point);
	
MAYINLINE
G4double
G4VoxelNavigation_ComputeStep(
			G4VoxelNavigation *This,
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
#ifdef CHECK 
			,GEOMETRYLOC G4double * Result
#endif		
// Perhaps this should be if CHECK == 3 or 4	
			);
	
MAYINLINE G4double G4VoxelNavigation_ComputeSafety(
	G4VoxelNavigation *This,
	G4ThreeVector localPoint,
	const G4NavigationHistory *history);

#endif

#endif
