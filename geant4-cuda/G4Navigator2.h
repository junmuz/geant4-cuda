
/**
 * Modified navigator (tracked particle added) header
 * based on G4Navigator.hh of Geant 4.9.3
 */

#ifndef G4NAVIGATOR_H
#define G4NAVIGATOR_H

#include "everything.h"

#include "G4NavigationHistory.h"

#ifdef ENABLE_VOXEL_NAVIGATION
#include "G4VoxelNavigation.h"
#endif

typedef struct
{
	// Helpers/Utility classes (?-alignment)
	
	// The current particle
	Particle particle;

	G4NavigationHistory fHistory;
	// Transformation and history of the current path
	// through the geometrical hierarchy.
	
#ifdef ENABLE_VOXEL_NAVIGATION
	G4VoxelNavigation fVoxelNav;
#endif
	
	G4ThreeVector fStepEndPoint;
	//  Endpoint of last ComputeStep 
	//  - can be used for optimisation (eg when computing safety)
	
	G4ThreeVector fLastLocatedPointLocal;
	// Position of the last located point relative to its containing volume.

	G4ThreeVector fExitNormal;  // Leaving volume normal, in the
							  // volume containing the exited
							  // volume's coordinate system
	G4ThreeVector fGrandMotherExitNormal;  // Leaving volume normal, in its 
										 // own coordinate system
								
	//G4ThreeVector  fPreviousSftOrigin;
	
	// --------- 4-byte alignment
	
	G4bool fEnteredDaughter;
	// A memory of whether in This Step a daughter volume is entered 
	// (set in Compute & Locate).
	//  After Compute: it expects to enter a daughter
	//  After Locate:  it has entered a daughter

	G4bool fExitedMother;
	// A similar memory whether the Step exited current "mother" volume
	// completely, not entering daughter.

	G4bool fWasLimitedByGeometry;
	// Set true if last Step was limited by geometry.

	G4bool fEntering;
	G4bool fExiting;
	// Entering/Exiting volumes blocking/setup
	// o If exiting
	//      volume ptr & replica number (set & used by Locate..())
	//      used for blocking on redescent of geometry
	// o If entering
	//      volume ptr & replica number (set by ComputeStep(),used by
	//      Locate..()) of volume for `automatic' entry
	
	// Count zero steps - as one or two can occur due to changing momentum at
	//                    a boundary or at an edge common between volumes
	//                  - several are likely a problem in the geometry
	//                    description or in the navigation
	//
	G4bool fLastStepWasZero;
	// Whether the last ComputeStep moved Zero. Used to check for edges.

	G4bool fLocatedOnEdge;       
	// Whether the Navigator has detected an edge

	G4bool fLocatedOutsideWorld;
	// Whether the last call to Locate methods left the world
	// G4PhysicalVolume* fLastVolumeLocated; 
	
	G4bool fValidExitNormal;    // Set true if have leaving volume normal
	
		// Utility information
	//
	//G4bool fCheck;
	// Check-mode flag  [if true, more strict checks are performed].
	G4bool fPushed;
	// Push flag  [if true, means a stuck particle has been pushed].
	
	G4int fNumberZeroSteps;
	// Number of preceding moves that were Zero. Reset to 0 after finite step
	
	#ifndef DOUBLE_PRECISION
	int align1;
	#endif
	
	// --------- 8-byte alignment
	
	G4double       fPreviousSafety; 
	// Memory of last safety origin & value. Used in ComputeStep to ensure
	// that origin of current Step is in the same volume as the point of the
	// last relocation
	
	GEOMETRYLOC G4VPhysicalVolume *fBlockedPhysicalVolume;
	
	GEOMETRYLOC G4VPhysicalVolume  *fTopPhysical;
	// A link to the topmost physical volume in the detector.
	// Must be positioned at the origin and unrotated.

}
G4Navigator;

#ifndef HOST_CODE

// FWDEC
MAYINLINE void G4Navigator_ctor( G4Navigator *This );
	
MAYINLINE void G4Navigator_SetWorldVolume(
	G4Navigator *This,
	GEOMETRYLOC G4VPhysicalVolume* pWorld );
	
MAYINLINE GEOMETRYLOC G4VPhysicalVolume* G4Navigator_LocateGlobalPointAndSetup(
		G4Navigator *This,
		G4ThreeVector globalPoint,
		const G4ThreeVector* pGlobalDirection,
		G4bool relativeSearch,
		G4bool ignoreDirection );
		
MAYINLINE G4double G4Navigator_ComputeStep(
		G4Navigator *This, 
		G4ThreeVector pGlobalpoint,
		G4ThreeVector pDirection,
		const G4double pCurrentProposedStepLength,
		G4double *pNewSafety);
		
MAYINLINE void G4Navigator_SetGeometricallyLimitedStep( G4Navigator *This );

MAYINLINE G4double G4NormalNavigation_ComputeStep(
	G4ThreeVector localPoint,
	G4ThreeVector localDirection,
	const G4double currentProposedStepLength,
	G4double *newSafety,
	G4NavigationHistory *history,
	G4bool *validExitNormal,
	G4ThreeVector *exitNormal,
	G4bool *exiting,
	G4bool *entering,
	GEOMETRYLOC G4VPhysicalVolume *(*pBlockedPhysical));
	
MAYINLINE G4double G4NormalNavigation_ComputeSafety(
	G4ThreeVector localPoint,
	const G4NavigationHistory *history );
	
MAYINLINE G4bool G4NormalNavigation_LevelLocate(
	G4NavigationHistory *history,
	GEOMETRYLOC const G4VPhysicalVolume *blockedVol,
	G4ThreeVector* globalPoint,
	const G4ThreeVector* globalDirection,
	G4bool pLocatedOnEdge, 
	G4ThreeVector* localPoint );
	
#ifdef ENABLE_VOXEL_NAVIGATION
MAYINLINE void G4VoxelNavigation_ctor( G4VoxelNavigation *This );
#endif

#ifdef ENABLE_COMBINED_NAVIGATION
MAYINLINE G4double G4CombinedNavigation_ComputeStep(
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
	GEOMETRYLOC G4VPhysicalVolume *(*pBlockedPhysical));
	
INLINEFUNC G4bool
G4CombinedNavigation_LevelLocate(
			G4VoxelNavigation *This,
			G4NavigationHistory* history,
			GEOMETRYLOC const G4VPhysicalVolume* blockedVol,
			G4ThreeVector globalPoint,
			const G4ThreeVector* globalDirection,
			const G4bool pLocatedOnEdge, 
			G4ThreeVector *localPoint );
			
MAYINLINE G4double
G4CombinedNavigation_ComputeSafety(
			G4VoxelNavigation *This,
			G4ThreeVector localPoint,
			const G4NavigationHistory *history );
#endif

#ifdef INLINE_EVERYTHING
#include "G4Navigator.c"
#include "G4NormalNavigation.c"
#ifdef ENABLE_VOXEL_NAVIGATION
#include "G4VoxelNavigation.c"
#endif
#ifdef ENABLE_COMBINED_NAVIGATION
#include "G4CombinedNavigation.c"
#endif
#endif

#endif // !defined( HOST_CODE )

#endif
