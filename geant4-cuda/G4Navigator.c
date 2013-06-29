
/**
 * G4Navigator implementation,
 * based on G4Navigator.(i)cc of Geant 4.9.3
 */

#include "G4Navigator.h"
#include "G4VPhysicalVolume.h"

#include "G4VSolid.h"

// ********************************************************************
// ResetState
//
// Resets stack and minimum of navigator state `machine'
// ********************************************************************
//
MAYINLINE void G4Navigator_ResetState( G4Navigator *This )
{
  This->fWasLimitedByGeometry  = false;
  This->fEntering              = false;
  This->fExiting               = false;
  This->fLocatedOnEdge         = false;
  This->fLastStepWasZero       = false;
  This->fEnteredDaughter       = false;
  This->fExitedMother          = false;
  This->fPushed                = false;

  This->fValidExitNormal       = false;
  This->fExitNormal            = G4ThreeVector_create(0,0,0);

  This->fPreviousSafety        = 0.0; 

  This->fNumberZeroSteps       = 0;
    
  This->fBlockedPhysicalVolume = GEOMETRYNULL;

  This->fLastLocatedPointLocal = G4ThreeVector_create( DBL_MAX, -DBL_MAX, 0.0 ); 
  This->fLocatedOutsideWorld   = false;
}

// class G4Navigator Inline implementation
//
// ********************************************************************

// ********************************************************************
// GetCurrentLocalCoordinate
//
// Returns the local coordinate of the current track
// ********************************************************************
//
/*INLINEFUNC
G4ThreeVector G4Navigator_GetCurrentLocalCoordinate( const G4Navigator *This )
{
	return This->fLastLocatedPointLocal;
}*/

// ********************************************************************
// ComputeLocalAxis
//
// Returns local direction of vector direction in world coord system
// ********************************************************************
//
INLINEFUNC
G4ThreeVector G4Navigator_ComputeLocalAxis( const G4Navigator *This, G4ThreeVector pVec)
{
	G4AffineTransform t =
		G4NavigationHistory_GetTopTransform( &(This->fHistory) );
	return G4AffineTransform_TransformAxis(&t, pVec);
}

// ********************************************************************
// ComputeLocalPoint
//
// Returns local coordinates of a point in the world coord system
// ********************************************************************
//
INLINEFUNC G4ThreeVector
G4Navigator_ComputeLocalPoint( const G4Navigator *This, G4ThreeVector pGlobalPoint)
{
	G4AffineTransform t =
		G4NavigationHistory_GetTopTransform( &(This->fHistory) );
	return G4AffineTransform_TransformPoint(&t, pGlobalPoint);
}

// ********************************************************************
// SetWorldVolume
//
// Sets the world (`topmost') volume
// ********************************************************************
//
MAYINLINE void G4Navigator_SetWorldVolume( G4Navigator *This, GEOMETRYLOC G4VPhysicalVolume* pWorld )
{
	This->fTopPhysical = pWorld;
	G4NavigationHistory_SetFirstEntry( &(This->fHistory), pWorld );
}

// ********************************************************************
// SetGeometrycallyLimitedStep
//
// Informs the navigator that the previous Step calculated
// by the geometry was taken in its entirety
// ********************************************************************
//
MAYINLINE void G4Navigator_SetGeometricallyLimitedStep( G4Navigator *This )
{
	This->fWasLimitedByGeometry = true;
}

// ********************************************************************
// ResetStackAndState
//
// Resets stack and minimum of navigator state `machine'
// ********************************************************************
//
INLINEFUNC
void G4Navigator_ResetStackAndState( G4Navigator *This )
{
	G4NavigationHistory_Reset( &(This->fHistory) );
	G4Navigator_ResetState( This );
}

// ********************************************************************
// VolumeType
// ********************************************************************
//
INLINEFUNC
EVolume G4Navigator_VolumeType( const G4Navigator *This, GEOMETRYLOC const G4VPhysicalVolume *pVol )
{
	(void)This;
	(void)pVol;
	
	return kNormal;
}

// ********************************************************************
// EnteredDaughterVolume
//
// To inform the caller if the track is entering a daughter volume
// ********************************************************************
//
/*INLINEFUNC
G4bool G4Navigator_EnteredDaughterVolume( const G4Navigator *This)
{
  return This->fEnteredDaughter;
}

// ********************************************************************
// ExitedMotherVolume
// ********************************************************************
//
INLINEFUNC
G4bool G4Navigator_ExitedMotherVolume( const G4Navigator *This)
{
  return This->fExitedMother;
}*/

// -------------------------- END "INLINES"

// ********************************************************************
// Constructor
// ********************************************************************
//
MAYINLINE void G4Navigator_ctor( G4Navigator *This )
{
	G4NavigationHistory_ctor( &(This->fHistory) );
#ifdef ENABLE_VOXEL_NAVIGATION
#ifdef NEW_NAVIGATION
	G4NewNavigation_ctor( &(This->fVoxelNav ) );
#else
	G4VoxelNavigation_ctor( &(This->fVoxelNav ) );
#endif
#endif
	
	G4Navigator_ResetStackAndState( This );
	
	This->fWasLimitedByGeometry = false;
	This->fTopPhysical = GEOMETRYNULL;
	This->fPushed = false;
	
	This->fStepEndPoint = G4ThreeVector_create( kInfinity, kInfinity, kInfinity );
}

// ********************************************************************
// LocateGlobalPointAndSetup
//
// Locate the point in the hierarchy return 0 if outside
// The direction is required 
//    - if on an edge shared by more than two surfaces 
//      (to resolve likely looping in tracking)
//    - at initial location of a particle
//      (to resolve potential ambiguity at boundary)
// 
// Flags on exit: (comments to be completed)
// fEntering         - True if entering `daughter' volume (or replica)
//                     whether daughter of last mother directly 
//                     or daughter of that volume's ancestor.
// ********************************************************************
//
MAYINLINE 
GEOMETRYLOC G4VPhysicalVolume* 
G4Navigator_LocateGlobalPointAndSetup(
		G4Navigator *This,
		G4ThreeVector globalPoint,
		const G4ThreeVector* pGlobalDirection,
		G4bool relativeSearch,
		G4bool ignoreDirection,
		GEOMETRYLOC float * Result
		)
{
  G4bool notKnownContained=true, noResult;
  GEOMETRYLOC G4VPhysicalVolume *targetPhysical;
  GEOMETRYLOC G4VSolid *targetSolid = GEOMETRYNULL;
  G4ThreeVector localPoint = G4ThreeVector_create(0,0,0);
  G4ThreeVector globalDirection = G4ThreeVector_create(0,0,0);
  EInside insideCode;
  G4bool considerDirection = (!ignoreDirection) || This->fLocatedOnEdge;

  if( considerDirection && pGlobalDirection != 0 )
  {
    globalDirection=*pGlobalDirection;
  }
  if ( 1 )
  { 	
     G4Navigator_ResetStackAndState( This );
  }
  else
  { 
    if ( This->fWasLimitedByGeometry )
    {
      This->fWasLimitedByGeometry = false;
      This->fEnteredDaughter = This->fEntering;   // Remember
      This->fExitedMother = This->fExiting;       // Remember
      if ( This->fExiting )
      { 
        if ( G4NavigationHistory_GetDepth( &(This->fHistory) ) )
        {  
          This->fBlockedPhysicalVolume = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
          G4NavigationHistory_BackLevel( &(This->fHistory) );
        }
        else
        {  
          This->fLastLocatedPointLocal = localPoint;
          This->fLocatedOutsideWorld = true;
          return GEOMETRYNULL;           // Have exited world volume
        }
        // A fix for the case where a volume is "entered" at an edge
        // and a coincident surface exists outside it.
        //  - This stops it from exiting further volumes and cycling
        //  - However ReplicaNavigator treats This case itself
        //
        if ( This->fLocatedOnEdge )
        { 
          This->fExiting= false;
        }
      }
      else
        if ( This->fEntering )
        {
		  G4NavigationHistory_NewLevel( &(This->fHistory), This->fBlockedPhysicalVolume, kNormal);
							
          This->fEntering = false;
          This->fBlockedPhysicalVolume = GEOMETRYNULL;
          
          G4AffineTransform t = G4NavigationHistory_GetTopTransform( &(This->fHistory) );
          localPoint = G4AffineTransform_TransformPoint(&t,globalPoint);
          
          notKnownContained = false;
        }
    }
    else
    {
      This->fBlockedPhysicalVolume = GEOMETRYNULL;
      This->fEntering = false;
      This->fEnteredDaughter = false;  // Full Step was not taken, did not enter
      This->fExiting = false;
      This->fExitedMother = false;     // Full Step was not taken, did not exit
    }

  }
  
  //
  // Search from top of history up through geometry until
  // containing volume found:
  // If on 
  // o OUTSIDE - Back up level, not/no longer exiting volumes
  // o SURFACE and EXITING - Back up level, setting new blocking no.s
  // else
  // o containing volume found
  //
  while (notKnownContained)
  {
	targetSolid =
	  G4LogicalVolume_GetSolid(
		  G4VPhysicalVolume_GetLogicalVolume(
				G4NavigationHistory_GetTopVolume(&(This->fHistory))));
	//targetSolid = fHistory.GetTopVolume()->GetLogicalVolume()->GetSolid();
	
	G4AffineTransform t = G4NavigationHistory_GetTopTransform( &(This->fHistory) );
	localPoint = G4AffineTransform_TransformPoint(&t,globalPoint);
	
	insideCode = G4VSolid_Inside(targetSolid,localPoint);
	
    if ( insideCode==kOutside )
    {
      if ( G4NavigationHistory_GetDepth( &(This->fHistory) ) )
      {
        This->fBlockedPhysicalVolume = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
        G4NavigationHistory_BackLevel( &(This->fHistory) );
        This->fExiting = false;
      }
      else
      {
        This->fLastLocatedPointLocal = localPoint;
        This->fLocatedOutsideWorld = true;
        return GEOMETRYNULL;         // Have exited world volume
      }
    }
    else
      if ( insideCode==kSurface )
      {
        G4bool isExiting = This->fExiting;
        if( (!This->fExiting)&&considerDirection )
        {
			// Figure out whether we are exiting This level's volume
			// by using the direction
			//
			G4bool directionExiting = false;

			G4AffineTransform t = G4NavigationHistory_GetTopTransform( &(This->fHistory) );
			G4ThreeVector localDirection =G4AffineTransform_TransformAxis(&t,globalDirection);

			G4ThreeVector normal = G4VSolid_SurfaceNormal(targetSolid, localPoint);
			directionExiting = G4ThreeVector_dot(normal,localDirection) > 0.0;
			isExiting = isExiting || directionExiting;
        }
        if( isExiting )
        {
          if ( G4NavigationHistory_GetDepth( &(This->fHistory) ) )
          {
            This->fBlockedPhysicalVolume = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
            G4NavigationHistory_BackLevel( &(This->fHistory) );
            //
            // Still on surface but exited volume not necessarily convex
            //
            This->fValidExitNormal = false;
          }
          else
          {
            This->fLastLocatedPointLocal = localPoint;
            This->fLocatedOutsideWorld = true;
            return GEOMETRYNULL;          // Have exited world volume
          }
        }
        else
        {
          notKnownContained=false;
        }
      }
      else
      {
        notKnownContained=false;
      }
  }  // END while (notKnownContained)
  //
  // Search downwards until deepest containing volume found,
  // blocking fBlockedPhysicalVolume/BlockedReplicaNum
  //
  // 3 Cases:
  //
  // o Parameterised daughters
  //   =>Must be one G4PVParameterised daughter & voxels
  // o Positioned daughters & voxels
  // o Positioned daughters & no voxels

  noResult = true;  // noResult should be renamed to 
                    // something like enteredLevel, as that is its meaning.
  do
  {
	  
	  
    // Determine `type' of current mother volume
    //
    targetPhysical = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
    
#ifdef ENABLE_COMBINED_NAVIGATION

	noResult = G4CombinedNavigation_LevelLocate(
				&(This->fVoxelNav),

				&(This->fHistory),
				This->fBlockedPhysicalVolume,
				globalPoint,
				pGlobalDirection,
				considerDirection,
				&localPoint);

#else
    
#ifdef ENABLE_VOXEL_NAVIGATION
		GEOMETRYLOC G4LogicalVolume *targetLogical = G4VPhysicalVolume_GetLogicalVolume(targetPhysical);
#ifdef NEW_NAVIGATION
		if ( G4LogicalVolume_GetVoxelHeader( targetLogical ) != GEOMETRYNULL )
		{
			// New navigation
			noResult =
				G4NewNavigation_LevelLocate(
					&(This->fVoxelNav),
					&(This->fHistory),
					This->fBlockedPhysicalVolume,
					globalPoint,
					pGlobalDirection,
					considerDirection,
					&localPoint);
		}



#else
    
    if ( G4LogicalVolume_GetVoxelHeader( targetLogical ) != GEOMETRYNULL )
    {
		// Voxel navigation
		noResult =
			G4VoxelNavigation_LevelLocate(
				&(This->fVoxelNav),
				&(This->fHistory),
				This->fBlockedPhysicalVolume,
				globalPoint,
				pGlobalDirection,
				considerDirection,
				&localPoint);
	}

	
#endif
	else
	{
	#endif
		// Normal navigation
    
		noResult = G4NormalNavigation_LevelLocate(
				&(This->fHistory),
				This->fBlockedPhysicalVolume,
				&globalPoint,
				pGlobalDirection,
				considerDirection,
				&localPoint);
				
	#ifdef ENABLE_VOXEL_NAVIGATION
	}
	#endif
	
#endif
			
    // LevelLocate returns true if it finds a daughter volume 
    // in which globalPoint is inside (or on the surface).

    if ( noResult )
    {
      // Entering a daughter after ascending
      //
      // The blocked volume is no longer valid - it was for another level
      //
      This->fBlockedPhysicalVolume = GEOMETRYNULL;

      // fEntering should be false -- else blockedVolume is assumed good.
      // fEnteredDaughter is used for ExitNormal
      //
      This->fEntering = false;
      This->fEnteredDaughter = true;
    }
  } while (noResult);

  This->fLastLocatedPointLocal = localPoint;

  This->fLocatedOutsideWorld= false;

  return targetPhysical;
}

// ********************************************************************
// LocateGlobalPointWithinVolume
//
// -> the state information of This Navigator and its subNavigators
//    is updated in order to start the next step at pGlobalpoint
// -> no check is performed whether pGlobalpoint is inside the 
//    original volume (This must be the case).
//
// Note: a direction could be added to the arguments, to aid in future
//       optional checking (via the old code below, flagged by OLD_LOCATE). 
//       [ This would be done only in verbose mode ]
// ********************************************************************
//
MAYINLINE void
G4Navigator_LocateGlobalPointWithinVolume( G4Navigator *This, G4ThreeVector pGlobalpoint)
{
	This->fLastLocatedPointLocal = G4Navigator_ComputeLocalPoint( This, pGlobalpoint );

#ifdef ENABLE_VOXEL_NAVIGATION

#ifdef NEW_NAVIGATION
	GEOMETRYLOC G4VPhysicalVolume* motherPhysical = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
	GEOMETRYLOC G4LogicalVolume* motherLogical = G4VPhysicalVolume_GetLogicalVolume( motherPhysical );
	GEOMETRYLOC G4SmartVoxelHeader* pVoxelHeader = G4LogicalVolume_GetVoxelHeader( motherLogical );

	if ( pVoxelHeader )
	{
		G4NewNavigation_VoxelLocate( &(This->fVoxelNav), pVoxelHeader, This->fLastLocatedPointLocal );
	}
#else
	// Update the state of the Sub Navigators 
	// - in particular any voxel information they store/cache
	//
	GEOMETRYLOC G4VPhysicalVolume* motherPhysical = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
	GEOMETRYLOC G4LogicalVolume* motherLogical = G4VPhysicalVolume_GetLogicalVolume( motherPhysical );
	GEOMETRYLOC G4SmartVoxelHeader* pVoxelHeader = G4LogicalVolume_GetVoxelHeader( motherLogical );

	if ( pVoxelHeader )
	{
		G4VoxelNavigation_VoxelLocate( &(This->fVoxelNav), pVoxelHeader, This->fLastLocatedPointLocal );
	}
#endif
#endif

	// Reset the state variables 
	//   - which would have been affected
	//     by the 'equivalent' call to LocateGlobalPointAndSetup
	//   - who's values have been invalidated by the 'move'.
	//
	This->fBlockedPhysicalVolume = GEOMETRYNULL; 
	This->fEntering = false;
	This->fEnteredDaughter = false;  // Boundary not encountered, did not enter
	This->fExiting = false;
	This->fExitedMother = false;     // Boundary not encountered, did not exit
}

// ********************************************************************
// ComputeStep
//
// Computes the next geometric Step: intersections with current
// mother and `daughter' volumes.
//
// NOTE:
//
// Flags on entry:
// --------------
// fValidExitNormal  - Normal of exited volume is valid (convex, not a 
//                     coincident boundary)
// fExitNormal       - Surface normal of exited volume
// fExiting          - True if have exited solid
//
// fBlockedPhysicalVolume - Ptr to exited volume (or 0)
// fBlockedReplicaNo - Replication no of exited volume
// fLastStepWasZero  - True if last Step size was zero.
//
// Flags on exit:
// -------------
// fValidExitNormal  - True if surface normal of exited volume is valid
// fExitNormal       - Surface normal of exited volume rotated to mothers
//                    reference system
// fExiting          - True if exiting mother
// fEntering         - True if entering `daughter' volume (or replica)
// fBlockedPhysicalVolume - Ptr to candidate (entered) volume
// fBlockedReplicaNo - Replication no of candidate (entered) volume
// fLastStepWasZero  - True if This Step size was zero.
// ********************************************************************
//
MAYINLINE 
G4double G4Navigator_ComputeStep(
		G4Navigator *This, 
		G4ThreeVector pGlobalpoint,
		G4ThreeVector pDirection,
		const G4double pCurrentProposedStepLength,
		G4double *pNewSafety
#ifdef NEW_NAVIGATION
		,SHAREDMEM int * Numbers_Of_Solid,

		SHAREDMEM int * Sum_Of_Solids,
		
		     SHAREDTYPE SolidInfo  * Solids,
		
			 SHAREDTYPE ResultInfo * Result_For_Current_Solid,
		
		     SHAREDTYPE FinalResult * Compacter_Result,

			 SHAREDMEM bool * noStepArray,

			 SHAREDMEM PointInformation * LocationArray,
			 GEOMETRYLOC G4SmartVoxelNode * nullVNode
			 
			 
			
#endif
			, G4bool cur_vol_local
#ifdef CHECK 
			,GEOMETRYLOC G4double * Result
#endif
		
		)
{ 
  G4ThreeVector localDirection = G4Navigator_ComputeLocalAxis(This,pDirection);
  G4double Step = DBL_MAX;
  GEOMETRYLOC G4VPhysicalVolume  *motherPhysical = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
  const G4double kCarTolerance = K_GEOMETRY_CAR_TOLERANCE;
  
  GEOMETRYLOC G4LogicalVolume *motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);

  G4ThreeVector newLocalPoint = G4Navigator_ComputeLocalPoint( This, pGlobalpoint);
  
  if( !G4ThreeVector_equal(newLocalPoint, This->fLastLocatedPointLocal) )
  {
    // Check whether the relocation is within safety
    //
    G4ThreeVector oldLocalPoint = This->fLastLocatedPointLocal;
    
    G4double moveLenSq = G4ThreeVector_diff2(newLocalPoint,oldLocalPoint);

    if ( moveLenSq >= kCarTolerance*kCarTolerance )
    {
      // Relocate the point within the same volume
      //
      G4Navigator_LocateGlobalPointWithinVolume( This, pGlobalpoint );
    }
  }
  
#ifdef ENABLE_COMBINED_NAVIGATION

	Step = G4CombinedNavigation_ComputeStep(
			&(This->fVoxelNav),
			This->fLastLocatedPointLocal,
			localDirection,
			pCurrentProposedStepLength,
			pNewSafety, // pointer
			&(This->fHistory),
			&(This->fValidExitNormal),
			&(This->fExitNormal),
			&(This->fExiting),
			&(This->fEntering),
			&(This->fBlockedPhysicalVolume));

#else
  
#ifdef ENABLE_VOXEL_NAVIGATION
#ifdef NEW_NAVIGATION
   if ( G4LogicalVolume_GetVoxelHeader(motherLogical) != GEOMETRYNULL )
  {
	Step = G4NewNavigation_ComputeStep(
			&(This->fVoxelNav),
			This->fLastLocatedPointLocal,
			localDirection,
			pCurrentProposedStepLength,
			pNewSafety, // pointer
			&(This->fHistory),
			&(This->fValidExitNormal),
			&(This->fExitNormal),
			&(This->fExiting),
			&(This->fEntering),
			&(This->fBlockedPhysicalVolume)
			 ,Numbers_Of_Solid,

			 Sum_Of_Solids,
		
		     Solids,
	
			Result_For_Current_Solid,
		
		    Compacter_Result,

			noStepArray,

			LocationArray,
			nullVNode,
			cur_vol_local
			#ifdef CHECK 
			, Result
			#endif
			);
#else
  if ( G4LogicalVolume_GetVoxelHeader(motherLogical) != GEOMETRYNULL )
  {
	if( cur_vol_local )
	Step = G4VoxelNavigation_ComputeStep(
			&(This->fVoxelNav),
			This->fLastLocatedPointLocal,
			localDirection,
			pCurrentProposedStepLength,
			pNewSafety, // pointer
			&(This->fHistory),
			&(This->fValidExitNormal),
			&(This->fExitNormal),
			&(This->fExiting),
			&(This->fEntering),
			&(This->fBlockedPhysicalVolume)
	
#ifdef CHECK 
			, Result
#endif			
			);
	else
	    return 0;
#endif
  }
  else
  {
#endif

	Step = G4NormalNavigation_ComputeStep(
			This->fLastLocatedPointLocal,
			localDirection,
			pCurrentProposedStepLength,
			pNewSafety, // pointer
			&(This->fHistory),
			&(This->fValidExitNormal),
			&(This->fExitNormal),
			&(This->fExiting),
			&(This->fEntering),
			&(This->fBlockedPhysicalVolume));
			
#ifdef ENABLE_VOXEL_NAVIGATION || NEW_NAVIGATION
  }
#endif

#endif

  // Remember last safety origin & value.
  //
  //This->fPreviousSftOrigin = pGlobalpoint;
  This->fPreviousSafety = *pNewSafety;

  // Count zero steps - one can occur due to changing momentum at a boundary
  //                  - one, two (or a few) can occur at common edges between
  //                    volumes
  //                  - more than two is likely a problem in the geometry
  //                    description or the Navigation 

  // Rule of thumb: likely at an Edge if two consecutive steps are zero,
  //                because at least two candidate volumes must have been
  //                checked
  //
  This->fLocatedOnEdge   = This->fLastStepWasZero && (Step==0.0);
  This->fLastStepWasZero = (Step==0.0);
  if (This->fPushed)  This->fPushed = This->fLastStepWasZero;

  // Handle large number of consecutive zero steps
  //
  if ( This->fLastStepWasZero )
  {
    This->fNumberZeroSteps++;
    if( This->fNumberZeroSteps > K_NAVIGATOR_ACTION_THRESHOLD_NOZEROSTEPS-1 )
    {
       // Act to recover This stuck track. Pushing it along direction
       //
       Step += 0.9*kCarTolerance;
       This->fPushed = true;
    }
    if( This->fNumberZeroSteps > K_NAVIGATOR_ABANDON_THRESHOLD_NOZEROSTEPS-1 )
    {
      // Must kill This stuck track
      //
	  myAbort();
    }
  }
  else
  {
    if (!This->fPushed)  This->fNumberZeroSteps = 0;
  }

  This->fEnteredDaughter = This->fEntering;   // I expect to enter a volume in This Step
  This->fExitedMother = This->fExiting;

  if( This->fExiting )
  {
    if(This->fValidExitNormal)
    {
      // Convention: fExitNormal is in the 'grand-mother' coordinate system
      //
      This->fGrandMotherExitNormal= This->fExitNormal;
    }
    else
    {
      // We must calculate the normal anyway (in order to have it if requested)
      //
      
      G4ThreeVector finalLocalPoint =
		G4ThreeVector_saxpy( Step, localDirection, This->fLastLocatedPointLocal );

      // Now fGrandMotherExitNormal is in the 'grand-mother' coordinate system
      //
      
      This->fGrandMotherExitNormal =
		G4VSolid_SurfaceNormal(
			G4LogicalVolume_GetSolid(motherLogical),finalLocalPoint);

      G4RotationMatrix mRot = G4VPhysicalVolume_GetObjectRotationValue(motherPhysical);
      G4RotationMatrix inv = G4RotationMatrix_inverse(&mRot);

      This->fGrandMotherExitNormal 
       = G4RotationMatrix_apply(&inv,This->fGrandMotherExitNormal);
      //*= (*mRot).inverse();
    }
  }
  
  
  This->fStepEndPoint = 
	G4ThreeVector_saxpy(Step, pDirection, pGlobalpoint );

  if( (Step == pCurrentProposedStepLength) && (!This->fExiting) && (!This->fEntering) )
  {
    // This if Step is not really limited by the geometry.
    // The Navigator is obliged to return "infinity"
    //
    Step = kInfinity;
  }

  return Step;
}

// ********************************************************************
// GetLocalExitNormal
//
// Obtains the Normal vector to a surface (in local coordinates)
// pointing out of previous volume and into current volume
// ********************************************************************
//
/*INLINEFUNC
G4ThreeVector G4Navigator_GetLocalExitNormal( G4Navigator *This, G4bool* valid )
{
  G4ThreeVector ExitNormal = G4ThreeVector_create(0.,0.,0.);

  if ( G4Navigator_EnteredDaughterVolume(This) )
  {
    ExitNormal=
      G4ThreeVector_negation(
	   G4VSolid_SurfaceNormal(
		G4LogicalVolume_GetSolid(
		 G4VPhysicalVolume_GetLogicalVolume(
		  G4NavigationHistory_GetTopVolume( &(This->fHistory) ))),
		This->fLastLocatedPointLocal));
		
    //-(G4NavigationHistory_GetTopVolume( &(This->fHistory) )->
	// GetLogicalVolume()->GetSolid()->SurfaceNormal(This->fLastLocatedPointLocal));
    *valid = true;
  }
  else
  {
    if( This->fExitedMother )
    {
      ExitNormal = This->fGrandMotherExitNormal;
      *valid = true;
    }
    else
    {
      // We are not at a boundary.
      // ExitNormal remains (0,0,0)
      //
      *valid = false;
    }
  }
  return ExitNormal;
}*/

// ********************************************************************
// ComputeSafety
//
// It assumes that it will be 
//  i) called at the Point in the same volume as the EndPoint of the
//     ComputeStep.
// ii) after (or at the end of) ComputeStep OR after the relocation.
// ********************************************************************
//
/*MAYINLINE 
G4double G4Navigator_ComputeSafety(
		G4Navigator *This,
		G4ThreeVector pGlobalpoint )
{
  G4double newSafety = 0.0;

  const G4double kCarTolerance = K_GEOMETRY_CAR_TOLERANCE;

  G4double distEndpointSq = G4ThreeVector_diff2(pGlobalpoint, This->fStepEndPoint); 
  G4bool   stayedOnEndpoint  = distEndpointSq < kCarTolerance*kCarTolerance; 
  G4bool   endpointOnSurface = This->fEnteredDaughter || This->fExitedMother;

  if( !(endpointOnSurface && stayedOnEndpoint) )
  {
    // Pseudo-relocate to This point (updates voxel information only)
    //
    G4Navigator_LocateGlobalPointWithinVolume( This, pGlobalpoint );
      // --->> Danger: Side effects on sub-navigator voxel information <<---
      //       Could be replaced again by 'granular' calls to sub-navigator
      //       locates (similar side-effects, but faster.  
      //       Solutions:
      //        1) Re-locate (to where?)
      //        2) Insure that the methods using (G4ComputeStep?)
      //           does a relocation (if information is disturbed only ?)

    G4ThreeVector localPoint = G4Navigator_ComputeLocalPoint(This,pGlobalpoint);
    
#ifdef ENABLE_COMBINED_NAVIGATION

	newSafety = G4CombinedNavigation_ComputeSafety(
			&(This->fVoxelNav),localPoint,&(This->fHistory));

#else
    
#ifdef ENABLE_VOXEL_NAVIGATION
    GEOMETRYLOC G4VPhysicalVolume *motherPhysical = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
    GEOMETRYLOC G4LogicalVolume *motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
    GEOMETRYLOC G4SmartVoxelHeader* pVoxelHeader = G4LogicalVolume_GetVoxelHeader(motherLogical);

	if ( pVoxelHeader != GEOMETRYNULL )
	{
		// Voxel navigation
		newSafety = G4VoxelNavigation_ComputeSafety(
			&(This->fVoxelNav),localPoint,&(This->fHistory));
	}
	else
	{
#endif
		// Normal navigation
		newSafety = G4NormalNavigation_ComputeSafety(
			localPoint,&(This->fHistory));
			
#ifdef ENABLE_VOXEL_NAVIGATION
	}
#endif
#endif
  }
  else // if( endpointOnSurface && stayedOnEndpoint )
  {
    newSafety = 0.0; 
  }

  // Remember last safety origin & value
  //
  //This->fPreviousSftOrigin = pGlobalpoint;
  This->fPreviousSafety = newSafety; 

  return newSafety;
}*/
