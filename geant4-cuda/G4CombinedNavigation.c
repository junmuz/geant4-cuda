
/** These functions combine G4NormalNavigation and G4VoxelNavigation */

#ifndef G4COMBINEDNAVIGATION_C
#define G4COMBINEDNAVIGATION_C

// ********************************************************************
// ComputeStep
// ********************************************************************
//
MAYINLINE G4double
G4CombinedNavigation_ComputeStep(
			G4VoxelNavigation *vox,
			G4ThreeVector localPoint,
			G4ThreeVector localDirection,
			const G4double currentProposedStepLength,
			G4double *newSafety,
			G4NavigationHistory *history,
			G4bool *validExitNormal,
			G4ThreeVector *exitNormal,
			G4bool *exiting,
			G4bool *entering,
			GEOMETRYLOC G4VPhysicalVolume *(*pBlockedPhysical) )
{
  GEOMETRYLOC G4VPhysicalVolume *motherPhysical, *samplePhysical,
	*blockedExitedVol = GEOMETRYNULL;
  GEOMETRYLOC G4LogicalVolume *motherLogical;
  GEOMETRYLOC G4VSolid *motherSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int sampleNo;

  G4bool initialNode, noStep;
  GEOMETRYLOC const G4SmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;

  motherPhysical = G4NavigationHistory_GetTopVolume( history );
  motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = G4LogicalVolume_GetSolid(motherLogical);
  
  const G4bool voxelized = G4LogicalVolume_GetVoxelHeader(motherLogical) != GEOMETRYNULL;

  //
  // Compute mother safety
  //

  motherSafety = G4VSolid_DistanceToOut(motherSolid, localPoint);
  ourSafety = motherSafety;                 // Working isotropic safety

  //
  // Compute daughter safeties & intersections
  //

  // Exiting normal optimisation
  //
  if ( *exiting && *validExitNormal )
  {
    if ( G4ThreeVector_dot(localDirection,*exitNormal)>=kMinExitingNormalCosine )
    {
      // Block exited daughter volume
      //
      blockedExitedVol = *pBlockedPhysical;
      ourSafety = 0;
    }
  }
  *exiting = false;
  *entering = false;

  G4int localNoDaughters = G4LogicalVolume_GetNoDaughters(motherLogical);
  
  initialNode = true;
  noStep = true;

  while (noStep)
  {
	if ( voxelized )
	{
		curVoxelNode = vox->fVoxelNode;
		curNoVolumes = G4VoxelNode_GetNoContained(curVoxelNode);
	}
	else
	{
		curNoVolumes = localNoDaughters;
		noStep = false;
	}
	
    for (contentNo=curNoVolumes-1; contentNo>=0; contentNo--)
    {
		if (voxelized)
		{
			sampleNo = G4VoxelNode_GetVolume(curVoxelNode,contentNo);
		}
		else
		{
			sampleNo = contentNo;
		}
		
		samplePhysical = G4LogicalVolume_GetDaughter(motherLogical,sampleNo);
        if ( samplePhysical!=blockedExitedVol )
        {
		  G4AffineTransform sampleTf =
			G4AffineTransform_create_full(
				G4VPhysicalVolume_GetObjectRotationValue(samplePhysical),
				G4VPhysicalVolume_GetTranslation(samplePhysical));
				
          G4AffineTransform_Invert(&sampleTf);
          
          const G4ThreeVector samplePoint =
				G4AffineTransform_TransformPoint(&sampleTf,localPoint);
                     
          GEOMETRYLOC const G4VSolid *sampleSolid =
				G4LogicalVolume_GetSolid(
					G4VPhysicalVolume_GetLogicalVolume(
						samplePhysical ));
          
          const G4double sampleSafety =
			G4VSolid_DistanceToIn(sampleSolid,samplePoint);

          if ( sampleSafety<ourSafety )
          {
            ourSafety = sampleSafety;
          }
          if ( sampleSafety<=ourStep )
          {
            sampleDirection =
				G4AffineTransform_TransformAxis( &sampleTf, localDirection );
            
            G4double sampleStep =
				G4VSolid_DistanceToIn_full(sampleSolid, samplePoint, sampleDirection);

            if ( sampleStep<=ourStep )
            {
              ourStep = sampleStep;
              *entering = true;
              *exiting = false;
              *pBlockedPhysical = samplePhysical;

            }
        }
      }
    }
    if (initialNode)
    {
	  if (voxelized)
	  {
		  initialNode = false;
		  voxelSafety = G4VoxelNavigation_ComputeVoxelSafety(vox,localPoint);
		  if ( voxelSafety<ourSafety )
		  {
			ourSafety = voxelSafety;
		  }
	  }
      if ( currentProposedStepLength<ourSafety )
      {
        // Guaranteed physics limited
        //      
        noStep = false;
        *entering = false;
        *exiting = false;
        *pBlockedPhysical = GEOMETRYNULL;
        ourStep = kInfinity;
      }
      else
      {
        //
        // Compute mother intersection if required
        //
        if ( motherSafety<=ourStep )
        {
          G4double motherStep =
			G4VSolid_DistanceToOut_full( motherSolid, localPoint, localDirection,
                                         true, validExitNormal, exitNormal);

          if ( motherStep<=ourStep )
          {
            ourStep = motherStep;
            *exiting = true;
            *entering = false;
            if ( *validExitNormal )
            {
				G4RotationMatrix rot = G4VPhysicalVolume_GetObjectRotationValue(motherPhysical);
				G4RotationMatrix inv = G4RotationMatrix_inverse(&rot);
				*exitNormal = G4RotationMatrix_apply( &inv, *exitNormal );
            }
          }
          else
          {
            *validExitNormal = false;
          }
        }
      }
      *newSafety = ourSafety;
    }
    if ( noStep ) // not voxelized => noStep == false
    {
      noStep = G4VoxelNavigation_LocateNextVoxel(vox, localPoint, localDirection, ourStep);
    }
  }  // end -while (noStep)- loop

  return ourStep;
}


// ********************************************************************
// LevelLocate
// ********************************************************************
//
INLINEFUNC G4bool
G4CombinedNavigation_LevelLocate(
			G4VoxelNavigation *vox,
			G4NavigationHistory* history,
			GEOMETRYLOC const G4VPhysicalVolume* blockedVol,
			G4ThreeVector globalPoint,
			const G4ThreeVector* globalDirection,
			const G4bool pLocatedOnEdge, 
			G4ThreeVector *localPoint )
{
  GEOMETRYLOC G4SmartVoxelHeader *targetVoxelHeader;
  GEOMETRYLOC G4SmartVoxelNode *targetVoxelNode;
  GEOMETRYLOC G4VPhysicalVolume *targetPhysical, *samplePhysical;
  GEOMETRYLOC G4LogicalVolume *targetLogical;
  GEOMETRYLOC G4VSolid *sampleSolid;
  G4ThreeVector samplePoint;
  G4int targetNoDaughters, sampleLogicalNo;
  
  targetPhysical = G4NavigationHistory_GetTopVolume(history);
  targetLogical = G4VPhysicalVolume_GetLogicalVolume(targetPhysical);
  targetVoxelHeader = G4LogicalVolume_GetVoxelHeader(targetLogical);
  
  const G4bool voxelized = targetVoxelHeader != GEOMETRYNULL;

  // Find the voxel containing the point
  //
  if ( voxelized )
  {
	targetVoxelNode =
	  G4VoxelNavigation_VoxelLocate(vox,targetVoxelHeader,*localPoint);

	targetNoDaughters=G4VoxelNode_GetNoContained(targetVoxelNode);
  }
  else
  {
	targetNoDaughters = G4LogicalVolume_GetNoDaughters(targetLogical);
  }
  
  if ( targetNoDaughters==0 ) return false;

  //
  // Search daughters in volume
  //
  for ( int sampleNo=targetNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    if ( voxelized )
    {
    	sampleLogicalNo = G4VoxelNode_GetVolume(targetVoxelNode,sampleNo);
    }
    else
    {
    	sampleLogicalNo = sampleNo;
    }
	
    samplePhysical =
    	G4LogicalVolume_GetDaughter( targetLogical, sampleLogicalNo );
                     
    if ( samplePhysical!=blockedVol )
    {
      // Setup history
      //
      G4NavigationHistory_NewLevel(history, samplePhysical, kNormal);
      
      sampleSolid =
		G4LogicalVolume_GetSolid(
			G4VPhysicalVolume_GetLogicalVolume( samplePhysical ));
			
      G4AffineTransform tf =
			G4NavigationHistory_GetTopTransform( history );
			
      samplePoint =
		G4AffineTransform_TransformPoint( &tf, globalPoint );

      if( G4AuxiliaryNavServices_CheckPointOnSurface(
			sampleSolid, samplePoint, globalDirection, 
			tf, pLocatedOnEdge) )
      {
        // Enter this daughter
        //
        *localPoint = samplePoint;
        return true;
      }
      else
      {
		  G4NavigationHistory_BackLevel( history );
      }
    }
  }
  return false;
}

// ********************************************************************
// ComputeSafety
//
// Calculates the isotropic distance to the nearest boundary from the
// specified point in the local coordinate system. 
// The localpoint utilised must be within the current volume.
// ********************************************************************
//
/*MAYINLINE G4double
G4CombinedNavigation_ComputeSafety(
			G4VoxelNavigation *vox,
			G4ThreeVector localPoint,
			const G4NavigationHistory *history )
{	
  GEOMETRYLOC G4VPhysicalVolume *motherPhysical, *samplePhysical;
  GEOMETRYLOC G4LogicalVolume *motherLogical;
  GEOMETRYLOC G4VSolid *motherSolid;
  G4double motherSafety, ourSafety;
  G4int sampleNo;
  GEOMETRYLOC const G4SmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;

  motherPhysical = G4NavigationHistory_GetTopVolume(history);
  motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = G4LogicalVolume_GetSolid(motherLogical);
  const G4bool voxelized = G4LogicalVolume_GetVoxelHeader(motherLogical) != GEOMETRYNULL;

  //
  // Compute mother safety
  //

  motherSafety = G4VSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety;                 // Working isotropic safety

  //
  // Compute daughter safeties 
  //

  if (voxelized)
  {
	//  Look only inside the current Voxel only (in the first version).
	//
	curVoxelNode = vox->fVoxelNode;
	curNoVolumes = G4VoxelNode_GetNoContained(curVoxelNode);
  }
  else
  {
	curNoVolumes = G4LogicalVolume_GetNoDaughters(motherLogical);
  }

  for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
  {
    if ( voxelized )
		sampleNo = G4VoxelNode_GetVolume(curVoxelNode,contentNo);
	else
		sampleNo = contentNo;
		
    samplePhysical = G4LogicalVolume_GetDaughter(motherLogical,sampleNo);

    G4AffineTransform sampleTf =
		G4AffineTransform_create_full(
				G4VPhysicalVolume_GetObjectRotationValue(samplePhysical),
				G4VPhysicalVolume_GetTranslation(samplePhysical));
	
	G4AffineTransform_Invert( &sampleTf );
	
    const G4ThreeVector samplePoint =
		G4AffineTransform_TransformPoint(&sampleTf, localPoint);
	
    GEOMETRYLOC const G4VSolid *sampleSolid =
		G4LogicalVolume_GetSolid(
			G4VPhysicalVolume_GetLogicalVolume( samplePhysical ));
	
    G4double sampleSafety =
		G4VSolid_DistanceToIn(sampleSolid, samplePoint);
    
    if ( sampleSafety<ourSafety )
    {
      ourSafety = sampleSafety;
    }

  }
  if ( voxelized )
  {
	voxelSafety = G4VoxelNavigation_ComputeVoxelSafety(vox,localPoint);
	if ( voxelSafety<ourSafety )
	{
	ourSafety = voxelSafety;
	}
  }
  return ourSafety;
}*/

#endif
