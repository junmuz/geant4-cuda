
/**
 * G4VoxelNavigation implementation
 * based on G4VoxelNavigation.(i)cc of Geant 4.9.3
 */

#include "G4VoxelNavigation.h"
#include "G4Navigator.h"

// FWD, defined in normal navigation
INLINEFUNC G4bool
G4AuxiliaryNavServices_CheckPointOnSurface(
				     GEOMETRYLOC const G4VSolid* sampleSolid, 
                     G4ThreeVector localPoint, 
                     const G4ThreeVector* globalDirection, 
                     G4AffineTransform sampleTransform,
                     const G4bool locatedOnEdge);
                     
INLINEFUNC G4bool
G4AuxiliaryNavServices_CheckPointExiting(
				   GEOMETRYLOC const G4VSolid* sampleSolid, 
                   G4ThreeVector localPoint, 
                   const G4ThreeVector* globalDirection, 
                   G4AffineTransform sampleTransform );

// class G4VoxelNavigation Inline implementation
//
// --------------------------------------------------------------------

// ********************************************************************
// VoxelLocate
// ********************************************************************
//
INLINEFUNC
GEOMETRYLOC G4SmartVoxelNode*
G4VoxelNavigation_VoxelLocate(
			G4VoxelNavigation *This,
			GEOMETRYLOC G4SmartVoxelHeader* pHead,
			G4ThreeVector localPoint )
{
  GEOMETRYLOC G4SmartVoxelHeader *targetVoxelHeader=pHead;
  GEOMETRYLOC G4SmartVoxelNode *targetVoxelNode = GEOMETRYNULL;
  GEOMETRYLOC const G4SmartVoxelProxy *sampleProxy;
  EAxis targetHeaderAxis;
  G4double targetHeaderMin, targetHeaderNodeWidth;
  G4int targetHeaderNoSlices, targetNodeNo;

  This->fVoxelDepth = 0;

  while ( targetVoxelNode == GEOMETRYNULL )
  {
    targetHeaderAxis = G4VoxelHeader_GetAxis(targetVoxelHeader);
    targetHeaderNoSlices = G4VoxelHeader_GetNoSlices(targetVoxelHeader);
    targetHeaderMin = G4VoxelHeader_GetMinExtent(targetVoxelHeader);
    targetHeaderNodeWidth =
		(G4VoxelHeader_GetMaxExtent(targetVoxelHeader)-targetHeaderMin)
                          / targetHeaderNoSlices;
    targetNodeNo = (G4int)(
		(G4ThreeVector_coord(localPoint,targetHeaderAxis)-targetHeaderMin)
                          / targetHeaderNodeWidth);
                          
    // Rounding protection
    //
    if ( targetNodeNo<0 )
    {
		targetNodeNo = 0;
    }
    else if ( targetNodeNo>=targetHeaderNoSlices )
	{
		targetNodeNo = targetHeaderNoSlices-1;
	}
         
    // Stack info for stepping
    //
    
    This->fVoxelAxisStack[This->fVoxelDepth] = targetHeaderAxis;
    This->fVoxelNoSlicesStack[This->fVoxelDepth] = targetHeaderNoSlices;
    This->fVoxelSliceWidthStack[This->fVoxelDepth] = targetHeaderNodeWidth;
    This->fVoxelNodeNoStack[This->fVoxelDepth] = targetNodeNo;
    This->fVoxelHeaderStack[This->fVoxelDepth] = targetVoxelHeader;
    sampleProxy = G4VoxelHeader_GetSlice(targetVoxelHeader, targetNodeNo);

    if ( G4VoxelProxy_IsNode(sampleProxy) )
    {
      targetVoxelNode = G4VoxelProxy_GetNode(sampleProxy);
    }
    else
    {
      targetVoxelHeader = G4VoxelProxy_GetHeader(sampleProxy);
      This->fVoxelDepth++;
      myAssert(This->fVoxelDepth < K_MAX_VOXEL_STACK_DEPTH);
    }
  }
  
  This->fVoxelNode = targetVoxelNode;
  return targetVoxelNode;
}


// ********************************************************************
// LocateNextVoxel
//
// Finds the next voxel from the current voxel and point
// in the specified direction
//
// Returns false if all voxels considered
//              [current Step ends inside same voxel or leaves all voxels]
//         true  otherwise
//              [the information on the next voxel is put into the set of
//               fVoxel* variables & "stacks"] 
// ********************************************************************
// 
MAYINLINE
G4bool
G4VoxelNavigation_LocateNextVoxel(
			G4VoxelNavigation *This,
			G4ThreeVector localPoint,
			G4ThreeVector localDirection,
			const G4double currentStep )
{
  GEOMETRYLOC G4SmartVoxelHeader *workHeader=GEOMETRYNULL, *newHeader=GEOMETRYNULL;
  GEOMETRYLOC G4SmartVoxelProxy *newProxy=GEOMETRYNULL;
  GEOMETRYLOC G4SmartVoxelNode *newVoxelNode= GEOMETRYNULL;
  G4ThreeVector targetPoint, voxelPoint;
  G4double workNodeWidth, workMinExtent, workCoord;
  G4double minVal, maxVal, newDistance=0.;
  G4double newHeaderMin, newHeaderNodeWidth;
  G4int depth=0, newDepth=0, workNodeNo=0, newNodeNo=0, newHeaderNoSlices=0;
  EAxis workHeaderAxis, newHeaderAxis;
  G4bool isNewVoxel=false;
  
  G4double currentDistance = currentStep;

  // Determine if end of Step within current voxel
  //
  for (depth=0; depth<This->fVoxelDepth; depth++)
  {
    targetPoint =
		G4ThreeVector_saxpy(currentDistance,localDirection,localPoint);
			
    newDistance = currentDistance;
    workHeader = This->fVoxelHeaderStack[depth];
    workHeaderAxis = This->fVoxelAxisStack[depth];
    workNodeNo = This->fVoxelNodeNoStack[depth];
    workNodeWidth = This->fVoxelSliceWidthStack[depth];
    workMinExtent = G4VoxelHeader_GetMinExtent(workHeader);
    workCoord = G4ThreeVector_coord(targetPoint,workHeaderAxis);
    minVal = workMinExtent+workNodeNo*workNodeWidth;

    if ( minVal<=workCoord+K_GEOMETRY_CAR_TOLERANCE*0.5 )
    {
      maxVal = minVal+workNodeWidth;
      if ( maxVal<=workCoord-K_GEOMETRY_CAR_TOLERANCE*0.5 )
      {
        // Must consider next voxel
        //
        newNodeNo = workNodeNo+1;
        newHeader = workHeader;
        newDistance = (maxVal-G4ThreeVector_coord(localPoint,workHeaderAxis))
                    / G4ThreeVector_coord(localDirection,workHeaderAxis);
        isNewVoxel = true;
        newDepth = depth;
      }
    }
    else
    {
      newNodeNo = workNodeNo-1;
      newHeader = workHeader;
      newDistance = (minVal-G4ThreeVector_coord(localPoint,workHeaderAxis))
                  / G4ThreeVector_coord(localDirection,workHeaderAxis);
      isNewVoxel = true;
      newDepth = depth;
    }
    currentDistance = newDistance;
  }
  targetPoint = 
	G4ThreeVector_saxpy(currentDistance,localDirection,localPoint);

  // Check if end of Step within collected boundaries of current voxel
  //
  depth = This->fVoxelDepth;
  {
    workHeader = This->fVoxelHeaderStack[depth];
    workHeaderAxis = This->fVoxelAxisStack[depth];
    workNodeNo = This->fVoxelNodeNoStack[depth];
    workNodeWidth = This->fVoxelSliceWidthStack[depth];
    workMinExtent = G4VoxelHeader_GetMinExtent(workHeader);
    workCoord = G4ThreeVector_coord(targetPoint,workHeaderAxis);
    minVal = workMinExtent+G4VoxelNode_GetMinEquivalentSliceNo(This->fVoxelNode)*workNodeWidth;

    if ( minVal<=workCoord+K_GEOMETRY_CAR_TOLERANCE*0.5 )
    {
      maxVal = workMinExtent+(G4VoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)+1)
                            *workNodeWidth;
      if ( maxVal<=workCoord-K_GEOMETRY_CAR_TOLERANCE*0.5 )
      {
        newNodeNo = G4VoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)+1;
        newHeader = workHeader;
        newDistance = (maxVal-G4ThreeVector_coord(localPoint,workHeaderAxis))
                    / G4ThreeVector_coord(localDirection,workHeaderAxis);
        isNewVoxel = true;
        newDepth = depth;
      }
    }
    else
    {
      newNodeNo = G4VoxelNode_GetMinEquivalentSliceNo(This->fVoxelNode)-1;
      newHeader = workHeader;
      newDistance = (minVal-G4ThreeVector_coord(localPoint,workHeaderAxis))
                  / G4ThreeVector_coord(localDirection,workHeaderAxis);
      isNewVoxel = true;
      newDepth = depth;
    }
    currentDistance = newDistance;
  }
  if (isNewVoxel)
  {
    // Compute new voxel & adjust voxel stack
    //
    // newNodeNo=Candidate node no at 
    // newDepth =refinement depth of crossed voxel boundary
    // newHeader=Header for crossed voxel
    // newDistance=distance to crossed voxel boundary (along the track)
    //
    if ( (newNodeNo<0) || (newNodeNo>=G4VoxelHeader_GetNoSlices(newHeader)))
    {
      // Leaving mother volume
      //
      isNewVoxel = false;
    }
    else
    {
      // Compute intersection point on the least refined
      // voxel boundary that is hit
      //
      voxelPoint = G4ThreeVector_saxpy(newDistance,localDirection,localPoint);
      myAssert(newDepth < K_MAX_VOXEL_STACK_DEPTH);
      This->fVoxelNodeNoStack[newDepth] = newNodeNo;
      This->fVoxelDepth = newDepth;
      newVoxelNode = 0;
      while ( newVoxelNode == GEOMETRYNULL )
      {
        newProxy = G4VoxelHeader_GetSlice(newHeader,newNodeNo);
        if ( G4VoxelProxy_IsNode(newProxy) )
        {
          newVoxelNode = G4VoxelProxy_GetNode(newProxy);
        }
        else
        {
          This->fVoxelDepth++;
          myAssert(This->fVoxelDepth < K_MAX_VOXEL_STACK_DEPTH);
          newHeader = G4VoxelProxy_GetHeader(newProxy);
          newHeaderAxis = G4VoxelHeader_GetAxis(newHeader);
          newHeaderNoSlices = G4VoxelHeader_GetNoSlices(newHeader);
          newHeaderMin = G4VoxelHeader_GetMinExtent(newHeader);
          newHeaderNodeWidth =
			(G4VoxelHeader_GetMaxExtent(newHeader)-newHeaderMin)
                             / newHeaderNoSlices;
          newNodeNo = (G4int)(
			(G4ThreeVector_coord(voxelPoint,newHeaderAxis)-newHeaderMin)
                             / newHeaderNodeWidth );
          // Rounding protection
          //
          if ( newNodeNo<0 )
          {
            newNodeNo=0;
          }
          else if ( newNodeNo>=newHeaderNoSlices )
               {
                 newNodeNo = newHeaderNoSlices-1;
               }
          // Stack info for stepping
          //
          This->fVoxelAxisStack[This->fVoxelDepth] = newHeaderAxis;
          This->fVoxelNoSlicesStack[This->fVoxelDepth] = newHeaderNoSlices;
          This->fVoxelSliceWidthStack[This->fVoxelDepth] = newHeaderNodeWidth;
          This->fVoxelNodeNoStack[This->fVoxelDepth] = newNodeNo;
          This->fVoxelHeaderStack[This->fVoxelDepth] = newHeader;
        }
      }
      This->fVoxelNode = newVoxelNode;
    }
  }
  return isNewVoxel;        
}

#define std_min(a,b) (((a)<(b))?(a):(b))

// ********************************************************************
// ComputeVoxelSafety
//
// Computes safety from specified point to voxel boundaries
// using already located point
// o collected boundaries for most derived level
// o adjacent boundaries for previous levels
// ********************************************************************
//
MAYINLINE
G4double
G4VoxelNavigation_ComputeVoxelSafety(
			const G4VoxelNavigation *This,
			G4ThreeVector localPoint)
{
  GEOMETRYLOC G4SmartVoxelHeader *curHeader;
  G4double voxelSafety, curNodeWidth;
  G4double curNodeOffset, minCurCommonDelta, maxCurCommonDelta;
  G4int minCurNodeNoDelta, maxCurNodeNoDelta;
  G4int localVoxelDepth, curNodeNo;
  EAxis curHeaderAxis;

  localVoxelDepth = This->fVoxelDepth;

  curHeader = This->fVoxelHeaderStack[localVoxelDepth];
  curHeaderAxis = This->fVoxelAxisStack[localVoxelDepth];
  curNodeNo = This->fVoxelNodeNoStack[localVoxelDepth];
  curNodeWidth = This->fVoxelSliceWidthStack[localVoxelDepth];
  
  // Compute linear intersection distance to boundaries of max/min
  // to collected nodes at current level
  //
  curNodeOffset = curNodeNo*curNodeWidth;
  
  maxCurNodeNoDelta = G4VoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)-curNodeNo;
  minCurNodeNoDelta = curNodeNo-G4VoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode);
  minCurCommonDelta = G4ThreeVector_coord(localPoint,curHeaderAxis)
	- G4VoxelHeader_GetMinExtent(curHeader) - curNodeOffset;
	
  maxCurCommonDelta = curNodeWidth-minCurCommonDelta;

  if ( minCurNodeNoDelta<maxCurNodeNoDelta )
  {
    voxelSafety = minCurNodeNoDelta*curNodeWidth;
    voxelSafety += minCurCommonDelta;
  }
  else if (maxCurNodeNoDelta < minCurNodeNoDelta)
       {
         voxelSafety = maxCurNodeNoDelta*curNodeWidth;
         voxelSafety += maxCurCommonDelta;
        }
        else    // (maxCurNodeNoDelta == minCurNodeNoDelta)
        {
          voxelSafety = minCurNodeNoDelta*curNodeWidth;
          voxelSafety += std_min(minCurCommonDelta,maxCurCommonDelta);
        }

  // Compute isotropic safety to boundaries of previous levels
  // [NOT to collected boundaries]
  //
  while ( (localVoxelDepth>0) && (voxelSafety>0) )
  {
    localVoxelDepth--;
    curHeader = This->fVoxelHeaderStack[localVoxelDepth];
    curHeaderAxis = This->fVoxelAxisStack[localVoxelDepth];
    curNodeNo = This->fVoxelNodeNoStack[localVoxelDepth];
    curNodeWidth = This->fVoxelSliceWidthStack[localVoxelDepth];
    curNodeOffset = curNodeNo*curNodeWidth;
    minCurCommonDelta = G4ThreeVector_coord(localPoint,curHeaderAxis)
                        - G4VoxelHeader_GetMinExtent(curHeader) - curNodeOffset;
                        
    maxCurCommonDelta = curNodeWidth-minCurCommonDelta;
    
    if ( minCurCommonDelta<voxelSafety )
    {
      voxelSafety = minCurCommonDelta;
    }
    if ( maxCurCommonDelta<voxelSafety )
    {
      voxelSafety = maxCurCommonDelta;
    }
  }
  if ( voxelSafety<0 )
  {
    voxelSafety = 0;
  }

  return voxelSafety;
}


// ********************************************************************
// Constructor
// ********************************************************************
//
MAYINLINE
void G4VoxelNavigation_ctor( G4VoxelNavigation *This )
{
	This->fVoxelDepth = -1;
	This->fVoxelNode = GEOMETRYNULL;
	
#ifdef USE_BLIST
	This->fBlist = NULL;
	This->fBlistSz = 0;
#endif
}

#ifdef USE_BLIST
INLINEFUNC
void G4VoxelNavigation_EnlargeAndResetBlist( G4VoxelNavigation *This, G4int n )
{
	// TODO: knowingly leaking memory
	if ( This->fBlistSz < n )
	{
		This->fBlist = realloc( This->fBlist, sizeof(char)*n );
		This->fBlistSz = n;
	}
	memset( This->fBlist, 0, sizeof(char)*n );
}
#endif

#ifndef ENABLE_COMBINED_NAVIGATION

// ********************************************************************
// LevelLocate
// ********************************************************************
//
INLINEFUNC
G4bool
G4VoxelNavigation_LevelLocate(
			G4VoxelNavigation *This,
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
  G4int targetNoDaughters;
  
  targetPhysical = G4NavigationHistory_GetTopVolume(history);
  targetLogical = G4VPhysicalVolume_GetLogicalVolume(targetPhysical);
  targetVoxelHeader = G4LogicalVolume_GetVoxelHeader(targetLogical);

  // Find the voxel containing the point
  //
  targetVoxelNode =
	G4VoxelNavigation_VoxelLocate(This,targetVoxelHeader,*localPoint);

  targetNoDaughters=G4VoxelNode_GetNoContained(targetVoxelNode);
  if ( targetNoDaughters==0 ) return false;

  //
  // Search daughters in volume
  //
  for ( int sampleNo=targetNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    samplePhysical =
		G4LogicalVolume_GetDaughter( targetLogical, 
			G4VoxelNode_GetVolume(targetVoxelNode,sampleNo));
                     
    if ( samplePhysical!=blockedVol )
    {
      // Setup history
      //
      G4NavigationHistory_NewLevel(history, samplePhysical, kNormal);
      
      sampleSolid =
		G4LogicalVolume_GetSolid(
			G4VPhysicalVolume_GetLogicalVolume( samplePhysical ));
			
	  G4AffineTransform tf = G4NavigationHistory_GetTopTransform( history );
			
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
  
  //return G4NormalNavigation_LevelLocate( history, blockedVol, &globalPoint, globalDirection, pLocatedOnEdge, localPoint );
}

// ********************************************************************
// ComputeStep
// ********************************************************************
//
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
			)
{
  /*return G4NormalNavigation_ComputeStep(
		localPoint, localDirection, currentProposedStepLength,
		newSafety, history, validExitNormal, exitNormal,
		exiting, entering, pBlockedPhysical );*/
		
  GEOMETRYLOC G4VPhysicalVolume *motherPhysical, *samplePhysical,
	*blockedExitedVol = GEOMETRYNULL;
  GEOMETRYLOC G4LogicalVolume *motherLogical;
  GEOMETRYLOC G4VSolid *motherSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int sampleNo; // , localNoDaughters;

  G4bool initialNode, noStep;
  GEOMETRYLOC const G4SmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;

  motherPhysical = G4NavigationHistory_GetTopVolume( history );
  motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = G4LogicalVolume_GetSolid(motherLogical);

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

#ifdef USE_BLIST
  G4VoxelNavigation_EnlargeAndResetBlist( This, G4LogicalVolume_GetNoDaughters(motherLogical) );
#endif

  initialNode = true;
  noStep = true;

  while ( noStep )
  {
    curVoxelNode = This->fVoxelNode;
    curNoVolumes = G4VoxelNode_GetNoContained(curVoxelNode);
    for (contentNo=curNoVolumes-1; contentNo>=0; contentNo--)
    {
      sampleNo = G4VoxelNode_GetVolume( curVoxelNode, contentNo);
      
#ifdef USE_BLIST
      if (!This->fBlist[sampleNo])
      {
		This->fBlist[sampleNo] = 1;
#endif
		
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
#ifdef USE_BLIST
          } // -- FBLIST
#endif
        }
      }
    }
	
	
    if (initialNode)
    {	
      initialNode = false;
      voxelSafety = G4VoxelNavigation_ComputeVoxelSafety(This,localPoint);
      if ( voxelSafety<ourSafety )
      {
        ourSafety = voxelSafety;
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
    if (noStep)
    {
      noStep = G4VoxelNavigation_LocateNextVoxel(This, localPoint, localDirection, ourStep);
    }
  }  // end -while (noStep)- loop
  int locationId = get_global_id(0);
	//Result [ locationId ] = ourStep;
  return ourStep;
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
G4VoxelNavigation_ComputeSafety(
			G4VoxelNavigation *This,
			G4ThreeVector localPoint,
			const G4NavigationHistory *history )
{
  GEOMETRYLOC G4VPhysicalVolume *motherPhysical, *samplePhysical;
  GEOMETRYLOC G4LogicalVolume *motherLogical;
  GEOMETRYLOC G4VSolid *motherSolid;
  G4double motherSafety, ourSafety;
  G4int sampleNo; // , localNoDaughters;
  GEOMETRYLOC const G4SmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;

  motherPhysical = G4NavigationHistory_GetTopVolume(history);
  motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = G4LogicalVolume_GetSolid(motherLogical);

  //
  // Compute mother safety
  //

  motherSafety = G4VSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety;                 // Working isotropic safety

  //
  // Compute daughter safeties 
  //

  //localNoDaughters = G4LogicalVolume_GetNoDaughters(motherLogical);

  //  Look only inside the current Voxel only (in the first version).
  //
  curVoxelNode = This->fVoxelNode;
  curNoVolumes = G4VoxelNode_GetNoContained(curVoxelNode);

  for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
  {
    sampleNo = G4VoxelNode_GetVolume(curVoxelNode,contentNo);
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
  voxelSafety = G4VoxelNavigation_ComputeVoxelSafety(This,localPoint);
  if ( voxelSafety<ourSafety )
  {
    ourSafety = voxelSafety;
  }
  return ourSafety;
}*/

#endif
