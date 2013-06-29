
/**
 * G4NewNavigation implementation
 * New navigation is a rewrite of some of the core parts of VoxelNavigation for better efficiency on GPUs. 
 * The algorithm of ComputeStep primarily, has been altered to get a better performance on the GPU.
 *
 * NOTE: NewNavigaion is not to be confused with COMBINED_NAVIGATION which was created by Otto Seiskari. The idea behind that was to use some  things from Voxel Navigation
 * and some things from Normal Navigation and come up with a solution midway.
 * On the other hand, NEW_NAVIGATION is a revamped algorithm for Voxel Navigation which is aimed specifically for GPUs.
 * The rest of the code is based on G4VoxelNavigation.(i)cc of Geant 4.9.3
 */


#include "G4NewNavigation.h"
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
G4NewNavigation_VoxelLocate(
			G4NewNavigation *This,
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
G4NewNavigation_LocateNextVoxel(
			G4NewNavigation *This,
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
G4NewNavigation_ComputeVoxelSafety(
			const G4NewNavigation *This,
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
void G4NewNavigation_ctor( G4NewNavigation *This )
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
void G4NewNavigation_EnlargeAndResetBlist( G4NewNavigation *This, G4int n )
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
G4NewNavigation_LevelLocate(
			G4NewNavigation *This,
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
	G4NewNavigation_VoxelLocate(This,targetVoxelHeader,*localPoint);

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










//------------------------------------ 

// ********************************************************************
// ComputeStep
// ********************************************************************
//


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
			 GEOMETRYLOC G4SmartVoxelNode  *nullVNode,
			 G4bool cur_vol_local
#ifdef CHECK 
			,GEOMETRYLOC float * Result
#endif			
			)
{
	
  GEOMETRYLOC G4VPhysicalVolume *motherPhysical, *samplePhysical,
	*blockedExitedVol = GEOMETRYNULL;
  GEOMETRYLOC G4LogicalVolume *motherLogical;
  GEOMETRYLOC G4VSolid *motherSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int sampleNo; // , localNoDaughters;
    
  
  //EDIT: Defining these here for global scope.
  GEOMETRYLOC  G4VSolid *sampleSolid;
  G4AffineTransform sampleTf;
  const G4ThreeVector samplePoint;
  G4double sampleSafety;

  int PrevSum;
	//  PrevSum stores the sum of all solids of one type so far. 
  int Count_of_Solid_type;
	// This integer is to fill up the shared mem array for solid types

// For definition of shared memory arrays see kernel trace in gpu.c 
// _--------------------------------__   

  G4bool initialNode;
  GEOMETRYLOC const G4SmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;
  G4double sampleStep;
  	
#if ( GLOBAL_MODE  == 1)
  int locationId = get_global_id(0);
#else  
  int locationId = get_local_id(0);
#endif

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
  G4NewNavigation_EnlargeAndResetBlist( This, G4LogicalVolume_GetNoDaughters(motherLogical) );
#endif
  
  // noStepLocal -> Local boolean storing the noStep value for the current track.
  // noStepAll -> Value returned after a reduction of all noStepLocals

  initialNode = true;
  G4bool noStepLocal = true;
  G4bool noStepAll = true;
 // GEOMETRYLOC G4SmartVoxelNode * nullVNodePointer;
  //G4SmartVoxelNode nullNode;
  //nullVNodePointer = &nullNode;
  int counter = -1;
  while ( noStepAll )
  {	 counter ++;
        if( !noStepLocal )
        {   
             curVoxelNode = nullVNode;
             curNoVolumes = 0;
        }
		else
		{
		    curVoxelNode = This->fVoxelNode;
			curNoVolumes = G4VoxelNode_GetNoContained( curVoxelNode );
		}
		int number_of_volume_types = Solidcount, Curr_vol_type=0;
			// Volume type is from the enum ESolids in everything.h

		int current_solid_sum = 0;
		// if (GLOBAL_MODE ==1)
		//     number_of_threads = get_global_size(0);
		int number_of_threads = get_local_size(0);
	 
		// STEP 1:  Iteration through all volume types. Min is calculated one type at a time.
		//REMEMBER:EDIT: Changed starting value of Curr_vol_type to 1	
		for (Curr_vol_type = 1; Curr_vol_type < number_of_volume_types; (Curr_vol_type)++ )
		{	
			
			current_solid_sum = 0;
			int number_of_solids;
			number_of_solids = curVoxelNode->SolidType[ Curr_vol_type ];
				// The number of solids of this type in the current voxel for the current track
				// NOTE: There may be an overhead here. If more than one thread is in the same voxel and reads this data, mem access not coalesced
		
            Numbers_Of_Solid[ locationId ] = number_of_solids;
                // At this stage a parallel scan sum has to be called which updates the sums for all the solids. The function is inlined.
			
			
		    BARRIER_FLEXIBLE;
		        // Barrier to ensure that all threads have filled up the shared mem array. See everything.h for definition
			
			//int all_threads = get_global_size(1);
			
			Prefix_Sum( Numbers_Of_Solid, Sum_Of_Solids, number_of_threads );
				// The array Sum_Of_Solids stores the final result after the Prefix sum scan. 
			
				//Result[PrevSum + current_solid_sum ] = locationId;

			PrevSum = Sum_Of_Solids[ locationId ];
			
			
			
			// NOTE: The use of contentNo below is to maintain some resemblance to Otto's definition in VoxelNavigation.c 
			if ( noStepLocal )
			{
			  for( contentNo = curNoVolumes-1; contentNo>=0 ; contentNo--)
			  {	
				sampleNo = G4VoxelNode_GetVolume(curVoxelNode,contentNo);
				samplePhysical = G4LogicalVolume_GetDaughter(motherLogical,sampleNo);
				//if ( samplePhysical != blockedExitedVol )
				{	
					// NOTE: blockedExitedVol check makes sense for the serial version. Does not make as much sense to keep it on in the parallel version as well.
					sampleSolid = G4LogicalVolume_GetSolid( G4VPhysicalVolume_GetLogicalVolume( samplePhysical ));
					
					// NOTE: We iterate over all solids and compare if the solid type is what we are looking for. If it is then it is added
					// to the array Solids in shared mem. However, a better implementation is to have a way such that the solids are returned in sorted order 
					// the first place.

					if( sampleSolid->type == Curr_vol_type)
					{	
						// Should the solid be stored or the Physical Volume??
						SolidInfo Info = { samplePhysical, locationId };
						Solids[ PrevSum + current_solid_sum ] = Info;
						
						// PrevSum is the sum of all solids of that type for all threads upto this element not including the current thread.
						// current_solid_sum at the end of all iterations should be equal to the number_of_solids.
		      			current_solid_sum++;
					}
				  }
			   } 
			}
			//EDIT: Change to BARRIER_LOCAL when not testing GLOBAL_MODE 
            BARRIER_FLEXIBLE;
			
			//MODIFY : Perhaps a check is in order here.
			// One check at this point could be whether current_solid_sum is equal to the number_of_solids.
			// Before proceeding to call the kernel for calculating the min, check if code works up to this point
			
			
			
			int Total_solids_of_this_type;
			Total_solids_of_this_type = Numbers_Of_Solid[ number_of_threads - 1 ] + Sum_Of_Solids[ number_of_threads - 1 ];
			
			if ( Total_solids_of_this_type > BlockSize * Multiplier)
		   {
			   Total_solids_of_this_type = BlockSize * Multiplier;
			   // Update Errors or return error from here?
			}
		   
			
				// Checking a candidate solid of current type.
				int k=0;
				int iterations = (Total_solids_of_this_type / number_of_threads);
				for( k = 0; k < iterations ; k++ )
				{	
					int Work_id = locationId + number_of_threads*k;
					if( ( Work_id ) < Total_solids_of_this_type)
					{	
						GEOMETRYLOC G4VPhysicalVolume * candPhysical = ( Solids[ Work_id  ].PVolume);
						GEOMETRYLOC G4VSolid * candSolid = ((candPhysical)->flogical)->fSolid;
						int candId = Solids[ Work_id].trackId;
						G4ThreeVector candGlobalPoint = LocationArray[candId].Point;
						G4ThreeVector candGlobalDirection = LocationArray[candId].Direction;
						G4AffineTransform candTf = G4AffineTransform_create_full(
					                                G4VPhysicalVolume_GetObjectRotationValue(candPhysical),
													G4VPhysicalVolume_GetTranslation(candPhysical));

						G4AffineTransform_Invert( &candTf );
				     	const G4ThreeVector candPoint=
						           G4AffineTransform_TransformPoint(&candTf,candGlobalPoint);
						const G4double candSafety = G4VSolid_DistanceToIn( candSolid, candPoint );
						const G4ThreeVector candDirection = G4AffineTransform_TransformAxis( &candTf, candGlobalDirection );
						const G4double candStep = G4VSolid_DistanceToIn_full( candSolid, candPoint, candDirection );
							// Result_For_Current_Solid should hold data for the sampleSafety and the sampleStep
						// NOTE: In this version of navigation the step is calculated along with the safety. It is not checked if
						// step<safety; the safety is only calculated here because the physics may require it later.
						ResultInfo Result_of_Solid = { candSafety, candStep, candId, candPhysical } ;
						Result_For_Current_Solid[ Work_id ] = Result_of_Solid ;					
					}
				    BARRIER_FLEXIBLE;
					// REMOVE: Barrier used in debugging stage. Remove when done.
				}
			
			    	
				// At this point all the safeties and steps have been calculated, now to find the minimum step for the current solid per track.
				// find_minimum...
				// Most basic implementation for finding minimum, a proper min finding algorithm that takes into account threadIds is hard to find.
				
				Find_minimum ( Result_For_Current_Solid, Compacter_Result, PrevSum, number_of_solids );
					// Minimum finding function that finds the minimum step per track in Result_For_Current_Solid and
					// stores the minimum step, the safety and the pointer to the sampleSolid ( physical? ) 
					// The minimum is basically compared to the existing value and stored if smaller. This way find_minimum does not
					// care about which solid type is currently being processed.
					// See gpu.c where Compacter Results initial step value is set to kInfinfity.
					//if( Curr_vol_type ==1)	

		}
		
			//if( Curr_vol_type == 1)
			
	    if (  Compacter_Result[ locationId ].step <= ourStep )	
       {
	     ourStep = Compacter_Result[ locationId ].step;
	 	 *pBlockedPhysical = Compacter_Result[ locationId ].PVolume;
	     *entering = true;
		 *exiting = false;
		 
		}
		
	    if (initialNode)
	    {	
			initialNode = false;
			voxelSafety = G4NewNavigation_ComputeVoxelSafety(This,localPoint);
			if ( voxelSafety<ourSafety )
			{
				ourSafety = voxelSafety;
			}

		    if ( currentProposedStepLength<ourSafety )
			{
				// Guaranteed physics limited
				//      
				noStepLocal = false;
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
	
		if (noStepLocal)
		{
			noStepLocal = G4NewNavigation_LocateNextVoxel(This, localPoint, localDirection, ourStep);
		}

		noStepArray[ locationId ] = noStepLocal;
		BARRIER_FLEXIBLE;

		noStepAll = NoStepReduction( noStepArray, number_of_threads);
		
		//Prefix_Sum ( noStepArray, noStepArray, number_of_threads);
		
 }  // end -while (noStep)- loop

 // Double check to make sure all threads have reached here before exiting

BARRIER_FLEXIBLE;

	//Result[ locationId ] = ourStep;	
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
G4NewNavigation_ComputeSafety(
			G4NewNavigation *This,
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
  voxelSafety = G4NewNavigation_ComputeVoxelSafety(This,localPoint);
  if ( voxelSafety<ourSafety )
  {
    ourSafety = voxelSafety;
  }
  return ourSafety;
}*/

#endif
