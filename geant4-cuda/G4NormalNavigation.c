
/**
 * G4NormalNavigation implementation
 * based on G4NormalNavigation.* of Geant 4.9.3
 */

#include "G4Navigator.h"
#include "G4AffineTransform.h"
#include "G4LogicalVolume.h"

INLINEFUNC G4bool
G4AuxiliaryNavServices_CheckPointOnSurface(
				     GEOMETRYLOC const G4VSolid* sampleSolid, 
                     G4ThreeVector localPoint, 
                     const G4ThreeVector* globalDirection, 
                     G4AffineTransform sampleTransform,
                     const G4bool locatedOnEdge)
{
  G4ThreeVector localDirection, sampleNormal;
  G4bool        enter = false;

  EInside insideSolid = 
	G4VSolid_Inside(sampleSolid, localPoint);
	
  if ( insideSolid!=kOutside ) 
  {
    G4bool checkDirection= locatedOnEdge && (globalDirection!=0);
    if( (insideSolid==kSurface) && checkDirection)
    {
      // We are probably located on an edge.
      //
      localDirection= G4AffineTransform_TransformAxis(&sampleTransform,*globalDirection); 

      // Check whether we enter the volume
      // 
      sampleNormal = G4VSolid_SurfaceNormal(sampleSolid,localPoint);
      if ( G4ThreeVector_dot(sampleNormal,localDirection) <= 0 )
      {
        if( G4ThreeVector_dot(sampleNormal,localDirection) == 0 )
        {
          // We can't decide yet, let's make sure we're entering the solid.
          // If by a confusion we entered the next solid we find out now
          // whether to leave or to enter.
          // This happens when we're on the surface or edge shared by two
          // solids
          //
          
          G4double distanceToIn =
			G4VSolid_DistanceToIn_full( sampleSolid, localPoint, localDirection );
			
          if( distanceToIn != kInfinity )
          {
            enter = true;
          } 
        }
        else
        {
          enter = true;
        }
      }
    }
    else
    {
      enter = true;
    }
  }
  return enter;
}

// --------------------------------------------------------------------

/*INLINEFUNC G4bool
G4AuxiliaryNavServices_CheckPointExiting(
				   GEOMETRYLOC const G4VSolid* sampleSolid, 
                   G4ThreeVector localPoint, 
                   const G4ThreeVector* globalDirection, 
                   G4AffineTransform sampleTransform )
{
  if( !globalDirection )  { return false; }

  G4ThreeVector localDirection, sampleNormal;
  G4bool        exiting = false;

  EInside insideSolid =
	G4VSolid_Inside(sampleSolid, localPoint);
  
  if( (insideSolid==kSurface) )
  {
    localDirection= G4AffineTransform_TransformAxis(&sampleTransform,*globalDirection); 

    // Check whether we are exiting the volume
    // 
    sampleNormal = G4VSolid_SurfaceNormal(sampleSolid,localPoint);
    
    if ( G4ThreeVector_dot(sampleNormal,localDirection) >= 0 )
    {
      if( G4ThreeVector_dot(sampleNormal,localDirection) == 0 )
      {
        // We can't decide yet, let's make sure we're entering the solid.
        // If by a confusion we entered the next solid we find out now
        // whether to leave or to exiting.
        // This happens when we're on the surface or edge shared by two
        // solids
        //
        
        G4double distanceToIn =
			G4VSolid_DistanceToIn_full( sampleSolid, localPoint, localDirection );
                 
        if( distanceToIn != kInfinity )
        {
          exiting = true;
        } 
      }
      else
      {
        exiting = true;
      }
    }
  }
  return exiting;
}*/

#ifndef ENABLE_COMBINED_NAVIGATION

// --------------------------------------------------------------------

MAYINLINE G4bool
G4NormalNavigation_LevelLocate(
	G4NavigationHistory *history,
	GEOMETRYLOC const G4VPhysicalVolume *blockedVol,
	G4ThreeVector* globalPoint,
	const G4ThreeVector* globalDirection,
	G4bool pLocatedOnEdge, 
	G4ThreeVector* localPoint )
{
  GEOMETRYLOC G4VPhysicalVolume *targetPhysical, *samplePhysical;
  GEOMETRYLOC G4LogicalVolume *targetLogical;
  GEOMETRYLOC G4VSolid *sampleSolid;
  G4ThreeVector samplePoint;
  G4int targetNoDaughters;
  
  targetPhysical = G4NavigationHistory_GetTopVolume(history);
  targetLogical = G4VPhysicalVolume_GetLogicalVolume(targetPhysical);
  targetNoDaughters = G4LogicalVolume_GetNoDaughters(targetLogical);
  
  if (targetNoDaughters == 0) return false;
    //
	// Search daughters in volume
	//
  for ( int sampleNo=targetNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
	  samplePhysical =
		G4LogicalVolume_GetDaughter(targetLogical,sampleNo);

	  if ( samplePhysical!=blockedVol )
	  {
		// Setup history
		//
		G4NavigationHistory_NewLevel(history, samplePhysical, kNormal );
		
		sampleSolid =
			G4LogicalVolume_GetSolid(
				G4VPhysicalVolume_GetLogicalVolume(samplePhysical));
				
		G4AffineTransform tf =
			G4NavigationHistory_GetTopTransform(history);
		
		samplePoint =
			G4AffineTransform_TransformPoint( &tf, *globalPoint );
			
		if( G4AuxiliaryNavServices_CheckPointOnSurface(
			sampleSolid, samplePoint, globalDirection, 
			tf, pLocatedOnEdge) )
		{
		  // Enter This daughter
		  //
		  *localPoint = samplePoint;
		  return true;
		}
		else
		{
			G4NavigationHistory_BackLevel(history);
		}
	  }
  }
  
  return false;
}


// ********************************************************************
// ComputeStep
// ********************************************************************
//
MAYINLINE 
G4double
G4NormalNavigation_ComputeStep(
	G4ThreeVector localPoint,
	G4ThreeVector localDirection,
	const G4double currentProposedStepLength,
	G4double *newSafety,
	G4NavigationHistory *history,
	G4bool *validExitNormal,
	G4ThreeVector *exitNormal,
	G4bool *exiting,
	G4bool *entering,
	GEOMETRYLOC G4VPhysicalVolume *(*pBlockedPhysical))
{
  GEOMETRYLOC G4VPhysicalVolume *motherPhysical, *samplePhysical, *blockedExitedVol=0;
  GEOMETRYLOC G4LogicalVolume *motherLogical;
  GEOMETRYLOC G4VSolid *motherSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int localNoDaughters, sampleNo;

  motherPhysical = G4NavigationHistory_GetTopVolume(history);
  
  motherLogical  = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid    = G4LogicalVolume_GetSolid(motherLogical);

  // Compute mother safety
  //
 
  motherSafety = G4VSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety; // Working isotropic safety

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
      blockedExitedVol =* pBlockedPhysical;
      ourSafety = 0;
    }
  }
  *exiting  = false;
  *entering = false;

  localNoDaughters = G4LogicalVolume_GetNoDaughters(motherLogical);
  
  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo--)
  {
    samplePhysical = G4LogicalVolume_GetDaughter(motherLogical,sampleNo);
    
    if ( samplePhysical!=blockedExitedVol )
    {
      G4AffineTransform sampleTf =
       G4AffineTransform_create_full(
		G4VPhysicalVolume_GetObjectRotationValue(samplePhysical),
		G4VPhysicalVolume_GetTranslation(samplePhysical));
		
	  G4AffineTransform_Invert(&sampleTf);
	  
      const G4ThreeVector samplePoint =
			G4AffineTransform_TransformPoint(&sampleTf, localPoint);
			
      GEOMETRYLOC const G4VSolid *sampleSolid =
		G4LogicalVolume_GetSolid(
			G4VPhysicalVolume_GetLogicalVolume( samplePhysical ));
              
      const G4double sampleSafety =
		G4VSolid_DistanceToIn(sampleSolid,samplePoint);
             
      if ( sampleSafety<ourSafety )
      {
        ourSafety=sampleSafety;
      }
      if ( sampleSafety<=ourStep )
      {
		
        sampleDirection = G4AffineTransform_TransformAxis(&sampleTf, localDirection);
        
        const G4double sampleStep =
			G4VSolid_DistanceToIn_full(sampleSolid,samplePoint,sampleDirection);

        if ( sampleStep<=ourStep )
        {
          ourStep  = sampleStep;
          *entering = true;
          *exiting  = false;
          *pBlockedPhysical = samplePhysical;
        }
      }
    }
  }
  
  if ( currentProposedStepLength<ourSafety )
  {
    // Guaranteed physics limited
    //
    *entering = false;
    *exiting  = false;
    *pBlockedPhysical = GEOMETRYNULL;
    ourStep = kInfinity;
  }
  else
  {
    // Compute mother intersection if required
    //
    if ( motherSafety<=ourStep )
    {
      G4double motherStep =
		G4VSolid_DistanceToOut_full(
			motherSolid,
			localPoint,
			localDirection,
			true,
			validExitNormal,
			exitNormal);

      if ( motherStep<=ourStep )
      {
        ourStep  = motherStep;
        *exiting  = true;
        *entering = false;
        if ( *validExitNormal )
        {
          G4RotationMatrix rot = G4VPhysicalVolume_GetObjectRotationValue(motherPhysical);
		  G4RotationMatrix inv = G4RotationMatrix_inverse(&rot);
          *exitNormal = G4RotationMatrix_apply(&inv, *exitNormal);
        }
      }
      else
      {
        *validExitNormal = false;
      }
    }
  }
  *newSafety = ourSafety;
  return ourStep;
}

// ********************************************************************
// ComputeSafety
// ********************************************************************
//
/*MAYINLINE 
G4double G4NormalNavigation_ComputeSafety(
	G4ThreeVector localPoint,
	const G4NavigationHistory *history )
{
  GEOMETRYLOC G4VPhysicalVolume *motherPhysical, *samplePhysical;
  GEOMETRYLOC G4LogicalVolume *motherLogical;
  GEOMETRYLOC G4VSolid *motherSolid;
  G4double motherSafety, ourSafety;
  G4int localNoDaughters, sampleNo;
  
  motherPhysical = G4NavigationHistory_GetTopVolume(history);
  
  motherLogical  = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid    = G4LogicalVolume_GetSolid(motherLogical);

  // Compute mother safety
  //
  
  motherSafety = G4VSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety; // Working isotropic safety


  // Compute daughter safeties 
  //
  
  localNoDaughters = G4LogicalVolume_GetNoDaughters(motherLogical);
  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    samplePhysical = G4LogicalVolume_GetDaughter(motherLogical,sampleNo);
    
    G4AffineTransform sampleTf =
       G4AffineTransform_create_full(
		G4VPhysicalVolume_GetObjectRotationValue(samplePhysical),
		G4VPhysicalVolume_GetTranslation(samplePhysical));
		
	G4AffineTransform_Invert(&sampleTf);
	
    const G4ThreeVector samplePoint =
			G4AffineTransform_TransformPoint(&sampleTf, localPoint);
			
    GEOMETRYLOC const G4VSolid *sampleSolid =
		G4LogicalVolume_GetSolid(
			G4VPhysicalVolume_GetLogicalVolume( samplePhysical ));
            
    const G4double sampleSafety =
			G4VSolid_DistanceToIn(sampleSolid,samplePoint);
            
    if ( sampleSafety<ourSafety )
    {
      ourSafety = sampleSafety;
    }

  }
  return ourSafety;
}*/

#endif
