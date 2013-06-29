
/**
 * G4VoxelLimits inline implementation
 * based on G4VoxelLimits.(i)cc of Geant 4.9.3
 */

#ifndef G4VOXEL_LIMITS_C
#define G4VOXEL_LIMITS_C

#include "G4BuildVoxels.h"

INLINEFUNC
G4double G4VoxelLimits_GetMaxXExtent( const G4VoxelLimits *This )
{
  return This->fxAxisMax;
}

INLINEFUNC
G4double G4VoxelLimits_GetMaxYExtent( const G4VoxelLimits *This )
{
  return This->fyAxisMax;
}

INLINEFUNC
G4double G4VoxelLimits_GetMaxZExtent( const G4VoxelLimits *This )
{
  return This->fzAxisMax;
}

INLINEFUNC
G4double G4VoxelLimits_GetMinXExtent( const G4VoxelLimits *This )
{
  return This->fxAxisMin;
}

INLINEFUNC
G4double G4VoxelLimits_GetMinYExtent( const G4VoxelLimits *This )
{
  return This->fyAxisMin;
}

INLINEFUNC
G4double G4VoxelLimits_GetMinZExtent( const G4VoxelLimits *This )
{
  return This->fzAxisMin;
}

INLINEFUNC
G4double G4VoxelLimits_GetMaxExtent(const G4VoxelLimits *This, EAxis pAxis)
{
  if (pAxis==kXAxis)
  {
    return G4VoxelLimits_GetMaxXExtent(This);
  }
  else if (pAxis==kYAxis)
  {
    return G4VoxelLimits_GetMaxYExtent(This);
  }
  else 
  {
    assert(pAxis==kZAxis);
    return G4VoxelLimits_GetMaxZExtent(This);
  }
}

INLINEFUNC
G4double G4VoxelLimits_GetMinExtent(const G4VoxelLimits *This, EAxis pAxis)
{
  if (pAxis==kXAxis)
  {
    return G4VoxelLimits_GetMinXExtent(This);
  }
  else if (pAxis==kYAxis)
  {
    return G4VoxelLimits_GetMinYExtent(This);
  }
  else 
  {
    assert(pAxis==kZAxis);
    return G4VoxelLimits_GetMinZExtent(This);
  }
}

INLINEFUNC
G4bool G4VoxelLimits_IsXLimited( const G4VoxelLimits *This )
{
  return (This->fxAxisMin==-kInfinity&&This->fxAxisMax==kInfinity) ? false : true;
}

INLINEFUNC
G4bool G4VoxelLimits_IsYLimited( const G4VoxelLimits *This )
{
  return (This->fyAxisMin==-kInfinity&&This->fyAxisMax==kInfinity) ? false : true;
}

INLINEFUNC
G4bool G4VoxelLimits_IsZLimited( const G4VoxelLimits *This )
{
  return (This->fzAxisMin==-kInfinity&&This->fzAxisMax==kInfinity) ? false : true;
}

INLINEFUNC
G4bool G4VoxelLimits_IsLimited( const G4VoxelLimits *This )
{
  return (G4VoxelLimits_IsXLimited(This)||G4VoxelLimits_IsYLimited(This)||G4VoxelLimits_IsZLimited(This));
}

INLINEFUNC
G4bool G4VoxelLimits_IsLimited_axis(const G4VoxelLimits *This, EAxis pAxis)
{
  if (pAxis==kXAxis)
  {
    return G4VoxelLimits_IsXLimited(This);
  }
  else if (pAxis==kYAxis)
  {
    return G4VoxelLimits_IsYLimited(This);
  }
  else 
  {
    assert(pAxis==kZAxis);
    return G4VoxelLimits_IsZLimited(This);
  }
}

INLINEFUNC
G4bool G4VoxelLimits_Inside(const G4VoxelLimits *This, G4ThreeVector pVec)
{
  return ((G4VoxelLimits_GetMinXExtent(This)<=pVec.x) &&
	  (G4VoxelLimits_GetMaxXExtent(This)>=pVec.x) &&
	  (G4VoxelLimits_GetMinYExtent(This)<=pVec.y) &&
	  (G4VoxelLimits_GetMaxYExtent(This)>=pVec.y) &&
	  (G4VoxelLimits_GetMinZExtent(This)<=pVec.z) &&
	  (G4VoxelLimits_GetMaxZExtent(This)>=pVec.z) ) ? true : false;
}

///////////////////////////////////////////////////////////////////////////
//
// Empty constructor and destructor
//

MAYINLINE
void G4VoxelLimits_ctor( G4VoxelLimits *This )
{
	This->fxAxisMin = -kInfinity;
	This->fyAxisMin = -kInfinity;
	This->fzAxisMin = -kInfinity;
	This->fxAxisMax = kInfinity;
	This->fyAxisMax = kInfinity;
	This->fzAxisMax = kInfinity;
}


///////////////////////////////////////////////////////////////////////////
//
// Further restrict limits
// No checks for illegal restrictions
//

MAYINLINE
void G4VoxelLimits_AddLimit(
			G4VoxelLimits *This,
			const EAxis pAxis, 
			const G4double pMin,
			const G4double pMax )
{
  if ( pAxis == kXAxis )
  {
    if ( pMin > This->fxAxisMin ) This->fxAxisMin = pMin ;    
    if ( pMax < This->fxAxisMax ) This->fxAxisMax = pMax ;    
  }
  else if ( pAxis == kYAxis )
  {
    if ( pMin > This->fyAxisMin ) This->fyAxisMin = pMin ;    
    if ( pMax < This->fyAxisMax ) This->fyAxisMax = pMax ;
  }
  else
  { 
    assert( pAxis == kZAxis ) ;

    if ( pMin > This->fzAxisMin ) This->fzAxisMin = pMin ;
    if ( pMax < This->fzAxisMax ) This->fzAxisMax = pMax ;
  }
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate the `outcode' for the specified vector:
// The following bits are set:
//   0      pVec.x()<fxAxisMin && IsXLimited()
//   1      pVec.x()>fxAxisMax && IsXLimited()
//   2      pVec.y()<fyAxisMin && IsYLimited()
//   3      pVec.y()>fyAxisMax && IsYLimited()
//   4      pVec.z()<fzAxisMin && IsZLimited()
//   5      pVec.z()>fzAxisMax && IsZLimited()
//

MAYINLINE
G4int G4VoxelLimits_OutCode( const G4VoxelLimits *This, G4ThreeVector pVec )
{
  G4int code = 0 ;                // The outcode

  if ( G4VoxelLimits_IsXLimited(This) )
  {
    if ( pVec.x < This->fxAxisMin ) code |= 0x01 ;
    if ( pVec.x > This->fxAxisMax ) code |= 0x02 ;
  }
  if ( G4VoxelLimits_IsYLimited(This) )
  {
    if ( pVec.y < This->fyAxisMin ) code |= 0x04 ;
    if ( pVec.y > This->fyAxisMax ) code |= 0x08 ;
  }
  if (G4VoxelLimits_IsZLimited(This))
  {
    if ( pVec.z < This->fzAxisMin ) code |= 0x10 ;
    if ( pVec.z > This->fzAxisMax ) code |= 0x20 ;
  }
  return code;
}


///////////////////////////////////////////////////////////////////////////
//
// ClipToLimits
//
// Clip the line segment pStart->pEnd to the volume described by the
// current limits. Return true if the line remains after clipping,
// else false, and leave the vectors in an undefined state.
//
// Process:
//
// Use Cohen-Sutherland clipping in 3D
// [Fundamentals of Interactive Computer Graphics,Foley & Van Dam]
//

MAYINLINE
G4bool G4VoxelLimits_ClipToLimits(
		const G4VoxelLimits *This,
		G4ThreeVector* pStart,
		G4ThreeVector* pEnd )
{
  G4int sCode, eCode ;
  G4bool remainsAfterClip ;
    
  // Determine if line is trivially inside (both outcodes==0) or outside
  // (logical AND of outcodes !=0)

  sCode = G4VoxelLimits_OutCode(This, *pStart);
  eCode = G4VoxelLimits_OutCode(This, *pEnd);

  if ( sCode & eCode )
  {
    // Trivially outside, no intersection with region

    remainsAfterClip = false;
  }
  else if ( sCode == 0 && eCode == 0 )
  {
    // Trivially inside, no intersections

    remainsAfterClip = true ;
  }
  else
  {
    // Line segment *may* cut volume boundaries
    // At most, one end point is inside

    G4double x1, y1, z1, x2, y2, z2 ;

    x1 = pStart->x ;
    y1 = pStart->y ;
    z1 = pStart->z ;

    x2 = pEnd->x ;
    y2 = pEnd->y ;
    z2 = pEnd->z ;
     
    while ( sCode != eCode )
    {
      // Copy vectors to work variables x1-z1,x2-z2
      // Ensure x1-z1 lies outside volume, swapping vectors and outcodes
      // if necessary

      if ( sCode )
      {
        if ( sCode & 0x01 )  // Clip against fxAxisMin
        {
          z1 += (This->fxAxisMin-x1)*(z2-z1)/(x2-x1);
          y1 += (This->fxAxisMin-x1)*(y2-y1)/(x2-x1);
          x1  = This->fxAxisMin;
        }
        else if ( sCode & 0x02 ) // Clip against fxAxisMax
        {
          z1 += (This->fxAxisMax-x1)*(z2-z1)/(x2-x1);
          y1 += (This->fxAxisMax-x1)*(y2-y1)/(x2-x1);
          x1  = This->fxAxisMax ;
        }
        else if ( sCode & 0x04 )  // Clip against fyAxisMin
        {
          x1 += (This->fyAxisMin-y1)*(x2-x1)/(y2-y1);
          z1 += (This->fyAxisMin-y1)*(z2-z1)/(y2-y1);
          y1  = This->fyAxisMin;
        }
        else if ( sCode & 0x08 )  // Clip against fyAxisMax
        {
          x1 += (This->fyAxisMax-y1)*(x2-x1)/(y2-y1);
          z1 += (This->fyAxisMax-y1)*(z2-z1)/(y2-y1);
          y1  = This->fyAxisMax;
        }
        else if ( sCode & 0x10 )  // Clip against fzAxisMin
        {
          x1 += (This->fzAxisMin-z1)*(x2-x1)/(z2-z1);
          y1 += (This->fzAxisMin-z1)*(y2-y1)/(z2-z1);
          z1  = This->fzAxisMin;
        }
        else if ( sCode & 0x20 )  // Clip against fzAxisMax
        {
          x1 += (This->fzAxisMax-z1)*(x2-x1)/(z2-z1);
          y1 += (This->fzAxisMax-z1)*(y2-y1)/(z2-z1);
          z1  = This->fzAxisMax;
        }
      }
      if ( eCode )  // Clip 2nd end: repeat of 1st, but 1<>2
      {
        if ( eCode & 0x01 )  // Clip against fxAxisMin
        {
          z2 += (This->fxAxisMin-x2)*(z1-z2)/(x1-x2);
          y2 += (This->fxAxisMin-x2)*(y1-y2)/(x1-x2);
          x2  = This->fxAxisMin;
        }
        else if ( eCode & 0x02 )  // Clip against fxAxisMax
        {
          z2 += (This->fxAxisMax-x2)*(z1-z2)/(x1-x2);
          y2 += (This->fxAxisMax-x2)*(y1-y2)/(x1-x2);
          x2  = This->fxAxisMax;
        }
        else if ( eCode & 0x04 )  // Clip against fyAxisMin
        {
          x2 += (This->fyAxisMin-y2)*(x1-x2)/(y1-y2);
          z2 += (This->fyAxisMin-y2)*(z1-z2)/(y1-y2);
          y2  = This->fyAxisMin;
        }
        else if (eCode&0x08)  // Clip against fyAxisMax
        {
          x2 += (This->fyAxisMax-y2)*(x1-x2)/(y1-y2);
          z2 += (This->fyAxisMax-y2)*(z1-z2)/(y1-y2);
          y2  = This->fyAxisMax;
        }
        else if ( eCode & 0x10 )  // Clip against fzAxisMin
        {
          x2 += (This->fzAxisMin-z2)*(x1-x2)/(z1-z2);
          y2 += (This->fzAxisMin-z2)*(y1-y2)/(z1-z2);
          z2  = This->fzAxisMin;
        }
        else if ( eCode & 0x20 )  // Clip against fzAxisMax
        {
          x2 += (This->fzAxisMax-z2)*(x1-x2)/(z1-z2);
          y2 += (This->fzAxisMax-z2)*(y1-y2)/(z1-z2);
          z2  = This->fzAxisMax;
        }
      }
      //  G4endl; G4cout<<"x1 = "<<x1<<"\t"<<"x2 = "<<x2<<G4endl<<G4endl;
      *pStart = G4ThreeVector_create(x1,y1,z1);
      *pEnd   = G4ThreeVector_create(x2,y2,z2);
      sCode  = G4VoxelLimits_OutCode(This,*pStart);
      eCode  = G4VoxelLimits_OutCode(This,*pEnd);
    }
    if ( sCode == 0 && eCode == 0 ) remainsAfterClip = true;
    else                            remainsAfterClip = false;
  }
  return remainsAfterClip;
}

#endif
