
/** G4Box implementation, based on G4Box.cc of Geant 4.9.3 */

// --------------------------------------------------------------------

#include "G4Box.h"

#ifndef HOST_CODE

//////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

INLINEFUNC
G4ThreeVector G4Box_ApproxSurfaceNormal( GEOMETRYLOC const G4Box *This, G4ThreeVector p )
{
  G4double distx, disty, distz ;
  G4ThreeVector norm ;

  // Calculate distances as if in 1st octant

  distx = fabs(fabs(p.x) - This->fDx) ;
  disty = fabs(fabs(p.y) - This->fDy) ;
  distz = fabs(fabs(p.z) - This->fDz) ;

  if ( distx <= disty )
  {
    if ( distx <= distz )     // Closest to X
    {
      if ( p.x < 0 ) norm = G4ThreeVector_create(-1.0,0,0) ;
      else             norm = G4ThreeVector_create( 1.0,0,0) ;
    }
    else                      // Closest to Z
    {
      if ( p.z < 0 ) norm = G4ThreeVector_create(0,0,-1.0) ;
      else             norm = G4ThreeVector_create(0,0, 1.0) ;
    }
  }
  else
  {
    if ( disty <= distz )      // Closest to Y
    {
      if ( p.y < 0 ) norm = G4ThreeVector_create(0,-1.0,0) ;
      else             norm = G4ThreeVector_create(0, 1.0,0) ;
    }
    else                       // Closest to Z
    {
      if ( p.z < 0 ) norm = G4ThreeVector_create(0,0,-1.0) ;
      else             norm = G4ThreeVector_create(0,0, 1.0) ;
    }
  }
  return norm;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned
SOLIDINLINE
G4ThreeVector G4Box_SurfaceNormal(GEOMETRYLOC const G4Box *This, G4ThreeVector p)
{
  G4double distx, disty, distz ;
  G4ThreeVector norm ;
  
  const G4double kCarTolerance = K_GEOMETRY_CAR_TOLERANCE;

  // Calculate distances as if in 1st octant

  distx = fabs(fabs(p.x) - This->fDx) ;
  disty = fabs(fabs(p.y) - This->fDy) ;
  distz = fabs(fabs(p.z) - This->fDz) ;

  // New code for particle on surface including edges and corners with specific
  // normals

  const G4double delta    = 0.5*kCarTolerance;
  const G4ThreeVector nX  = G4ThreeVector_create( 1.0, 0,0  );
  const G4ThreeVector nmX = G4ThreeVector_create(-1.0, 0,0  );
  const G4ThreeVector nY  = G4ThreeVector_create( 0, 1.0,0  );
  const G4ThreeVector nmY = G4ThreeVector_create( 0,-1.0,0  );
  const G4ThreeVector nZ  = G4ThreeVector_create( 0, 0,  1.0);
  const G4ThreeVector nmZ = G4ThreeVector_create( 0, 0,- 1.0);

  G4ThreeVector
	normX = G4ThreeVector_create(0.,0.,0.),
	normY = G4ThreeVector_create(0.,0.,0.),
	normZ = G4ThreeVector_create(0.,0.,0.);
	
  G4ThreeVector sumnorm = G4ThreeVector_create(0., 0., 0.);
  G4int noSurfaces=0; 

  if (distx <= delta)         // on X/mX surface and around
  {
    noSurfaces ++; 
    if ( p.x >= 0.){        // on +X surface
      normX= nX ;    // G4ThreeVector( 1.0, 0., 0. );
    }else{
      normX= nmX;    // G4ThreeVector(-1.0, 0., 0. ); 
    }
    sumnorm= normX; 
  }

  if (disty <= delta)    // on one of the +Y or -Y surfaces
  {
    noSurfaces ++; 
    if ( p.y >= 0.){        // on +Y surface
      normY= nY;
    }else{
      normY = nmY; 
    }
    G4ThreeVector_sum_assign( &sumnorm, normY ); // sumnorm += normY; 
  }

  if (distz <= delta)    // on one of the +Z or -Z surfaces
  {
    noSurfaces ++; 
    if ( p.z >= 0.){        // on +Z surface
      normZ= nZ;
    }else{
      normZ = nmZ; 
    }
    G4ThreeVector_sum_assign( &sumnorm, normZ ); // sumnorm += normZ; 
  }

  // sumnorm= normX + normY + normZ; 
  const G4double invSqrt2 = 1.0 / sqrt( 2.0); // TODO rsqrt
  const G4double invSqrt3 = 1.0 / sqrt( 3.0); 

  norm= G4ThreeVector_create( 0., 0., 0.); 
  if( noSurfaces > 0 )
  { 
    if( noSurfaces == 1 ){ 
      norm= sumnorm; 
    }else{
      // norm = sumnorm . unit(); 
      if( noSurfaces == 2 ) { 
        // 2 surfaces -> on edge 
        norm = G4ThreeVector_mult(sumnorm, invSqrt2); 
      } else { 
        // 3 surfaces (on corner)
        norm = G4ThreeVector_mult(sumnorm, invSqrt3); 
      }
    }
  }else{
     norm = G4Box_ApproxSurfaceNormal(This, p);
  }
  
  return norm;
}



///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to box from an outside point
// - return kInfinity if no intersection.
//
// ALGORITHM:
//
// Check that if point lies outside x/y/z extent of box, travel is towards
// the box (ie. there is a possibility of an intersection)
//
// Calculate pairs of minimum and maximum distances for x/y/z travel for
// intersection with the box's x/y/z extent.
// If there is a valid intersection, it is given by the maximum min distance
// (ie. distance to satisfy x/y/z intersections) *if* <= minimum max distance
// (ie. distance after which 1+ of x/y/z intersections not satisfied)
//
// NOTE:
//
// `Inside' safe - meaningful answers given if point is inside the exact
// shape.
SOLIDINLINE
G4double G4Box_DistanceToIn_full(GEOMETRYLOC const G4Box *This, G4ThreeVector p,G4ThreeVector v)
{
  G4double safx, safy, safz ;
  G4double smin=0.0, sminy, sminz ; // , sminx ;
  G4double smax=kInfinity, smaxy, smaxz ; // , smaxx ;  // they always > 0
  G4double stmp ;
  G4double sOut=kInfinity, sOuty=kInfinity, sOutz=kInfinity ;
  const G4double kCarTolerance = K_GEOMETRY_CAR_TOLERANCE;

  safx = fabs(p.x) - This->fDx ;     // minimum distance to x surface of shape
  safy = fabs(p.y) - This->fDy ;
  safz = fabs(p.z) - This->fDz ;

  // Will we intersect?
  // If safx/y/z is >-tol/2 the point is outside/on the box's x/y/z extent.
  // If both p.x/y/z and v.x/y/z repectively are both positive/negative,
  // travel is in a direction away from the shape.

  if (    ((p.x*v.x >= 0.0) && safx > -kCarTolerance*0.5) 
       || ((p.y*v.y >= 0.0) && safy > -kCarTolerance*0.5)
       || ((p.z*v.z >= 0.0) && safz > -kCarTolerance*0.5)   ) 
  {
    return kInfinity ;  // travel away or parallel within tolerance
  }

  // Compute min / max distances for x/y/z travel:
  // X Planes

  if ( v.x)
  {
    stmp = 1.0/fabs(v.x) ;

    if (safx >= 0.0)
    {
      smin = safx*stmp ;
      smax = (This->fDx+fabs(p.x))*stmp ;
    }
    else
    {
      if (v.x > 0)  sOut = (This->fDx - p.x)*stmp ;
      if (v.x < 0)  sOut = (This->fDx + p.x)*stmp ;
    }
  }

  // Y Planes

  if ( v.y) 
  {
    stmp = 1.0/fabs(v.y) ;

    if (safy >= 0.0)
    {
      sminy = safy*stmp ;
      smaxy = (This->fDy+fabs(p.y))*stmp ;

      if (sminy > smin) smin=sminy ;
      if (smaxy < smax) smax=smaxy ;

      if (smin >= smax-kCarTolerance*0.5)
      {
        return kInfinity ;  // touch XY corner
      }
    }
    else
    {
      if (v.y > 0)  sOuty = (This->fDy - p.y)*stmp ;
      if (v.y < 0)  sOuty = (This->fDy + p.y)*stmp ;
      if( sOuty < sOut ) sOut = sOuty ;
    }     
  }

  // Z planes

  if ( v.z )
  {
    stmp = 1.0/fabs(v.z) ;

    if ( safz >= 0.0)
    {
      sminz = safz*stmp ;
      smaxz = (This->fDz+fabs(p.z))*stmp ;

      if (sminz > smin) smin = sminz ;
      if (smaxz < smax) smax = smaxz ;

      if (smin >= smax-kCarTolerance*0.5)
      { 
        return kInfinity ;    // touch ZX or ZY corners
      }
    }
    else
    {
      if (v.z > 0)  sOutz = (This->fDz - p.z)*stmp ;
      if (v.z < 0)  sOutz = (This->fDz + p.z)*stmp ;
      if( sOutz < sOut ) sOut = sOutz ;
    }
  }

  if ( sOut <= smin + 0.5*kCarTolerance) // travel over edge
  {
    return kInfinity ;
  }
  if (smin < 0.5*kCarTolerance)  smin = 0.0 ;

  return smin ;
}

//////////////////////////////////////////////////////////////////////////
// 
// Appoximate distance to box.
// Returns largest perpendicular distance to the closest x/y/z sides of
// the box, which is the most fast estimation of the shortest distance to box
// - If inside return 0
SOLIDINLINE
G4double G4Box_DistanceToIn(GEOMETRYLOC const G4Box *This, G4ThreeVector p)
{
  G4double safex, safey, safez, safe = 0.0 ;

  safex = fabs(p.x) - This->fDx ;
  safey = fabs(p.y) - This->fDy ;
  safez = fabs(p.z) - This->fDz ;

  if (safex > safe) safe = safex ;
  if (safey > safe) safe = safey ;
  if (safez > safe) safe = safez ;

  return safe ;
}

/////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface of box from inside
// by calculating distances to box's x/y/z planes.
// Smallest distance is exact distance to exiting.
// - Eliminate one side of each pair by considering direction of v
// - when leaving a surface & v.close, return 0
SOLIDINLINE
G4double G4Box_DistanceToOut_full( GEOMETRYLOC const G4Box *This, G4ThreeVector p,G4ThreeVector v,
                               const G4bool calcNorm,
                                G4bool *validNorm,G4ThreeVector *n)
{
  const G4double kCarTolerance = K_GEOMETRY_CAR_TOLERANCE;
  enum {kBoxUndefined,kPX,kMX,kPY,kMY,kPZ,kMZ} side = kBoxUndefined ;
  G4double pdist,stmp,snxt;

  if (calcNorm) *validNorm = true ; // All normals are valid

  if (v.x > 0)   // X planes  
  {
    pdist = This->fDx - p.x ;

    if (pdist > kCarTolerance*0.5)
    {
      snxt = pdist/v.x ;
      side = kPX ;
    }
    else
    {
      if (calcNorm) *n    = G4ThreeVector_create(1,0,0) ;
      return         snxt = 0 ;
    }
  }
  else if (v.x < 0) 
  {
    pdist = This->fDx + p.x ;

    if (pdist > kCarTolerance*0.5)
    {
      snxt = -pdist/v.x ;
      side = kMX ;
    }
    else
    {
      if (calcNorm) *n   = G4ThreeVector_create(-1,0,0) ;
      return        snxt = 0 ;
    }
  }
  else snxt = kInfinity ;

  if ( v.y > 0 )   // Y planes  
  {
    pdist=This->fDy-p.y;

    if (pdist>kCarTolerance*0.5)
    {
      stmp=pdist/v.y;

      if (stmp<snxt)
      {
        snxt=stmp;
        side=kPY;
      }
    }
    else
    {
      if (calcNorm) *n    = G4ThreeVector_create(0,1,0) ;
      return         snxt = 0 ;
    }
  }
  else if ( v.y < 0 ) 
  {
    pdist = This->fDy + p.y ;

    if (pdist > kCarTolerance*0.5)
    {
      stmp=-pdist/v.y;

      if (stmp<snxt)
      {
        snxt=stmp;
        side=kMY;
      }
    }
    else
    {
      if (calcNorm) *n    = G4ThreeVector_create(0,-1,0) ;
      return         snxt = 0 ;
    }
  }
  if (v.z>0)        // Z planes 
  {
    pdist=This->fDz-p.z;

    if (pdist > kCarTolerance*0.5)
    {
      stmp=pdist/v.z;

      if (stmp < snxt)
      {
        snxt=stmp;
        side=kPZ;
      }
    }
    else
    {
      if (calcNorm) *n    = G4ThreeVector_create(0,0,1) ;
      return         snxt = 0 ;
    }
  }
  else if (v.z<0) 
  {
    pdist = This->fDz + p.z ;

    if (pdist > kCarTolerance*0.5)
    {
      stmp=-pdist/v.z;

      if (stmp < snxt)
      {
        snxt=stmp;
        side=kMZ;
      }
    }
    else
    {
      if (calcNorm) *n    = G4ThreeVector_create(0,0,-1) ;
      return         snxt = 0 ;
    }
  }
  if (calcNorm)
  {      
    switch (side)
    {
      case kPX:
        *n=G4ThreeVector_create(1,0,0);
        break;
      case kMX:
        *n=G4ThreeVector_create(-1,0,0);
        break;
      case kPY:
        *n=G4ThreeVector_create(0,1,0);
        break;
      case kMY:
        *n=G4ThreeVector_create(0,-1,0);
        break;
      case kPZ:
        *n=G4ThreeVector_create(0,0,1);
        break;
      case kMZ:
        *n=G4ThreeVector_create(0,0,-1);
        break;
      default:
        break;
    }
  }
  return snxt;
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - If outside return 0
SOLIDINLINE
G4double G4Box_DistanceToOut( GEOMETRYLOC const G4Box *This, G4ThreeVector p )
{
  G4double safx1,safx2,safy1,safy2,safz1,safz2,safe=0.0;


  safx1 = This->fDx - p.x ;
  safx2 = This->fDx + p.x ;
  safy1 = This->fDy - p.y ;
  safy2 = This->fDy + p.y ;
  safz1 = This->fDz - p.z ;
  safz2 = This->fDz + p.z ;  
  
  // shortest Dist to any boundary now MIN(safx1,safx2,safy1..)

  if (safx2 < safx1) safe = safx2 ;
  else               safe = safx1 ;
  if (safy1 < safe)  safe = safy1 ;
  if (safy2 < safe)  safe = safy2 ;
  if (safz1 < safe)  safe = safz1 ;
  if (safz2 < safe)  safe = safz2 ;

  if (safe < 0) safe = 0 ;
  return safe ;  
}

#endif // not defined host code

/////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

SOLIDINLINE
EInside G4Box_Inside(GEOMETRYLOC const G4Box *This, G4ThreeVector p)
{
  const G4double kCarTolerance = K_GEOMETRY_CAR_TOLERANCE;
  EInside in = kOutside ;

  if ( fabs(p.x) <= This->fDx - kCarTolerance*0.5 )
  {
    if (fabs(p.y) <= This->fDy - kCarTolerance*0.5 )
    {
      if      (fabs(p.z) <= This->fDz - kCarTolerance*0.5 ) in = kInside ;
      else if (fabs(p.z) <= This->fDz + kCarTolerance*0.5 ) in = kSurface ;
    }
    else if (fabs(p.y) <= This->fDy + kCarTolerance*0.5 )
    {
      if (fabs(p.z) <= This->fDz + kCarTolerance*0.5 ) in = kSurface ;
    }
  }
  else if (fabs(p.x) <= This->fDx + kCarTolerance*0.5 )
  {
    if (fabs(p.y) <= This->fDy + kCarTolerance*0.5 )
    {
      if (fabs(p.z) <= This->fDz + kCarTolerance*0.5) in = kSurface ;
    }
  }
  return in ;
}

#ifdef HOST_CODE

SOLIDINLINE
void G4Box_ctor(G4Box *This,G4double pX,G4double pY,G4double pZ)
{
	This->solid.type = kBox;
	
	This->fDx = pX;
	This->fDy = pY;
	This->fDz = pZ;
}


//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

#ifdef ENABLE_VOXEL_NAVIGATION

#include "G4AffineTransform.h"
#include "G4VoxelLimits.c"
#include "G4VSolid.c"

INLINEFUNC 
G4ThreeVectorList* G4Box_CreateRotatedVertices(const G4Box *This, G4AffineTransform pTransform)
{
	G4ThreeVectorList* vertices = new G4ThreeVectorList;
	vertices->reserve(8);

	G4ThreeVector vertex0 = G4ThreeVector_create(-This->fDx,-This->fDy,-This->fDz) ;
	G4ThreeVector vertex1 = G4ThreeVector_create(This->fDx,-This->fDy,-This->fDz) ;
	G4ThreeVector vertex2 = G4ThreeVector_create(This->fDx,This->fDy,-This->fDz) ;
	G4ThreeVector vertex3 = G4ThreeVector_create(-This->fDx,This->fDy,-This->fDz) ;
	G4ThreeVector vertex4 = G4ThreeVector_create(-This->fDx,-This->fDy,This->fDz) ;
	G4ThreeVector vertex5 = G4ThreeVector_create(This->fDx,-This->fDy,This->fDz) ;
	G4ThreeVector vertex6 = G4ThreeVector_create(This->fDx,This->fDy,This->fDz) ;
	G4ThreeVector vertex7 = G4ThreeVector_create(-This->fDx,This->fDy,This->fDz) ;

	vertices->push_back(G4AffineTransform_TransformPoint(&pTransform, vertex0));
	vertices->push_back(G4AffineTransform_TransformPoint(&pTransform, vertex1));
	vertices->push_back(G4AffineTransform_TransformPoint(&pTransform, vertex2));
	vertices->push_back(G4AffineTransform_TransformPoint(&pTransform, vertex3));
	vertices->push_back(G4AffineTransform_TransformPoint(&pTransform, vertex4));
	vertices->push_back(G4AffineTransform_TransformPoint(&pTransform, vertex5));
	vertices->push_back(G4AffineTransform_TransformPoint(&pTransform, vertex6));
	vertices->push_back(G4AffineTransform_TransformPoint(&pTransform, vertex7));
	return vertices;
}

SOLIDINLINE
G4bool G4Box_CalculateExtent(
			   const G4Box *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax)
{
  const G4double kCarTolerance = K_GEOMETRY_CAR_TOLERANCE;
    
  if (!G4AffineTransform_IsRotated(&pTransform))
  {
    // Special case handling for unrotated boxes
    // Compute x/y/z mins and maxs respecting limits, with early returns
    // if outside limits. Then switch() on pAxis

    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;

    xoffset = G4AffineTransform_NetTranslation(&pTransform).x;
    xMin    = xoffset - This->fDx ;
    xMax    = xoffset + This->fDx ;

    if (G4VoxelLimits_IsXLimited(&pVoxelLimit))
    {
      if ( xMin > G4VoxelLimits_GetMaxXExtent(&pVoxelLimit)+kCarTolerance || 
           xMax <G4VoxelLimits_GetMinXExtent(&pVoxelLimit)-kCarTolerance    ) return false ;
      else
      {
        if (xMin < G4VoxelLimits_GetMinXExtent(&pVoxelLimit))
        {
          xMin = G4VoxelLimits_GetMinXExtent(&pVoxelLimit) ;
        }
        if (xMax > G4VoxelLimits_GetMaxXExtent(&pVoxelLimit))
        {
          xMax = G4VoxelLimits_GetMaxXExtent(&pVoxelLimit) ;
        }
      }
    }
    yoffset = G4AffineTransform_NetTranslation(&pTransform).y ;
    yMin    = yoffset - This->fDy ;
    yMax    = yoffset + This->fDy ;

    if (G4VoxelLimits_IsYLimited(&pVoxelLimit))
    {
      if ( yMin > G4VoxelLimits_GetMaxYExtent(&pVoxelLimit)+kCarTolerance ||
           yMax < G4VoxelLimits_GetMinYExtent(&pVoxelLimit)-kCarTolerance   ) return false ;
      else
      {
        if (yMin < G4VoxelLimits_GetMinYExtent(&pVoxelLimit))
        {
          yMin = G4VoxelLimits_GetMinYExtent(&pVoxelLimit) ;
        }
        if (yMax > G4VoxelLimits_GetMaxYExtent(&pVoxelLimit))
        {
          yMax = G4VoxelLimits_GetMaxYExtent(&pVoxelLimit) ;
        }
      }
    }
    zoffset = G4AffineTransform_NetTranslation(&pTransform).z;
    zMin    = zoffset - This->fDz ;
    zMax    = zoffset + This->fDz ;

    if (G4VoxelLimits_IsZLimited(&pVoxelLimit))
    {
      if ( zMin > G4VoxelLimits_GetMaxZExtent(&pVoxelLimit)+kCarTolerance ||
           zMax < G4VoxelLimits_GetMinZExtent(&pVoxelLimit)-kCarTolerance   ) return false ;
      else
      {
        if (zMin < G4VoxelLimits_GetMinZExtent(&pVoxelLimit))
        {
          zMin = G4VoxelLimits_GetMinZExtent(&pVoxelLimit) ;
        }
        if (zMax > G4VoxelLimits_GetMaxZExtent(&pVoxelLimit))
        {
          zMax = G4VoxelLimits_GetMaxZExtent(&pVoxelLimit) ;
        }
      }
    }
    switch (pAxis)
    {
      case kXAxis:
        *pMin = xMin ;
        *pMax = xMax ;
        break ;
      case kYAxis:
        *pMin=yMin;
        *pMax=yMax;
        break;
      case kZAxis:
        *pMin=zMin;
        *pMax=zMax;
        break;
      default:
        break;
    }
    *pMin -= kCarTolerance ;
    *pMax += kCarTolerance ;

    return true;
  }
  else  // General rotated case - create and clip mesh to boundaries
  {
	const G4AffineTransform invt = G4AffineTransform_Inverse( &pTransform );
		
    G4bool existsAfterClip = false ;
    G4ThreeVectorList* vertices = G4Box_CreateRotatedVertices(This,pTransform);
    
    *pMin = +kInfinity ;
    *pMax = -kInfinity ;

    // Calculate rotated vertex coordinates

    G4VSolid_ClipCrossSection(vertices,0,&pVoxelLimit,pAxis,*pMin,*pMax) ;
    G4VSolid_ClipCrossSection(vertices,4,&pVoxelLimit,pAxis,*pMin,*pMax) ;
    G4VSolid_ClipBetweenSections(vertices,0,&pVoxelLimit,pAxis,*pMin,*pMax) ;

    if (G4VoxelLimits_IsLimited_axis(&pVoxelLimit,pAxis) == false) 
    {  
      if ( *pMin != kInfinity || *pMax != -kInfinity ) 
      {
        existsAfterClip = true ;

        // Add 2*tolerance to avoid precision troubles

        *pMin           -= kCarTolerance;
        *pMax           += kCarTolerance;
      }
    }      
    else
    {
      G4ThreeVector clipCentre = G4ThreeVector_create(
       ( G4VoxelLimits_GetMinXExtent(&pVoxelLimit)+G4VoxelLimits_GetMaxXExtent(&pVoxelLimit))*0.5,
       ( G4VoxelLimits_GetMinYExtent(&pVoxelLimit)+G4VoxelLimits_GetMaxYExtent(&pVoxelLimit))*0.5,
       ( G4VoxelLimits_GetMinZExtent(&pVoxelLimit)+G4VoxelLimits_GetMaxZExtent(&pVoxelLimit))*0.5);

      if ( *pMin != kInfinity || *pMax != -kInfinity )
      {
        existsAfterClip = true ;
  

        // Check to see if endpoints are in the solid

		G4ThreeVector_set_coord(&clipCentre, pAxis, G4VoxelLimits_GetMinExtent(&pVoxelLimit,pAxis));

		{
		const G4ThreeVector tp = G4AffineTransform_TransformPoint(&invt, clipCentre);

        if (G4Box_Inside(This, tp) != kOutside)
        {
          *pMin = G4VoxelLimits_GetMinExtent(&pVoxelLimit,pAxis);
        }
        else
        {
          *pMin -= kCarTolerance;
        }
		}
		
        G4ThreeVector_set_coord(&clipCentre, pAxis, G4VoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis));

		{
		const G4ThreeVector tp = G4AffineTransform_TransformPoint(&invt, clipCentre);

        if (G4Box_Inside(This, tp) != kOutside)
        {
          *pMax = G4VoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis);
        }
        else
        {
          *pMax += kCarTolerance;
        }
        
		}
      }

      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.
      
      else
      {
	      const G4ThreeVector tp = G4AffineTransform_TransformPoint(&invt, clipCentre);
	      
		  if (G4Box_Inside(This, tp) != kOutside)
		  {
			existsAfterClip = true ;
			*pMin            = G4VoxelLimits_GetMinExtent(&pVoxelLimit,pAxis);
			*pMax            = G4VoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis);
		  }
	  }
    }
    delete vertices;
    return existsAfterClip;
  } 
}

#endif // voxel nav
#endif // host code
