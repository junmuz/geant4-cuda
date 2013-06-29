
/** G4Cons implementation, based on G4Cons.cc of Geant 4.9.3 */

// --------------------------------------------------------------------

#include "G4Cons.h"

#define GET_hDPhi (0.5*This->fDPhi) // half delta phi
#define GET_cPhi (This->fSPhi + GET_hDPhi)
#define GET_ePhi (This->fSPhi + This->fDPhi)
  
#ifdef DISABLE_CACHE_MEMBERS
  #define GET_sinCPhi sin(GET_cPhi)
  #define GET_cosCPhi cos(GET_cPhi)
  #define GET_cosHDPhiIT cos(GET_hDPhi - 0.5*K_GEOMETRY_ANG_TOLERANCE)
  #define GET_cosHDPhiOT cos(GET_hDPhi + 0.5*K_GEOMETRY_ANG_TOLERANCE)
  #define GET_sinSPhi sin(This->fSPhi)
  #define GET_cosSPhi cos(This->fSPhi)
  #define GET_sinEPhi sin(GET_ePhi)
  #define GET_cosEPhi cos(GET_ePhi)
#else
  #define GET_sinCPhi (This->sinCPhi)  
  #define GET_cosCPhi (This->cosCPhi)
  #define GET_cosHDPhiIT (This->cosHDPhiIT)
  #define GET_cosHDPhiOT (This->cosHDPhiOT)
  #define GET_sinSPhi (This->sinSPhi)
  #define GET_cosSPhi (This->cosSPhi)
  #define GET_sinEPhi (This->sinEPhi)
  #define GET_cosEPhi (This->cosEPhi)
#endif

#ifdef HOST_CODE

INLINEFUNC
G4double G4Cons_GetInnerRadiusMinusZ( const G4Cons *This )
{
  return This->fRmin1 ;
}

INLINEFUNC
G4double G4Cons_GetOuterRadiusMinusZ( const G4Cons *This )
{
  return This->fRmax1 ;
}

INLINEFUNC
G4double G4Cons_GetInnerRadiusPlusZ( const G4Cons *This )
{
  return This->fRmin2 ;
}

INLINEFUNC
G4double G4Cons_GetOuterRadiusPlusZ( const G4Cons *This )
{
  return This->fRmax2 ;
}

INLINEFUNC
G4double G4Cons_GetZHalfLength( const G4Cons *This )
{
  return This->fDz ;
}

#ifdef ENABLE_SLICED_CONS
INLINEFUNC  
G4double G4Cons_GetStartPhiAngle( const G4Cons *This )
{
  return This->fSPhi ;
}

INLINEFUNC
G4double G4Cons_GetDeltaPhiAngle( const G4Cons *This )
{
  return This->fDPhi;
}

INLINEFUNC 
void G4Cons_InitializeTrigonometry( G4Cons *This )
{
#ifndef DISABLE_CACHE_MEMBERS
  This->sinCPhi    = std::sin(GET_cPhi);
  This->cosCPhi    = std::cos(GET_cPhi);
  This->cosHDPhiIT = std::cos(GET_hDPhi - 0.5*K_GEOMETRY_ANG_TOLERANCE); // inner/outer tol half dphi
  This->cosHDPhiOT = std::cos(GET_hDPhi + 0.5*K_GEOMETRY_ANG_TOLERANCE);
  This->sinSPhi = std::sin(This->fSPhi);
  This->cosSPhi = std::cos(This->fSPhi);
  This->sinEPhi = std::sin(GET_ePhi);
  This->cosEPhi = std::cos(GET_ePhi);
#else
  (void)This;
#endif
}

INLINEFUNC void G4Cons_CheckSPhiAngle( G4Cons *This, G4double sPhi )
{
  // Ensure fSphi in 0-2PI or -2PI-0 range if shape crosses 0

  if ( sPhi < 0 )
  {
    This->fSPhi = twopi - std::fmod(std::fabs(sPhi),twopi);
  }
  else
  {
    This->fSPhi = std::fmod(sPhi,twopi) ;
  }
  if ( This->fSPhi+This->fDPhi > twopi )
  {
    This->fSPhi -= twopi ;
  }
}

INLINEFUNC void G4Cons_CheckDPhiAngle( G4Cons *This, G4double dPhi )
{
  This->fPhiFullCone = true;
  if ( dPhi >= twopi-K_GEOMETRY_ANG_TOLERANCE*0.5 )
  {
    This->fDPhi=twopi;
    This->fSPhi=0;
  }
  else
  {
    This->fPhiFullCone = false;
    myAssert( dPhi > 0 );
    This->fDPhi = dPhi;
  }
}

INLINEFUNC void G4Cons_CheckPhiAngles( G4Cons *This, G4double sPhi, G4double dPhi )
{
  G4Cons_CheckDPhiAngle(This,dPhi);
  if ( (This->fDPhi<twopi) && (sPhi) ) { G4Cons_CheckSPhiAngle(This,sPhi); }
  G4Cons_InitializeTrigonometry(This);
}
#endif

// --------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////
//
// Private enum: Not for external use - used by distanceToOut

//enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};

//////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//               - note if pDPhi>2PI then reset to 2PI

SOLIDINLINE void G4Cons_ctor( G4Cons *This,
                      G4double  pRmin1, G4double pRmax1,
                      G4double  pRmin2, G4double pRmax2,
                      G4double pDz,
                      G4double pSPhi, G4double pDPhi )
{
	This->solid.type = kCons;
	
#ifdef ENABLE_SLICED_CONS
	This->fSPhi = 0;
	This->fDPhi = 0;
#else
	(void)pSPhi;
	(void)pDPhi;
#endif

	// Check z-len
	//
	myAssert( pDz > 0 );
	This->fDz = pDz;

	// Check radii
	//
	myAssert ( (pRmin1<=pRmax1) && (pRmin2<=pRmax2) && (pRmin1>=0) && (pRmin2>=0) );
  
	This->fRmin1 = pRmin1 ; 
	This->fRmax1 = pRmax1 ;
	This->fRmin2 = pRmin2 ; 
	This->fRmax2 = pRmax2 ;
	if( (pRmin1 == 0.0) && (pRmin2 > 0.0) ) { This->fRmin1 = 1e3*K_GEOMETRY_RAD_TOLERANCE ; }
	if( (pRmin2 == 0.0) && (pRmin1 > 0.0) ) { This->fRmin2 = 1e3*K_GEOMETRY_RAD_TOLERANCE ; }

#ifdef ENABLE_SLICED_CONS
	// Check angles
	//
	G4Cons_CheckPhiAngles(This, pSPhi, pDPhi);
#endif
}

#endif // host code

/////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

SOLIDINLINE EInside G4Cons_Inside(GEOMETRYLOC const G4Cons *This, G4ThreeVector p)
{
  G4double r2, rl, rh, tolRMin, tolRMax; // rh2, rl2 ;
  EInside in;
  const G4double halfCarTolerance=K_GEOMETRY_CAR_TOLERANCE*0.5;
  const G4double halfRadTolerance=K_GEOMETRY_RAD_TOLERANCE*0.5;
#ifdef ENABLE_SLICED_CONS
  const G4double halfAngTolerance=K_GEOMETRY_ANG_TOLERANCE*0.5;
  G4double pPhi;
#endif

  if (fabs(p.z) > This->fDz + halfCarTolerance )  { return in = kOutside; }
  else if(fabs(p.z) >= This->fDz - halfCarTolerance )    { in = kSurface; }
  else                                                    { in = kInside;  }

  r2 = p.x*p.x + p.y*p.y ;
  rl = 0.5*(This->fRmin2*(p.z + This->fDz) + This->fRmin1*(This->fDz - p.z))/This->fDz ;
  rh = 0.5*(This->fRmax2*(p.z+This->fDz)+This->fRmax1*(This->fDz-p.z))/This->fDz;

  // rh2 = rh*rh;

  tolRMin = rl - halfRadTolerance;
  if ( tolRMin < 0 )  { tolRMin = 0; }
  tolRMax = rh + halfRadTolerance;

  if ( (r2<tolRMin*tolRMin) || (r2>tolRMax*tolRMax) ) { return in = kOutside; }

  if (rl) { tolRMin = rl + halfRadTolerance; }
  else    { tolRMin = 0.0; }
  tolRMax = rh - halfRadTolerance;
      
  if (in == kInside) // else it's kSurface already
  {
     if ( (r2 < tolRMin*tolRMin) || (r2 >= tolRMax*tolRMax) ) { in = kSurface; }
  }
  
#ifdef ENABLE_SLICED_CONS
  if ( !This->fPhiFullCone && ((p.x != 0.0) || (p.y != 0.0)) )
  {
    pPhi = atan2(p.y,p.x) ;

    if ( pPhi < This->fSPhi - halfAngTolerance  )             { pPhi += twopi; }
    else if ( pPhi > This->fSPhi + This->fDPhi + halfAngTolerance ) { pPhi -= twopi; }
    
    if ( (pPhi < This->fSPhi - halfAngTolerance) ||          
         (pPhi > This->fSPhi + This->fDPhi + halfAngTolerance) )  { return in = kOutside; }
      
    else if (in == kInside)  // else it's kSurface anyway already
    {
       if ( (pPhi < This->fSPhi + halfAngTolerance) || 
            (pPhi > This->fSPhi + This->fDPhi - halfAngTolerance) )  { in = kSurface; }
    }
  }
  else if ( !This->fPhiFullCone ) { in = kSurface; }
#endif

  return in ;
}

#if defined(HOST_CODE) && defined(ENABLE_VOXEL_NAVIGATION)

#include "G4VSolid.c"

////////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -This->fDz cross section
//          [4-7] +This->fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

INLINEFUNC G4ThreeVectorList*
G4Cons_CreateRotatedVertices(const G4Cons *This, const G4AffineTransform* pTransform)
{
  G4ThreeVectorList* vertices ;
  G4ThreeVector vertex0, vertex1, vertex2, vertex3 ;
  G4double meshAngle, meshRMax1, meshRMax2, crossAngle;
  G4double cosCrossAngle, sinCrossAngle, sAngle ;
  G4double rMaxX1, rMaxX2, rMaxY1, rMaxY2, rMinX1, rMinX2, rMinY1, rMinY2 ;
  G4int crossSection, noCrossSections ;

  // Compute no of cross-sections necessary to mesh cone
    
#ifdef ENABLE_SLICED_CONS
  noCrossSections = G4int(This->fDPhi/K_GEOMETRY_MESH_ANGLE_DEFAULT) + 1 ;
#else
  noCrossSections = G4int(twopi/K_GEOMETRY_MESH_ANGLE_DEFAULT) + 1 ;
#endif

  if (noCrossSections < K_GEOMETRY_MIN_MESH_SECTIONS)
  {
    noCrossSections = K_GEOMETRY_MIN_MESH_SECTIONS ;
  }
  else if (noCrossSections > K_GEOMETRY_MAX_MESH_SECTIONS)
  {
    noCrossSections = K_GEOMETRY_MAX_MESH_SECTIONS ;
  }
#ifdef ENABLE_SLICED_CONS
  meshAngle = This->fDPhi/(noCrossSections - 1) ;
#else
  meshAngle = twopi/(noCrossSections - 1) ;
#endif

  meshRMax1 = This->fRmax1/std::cos(meshAngle*0.5) ;
  meshRMax2 = This->fRmax2/std::cos(meshAngle*0.5) ;

  // If complete in phi, set start angle such that mesh will be at RMax
  // on the x axis. Will give better extent calculations when not rotated.

#ifdef ENABLE_SLICED_CONS
  if ( This->fPhiFullCone && (This->fSPhi == 0.0) )
  {
    sAngle = -meshAngle*0.5 ;
  }
  else
  {
    sAngle = This->fSPhi ;
  } 
#else
  sAngle = -meshAngle*0.5 ;
#endif
  vertices = new G4ThreeVectorList();
  vertices->reserve(noCrossSections*4) ;

  if (vertices)
  {
    for (crossSection = 0 ; crossSection < noCrossSections ; crossSection++)
    {
      // Compute coordinates of cross section at section crossSection

      crossAngle    = sAngle + crossSection*meshAngle ;
      cosCrossAngle = std::cos(crossAngle) ;
      sinCrossAngle = std::sin(crossAngle) ;

      rMaxX1 = meshRMax1*cosCrossAngle ;
      rMaxY1 = meshRMax1*sinCrossAngle ;
      rMaxX2 = meshRMax2*cosCrossAngle ;
      rMaxY2 = meshRMax2*sinCrossAngle ;
        
      rMinX1 = This->fRmin1*cosCrossAngle ;
      rMinY1 = This->fRmin1*sinCrossAngle ;
      rMinX2 = This->fRmin2*cosCrossAngle ;
      rMinY2 = This->fRmin2*sinCrossAngle ;
        
      vertex0 = G4ThreeVector_create(rMinX1,rMinY1,-This->fDz) ;
      vertex1 = G4ThreeVector_create(rMaxX1,rMaxY1,-This->fDz) ;
      vertex2 = G4ThreeVector_create(rMaxX2,rMaxY2,+This->fDz) ;
      vertex3 = G4ThreeVector_create(rMinX2,rMinY2,+This->fDz) ;

      vertices->push_back(G4AffineTransform_TransformPoint(pTransform,vertex0)) ;
      vertices->push_back(G4AffineTransform_TransformPoint(pTransform,vertex1)) ;
      vertices->push_back(G4AffineTransform_TransformPoint(pTransform,vertex2)) ;
      vertices->push_back(G4AffineTransform_TransformPoint(pTransform,vertex3)) ;
    }
  }
  else
  {
	myAssert(false);
  }

  return vertices ;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

SOLIDINLINE
G4bool G4Cons_CalculateExtent(
			   const G4Cons *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax)
{
  if ( !G4AffineTransform_IsRotated(&pTransform)
#ifdef ENABLE_SLICED_CONS
	&& (This->fDPhi == twopi)
#endif
    && (This->fRmin1 == 0) && (This->fRmin2 == 0) )
  {
    // Special case handling for unrotated solid cones
    // Compute z/x/y mins and maxs for bounding box respecting limits,
    // with early returns if outside limits. Then switch() on pAxis,
    // and compute exact x and y limit for x/y case
      
    G4double xoffset, xMin, xMax ;
    G4double yoffset, yMin, yMax ;
    G4double zoffset, zMin, zMax ;

    G4double diff1, diff2, maxDiff, newMin, newMax, RMax ;
    G4double xoff1, xoff2, yoff1, yoff2 ;
      
    zoffset = G4AffineTransform_NetTranslation(&pTransform).z;
    zMin    = zoffset - This->fDz ;
    zMax    = zoffset + This->fDz ;

    if (G4VoxelLimits_IsZLimited(&pVoxelLimit))
    {
      if( (zMin > G4VoxelLimits_GetMaxZExtent(&pVoxelLimit) + K_GEOMETRY_CAR_TOLERANCE) || 
          (zMax < G4VoxelLimits_GetMinZExtent(&pVoxelLimit) - K_GEOMETRY_CAR_TOLERANCE)  )
      {
        return false ;
      }
      else
      {
        if ( zMin < G4VoxelLimits_GetMinZExtent(&pVoxelLimit) )
        {
          zMin = G4VoxelLimits_GetMinZExtent(&pVoxelLimit) ;
        }
        if ( zMax > G4VoxelLimits_GetMaxZExtent(&pVoxelLimit) )
        {
          zMax = G4VoxelLimits_GetMaxZExtent(&pVoxelLimit) ;
        }
      }
    }
    xoffset = G4AffineTransform_NetTranslation(&pTransform).x ;
    RMax    = (This->fRmax2 >= This->fRmax1) ?  zMax : zMin  ;                  
    xMax    = xoffset + (This->fRmax1 + This->fRmax2)*0.5 + 
              (RMax - zoffset)*(This->fRmax2 - This->fRmax1)/(2*This->fDz) ;
    xMin    = 2*xoffset-xMax ;

    if (G4VoxelLimits_IsXLimited(&pVoxelLimit))
    {
      if ( (xMin > G4VoxelLimits_GetMaxXExtent(&pVoxelLimit) + K_GEOMETRY_CAR_TOLERANCE) || 
           (xMax < G4VoxelLimits_GetMinXExtent(&pVoxelLimit) - K_GEOMETRY_CAR_TOLERANCE)  )
      {
        return false ;
      }
      else
      {
        if ( xMin < G4VoxelLimits_GetMinXExtent(&pVoxelLimit) )
        {
          xMin = G4VoxelLimits_GetMinXExtent(&pVoxelLimit) ;
        }
        if ( xMax > G4VoxelLimits_GetMaxXExtent(&pVoxelLimit) )
        {
          xMax=G4VoxelLimits_GetMaxXExtent(&pVoxelLimit) ;
        }
      }
    }
    yoffset = G4AffineTransform_NetTranslation(&pTransform).y ;
    yMax    = yoffset + (This->fRmax1 + This->fRmax2)*0.5 + 
              (RMax - zoffset)*(This->fRmax2 - This->fRmax1)/(2*This->fDz) ;
    yMin    = 2*yoffset-yMax ;
    RMax    = yMax - yoffset ;  // = max radius due to Zmax/Zmin cuttings

    if (G4VoxelLimits_IsYLimited(&pVoxelLimit))
    {
      if ( (yMin > G4VoxelLimits_GetMaxYExtent(&pVoxelLimit) + K_GEOMETRY_CAR_TOLERANCE) || 
           (yMax < G4VoxelLimits_GetMinYExtent(&pVoxelLimit) - K_GEOMETRY_CAR_TOLERANCE)  )
      {
        return false ;
      }
      else
      {
        if ( yMin < G4VoxelLimits_GetMinYExtent(&pVoxelLimit) )
        {
          yMin = G4VoxelLimits_GetMinYExtent(&pVoxelLimit) ;
        }
        if ( yMax > G4VoxelLimits_GetMaxYExtent(&pVoxelLimit) )
        {
          yMax = G4VoxelLimits_GetMaxYExtent(&pVoxelLimit) ;
        }
      }
    }    
    switch (pAxis) // Known to cut cones
    {
      case kXAxis:
        yoff1 = yoffset - yMin ;
        yoff2 = yMax - yoffset ;

        if ((yoff1 >= 0) && (yoff2 >= 0)) // Y limits cross max/min x
        {                                 // => no change
          *pMin = xMin ;
          *pMax = xMax ;
        }
        else
        {
          // Y limits don't cross max/min x => compute max delta x,
          // hence new mins/maxs
         
          diff1   = std::sqrt(RMax*RMax - yoff1*yoff1) ;
          diff2   = std::sqrt(RMax*RMax - yoff2*yoff2) ;
          maxDiff = (diff1>diff2) ? diff1:diff2 ;
          newMin  = xoffset - maxDiff ;
          newMax  = xoffset + maxDiff ;
          *pMin    = ( newMin < xMin ) ? xMin : newMin  ;
          *pMax    = ( newMax > xMax) ? xMax : newMax ;
        } 
      break ;

      case kYAxis:
        xoff1 = xoffset - xMin ;
        xoff2 = xMax - xoffset ;

        if ((xoff1 >= 0) && (xoff2 >= 0) ) // X limits cross max/min y
        {                                  // => no change
          *pMin = yMin ;
          *pMax = yMax ;
        }
        else
        {
          // X limits don't cross max/min y => compute max delta y,
          // hence new mins/maxs

          diff1   = std::sqrt(RMax*RMax - xoff1*xoff1) ;
          diff2   = std::sqrt(RMax*RMax-xoff2*xoff2) ;
          maxDiff = (diff1 > diff2) ? diff1:diff2 ;
          newMin  = yoffset - maxDiff ;
          newMax  = yoffset + maxDiff ;
          *pMin    = (newMin < yMin) ? yMin : newMin ;
          *pMax    = (newMax > yMax) ? yMax : newMax ;
        }
      break ;

      case kZAxis:
        *pMin = zMin ;
        *pMax = zMax ;
      break ;
      
      default:
      break ;
    }
    *pMin -= K_GEOMETRY_CAR_TOLERANCE ;
    *pMax += K_GEOMETRY_CAR_TOLERANCE ;

    return true ;
  }
  else   // Calculate rotated vertex coordinates
  {
    G4int i, noEntries, noBetweenSections4 ;
    G4bool existsAfterClip = false ;
    G4ThreeVectorList* vertices = G4Cons_CreateRotatedVertices(This, &pTransform);

    *pMin = +kInfinity ;
    *pMax = -kInfinity ;

    noEntries          = vertices->size() ;
    noBetweenSections4 = noEntries-4 ;
      
    for ( i = 0 ; i < noEntries ; i += 4 )
    {
      G4VSolid_ClipCrossSection(vertices, i, &pVoxelLimit, pAxis, *pMin, *pMax) ;
    }
    for ( i = 0 ; i < noBetweenSections4 ; i += 4 )
    {
      G4VSolid_ClipBetweenSections(vertices, i, &pVoxelLimit, pAxis, *pMin, *pMax) ;
    }    
    if ( (*pMin != kInfinity) || (*pMax != -kInfinity) )
    {
      existsAfterClip = true ;
        
      // Add 2*tolerance to avoid precision troubles

      *pMin -= K_GEOMETRY_CAR_TOLERANCE ;
      *pMax += K_GEOMETRY_CAR_TOLERANCE ;
    }
    else
    {
      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.
       
      G4ThreeVector clipCentre = G4ThreeVector_create(
      (G4VoxelLimits_GetMinXExtent(&pVoxelLimit) + G4VoxelLimits_GetMaxXExtent(&pVoxelLimit))*0.5,
      (G4VoxelLimits_GetMinYExtent(&pVoxelLimit) + G4VoxelLimits_GetMaxYExtent(&pVoxelLimit))*0.5,
      (G4VoxelLimits_GetMinZExtent(&pVoxelLimit) + G4VoxelLimits_GetMaxZExtent(&pVoxelLimit))*0.5  ) ;
        
      
      const G4AffineTransform invt = G4AffineTransform_Inverse(&pTransform);
      if (G4Cons_Inside(This,G4AffineTransform_TransformPoint(&invt,clipCentre)) != kOutside)
      {
        existsAfterClip = true ;
        *pMin            = G4VoxelLimits_GetMinExtent(&pVoxelLimit,pAxis) ;
        *pMax            = G4VoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis) ;
      }
    }
    delete vertices ;
    return existsAfterClip ;
  }
}

#endif // voxel nav & host code

////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

INLINEFUNC
G4ThreeVector G4Cons_ApproxSurfaceNormal( GEOMETRYLOC const G4Cons *This, G4ThreeVector p )
{
  enum {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ} side ;
  G4ThreeVector norm ;
  G4double rho;
  G4double distZ, distRMin, distRMax, distMin ;
  G4double tanRMin, secRMin, pRMin, widRMin ;
  G4double tanRMax, secRMax, pRMax, widRMax ;
#ifdef ENABLE_SLICED_CONS
  G4double distSPhi, distEPhi, phi;
#endif

  distZ = fabs(fabs(p.z) - This->fDz) ;
  rho   = sqrt(p.x*p.x + p.y*p.y) ;

  tanRMin  = (This->fRmin2 - This->fRmin1)*0.5/This->fDz ;
  secRMin  = sqrt(1 + tanRMin*tanRMin) ;
  pRMin    = rho - p.z*tanRMin ;
  widRMin  = This->fRmin2 - This->fDz*tanRMin ;
  distRMin = fabs(pRMin - widRMin)/secRMin ;

  tanRMax  = (This->fRmax2 - This->fRmax1)*0.5/This->fDz ;
  secRMax  = sqrt(1+tanRMax*tanRMax) ;
  pRMax    = rho - p.z*tanRMax ;
  widRMax  = This->fRmax2 - This->fDz*tanRMax ;
  distRMax = fabs(pRMax - widRMax)/secRMax ;
  
  if (distRMin < distRMax)  // First minimum
  {
    if (distZ < distRMin)
    {
      distMin = distZ ;
      side    = kNZ ;
    }
    else
    {
      distMin = distRMin ;
      side    = kNRMin ;
    }
  }
  else
  {
    if (distZ < distRMax)
    {
      distMin = distZ ;
      side    = kNZ ;
    }
    else
    {
      distMin = distRMax ;
      side    = kNRMax ;
    }
  }
#ifdef ENABLE_SLICED_CONS
  if ( !This->fPhiFullCone && rho )  // Protected against (0,0,z) 
  {
    phi = atan2(p.y,p.x) ;

    if (phi < 0)  { phi += twopi; }

    if (This->fSPhi < 0)  { distSPhi = fabs(phi - (This->fSPhi + twopi))*rho; }
    else            { distSPhi = fabs(phi - This->fSPhi)*rho; }

    distEPhi = fabs(phi - This->fSPhi - This->fDPhi)*rho ;

    // Find new minimum

    if (distSPhi < distEPhi)
    {
      if (distSPhi < distMin)  { side = kNSPhi; }
    }
    else 
    {
      if (distEPhi < distMin)  { side = kNEPhi; }
    }
  }
#endif 
  switch (side)
  {
    case kNRMin:      // Inner radius
      rho *= secRMin ;
      norm = G4ThreeVector_create(-p.x/rho, -p.y/rho, tanRMin/secRMin) ;
      break ;
    case kNRMax:      // Outer radius
      rho *= secRMax ;
      norm = G4ThreeVector_create(p.x/rho, p.y/rho, -tanRMax/secRMax) ;
      break ;
    case kNZ:      // +/- dz
      if (p.z > 0)  { norm = G4ThreeVector_create(0,0,1);  }
      else            { norm = G4ThreeVector_create(0,0,-1); }
      break ;
#ifdef ENABLE_SLICED_CONS
    case kNSPhi:
      norm = G4ThreeVector_create(sin(This->fSPhi), -cos(This->fSPhi), 0) ;
      break ;
    case kNEPhi:
      norm=G4ThreeVector_create(-sin(This->fSPhi+This->fDPhi), cos(This->fSPhi+This->fDPhi), 0) ;
      break ;
#endif
    default:
      myAssert(false);
      break ;    
  }
  return norm ;
}


////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

SOLIDINLINE
G4ThreeVector G4Cons_SurfaceNormal(GEOMETRYLOC const G4Cons *This, G4ThreeVector p)
{
  G4int noSurfaces = 0;
  G4double rho;
  G4double distZ, distRMin, distRMax;
  G4double tanRMin, secRMin, pRMin, widRMin;
  G4double tanRMax, secRMax, pRMax, widRMax;

  const G4double delta  = 0.5*K_GEOMETRY_CAR_TOLERANCE;
  
  G4ThreeVector norm, sumnorm = G4ThreeVector_create(0.,0.,0.), nZ = G4ThreeVector_create(0.,0.,1.);
  G4ThreeVector nR, nr = G4ThreeVector_create(0.,0.,0.);
#ifdef ENABLE_SLICED_CONS
  const G4double dAngle = 0.5*K_GEOMETRY_ANG_TOLERANCE;
  G4double distSPhi = kInfinity, distEPhi = kInfinity, pPhi;
  G4ThreeVector nPs, nPe;
#endif

  distZ = fabs(fabs(p.z) - This->fDz);
  rho   = sqrt(p.x*p.x + p.y*p.y);

  tanRMin  = (This->fRmin2 - This->fRmin1)*0.5/This->fDz;
  secRMin  = sqrt(1 + tanRMin*tanRMin);
  pRMin    = rho - p.z*tanRMin;
  widRMin  = This->fRmin2 - This->fDz*tanRMin;
  distRMin = fabs(pRMin - widRMin)/secRMin;

  tanRMax  = (This->fRmax2 - This->fRmax1)*0.5/This->fDz;
  secRMax  = sqrt(1+tanRMax*tanRMax);
  pRMax    = rho - p.z*tanRMax;
  widRMax  = This->fRmax2 - This->fDz*tanRMax;
  distRMax = fabs(pRMax - widRMax)/secRMax;

#ifdef ENABLE_SLICED_CONS
  if (!This->fPhiFullCone)   // Protected against (0,0,z) 
  {
    if ( rho )
    {
      pPhi = atan2(p.y,p.x);

      if (pPhi  < This->fSPhi-delta)            { pPhi += twopi; }
      else if (pPhi > This->fSPhi+This->fDPhi+delta)  { pPhi -= twopi; }

      distSPhi = fabs( pPhi - This->fSPhi ); 
      distEPhi = fabs( pPhi - This->fSPhi - This->fDPhi ); 
    }
    else if( !(This->fRmin1) || !(This->fRmin2) )
    {
      distSPhi = 0.; 
      distEPhi = 0.; 
    }
    nPs = G4ThreeVector_create(sin(This->fSPhi), -cos(This->fSPhi), 0);
    nPe = G4ThreeVector_create(-sin(This->fSPhi+This->fDPhi), cos(This->fSPhi+This->fDPhi), 0);
  }
#endif
  if ( rho > delta )   
  {
    nR = G4ThreeVector_create(p.x/rho/secRMax, p.y/rho/secRMax, -tanRMax/secRMax);
    if (This->fRmin1 || This->fRmin2)
    {
      nr = G4ThreeVector_create(-p.x/rho/secRMin,-p.y/rho/secRMin,tanRMin/secRMin);
    }
  }

  if( distRMax <= delta )
  {
    noSurfaces ++;
    G4ThreeVector_sum_assign(&sumnorm,nR); //sumnorm += nR;
  }
  if( (This->fRmin1 || This->fRmin2) && (distRMin <= delta) )
  {
    noSurfaces ++;
    G4ThreeVector_sum_assign(&sumnorm,nr); //sumnorm += nr;
  }
#ifdef ENABLE_SLICED_CONS
  if( !This->fPhiFullCone )   
  {
    if (distSPhi <= dAngle)
    {
      noSurfaces ++;
      G4ThreeVector_sum_assign(&sumnorm,nPs); //sumnorm += nPs;
    }
    if (distEPhi <= dAngle) 
    {
      noSurfaces ++;
      G4ThreeVector_sum_assign(&sumnorm,nPe); //sumnorm += nPe;
    }
  }
#endif
  if (distZ <= delta)  
  {
    noSurfaces ++;
    if ( p.z >= 0.)  { G4ThreeVector_sum_assign(&sumnorm,nZ); } //sumnorm += nZ;
    else               { G4ThreeVector_subtract_assign(&sumnorm,nZ); } //sumnorm -= nZ;
  }
  if ( noSurfaces == 0 )
  {
	 myAssert(false);
     norm = G4Cons_ApproxSurfaceNormal(This,p);
  }
  else if ( noSurfaces == 1 )  { norm = sumnorm; }
  else                         { norm = G4ThreeVector_unit(sumnorm); }

  return norm ;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outside cone, compute intersection with rmax1*0.5
//        - if at valid phi,z return
//        - if inside outer cone, handle case when on tolerant outer cone
//          boundary and heading inwards(->0 to in)
//
// -> Compute intersection with inner cone, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - `if valid' implies tolerant checking of intersection points
// - z, phi intersection from Tubs

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

SOLIDINLINE G4double G4Cons_DistanceToIn_full(
				GEOMETRYLOC const G4Cons *This,
				G4ThreeVector p,
				G4ThreeVector v)
{
  G4double snxt = kInfinity ;      // snxt = default return value
  //const G4double dRmax = 100*MIN(This->fRmax1,This->fRmax2);
  const G4double halfCarTolerance=K_GEOMETRY_CAR_TOLERANCE*0.5;
  const G4double halfRadTolerance=K_GEOMETRY_RAD_TOLERANCE*0.5;

  G4double tanRMax,secRMax,rMaxAv,rMaxOAv ;  // Data for cones
  G4double tanRMin,secRMin,rMinAv,rMinOAv ;
  G4double rout,rin ;

  G4double tolORMin,tolORMin2,tolIRMin,tolIRMin2 ; // `generous' radii squared
  G4double tolORMax2,tolIRMax,tolIRMax2 ;
  G4double tolODz,tolIDz ;

  G4double s,xi,yi,zi,ri=0.,risec,rhoi2 ; // Intersection point vars

  G4double t1,t2,t3,b,c,d ;    // Quadratic solver variables 
  G4double nt1,nt2,nt3 ;
#ifdef ENABLE_SLICED_CONS
  G4double Dist,Comp,cosPsi ;
#endif

  G4ThreeVector Normal;

  // Cone Precalcs

  tanRMin = (This->fRmin2 - This->fRmin1)*0.5/This->fDz ;
  secRMin = sqrt(1.0 + tanRMin*tanRMin) ;
  rMinAv  = (This->fRmin1 + This->fRmin2)*0.5 ;

  if (rMinAv > halfRadTolerance)
  {
    rMinOAv = rMinAv - halfRadTolerance ;
  }
  else
  {
    rMinOAv = 0.0 ;
  }  
  tanRMax = (This->fRmax2 - This->fRmax1)*0.5/This->fDz ;
  secRMax = sqrt(1.0 + tanRMax*tanRMax) ;
  rMaxAv  = (This->fRmax1 + This->fRmax2)*0.5 ;
  rMaxOAv = rMaxAv + halfRadTolerance ;
   
  // Intersection with z-surfaces

  tolIDz = This->fDz - halfCarTolerance ;
  tolODz = This->fDz + halfCarTolerance ;

  if (fabs(p.z) >= tolIDz)
  {
    if ( p.z*v.z < 0 )    // at +Z going in -Z or visa versa
    {
      s = (fabs(p.z) - This->fDz)/fabs(v.z) ;  // Z intersect distance

      if( s < 0.0 )  { s = 0.0; }                      // negative dist -> zero

      xi   = p.x + s*v.x ;  // Intersection coords
      yi   = p.y + s*v.y ;
      rhoi2 = xi*xi + yi*yi  ;

      // Check validity of intersection
      // Calculate (outer) tolerant radi^2 at intersecion

      if (v.z > 0)
      {
        tolORMin  = This->fRmin1 - halfRadTolerance*secRMin ;
        tolIRMin  = This->fRmin1 + halfRadTolerance*secRMin ;
        tolIRMax  = This->fRmax1 - halfRadTolerance*secRMin ;
        tolORMax2 = (This->fRmax1 + halfRadTolerance*secRMax)*
                    (This->fRmax1 + halfRadTolerance*secRMax) ;
      }
      else
      {
        tolORMin  = This->fRmin2 - halfRadTolerance*secRMin ;
        tolIRMin  = This->fRmin2 + halfRadTolerance*secRMin ;
        tolIRMax  = This->fRmax2 - halfRadTolerance*secRMin ;
        tolORMax2 = (This->fRmax2 + halfRadTolerance*secRMax)*
                    (This->fRmax2 + halfRadTolerance*secRMax) ;
      }
      if ( tolORMin > 0 ) 
      {
        tolORMin2 = tolORMin*tolORMin ;
        tolIRMin2 = tolIRMin*tolIRMin ;
      }
      else                
      {
        tolORMin2 = 0.0 ;
        tolIRMin2 = 0.0 ;
      }
      if ( tolIRMax > 0 )  { tolIRMax2 = tolIRMax*tolIRMax; }     
      else                 { tolIRMax2 = 0.0; }
      
      if ( (tolIRMin2 <= rhoi2) && (rhoi2 <= tolIRMax2) )
      {
		  
#ifdef ENABLE_SLICED_CONS
        if ( !This->fPhiFullCone && rhoi2 )
        {
          // Psi = angle made with central (average) phi of shape

          cosPsi = (xi*GET_cosCPhi + yi*GET_sinCPhi)/sqrt(rhoi2) ;

          if (cosPsi >= GET_cosHDPhiIT)  { return s; }
        }
        else
        {
          return s;
        }
#else
		return s;
#endif
      }
    }
    else  // On/outside extent, and heading away  -> cannot intersect
    {
      return snxt ;  
    }
  }
    
// ----> Can not intersect z surfaces


// Intersection with outer cone (possible return) and
//                   inner cone (must also check phi)
//
// Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
//
// Intersects with x^2+y^2=(a*z+b)^2
//
// where a=tanRMax or tanRMin
//       b=rMaxAv  or rMinAv
//
// (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0 ;
//     t1                        t2                      t3  
//
//  \--------u-------/       \-----------v----------/ \---------w--------/
//

  t1   = 1.0 - v.z*v.z ;
  t2   = p.x*v.x + p.y*v.y ;
  t3   = p.x*p.x + p.y*p.y ;
  rin  = tanRMin*p.z + rMinAv ;
  rout = tanRMax*p.z + rMaxAv ;

  // Outer Cone Intersection
  // Must be outside/on outer cone for valid intersection

  nt1 = t1 - (tanRMax*v.z)*(tanRMax*v.z) ;
  nt2 = t2 - tanRMax*v.z*rout ;
  nt3 = t3 - rout*rout ;

  if (fabs(nt1) > K_GEOMETRY_RAD_TOLERANCE)  // Equation quadratic => 2 roots
  {
    b = nt2/nt1;
    c = nt3/nt1;
    d = b*b-c  ;
    if ( (nt3 > rout*K_GEOMETRY_RAD_TOLERANCE*secRMax) || (rout < 0) )
    {
      // If outside real cone (should be rho-rout>K_GEOMETRY_RAD_TOLERANCE*0.5
      // NOT rho^2 etc) saves a sqrt() at expense of accuracy

      if (d >= 0)
      {
          
        if ((rout < 0) && (nt3 <= 0))
        {
          // Inside `shadow cone' with -ve radius
          // -> 2nd root could be on real cone

          s = -b + sqrt(d) ;
        }
        else
        {
          if ((b <= 0) && (c >= 0)) // both >=0, try smaller root
          {
            s = -b - sqrt(d) ;
          }
          else
          {
            if ( c <= 0 ) // second >=0
            {
              s = -b + sqrt(d) ;
            }
            else  // both negative, travel away
            {
              return kInfinity ;
            }
          }
        }
        if ( s > 0 )  // If 'forwards'. Check z intersection
        {
		  // myAssert( s <= fRmax ); // TODO?
          /*if ( s>dRmax ) // Avoid rounding errors due to precision issues on
          {              // 64 bits systems. Split long distances and recompute
            G4double fTerm = s-fmod(s,dRmax);
            s = fTerm + DistanceToIn(p+fTerm*v,v);
          } */
          
          zi = p.z + s*v.z ;

          if (fabs(zi) <= tolODz)
          {
            // Z ok. Check phi intersection if reqd

#ifdef ENABLE_SLICED_CONS
            if ( This->fPhiFullCone )  { return s; }
            else
            {
              xi     = p.x + s*v.x ;
              yi     = p.y + s*v.y ;
              ri     = rMaxAv + zi*tanRMax ;
              cosPsi = (xi*GET_cosCPhi + yi*GET_sinCPhi)/ri ;

              if ( cosPsi >= GET_cosHDPhiIT )  { return s; }
            }
#else
			return s;
#endif
          }
        }                // end if (s>0)
      }
    }
    else
    {
      // Inside outer cone
      // check not inside, and heading through G4Cons (-> 0 to in)

      if ( ( t3  > (rin + halfRadTolerance*secRMin)*
                   (rin + halfRadTolerance*secRMin) )
        && (nt2 < 0) && (d >= 0) && (fabs(p.z) <= tolIDz) )
      {
        // Inside cones, delta r -ve, inside z extent
        // Point is on the Surface => check Direction using  Normal.dot(v)

        xi     = p.x ;
        yi     = p.y  ;
        risec  = sqrt(xi*xi + yi*yi)*secRMax ;
        Normal = G4ThreeVector_create(xi/risec,yi/risec,-tanRMax/secRMax) ;
#ifdef ENABLE_SLICED_CONS
        if ( !This->fPhiFullCone )
        {
          cosPsi = (p.x*GET_cosCPhi + p.y*GET_sinCPhi)/sqrt(t3) ;
          if ( cosPsi >= GET_cosHDPhiIT )
          {
            if ( G4ThreeVector_dot(Normal,v) <= 0 )  { return 0.0; }
          }
        }
        else
        {             
          if ( G4ThreeVector_dot(Normal,v) <= 0 )  { return 0.0; }
        }
#else
          if ( G4ThreeVector_dot(Normal,v) <= 0 )  { return 0.0; }
#endif
      }
    }
  }
  else  //  Single root case 
  {
    if ( fabs(nt2) > K_GEOMETRY_RAD_TOLERANCE )
    {
      s = -0.5*nt3/nt2 ;

      if ( s < 0 )  { return kInfinity; }   // travel away
      else  // s >= 0,  If 'forwards'. Check z intersection
      {
        zi = p.z + s*v.z ;

        if ((fabs(zi) <= tolODz) && (nt2 < 0))
        {
          // Z ok. Check phi intersection if reqd
#ifdef ENABLE_SLICED_CONS
          if ( This->fPhiFullCone )  { return s; }
          else
          {
            xi     = p.x + s*v.x ;
            yi     = p.y + s*v.y ;
            ri     = rMaxAv + zi*tanRMax ;
            cosPsi = (xi*GET_cosCPhi + yi*GET_sinCPhi)/ri ;

            if (cosPsi >= GET_cosHDPhiIT)  { return s; }
          }
#else
			return s;
#endif
        }
      }
    }
    else  //    travel || cone surface from its origin
    {
      s = kInfinity ;
    }
  }

  // Inner Cone Intersection
  // o Space is divided into 3 areas:
  //   1) Radius greater than real inner cone & imaginary cone & outside
  //      tolerance
  //   2) Radius less than inner or imaginary cone & outside tolarance
  //   3) Within tolerance of real or imaginary cones
  //      - Extra checks needed for 3's intersections
  //        => lots of duplicated code

  if (rMinAv)
  { 
    nt1 = t1 - (tanRMin*v.z)*(tanRMin*v.z) ;
    nt2 = t2 - tanRMin*v.z*rin ;
    nt3 = t3 - rin*rin ;
 
    if ( nt1 )
    {
      if ( nt3 > rin*K_GEOMETRY_RAD_TOLERANCE*secRMin )
      {
        // At radius greater than real & imaginary cones
        // -> 2nd root, with zi check

        b = nt2/nt1 ;
        c = nt3/nt1 ;
        d = b*b-c ;
        if (d >= 0)   // > 0
        {
          s = -b + sqrt(d) ;

          if ( s >= 0 )   // > 0
          {
			//myAssert( s <= dRmax ); //FIXME
            /*if ( s>dRmax ) // Avoid rounding errors due to precision issues on
            {              // 64 bits systems. Split long distance and recompute
              G4double fTerm = s-fmod(s,dRmax);
              s = fTerm + DistanceToIn(p+fTerm*v,v);
            } */
            zi = p.z + s*v.z ;

            if ( fabs(zi) <= tolODz )
            {
#ifdef ENABLE_SLICED_CONS
              if ( !This->fPhiFullCone )
              {
                xi     = p.x + s*v.x ;
                yi     = p.y + s*v.y ;
                ri     = rMinAv + zi*tanRMin ;
                cosPsi = (xi*GET_cosCPhi + yi*GET_sinCPhi)/ri ;

                if (cosPsi >= GET_cosHDPhiIT)
                { 
                  if ( s > halfRadTolerance )  { snxt=s; }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec  = sqrt(xi*xi + yi*yi)*secRMin ;
                    Normal = G4ThreeVector_create(-xi/risec,-yi/risec,tanRMin/secRMin);
                    if ( G4ThreeVector_dot(Normal,v) <= 0 )  { snxt = s; }
                  } 
                }
              }
              else
              {
#endif
                if ( s > halfRadTolerance )  { return s; }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x + s*v.x ;
                  yi     = p.y + s*v.y ;
                  risec  = sqrt(xi*xi + yi*yi)*secRMin ;
                  Normal = G4ThreeVector_create(-xi/risec,-yi/risec,tanRMin/secRMin) ;
                  if ( G4ThreeVector_dot(Normal,v) <= 0 )  { return s; }
                }
#ifdef ENABLE_SLICED_CONS
              }
#endif
            }
          }
        }
      }
      else  if ( nt3 < -rin*K_GEOMETRY_RAD_TOLERANCE*secRMin )
      {
        // Within radius of inner cone (real or imaginary)
        // -> Try 2nd root, with checking intersection is with real cone
        // -> If check fails, try 1st root, also checking intersection is
        //    on real cone

        b = nt2/nt1 ;
        c = nt3/nt1 ;
        d = b*b - c ;

        if ( d >= 0 )  // > 0
        {
          s  = -b + sqrt(d) ;
          zi = p.z + s*v.z ;
          ri = rMinAv + zi*tanRMin ;

          if ( ri > 0 )
          {
            if ( (s >= 0) && (fabs(zi) <= tolODz) )  // s > 0
            {
			  // myAssert( s <= dRmax ); // FIXME
              /*if ( s>dRmax ) // Avoid rounding errors due to precision issues
              {              // seen on 64 bits systems. Split and recompute
                G4double fTerm = s-fmod(s,dRmax);
                s = fTerm + DistanceToIn(p+fTerm*v,v);
              } */
#ifdef ENABLE_SLICED_CONS
              if ( !This->fPhiFullCone )
              {
                xi     = p.x + s*v.x ;
                yi     = p.y + s*v.y ;
                cosPsi = (xi*GET_cosCPhi + yi*GET_sinCPhi)/ri ;

                if (cosPsi >= GET_cosHDPhiOT)
                {
                  if ( s > halfRadTolerance )  { snxt=s; }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec  = sqrt(xi*xi + yi*yi)*secRMin ;
                    Normal = G4ThreeVector_create(-xi/risec,-yi/risec,tanRMin/secRMin);
                    if ( G4ThreeVector_dot(Normal,v) <= 0 )  { snxt = s; } 
                  }
                }
              }
              else
              {
#endif
                if( s > halfRadTolerance )  { return s; }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x + s*v.x ;
                  yi     = p.y + s*v.y ;
                  risec  = sqrt(xi*xi + yi*yi)*secRMin ;
                  Normal = G4ThreeVector_create(-xi/risec,-yi/risec,tanRMin/secRMin) ;
                  if ( G4ThreeVector_dot(Normal,v) <= 0 )  { return s; }
                } 
#ifdef ENABLE_SLICED_CONS
              }
#endif
            }
          }
          else
          {
            s  = -b - sqrt(d) ;
            zi = p.z + s*v.z ;
            ri = rMinAv + zi*tanRMin ;

            if ( (s >= 0) && (ri > 0) && (fabs(zi) <= tolODz) ) // s>0
            {
			  // myAssert( s <= dRmax ); // FIXME
              /*if ( s>dRmax ) // Avoid rounding errors due to precision issues
              {              // seen on 64 bits systems. Split and recompute
                G4double fTerm = s-fmod(s,dRmax);
                s = fTerm + DistanceToIn(p+fTerm*v,v);
              } */
#ifdef ENABLE_SLICED_CONS
              if ( !This->fPhiFullCone )
              {
                xi     = p.x + s*v.x ;
                yi     = p.y + s*v.y ;
                cosPsi = (xi*GET_cosCPhi + yi*GET_sinCPhi)/ri ;

                if (cosPsi >= GET_cosHDPhiIT)
                {
                  if ( s > halfRadTolerance )  { snxt=s; }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec  = sqrt(xi*xi + yi*yi)*secRMin ;
                    Normal = G4ThreeVector_create(-xi/risec,-yi/risec,tanRMin/secRMin);
                    if ( G4ThreeVector_dot(Normal,v) <= 0 )  { snxt = s; } 
                  }
                }
              }
              else
              {
#endif
                if ( s > halfRadTolerance )  { return s; }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x + s*v.x ;
                  yi     = p.y + s*v.y ;
                  risec  = sqrt(xi*xi + yi*yi)*secRMin ;
                  Normal = G4ThreeVector_create(-xi/risec,-yi/risec,tanRMin/secRMin) ;
                  if ( G4ThreeVector_dot(Normal,v) <= 0 )  { return s; }
                } 
#ifdef ENABLE_SLICED_CONS
              }
#endif
            }
          }
        }
      }
      else
      {
        // Within kRadTol*0.5 of inner cone (real OR imaginary)
        // ----> Check not travelling through (=>0 to in)
        // ----> if not:
        //    -2nd root with validity check

        if ( fabs(p.z) <= tolODz )
        {
          if ( nt2 > 0 )
          {
            // Inside inner real cone, heading outwards, inside z range

#ifdef ENABLE_SLICED_CONS
            if ( !This->fPhiFullCone )
            {
              cosPsi = (p.x*GET_cosCPhi + p.y*GET_sinCPhi)/sqrt(t3) ;

              if (cosPsi >= GET_cosHDPhiIT)  { return 0.0; }
            }
            else  { return 0.0; }
#else
			return 0.0;
#endif
          }
          else
          {
            // Within z extent, but not travelling through
            // -> 2nd root or kInfinity if 1st root on imaginary cone

            b = nt2/nt1 ;
            c = nt3/nt1 ;
            d = b*b - c ;

            if ( d >= 0 )   // > 0
            {
              s  = -b - sqrt(d) ;
              zi = p.z + s*v.z ;
              ri = rMinAv + zi*tanRMin ;
              
              if ( ri > 0 )   // 2nd root
              {
                s  = -b + sqrt(d) ;
                zi = p.z + s*v.z ;

                if ( (s >= 0) && (fabs(zi) <= tolODz) )  // s>0
                {
				  // myAssert( s <= dRmax ); // FIXME
                  /*if ( s>dRmax ) // Avoid rounding errors due to precision issue
                  {              // seen on 64 bits systems. Split and recompute
                    G4double fTerm = s-fmod(s,dRmax);
                    s = fTerm + DistanceToIn(p+fTerm*v,v);
                  } */
#ifdef ENABLE_SLICED_CONS
                  if ( !This->fPhiFullCone )
                  {
                    xi     = p.x + s*v.x ;
                    yi     = p.y + s*v.y ;
                    ri     = rMinAv + zi*tanRMin ;
                    cosPsi = (xi*GET_cosCPhi + yi*GET_sinCPhi)/ri ;

                    if ( cosPsi >= GET_cosHDPhiIT )  { snxt = s; }
                  }
                  else  { return s; }
#else
					return s;
#endif
                }
              }
              else  { return kInfinity; }
            }
          }
        }
        else   // 2nd root
        {
          b = nt2/nt1 ;
          c = nt3/nt1 ;
          d = b*b - c ;

          if ( d > 0 )
          {  
            s  = -b + sqrt(d) ;
            zi = p.z + s*v.z ;

            if ( (s >= 0) && (fabs(zi) <= tolODz) )  // s>0
            {
			  // myAssert( s <= dRmax ); // FIXME
              /*if ( s>dRmax ) // Avoid rounding errors due to precision issues
              {              // seen on 64 bits systems. Split and recompute
                G4double fTerm = s-fmod(s,dRmax);
                s = fTerm + DistanceToIn(p+fTerm*v,v);
              } */
#ifdef ENABLE_SLICED_CONS
              if ( !This->fPhiFullCone )
              {
                xi     = p.x + s*v.x;
                yi     = p.y + s*v.y;
                ri     = rMinAv + zi*tanRMin ;
                cosPsi = (xi*GET_cosCPhi + yi*GET_sinCPhi)/ri;

                if (cosPsi >= GET_cosHDPhiIT)  { snxt = s; }
              }
              else  { return s; }
#else
				return s;
#endif
            }
          }
        }
      }
    }
  }

  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to K_GEOMETRY_CAR_TOLERANCE*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> Should use some form of loop Construct

#ifdef ENABLE_SLICED_CONS
  if ( !This->fPhiFullCone )
  {
    // First phi surface (starting phi)

    Comp    = v.x*GET_sinSPhi - v.y*GET_cosSPhi ;
                    
    if ( Comp < 0 )    // Component in outwards normal dirn
    {
      Dist = (p.y*GET_cosSPhi - p.x*GET_sinSPhi) ;

      if (Dist < halfCarTolerance)
      {
        s = Dist/Comp ;

        if ( s < snxt )
        {
          if ( s < 0 )  { s = 0.0; }

          zi = p.z + s*v.z ;

          if ( fabs(zi) <= tolODz )
          {
            xi        = p.x + s*v.x ;
            yi        = p.y + s*v.y ;
            rhoi2     = xi*xi + yi*yi ;
            tolORMin2 = (rMinOAv + zi*tanRMin)*(rMinOAv + zi*tanRMin) ;
            tolORMax2 = (rMaxOAv + zi*tanRMax)*(rMaxOAv + zi*tanRMax) ;

            if ( (rhoi2 >= tolORMin2) && (rhoi2 <= tolORMax2) )
            {
              // z and r intersections good - check intersecting with
              // correct half-plane

              if ((yi*GET_cosCPhi - xi*GET_sinCPhi) <= 0 )  { snxt = s; }
            }
          }
        }
      }
    }

    // Second phi surface (Ending phi)

    Comp    = -(v.x*GET_sinEPhi - v.y*GET_cosEPhi) ;
        
    if ( Comp < 0 )   // Component in outwards normal dirn
    {
      Dist = -(p.y*GET_cosEPhi - p.x*GET_sinEPhi) ;
      if (Dist < halfCarTolerance)
      {
        s = Dist/Comp ;

        if ( s < snxt )
        {
          if ( s < 0 )  { s = 0.0; }

          zi = p.z + s*v.z ;

          if (fabs(zi) <= tolODz)
          {
            xi        = p.x + s*v.x ;
            yi        = p.y + s*v.y ;
            rhoi2     = xi*xi + yi*yi ;
            tolORMin2 = (rMinOAv + zi*tanRMin)*(rMinOAv + zi*tanRMin) ;
            tolORMax2 = (rMaxOAv + zi*tanRMax)*(rMaxOAv + zi*tanRMax) ;

            if ( (rhoi2 >= tolORMin2) && (rhoi2 <= tolORMax2) )
            {
              // z and r intersections good - check intersecting with
              // correct half-plane

              if ( (yi*GET_cosCPhi - xi*GET_sinCPhi) >= 0.0 )  { snxt = s; }
            }
          }
        }
      }
    }
  }
#endif
  if (snxt < halfCarTolerance)  { snxt = 0.; }

  return snxt ;
}

//////////////////////////////////////////////////////////////////////////////
// 
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

SOLIDINLINE
G4double G4Cons_DistanceToIn(GEOMETRYLOC const G4Cons *This, G4ThreeVector p)
{
  G4double safe=0.0, rho, safeR1, safeR2, safeZ;
  G4double tanRMin, secRMin, pRMin ;
  G4double tanRMax, secRMax, pRMax ;
#ifdef ENABLE_SLICED_CONS
  G4double safePhi, cosPsi;
#endif

  rho   = sqrt(p.x*p.x + p.y*p.y) ;
  safeZ = fabs(p.z) - This->fDz ;

  if ( This->fRmin1 || This->fRmin2 )
  {
    tanRMin = (This->fRmin2 - This->fRmin1)*0.5/This->fDz ;
    secRMin = sqrt(1.0 + tanRMin*tanRMin) ;
    pRMin   = tanRMin*p.z + (This->fRmin1 + This->fRmin2)*0.5 ;
    safeR1  = (pRMin - rho)/secRMin ;

    tanRMax = (This->fRmax2 - This->fRmax1)*0.5/This->fDz ;
    secRMax = sqrt(1.0 + tanRMax*tanRMax) ;
    pRMax   = tanRMax*p.z + (This->fRmax1 + This->fRmax2)*0.5 ;
    safeR2  = (rho - pRMax)/secRMax ;

    if ( safeR1 > safeR2) { safe = safeR1; }
    else                  { safe = safeR2; }
  }
  else
  {
    tanRMax = (This->fRmax2 - This->fRmax1)*0.5/This->fDz ;
    secRMax = sqrt(1.0 + tanRMax*tanRMax) ;
    pRMax   = tanRMax*p.z + (This->fRmax1 + This->fRmax2)*0.5 ;
    safe    = (rho - pRMax)/secRMax ;
  }
  if ( safeZ > safe )  { safe = safeZ; }

#ifdef ENABLE_SLICED_CONS
  if ( !This->fPhiFullCone && rho )
  {
    // Psi=angle from central phi to point

    cosPsi = (p.x*GET_cosCPhi + p.y*GET_sinCPhi)/rho ;

    if ( cosPsi < cos(This->fDPhi*0.5) ) // Point lies outside phi range
    {
      if ( (p.y*GET_cosCPhi - p.x*GET_sinCPhi) <= 0.0 )
      {
        safePhi = fabs(p.x*sin(This->fSPhi)-p.y*cos(This->fSPhi));
      }
      else
      {
        safePhi = fabs(p.x*GET_sinEPhi-p.y*GET_cosEPhi);
      }
      if ( safePhi > safe )  { safe = safePhi; }
    }
  }
#endif
  if ( safe < 0.0 )  { safe = 0.0; }

  return safe ;
}

///////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from 'inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

SOLIDINLINE G4double G4Cons_DistanceToOut_full(
			   GEOMETRYLOC const G4Cons *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n)
{
  typedef enum {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ} ESide;
  ESide side = kNull, sider = kNull, sidetol = kNull ;

  const G4double halfCarTolerance=K_GEOMETRY_CAR_TOLERANCE*0.5;
  const G4double halfRadTolerance=K_GEOMETRY_RAD_TOLERANCE*0.5;
  
#ifdef ENABLE_SLICED_CONS
  const G4double halfAngTolerance=K_GEOMETRY_ANG_TOLERANCE*0.5;
  const G4double pi = M_PI;
  
  G4double sphi, pDistS, compS, pDistE, compE, sphi2, vphi;
  ESide sidephi = kNull;
#endif

  G4double snxt,sr,pdist ;

  G4double tanRMax, secRMax, rMaxAv ;  // Data for outer cone
  G4double tanRMin, secRMin, rMinAv ;  // Data for inner cone

  G4double t1, t2, t3, rout, rin, nt1, nt2, nt3 ;
  G4double b, c, d, sr2, sr3 ;

  // Vars for intersection within tolerance

  G4double slentol = kInfinity ;

  // Vars for phi intersection:

  G4double xi, yi, risec ;
  G4double zi, ri, deltaRoi2 ;

  // Z plane intersection

  if ( v.z > 0.0 )
  {
    pdist = This->fDz - p.z ;

    if (pdist > halfCarTolerance)
    {
      snxt = pdist/v.z ;
      side = kPZ ;
    }
    else
    {
      if (calcNorm)
      {
        *n         = G4ThreeVector_create(0,0,1) ;
        *validNorm = true ;
      }
      return  snxt = 0.0;
    }
  }
  else if ( v.z < 0.0 )
  {
    pdist = This->fDz + p.z ;

    if ( pdist > halfCarTolerance)
    {
      snxt = -pdist/v.z ;
      side = kMZ ;
    }
    else
    {
      if ( calcNorm )
      {
        *n         = G4ThreeVector_create(0,0,-1) ;
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
  }
  else     // Travel perpendicular to z axis
  {
    snxt = kInfinity ;    
    side = kNull ;
  }

  // Radial Intersections
  //
  // Intersection with outer cone (possible return) and
  //                   inner cone (must also check phi)
  //
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=(a*z+b)^2
  //
  // where a=tanRMax or tanRMin
  //       b=rMaxAv  or rMinAv
  //
  // (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0 ;
  //     t1                        t2                      t3  
  //
  //  \--------u-------/       \-----------v----------/ \---------w--------/

  tanRMax = (This->fRmax2 - This->fRmax1)*0.5/This->fDz ;
  secRMax = sqrt(1.0 + tanRMax*tanRMax) ;
  rMaxAv  = (This->fRmax1 + This->fRmax2)*0.5 ;


  t1   = 1.0 - v.z*v.z ;      // since v normalised
  t2   = p.x*v.x + p.y*v.y ;
  t3   = p.x*p.x + p.y*p.y ;
  rout = tanRMax*p.z + rMaxAv ;

  nt1 = t1 - (tanRMax*v.z)*(tanRMax*v.z) ;
  nt2 = t2 - tanRMax*v.z*rout ;
  nt3 = t3 - rout*rout ;

  if (v.z > 0.0)
  {
    deltaRoi2 = snxt*snxt*t1 + 2*snxt*t2 + t3
                - This->fRmax2*(This->fRmax2 + K_GEOMETRY_RAD_TOLERANCE*secRMax);
  }
  else if ( v.z < 0.0 )
  {
    deltaRoi2 = snxt*snxt*t1 + 2*snxt*t2 + t3
                - This->fRmax1*(This->fRmax1 + K_GEOMETRY_RAD_TOLERANCE*secRMax);
  }
  else
  {
    deltaRoi2 = 1.0;
  }

  if ( nt1 && (deltaRoi2 > 0.0) )  
  {
    // Equation quadratic => 2 roots : second root must be leaving

    b = nt2/nt1 ;
    c = nt3/nt1 ;
    d = b*b - c ;

    if ( d >= 0 )
    {
      // Check if on outer cone & heading outwards
      // NOTE: Should use rho-rout>-K_GEOMETRY_RAD_TOLERANCE*0.5
        
      if (nt3 > -halfRadTolerance && nt2 >= 0 )
      {
        if (calcNorm)
        {
          risec      = sqrt(t3)*secRMax ;
          *validNorm = true ;
          *n         = G4ThreeVector_create(p.x/risec,p.y/risec,-tanRMax/secRMax);
        }
        return snxt=0 ;
      }
      else
      {
        sider = kRMax  ;
        sr    = -b - sqrt(d) ; // was +srqrt(d), vmg 28.04.99
        zi    = p.z + sr*v.z ;
        ri    = tanRMax*zi + rMaxAv ;
          
        if ((ri >= 0) && (-halfRadTolerance <= sr) && (sr <= halfRadTolerance))
        {
          // An intersection within the tolerance
          //   we will Store it in case it is good -
          // 
          slentol = sr ;
          sidetol = kRMax ;
        }            
        if ( (ri < 0) || (sr < halfRadTolerance) )
        {
          // Safety: if both roots -ve ensure that sr cannot `win'
          //         distance to out

          sr2 = -b + sqrt(d) ;
          zi  = p.z + sr2*v.z ;
          ri  = tanRMax*zi + rMaxAv ;

          if ((ri >= 0) && (sr2 > halfRadTolerance))
          {
            sr = sr2;
          }
          else
          {
            sr = kInfinity ;

            if( (-halfRadTolerance <= sr2) && ( sr2 <= halfRadTolerance) )
            {
              // An intersection within the tolerance.
              // Storing it in case it is good.

              slentol = sr2 ;
              sidetol = kRMax ;
            }
          }
        }
      }
    }
    else
    {
      // No intersection with outer cone & not parallel
      // -> already outside, no intersection

      if ( calcNorm )
      {
        risec      = sqrt(t3)*secRMax;
        *validNorm = true;
        *n         = G4ThreeVector_create(p.x/risec,p.y/risec,-tanRMax/secRMax);
      }
      return snxt = 0.0 ;
    }
  }
  else if ( nt2 && (deltaRoi2 > 0.0) )
  {
    // Linear case (only one intersection) => point outside outer cone

    if ( calcNorm )
    {
      risec      = sqrt(t3)*secRMax;
      *validNorm = true;
      *n         = G4ThreeVector_create(p.x/risec,p.y/risec,-tanRMax/secRMax);
    }
    return snxt = 0.0 ;
  }
  else
  {
    // No intersection -> parallel to outer cone
    // => Z or inner cone intersection

    sr = kInfinity ;
  }

  // Check possible intersection within tolerance

  if ( slentol <= halfCarTolerance )
  {
    // An intersection within the tolerance was found.  
    // We must accept it only if the momentum points outwards.  
    //
    // G4ThreeVector ptTol ;  // The point of the intersection  
    // ptTol= p + slentol*v ;
    // ri=tanRMax*zi+rMaxAv ;
    //
    // Calculate a normal vector,  as below

    xi    = p.x + slentol*v.x;
    yi    = p.y + slentol*v.y;
    risec = sqrt(xi*xi + yi*yi)*secRMax;
    G4ThreeVector Normal = G4ThreeVector_create(xi/risec,yi/risec,-tanRMax/secRMax);

    if ( G4ThreeVector_dot(Normal,v) > 0 )    // We will leave the Cone immediatelly
    {
      if ( calcNorm ) 
      {
        *n         = G4ThreeVector_unit(Normal);
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
    else // On the surface, but not heading out so we ignore this intersection
    {    //                                        (as it is within tolerance).
      slentol = kInfinity ;
    }
  }

  // Inner Cone intersection

  if ( This->fRmin1 || This->fRmin2 )
  {
    tanRMin = (This->fRmin2 - This->fRmin1)*0.5/This->fDz ;
    nt1     = t1 - (tanRMin*v.z)*(tanRMin*v.z) ;

    if ( nt1 )
    {
      secRMin = sqrt(1.0 + tanRMin*tanRMin) ;
      rMinAv  = (This->fRmin1 + This->fRmin2)*0.5 ;    
      rin     = tanRMin*p.z + rMinAv ;
      nt2     = t2 - tanRMin*v.z*rin ;
      nt3     = t3 - rin*rin ;
      
      // Equation quadratic => 2 roots : first root must be leaving

      b = nt2/nt1 ;
      c = nt3/nt1 ;
      d = b*b - c ;

      if ( d >= 0.0 )
      {
        // NOTE: should be rho-rin<K_GEOMETRY_RAD_TOLERANCE*0.5,
        //       but using squared versions for efficiency

        if (nt3 < K_GEOMETRY_RAD_TOLERANCE*(rin + K_GEOMETRY_RAD_TOLERANCE*0.25)) 
        {
          if ( nt2 < 0.0 )
          {
            if (calcNorm)  { *validNorm = false; }
            return          snxt      = 0.0;
          }
        }
        else
        {
          sr2 = -b - sqrt(d) ;
          zi  = p.z + sr2*v.z ;
          ri  = tanRMin*zi + rMinAv ;

          if( (ri>=0.0)&&(-halfRadTolerance<=sr2)&&(sr2<=halfRadTolerance) )
          {
            // An intersection within the tolerance
            // storing it in case it is good.

            slentol = sr2 ;
            sidetol = kRMax ;
          }
          if( (ri<0) || (sr2 < halfRadTolerance) )
          {
            sr3 = -b + sqrt(d) ;

            // Safety: if both roots -ve ensure that sr cannot `win'
            //         distancetoout

            if  ( sr3 > halfRadTolerance )
            {
              if( sr3 < sr )
              {
                zi = p.z + sr3*v.z ;
                ri = tanRMin*zi + rMinAv ;

                if ( ri >= 0.0 )
                {
                  sr=sr3 ;
                  sider=kRMin ;
                }
              } 
            }
            else if ( sr3 > -halfRadTolerance )
            {
              // Intersection in tolerance. Store to check if it's good

              slentol = sr3 ;
              sidetol = kRMin ;
            }
          }
          else if ( (sr2 < sr) && (sr2 > halfCarTolerance) )
          {
            sr    = sr2 ;
            sider = kRMin ;
          }
          else if (sr2 > -halfCarTolerance)
          {
            // Intersection in tolerance. Store to check if it's good

            slentol = sr2 ;
            sidetol = kRMin ;
          }    
          if( slentol <= halfCarTolerance  )
          {
            // An intersection within the tolerance was found. 
            // We must accept it only if  the momentum points outwards. 

            G4ThreeVector Normal ; 
            
            // Calculate a normal vector,  as below

            xi     = p.x + slentol*v.x ;
            yi     = p.y + slentol*v.y ;
            if( sidetol==kRMax )
            {
              risec  = sqrt(xi*xi + yi*yi)*secRMax ;
              Normal = G4ThreeVector_create(xi/risec,yi/risec,-tanRMax/secRMax) ;
            }
            else
            {
              risec  = sqrt(xi*xi + yi*yi)*secRMin ;
              Normal = G4ThreeVector_create(-xi/risec,-yi/risec,tanRMin/secRMin) ;
            }
            if( G4ThreeVector_dot(Normal,v) > 0 )
            {
              // We will leave the cone immediately

              if( calcNorm ) 
              {
                *n         = G4ThreeVector_unit(Normal) ;
                *validNorm = true ;
              }
              return snxt = 0.0 ;
            }
            else 
            { 
              // On the surface, but not heading out so we ignore this
              // intersection (as it is within tolerance). 

              slentol = kInfinity ;
            }        
          }
        }
      }
    }
  }

  // Linear case => point outside inner cone ---> outer cone intersect
  //
  // Phi Intersection
  
#ifdef ENABLE_SLICED_CONS
  if ( !This->fPhiFullCone )
  {
    // add angle calculation with correction 
    // of the difference in domain of atan2 and Sphi

    vphi = atan2(v.y,v.x) ;

    if ( vphi < This->fSPhi - halfAngTolerance  )              { vphi += twopi; }
    else if ( vphi > This->fSPhi + This->fDPhi + halfAngTolerance )  { vphi -= twopi; }

    if ( p.x || p.y )   // Check if on z axis (rho not needed later)
    {
      // pDist -ve when inside

      pDistS = p.x*GET_sinSPhi - p.y*GET_cosSPhi ;
      pDistE = -p.x*GET_sinEPhi + p.y*GET_cosEPhi ;

      // Comp -ve when in direction of outwards normal

      compS = -GET_sinSPhi*v.x + GET_cosSPhi*v.y ;
      compE = GET_sinEPhi*v.x - GET_cosEPhi*v.y ;

      sidephi = kNull ;

      if( ( (This->fDPhi <= pi) && ( (pDistS <= halfCarTolerance)
                            && (pDistE <= halfCarTolerance) ) )
         || ( (This->fDPhi >  pi) && !((pDistS >  halfCarTolerance)
                              && (pDistE >  halfCarTolerance) ) )  )
      {
        // Inside both phi *full* planes
        if ( compS < 0 )
        {
          sphi = pDistS/compS ;
          if (sphi >= -halfCarTolerance)
          {
            xi = p.x + sphi*v.x ;
            yi = p.y + sphi*v.y ;

            // Check intersecting with correct half-plane
            // (if not -> no intersect)
            //
            if ( (fabs(xi)<=K_GEOMETRY_CAR_TOLERANCE)
              && (fabs(yi)<=K_GEOMETRY_CAR_TOLERANCE) )
            {
              sidephi= kSPhi;
              if ( ( This->fSPhi-halfAngTolerance <= vphi )
                && ( This->fSPhi+This->fDPhi+halfAngTolerance >=vphi ) )
              {
                sphi = kInfinity;
              }
            }
            else
            if ( (yi*GET_cosCPhi-xi*GET_sinCPhi)>=0 )
            {
              sphi = kInfinity ;
            }
            else
            {
              sidephi = kSPhi ;
              if ( pDistS > -halfCarTolerance )
              {
                sphi = 0.0 ; // Leave by sphi immediately
              }    
            }       
          }
          else
          {
            sphi = kInfinity ;
          }
        }
        else
        {
          sphi = kInfinity ;
        }

        if ( compE < 0 )
        {
          sphi2 = pDistE/compE ;

          // Only check further if < starting phi intersection
          //
          if ( (sphi2 > -halfCarTolerance) && (sphi2 < sphi) )
          {
            xi = p.x + sphi2*v.x ;
            yi = p.y + sphi2*v.y ;

            // Check intersecting with correct half-plane

            if ( (fabs(xi)<=K_GEOMETRY_CAR_TOLERANCE)
              && (fabs(yi)<=K_GEOMETRY_CAR_TOLERANCE) )
            {
              // Leaving via ending phi

              if(!( (This->fSPhi-halfAngTolerance <= vphi)
                 && (This->fSPhi+This->fDPhi+halfAngTolerance >= vphi) ) )
              {
                sidephi = kEPhi ;
                if ( pDistE <= -halfCarTolerance )  { sphi = sphi2; }
                else                                { sphi = 0.0; }
              }
            }
            else // Check intersecting with correct half-plane
            if ( yi*GET_cosCPhi-xi*GET_sinCPhi >= 0 )
            {
              // Leaving via ending phi

              sidephi = kEPhi ;
              if ( pDistE <= -halfCarTolerance )  { sphi = sphi2; }
              else                                { sphi = 0.0; }
            }
          }
        }
      }
      else
      {
        sphi = kInfinity ;
      }
    }
    else
    {
      // On z axis + travel not || to z axis -> if phi of vector direction
      // within phi of shape, Step limited by rmax, else Step =0

      if ( (This->fSPhi-halfAngTolerance <= vphi)
        && (vphi <= This->fSPhi+This->fDPhi+halfAngTolerance) )
      {
        sphi = kInfinity ;
      }
      else
      {
        sidephi = kSPhi  ;   // arbitrary 
        sphi    = 0.0 ;
      }
    }      
    if ( sphi < snxt )  // Order intersecttions
    {
      snxt=sphi ;
      side=sidephi ;
    }
  }
#endif
  if ( sr < snxt )  // Order intersections
  {
    snxt = sr    ;
    side = sider ;
  }
  if (calcNorm)
  {
    switch(side)
    {                     // Note: returned vector not normalised
      case kRMax:         // (divide by frmax for unit vector)
        xi         = p.x + snxt*v.x ;
        yi         = p.y + snxt*v.y ;
        risec      = sqrt(xi*xi + yi*yi)*secRMax ;
        *n         = G4ThreeVector_create(xi/risec,yi/risec,-tanRMax/secRMax) ;
        *validNorm = true ;
        break ;
      case kRMin:
        *validNorm = false ;  // Rmin is inconvex
        break ;
#ifdef ENABLE_SLICED_CONS
      case kSPhi:
        if ( This->fDPhi <= pi )
        {
          *n         = G4ThreeVector_create(GET_sinSPhi, -GET_cosSPhi, 0);
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;
      case kEPhi:
        if ( This->fDPhi <= pi )
        {
          *n = G4ThreeVector_create(-GET_sinEPhi, GET_cosEPhi, 0);
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;
#endif
      case kPZ:
        *n         = G4ThreeVector_create(0,0,1) ;
        *validNorm = true ;
        break ;
      case kMZ:
        *n         = G4ThreeVector_create(0,0,-1) ;
        *validNorm = true ;
        break ;
      default:
        myAssert(false);
        break ;
    }
  }
  if (snxt < halfCarTolerance)  { snxt = 0.; }

  return snxt ;
}

//////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

SOLIDINLINE G4double G4Cons_DistanceToOut(GEOMETRYLOC const G4Cons *This, G4ThreeVector p)
{
  G4double safe=0.0, rho, safeR1, safeR2, safeZ;
  G4double tanRMin, secRMin, pRMin;
  G4double tanRMax, secRMax, pRMax;

  rho = sqrt(p.x*p.x + p.y*p.y) ;
  safeZ = This->fDz - fabs(p.z) ;

  if (This->fRmin1 || This->fRmin2)
  {
    tanRMin = (This->fRmin2 - This->fRmin1)*0.5/This->fDz ;
    secRMin = sqrt(1.0 + tanRMin*tanRMin) ;
    pRMin   = tanRMin*p.z + (This->fRmin1 + This->fRmin2)*0.5 ;
    safeR1  = (rho - pRMin)/secRMin ;
  }
  else
  {
    safeR1 = kInfinity ;
  }

  tanRMax = (This->fRmax2 - This->fRmax1)*0.5/This->fDz ;
  secRMax = sqrt(1.0 + tanRMax*tanRMax) ;
  pRMax   = tanRMax*p.z + (This->fRmax1+This->fRmax2)*0.5 ;
  safeR2  = (pRMax - rho)/secRMax ;

  if (safeR1 < safeR2)  { safe = safeR1; }
  else                  { safe = safeR2; }
  if (safeZ < safe)     { safe = safeZ ; }

  // Check if phi divided, Calc distances closest phi plane

#ifdef ENABLE_SLICED_CONS
  if (!This->fPhiFullCone)
  {
    // Above/below central phi of G4Cons
    
    G4double safePhi;

    if ( (p.y*GET_cosCPhi - p.x*GET_sinCPhi) <= 0 )
    {
      safePhi = -(p.x*GET_sinSPhi - p.y*GET_cosSPhi) ;
    }
    else
    {
      safePhi = (p.x*GET_sinEPhi - p.y*GET_cosEPhi) ;
    }
    if (safePhi < safe)  { safe = safePhi; }
  }
#endif
  if ( safe < 0 )  { safe = 0; }

  return safe ;
}


#ifdef ENABLE_G4POLYCONE

#define SQ(x) ((x)*(x))

///////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from 'inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

SOLIDINLINE G4double G4PolyConeCons_DistanceToOut_full(
			   const G4Cons *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n,
			   G4bool end1, G4bool end2)
{
  if ( This->fDz < K_GEOMETRY_CAR_TOLERANCE )
  {
	const G4double r2 = p.x*p.x + p.y*p.y;
	
	if ( ( r2 >= SQ(This->fRmin1) || r2 >= SQ(This->fRmin2) ) &&
		 ( r2 <= SQ(This->fRmax1) || r2 <= SQ(This->fRmax2) ) )
	{
		if ( v.z < 0 )
		{
			if ( (r2 <= SQ(This->fRmin1) && r2 >= SQ(This->fRmin2)) ||
				 (r2 >= SQ(This->fRmax1) && r2 <= SQ(This->fRmax2) ) || end1 )
			{
				if ( calcNorm )
				{
					*validNorm = true;
					*n = G4ThreeVector_create(0,0,-1);
				}
				return 0;
			}
		}
		else
		{
			if ( (r2 >= SQ(This->fRmin1) && r2 <= SQ(This->fRmin2)) ||
				 (r2 <= SQ(This->fRmax1) && r2 >= SQ(This->fRmax2) ) || end2 )
			{
				if ( calcNorm )
				{
					*validNorm = true;
					*n = G4ThreeVector_create(0,0,1);
				}
				return 0;
			}
		}
	}
	
	return kInfinity;
  }
  else
  {
	
  typedef enum {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ} ESide;
  ESide side = kNull, sider = kNull, sidetol = kNull ;

  const G4double halfCarTolerance=K_GEOMETRY_CAR_TOLERANCE*0.5;
  const G4double halfRadTolerance=K_GEOMETRY_RAD_TOLERANCE*0.5;
  
#ifdef ENABLE_SLICED_POLYCONS
  const G4double halfAngTolerance=K_GEOMETRY_ANG_TOLERANCE*0.5;
  const G4double pi = M_PI;
  
  G4double sphi, pDistS, compS, pDistE, compE, sphi2, vphi;
  ESide sidephi = kNull;
#endif

  G4double snxt,sr,pdist ;

  G4double tanRMax, secRMax, rMaxAv ;  // Data for outer cone
  G4double tanRMin, secRMin, rMinAv ;  // Data for inner cone

  G4double t1, t2, t3, rout, rin, nt1, nt2, nt3 ;
  G4double b, c, d, sr2, sr3 ;

  // Vars for intersection within tolerance

  G4double slentol = kInfinity ;

  // Vars for phi intersection:

  G4double xi, yi, risec ;
  G4double zi, ri, deltaRoi2 ;

  // Z plane intersection

  if ( v.z > 0.0 )
  {
    pdist = This->fDz - p.z ;

    if (pdist > halfCarTolerance)
    {
      snxt = pdist/v.z ;
      side = kPZ ;
    }
    else
    {
      if (calcNorm)
      {
        *n         = G4ThreeVector_create(0,0,1) ;
        *validNorm = true ;
      }
      if (!end2) return kInfinity;
      return  snxt = 0.0;
    }
  }
  else if ( v.z < 0.0 )
  {
    pdist = This->fDz + p.z ;

    if ( pdist > halfCarTolerance)
    {
      snxt = -pdist/v.z ;
      side = kMZ ;
    }
    else
    {
      if ( calcNorm )
      {
        *n         = G4ThreeVector_create(0,0,-1) ;
        *validNorm = true ;
      }
      if (!end1) return kInfinity;
      return snxt = 0.0 ;
    }
  }
  else     // Travel perpendicular to z axis
  {
    snxt = kInfinity ;    
    side = kNull ;
  }

  // Radial Intersections
  //
  // Intersection with outer cone (possible return) and
  //                   inner cone (must also check phi)
  //
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=(a*z+b)^2
  //
  // where a=tanRMax or tanRMin
  //       b=rMaxAv  or rMinAv
  //
  // (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0 ;
  //     t1                        t2                      t3  
  //
  //  \--------u-------/       \-----------v----------/ \---------w--------/

  tanRMax = (This->fRmax2 - This->fRmax1)*0.5/This->fDz ;
  secRMax = sqrt(1.0 + tanRMax*tanRMax) ;
  rMaxAv  = (This->fRmax1 + This->fRmax2)*0.5 ;


  t1   = 1.0 - v.z*v.z ;      // since v normalised
  t2   = p.x*v.x + p.y*v.y ;
  t3   = p.x*p.x + p.y*p.y ;
  rout = tanRMax*p.z + rMaxAv ;

  nt1 = t1 - (tanRMax*v.z)*(tanRMax*v.z) ;
  nt2 = t2 - tanRMax*v.z*rout ;
  nt3 = t3 - rout*rout ;

  if (v.z > 0.0)
  {
    deltaRoi2 = snxt*snxt*t1 + 2*snxt*t2 + t3
                - This->fRmax2*(This->fRmax2 + K_GEOMETRY_RAD_TOLERANCE*secRMax);
  }
  else if ( v.z < 0.0 )
  {
    deltaRoi2 = snxt*snxt*t1 + 2*snxt*t2 + t3
                - This->fRmax1*(This->fRmax1 + K_GEOMETRY_RAD_TOLERANCE*secRMax);
  }
  else
  {
    deltaRoi2 = 1.0;
  }

  if ( nt1 && (deltaRoi2 > 0.0) )  
  {
    // Equation quadratic => 2 roots : second root must be leaving

    b = nt2/nt1 ;
    c = nt3/nt1 ;
    d = b*b - c ;

    if ( d >= 0 )
    {
      // Check if on outer cone & heading outwards
      // NOTE: Should use rho-rout>-K_GEOMETRY_RAD_TOLERANCE*0.5
        
      if (nt3 > -halfRadTolerance && nt2 >= 0 )
      {
        if (calcNorm)
        {
          risec      = sqrt(t3)*secRMax ;
          *validNorm = true ;
          *n         = G4ThreeVector_create(p.x/risec,p.y/risec,-tanRMax/secRMax);
        }
        return snxt=0 ;
      }
      else
      {
        sider = kRMax  ;
        sr    = -b - sqrt(d) ; // was +srqrt(d), vmg 28.04.99
        zi    = p.z + sr*v.z ;
        ri    = tanRMax*zi + rMaxAv ;
          
        if ((ri >= 0) && (-halfRadTolerance <= sr) && (sr <= halfRadTolerance))
        {
          // An intersection within the tolerance
          //   we will Store it in case it is good -
          // 
          slentol = sr ;
          sidetol = kRMax ;
        }            
        if ( (ri < 0) || (sr < halfRadTolerance) )
        {
          // Safety: if both roots -ve ensure that sr cannot `win'
          //         distance to out

          sr2 = -b + sqrt(d) ;
          zi  = p.z + sr2*v.z ;
          ri  = tanRMax*zi + rMaxAv ;

          if ((ri >= 0) && (sr2 > halfRadTolerance))
          {
            sr = sr2;
          }
          else
          {
            sr = kInfinity ;

            if( (-halfRadTolerance <= sr2) && ( sr2 <= halfRadTolerance) )
            {
              // An intersection within the tolerance.
              // Storing it in case it is good.

              slentol = sr2 ;
              sidetol = kRMax ;
            }
          }
        }
      }
    }
    else
    {
      // No intersection with outer cone & not parallel
      // -> already outside, no intersection

      if ( calcNorm )
      {
        risec      = sqrt(t3)*secRMax;
        *validNorm = true;
        *n         = G4ThreeVector_create(p.x/risec,p.y/risec,-tanRMax/secRMax);
      }
      return snxt = 0.0 ;
    }
  }
  else if ( nt2 && (deltaRoi2 > 0.0) )
  {
    // Linear case (only one intersection) => point outside outer cone

    if ( calcNorm )
    {
      risec      = sqrt(t3)*secRMax;
      *validNorm = true;
      *n         = G4ThreeVector_create(p.x/risec,p.y/risec,-tanRMax/secRMax);
    }
    return snxt = 0.0 ;
  }
  else
  {
    // No intersection -> parallel to outer cone
    // => Z or inner cone intersection

    sr = kInfinity ;
  }

  // Check possible intersection within tolerance

  if ( slentol <= halfCarTolerance )
  {
    // An intersection within the tolerance was found.  
    // We must accept it only if the momentum points outwards.  
    //
    // G4ThreeVector ptTol ;  // The point of the intersection  
    // ptTol= p + slentol*v ;
    // ri=tanRMax*zi+rMaxAv ;
    //
    // Calculate a normal vector,  as below

    xi    = p.x + slentol*v.x;
    yi    = p.y + slentol*v.y;
    risec = sqrt(xi*xi + yi*yi)*secRMax;
    G4ThreeVector Normal = G4ThreeVector_create(xi/risec,yi/risec,-tanRMax/secRMax);

    if ( G4ThreeVector_dot(Normal,v) > 0 )    // We will leave the Cone immediatelly
    {
      if ( calcNorm ) 
      {
        *n         = G4ThreeVector_unit(Normal);
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
    else // On the surface, but not heading out so we ignore this intersection
    {    //                                        (as it is within tolerance).
      slentol = kInfinity ;
    }
  }

  // Inner Cone intersection

  if ( This->fRmin1 || This->fRmin2 )
  {
    tanRMin = (This->fRmin2 - This->fRmin1)*0.5/This->fDz ;
    nt1     = t1 - (tanRMin*v.z)*(tanRMin*v.z) ;

    if ( nt1 )
    {
      secRMin = sqrt(1.0 + tanRMin*tanRMin) ;
      rMinAv  = (This->fRmin1 + This->fRmin2)*0.5 ;    
      rin     = tanRMin*p.z + rMinAv ;
      nt2     = t2 - tanRMin*v.z*rin ;
      nt3     = t3 - rin*rin ;
      
      // Equation quadratic => 2 roots : first root must be leaving

      b = nt2/nt1 ;
      c = nt3/nt1 ;
      d = b*b - c ;

      if ( d >= 0.0 )
      {
        // NOTE: should be rho-rin<K_GEOMETRY_RAD_TOLERANCE*0.5,
        //       but using squared versions for efficiency

        if (nt3 < K_GEOMETRY_RAD_TOLERANCE*(rin + K_GEOMETRY_RAD_TOLERANCE*0.25)) 
        {
          if ( nt2 < 0.0 )
          {
            if (calcNorm)  { *validNorm = false; }
            return          snxt      = 0.0;
          }
        }
        else
        {
          sr2 = -b - sqrt(d) ;
          zi  = p.z + sr2*v.z ;
          ri  = tanRMin*zi + rMinAv ;

          if( (ri>=0.0)&&(-halfRadTolerance<=sr2)&&(sr2<=halfRadTolerance) )
          {
            // An intersection within the tolerance
            // storing it in case it is good.

            slentol = sr2 ;
            sidetol = kRMax ;
          }
          if( (ri<0) || (sr2 < halfRadTolerance) )
          {
            sr3 = -b + sqrt(d) ;

            // Safety: if both roots -ve ensure that sr cannot `win'
            //         distancetoout

            if  ( sr3 > halfRadTolerance )
            {
              if( sr3 < sr )
              {
                zi = p.z + sr3*v.z ;
                ri = tanRMin*zi + rMinAv ;

                if ( ri >= 0.0 )
                {
                  sr=sr3 ;
                  sider=kRMin ;
                }
              } 
            }
            else if ( sr3 > -halfRadTolerance )
            {
              // Intersection in tolerance. Store to check if it's good

              slentol = sr3 ;
              sidetol = kRMin ;
            }
          }
          else if ( (sr2 < sr) && (sr2 > halfCarTolerance) )
          {
            sr    = sr2 ;
            sider = kRMin ;
          }
          else if (sr2 > -halfCarTolerance)
          {
            // Intersection in tolerance. Store to check if it's good

            slentol = sr2 ;
            sidetol = kRMin ;
          }    
          if( slentol <= halfCarTolerance  )
          {
            // An intersection within the tolerance was found. 
            // We must accept it only if  the momentum points outwards. 

            G4ThreeVector Normal ; 
            
            // Calculate a normal vector,  as below

            xi     = p.x + slentol*v.x ;
            yi     = p.y + slentol*v.y ;
            if( sidetol==kRMax )
            {
              risec  = sqrt(xi*xi + yi*yi)*secRMax ;
              Normal = G4ThreeVector_create(xi/risec,yi/risec,-tanRMax/secRMax) ;
            }
            else
            {
              risec  = sqrt(xi*xi + yi*yi)*secRMin ;
              Normal = G4ThreeVector_create(-xi/risec,-yi/risec,tanRMin/secRMin) ;
            }
            if( G4ThreeVector_dot(Normal,v) > 0 )
            {
              // We will leave the cone immediately

              if( calcNorm ) 
              {
                *n         = G4ThreeVector_unit(Normal) ;
                *validNorm = true ;
              }
              return snxt = 0.0 ;
            }
            else 
            { 
              // On the surface, but not heading out so we ignore this
              // intersection (as it is within tolerance). 

              slentol = kInfinity ;
            }        
          }
        }
      }
    }
  }

  // Linear case => point outside inner cone ---> outer cone intersect
  //
  // Phi Intersection
  
#ifdef ENABLE_SLICED_POLYCONS
  if ( !This->fPhiFullCone )
  {
    // add angle calculation with correction 
    // of the difference in domain of atan2 and Sphi

    vphi = atan2(v.y,v.x) ;

    if ( vphi < This->fSPhi - halfAngTolerance  )              { vphi += twopi; }
    else if ( vphi > This->fSPhi + This->fDPhi + halfAngTolerance )  { vphi -= twopi; }

    if ( p.x || p.y )   // Check if on z axis (rho not needed later)
    {
      // pDist -ve when inside

      pDistS = p.x*GET_sinSPhi - p.y*GET_cosSPhi ;
      pDistE = -p.x*GET_sinEPhi + p.y*GET_cosEPhi ;

      // Comp -ve when in direction of outwards normal

      compS = -GET_sinSPhi*v.x + GET_cosSPhi*v.y ;
      compE = GET_sinEPhi*v.x - GET_cosEPhi*v.y ;

      sidephi = kNull ;

      if( ( (This->fDPhi <= pi) && ( (pDistS <= halfCarTolerance)
                            && (pDistE <= halfCarTolerance) ) )
         || ( (This->fDPhi >  pi) && !((pDistS >  halfCarTolerance)
                              && (pDistE >  halfCarTolerance) ) )  )
      {
        // Inside both phi *full* planes
        if ( compS < 0 )
        {
          sphi = pDistS/compS ;
          if (sphi >= -halfCarTolerance)
          {
            xi = p.x + sphi*v.x ;
            yi = p.y + sphi*v.y ;

            // Check intersecting with correct half-plane
            // (if not -> no intersect)
            //
            if ( (fabs(xi)<=K_GEOMETRY_CAR_TOLERANCE)
              && (fabs(yi)<=K_GEOMETRY_CAR_TOLERANCE) )
            {
              sidephi= kSPhi;
              if ( ( This->fSPhi-halfAngTolerance <= vphi )
                && ( This->fSPhi+This->fDPhi+halfAngTolerance >=vphi ) )
              {
                sphi = kInfinity;
              }
            }
            else
            if ( (yi*GET_cosCPhi-xi*GET_sinCPhi)>=0 )
            {
              sphi = kInfinity ;
            }
            else
            {
              sidephi = kSPhi ;
              if ( pDistS > -halfCarTolerance )
              {
                sphi = 0.0 ; // Leave by sphi immediately
              }    
            }       
          }
          else
          {
            sphi = kInfinity ;
          }
        }
        else
        {
          sphi = kInfinity ;
        }

        if ( compE < 0 )
        {
          sphi2 = pDistE/compE ;

          // Only check further if < starting phi intersection
          //
          if ( (sphi2 > -halfCarTolerance) && (sphi2 < sphi) )
          {
            xi = p.x + sphi2*v.x ;
            yi = p.y + sphi2*v.y ;

            // Check intersecting with correct half-plane

            if ( (fabs(xi)<=K_GEOMETRY_CAR_TOLERANCE)
              && (fabs(yi)<=K_GEOMETRY_CAR_TOLERANCE) )
            {
              // Leaving via ending phi

              if(!( (This->fSPhi-halfAngTolerance <= vphi)
                 && (This->fSPhi+This->fDPhi+halfAngTolerance >= vphi) ) )
              {
                sidephi = kEPhi ;
                if ( pDistE <= -halfCarTolerance )  { sphi = sphi2; }
                else                                { sphi = 0.0; }
              }
            }
            else // Check intersecting with correct half-plane
            if ( yi*GET_cosCPhi-xi*GET_sinCPhi >= 0 )
            {
              // Leaving via ending phi

              sidephi = kEPhi ;
              if ( pDistE <= -halfCarTolerance )  { sphi = sphi2; }
              else                                { sphi = 0.0; }
            }
          }
        }
      }
      else
      {
        sphi = kInfinity ;
      }
    }
    else
    {
      // On z axis + travel not || to z axis -> if phi of vector direction
      // within phi of shape, Step limited by rmax, else Step =0

      if ( (This->fSPhi-halfAngTolerance <= vphi)
        && (vphi <= This->fSPhi+This->fDPhi+halfAngTolerance) )
      {
        sphi = kInfinity ;
      }
      else
      {
        sidephi = kSPhi  ;   // arbitrary 
        sphi    = 0.0 ;
      }
    }      
    if ( sphi < snxt )  // Order intersecttions
    {
      snxt=sphi ;
      side=sidephi ;
    }
  }
#endif
  if ( sr < snxt )  // Order intersections
  {
    snxt = sr    ;
    side = sider ;
  }
  //if (calcNorm)
  {
    switch(side)
    {                     // Note: returned vector not normalised
      case kRMax:         // (divide by frmax for unit vector)
        xi         = p.x + snxt*v.x ;
        yi         = p.y + snxt*v.y ;
        risec      = sqrt(xi*xi + yi*yi)*secRMax ;
        *n         = G4ThreeVector_create(xi/risec,yi/risec,-tanRMax/secRMax) ;
        *validNorm = true ;
        break ;
      case kRMin:
        *validNorm = false ;  // Rmin is inconvex
        break ;
#ifdef ENABLE_SLICED_POLYCONS
      case kSPhi:
        if ( This->fDPhi <= pi )
        {
          *n         = G4ThreeVector_create(GET_sinSPhi, -GET_cosSPhi, 0);
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;
      case kEPhi:
        if ( This->fDPhi <= pi )
        {
          *n = G4ThreeVector_create(-GET_sinEPhi, GET_cosEPhi, 0);
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;
#endif
      case kPZ:
        *n         = G4ThreeVector_create(0,0,1) ;
        *validNorm = true ;
        break ;
      case kMZ:
        *n         = G4ThreeVector_create(0,0,-1) ;
        *validNorm = true ;
        break ;
      default:
        myAssert(false);
        break ;
    }
  }
  if (snxt < halfCarTolerance)  { snxt = 0.; }

  if ( (side == kPZ && !end2) || (side == kMZ && !end1) )
  {
	if (calcNorm) *validNorm = false;
	return kInfinity;
  }
  return snxt ;
  
  }
}

/////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

SOLIDINLINE EInside G4PolyConeCons_Inside(GEOMETRYLOC const G4Cons *This, G4ThreeVector p, G4bool end1, G4bool end2)
{
  const G4double r2 = p.x*p.x + p.y*p.y;
		
  if ( This->fDz < K_GEOMETRY_CAR_TOLERANCE )
  {
		if (( r2 >= SQ(This->fRmax2) &&
			  r2 >= SQ(This->fRmax1) ) || 
		    ( r2 <= SQ(This->fRmin1) &&
			  r2 <= SQ(This->fRmin2) )) return kOutside;
			  
		if ( r2 >= SQ(This->fRmin1) && r2 >= SQ(This->fRmin2 ) &&
			 r2 <= SQ(This->fRmax1) && r2 <= SQ(This->fRmax2 ) )
		{
			if ( end1 || end2 ) return kSurface;
			else return kInside;
		}
			
		return kSurface;
  }
  else
  {
	
  G4double rl, rh, tolRMin, tolRMax; // rh2, rl2 ;
  EInside in;
  const G4double halfCarTolerance=K_GEOMETRY_CAR_TOLERANCE*0.5;
  const G4double halfRadTolerance=K_GEOMETRY_RAD_TOLERANCE*0.5;
#ifdef ENABLE_SLICED_CONS
  const G4double halfAngTolerance=K_GEOMETRY_ANG_TOLERANCE*0.5;
  G4double pPhi;
#endif

  if (fabs(p.z) > This->fDz + halfCarTolerance )  { return in = kOutside; }
  else if(fabs(p.z) >= This->fDz - halfCarTolerance )
  {
	  if ( (p.z < 0 && end1) || (p.z > 0 && end2) ) return in = kSurface;
	  else return in = kInside;
  }
  else in = kInside;

  rl = 0.5*(This->fRmin2*(p.z + This->fDz) + This->fRmin1*(This->fDz - p.z))/This->fDz ;
  rh = 0.5*(This->fRmax2*(p.z+This->fDz)+This->fRmax1*(This->fDz-p.z))/This->fDz;

  // rh2 = rh*rh;

  tolRMin = rl - halfRadTolerance;
  if ( tolRMin < 0 )  { tolRMin = 0; }
  tolRMax = rh + halfRadTolerance;

  if ( (r2<tolRMin*tolRMin) || (r2>tolRMax*tolRMax) ) { return in = kOutside; }

  if (rl) { tolRMin = rl + halfRadTolerance; }
  else    { tolRMin = 0.0; }
  tolRMax = rh - halfRadTolerance;
      
  if (in == kInside) // else it's kSurface already
  {
     if ( (r2 < tolRMin*tolRMin) || (r2 >= tolRMax*tolRMax) ) { in = kSurface; }
  }
  
#ifdef ENABLE_SLICED_CONS
  if ( !This->fPhiFullCone && ((p.x != 0.0) || (p.y != 0.0)) )
  {
    pPhi = atan2(p.y,p.x) ;

    if ( pPhi < This->fSPhi - halfAngTolerance  )             { pPhi += twopi; }
    else if ( pPhi > This->fSPhi + This->fDPhi + halfAngTolerance ) { pPhi -= twopi; }
    
    if ( (pPhi < This->fSPhi - halfAngTolerance) ||          
         (pPhi > This->fSPhi + This->fDPhi + halfAngTolerance) )  { return in = kOutside; }
      
    else if (in == kInside)  // else it's kSurface anyway already
    {
       if ( (pPhi < This->fSPhi + halfAngTolerance) || 
            (pPhi > This->fSPhi + This->fDPhi - halfAngTolerance) )  { in = kSurface; }
    }
  }
  else if ( !This->fPhiFullCone ) { in = kSurface; }
#endif

  return in ;
  
  }
}

//////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

/*SOLIDINLINE G4double G4PolyConeCons_DistanceToOut(
	const G4Cons *This, G4ThreeVector p, G4bool end1, G4bool end2)
{
  G4double safe=0.0, rho, safeR1, safeR2, safeZ = kInfinity;
  G4double tanRMin, secRMin, pRMin;
  G4double tanRMax, secRMax, pRMax;

  rho = sqrt(p.x*p.x + p.y*p.y) ;

  if (This->fRmin1 || This->fRmin2)
  {
    tanRMin = (This->fRmin2 - This->fRmin1)*0.5/This->fDz ;
    secRMin = sqrt(1.0 + tanRMin*tanRMin) ;
    pRMin   = tanRMin*p.z + (This->fRmin1 + This->fRmin2)*0.5 ;
    safeR1  = (rho - pRMin)/secRMin ;
  }
  else
  {
    safeR1 = kInfinity ;
  }

  tanRMax = (This->fRmax2 - This->fRmax1)*0.5/This->fDz ;
  secRMax = sqrt(1.0 + tanRMax*tanRMax) ;
  pRMax   = tanRMax*p.z + (This->fRmax1+This->fRmax2)*0.5 ;
  safeR2  = (pRMax - rho)/secRMax ;

  if (safeR1 < safeR2)  { safe = safeR1; }
  else                  { safe = safeR2; }
  
  if ( end1 )
  {
	  safeZ = p.z - (-This->fDz);
	  if ( safeZ < safe ) safe = safeZ;
  }
  
  if ( end2 )
  {
	  safeZ = This->fDz - p.z;
	  if ( safeZ < safe ) safe = safeZ;
  }

  // Check if phi divided, Calc distances closest phi plane

#ifdef ENABLE_SLICED_POLYCONS
  if (!This->fPhiFullCone)
  {
    // Above/below central phi of G4Cons
    
    G4double safePhi;

    if ( (p.y*GET_cosCPhi - p.x*GET_sinCPhi) <= 0 )
    {
      safePhi = -(p.x*GET_sinSPhi - p.y*GET_cosSPhi) ;
    }
    else
    {
      safePhi = (p.x*GET_sinEPhi - p.y*GET_cosEPhi) ;
    }
    if (safePhi < safe)  { safe = safePhi; }
  }
#endif
  if ( safe < 0 )  { safe = 0; }

  return safe ;
}*/

#undef GET_sinCPhi
#undef GET_cosCPhi
#undef GET_cosHDPhiIT
#undef GET_cosHDPhiOT
#undef GET_sinSPhi
#undef GET_cosSPhi
#undef GET_sinEPhi
#undef GET_cosEPhi
#undef GET_hDPhi
#undef GET_cPhi
#undef GET_ePhi
  

#endif
