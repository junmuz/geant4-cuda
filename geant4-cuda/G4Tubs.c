
/**
 * G4Tubs implementation
 * based on G4Tubs.cc of Geant 4.9.3
 */

#include "G4Tubs.h"

#ifdef HOST_CODE

INLINEFUNC
G4double G4Tubs_GetInnerRadius ( const G4Tubs *This )
{
  return This->fRMin;
}

INLINEFUNC
G4double G4Tubs_GetOuterRadius ( const G4Tubs *This )
{
  return This->fRMax;
}

INLINEFUNC
G4double G4Tubs_GetZHalfLength ( const G4Tubs *This )
{
  return This->fDz;
}

#ifdef ENABLE_SLICED_TUBS
INLINEFUNC
G4double G4Tubs_GetStartPhiAngle ( const G4Tubs *This )
{
  return This->fSPhi;
}

INLINEFUNC
G4double G4Tubs_GetDeltaPhiAngle ( const G4Tubs *This )
{
  return This->fDPhi;
}

INLINEFUNC 
void G4Tubs_InitializeTrigonometry( G4Tubs *This )
{
  G4double hDPhi = 0.5*This->fDPhi;                       // half delta phi
  G4double cPhi  = This->fSPhi + hDPhi; 
  G4double ePhi  = This->fSPhi + This->fDPhi;

  This->sinCPhi    = std::sin(cPhi);
  This->cosCPhi    = std::cos(cPhi);
  This->cosHDPhiIT = std::cos(hDPhi - 0.5*K_GEOMETRY_ANG_TOLERANCE); // inner/outer tol half dphi
  This->cosHDPhiOT = std::cos(hDPhi + 0.5*K_GEOMETRY_ANG_TOLERANCE);
  This->sinSPhi = std::sin(This->fSPhi);
  This->cosSPhi = std::cos(This->fSPhi);
  This->sinEPhi = std::sin(ePhi);
  This->cosEPhi = std::cos(ePhi);
}

INLINEFUNC void G4Tubs_CheckSPhiAngle(G4Tubs *This, G4double sPhi)
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

INLINEFUNC void G4Tubs_CheckDPhiAngle(G4Tubs *This, G4double dPhi)
{
  This->fPhiFullTube = true;
  if ( dPhi >= twopi-K_GEOMETRY_ANG_TOLERANCE*0.5 )
  {
    This->fDPhi=twopi;
    This->fSPhi=0;
  }
  else
  {
    This->fPhiFullTube = false;
    This->fDPhi = dPhi;
    myAssert( dPhi > 0 );
  }
}

INLINEFUNC void G4Tubs_CheckPhiAngles(G4Tubs *This, G4double sPhi, G4double dPhi)
{
  G4Tubs_CheckDPhiAngle(This,dPhi);
  if ( (This->fDPhi<twopi) && (sPhi) ) { G4Tubs_CheckSPhiAngle(This,sPhi); }
  G4Tubs_InitializeTrigonometry(This);
}
#endif

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

SOLIDINLINE
void G4Tubs_ctor( G4Tubs *This, G4double pRMin, G4double pRMax,
                  G4double pDz, G4double pSPhi, G4double pDPhi )
{
	This->solid.type = kTubs;
	
#ifdef ENABLE_SLICED_TUBS
	This->fSPhi = 0;
	This->fDPhi = 0;
#else
	(void)pSPhi;
	(void)pDPhi;
#endif
	
	myAssert( pDz > 0 );
	
	This->fDz = pDz;
	
	myAssert( (pRMin < pRMax) && (pRMin >= 0) ); // Check radii
	This->fRMin = pRMin;
	This->fRMax = pRMax;

#ifdef ENABLE_SLICED_TUBS
	// Check angles
	G4Tubs_CheckPhiAngles(This, pSPhi, pDPhi);
#endif
}

#else // host code

/////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

INLINEFUNC
G4ThreeVector G4Tubs_ApproxSurfaceNormal( GEOMETRYLOC const G4Tubs *This, G4ThreeVector p )
{
  enum {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ} side;
  G4ThreeVector norm ;
  G4double rho;
  G4double distZ, distRMin, distRMax ;
#ifdef ENABLE_SLICED_TUBS
  G4double phi, distSPhi, distEPhi, distMin;
#endif

  rho = sqrt(p.x*p.x + p.y*p.y) ;

  distRMin = fabs(rho - This->fRMin) ;
  distRMax = fabs(rho - This->fRMax) ;
  distZ    = fabs(fabs(p.z) - This->fDz) ;

  if (distRMin < distRMax) // First minimum
  {
    if ( distZ < distRMin )
    {
#ifdef ENABLE_SLICED_TUBS
       distMin = distZ ;
#endif
       side    = kNZ ;
    }
    else
    {
#ifdef ENABLE_SLICED_TUBS
      distMin = distRMin ;
#endif
      side    = kNRMin   ;
    }
  }
  else
  {
    if ( distZ < distRMax )
    {
#ifdef ENABLE_SLICED_TUBS
      distMin = distZ ;
#endif
      side    = kNZ   ;
    }
    else
    {
#ifdef ENABLE_SLICED_TUBS
      distMin = distRMax ;
#endif
      side    = kNRMax   ;
    }
  }   
#ifdef ENABLE_SLICED_TUBS
  if (!This->fPhiFullTube  &&  rho ) // Protected against (0,0,z) 
  {
    phi = atan2(p.y,p.x) ;

    if ( phi < 0 )  { phi += twopi; }

    if ( This->fSPhi < 0 )
    {
      distSPhi = fabs(phi - (This->fSPhi + twopi))*rho ;
    }
    else
    {
      distSPhi = fabs(phi - This->fSPhi)*rho ;
    }
    distEPhi = fabs(phi - This->fSPhi - This->fDPhi)*rho ;
                                      
    if (distSPhi < distEPhi) // Find new minimum
    {
      if ( distSPhi < distMin )
      {
        side = kNSPhi ;
      }
    }
    else
    {
      if ( distEPhi < distMin )
      {
        side = kNEPhi ;
      }
    }
  }
#endif  
  switch ( side )
  {
    case kNRMin : // Inner radius
    {                      
      norm = G4ThreeVector_create(-p.x/rho, -p.y/rho, 0) ;
      break ;
    }
    case kNRMax : // Outer radius
    {                  
      norm = G4ThreeVector_create(p.x/rho, p.y/rho, 0) ;
      break ;
    }
    case kNZ : //    + or - dz
    {                              
      if ( p.z > 0 )  { norm = G4ThreeVector_create(0,0,1) ; }
      else              { norm = G4ThreeVector_create(0,0,-1); }
      break ;
    }
#ifdef ENABLE_SLICED_TUBS
    case kNSPhi:
    {
      norm = G4ThreeVector_create(sin(This->fSPhi), -cos(This->fSPhi), 0) ;
      break ;
    }
    case kNEPhi:
    {
      norm = G4ThreeVector_create(-sin(This->fSPhi+This->fDPhi), cos(This->fSPhi+This->fDPhi), 0) ;
      break;
    }
#endif
    default:
    {
      myAssert(false);
      break ;
    }    
  }                
  return norm;
}

///////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

SOLIDINLINE G4ThreeVector G4Tubs_SurfaceNormal(GEOMETRYLOC const G4Tubs *This, G4ThreeVector p)
{
  G4int noSurfaces = 0;
  G4double rho;
  G4double distZ, distRMin, distRMax;

  const G4double halfCarTolerance = 0.5*K_GEOMETRY_CAR_TOLERANCE;
  
#ifdef ENABLE_SLICED_TUBS
  const G4double halfAngTolerance = 0.5*K_GEOMETRY_ANG_TOLERANCE;
  G4double distSPhi = kInfinity, distEPhi = kInfinity, pPhi;
  G4ThreeVector nPs, nPe;
#endif

  G4ThreeVector norm, sumnorm = G4ThreeVector_create(0.,0.,0.);
  G4ThreeVector nZ = G4ThreeVector_create(0, 0, 1.0);
  G4ThreeVector nR;

  rho = sqrt(p.x*p.x + p.y*p.y);

  distRMin = fabs(rho - This->fRMin);
  distRMax = fabs(rho - This->fRMax);
  distZ    = fabs(fabs(p.z) - This->fDz);

#ifdef ENABLE_SLICED_TUBS
  if (!This->fPhiFullTube)    // Protected against (0,0,z) 
  {
    if ( rho > halfCarTolerance )
    {
      pPhi = atan2(p.y,p.x);
    
      if(pPhi  < This->fSPhi- halfCarTolerance)           { pPhi += twopi; }
      else if(pPhi > This->fSPhi+This->fDPhi+ halfCarTolerance) { pPhi -= twopi; }

      distSPhi = fabs(pPhi - This->fSPhi);       
      distEPhi = fabs(pPhi - This->fSPhi - This->fDPhi); 
    }
    else if( !This->fRMin )
    {
      distSPhi = 0.; 
      distEPhi = 0.; 
    }
    nPs = G4ThreeVector_create(sin(This->fSPhi),-cos(This->fSPhi),0);
    nPe = G4ThreeVector_create(-sin(This->fSPhi+This->fDPhi),cos(This->fSPhi+This->fDPhi),0);
  }
#endif
  if ( rho > halfCarTolerance ) { nR = G4ThreeVector_create(p.x/rho,p.y/rho,0); }

  if( distRMax <= halfCarTolerance )
  {
    noSurfaces ++;
    G4ThreeVector_sum_assign( &sumnorm, nR );
  }
  if( This->fRMin && (distRMin <= halfCarTolerance) )
  {
    noSurfaces ++;
    G4ThreeVector_subtract_assign( &sumnorm, nR );
  }
#ifdef ENABLE_SLICED_TUBS
  if( This->fDPhi < twopi )   
  {
    if (distSPhi <= halfAngTolerance)  
    {
      noSurfaces ++;
      G4ThreeVector_sum_assign( &sumnorm, nPs );
    }
    if (distEPhi <= halfAngTolerance)  
    {
      noSurfaces ++;
      G4ThreeVector_sum_assign( &sumnorm, nPe );
    }
  }
#endif
  
  if (distZ <= halfCarTolerance)  
  {
    noSurfaces ++;
    if ( p.z >= 0.)  { G4ThreeVector_sum_assign( &sumnorm, nZ ); }
    else               { G4ThreeVector_subtract_assign( &sumnorm, nZ ); }
  }
  if ( noSurfaces == 0 )
  {
     norm = G4Tubs_ApproxSurfaceNormal(This,p);
  }
  else if ( noSurfaces == 1 )  { norm = sumnorm; }
  else                         { norm = G4ThreeVector_unit(sumnorm); }

  return norm;
}

////////////////////////////////////////////////////////////////////
//
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - 'if valid' implies tolerant checking of intersection points

SOLIDINLINE G4double G4Tubs_DistanceToIn_full(
				GEOMETRYLOC const G4Tubs *This,
				G4ThreeVector p,
				G4ThreeVector v)
{
  G4double snxt = kInfinity ;      // snxt = default return value
  G4double tolIRMax2 ;  // 'generous' radii squared
  G4double tolORMax2, tolIRMin2, tolODz, tolIDz ;
  //const G4double dRmax = 100.*This->fRMax;

  const G4double halfCarTolerance = 0.5*K_GEOMETRY_CAR_TOLERANCE;
  const G4double halfRadTolerance = 0.5*K_GEOMETRY_RAD_TOLERANCE;

  // Intersection point variables
  //
  G4double s, xi, yi, zi, rho2;
  G4double t1, t2, t3, b, c, d ;     // Quadratic solver variables 
  
#ifdef ENABLE_SLICED_TUBS
  G4double Dist, cosPsi, Comp, inum, iden, tolORMin2;
#endif
  
  // Calculate tolerant rmin and rmax

  if (This->fRMin > K_GEOMETRY_RAD_TOLERANCE)
  {
#ifdef ENABLE_SLICED_TUBS
    tolORMin2 = (This->fRMin - halfRadTolerance)*(This->fRMin - halfRadTolerance) ;
#endif
    tolIRMin2 = (This->fRMin + halfRadTolerance)*(This->fRMin + halfRadTolerance) ;
  }
  else
  {
#ifdef ENABLE_SLICED_TUBS
    tolORMin2 = 0.0 ;
#endif
    tolIRMin2 = 0.0 ;
  }
  tolORMax2 = (This->fRMax + halfRadTolerance)*(This->fRMax + halfRadTolerance) ;
  tolIRMax2 = (This->fRMax - halfRadTolerance)*(This->fRMax - halfRadTolerance) ;

  // Intersection with Z surfaces

  tolIDz = This->fDz - halfCarTolerance ;
  tolODz = This->fDz + halfCarTolerance ;

  if (fabs(p.z) >= tolIDz)
  {
    if ( p.z*v.z < 0 )    // at +Z going in -Z or visa versa
    {
      s = (fabs(p.z) - This->fDz)/fabs(v.z) ;   // Z intersect distance

      if(s < 0.0)  { s = 0.0; }

      xi   = p.x + s*v.x ;                // Intersection coords
      yi   = p.y + s*v.y ;
      rho2 = xi*xi + yi*yi ;

      // Check validity of intersection

      if ((tolIRMin2 <= rho2) && (rho2 <= tolIRMax2))
      {
#ifdef ENABLE_SLICED_TUBS
        if (!This->fPhiFullTube && rho2)
        {
          // Psi = angle made with central (average) phi of shape
          //
          inum   = xi*This->cosCPhi + yi*This->sinCPhi ;
          iden   = sqrt(rho2) ;
          cosPsi = inum/iden ;
          if (cosPsi >= This->cosHDPhiIT)  { return s ; }
        }
        else
        {
          return s ;
        }
#else
		return s;
#endif
      }
    }
    else
    {
      if ( snxt<halfCarTolerance )  { snxt=0; }
      return snxt ;  // On/outside extent, and heading away
                     // -> cannot intersect
    }
  }

  // -> Can not intersect z surfaces
  //
  // Intersection with rmax (possible return) and rmin (must also check phi)
  //
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //            t1                t2                t3

  t1 = 1.0 - v.z*v.z ;
  t2 = p.x*v.x + p.y*v.y ;
  t3 = p.x*p.x + p.y*p.y ;

  if ( t1 > 0 )        // Check not || to z axis
  {
    b = t2/t1 ;
    c = t3 - This->fRMax*This->fRMax ;
    if ((t3 >= tolORMax2) && (t2<0))   // This also handles the tangent case
    {
      // Try outer cylinder intersection
      //          c=(t3-fRMax*fRMax)/t1;

      c /= t1 ;
      d = b*b - c ;

      if (d >= 0)  // If real root
      {
        s = -b - sqrt(d) ;
        if (s >= 0)  // If 'forwards'
        {
		  //  myAssert( s <= dRmax ); // TODO: fix!
          /*if ( s>dRmax ) // Avoid rounding errors due to precision issues on
          {              // 64 bits systems. Split long distances and recompute
         
            G4double fTerm = s-fmod(s,dRmax);
           s = fTerm + DistanceToIn(p+fTerm*v,v);
          } */
          // Check z intersection
          //
          zi = p.z + s*v.z ;
          if (fabs(zi)<=tolODz)
          {
            // Z ok. Check phi intersection if reqd
            //
#ifdef ENABLE_SLICED_TUBS
            if (This->fPhiFullTube)
            {
              return s ;
            }
            else
            {
              xi     = p.x + s*v.x ;
              yi     = p.y + s*v.y ;
              cosPsi = (xi*This->cosCPhi + yi*This->sinCPhi)/This->fRMax ;
              if (cosPsi >= This->cosHDPhiIT)  { return s ; }
            }
#else
			return s;
#endif
          }  //  end if std::fabs(zi)
        }    //  end if (s>=0)
      }      //  end if (d>=0)
    }        //  end if (r>=fRMax)
    else 
    {
      // Inside outer radius :
      // check not inside, and heading through tubs (-> 0 to in)

      if ((t3 > tolIRMin2) && (t2 < 0) && (fabs(p.z) <= tolIDz))
      {
        // Inside both radii, delta r -ve, inside z extent

#ifdef ENABLE_SLICED_TUBS
        if (!This->fPhiFullTube)
        {
          inum   = p.x*This->cosCPhi + p.y*This->sinCPhi ;
          iden   = sqrt(t3) ;
          cosPsi = inum/iden ;
          if (cosPsi >= This->cosHDPhiIT)
          {
            // In the old version, the small negative tangent for the point
            // on surface was not taken in account, and returning 0.0 ...
            // New version: check the tangent for the point on surface and 
            // if no intersection, return kInfinity, if intersection instead
            // return s.
            //
            c = t3-This->fRMax*This->fRMax; 
            if ( c<=0.0 )
            {
              return 0.0;
            }
            else
            {
              c = c/t1 ;
              d = b*b-c;
              if ( d>=0.0 )
              {
                snxt = c/(-b+sqrt(d)); // using safe solution
                                            // for quadratic equation 
                if ( snxt < halfCarTolerance ) { snxt=0; }
                return snxt ;
              }      
              else
              {
                return kInfinity;
              }
            }
          } 
        }
        else
        {
#endif
          // In the old version, the small negative tangent for the point
          // on surface was not taken in account, and returning 0.0 ...
          // New version: check the tangent for the point on surface and 
          // if no intersection, return kInfinity, if intersection instead
          // return s.
          //
          c = t3 - This->fRMax*This->fRMax; 
          if ( c<=0.0 )
          {
            return 0.0;
          }
          else
          {
            c = c/t1 ;
            d = b*b-c;
            if ( d>=0.0 )
            {
              snxt= c/(-b+sqrt(d)); // using safe solution
                                         // for quadratic equation 
              if ( snxt < halfCarTolerance ) { snxt=0; }
              return snxt ;
            }      
            else
            {
              return kInfinity;
            }
          }
#ifdef ENABLE_SLICED_TUBS
        } // end if   (!fPhiFullTube)
#endif
      }   // end if   (t3>tolIRMin2)
    }     // end if   (Inside Outer Radius) 
    if ( This->fRMin )    // Try inner cylinder intersection
    {
      c = (t3 - This->fRMin*This->fRMin)/t1 ;
      d = b*b - c ;
      if ( d >= 0.0 )  // If real root
      {
        // Always want 2nd root - we are outside and know rmax Hit was bad
        // - If on surface of rmin also need farthest root

        s = -b + sqrt(d) ;
        if (s >= -halfCarTolerance)  // check forwards
        {
          // Check z intersection
          //
          if(s < 0.0)  { s = 0.0; }
          
		  //myAssert( s <= dRmax ); // TODO!
          
          /*if ( s>dRmax ) // Avoid rounding errors due to precision issues seen
          {              // 64 bits systems. Split long distances and recompute
          
            G4double fTerm = s-fmod(s,dRmax);
            s = This->fTerm + G4Tubs_DistanceToIn(This,p+This->fTerm*v,v);
          } */
          
          zi = p.z + s*v.z ;
          if (fabs(zi) <= tolODz)
          {
            // Z ok. Check phi
            //
#ifdef ENABLE_SLICED_TUBS
            if ( This->fPhiFullTube )
            {
              return s ; 
            }
            else
            {
              xi     = p.x + s*v.x ;
              yi     = p.y + s*v.y ;
              cosPsi = (xi*This->cosCPhi + yi*This->sinCPhi)/This->fRMin ;
              if (cosPsi >= This->cosHDPhiIT)
              {
                // Good inner radius isect
                // - but earlier phi isect still possible

                snxt = s ;
              }
            }
#else
			return s;
#endif
          }        //    end if std::fabs(zi)
        }          //    end if (s>=0)
      }            //    end if (d>=0)
    }              //    end if (fRMin)
  }

  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to K_GEOMETRY_CAR_TOLERANCE*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> use some form of loop Construct ?
  //
#ifdef ENABLE_SLICED_TUBS
  if ( !This->fPhiFullTube )
  {
    // First phi surface (Starting phi)
    //
    Comp    = v.x*This->sinSPhi - v.y*This->cosSPhi ;
                    
    if ( Comp < 0 )  // Component in outwards normal dirn
    {
      Dist = (p.y*This->cosSPhi - p.x*This->sinSPhi) ;

      if ( Dist < halfCarTolerance )
      {
        s = Dist/Comp ;

        if (s < snxt)
        {
          if ( s < 0 )  { s = 0.0; }
          zi = p.z + s*v.z ;
          if ( fabs(zi) <= tolODz )
          {
            xi   = p.x + s*v.x ;
            yi   = p.y + s*v.y ;
            rho2 = xi*xi + yi*yi ;

            if ( ( (rho2 >= tolIRMin2) && (rho2 <= tolIRMax2) )
              || ( (rho2 >  tolORMin2) && (rho2 <  tolIRMin2)
                && ( v.y*This->cosSPhi - v.x*This->sinSPhi >  0 )
                && ( v.x*This->cosSPhi + v.y*This->sinSPhi >= 0 )     )
              || ( (rho2 > tolIRMax2) && (rho2 < tolORMax2)
                && (v.y*This->cosSPhi - v.x*This->sinSPhi > 0)
                && (v.x*This->cosSPhi + v.y*This->sinSPhi < 0) )    )
            {
              // z and r intersections good
              // - check intersecting with correct half-plane
              //
              if ((yi*This->cosCPhi-xi*This->sinCPhi) <= halfCarTolerance) { snxt = s; }
            }
          }
        }
      }    
    }
      
    // Second phi surface (Ending phi)

    Comp    = -(v.x*This->sinEPhi - v.y*This->cosEPhi) ;
        
    if (Comp < 0 )  // Component in outwards normal dirn
    {
      Dist = -(p.y*This->cosEPhi - p.x*This->sinEPhi) ;

      if ( Dist < halfCarTolerance )
      {
        s = Dist/Comp ;

        if (s < snxt)
        {
          if ( s < 0 )  { s = 0; }
          zi = p.z + s*v.z ;
          if ( fabs(zi) <= tolODz )
          {
            xi   = p.x + s*v.x ;
            yi   = p.y + s*v.y ;
            rho2 = xi*xi + yi*yi ;
            if ( ( (rho2 >= tolIRMin2) && (rho2 <= tolIRMax2) )
                || ( (rho2 > tolORMin2)  && (rho2 < tolIRMin2)
                  && (v.x*This->sinEPhi - v.y*This->cosEPhi >  0)
                  && (v.x*This->cosEPhi + v.y*This->sinEPhi >= 0) )
                || ( (rho2 > tolIRMax2) && (rho2 < tolORMax2)
                  && (v.x*This->sinEPhi - v.y*This->cosEPhi > 0)
                  && (v.x*This->cosEPhi + v.y*This->sinEPhi < 0) ) )
            {
              // z and r intersections good
              // - check intersecting with correct half-plane
              //
              if ( (yi*This->cosCPhi-xi*This->sinCPhi) >= 0 ) { snxt = s; }
            }                         //?? >=-halfCarTolerance
          }
        }
      }
    }         //  Comp < 0
  }           //  !fPhiFullTube 
#endif
  if ( snxt<halfCarTolerance )  { snxt=0; }
  return snxt ;
}
 
//////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points
//   Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

SOLIDINLINE G4double G4Tubs_DistanceToIn(GEOMETRYLOC const G4Tubs *This, G4ThreeVector p)
{
  G4double safe=0.0, rho, safe1, safe2, safe3 ;
  
#ifdef ENABLE_SLICED_TUBS
  G4double safePhi, cosPsi ;
#endif

  rho   = sqrt(p.x*p.x + p.y*p.y) ;
  safe1 = This->fRMin - rho ;
  safe2 = rho - This->fRMax ;
  safe3 = fabs(p.z) - This->fDz ;

  if ( safe1 > safe2 ) { safe = safe1; }
  else                 { safe = safe2; }
  if ( safe3 > safe )  { safe = safe3; }

#ifdef ENABLE_SLICED_TUBS
  if ( (!This->fPhiFullTube) && (rho) )
  {
    // Psi=angle from central phi to point
    //
    cosPsi = (p.x*This->cosCPhi + p.y*This->sinCPhi)/rho ;
    
    if ( cosPsi < cos(This->fDPhi*0.5) )
    {
      // Point lies outside phi range

      if ( (p.y*This->cosCPhi - p.x*This->sinCPhi) <= 0 )
      {
        safePhi = fabs(p.x*This->sinSPhi - p.y*This->cosSPhi) ;
      }
      else
      {
        safePhi = fabs(p.x*This->sinEPhi - p.y*This->cosEPhi) ;
      }
      if ( safePhi > safe )  { safe = safePhi; }
    }
  }
#endif

  if ( safe < 0 )  { safe = 0; }
  return safe ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

SOLIDINLINE G4double G4Tubs_DistanceToOut_full(
			   GEOMETRYLOC const G4Tubs *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n)
{  
  typedef enum {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ} ESide;
  ESide side=kNull , sider=kNull;
  G4double snxt, sr=kInfinity, pdist ;
  G4double deltaR, t1, t2, t3, b, c, d2, roMin2 ;

  const G4double halfCarTolerance = K_GEOMETRY_CAR_TOLERANCE*0.5;
 
  // Vars for phi intersection:

  G4double xi, yi, roi2 ;
 
 #ifdef ENABLE_SLICED_TUBS
  const G4double pi = M_PI;
  const G4double halfAngTolerance = K_GEOMETRY_ANG_TOLERANCE*0.5;
  
  G4double sphi=kInfinity, sphi2, pDistS, compS, pDistE, compE, vphi;
  ESide sidephi=kNull;
#endif
 
  // Z plane intersection

  if (v.z > 0 )
  {
    pdist = This->fDz - p.z ;
    if ( pdist > halfCarTolerance )
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
      return snxt = 0 ;
    }
  }
  else if ( v.z < 0 )
  {
    pdist = This->fDz + p.z ;

    if ( pdist > halfCarTolerance )
    {
      snxt = -pdist/v.z ;
      side = kMZ ;
    }
    else
    {
      if (calcNorm)
      {
        *n         = G4ThreeVector_create(0,0,-1) ;
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
  }
  else
  {
    snxt = kInfinity ;    // Travel perpendicular to z axis
    side = kNull;
  }

  // Radial Intersections
  //
  // Find intersection with cylinders at rmax/rmin
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //
  //            t1                t2                    t3

  t1   = 1.0 - v.z*v.z ;      // since v normalised
  t2   = p.x*v.x + p.y*v.y ;
  t3   = p.x*p.x + p.y*p.y ;

  if ( snxt > 10*(This->fDz+This->fRMax) )  { roi2 = 2*This->fRMax*This->fRMax; }
  else  { roi2 = snxt*snxt*t1 + 2*snxt*t2 + t3; }        // radius^2 on +-fDz

  if ( t1 > 0 ) // Check not parallel
  {
    // Calculate sr, r exit distance
     
    if ( (t2 >= 0.0) && (roi2 > This->fRMax*(This->fRMax + K_GEOMETRY_RAD_TOLERANCE)) )
    {
      // Delta r not negative => leaving via rmax

      deltaR = t3 - This->fRMax*This->fRMax ;

      // NOTE: Should use rho-fRMax<-K_GEOMETRY_RAD_TOLERANCE*0.5
      // - avoid sqrt for efficiency

      if ( deltaR < -K_GEOMETRY_RAD_TOLERANCE*This->fRMax )
      {
        b     = t2/t1 ;
        c     = deltaR/t1 ;
        d2    = b*b-c;
        if( d2 >= 0 ) { sr = -b + sqrt(d2); }
        else          { sr = 0.; }
        sider = kRMax ;
      }
      else
      {
        // On tolerant boundary & heading outwards (or perpendicular to)
        // outer radial surface -> leaving immediately

        if ( calcNorm ) 
        {
          *n         = G4ThreeVector_create(p.x/This->fRMax,p.y/This->fRMax,0) ;
          *validNorm = true ;
        }
        return snxt = 0 ; // Leaving by rmax immediately
      }
    }             
    else if ( t2 < 0. ) // i.e.  t2 < 0; Possible rmin intersection
    {
      roMin2 = t3 - t2*t2/t1 ; // min ro2 of the plane of movement 

      if ( This->fRMin && (roMin2 < This->fRMin*(This->fRMin - K_GEOMETRY_RAD_TOLERANCE)) )
      {
        deltaR = t3 - This->fRMin*This->fRMin ;
        b      = t2/t1 ;
        c      = deltaR/t1 ;
        d2     = b*b - c ;

        if ( d2 >= 0 )   // Leaving via rmin
        {
          // NOTE: SHould use rho-rmin>K_GEOMETRY_RAD_TOLERANCE*0.5
          // - avoid sqrt for efficiency

          if (deltaR > K_GEOMETRY_RAD_TOLERANCE*This->fRMin)
          {
            sr    = -b-sqrt(d2) ;
            sider = kRMin ;
          }
          else
          {
            if ( calcNorm ) { *validNorm = false; }  // Concave side
            return snxt = 0.0;
          }
        }
        else    // No rmin intersect -> must be rmax intersect
        {
          deltaR = t3 - This->fRMax*This->fRMax ;
          c     = deltaR/t1 ;
          d2    = b*b-c;
          if( d2 >=0. )
          {
            sr     = -b + sqrt(d2) ;
            sider  = kRMax ;
          }
          else // Case: On the border+t2<K_GEOMETRY_RAD_TOLERANCE
               //       (v is perpendicular to the surface)
          {
            if (calcNorm)
            {
              *n = G4ThreeVector_create(p.x/This->fRMax,p.y/This->fRMax,0) ;
              *validNorm = true ;
            }
            return snxt = 0.0;
          }
        }
      }
      else if ( roi2 > This->fRMax*(This->fRMax + K_GEOMETRY_RAD_TOLERANCE) )
           // No rmin intersect -> must be rmax intersect
      {
        deltaR = t3 - This->fRMax*This->fRMax ;
        b      = t2/t1 ;
        c      = deltaR/t1;
        d2     = b*b-c;
        if( d2 >= 0 )
        {
          sr     = -b + sqrt(d2) ;
          sider  = kRMax ;
        }
        else // Case: On the border+t2<K_GEOMETRY_RAD_TOLERANCE
             //       (v is perpendicular to the surface)
        {
          if (calcNorm)
          {
            *n = G4ThreeVector_create(p.x/This->fRMax,p.y/This->fRMax,0) ;
            *validNorm = true ;
          }
          return snxt = 0.0;
        }
      }
    }
    
    // Phi Intersection

#ifdef ENABLE_SLICED_TUBS
    if ( !This->fPhiFullTube )
    {
      // add angle calculation with correction 
      // of the difference in domain of atan2 and Sphi
      //
      vphi = atan2(v.y,v.x) ;
     
      if ( vphi < This->fSPhi - halfAngTolerance  )             { vphi += twopi; }
      else if ( vphi > This->fSPhi + This->fDPhi + halfAngTolerance ) { vphi -= twopi; }


      if ( p.x || p.y )  // Check if on z axis (rho not needed later)
      {
        // pDist -ve when inside

        pDistS = p.x*This->sinSPhi - p.y*This->cosSPhi ;
        pDistE = -p.x*This->sinEPhi + p.y*This->cosEPhi ;

        // Comp -ve when in direction of outwards normal

        compS   = -This->sinSPhi*v.x + This->cosSPhi*v.y ;
        compE   =  This->sinEPhi*v.x - This->cosEPhi*v.y ;
       
        sidephi = kNull;
        
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
              if( (fabs(xi)<=K_GEOMETRY_CAR_TOLERANCE)&&(fabs(yi)<=K_GEOMETRY_CAR_TOLERANCE) )
              {
                sidephi = kSPhi;
                if (((This->fSPhi-halfAngTolerance)<=vphi)
                   &&((This->fSPhi+This->fDPhi+halfAngTolerance)>=vphi))
                {
                  sphi = kInfinity;
                }
              }
              else if ( yi*This->cosCPhi-xi*This->sinCPhi >=0 )
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
              
              if ((fabs(xi)<=K_GEOMETRY_CAR_TOLERANCE)&&(fabs(yi)<=K_GEOMETRY_CAR_TOLERANCE))
              {
                // Leaving via ending phi
                //
                if( !((This->fSPhi-halfAngTolerance <= vphi)
                     &&(This->fSPhi+This->fDPhi+halfAngTolerance >= vphi)) )
                {
                  sidephi = kEPhi ;
                  if ( pDistE <= -halfCarTolerance )  { sphi = sphi2 ; }
                  else                                { sphi = 0.0 ;   }
                }
              } 
              else    // Check intersecting with correct half-plane 

              if ( (yi*This->cosCPhi-xi*This->sinCPhi) >= 0)
              {
                // Leaving via ending phi
                //
                sidephi = kEPhi ;
                if ( pDistE <= -halfCarTolerance ) { sphi = sphi2 ; }
                else                               { sphi = 0.0 ;   }
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
               
        if ( (This->fSPhi - halfAngTolerance <= vphi)
           && (vphi <= This->fSPhi + This->fDPhi + halfAngTolerance ) )
        {
          sphi = kInfinity ;
        }
        else
        {
          sidephi = kSPhi ; // arbitrary 
          sphi    = 0.0 ;
        }
      }
      if (sphi < snxt)  // Order intersecttions
      {
        snxt = sphi ;
        side = sidephi ;
      }
    }
#endif
    if (sr < snxt)  // Order intersections
    {
      snxt = sr ;
      side = sider ;
    }
  }
  if (calcNorm)
  {
    switch(side)
    {
      case kRMax:
        // Note: returned vector not normalised
        // (divide by fRMax for unit vector)
        //
        xi = p.x + snxt*v.x ;
        yi = p.y + snxt*v.y ;
        *n = G4ThreeVector_create(xi/This->fRMax,yi/This->fRMax,0) ;
        *validNorm = true ;
        break ;

      case kRMin:
        *validNorm = false ;  // Rmin is inconvex
        break ;

#ifdef ENABLE_SLICED_TUBS
      case kSPhi:
        if ( This->fDPhi <= pi )
        {
          *n         = G4ThreeVector_create(This->sinSPhi,-This->cosSPhi,0) ;
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;

      case kEPhi:
        if (This->fDPhi <= pi)
        {
          *n = G4ThreeVector_create(-This->sinEPhi,This->cosEPhi,0) ;
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
  if ( snxt<halfCarTolerance )  { snxt=0 ; }

  return snxt ;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

SOLIDINLINE G4double G4Tubs_DistanceToOut(GEOMETRYLOC const G4Tubs *This, G4ThreeVector p)
{
  G4double safe=0.0, rho, safeR1, safeR2, safeZ;
  rho = sqrt(p.x*p.x + p.y*p.y) ;
  
#ifdef ENABLE_SLICED_TUBS
  G4double safePhi;
#endif

  if ( This->fRMin )
  {
    safeR1 = rho   - This->fRMin ;
    safeR2 = This->fRMax - rho ;
 
    if ( safeR1 < safeR2 ) { safe = safeR1 ; }
    else                   { safe = safeR2 ; }
  }
  else
  {
    safe = This->fRMax - rho ;
  }
  safeZ = This->fDz - fabs(p.z) ;

  if ( safeZ < safe )  { safe = safeZ ; }

  // Check if phi divided, Calc distances closest phi plane
  //
#ifdef ENABLE_SLICED_TUBS
  if ( !This->fPhiFullTube )
  {
    if ( p.y*This->cosCPhi-p.x*This->sinCPhi <= 0 )
    {
      safePhi = -(p.x*This->sinSPhi - p.y*This->cosSPhi) ;
    }
    else
    {
      safePhi = (p.x*This->sinEPhi - p.y*This->cosEPhi) ;
    }
    if (safePhi < safe)  { safe = safePhi ; }
  }
#endif
  if ( safe < 0 )  { safe = 0 ; }

  return safe ;  
}


#endif // not host code

///////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

SOLIDINLINE EInside G4Tubs_Inside(GEOMETRYLOC const G4Tubs *This, G4ThreeVector p)
{
  G4double r2,tolRMin,tolRMax;
  EInside in = kOutside ;
  const G4double halfCarTolerance=K_GEOMETRY_CAR_TOLERANCE*0.5;
  const G4double halfRadTolerance=K_GEOMETRY_RAD_TOLERANCE*0.5;
  
#ifdef ENABLE_SLICED_TUBS
  const G4double halfAngTolerance=K_GEOMETRY_ANG_TOLERANCE*0.5;
  G4double pPhi;
#endif

  if (fabs(p.z) <= This->fDz - halfCarTolerance)
  {
    r2 = p.x*p.x + p.y*p.y ;

    if (This->fRMin) { tolRMin = This->fRMin + halfRadTolerance ; }
    else       { tolRMin = 0 ; }

    tolRMax = This->fRMax - halfRadTolerance ;
      
    if ((r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax))
    {
#ifndef ENABLE_SLICED_TUBS
		in = kInside;
#else
      if ( This->fPhiFullTube )
      {
        in = kInside ;
      }
      else
      {
        // Try inner tolerant phi boundaries (=>inside)
        // if not inside, try outer tolerant phi boundaries

        if ((tolRMin==0)&&(p.x<=halfCarTolerance)&&(p.y<=halfCarTolerance))
        {
          in=kSurface;
        }
        else
        {
          pPhi = atan2(p.y,p.x) ;
          if ( pPhi < -halfAngTolerance )  { pPhi += twopi; } // 0<=pPhi<2pi

          if ( This->fSPhi >= 0 )
          {
            if ( (fabs(pPhi) < halfAngTolerance)
              && (fabs(This->fSPhi + This->fDPhi - twopi) < halfAngTolerance) )
            { 
              pPhi += twopi ; // 0 <= pPhi < 2pi
            }
            if ( (pPhi >= This->fSPhi + halfAngTolerance)
              && (pPhi <= This->fSPhi + This->fDPhi - halfAngTolerance) )
            {
              in = kInside ;
            }
            else if ( (pPhi >= This->fSPhi - halfAngTolerance)
                   && (pPhi <= This->fSPhi + This->fDPhi + halfAngTolerance) )
            {
              in = kSurface ;
            }
          }
          else  // fSPhi < 0
          {
            if ( (pPhi <= This->fSPhi + twopi - halfAngTolerance)
              && (pPhi >= This->fSPhi + This->fDPhi  + halfAngTolerance) ) {;} //kOutside
            else if ( (pPhi <= This->fSPhi + twopi + halfAngTolerance)
                   && (pPhi >= This->fSPhi + This->fDPhi  - halfAngTolerance) )
            {
              in = kSurface ;
            }
            else
            {
              in = kInside ;
            }
          }
        }                    
      }
#endif
    }
    else  // Try generous boundaries
    {
      tolRMin = This->fRMin - halfRadTolerance ;
      tolRMax = This->fRMax + halfRadTolerance ;

      if ( tolRMin < 0 )  { tolRMin = 0; }

      if ( (r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax) )
      {
#ifndef ENABLE_SLICED_TUBS
		in = kSurface;
#else
        if (This->fPhiFullTube || (r2 <=halfRadTolerance*halfRadTolerance) )
        {                        // Continuous in phi or on z-axis
          in = kSurface ;
        }
        else // Try outer tolerant phi boundaries only
        {
          pPhi = atan2(p.y,p.x) ;

          if ( pPhi < -halfAngTolerance)  { pPhi += twopi; } // 0<=pPhi<2pi
          if ( This->fSPhi >= 0 )
          {
            if ( (fabs(pPhi) < halfAngTolerance)
              && (fabs(This->fSPhi + This->fDPhi - twopi) < halfAngTolerance) )
            { 
              pPhi += twopi ; // 0 <= pPhi < 2pi
            }
            if ( (pPhi >= This->fSPhi - halfAngTolerance)
              && (pPhi <= This->fSPhi + This->fDPhi + halfAngTolerance) )
            {
              in = kSurface ;
            }
          }
          else  // fSPhi < 0
          {
            if ( (pPhi <= This->fSPhi + twopi - halfAngTolerance)
              && (pPhi >= This->fSPhi + This->fDPhi + halfAngTolerance) ) {;} // kOutside
            else
            {
              in = kSurface ;
            }
          }
        }
#endif
      }
    }
  }
  else if (fabs(p.z) <= This->fDz + halfCarTolerance)
  {                                          // Check within tolerant r limits
    r2      = p.x*p.x + p.y*p.y ;
    tolRMin = This->fRMin - halfRadTolerance ;
    tolRMax = This->fRMax + halfRadTolerance ;

    if ( tolRMin < 0 )  { tolRMin = 0; }

    if ( (r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax) )
    {
#ifndef ENABLE_SLICED_TUBS
		in = kSurface;
#else
      if (This->fPhiFullTube || (r2 <=halfRadTolerance*halfRadTolerance))
      {                        // Continuous in phi or on z-axis
        in = kSurface ;
      }
      else // Try outer tolerant phi boundaries
      {
        pPhi = atan2(p.y,p.x) ;

        if ( pPhi < -halfAngTolerance )  { pPhi += twopi; }  // 0<=pPhi<2pi
        if ( This->fSPhi >= 0 )
        {
          if ( (fabs(pPhi) < halfAngTolerance)
            && (fabs(This->fSPhi + This->fDPhi - twopi) < halfAngTolerance) )
          { 
            pPhi += twopi ; // 0 <= pPhi < 2pi
          }
          if ( (pPhi >= This->fSPhi - halfAngTolerance)
            && (pPhi <= This->fSPhi + This->fDPhi + halfAngTolerance) )
          {
            in = kSurface;
          }
        }
        else  // fSPhi < 0
        {
          if ( (pPhi <= This->fSPhi + twopi - halfAngTolerance)
            && (pPhi >= This->fSPhi + This->fDPhi  + halfAngTolerance) ) {;}
          else
          {
            in = kSurface ;
          }
        }      
      }
#endif
    }
  }
  return in;
}

#if defined( HOST_CODE ) && defined( ENABLE_VOXEL_NAVIGATION )

#include "G4AffineTransform.h"
#include "G4VoxelLimits.c"
#include "G4VSolid.c"

/////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -fDz cross section
//          [4-7] +fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

G4ThreeVectorList*
G4Tubs_CreateRotatedVertices( const G4Tubs *This, G4AffineTransform *pTransform )
{
  G4ThreeVectorList* vertices ;
  G4ThreeVector vertex0, vertex1, vertex2, vertex3 ;
  G4double meshAngle, meshRMax, crossAngle,
           cosCrossAngle, sinCrossAngle, sAngle;
  G4double rMaxX, rMaxY, rMinX, rMinY, meshRMin ;
  G4int crossSection, noCrossSections;

#ifdef ENABLE_SLICED_TUBS
  // Compute no of cross-sections necessary to mesh tube
  //
  noCrossSections = G4int(This->fDPhi/K_GEOMETRY_MESH_ANGLE_DEFAULT) + 1 ;
#else
  noCrossSections = G4int(twopi/K_GEOMETRY_MESH_ANGLE_DEFAULT) + 1 ;
#endif

  if ( noCrossSections < K_GEOMETRY_MIN_MESH_SECTIONS )
  {
    noCrossSections = K_GEOMETRY_MIN_MESH_SECTIONS ;
  }
  else if (noCrossSections > K_GEOMETRY_MAX_MESH_SECTIONS)
  {
    noCrossSections = K_GEOMETRY_MAX_MESH_SECTIONS ;
  }
  // noCrossSections = 4 ;

#ifdef ENABLE_SLICED_TUBS
  meshAngle = This->fDPhi/(noCrossSections - 1) ;
  // meshAngle = fDPhi/(noCrossSections) ;
#else
  meshAngle = twopi/(noCrossSections - 1) ;
#endif

  meshRMax  = (This->fRMax+100*K_GEOMETRY_CAR_TOLERANCE)/std::cos(meshAngle*0.5) ;
  meshRMin = This->fRMin - 100*K_GEOMETRY_CAR_TOLERANCE ; 
 
  // If complete in phi, set start angle such that mesh will be at fRMax
  // on the x axis. Will give better extent calculations when not rotated.

#ifdef ENABLE_SLICED_TUBS
  if ((This->fPhiFullTube) && (This->fSPhi == 0) )  { sAngle = -meshAngle*0.5 ; }
  else                                { sAngle =  This->fSPhi ; }
#else
	sAngle = -meshAngle*0.5;
#endif
    
  vertices = new G4ThreeVectorList;
  vertices->reserve(noCrossSections*4);
    
  if ( vertices )
  {
    for (crossSection = 0 ; crossSection < noCrossSections ; crossSection++ )
    {
      // Compute coordinates of cross section at section crossSection

      crossAngle    = sAngle + crossSection*meshAngle ;
      cosCrossAngle = std::cos(crossAngle) ;
      sinCrossAngle = std::sin(crossAngle) ;

      rMaxX = meshRMax*cosCrossAngle ;
      rMaxY = meshRMax*sinCrossAngle ;

      if(meshRMin <= 0.0)
      {
        rMinX = 0.0 ;
        rMinY = 0.0 ;
      }
      else
      {
        rMinX = meshRMin*cosCrossAngle ;
        rMinY = meshRMin*sinCrossAngle ;
      }
      vertex0 = G4ThreeVector_create(rMinX,rMinY,-This->fDz) ;
      vertex1 = G4ThreeVector_create(rMaxX,rMaxY,-This->fDz) ;
      vertex2 = G4ThreeVector_create(rMaxX,rMaxY,+This->fDz) ;
      vertex3 = G4ThreeVector_create(rMinX,rMinY,+This->fDz) ;

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


////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
SOLIDINLINE
G4bool G4Tubs_CalculateExtent(
			   const G4Tubs *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax )
{
  if ( (!G4AffineTransform_IsRotated(&pTransform)) &&
#ifdef ENABLE_SLICED_TUBS
	(This->fDPhi == twopi) &&
#endif
	(This->fRMin == 0) )
  {
    // Special case handling for unrotated solid tubes
    // Compute x/y/z mins and maxs fro bounding box respecting limits,
    // with early returns if outside limits. Then switch() on pAxis,
    // and compute exact x and y limit for x/y case
      
    G4double xoffset, xMin, xMax;
    G4double yoffset, yMin, yMax;
    G4double zoffset, zMin, zMax;

    G4double diff1, diff2, maxDiff, newMin, newMax;
    G4double xoff1, xoff2, yoff1, yoff2, delta;

    xoffset = G4AffineTransform_NetTranslation(&pTransform).x;
    xMin = xoffset - This->fRMax;
    xMax = xoffset + This->fRMax;

    if (G4VoxelLimits_IsXLimited(&pVoxelLimit))
    {
      if ( (xMin > G4VoxelLimits_GetMaxXExtent(&pVoxelLimit))
        || (xMax < G4VoxelLimits_GetMinXExtent(&pVoxelLimit)) )
      {
        return false;
      }
      else
      {
        if (xMin < G4VoxelLimits_GetMinXExtent(&pVoxelLimit))
        {
          xMin = G4VoxelLimits_GetMinXExtent(&pVoxelLimit);
        }
        if (xMax > G4VoxelLimits_GetMaxXExtent(&pVoxelLimit))
        {
          xMax = G4VoxelLimits_GetMaxXExtent(&pVoxelLimit);
        }
      }
    }
    yoffset = G4AffineTransform_NetTranslation(&pTransform).y;
    yMin    = yoffset - This->fRMax;
    yMax    = yoffset + This->fRMax;

    if ( G4VoxelLimits_IsYLimited(&pVoxelLimit) )
    {
      if ( (yMin > G4VoxelLimits_GetMaxYExtent(&pVoxelLimit))
        || (yMax < G4VoxelLimits_GetMinYExtent(&pVoxelLimit)) )
      {
        return false;
      }
      else
      {
        if (yMin < G4VoxelLimits_GetMinYExtent(&pVoxelLimit))
        {
          yMin = G4VoxelLimits_GetMinYExtent(&pVoxelLimit);
        }
        if (yMax > G4VoxelLimits_GetMaxYExtent(&pVoxelLimit))
        {
          yMax=G4VoxelLimits_GetMaxYExtent(&pVoxelLimit);
        }
      }
    }
    zoffset = G4AffineTransform_NetTranslation(&pTransform).z;
    zMin    = zoffset - This->fDz;
    zMax    = zoffset + This->fDz;

    if ( G4VoxelLimits_IsZLimited(&pVoxelLimit) )
    {
      if ( (zMin > G4VoxelLimits_GetMaxZExtent(&pVoxelLimit))
        || (zMax < G4VoxelLimits_GetMinZExtent(&pVoxelLimit)) )
      {
        return false;
      }
      else
      {
        if (zMin < G4VoxelLimits_GetMinZExtent(&pVoxelLimit))
        {
          zMin = G4VoxelLimits_GetMinZExtent(&pVoxelLimit);
        }
        if (zMax > G4VoxelLimits_GetMaxZExtent(&pVoxelLimit))
        {
          zMax = G4VoxelLimits_GetMaxZExtent(&pVoxelLimit);
        }
      }
    }
    switch ( pAxis )  // Known to cut cylinder
    {
      case kXAxis :
      {
        yoff1 = yoffset - yMin;
        yoff2 = yMax    - yoffset;

        if ( (yoff1 >= 0) && (yoff2 >= 0) ) // Y limits cross max/min x
        {                                   // => no change
          *pMin = xMin;
          *pMax = xMax;
        }
        else
        {
          // Y limits don't cross max/min x => compute max delta x,
          // hence new mins/maxs

          delta   = This->fRMax*This->fRMax - yoff1*yoff1;
          diff1   = (delta>0.) ? std::sqrt(delta) : 0.;
          delta   = This->fRMax*This->fRMax - yoff2*yoff2;
          diff2   = (delta>0.) ? std::sqrt(delta) : 0.;
          maxDiff = (diff1 > diff2) ? diff1:diff2;
          newMin  = xoffset - maxDiff;
          newMax  = xoffset + maxDiff;
          *pMin    = (newMin < xMin) ? xMin : newMin;
          *pMax    = (newMax > xMax) ? xMax : newMax;
        }
        break;
      }
      case kYAxis :
      {
        xoff1 = xoffset - xMin;
        xoff2 = xMax - xoffset;

        if ( (xoff1 >= 0) && (xoff2 >= 0) ) // X limits cross max/min y
        {                                   // => no change
          *pMin = yMin;
          *pMax = yMax;
        }
        else
        {
          // X limits don't cross max/min y => compute max delta y,
          // hence new mins/maxs

          delta   = This->fRMax*This->fRMax - xoff1*xoff1;
          diff1   = (delta>0.) ? std::sqrt(delta) : 0.;
          delta   = This->fRMax*This->fRMax - xoff2*xoff2;
          diff2   = (delta>0.) ? std::sqrt(delta) : 0.;
          maxDiff = (diff1 > diff2) ? diff1 : diff2;
          newMin  = yoffset - maxDiff;
          newMax  = yoffset + maxDiff;
          *pMin    = (newMin < yMin) ? yMin : newMin;
          *pMax    = (newMax > yMax) ? yMax : newMax;
        }
        break;
      }
      case kZAxis:
      {
        *pMin = zMin;
        *pMax = zMax;
        break;
      }
      default:
        break;
    }
    *pMin -= K_GEOMETRY_CAR_TOLERANCE;
    *pMax += K_GEOMETRY_CAR_TOLERANCE;
    return true;
  }

  else // Calculate rotated vertex coordinates
  {
    G4int i, noEntries, noBetweenSections4;
    G4bool existsAfterClip = false;
    G4ThreeVectorList* vertices = G4Tubs_CreateRotatedVertices(This,&pTransform);

    *pMin =  kInfinity;
    *pMax = -kInfinity;

    noEntries = vertices->size();
    noBetweenSections4 = noEntries - 4;
    
    for ( i = 0 ; i < noEntries ; i += 4 )
    {
      G4VSolid_ClipCrossSection(vertices, i, &pVoxelLimit, pAxis, *pMin, *pMax);
    }
    for ( i = 0 ; i < noBetweenSections4 ; i += 4 )
    {
      G4VSolid_ClipBetweenSections(vertices, i, &pVoxelLimit, pAxis, *pMin, *pMax);
    }
    if ( (*pMin != kInfinity) || (*pMax != -kInfinity) )
    {
      existsAfterClip = true;
      *pMin -= K_GEOMETRY_CAR_TOLERANCE; // Add 2*tolerance to avoid precision troubles
      *pMax += K_GEOMETRY_CAR_TOLERANCE;
    }
    else
    {
      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.

      G4ThreeVector clipCentre = G4ThreeVector_create(
             (G4VoxelLimits_GetMinXExtent(&pVoxelLimit)+G4VoxelLimits_GetMaxXExtent(&pVoxelLimit))*0.5,
             (G4VoxelLimits_GetMinYExtent(&pVoxelLimit)+G4VoxelLimits_GetMaxYExtent(&pVoxelLimit))*0.5,
             (G4VoxelLimits_GetMinZExtent(&pVoxelLimit)+G4VoxelLimits_GetMaxZExtent(&pVoxelLimit))*0.5 );
        
      const G4AffineTransform invt = G4AffineTransform_Inverse(&pTransform);
        
      if ( G4Tubs_Inside(This,G4AffineTransform_TransformPoint(&invt,clipCentre)) != kOutside )
      {
        existsAfterClip = true;
        *pMin            = G4VoxelLimits_GetMinExtent(&pVoxelLimit, pAxis);
        *pMax            = G4VoxelLimits_GetMaxExtent(&pVoxelLimit, pAxis);
      }
    }
    delete vertices;
    return existsAfterClip;
  }
}

#endif // voxel nav & host code
