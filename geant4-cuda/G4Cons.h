
/** G4Cons header, based on G4Cons.hh of Geant 4.9.3 */

#ifndef G4Cons_HH
#define G4Cons_HH

#include "G4VSolid.h"

#define ENABLE_SLICED_CONS
//#define ENABLE_SLICED_POLYCONS
//#include <stdio.h>

//#define DISABLE_CACHE_MEMBERS

typedef struct
{
	G4VSolid solid;
	
    G4double fRmin1, fRmin2, fRmax1, fRmax2, fDz;
    
#ifdef ENABLE_SLICED_CONS
    G4double fSPhi, fDPhi;
      //
      // Radial and angular dimensions

#ifndef DISABLE_CACHE_MEMBERS
    G4double sinCPhi, cosCPhi, cosHDPhiOT, cosHDPhiIT,
             sinSPhi, cosSPhi, sinEPhi, cosEPhi;
      //
      // Cached trigonometric values
#endif

    G4bool fPhiFullCone;
      //
      // Flag for identification of section or full cone
      
     #ifdef DOUBLE_PRECISION
     int align1;
     #endif
#endif
}
G4Cons;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HOST_CODE
SOLIDINLINE void G4Cons_ctor(G4Cons *This, 
                      G4double  pRmin1, G4double pRmax1,
                      G4double  pRmin2, G4double pRmax2,
                      G4double pDz,
                      G4double pSPhi, G4double pDPhi);

#ifdef ENABLE_VOXEL_NAVIGATION
SOLIDINLINE
G4bool G4Cons_CalculateExtent(
			   const G4Cons *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax);
#endif
			   
#else
SOLIDINLINE EInside G4Cons_Inside(GEOMETRYLOC const G4Cons *This, G4ThreeVector p);

SOLIDINLINE G4ThreeVector G4Cons_SurfaceNormal(GEOMETRYLOC const G4Cons *This, G4ThreeVector p);


SOLIDINLINE G4double G4Cons_DistanceToIn_full(
				GEOMETRYLOC const G4Cons *This,
				G4ThreeVector p,
				G4ThreeVector v);

SOLIDINLINE G4double G4Cons_DistanceToIn(GEOMETRYLOC const G4Cons *This, G4ThreeVector p);

SOLIDINLINE G4double G4Cons_DistanceToOut_full(
			   GEOMETRYLOC const G4Cons *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n);

SOLIDINLINE G4double G4Cons_DistanceToOut(GEOMETRYLOC const G4Cons *This, G4ThreeVector p);
#endif

#if defined(INLINE_EVERYTHING) || defined(HOST_CODE)
#include "G4Cons.c"
#endif

#ifdef __cplusplus
}
#endif

#endif
