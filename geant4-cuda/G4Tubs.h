
/**
 * G4Tubs header
 * based on G4Tubs.hh of Geant 4.9.3
 */

#ifndef G4TUBS_HH
#define G4TUBS_HH

#include "G4VSolid.h"

//#define ENABLE_SLICED_TUBS

typedef struct
{
	G4VSolid solid;
	
    G4double fRMin, fRMax, fDz;
      //
      // Radial and angular dimensions
   
#ifdef ENABLE_SLICED_TUBS
	G4double fSPhi, fDPhi;

	// TODO: it may not be smart to cache these on the GPU
	
    G4double sinCPhi, cosCPhi, cosHDPhiOT, cosHDPhiIT,
             sinSPhi, cosSPhi, sinEPhi, cosEPhi;
      //
      // Cached trigonometric values

    G4bool fPhiFullTube;
      //
      // Flag for identification of section or full tube
      
     #ifdef DOUBLE_PRECISION
     int align1;
     #endif
#endif
}
G4Tubs;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HOST_CODE
SOLIDINLINE void G4Tubs_ctor(G4Tubs *This, G4double pRMin, G4double pRMax,
                  G4double pDz, G4double pSPhi, G4double pDPhi );

#ifdef ENABLE_VOXEL_NAVIGATION
SOLIDINLINE
G4bool G4Tubs_CalculateExtent(
			   const G4Tubs *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax);
#endif

#else

SOLIDINLINE EInside G4Tubs_Inside(GEOMETRYLOC const G4Tubs *This, G4ThreeVector p);

SOLIDINLINE G4ThreeVector G4Tubs_SurfaceNormal(GEOMETRYLOC const G4Tubs *This, G4ThreeVector p);


SOLIDINLINE G4double G4Tubs_DistanceToIn_full(
				GEOMETRYLOC const G4Tubs *This,
				G4ThreeVector p,
				G4ThreeVector v);

SOLIDINLINE G4double G4Tubs_DistanceToIn(GEOMETRYLOC const G4Tubs *This, G4ThreeVector p);

SOLIDINLINE G4double G4Tubs_DistanceToOut_full(
			   GEOMETRYLOC const G4Tubs *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n);

SOLIDINLINE G4double G4Tubs_DistanceToOut(GEOMETRYLOC const G4Tubs *This, G4ThreeVector p);
#endif

#if defined(INLINE_EVERYTHING) || defined(HOST_CODE)
#include "G4Tubs.c"
#endif

#ifdef __cplusplus
}
#endif

#endif
