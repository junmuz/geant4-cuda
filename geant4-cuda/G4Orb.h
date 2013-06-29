
/**
 * G4Orb header
 * based on G4Orb.hh of Geant 4.9.3
 */

#ifndef G4Orb_HH
#define G4Orb_HH

#include "G4VSolid.h"

typedef struct
{
	G4VSolid solid;
	
    G4double fRmax;
    G4double fRmaxTolerance;
#ifndef DOUBLE_PRECISION
	G4double align;
#endif
}
G4Orb;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HOST_CODE
SOLIDINLINE void G4Orb_ctor(G4Orb *This, G4double radius);

#ifdef ENABLE_VOXEL_NAVIGATION
SOLIDINLINE
G4bool G4Orb_CalculateExtent(
			   const G4Orb *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax);
#endif

#else

SOLIDINLINE EInside G4Orb_Inside(GEOMETRYLOC const G4Orb *This, G4ThreeVector p);

SOLIDINLINE G4ThreeVector G4Orb_SurfaceNormal(GEOMETRYLOC const G4Orb *This, G4ThreeVector p);


SOLIDINLINE G4double G4Orb_DistanceToIn_full(
				GEOMETRYLOC const G4Orb *This,
				G4ThreeVector p,
				G4ThreeVector v);

SOLIDINLINE G4double G4Orb_DistanceToIn(GEOMETRYLOC const G4Orb *This, G4ThreeVector p);

SOLIDINLINE G4double G4Orb_DistanceToOut_full(
			   GEOMETRYLOC const G4Orb *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n);

SOLIDINLINE G4double G4Orb_DistanceToOut(GEOMETRYLOC const G4Orb *This, G4ThreeVector p);
#endif

#if defined(INLINE_EVERYTHING) || defined(HOST_CODE)
#include "G4Orb.c"
#endif

#ifdef __cplusplus
}
#endif

#endif
