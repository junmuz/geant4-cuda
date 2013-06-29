
/** G4Box header, based on G4Box.hh of Geant 4.9.3 */

// --------------------------------------------------------------------
#ifndef G4BOX_HH
#define G4BOX_HH

#include "G4VSolid.h"

typedef struct
{
	G4VSolid solid;
    G4double fDx,fDy,fDz;
}
G4Box;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HOST_CODE
SOLIDINLINE void G4Box_ctor(G4Box *This, G4double x, G4double y, G4double z);

#ifdef ENABLE_VOXEL_NAVIGATION
SOLIDINLINE
G4bool G4Box_CalculateExtent(
			   const G4Box *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax);
#endif
			   
#else
SOLIDINLINE EInside G4Box_Inside(GEOMETRYLOC const G4Box *This, G4ThreeVector p);

SOLIDINLINE G4ThreeVector G4Box_SurfaceNormal(GEOMETRYLOC const G4Box *This, G4ThreeVector p);


SOLIDINLINE G4double G4Box_DistanceToIn_full(
				GEOMETRYLOC const G4Box *This,
				G4ThreeVector p,
				G4ThreeVector v);

SOLIDINLINE G4double G4Box_DistanceToIn(GEOMETRYLOC const G4Box *This, G4ThreeVector p);

SOLIDINLINE G4double G4Box_DistanceToOut_full(
			   GEOMETRYLOC const G4Box *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n);

SOLIDINLINE G4double G4Box_DistanceToOut(GEOMETRYLOC const G4Box *This, G4ThreeVector p);
#endif

#if defined(INLINE_EVERYTHING) || defined(HOST_CODE)
#include "G4Box.c"
#endif

#ifdef __cplusplus
}
#endif

#endif
