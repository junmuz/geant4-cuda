
/**
 * Switch functions to replace virtual function calls in G4VSolid
 */

#ifndef G4VSOLID_INLINE_C
#define G4VSOLID_INLINE_C

#define ENABLE_G4BOX
#ifndef NO_ORB
#define ENABLE_G4ORB
#endif
#ifndef ONLY_BOX_AND_ORB
#define ENABLE_G4TUBS
#define ENABLE_G4CONS
//#define ENABLE_G4POLYCONE
#endif

#ifdef ENABLE_G4BOX
	#include "G4Box.h"
#endif
#ifdef ENABLE_G4ORB
	#include "G4Orb.h"
#endif
#ifdef ENABLE_G4TUBS
	#include "G4Tubs.h"
#endif
#ifdef ENABLE_G4CONS
	#include "G4Cons.h"
#endif
#ifdef ENABLE_G4POLYCONE
	#include "G4PolyCone.h"
#endif

#ifndef HOST_CODE

INLINEFUNC
EInside G4VSolid_Inside(GEOMETRYLOC const G4VSolid *This, G4ThreeVector p)
{
	switch(This->type)
	{
#ifdef ENABLE_G4BOX
		case kBox:
			return G4Box_Inside((GEOMETRYLOC const G4Box*)This,p);
#endif
#ifdef ENABLE_G4ORB
		case kOrb:
			return G4Orb_Inside((GEOMETRYLOC const G4Orb*)This,p);
#endif
#ifdef ENABLE_G4TUBS
		case kTubs:
			return G4Tubs_Inside((GEOMETRYLOC const G4Tubs*)This,p);
#endif
#ifdef ENABLE_G4CONS
		case kCons:
			return G4Cons_Inside((GEOMETRYLOC const G4Cons*)This,p);
#endif
#ifdef ENABLE_G4POLYCONE
		case kPolyCone:
			return G4PolyCone_Inside((GEOMETRYLOC const G4PolyCone*)This,p);
#endif

		default:
			myAssert(false);
			return kOutside;
	}
}

INLINEFUNC
G4ThreeVector G4VSolid_SurfaceNormal(GEOMETRYLOC const G4VSolid *This, G4ThreeVector p)
{
	switch(This->type)
	{
#ifdef ENABLE_G4BOX
		case kBox:
			return G4Box_SurfaceNormal((GEOMETRYLOC const G4Box*)This,p);
#endif
#ifdef ENABLE_G4ORB
		case kOrb:
			return G4Orb_SurfaceNormal((GEOMETRYLOC const G4Orb*)This,p);
#endif
#ifdef ENABLE_G4TUBS
		case kTubs:
			return G4Tubs_SurfaceNormal((GEOMETRYLOC const G4Tubs*)This,p);
#endif
#ifdef ENABLE_G4CONS
		case kCons:
			return G4Cons_SurfaceNormal((GEOMETRYLOC const G4Cons*)This,p);
#endif
#ifdef ENABLE_G4POLYCONE
		case kPolyCone:
			return G4PolyCone_SurfaceNormal((GEOMETRYLOC const G4PolyCone*)This,p);
#endif
		default:
			myAssert(false);
			return G4ThreeVector_create(0,0,0);
	}
}

INLINEFUNC
G4double G4VSolid_DistanceToIn_full(
				GEOMETRYLOC const G4VSolid *This,
				G4ThreeVector p,
				G4ThreeVector v)
{
	switch(This->type)
	{
#ifdef ENABLE_G4BOX
		case kBox:
			return G4Box_DistanceToIn_full((GEOMETRYLOC const G4Box*)This,p,v);
#endif
#ifdef ENABLE_G4ORB
		case kOrb:
			return G4Orb_DistanceToIn_full((GEOMETRYLOC const G4Orb*)This,p,v);
#endif
#ifdef ENABLE_G4TUBS
		case kTubs:
			return G4Tubs_DistanceToIn_full((GEOMETRYLOC const G4Tubs*)This,p,v);
#endif
#ifdef ENABLE_G4CONS
		case kCons:
			return G4Cons_DistanceToIn_full((GEOMETRYLOC const G4Cons*)This,p,v);
#endif
#ifdef ENABLE_G4POLYCONE
		case kPolyCone:
			return G4PolyCone_DistanceToIn_full((GEOMETRYLOC const G4PolyCone*)This,p,v);
#endif
		default:
			myAssert(false);
			return 0;
	}
}

INLINEFUNC
G4double G4VSolid_DistanceToIn(GEOMETRYLOC const G4VSolid *This, G4ThreeVector p)
{
	switch(This->type)
	{
#ifdef ENABLE_G4BOX
		case kBox:
			return G4Box_DistanceToIn((GEOMETRYLOC const G4Box*)This,p);
#endif
#ifdef ENABLE_G4ORB
		case kOrb:
			return G4Orb_DistanceToIn((GEOMETRYLOC const G4Orb*)This,p);
#endif
#ifdef ENABLE_G4TUBS
		case kTubs:
			return G4Tubs_DistanceToIn((GEOMETRYLOC const G4Tubs*)This,p);
#endif
#ifdef ENABLE_G4CONS
		case kCons:
			return G4Cons_DistanceToIn((GEOMETRYLOC const G4Cons*)This,p);
#endif
#ifdef ENABLE_G4POLYCONE
		case kPolyCone:
			return G4PolyCone_DistanceToIn((GEOMETRYLOC const G4PolyCone*)This,p);
#endif
		default:
			myAssert(false);
			return 0;
	}
}

INLINEFUNC
G4double G4VSolid_DistanceToOut_full(
			   GEOMETRYLOC const G4VSolid *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n)
{
	switch(This->type)
	{
#ifdef ENABLE_G4BOX
		case kBox:
			return G4Box_DistanceToOut_full((GEOMETRYLOC const G4Box*)This,p,v,calcNorm,validNorm,n);
#endif
#ifdef ENABLE_G4ORB
		case kOrb:
			return G4Orb_DistanceToOut_full((GEOMETRYLOC const G4Orb*)This,p,v,calcNorm,validNorm,n);
#endif
#ifdef ENABLE_G4TUBS
		case kTubs:
			return G4Tubs_DistanceToOut_full((GEOMETRYLOC const G4Tubs*)This,p,v,calcNorm,validNorm,n);
#endif
#ifdef ENABLE_G4CONS
		case kCons:
			return G4Cons_DistanceToOut_full((GEOMETRYLOC const G4Cons*)This,p,v,calcNorm,validNorm,n);
#endif
#ifdef ENABLE_G4POLYCONE
		case kPolyCone:
			return G4PolyCone_DistanceToOut_full((GEOMETRYLOC const G4PolyCone*)This,p,v,calcNorm,validNorm,n);
#endif
		default:
			myAssert(false);
			return 0;
	}
}

INLINEFUNC
G4double G4VSolid_DistanceToOut(GEOMETRYLOC const G4VSolid *This, G4ThreeVector p)
{
	switch(This->type)
	{
#ifdef ENABLE_G4BOX
		case kBox:
			return G4Box_DistanceToOut((GEOMETRYLOC const G4Box*)This,p);
#endif
#ifdef ENABLE_G4ORB
		case kOrb:
			return G4Orb_DistanceToOut((GEOMETRYLOC const G4Orb*)This,p);
#endif
#ifdef ENABLE_G4TUBS
		case kTubs:
			return G4Tubs_DistanceToOut((GEOMETRYLOC const G4Tubs*)This,p);
#endif
#ifdef ENABLE_G4CONS
		case kCons:
			return G4Cons_DistanceToOut((GEOMETRYLOC const G4Cons*)This,p);
#endif
#ifdef ENABLE_G4POLYCONE
		case kPolyCone:
			return G4PolyCone_DistanceToOut((GEOMETRYLOC const G4PolyCone*)This,p);
#endif
		default:
			myAssert(false);
			return 0;
	}
}

#else

#ifdef ENABLE_VOXEL_NAVIGATION
INLINEFUNC
G4bool G4VSolid_CalculateExtent(
			   const G4VSolid *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax)
{
	switch(This->type)
	{
#ifdef ENABLE_G4BOX
		case kBox:
			return G4Box_CalculateExtent(
				(GEOMETRYLOC const G4Box*)This, pAxis,
				pVoxelLimit, pTransform, pMin, pMax);
#endif
#ifdef ENABLE_G4ORB
		case kOrb:
			return G4Orb_CalculateExtent(
				(GEOMETRYLOC const G4Orb*)This, pAxis,
				pVoxelLimit, pTransform, pMin, pMax);
#endif
#ifdef ENABLE_G4TUBS
		case kTubs:
			return G4Tubs_CalculateExtent(
				(GEOMETRYLOC const G4Tubs*)This, pAxis,
				pVoxelLimit, pTransform, pMin, pMax);
#endif
#ifdef ENABLE_G4CONS
		case kCons:
			return G4Cons_CalculateExtent(
				(GEOMETRYLOC const G4Cons*)This, pAxis,
				pVoxelLimit, pTransform, pMin, pMax);
#endif
#ifdef ENABLE_G4POLYCONE
		case kPolyCone:
			return G4PolyCone_CalculateExtent(
				(GEOMETRYLOC const G4PolyCone*)This, pAxis,
				pVoxelLimit, pTransform, pMin, pMax);
#endif
		default:
			myAssert(false);
			return 0;
	}
}
#endif

#endif

#endif
