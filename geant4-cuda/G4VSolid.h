
/**
 * G4VSolid header
 * based on G4VSolid.hh of Geant 4.9.3
 */

#ifndef G4VSOLID_HH
#define G4VSOLID_HH

#include "everything.h"
#include "G4ThreeVector.h"

#define SOLIDINLINE INLINEFUNC


// EDIT: Moving definition of ESolid to everything.h
//typedef enum { kBox, kOrb, kTubs, kCons, kPolyCone } ESolid;


typedef struct G4VSolid
{
	ESolid type;
#ifdef DOUBLE_PRECISION
	int align;
#endif
}
G4VSolid;

#ifndef HOST_CODE

INLINEFUNC
EInside G4VSolid_Inside(GEOMETRYLOC const G4VSolid *This, G4ThreeVector p);
  // Returns kOutside if the point at offset p is outside the shapes
  // boundaries plus Tolerance/2, kSurface if the point is <= Tolerance/2
  // from a surface, otherwise kInside.

INLINEFUNC
G4ThreeVector G4VSolid_SurfaceNormal(GEOMETRYLOC const G4VSolid *This, G4ThreeVector p);
  // Returns the outwards pointing unit normal of the shape for the
  // surface closest to the point at offset p.

INLINEFUNC
G4double G4VSolid_DistanceToIn_full(
				GEOMETRYLOC const G4VSolid *This,
				G4ThreeVector p,
				G4ThreeVector v);
				
  // Return the distance along the normalised vector v to the shape,
  // from the point at offset p. If there is no intersection, return
  // kInfinity. The first intersection resulting from `leaving' a
  // surface/volume is discarded. Hence, it is tolerant of points on
  // the surface of the shape.

INLINEFUNC
G4double G4VSolid_DistanceToIn(GEOMETRYLOC const G4VSolid *This, G4ThreeVector p);
  // Calculate the distance to the nearest surface of a shape from an
  // outside point. The distance can be an underestimate.

INLINEFUNC
G4double G4VSolid_DistanceToOut_full(
			   GEOMETRYLOC const G4VSolid *This,
			   G4ThreeVector p,
			   G4ThreeVector v,
			   const G4bool calcNorm,
			   G4bool *validNorm,
			   G4ThreeVector *n);
  // Return the distance along the normalised vector v to the shape,
  // from a point at an offset p inside or on the surface of the shape.
  // Intersections with surfaces, when the point is < Tolerance/2 from a
  // surface must be ignored.
  // If calcNorm==true:
  //    validNorm set true if the solid lies entirely behind or on the
  //              exiting surface.
  //    n set to exiting outwards normal vector (undefined Magnitude).
  //    validNorm set to false if the solid does not lie entirely behind
  //              or on the exiting surface
  // If calcNorm==false:
  //    validNorm and n are unused.
  //
  // Must be called as solid.DistanceToOut(p,v) or by specifying all
  // the parameters.

INLINEFUNC
G4double G4VSolid_DistanceToOut(GEOMETRYLOC const G4VSolid *This, G4ThreeVector p);
  // Calculate the distance to the nearest surface of a shape from an
  // inside point. The distance can be an underestimate.
  
#else // not defined host code

#ifdef ENABLE_VOXEL_NAVIGATION

#include <vector>
typedef std::vector<G4ThreeVector> G4ThreeVectorList;

// Angle for mesh `wedges' in rads
// Works best when simple fraction of pi/2
#define K_GEOMETRY_MESH_ANGLE_DEFAULT (M_PI/4)

// Min wedges+1 to make
#define K_GEOMETRY_MIN_MESH_SECTIONS 3

// max wedges+1 to make
 // =>10 degrees/wedge for complete tube
#define K_GEOMETRY_MAX_MESH_SECTIONS 37

#include "G4Voxels.h"
#include "G4BuildVoxels.h"
#include "G4AffineTransform.h"

INLINEFUNC
G4bool G4VSolid_CalculateExtent(
			   const G4VSolid *This,
			   const EAxis pAxis,
			   G4VoxelLimits pVoxelLimit,
			   G4AffineTransform pTransform,
			   G4double* pMin, G4double* pMax);
  // Calculate the minimum and maximum extent of the solid, when under the
  // specified transform, and within the specified limits. If the solid
  // is not intersected by the region, return false, else return true.
  
#endif // voxel navigation
  
#ifdef INLINE_EVERYTHING
#include "G4VSolid.c"
#endif
  
#endif // host code

#endif // guard
