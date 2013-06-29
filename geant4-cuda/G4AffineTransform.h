
/**
 * Affine transform data structure based on a C++ class in CLHEP
 */

#ifndef G4AFFINETRANSFORM_H
#define G4AFFINETRANSFORM_H

#include "everything.h"
#include "G4ThreeVector.h"
#include "G4RotationMatrix.h"

typedef struct
{
  G4double rxx,rxy,rxz;
  G4double ryx,ryy,ryz;
  G4double rzx,rzy,rzz;
  G4double tx,ty,tz;
}
G4AffineTransform;

#include "G4AffineTransform_inline.c"

#endif
