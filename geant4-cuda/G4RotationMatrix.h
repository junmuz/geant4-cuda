
/**
 * Rotation matrix data structure based on a C++ class in CLHEP
 */

#ifndef G4ROTATIONMATRIX_H
#define G4ROTATIONMATRIX_H

#include "everything.h"

typedef struct
{
	G4double
		rxx, rxy, rxz, 
		ryx, ryy, ryz, 
		rzx, rzy, rzz;

//#ifndef DOUBLE_PRECISION
	G4double align;
//#endif
}
G4RotationMatrix;

#include "G4RotationMatrix_inline.c"

#endif
