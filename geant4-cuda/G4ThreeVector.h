
/**
 * Three-vector data structure
 */

#ifndef G4THREEVECTOR_HH
#define G4THREEVECTOR_HH

#include "everything.h"

#if defined(GPU)
	#define THREEVECTOR_EXTRAMEMBER
#endif

typedef struct
{
	G4double x,y,z;
	#ifdef THREEVECTOR_EXTRAMEMBER
	G4double w;
	#endif
}
G4ThreeVector;

#include "G4ThreeVector_inline.c"

#endif
