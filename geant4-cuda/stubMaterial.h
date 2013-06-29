
/** Dummy material */

#ifndef __STUBMATERIAL_HH__
#define __STUBMATERIAL_HH__

#include "everything.h"

typedef struct //__attribute__((packed))
{
	G4double property;
}
StubMaterial;

#define G4Material StubMaterial

#endif
