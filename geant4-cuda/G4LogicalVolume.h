
/** G4LogicalVolume header, based on G4LogicalVolume.hh of Geant 4.9.3 */

#ifndef G4LOGICALVOLUME_H
#define G4LOGICALVOLUME_H

#include "everything.h"

#include "stubMaterial.h"

#ifdef ENABLE_VOXEL_NAVIGATION
#include "G4Voxels.h"
#endif

struct G4VPhysicalVolume;
struct G4VSolid;

typedef struct
{
	G4int fNoDaughters;
	
	// ALLOCATE MANUALLY !!!!
	GEOMETRYLOC struct G4VPhysicalVolume * GEOMETRYLOC *fDaughters;

	  int check;
	
	  GEOMETRYLOC G4Material* fMaterial;
	// Pointer to material at this node.
	GEOMETRYLOC struct G4VSolid* fSolid;
	// Pointer to solid.
    
#ifdef ENABLE_VOXEL_NAVIGATION
	GEOMETRYLOC G4SmartVoxelHeader *fVoxel;
	int align;
#endif
}
G4LogicalVolume;

#include "G4VSolid.h"
#include "G4VSolid_inline.c"

#include "G4LogicalVolume_inline.c"

#endif
