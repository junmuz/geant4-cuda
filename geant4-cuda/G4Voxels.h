
/** Classes needed for voxel handling on the device */

#ifndef G4VOXELS_H
#define G4VOXELS_H

#include "everything.h"

struct G4SmartVoxelProxy;

typedef struct
{
	G4double fmaxExtent;
	G4double fminExtent;
	  // Max and min coordinate along faxis.
	  
	GEOMETRYLOC struct G4SmartVoxelProxy* GEOMETRYLOC * fslices;
	  // Slices along axis. 
	
	G4int fNumSlices;

	G4int fminEquivalent;
	G4int fmaxEquivalent;
	  // Min and max equivalent slice nos for previous level.
	  
	EAxis faxis;
	EAxis fparamAxis;
	  // Axis for slices.
}
G4SmartVoxelHeader;

#ifndef NEW_NAVIGATION
typedef struct
{
	GEOMETRYLOC G4int *fcontents;
	
	G4int fminEquivalent;
	G4int fmaxEquivalent;
	
	G4int fNumContents;
}
G4SmartVoxelNode;

#else
// Adding support for new navigation needs new type of voxels.

typedef struct
{
	GEOMETRYLOC G4int *fcontents;
	
	G4int fminEquivalent;
	G4int fmaxEquivalent;
	
	G4int fNumContents;
	
	G4int SolidType[ Solidcount ];
		// Array with each element storing the sum of solids of a certain type.

}
G4SmartVoxelNode;
#endif

typedef struct G4SmartVoxelProxy
{
	GEOMETRYLOC G4SmartVoxelHeader* fHeader;
    GEOMETRYLOC G4SmartVoxelNode* fNode;
}
G4SmartVoxelProxy;




#include "G4Voxels.c"

#endif
