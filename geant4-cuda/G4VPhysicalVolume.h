
/**
 * G4VPhysicalVolume de-abstracted to equal G4PVPlacement
 */

#ifndef G4VPHYSICALVOLUME_H
#define G4VPHYSICALVOLUME_H

#include "everything.h"
#include "G4LogicalVolume.h"
#include "G4RotationMatrix.h"

typedef struct G4VPhysicalVolume
{
//public:  // with description
	
    G4RotationMatrix frot;
    G4ThreeVector ftrans;

//  private:
	//EDIT:
	int guard1;
    GEOMETRYLOC G4LogicalVolume *flogical;   // The logical volume representing the
                                 // physical and tracking attributes of
                                 // the volume
	int guard2;
	GEOMETRYLOC G4LogicalVolume   *flmother; // The current mother logical volume
	int guard3;
	int count;

#ifdef OPENCL_CODE
	int counter_shadow;// NOTE: This is not used. Exists only to keep sizes on GPU and CPU consistent.
#else
	static int counter;
#endif
/*
NOTE: The reason for the strange definitions of counter and counter_shadow above
is that in OpenCL ( version 1.1) it is not alowed to define a static member within
a struct. In order to keep the sizes of the PhysicalVolumes consistent on teh GPU and 
the CPU I had to define the counter_shadow. 
This can be modified if future releases allow for it.

*/
}
G4VPhysicalVolume;

#ifndef OPENCL_CODE
int G4VPhysicalVolume::counter = 0;
#endif


#include "G4VPhysicalVolume_inline.c"

#endif
