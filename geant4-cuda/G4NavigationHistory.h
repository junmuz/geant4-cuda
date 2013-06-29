
/**
 * G4NavigationHistory header,
 * based on G4NavigationHistory.hh of Geant 4.9.3
 */

#ifndef G4NAVIGATIONHISTORY_H
#define G4NAVIGATIONHISTORY_H

#include "everything.h"

#include "G4AffineTransform.h"
#include "G4VPhysicalVolume.h"

#define K_NAVIGATION_HISTORY_DEPTH 16
// defined in Makefile
//#define K_NAVIGATION_HISTORY_DEPTH 16 // pick an even number

typedef struct
{
//private:

   G4AffineTransform  fTransform;
     // Compounded global->local transformation (takes a point in the 
     // global reference system to the system of the volume at this level)

   GEOMETRYLOC G4VPhysicalVolume* fPhysicalVolumePtr;
     // Physical volume ptrs, for this level's volume

   EVolume fVolumeType;
     // Volume `type' 
}
G4NavigationLevel;


typedef struct
{
	G4NavigationLevel fNavHistory[K_NAVIGATION_HISTORY_DEPTH];

	G4int fStackDepth;
	// Depth of stack: effectively depth in geometrical tree
	
	int align;
}
G4NavigationHistory;

#include "G4NavigationHistory_inline.c"

#endif
