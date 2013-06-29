
/** Contains things needed for building voxelizations on the host */

#ifndef BUILD_VOXELS_H
#define BUILD_VOXELS_H

#include "G4Voxels.h"

typedef struct
{
    G4double fxAxisMin,fxAxisMax;
    G4double fyAxisMin,fyAxisMax;
    G4double fzAxisMin,fzAxisMax;
}
G4VoxelLimits;

#ifdef INLINE_EVERYTHING
#include "G4VoxelLimits.c"
#include "G4VoxelHeader.cpp"
#endif

#endif
