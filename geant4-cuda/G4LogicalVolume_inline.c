
/**
 * G4LogicalVolume inline implementation,
 * based on G4LogicalVolume.icc of Geant 4.9.3 
 */

#ifndef G4LOGICALVOLUME_INLINE_C
#define G4LOGICALVOLUME_INLINE_C

#ifdef ENABLE_VOXEL_NAVIGATION

// ********************************************************************
// GetVoxelHeader
// ********************************************************************
//
INLINEFUNC
GEOMETRYLOC G4SmartVoxelHeader * G4LogicalVolume_GetVoxelHeader(GEOMETRYLOC const G4LogicalVolume* This)
{
	return This->fVoxel;
}

#endif

// ********************************************************************
// GetNoDaughters
// ********************************************************************
//
INLINEFUNC
G4int G4LogicalVolume_GetNoDaughters(GEOMETRYLOC const G4LogicalVolume* This)
{
  return This->fNoDaughters;
}

// ********************************************************************
// GetDaughter
// ********************************************************************
//
INLINEFUNC
GEOMETRYLOC struct G4VPhysicalVolume* G4LogicalVolume_GetDaughter(GEOMETRYLOC const G4LogicalVolume* This, const G4int i)
{
  return This->fDaughters[i];
}

// ********************************************************************
// GetSolid
// ********************************************************************
//
INLINEFUNC
GEOMETRYLOC struct G4VSolid* G4LogicalVolume_GetSolid(GEOMETRYLOC const G4LogicalVolume* This)
{
  return This->fSolid;
}

// ********************************************************************
// GetMaterial
// ********************************************************************
//
INLINEFUNC
GEOMETRYLOC G4Material* G4LogicalVolume_GetMaterial(GEOMETRYLOC const G4LogicalVolume* This)
{
  return This->fMaterial;
}

#ifdef HOST_CODE

#ifdef ENABLE_VOXEL_NAVIGATION
// ********************************************************************
// SetVoxelHeader
// ********************************************************************
//
INLINEFUNC
void G4LogicalVolume_SetVoxelHeader(G4LogicalVolume* This, G4SmartVoxelHeader * pVoxel)
{
  This->fVoxel = pVoxel;
}
#endif

INLINEFUNC
void G4LogicalVolume_AddDaughter(G4LogicalVolume* This, struct G4VPhysicalVolume* pNewDaughter)
{
	This->fNoDaughters++;
	
	// WARNING: the data structure must be allocated manually
	
	// This->fDaughters = realloc( This->fDaughters, sizeof(struct G4VPhysicalVolume*)*This->fNoDaughters );
	This->fDaughters[This->fNoDaughters-1] = pNewDaughter;
}

INLINEFUNC
void
G4LogicalVolume_ctor(
	G4LogicalVolume *This, struct G4VSolid* pSolid, G4Material* pMaterial )
{
	assert( pSolid != NULL );
	assert( pMaterial != NULL );
	
	This->fDaughters = NULL;
	This->fNoDaughters = 0;
	This->fSolid = pSolid;
	This->fMaterial = pMaterial;
#ifdef ENABLE_VOXEL_NAVIGATION
	This->fVoxel = NULL;
#endif
}

INLINEFUNC
void
G4LogicalVolume_dtor( G4LogicalVolume *This )
{
	(void)This;
	//free(This->fDaughters);
}

#endif

#endif
