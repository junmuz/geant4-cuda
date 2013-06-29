
/**
 * G4VPhysicalVolume inline implementation
 * based on G4VPhysicalVolume.icc of Geant 4.9.3
 */

INLINEFUNC
G4ThreeVector G4VPhysicalVolume_GetTranslation(GEOMETRYLOC const G4VPhysicalVolume *This)
{
  return This->ftrans;
}

/*INLINEFUNC
const G4RotationMatrix* G4VPhysicalVolume_GetRotation(GEOMETRYLOC const G4VPhysicalVolume *This)
{
  return &(This->frot);
}

INLINEFUNC
G4RotationMatrix* G4VPhysicalVolume_GetRotation_nonconst(GEOMETRYLOC G4VPhysicalVolume *This)
{
  return &(This->frot);
}*/

INLINEFUNC
GEOMETRYLOC G4LogicalVolume* G4VPhysicalVolume_GetLogicalVolume(GEOMETRYLOC const G4VPhysicalVolume *This)
{
  return This->flogical;
}

INLINEFUNC
GEOMETRYLOC G4LogicalVolume* G4VPhysicalVolume_GetMotherLogical(GEOMETRYLOC const G4VPhysicalVolume *This)
{
  return This->flmother;
}

INLINEFUNC
G4RotationMatrix G4VPhysicalVolume_GetObjectRotationValue(GEOMETRYLOC const G4VPhysicalVolume *This)
{
  return This->frot;
}

INLINEFUNC
G4ThreeVector G4VPhysicalVolume_GetObjectTranslation(GEOMETRYLOC const G4VPhysicalVolume *This)
{
	return This->ftrans;
}

#ifdef HOST_CODE

INLINEFUNC
void G4VPhysicalVolume_SetLogicalVolume(G4VPhysicalVolume *This, G4LogicalVolume *pLogical)
{
  This->flogical=pLogical;
}

INLINEFUNC
void G4VPhysicalVolume_SetMotherLogical(G4VPhysicalVolume *This, G4LogicalVolume *pMother)
{
  This->flmother=pMother;
}

// Constructor: init parameters
//

INLINEFUNC
void G4VPhysicalVolume_ctor(
		G4VPhysicalVolume *This,
		G4RotationMatrix pRot,
		G4ThreeVector tlate,
		GEOMETRYLOC G4LogicalVolume* pLogical )
{	
	This->frot = pRot;
	This->ftrans = tlate;
	This->flogical = pLogical;
	This->flmother = NULL;
	
	
}

#endif // defined(HOST_CODE)
