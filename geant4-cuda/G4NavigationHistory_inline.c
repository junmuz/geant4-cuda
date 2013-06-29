
/**
 * G4NavigationHistory (inline) implementation,
 * based on G4NavigationHistory.* and G4NavigationLevel.* of Geant 4.9.3 
 */

INLINEFUNC
void G4NavigationLevel_ctor(
			G4NavigationLevel *This,
			GEOMETRYLOC G4VPhysicalVolume* pPhysVol,
			G4AffineTransform afTransform,
			EVolume            volTp )
{
	This->fTransform = afTransform;
	This->fPhysicalVolumePtr = pPhysVol;
	This->fVolumeType = volTp;
}

INLINEFUNC
void G4NavigationLevel_ctor_relative(
			G4NavigationLevel *This,
			GEOMETRYLOC G4VPhysicalVolume* pPhysVol,
			G4AffineTransform levelAbove,
			G4AffineTransform relativeCurrent,
			EVolume            volTp )
{
	This->fPhysicalVolumePtr = pPhysVol;
	This->fVolumeType = volTp;
	G4AffineTransform_InverseProduct(&(This->fTransform), &levelAbove, &relativeCurrent );
}

INLINEFUNC
G4NavigationLevel G4NavigationLevel_create(
			GEOMETRYLOC G4VPhysicalVolume* pPhysVol,
			G4AffineTransform afTransform,
			EVolume            volTp )
{
	G4NavigationLevel lev;
	G4NavigationLevel_ctor( &lev, pPhysVol, afTransform, volTp );
	return lev;
}

INLINEFUNC
G4NavigationLevel G4NavigationLevel_create_relative(
	GEOMETRYLOC G4VPhysicalVolume* pPhysVol,
	G4AffineTransform levelAbove,
	G4AffineTransform relativeCurrent,
	EVolume            volTp)
{
	G4NavigationLevel lev;
	G4NavigationLevel_ctor_relative( &lev, pPhysVol, levelAbove, relativeCurrent, volTp );
	return lev;
}


INLINEFUNC
GEOMETRYLOC G4VPhysicalVolume* G4NavigationLevel_GetPhysicalVolume(
	const G4NavigationLevel *This ) 
{ 
  return This->fPhysicalVolumePtr;
}

INLINEFUNC
G4AffineTransform G4NavigationLevel_GetTransform(
	const G4NavigationLevel *This )
{ 
  return This->fTransform;
} 

INLINEFUNC
const G4AffineTransform* G4NavigationLevel_GetPtrTransform(
	const G4NavigationLevel *This )
{ 
  return &(This->fTransform);
} 

INLINEFUNC
EVolume G4NavigationLevel_GetVolumeType(
	const G4NavigationLevel *This )
{ 
  return This->fVolumeType;
}

// class G4NavigationHistory Inline implementation
//
// ----------------------------------------------------------------------

INLINEFUNC
void G4NavigationHistory_Reset( G4NavigationHistory *This )
{
	This->fStackDepth = 0;
}

INLINEFUNC
void G4NavigationHistory_Clear( G4NavigationHistory *This )
{
  G4AffineTransform origin = G4AffineTransform_create_vector(G4ThreeVector_create(0.,0.,0.));
  G4NavigationLevel tmpNavLevel = G4NavigationLevel_create(0, origin, kNormal) ;

  G4NavigationHistory_Reset( This );
  for (G4int ilev=K_NAVIGATION_HISTORY_DEPTH-1; ilev>=0; ilev--)
  {
     This->fNavHistory[ilev] = tmpNavLevel;
  }
}

INLINEFUNC
void G4NavigationHistory_ctor( G4NavigationHistory *This )
{
	This->fStackDepth = 0;
	G4NavigationHistory_Clear( This );
}

INLINEFUNC
void G4NavigationHistory_dtor( G4NavigationHistory *This )
{
	(void)This;
}



INLINEFUNC
void G4NavigationHistory_SetFirstEntry(
	G4NavigationHistory *This, GEOMETRYLOC G4VPhysicalVolume* pVol)
{
  G4ThreeVector translation = G4ThreeVector_create(0.,0.,0.);

  // Protection needed in case pVol=null 
  // so that a touchable-history can signal OutOfWorld 
  //
  if( pVol!=GEOMETRYNULL )
  {
    translation = G4VPhysicalVolume_GetTranslation( pVol );
  }
  This->fNavHistory[0] =
    G4NavigationLevel_create( pVol, G4AffineTransform_create_vector(translation), kNormal );
}

INLINEFUNC
const G4AffineTransform* G4NavigationHistory_GetPtrTopTransform(
	const G4NavigationHistory *This )
{
  return G4NavigationLevel_GetPtrTransform( &(This->fNavHistory[This->fStackDepth]) );
}

INLINEFUNC
G4AffineTransform G4NavigationHistory_GetTopTransform(
	const G4NavigationHistory *This )
{
  return G4NavigationLevel_GetTransform( &(This->fNavHistory[This->fStackDepth]) );
}

INLINEFUNC
EVolume G4NavigationHistory_GetTopVolumeType(
	const G4NavigationHistory *This )
{
  return G4NavigationLevel_GetVolumeType( &(This->fNavHistory[This->fStackDepth]) );
}

INLINEFUNC
GEOMETRYLOC G4VPhysicalVolume* G4NavigationHistory_GetTopVolume(
	const G4NavigationHistory *This )
{
  return G4NavigationLevel_GetPhysicalVolume( &(This->fNavHistory[This->fStackDepth]) );
}

INLINEFUNC
G4int G4NavigationHistory_GetDepth(
	const G4NavigationHistory *This )
{
  return This->fStackDepth;
}

INLINEFUNC
G4AffineTransform
G4NavigationHistory_GetTransform(
	const G4NavigationHistory *This, G4int n )
{
  return G4NavigationLevel_GetTransform( &(This->fNavHistory[n]) );
}

INLINEFUNC
EVolume G4NavigationHistory_GetVolumeType(
	const G4NavigationHistory *This, G4int n )
{
  return G4NavigationLevel_GetVolumeType( &(This->fNavHistory[n]) );
}

INLINEFUNC
GEOMETRYLOC G4VPhysicalVolume* G4NavigationHistory_GetVolume(
	const G4NavigationHistory *This, G4int n )
{
  return G4NavigationLevel_GetPhysicalVolume( &(This->fNavHistory[n]) );
}

INLINEFUNC
G4int G4NavigationHistory_GetMaxDepth(
	const G4NavigationHistory *This )
{
	(void)This;
	return K_NAVIGATION_HISTORY_DEPTH;
}

INLINEFUNC
void G4NavigationHistory_BackLevel( G4NavigationHistory *This )
{
  myAssert( This->fStackDepth>0 );

  // Tell  the  level  that I am forgetting it
  // delete fNavHistory(fStackDepth);
  //
  This->fStackDepth--;
}

/*INLINEFUNC
void G4NavigationHistory_BackLevel_n(G4NavigationHistory *This, G4int n)
{
  myAssert( n<=This->fStackDepth );
  This->fStackDepth-=n;
}*/

/*static
void G4NavigationHistory::EnlargeHistory()
{
  G4int len = fNavHistory.size();
  if ( len==fStackDepth )
  {
    // Note: Resize operation clears additional entries
    //
    G4int nlen = len+kHistoryStride;
    fNavHistory.resize(nlen);
  }  
}
*/

INLINEFUNC
void G4NavigationHistory_NewLevel(
		G4NavigationHistory *This,
		GEOMETRYLOC G4VPhysicalVolume *pNewMother,
		EVolume vType )
{
  This->fStackDepth++;
  
  myAssert(This->fStackDepth < K_NAVIGATION_HISTORY_DEPTH);
  
  This->fNavHistory[This->fStackDepth] =
    G4NavigationLevel_create_relative(
			pNewMother, 
			G4NavigationLevel_GetTransform( &(This->fNavHistory[This->fStackDepth-1]) ),
			G4AffineTransform_create_full(
				G4VPhysicalVolume_GetObjectRotationValue( pNewMother ),
				G4VPhysicalVolume_GetTranslation( pNewMother )),
			vType ); 
  // The constructor computes the new global->local transform
}

