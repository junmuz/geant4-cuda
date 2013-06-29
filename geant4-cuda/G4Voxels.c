
/**
 * Device voxel handling inline implementations
 * based on G4SmartVoxel*.(i)cc of Geant 4.9.3
 */

#include "G4Voxels.h"

#ifndef G4VOXELS_C
#define G4VOXELS_C

INLINEFUNC
void G4VoxelNode_ctor( GEOMETRYLOC G4SmartVoxelNode *This, G4int no )
{
	This->fmaxEquivalent = no;
	This->fminEquivalent = no;
	This->fcontents = NULL;
	This->fNumContents = 0;
#ifdef NEW_NAVIGATION	
	for ( int i=0; i<Solidcount;i++)
		This->SolidType[i] = 0;
#endif
}

INLINEFUNC G4int
G4VoxelNode_GetNoContained(GEOMETRYLOC const G4SmartVoxelNode *This)
{
	return This->fNumContents;
}
	
INLINEFUNC G4int
G4VoxelNode_GetVolume(
	GEOMETRYLOC const G4SmartVoxelNode *This, G4int contentNo)
{
	myAssert( contentNo >= 0 && contentNo < This->fNumContents );
	return This->fcontents[contentNo];
}
	
INLINEFUNC G4int
G4VoxelNode_GetMaxEquivalentSliceNo(
	GEOMETRYLOC const G4SmartVoxelNode *This )
{
	return This->fmaxEquivalent;
}
	
INLINEFUNC G4int
G4VoxelNode_GetMinEquivalentSliceNo(
	GEOMETRYLOC const G4SmartVoxelNode *This )
{
	return This->fminEquivalent;
}
	
INLINEFUNC G4int
G4VoxelHeader_GetMaxEquivalentSliceNo(
	GEOMETRYLOC const G4SmartVoxelHeader *This )
{
	return This->fmaxEquivalent;
}
	
INLINEFUNC G4int
G4VoxelHeader_GetMinEquivalentSliceNo(
	GEOMETRYLOC const G4SmartVoxelHeader *This )
{
	return This->fminEquivalent;
}
	
INLINEFUNC EAxis
G4VoxelHeader_GetAxis( GEOMETRYLOC const G4SmartVoxelHeader *This )
{
	return This->faxis;
}

INLINEFUNC G4int
G4VoxelHeader_GetNoSlices( GEOMETRYLOC const G4SmartVoxelHeader *This )
{
	return This->fNumSlices;
}

INLINEFUNC G4double
G4VoxelHeader_GetMinExtent( GEOMETRYLOC const G4SmartVoxelHeader *This )
{
	return This->fminExtent;
}

INLINEFUNC G4double
G4VoxelHeader_GetMaxExtent( GEOMETRYLOC const G4SmartVoxelHeader *This )
{
	return This->fmaxExtent;
}

INLINEFUNC GEOMETRYLOC G4SmartVoxelProxy*
G4VoxelHeader_GetSlice( GEOMETRYLOC const G4SmartVoxelHeader *This, G4int n )
{
	myAssert( n >= 0 && n < This->fNumSlices );
	return This->fslices[n];
}

INLINEFUNC G4bool
G4VoxelProxy_IsNode( GEOMETRYLOC const G4SmartVoxelProxy *This )
{
	return This->fNode != GEOMETRYNULL;
}

INLINEFUNC G4bool
G4VoxelProxy_IsHeader( GEOMETRYLOC const G4SmartVoxelProxy *This )
{
	return This->fHeader != GEOMETRYNULL;
}

INLINEFUNC GEOMETRYLOC G4SmartVoxelNode*
G4VoxelProxy_GetNode( GEOMETRYLOC const G4SmartVoxelProxy *This )
{
	return This->fNode;
}

INLINEFUNC GEOMETRYLOC G4SmartVoxelHeader*
G4VoxelProxy_GetHeader( GEOMETRYLOC const G4SmartVoxelProxy *This )
{
	return This->fHeader;
}

#ifdef HOST_CODE

#include <cstdlib>

INLINEFUNC void
G4VoxelNode_SetMaxEquivalentSliceNo( 
	G4SmartVoxelNode *This, G4int val )
{
	This->fmaxEquivalent = val;
}
	
INLINEFUNC void
G4VoxelNode_SetMinEquivalentSliceNo(
	G4SmartVoxelNode *This, G4int val )
{
	This->fminEquivalent = val;
}

INLINEFUNC void
G4VoxelHeader_SetMaxEquivalentSliceNo( 
	G4SmartVoxelHeader *This, G4int val )
{
	This->fmaxEquivalent = val;
}
	
INLINEFUNC void
G4VoxelHeader_SetMinEquivalentSliceNo(
	G4SmartVoxelHeader *This, G4int val )
{
	This->fminEquivalent = val;
}

INLINEFUNC void
G4VoxelHeader_SetMinExtent( G4SmartVoxelHeader *This, G4double val )
{
	This->fminExtent = val;
}

INLINEFUNC void
G4VoxelHeader_SetMaxExtent( G4SmartVoxelHeader *This, G4double val )
{
	This->fmaxExtent = val;
}

	
INLINEFUNC void
G4VoxelHeader_SetAxis( G4SmartVoxelHeader *This, EAxis val )
{
	This->faxis = val;
}


INLINEFUNC G4bool G4VoxelNode_operator_equal(
	const G4SmartVoxelNode *This,
	const G4SmartVoxelNode *v )
{
  G4int maxNode=G4VoxelNode_GetNoContained(This);
  if (maxNode==G4VoxelNode_GetNoContained(v))
  {
    for (G4int node=0;node<maxNode;node++)
    {
      if (G4VoxelNode_GetVolume(This,node)!=G4VoxelNode_GetVolume(v,node))
      {
        return false;
      }
    }
    return true;
  }
  return false;
}

INLINEFUNC
G4bool G4VoxelProxy_operator_equal (
		const G4SmartVoxelProxy *This,
		const G4SmartVoxelProxy *v)
{
  return (This==v) ? true : false;
}


#ifdef NEW_NAVIGATION
INLINEFUNC
void G4VoxelNode_Insert( G4SmartVoxelNode *This, G4int thing, int type)
{	static int i;
	This->fNumContents++;
	This->fcontents = (G4int*)std::realloc( This->fcontents, sizeof(G4int)*This->fNumContents );
	This->fcontents[This->fNumContents-1] = thing;

	(This->SolidType[type])++;
	//REMOVE
	//if(i<20)
		//std::cout<<" From insert voxel node, SolidType count is =" << This->SolidType[type]<<" for the solid of type = "<<type<< " and the fNumContents is "<<This->fNumContents<< " and this is "<<This<< std::endl;
	
	if(type==0)
		std::cout<<" Box\n";
	
	i++;
}
#else
void G4VoxelNode_Insert( G4SmartVoxelNode *This, G4int thing)
{
	This->fNumContents++;
	This->fcontents = (G4int*)std::realloc( This->fcontents, sizeof(G4int)*This->fNumContents );
	This->fcontents[This->fNumContents-1] = thing;
}
#endif



INLINEFUNC
void G4VoxelNode_dtor( G4SmartVoxelNode *This )
{
	std::free( This->fcontents );
}

INLINEFUNC
void G4VoxelProxy_ctor_node( G4SmartVoxelProxy *This, G4SmartVoxelNode *n )
{
	This->fNode = n;
	This->fHeader = NULL;
}

INLINEFUNC
void G4VoxelProxy_ctor_header( G4SmartVoxelProxy *This, G4SmartVoxelHeader *h )
{
	This->fNode = NULL;
	This->fHeader = h;
}

#endif

#endif
