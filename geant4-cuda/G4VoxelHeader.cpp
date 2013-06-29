
/**
 * G4VoxelHeader (host side) implementation
 * based on G4VoxelHeader.cc of Geant 4.9.3
 */

#ifndef G4VOXELHEADER_CPP
#define G4VOXELHEADER_CPP

#include "G4BuildVoxels.h"
#include "G4LogicalVolume.h"
#include "G4VPhysicalVolume.h"
#include "G4AffineTransform.h"

#include "G4LogicalVolume_inline.c"
#include "G4VSolid.h"
#include "G4VSolid_inline.c"

#include <vector>
#include <algorithm>
#include <cstdlib>

namespace 
{

typedef std::vector<G4SmartVoxelProxy*> G4ProxyVector;
typedef std::vector<G4SmartVoxelNode*> G4NodeVector;
typedef std::vector<G4int> G4VolumeNosVector;
typedef std::vector<G4double> G4VolumeExtentVector;

}

extern "C"
{

static int minVoxelVolumesLevel1 = K_MIN_VOXEL_VOLUMES_LEVEL_1;
static int minVoxelVolumesLevel2 = K_MIN_VOXEL_VOLUMES_LEVEL_2;
static int minVoxelVolumesLevel3 = K_MIN_VOXEL_VOLUMES_LEVEL_3;

void G4VoxelHeader_SetMinVoxelLimits( int lev1, int lev2, int lev3 )
{
	minVoxelVolumesLevel1 = lev1;
	minVoxelVolumesLevel2 = lev2;
	minVoxelVolumesLevel3 = lev3;
}

int G4VoxelHeader_GetMinVoxelVolumesLevel1()
{
	return minVoxelVolumesLevel1;
}

int G4VoxelHeader_GetMinVoxelVolumesLevel2()
{
	return minVoxelVolumesLevel2;
}

int G4VoxelHeader_GetMinVoxelVolumesLevel3()
{
	return minVoxelVolumesLevel3;
}
	
// FWD

MAYINLINE void
G4VoxelHeader_ctor_protected(
		G4SmartVoxelHeader *This,
		G4LogicalVolume* pVolume,
		G4VoxelLimits pLimits,
		const G4VolumeNosVector* pCandidates,
		G4int pSlice,
		G4double smartless);
	

MAYINLINE void
G4VoxelHeader_SetSlices(
	G4SmartVoxelHeader *This,
	const G4ProxyVector *slices )
{
	This->fNumSlices = slices->size();
	
	std::size_t len = sizeof(G4SmartVoxelProxy*)*This->fNumSlices;
	
	This->fslices = (G4SmartVoxelProxy**)std::realloc( This->fslices, len );
	std::memcpy( This->fslices, &((*slices)[0]), len );
}

// ***************************************************************************
// Equality operator: returns true if contents are equivalent.
// Implies a deep search through contained nodes/header.
// Compares headers' axes,sizes,extents. Returns false if different.
// For each contained proxy, determines whether node/header, compares and
// returns if different. Compares and returns if proxied nodes/headers
// are different.
// ***************************************************************************
//
MAYINLINE
G4bool G4VoxelHeader_operator_equal(
	const G4SmartVoxelHeader *This,
	G4SmartVoxelHeader* pHead)
{
  if ( (G4VoxelHeader_GetAxis(This)      == G4VoxelHeader_GetAxis(pHead))
    && (G4VoxelHeader_GetNoSlices(This)  == G4VoxelHeader_GetNoSlices(This))
    && (G4VoxelHeader_GetMinExtent(This) == G4VoxelHeader_GetMinExtent(This))
    && (G4VoxelHeader_GetMaxExtent(This) == G4VoxelHeader_GetMaxExtent(This)) )
  {
    G4int node, maxNode;
    G4SmartVoxelProxy *leftProxy, *rightProxy;
    G4SmartVoxelHeader *leftHeader, *rightHeader;
    G4SmartVoxelNode *leftNode, *rightNode;

    maxNode=G4VoxelHeader_GetNoSlices(This);
    for (node=0; node<maxNode; node++)
    {
      leftProxy  = G4VoxelHeader_GetSlice(This,node);
      rightProxy = G4VoxelHeader_GetSlice(pHead,node);
      if (G4VoxelProxy_IsHeader(leftProxy))
      {
        if (G4VoxelProxy_IsNode(rightProxy))
        {
          return false;
        }
        else
        {
          leftHeader  = G4VoxelProxy_GetHeader(leftProxy);
          rightHeader = G4VoxelProxy_GetHeader(rightProxy);
          if (!(G4VoxelHeader_operator_equal(leftHeader,rightHeader)))
          {
            return false;
          }
        }
      }
      else
      {
        if (G4VoxelProxy_IsHeader(rightProxy))
        {
          return false;
        }
        else
        {
          leftNode  = G4VoxelProxy_GetNode(leftProxy);
          rightNode = G4VoxelProxy_GetNode(rightProxy);
          if (!(G4VoxelNode_operator_equal(leftNode,rightNode)))
          {
            return false;
          }
        }
      }
    }
    return true;
  }
  else
  {
    return false;
  }
}

// ***************************************************************************
// Builds the nodes corresponding to slices between the specified limits
// and along the specified axis, using candidate volume no.s in the vector
// pCandidates. If the `daughters' are replicated volumes (ie. the logical
// volume has a single replicated/parameterised volume for a daughter)
// the candidate no.s are interpreted as PARAMETERISED volume no.s & 
// PARAMETERISATIONs are applied to compute transformations & solid
// dimensions appropriately. The volume must be parameterised - ie. has a
// parameterisation object & non-consuming) - in this case.
// 
// Returns pointer to built node "structure" (guaranteed non NULL) consisting
// of G4SmartVoxelNodeProxies referring to G4SmartVoxelNodes.
// ***************************************************************************
//
MAYINLINE
G4ProxyVector* G4VoxelHeader_BuildNodes(
		G4SmartVoxelHeader *This,
		G4LogicalVolume* pVolume,
		G4VoxelLimits pLimits,
		const G4VolumeNosVector* pCandidates,
		EAxis pAxis,
		G4double smartlessUser )
{
  (void)This;
	
  G4double motherMinExtent= kInfinity, motherMaxExtent= -kInfinity,
           targetMinExtent= kInfinity, targetMaxExtent= -kInfinity;
  G4VPhysicalVolume *pDaughter=0;
  G4VSolid *targetSolid;
  G4AffineTransform targetTransform = G4AffineTransform_create_id();
  G4int nCandidates = pCandidates->size();
  G4int nVol, nNode, targetVolNo;
  G4VoxelLimits noLimits; G4VoxelLimits_ctor(&noLimits);
    
  // Compute extent of logical volume's solid along this axis
  // NOTE: results stored locally and not preserved/reused
  //
  G4VSolid* outerSolid = G4LogicalVolume_GetSolid(pVolume);
  assert( outerSolid != NULL );
  
  const G4AffineTransform origin = G4AffineTransform_create_id();
  if( !G4VSolid_CalculateExtent(outerSolid, pAxis, pLimits, origin,
                                   &motherMinExtent, &motherMaxExtent) )
  {
    G4VSolid_CalculateExtent(outerSolid, pAxis, noLimits, origin,
                                &motherMinExtent, &motherMaxExtent);
  }
  G4VolumeExtentVector minExtents(nCandidates,0.);
  G4VolumeExtentVector maxExtents(nCandidates,0.);
  int type;
		//type will contain the 'type' of the solid ( see enum ESolids in everything.h) 
  // Compute extents
  //
  for (nVol=0; nVol<nCandidates; nVol++)
  {
	targetVolNo=(*pCandidates)[nVol];

	pDaughter=G4LogicalVolume_GetDaughter(pVolume,targetVolNo);

	// Setup daughter's transformations
	//
	targetTransform = G4AffineTransform_create_full(
		G4VPhysicalVolume_GetObjectRotationValue(pDaughter),
		G4VPhysicalVolume_GetTranslation(pDaughter));
		
	// Get underlying (and setup) solid
	//
	targetSolid = G4LogicalVolume_GetSolid(
		G4VPhysicalVolume_GetLogicalVolume(pDaughter) );
	

	type = targetSolid->type;
		// type gets the Solids type and is now passed in to voxel node insert.
	assert( targetSolid != NULL );

	// Calculate extents
	//
	if(!G4VSolid_CalculateExtent(targetSolid, pAxis, pLimits, targetTransform,
									 &targetMinExtent, &targetMaxExtent))
	{
	  G4VSolid_CalculateExtent(targetSolid, pAxis, noLimits, targetTransform,
								   &targetMinExtent, &targetMaxExtent);
	}
	minExtents[nVol] = targetMinExtent;
	maxExtents[nVol] = targetMaxExtent;

	// Check not entirely outside mother when processing toplevel nodes
	//
	if( (!G4VoxelLimits_IsLimited(&pLimits)) && ((targetMaxExtent<=motherMinExtent) ||(targetMinExtent>=motherMaxExtent)) )
	{
		return NULL;
	}

  }

  // Extents of all daughters known

  // Calculate minimum slice width, only including volumes inside the limits
  //
  G4double minWidth = kInfinity;
  G4double currentWidth;
  for (nVol=0; nVol<nCandidates; nVol++)
  {
    // currentWidth should -always- be a positive value. Inaccurate computed extent
    // from the solid or situations of malformed geometries (overlaps) may lead to
    // negative values and therefore unpredictable crashes !
    //
    currentWidth = std::abs(maxExtents[nVol]-minExtents[nVol]);
    if ( (currentWidth<minWidth)
      && (maxExtents[nVol]>=G4VoxelLimits_GetMinExtent(&pLimits,pAxis))
      && (minExtents[nVol]<=G4VoxelLimits_GetMaxExtent(&pLimits,pAxis)) )
    {
      minWidth = currentWidth;
    }
  }

  // No. of Nodes formula - nearest integer to
  // mother width/half min daughter width +1
  //
  G4double noNodesExactD = ((motherMaxExtent-motherMinExtent)*2.0/minWidth)+1.0;

  // Compare with "smartless quality", i.e. the average number of slices
  // used per contained volume.
  //
  G4double smartlessComputed = noNodesExactD / nCandidates;
  
  G4double smartless = (smartlessComputed <= smartlessUser)
                       ? smartlessComputed : smartlessUser;
  
  G4double noNodesSmart = smartless*nCandidates;
  G4int    noNodesExactI = G4int(noNodesSmart);
  G4int    noNodes = ((noNodesSmart-noNodesExactI)>=0.5)
                     ? noNodesExactI+1 : noNodesExactI;
  if( noNodes == 0 ) { noNodes=1; }
  
  const G4int kMaxVoxelNodes = K_MAX_VOXEL_NODES;

  if (noNodes > kMaxVoxelNodes)
  {
    noNodes=kMaxVoxelNodes; // TODO?
  }
  G4double nodeWidth = (motherMaxExtent-motherMinExtent)/noNodes;

// Create G4VoxelNodes. Will Add proxies before setting fslices
//
  G4NodeVector* nodeList = new G4NodeVector();
  nodeList->reserve(noNodes);
  
  for (nNode=0; nNode<noNodes; nNode++)
  {
    G4SmartVoxelNode *pNode;
    pNode = new G4SmartVoxelNode;
    G4VoxelNode_ctor(pNode,nNode);
    
    nodeList->push_back(pNode);
  }

  // All nodes created (empty)

  // Fill nodes: Step through extent lists
  //
  for (nVol=0; nVol<nCandidates; nVol++)
  {
    G4int nodeNo, minContainingNode, maxContainingNode;
    minContainingNode = G4int((minExtents[nVol]-motherMinExtent)/nodeWidth);
    maxContainingNode = G4int((maxExtents[nVol]-motherMinExtent)/nodeWidth);

    // Only add nodes that are inside the limits of the axis
    //
    if ( (maxContainingNode>=0) && (minContainingNode<noNodes) )
    {
      // If max extent is on max boundary => maxContainingNode=noNodes:
      // should be one less as nodeList has noNodes entries
      //
      if (maxContainingNode>=noNodes)
      {
        maxContainingNode = noNodes-1;
      }
      //
      // Protection against protruding volumes
      //
      if (minContainingNode<0)
      {
        minContainingNode=0;
      }
      for (nodeNo=minContainingNode; nodeNo<=maxContainingNode; nodeNo++)
      {
		 
#ifdef NEW_NAVIGATION
		  G4VoxelNode_Insert( (*nodeList)[nodeNo], (*pCandidates)[nVol], type );
			// Two definitions for Voxelnode_Insert with and without new navigation
			// Perhaps not the cleanest method but convenient for now.
#else
		  G4VoxelNode_Insert( (*nodeList)[nodeNo], (*pCandidates)[nVol] );
#endif
      }
    }
  }

  // All nodes filled

  // Create proxy List : caller has deletion responsibility
  // (but we must delete nodeList *itself* - not the contents)
  //
  G4ProxyVector* proxyList = new G4ProxyVector();
  proxyList->reserve(noNodes);
  
  //
  // Fill proxy List
  //
  for (nNode=0; nNode<noNodes; nNode++)
  {

    G4SmartVoxelProxy* pProxyNode = new G4SmartVoxelProxy;
    G4VoxelProxy_ctor_node( pProxyNode, ((*nodeList)[nNode]) );
    
    proxyList->push_back(pProxyNode);
  }
  delete nodeList;
  return proxyList;
}


// ***************************************************************************
// Collects common nodes at our level, deleting all but one to save
// memory, and adjusting stored slice pointers appropriately.
//
// Preconditions:
// o the slices have not previously be "collected"
// o all of the slices are nodes.
// ***************************************************************************
//
MAYINLINE
void G4VoxelHeader_CollectEquivalentNodes( G4SmartVoxelHeader *This )
{
  G4int sliceNo, maxNo, equivNo;
  G4int maxNode=G4VoxelHeader_GetNoSlices(This);
  G4SmartVoxelNode *equivNode;
  G4SmartVoxelProxy *equivProxy;

  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
  {
    equivProxy=This->fslices[sliceNo];

    // Assumption (see preconditions): all slices are nodes
    //
    equivNode = G4VoxelProxy_GetNode(equivProxy);
    maxNo = G4VoxelNode_GetMaxEquivalentSliceNo(equivNode);
    if (maxNo != sliceNo)
    {
      // Do collection between sliceNo and maxNo inclusive
      //
      for (equivNo=sliceNo+1; equivNo<=maxNo; equivNo++)
      {
		G4SmartVoxelNode *dyingNode = G4VoxelProxy_GetNode(This->fslices[equivNo]);
		G4VoxelNode_dtor( dyingNode );
        delete dyingNode;
        delete This->fslices[equivNo];
        This->fslices[equivNo] = equivProxy;
      }
      sliceNo = maxNo;
    }
  }
}

// ***************************************************************************
// Destructor:
// deletes all proxies and underlying objects.
// ***************************************************************************
//
MAYINLINE void
G4VoxelHeader_dtor(G4SmartVoxelHeader *This)
{
  // Manually destroy underlying nodes/headers
  // Delete collected headers and nodes once only
  //
  G4int node, proxy, maxNode=G4VoxelHeader_GetNoSlices(This);
  G4SmartVoxelProxy *lastProxy=0;
  G4SmartVoxelNode *dyingNode, *lastNode=0;
  G4SmartVoxelHeader *dyingHeader, *lastHeader=0;

  for (node=0; node<maxNode; node++)
  {
    if (G4VoxelProxy_IsHeader(This->fslices[node]))
    {
      dyingHeader = G4VoxelProxy_GetHeader(This->fslices[node]);
      if (lastHeader!=dyingHeader)
      {
        lastHeader = dyingHeader;
        lastNode = 0;
        G4VoxelHeader_dtor( dyingHeader );
        delete dyingHeader;
      }
    }
    else
    {
      dyingNode = G4VoxelProxy_GetNode(This->fslices[node]);
      if (dyingNode!=lastNode)
      {
        lastNode=dyingNode;
        lastHeader=0;
        
        G4VoxelNode_dtor( dyingNode );
        delete dyingNode;
      }
    }
  }
  // Delete proxies
  //
  for (proxy=0; proxy<maxNode; proxy++)
  {
    if (This->fslices[proxy]!=lastProxy)
    {
      lastProxy = This->fslices[proxy];
      delete lastProxy;
    }
  }
  
  std::free(This->fslices);
}

// ***************************************************************************
// Collects common headers at our level, deleting all but one to save
// memory, and adjusting stored slice pointers appropriately.
// 
// Preconditions:
// o if a header forms part of a range of equivalent slices
//   (ie. GetMaxEquivalentSliceNo()>GetMinEquivalentSliceNo()),
//   it is assumed that all slices in the range are headers.
// o this will be true if a constant Expression is used to evaluate
//   when to refine nodes.
// ***************************************************************************
//
MAYINLINE
void G4VoxelHeader_CollectEquivalentHeaders( G4SmartVoxelHeader *This )
{
  G4int sliceNo, maxNo, equivNo;
  G4int maxNode = G4VoxelHeader_GetNoSlices(This);
  G4SmartVoxelHeader *equivHeader, *sampleHeader;
  G4SmartVoxelProxy *equivProxy;

  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
  {
    equivProxy = This->fslices[sliceNo];
    if (G4VoxelProxy_IsHeader(equivProxy))
    {
      equivHeader = G4VoxelProxy_GetHeader(equivProxy);
      maxNo = G4VoxelHeader_GetMaxEquivalentSliceNo(equivHeader);
      if (maxNo != sliceNo)
      {
        // Attempt collection between sliceNo and maxNo inclusive:
        // look for common headers. All slices between sliceNo and maxNo
        // are guaranteed to be headers but may not have equal contents
        //

        for (equivNo=sliceNo+1; equivNo<=maxNo; equivNo++)
        {
          sampleHeader = G4VoxelProxy_GetHeader(This->fslices[equivNo]);
          if ( G4VoxelHeader_operator_equal(sampleHeader,equivHeader) )
          {

            // Delete sampleHeader + proxy and replace with equivHeader/Proxy
            //
            G4VoxelHeader_dtor( sampleHeader );
            delete sampleHeader;
            delete This->fslices[equivNo];
            This->fslices[equivNo] = equivProxy;
          }
          else
          {
            // Not equal. Set this header to be
            // the current header for comparisons
            //
            equivProxy  = This->fslices[equivNo];
            equivHeader = G4VoxelProxy_GetHeader(equivProxy);
          }

        }
        // Skip past examined slices
        //
        sliceNo = maxNo;
      }
    }
  }
}


// ***************************************************************************
// Calculate a "quality value" for the specified vector of voxels.
// The value returned should be >0 and such that the smaller the number
// the higher the quality of the slice.
//
// Preconditions: pSlice must consist of G4SmartVoxelNodeProxies only
// Process:
// o Examine each node in turn, summing:
//      no. of non-empty nodes
//      no. of volumes in each node
// o Calculate Quality=sigma(volumes in nod)/(no. of non-empty nodes)
//      if all nodes empty, return kInfinity
// o Call G4Exception on finding a G4SmartVoxelHeaderProxy
// ***************************************************************************
//
MAYINLINE
G4double G4VoxelHeader_CalculateQuality(
		G4SmartVoxelHeader *This,
		G4ProxyVector *pSlice)
{
  (void)This;
	
  G4double quality;
  G4int nNodes = pSlice->size();
  G4int noContained, maxContained=0, sumContained=0, sumNonEmptyNodes=0;
  G4SmartVoxelNode *node;

  for (G4int i=0; i<nNodes; i++)
  {
    if (G4VoxelProxy_IsNode((*pSlice)[i]))
    {
      // Definitely a node. Add info to running totals
      //
      node = G4VoxelProxy_GetNode((*pSlice)[i]);
      noContained = G4VoxelNode_GetNoContained(node);
      if (noContained)
      {
        sumNonEmptyNodes++;
        sumContained += noContained;
        //
        // Calc maxContained for statistics
        //
        if (noContained>maxContained)
        {
          maxContained = noContained;
        }
      }
    }
  }

  // Calculate quality with protection against no non-empty nodes
  //
  if (sumNonEmptyNodes)
  {
    quality = sumContained/sumNonEmptyNodes;
  }
  else
  {
    quality = kInfinity;
  }

  return quality;
}

// ***************************************************************************
// Calculates and stores the minimum and maximum equivalent neighbour
// values for all slices at our level.
//
// Precondition: all slices are nodes.
// For each potential start of a group of equivalent nodes:
// o searches forwards in fslices to find group end
// o loops from start to end setting start and end slices.
// ***************************************************************************
//
MAYINLINE
void G4VoxelHeader_BuildEquivalentSliceNos( G4SmartVoxelHeader *This )
{
  G4int sliceNo, minNo, maxNo, equivNo;
  G4int maxNode = G4VoxelHeader_GetNoSlices(This);
  G4SmartVoxelNode *startNode, *sampleNode;
  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
  {
    minNo = sliceNo;

    // Get first node (see preconditions - will throw exception if a header)
    //
    startNode = G4VoxelProxy_GetNode(This->fslices[minNo]);
    assert(startNode != NULL);

    // Find max equivalent
    //
    for (equivNo=minNo+1; equivNo<maxNode; equivNo++)
    {
      sampleNode = G4VoxelProxy_GetNode(This->fslices[equivNo]);
      assert( sampleNode != NULL );
      if (!(G4VoxelNode_operator_equal(startNode,sampleNode))) { break; }
    }
    maxNo = equivNo-1;
    if (maxNo != minNo)
    {
      // Set min and max nos
      //
      for (equivNo=minNo; equivNo<=maxNo; equivNo++)
      {
        sampleNode = G4VoxelProxy_GetNode(This->fslices[equivNo]);
        G4VoxelNode_SetMinEquivalentSliceNo(sampleNode,minNo);
        G4VoxelNode_SetMaxEquivalentSliceNo(sampleNode,maxNo);
      }
      // Advance outer loop to end of equivalent group
      //
      sliceNo = maxNo;
    }
  }
}

// ***************************************************************************
// Examined each contained node, refines (creates a replacement additional
// dimension of voxels) when there is more than one voxel in the slice.
// Does not refine further if already limited in two dimensions (=> this
// is the third level of limits)
//
// Preconditions: slices (nodes) have been built.
// ***************************************************************************
//
MAYINLINE
void G4VoxelHeader_RefineNodes(
		G4SmartVoxelHeader *This,
		G4LogicalVolume* pVolume,
		G4VoxelLimits pLimits,
		G4double smartless)
{
  G4int refinedDepth=0, minVolumes;
  G4int maxNode = G4VoxelHeader_GetNoSlices(This);

  if (G4VoxelLimits_IsXLimited(&pLimits))
  {
    refinedDepth++;
  }
  if (G4VoxelLimits_IsYLimited(&pLimits)) 
  {
    refinedDepth++;
  }
  if (G4VoxelLimits_IsZLimited(&pLimits))
  {
    refinedDepth++;
  }
  
  //const G4int kMinVoxelVolumesLevel2 = K_MIN_VOXEL_VOLUMES_LEVEL_2;
  //const G4int kMinVoxelVolumesLevel3 = K_MIN_VOXEL_VOLUMES_LEVEL_3;

  // Calculate minimum number of volumes necessary to refine
  //
  switch (refinedDepth)
  {
    case 0:
      minVolumes=minVoxelVolumesLevel2;
      break;
    case 1:
      minVolumes=minVoxelVolumesLevel3;
      break;
    default:
      minVolumes=10000;   // catch refinedDepth=3 and errors
      break;
  }

  if (refinedDepth<2)
  {
    G4int targetNo, noContainedDaughters, minNo, maxNo, replaceNo, i;
    G4double sliceWidth = (This->fmaxExtent-This->fminExtent)/maxNode;
    G4VoxelLimits newLimits; G4VoxelLimits_ctor(&newLimits);
    G4SmartVoxelNode* targetNode;
    G4SmartVoxelProxy* targetNodeProxy;
    G4SmartVoxelHeader* replaceHeader;
    G4SmartVoxelProxy* replaceHeaderProxy;
    G4VolumeNosVector* targetList;
    G4SmartVoxelProxy* lastProxy;
      
    for (targetNo=0; targetNo<maxNode; targetNo++)
    {
      // Assume all slices are nodes (see preconditions)
      //
      targetNodeProxy = This->fslices[targetNo];
      targetNode = G4VoxelProxy_GetNode(targetNodeProxy);

      if (G4VoxelNode_GetNoContained(targetNode) >= minVolumes)
      {
        noContainedDaughters = G4VoxelNode_GetNoContained(targetNode);
        targetList = new G4VolumeNosVector();
        targetList->reserve(noContainedDaughters);
        
        for (i=0; i<noContainedDaughters; i++)
        {
          targetList->push_back(G4VoxelNode_GetVolume(targetNode,i));
        }
        minNo = G4VoxelNode_GetMinEquivalentSliceNo(targetNode);
        maxNo = G4VoxelNode_GetMaxEquivalentSliceNo(targetNode);

        // Delete node proxies at start of collected sets of nodes/headers
        //
        lastProxy=0;
        for (replaceNo=minNo; replaceNo<=maxNo; replaceNo++)
        {
          if (lastProxy != This->fslices[replaceNo])
          {
            lastProxy=This->fslices[replaceNo];
            delete lastProxy;
          }
        }
        // Delete node to be replaced
        //
        G4VoxelNode_dtor( targetNode );
        delete targetNode;

        // Create new headers + proxies and replace in fslices
        //
        newLimits = pLimits;
        G4VoxelLimits_AddLimit(&newLimits,This->faxis,This->fminExtent+sliceWidth*minNo,
                           This->fminExtent+sliceWidth*(maxNo+1));
                           
        replaceHeader = new G4SmartVoxelHeader;
        G4VoxelHeader_ctor_protected(
			replaceHeader,pVolume,newLimits,targetList,replaceNo,smartless);
			
        G4VoxelHeader_SetMinEquivalentSliceNo(replaceHeader,minNo);
        G4VoxelHeader_SetMaxEquivalentSliceNo(replaceHeader,maxNo);
        replaceHeaderProxy = new G4SmartVoxelProxy;
        G4VoxelProxy_ctor_header(replaceHeaderProxy,replaceHeader);
        
        for (replaceNo=minNo; replaceNo<=maxNo; replaceNo++)
        {
          This->fslices[replaceNo] = replaceHeaderProxy;
        }
        // Finished replacing current `equivalent' group
        //
        delete targetList;
        targetNo=maxNo;
      }
    }
  }
}

// ***************************************************************************
// Builds and refines voxels between specified limits, considering only
// the physical volumes numbered `pCandidates'.
// o Chooses axis
// o Determines min and max extents (of mother solid) within limits.
// ***************************************************************************
//
MAYINLINE void
G4VoxelHeader_BuildVoxelsWithinLimits(
		G4SmartVoxelHeader *This, 
		G4LogicalVolume* pVolume,
		G4VoxelLimits pLimits,
		const G4VolumeNosVector* pCandidates,
		G4double smartless)
{
  // Choose best axis for slicing by:
  // 1. Trying all unlimited cartesian axes
  // 2. Select axis which gives greatest no slices

  G4ProxyVector *pGoodSlices = NULL, *pTestSlices, *tmpSlices;
  G4double goodSliceScore=kInfinity, testSliceScore;
  EAxis goodSliceAxis = kXAxis;
  EAxis testAxis      = kXAxis;
  G4int node, maxNode, iaxis;
  G4VoxelLimits noLimits; G4VoxelLimits_ctor(&noLimits);

  // Try all non-limited cartesian axes
  //
  for (iaxis=0; iaxis<3; iaxis++)
  {
    switch(iaxis)
    {
      case 0:
        testAxis = kXAxis;
        break;
      case 1:
        testAxis = kYAxis;
        break;
      case 2:
        testAxis = kZAxis;
        break;
    }
    if (!G4VoxelLimits_IsLimited_axis(&pLimits,testAxis))
    {
      pTestSlices = G4VoxelHeader_BuildNodes(This,pVolume,pLimits,pCandidates,testAxis,smartless);
      
      if (pTestSlices == NULL) return;
      
      testSliceScore = G4VoxelHeader_CalculateQuality(This,pTestSlices);
      if ( (!pGoodSlices) || (testSliceScore<goodSliceScore) )
      {
        goodSliceAxis  = testAxis;
        goodSliceScore = testSliceScore;
        tmpSlices      = pGoodSlices;
        pGoodSlices    = pTestSlices;
        pTestSlices    = tmpSlices;
      }
      if (pTestSlices)
      {
        // Destroy pTestSlices and all its contents
        //
        maxNode=pTestSlices->size();
        for (node=0; node<maxNode; node++)
        {
			G4SmartVoxelNode *dyingNode = G4VoxelProxy_GetNode((*pTestSlices)[node]);
			G4VoxelNode_dtor( dyingNode );
			delete dyingNode;
        }
        G4SmartVoxelProxy* tmpProx;
        while (pTestSlices->size()>0)
        {
          tmpProx = pTestSlices->back();
          pTestSlices->pop_back();
          for (G4ProxyVector::iterator i=pTestSlices->begin();
                                       i!=pTestSlices->end(); i++)
          {
            if (*i==tmpProx)
            {
              pTestSlices->erase(i); i--;
            }
          }
          if ( tmpProx ) { delete tmpProx; }
        } 
        delete pTestSlices;
      }
    }
  }
  // Check for error case.. when limits already 3d,
  // so cannot select a new axis
  //
  assert( pGoodSlices != NULL );

  // 
  // We have selected pGoodSlices, with a score testSliceScore
  //

  // Store chosen axis, slice ptr
  //
  
  // TODO!!
  G4VoxelHeader_SetSlices(This,pGoodSlices);
  
  delete pGoodSlices;   // Destroy slices vector, but not contained
                        // proxies or nodes
  This->faxis=goodSliceAxis;

  // Calculate and set min and max extents given our axis
  //
  
  G4VSolid* outerSolid = G4LogicalVolume_GetSolid(pVolume);
  const G4AffineTransform origin = G4AffineTransform_create_id();
  
  assert( outerSolid != NULL );
  
  if(!G4VSolid_CalculateExtent(outerSolid,This->faxis,pLimits,origin,&(This->fminExtent),&(This->fmaxExtent)))
  {
    G4VSolid_CalculateExtent(outerSolid,This->faxis,noLimits,origin,&(This->fminExtent),&(This->fmaxExtent));
  }

  // Calculate equivalent nos
  //
  G4VoxelHeader_BuildEquivalentSliceNos(This);
  G4VoxelHeader_CollectEquivalentNodes(This);     // Collect common nodes
  G4VoxelHeader_RefineNodes(This,pVolume,pLimits,smartless); // Refine nodes creating headers

  // No common headers can exist because collapsed by construction
}

// ***************************************************************************
// Builds voxels for daughters specified volume, in NON-REPLICATED case
// o Create List of target volume nos (all daughters; 0->noDaughters-1)
// o BuildWithinLimits does Build & also determines mother dimensions.
// ***************************************************************************
//
void G4VoxelHeader_BuildVoxels(G4SmartVoxelHeader *This, G4LogicalVolume* pVolume, G4double smartless)
{
	G4VoxelLimits limits;   // Create `unlimited' limits object
	G4VoxelLimits_ctor(&limits);
	G4int nDaughters = G4LogicalVolume_GetNoDaughters(pVolume);

	G4VolumeNosVector targetList;
	targetList.reserve(nDaughters);

	for (G4int i=0; i<nDaughters; i++)
	{
		targetList.push_back(i);
	}
	G4VoxelHeader_BuildVoxelsWithinLimits(This, pVolume, limits, &targetList, smartless);
}


// ***************************************************************************
// Returns true if all slices have equal contents.
// Preconditions: all equal slices have been collected.
// Procedure:
// o checks all slice proxy pointers are equal
// o returns true if only one slice or all slice proxies pointers equal.
// ***************************************************************************
//
MAYINLINE
G4bool G4VoxelHeader_AllSlicesEqual( const G4SmartVoxelHeader *This )
{
  G4int noSlices = G4VoxelHeader_GetNoSlices(This);
  G4SmartVoxelProxy* refProxy;

  if (noSlices>1)
  {
    refProxy=This->fslices[0];
    for (G4int i=1; i<noSlices; i++)
    {
      if (refProxy!=This->fslices[i])
      {
        return false;
      }
    }
  }
  return true;
}


// ***************************************************************************
// Constructor for topmost header, to begin voxel construction at a
// given logical volume.
// Constructs target List of volumes, calls "Build and refine" constructor.
// Assumes all daughters represent single volumes (ie. no divisions
// or parametric)
// ***************************************************************************
//
MAYINLINE void
G4VoxelHeader_ctor(
		G4SmartVoxelHeader *This,
		G4LogicalVolume* pVolume,
		G4int pSlice,
		G4double smartless )
{
	assert( pVolume != NULL );
		
	This->fNumSlices = 0;
	This->fslices = NULL;
	
	This->faxis = kUndefined;
	This->fminEquivalent = pSlice;
	This->fmaxEquivalent = pSlice;
	This->fparamAxis = kUndefined;

	if ( G4LogicalVolume_GetNoDaughters(pVolume) >= minVoxelVolumesLevel1 )
	{
		G4VoxelHeader_BuildVoxels(This,pVolume,smartless);
	}
}

// ***************************************************************************
// Prs and refines voxels between specified limits, considering only
// the physical volumes numbered `pCandidates'. `pSlice' is used to set max
// and min equivalent slice nos for the header - they apply to the level
// of the header, not its nodes.
// ***************************************************************************
//
MAYINLINE void
G4VoxelHeader_ctor_protected(
		G4SmartVoxelHeader *This,
		G4LogicalVolume* pVolume,
		G4VoxelLimits pLimits,
		const G4VolumeNosVector* pCandidates,
		G4int pSlice,
		G4double smartless)
{
	This->fNumSlices = 0;
	This->fslices = NULL;
	
	This->faxis = kUndefined;
	This->fminEquivalent = pSlice;
	This->fmaxEquivalent = pSlice;
	This->fparamAxis = kUndefined;

	G4VoxelHeader_BuildVoxelsWithinLimits(This,pVolume,pLimits,pCandidates, smartless);
}


}

#endif
