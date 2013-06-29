
/**
 * G4VSolid host side implementation
 * based on G4VSolid.cc of Geant 4.9.3
 */

#ifndef G4VSOLID_C
#define G4VSOLID_C

#ifdef HOST_CODE

#include "G4VSolid.h"

#include "G4VoxelLimits.c"

////////////////////////////////////////////////////////////////////////////
//
// pVoxelLimits must be only limited along one axis, and either the maximum
// along the axis must be +kInfinity, or the minimum -kInfinity

MAYINLINE void
G4VSolid_ClipPolygonToSimpleLimits( G4ThreeVectorList& pPolygon,
                                     G4ThreeVectorList& outputPolygon,
                               const G4VoxelLimits *pVoxelLimit )
{
  G4int i;
  G4int noVertices=pPolygon.size();
  G4ThreeVector vEnd,vStart;

  for (i = 0 ; i < noVertices ; i++ )
  {
    vStart = pPolygon[i];
    // G4cout << "i = " << i << G4endl;
    if ( i == noVertices-1 )    vEnd = pPolygon[0];
    else                        vEnd = pPolygon[i+1];

    if ( G4VoxelLimits_Inside(pVoxelLimit,vStart) )
    {
      if (G4VoxelLimits_Inside(pVoxelLimit,vEnd))
      {
        // vStart and vEnd inside -> output end point
        //
        outputPolygon.push_back(vEnd);
      }
      else
      {
        // vStart inside, vEnd outside -> output crossing point
        //
        // G4cout << "vStart inside, vEnd outside" << G4endl;
        G4VoxelLimits_ClipToLimits(pVoxelLimit,&vStart,&vEnd);
        outputPolygon.push_back(vEnd);
      }    
    }
    else
    {
      if (G4VoxelLimits_Inside(pVoxelLimit,vEnd))
      {
        // vStart outside, vEnd inside -> output inside section
        //
        // G4cout << "vStart outside, vEnd inside" << G4endl;
        G4VoxelLimits_ClipToLimits(pVoxelLimit,&vStart,&vEnd);
        outputPolygon.push_back(vStart);
        outputPolygon.push_back(vEnd);  
      }
      else  // Both point outside -> no output
      {
        // outputPolygon.push_back(vStart);
        // outputPolygon.push_back(vEnd);  
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Clip the convex polygon described by the vertices at
// pSectionIndex ->pSectionIndex+3 within pVertices to the limits pVoxelLimit
//
// Set pMin to the smallest
//
// Calculate the extent of the polygon along pAxis, when clipped to the
// limits pVoxelLimit. If the polygon exists after clippin, set pMin to
// the polygon's minimum extent along the axis if <pMin, and set pMax to
// the polygon's maximum extent along the axis if >pMax.
//
// The polygon is described by a set of vectors, where each vector represents
// a vertex, so that the polygon is described by the vertex sequence:
//   0th->1st 1st->2nd 2nd->... nth->0th
//
// Modifications to the polygon are made
//
// NOTE: Execessive copying during clipping

MAYINLINE void G4VSolid_ClipPolygon( G4ThreeVectorList& pPolygon,
                           const G4VoxelLimits *pVoxelLimit,
                           const EAxis )
{
  G4ThreeVectorList outputPolygon;

  if ( G4VoxelLimits_IsLimited(pVoxelLimit) )
  {
    if (G4VoxelLimits_IsXLimited(pVoxelLimit) ) // && pAxis != kXAxis)
    {
      G4VoxelLimits simpleLimit1; G4VoxelLimits_ctor(&simpleLimit1);
      G4VoxelLimits_AddLimit(&simpleLimit1, kXAxis,G4VoxelLimits_GetMinXExtent(pVoxelLimit),kInfinity);
      //  G4cout<<"MinXExtent()"<<G4endl;
      G4VSolid_ClipPolygonToSimpleLimits(pPolygon,outputPolygon,&simpleLimit1);
   
      pPolygon.clear();

      if ( !outputPolygon.size() )  return;

      G4VoxelLimits simpleLimit2; G4VoxelLimits_ctor(&simpleLimit2);
      //  G4cout<<"MaxXExtent()"<<G4endl;
      G4VoxelLimits_AddLimit(&simpleLimit2, kXAxis,-kInfinity,G4VoxelLimits_GetMaxXExtent(pVoxelLimit));
      G4VSolid_ClipPolygonToSimpleLimits(outputPolygon,pPolygon,&simpleLimit2);

      if ( !pPolygon.size() )       return;
      else                          outputPolygon.clear();
    }
    if ( G4VoxelLimits_IsYLimited(pVoxelLimit) ) // && pAxis != kYAxis)
    {
      G4VoxelLimits simpleLimit1; G4VoxelLimits_ctor(&simpleLimit1);
      G4VoxelLimits_AddLimit(&simpleLimit1, kYAxis,G4VoxelLimits_GetMinYExtent(pVoxelLimit),kInfinity);
      G4VSolid_ClipPolygonToSimpleLimits(pPolygon,outputPolygon,&simpleLimit1);

      // Must always clear pPolygon - for clip to simpleLimit2 and in case of
      // early exit

      pPolygon.clear();

      if ( !outputPolygon.size() )  return;

      G4VoxelLimits simpleLimit2; G4VoxelLimits_ctor(&simpleLimit2);
      G4VoxelLimits_AddLimit(&simpleLimit2, kYAxis,-kInfinity,G4VoxelLimits_GetMaxYExtent(pVoxelLimit));
      G4VSolid_ClipPolygonToSimpleLimits(outputPolygon,pPolygon,&simpleLimit2);

      if ( !pPolygon.size() )       return;
      else                          outputPolygon.clear();
    }
    if ( G4VoxelLimits_IsZLimited(pVoxelLimit) ) // && pAxis != kZAxis)
    {
      G4VoxelLimits simpleLimit1; G4VoxelLimits_ctor(&simpleLimit1);
      G4VoxelLimits_AddLimit(&simpleLimit1, kZAxis,G4VoxelLimits_GetMinZExtent(pVoxelLimit),kInfinity);
      G4VSolid_ClipPolygonToSimpleLimits(pPolygon,outputPolygon,&simpleLimit1);

      // Must always clear pPolygon - for clip to simpleLimit2 and in case of
      // early exit

      pPolygon.clear();

      if ( !outputPolygon.size() )  return;

      G4VoxelLimits simpleLimit2; G4VoxelLimits_ctor(&simpleLimit2);
      G4VoxelLimits_AddLimit(&simpleLimit2, kZAxis,-kInfinity,G4VoxelLimits_GetMaxZExtent(pVoxelLimit));
      G4VSolid_ClipPolygonToSimpleLimits(outputPolygon,pPolygon,&simpleLimit2);

      // Return after final clip - no cleanup
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// Calculate the maximum and minimum extents of the convex polygon pPolygon
// along the axis pAxis, within the limits pVoxelLimit
//

MAYINLINE void
G4VSolid_CalculateClippedPolygonExtent(G4ThreeVectorList& pPolygon,
                                  const G4VoxelLimits* pVoxelLimit,
                                  const EAxis pAxis, 
                                        G4double& pMin,
                                        G4double& pMax)
{
  G4int noLeft,i;
  G4double component;
  
  G4VSolid_ClipPolygon(pPolygon,pVoxelLimit,pAxis);
  noLeft = pPolygon.size();

  if ( noLeft )
  {
    for (i=0;i<noLeft;i++)
    {
      component = G4ThreeVector_coord(pPolygon[i], pAxis);
 
      if (component < pMin) 
      { 
        pMin = component;      
      }
      if (component > pMax)
      {  
        pMax = component;  
      }    
    }
  }
}

///////////////////////////////////////////////////////////////////////////
// 
// Calculate the maximum and minimum extents of the polygon described
// by the vertices: pSectionIndex->pSectionIndex+1->
//                   pSectionIndex+2->pSectionIndex+3->pSectionIndex
// in the List pVertices
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices
//

MAYINLINE void G4VSolid_ClipCrossSection(       G4ThreeVectorList* pVertices,
                                 const G4int pSectionIndex,
                                 const G4VoxelLimits *pVoxelLimit,
                                 const EAxis pAxis, 
                                       G4double& pMin, G4double& pMax)
{
  G4ThreeVectorList polygon;
  polygon.reserve(4);
  polygon.push_back((*pVertices)[pSectionIndex]);
  polygon.push_back((*pVertices)[pSectionIndex+1]);
  polygon.push_back((*pVertices)[pSectionIndex+2]);
  polygon.push_back((*pVertices)[pSectionIndex+3]);
  G4VSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
}

//////////////////////////////////////////////////////////////////////////////////
//
// Calculate the maximum and minimum extents of the polygons
// joining the CrossSections at pSectionIndex->pSectionIndex+3 and
//                              pSectionIndex+4->pSectionIndex7
//
// in the List pVertices, within the boundaries of the voxel limits pVoxelLimit
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices

MAYINLINE void G4VSolid_ClipBetweenSections(      G4ThreeVectorList* pVertices,
                                   const G4int pSectionIndex,
                                   const G4VoxelLimits* pVoxelLimit,
                                   const EAxis pAxis, 
                                         G4double& pMin, G4double& pMax)
{
  G4ThreeVectorList polygon;
  polygon.reserve(4);
  polygon.push_back((*pVertices)[pSectionIndex]);
  polygon.push_back((*pVertices)[pSectionIndex+4]);
  polygon.push_back((*pVertices)[pSectionIndex+5]);
  polygon.push_back((*pVertices)[pSectionIndex+1]);
  // G4cout<<"ClipBetweenSections: 0-4-5-1"<<G4endl;
  G4VSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+1]);
  polygon.push_back((*pVertices)[pSectionIndex+5]);
  polygon.push_back((*pVertices)[pSectionIndex+6]);
  polygon.push_back((*pVertices)[pSectionIndex+2]);
  // G4cout<<"ClipBetweenSections: 1-5-6-2"<<G4endl;
  G4VSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+2]);
  polygon.push_back((*pVertices)[pSectionIndex+6]);
  polygon.push_back((*pVertices)[pSectionIndex+7]);
  polygon.push_back((*pVertices)[pSectionIndex+3]);
  //  G4cout<<"ClipBetweenSections: 2-6-7-3"<<G4endl;
  G4VSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+3]);
  polygon.push_back((*pVertices)[pSectionIndex+7]);
  polygon.push_back((*pVertices)[pSectionIndex+4]);
  polygon.push_back((*pVertices)[pSectionIndex]);
  //  G4cout<<"ClipBetweenSections: 3-7-4-0"<<G4endl;
  G4VSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
}

#endif // host code

#endif // guard
