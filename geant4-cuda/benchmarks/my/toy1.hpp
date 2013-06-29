
/** A very simple test geometry with a sphere inside a box */

#ifndef TOYGEOMETRY1_HPP
#define TOYGEOMETRY1_HPP

#include "../../geometry_common.hpp"

class ToyGeometry1 : public BasicGeometry
{
public:
	ToyGeometry1( int maxsz = 1024 )
	 : BasicGeometry( maxsz ) {}
	
	double getScale() const
	{
		return 1.0;
	}

	CameraParameters getCamera() const
	{
		CameraParameters p;
		p.heading = 0;
		p.pitch = 0;
		p.dist = 0.1;
		p.yfov = 80;
		return p;
	}

private:

	void create_impl()
	{
		printf ("Entered in toy1.hpp\n");
		const G4double mybox_x = 2.;
		const G4double mybox_y = 5.;
		const G4double mybox_z = 2.1;
		const G4double myorb_rad = 0.5;
		G4RotationMatrix idRot = G4RotationMatrix_create_id();
		G4ThreeVector zeroTrans = G4ThreeVector_create( 0,0,0 );
		G4ThreeVector orbTrans = G4ThreeVector_create( 0.7, 2.7, 0.7 );
	
		// Allocate
		
		// World volume must be first
		G4VPhysicalVolume *phys_box = reserveThing<G4VPhysicalVolume>();
		G4VPhysicalVolume *phys_orb = reserveThing<G4VPhysicalVolume>();
		G4LogicalVolume *logical_box = reserveThing<G4LogicalVolume>();
		G4LogicalVolume *logical_orb = reserveThing<G4LogicalVolume>();
				
		G4Box *mybox = reserveThing<G4Box>();
		G4Orb *myorb = reserveThing<G4Orb>();
				
		StubMaterial *box_material = reserveThing<StubMaterial>();
		StubMaterial *orb_material = reserveThing<StubMaterial>();
		
		// Construct
		box_material->property = 500.0;
		orb_material->property = 2000.0;
		
		G4Box_ctor(mybox, mybox_x, mybox_y, mybox_z);
		G4Orb_ctor(myorb, myorb_rad);
		
		G4LogicalVolume_ctor(logical_box, (G4VSolid*)mybox, box_material);
		G4LogicalVolume_ctor(logical_orb, (G4VSolid*)myorb, orb_material);
		
		G4VPhysicalVolume_ctor(phys_box, idRot, zeroTrans, logical_box);
		G4VPhysicalVolume_ctor(phys_orb, idRot, orbTrans, logical_orb);
		
		addLogicalVolumeDaughter( logical_box, phys_orb );
		G4VPhysicalVolume_SetMotherLogical(phys_orb, logical_box);
		
		// Register pointers for relocation
		addLogicalVolumePointers( logical_orb );
		addLogicalVolumePointers( logical_box );
		addPhysicalVolumePointers( phys_orb );
		addPhysicalVolumePointers( phys_box );
	}
};

#endif
