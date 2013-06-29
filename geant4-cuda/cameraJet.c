
/** Raytracer camera "particle jet" generation */

#include "geometry.hpp"
#include "stubParticle.h"
#include <math.h>

#define DEG2RAD(x) ((x)/180.0*M_PI)

void createCameraJet( StubParticle *container, int xres, int yres, CameraParameters params )
{
	const G4ThreeVector target = G4ThreeVector_create(
		params.target_x,
		params.target_y,
		params.target_z );
		
	const G4double
		heading = params.heading,
		pitch = params.pitch,
		roll = params.roll,
		dist = params.dist;
		
	const G4double yfov = params.yfov;
	
	const G4double ymul = tan(DEG2RAD(yfov)/2)*2;
	const G4double xmul = ymul / yres * xres;
	
	for ( int y = 0; y<yres; ++y )
	{
		for ( int x = 0; x<xres; ++x )
		{
			G4double dy = (y/(G4double)(yres-1)-0.5)*ymul;
			G4double dx = (x/(G4double)(xres-1)-0.5)*xmul;
			
			StubParticle p;
			p.pos = G4ThreeVector_create(0,-dist,0);
			G4ThreeVector_rotate( &(p.pos), DEG2RAD(heading), DEG2RAD(pitch), DEG2RAD(roll) );
			
			G4ThreeVector_sum_assign( &(p.pos), target );
			
			p.dir = G4ThreeVector_unit(G4ThreeVector_create(dx,1,dy));
				G4ThreeVector_rotate( &(p.dir), DEG2RAD(heading), DEG2RAD(pitch), 0);
				
			
			container[y*xres+x] = p;
		}
	}
}
