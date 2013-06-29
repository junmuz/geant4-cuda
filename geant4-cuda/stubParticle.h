
/** Dummy particle */

#ifndef __STUBPARTICLE_HH__
#define __STUBPARTICLE_HH__

#include "G4ThreeVector.h"

typedef struct //__attribute__((packed))
{
	G4ThreeVector pos, dir;
} 
StubParticle;

typedef struct
{
	G4ThreeVector pos, dir;
	G4double t;
}
ParticleWithLifetime;

#ifdef PHYSICS
typedef ParticleWithLifetime Particle;
#else
typedef StubParticle Particle;
#endif

//#define G4Particle StubParticle

#endif
