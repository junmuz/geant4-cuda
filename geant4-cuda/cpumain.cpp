
/** Main CPU host file */

#include "myutils.h"

#include "everything.h"
#include "stubParticle.h"
#include "G4VPhysicalVolume.h"

#include <vector>
#include <cstdio>

#include "hostcommons.hpp"

// Test relocation on CPU to ensure it works properly
#define RELOCATE

extern "C" {
void cpuexec(
	int problemSz,
	const void *particles,
	G4VPhysicalVolume *geometryRoot,
	void *output,
	G4double phys_step);
}

int main(int argc, char *argv[])
{
#ifndef NOCATCH
  try
  {
#endif
	TestCase testCase( argc, argv );

#ifdef RELOCATE
	std::vector<Geometry::byte> altbuffer;
#endif
	
#ifdef RELOCATE
	altbuffer.resize(testCase.geom->size());
	testCase.geom->relocate( &(altbuffer[0]) );
	std::memcpy( &(altbuffer[0]), testCase.geom->getBuffer(), testCase.geom->size() );
	std::memset( testCase.geom->getBuffer(), 0, testCase.geom->size() );
#endif
	
#ifdef RELOCATE
	testCase.geom.reset();

	G4VPhysicalVolume *geomRoot =
		(G4VPhysicalVolume*)&(altbuffer[0]);
#else
	G4VPhysicalVolume *geomRoot =
		(G4VPhysicalVolume*)testCase.geom->getBuffer();
#endif
	
	const my_clock_t t0 = my_clock();
	BEGIN_ENERGY_MEASUREMENT;
	cpuexec(
		// number of rounds multiplies problem size, no separate rounds
		testCase.getSize() * testCase.getRounds(),
		&(testCase.input[0]),
		geomRoot,
		&(testCase.output[0]),
		testCase.phys_step );
	const my_clock_t t1 = my_clock();
	END_ENERGY_MEASUREMENT;
	
	std::cerr << "Elapsed: " << tdiffms(t0,t1) << " ms\n";
	
	testCase.outputData( "imgcpu.txt" );
	
	return EXIT_SUCCESS;
#ifndef NOCATCH
  }
  catch ( const std::runtime_error &e )
  {
	std::cerr << e.what() << std::endl;
	return EXIT_FAILURE;
  }
#endif
}
