
/** Test cases and host input & output */

#ifndef HOST_COMMONS_HPP
#define HOST_COMMONS_HPP

#include <iostream>
#include <sstream>
#include <memory>
#include <ctime>

#define PRINT_OUTPUT
#define PRINT_INPUT
//#define NOCATCH

#ifndef NO_ORB
#include "benchmarks/my/toy1.hpp"
#endif
#include "benchmarks/my/toy2.hpp"
#ifndef ONLY_BOX_AND_ORB
#ifndef NO_ORB
#include "benchmarks/Boyadzhiev/spheres.hpp"
#endif
#include "benchmarks/CMS/CMSModel.hpp"
#endif


#include "cameraJet.c"

// Use boost for random number generation
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/variate_generator.hpp"
#include "boost/version.hpp"
// VoxelHeader and Voxels are only if you need Voxel Navigation
#ifdef ENABLE_VOXEL_NAVIGATION
	#include "G4Voxels.h"
	#include "G4VoxelHeader.cpp"
#endif
#define MEASURE_ENERGY

#ifdef MEASURE_ENERGY
#define BEGIN_ENERGY_MEASUREMENT std::system("./measurebegin.sh")
#define END_ENERGY_MEASUREMENT std::system("./measureend.sh")
#else
#define BEGIN_ENERGY_MEASUREMENT (void)0
#define END_ENERGY_MEASUREMENT (void)0
#endif

class TestCaseBase
{
public:
	enum GeometryType { INVALID_TYPE, TOY1, TOY2, SPHERES, CMS };
	//EDIT: Making geom BasicGeometry
	std::auto_ptr<BasicGeometry> geom;
	G4double phys_step;
	
	/** Return problem size (number of particles per round) */
	virtual int getSize() const = 0;
	
	/** Return the number of rounds */
	virtual int getRounds() const = 0;
	
	/** Output results to file named fn */
	virtual void outputData( const char *fn ) const = 0;
	
	virtual ~TestCaseBase() {}
	
protected:

	bool tryInput( const std::string& src, std::string prefix, std::string &target )
	{
		prefix = "-"+prefix+"=";
		if ( src.find(prefix)==0 )
		{
			target = src.substr(prefix.size());
			return true;
		}
		return false;
	}
	
	int toInt( const std::string& str )
	{
		int i;
		std::istringstream ss( str );
		if ( !(ss >> i) )
			throw std::runtime_error("Invalid integer value \""+str+"\"");
		return i;
	}
	
	double toFloat( std::string str )
	{
		double d;
		std::istringstream ss( str );
		if ( !(ss >> d) )
			throw std::runtime_error("Invalid floating point value \""+str+"\"");
		return d;
	}

	TestCaseBase( int argc, char *argv[] )
	{
		double step = -1;
		int vox1 = K_MIN_VOXEL_VOLUMES_LEVEL_1;
		int vox2 = K_MIN_VOXEL_VOLUMES_LEVEL_2;
		int vox3 = K_MIN_VOXEL_VOLUMES_LEVEL_3;
		GeometryType type = TestCaseBase::INVALID_TYPE;

		for( int i=1; i<argc; ++i )
		{
			const std::string arg(argv[i]);
			std::string val;
			
			if ( tryInput( arg, "test", val ) )
			{
				const std::string name = val;
				if ( name == "toy1" ) type = TestCaseBase::TOY1;
				else if ( name == "toy2" ) type = TestCaseBase::TOY2;
				else if ( name == "spheres" ) type = TestCaseBase::SPHERES;
				else if ( name == "cms" ) type = TestCaseBase::CMS;
			}
			else if ( tryInput( arg, "step", val ) )
				step = toFloat(val);
				
			else if ( tryInput( arg, "vox1", val ) )
				vox1 = toInt( val );
				
			else if ( tryInput( arg, "vox2", val ) )
				vox2 = toInt( val );
				
			else if ( tryInput( arg, "vox3", val ) )
				vox3 = toInt( val );
		}
#ifdef ENABLE_VOXEL_NAVIGATION
		G4VoxelHeader_SetMinVoxelLimits( vox1, vox2, vox3 );
#endif		
		// Physical step scale (maximum relative size of G4 physical steps)
		if ( step <= 0 )
			throw std::runtime_error("Invalid or missing step size");
		
		// Geometry type
		switch( type )
		{
#ifndef NO_ORB
		case TestCaseBase::TOY1:
			//EDIT: Changing to BasicGeomtry to get access to ptrs. Should not change much though.
			geom = std::auto_ptr<BasicGeometry>( new ToyGeometry1 );
			std::cerr << "Test case: Toy Geometry 1" << std::endl;
			break;
#endif

		case TestCaseBase::TOY2:
			//EDIT: Changing to BasicGeomtry to get access to ptrs. Should not change much though.
			geom = std::auto_ptr<BasicGeometry>( new ToyGeometry2 );
			std::cerr << "Test case: Toy Geometry 2" << std::endl;
			break;

#ifndef ONLY_BOX_AND_ORB
#ifndef NO_ORB
		case TestCaseBase::SPHERES:
			//EDIT: Changing to BasicGeomtry to get access to ptrs. Should not change much though.
			geom = std::auto_ptr<BasicGeometry>( new BoyadzhievSpheres );
			std::cerr << "Test case: Spherical Shells" << std::endl;
			break;
#endif

		case TestCaseBase::CMS:
			//EDIT: Changing to BasicGeomtry to get access to ptrs. Should not change much though.
			geom = std::auto_ptr<BasicGeometry>( new CMSModel );
			std::cerr << "Test case: CMS" << std::endl;
			break;
#endif

		default:
			throw std::runtime_error("Invalid or missing geometry type");
		}
		
#ifdef ENABLE_VOXEL_NAVIGATION
		std::cerr
			<< "Physical step scale: " << step << std::endl
			<< "Min. voxels per level: "
			<< G4VoxelHeader_GetMinVoxelVolumesLevel1() << ", "
			<< G4VoxelHeader_GetMinVoxelVolumesLevel2() << ", "
			<< G4VoxelHeader_GetMinVoxelVolumesLevel3() << std::endl;
#endif
		printf ("Geometry Creation Test\n");
		geom->create();
		phys_step = step * geom->getScale();
		printf ("Geometry Creation Test Passed\n");
		std::cerr
			<< "\n\tGeometry buffer size: " << geom->size() << " bytes\n"
			<< "\tNumber of voxel nodes: " << geom->getNumVoxelNodes() << std::endl
			<< "\tPrecision: " <<
				#ifdef DOUBLE_PRECISION
					"double"
				#else
					"single"
				#endif
			<< " (" << sizeof(G4double)*8 << "-bit)\n\n";
	}
};

class RaytracingTestCase : public TestCaseBase
{
public:
	int xres, yres;
	std::vector<StubParticle> input;
	std::vector<G4double> output;
	
	int getSize() const { return xres*yres; }
	int getRounds() const { return 1; }

	RaytracingTestCase( int argc, char *argv[] )
	:
		TestCaseBase( argc, argv )
	{
		xres = yres = -1;
		
		for( int i=1; i<argc; ++i )
		{
			const std::string prefixXres = "-xres=";
			const std::string prefixYres = "-yres=";
			const std::string arg(argv[i]);
			if ( arg.find(prefixXres)==0 )
			{
				std::istringstream ss( arg.substr(prefixXres.size()) );
				if ( !(ss >> xres) )
					throw std::runtime_error("Invalid image dimension \'"+ss.str()+"\'");
			}
			else if ( arg.find(prefixYres)==0 )
			{
				std::istringstream ss( arg.substr(prefixYres.size()) );
				if ( !(ss >> yres) )
					throw std::runtime_error("Invalid image dimension \'"+ss.str()+"\'");
			}
		}
		
		if ( xres <= 0 || yres <= 0 )
			throw std::runtime_error("Invalid or missing image dimension(s)");
			
		std::cerr << "--- Raytracing test case, "
			<< xres << "*" << yres << "\n\n";
			
		input.resize(xres*yres);
		output.resize(xres*yres);

		createCameraJet( &(input[0]), xres, yres, geom->getCamera() );
	}

	// Image output
	void outputData( const char *fn ) const
	{
	#ifdef PRINT_OUTPUT
		std::ofstream outfile( fn );
	
		for ( int y=0; y<yres; ++y )
		{
			for ( int x=0; x<xres; ++x )
			{
				outfile << output[y*xres + x] << "\t";
			}
			outfile << "\n";
		}
	#else
		(void)fn;
	#endif
	}
};

class SimplePhysicsTestCase : public TestCaseBase
{
public:
	std::vector<ParticleWithLifetime> input;
	std::vector<G4double> output;
	int nParticles;
	int nRounds;
	int randomSeed;
	double attenuationPerDensity;
	
	int getSize() const { return nParticles; }
	int getRounds() const { return nRounds; }
	
	SimplePhysicsTestCase( int argc, char *argv[] )
	:
		TestCaseBase( argc, argv )
	{
		nParticles = -1;
		randomSeed = std::time(NULL) ^ std::clock();
		attenuationPerDensity = 1.0;
		nRounds = 1;
		bool raytracerLike = false;
		
		for( int i=1; i<argc; ++i )
		{
			const std::string arg(argv[i]);
			std::string val;
			
			if ( tryInput( arg, "n", val ) )
				nParticles = toInt( val );
			else if ( tryInput( arg, "seed", val ) )
				randomSeed = toInt( val );
			else if ( tryInput( arg, "mean", val ) )
				attenuationPerDensity = toFloat( val );
			else if ( tryInput( arg, "rounds", val ) )
				nRounds = toInt( val );
			else if ( tryInput( arg, "raytracer", val ) )
				raytracerLike = toInt( val ) != 0;
		}
		
		if ( nParticles <= 0 )
			throw std::runtime_error("Invalid or missing problem size");
			
		if ( nRounds <= 0 )
			throw std::runtime_error("Invalid number of rounds");
			
		if ( attenuationPerDensity <= 0.0 )
			throw std::runtime_error("Invalid attenuation coefficient");
			
		std::cerr << "--- Simple physics test, "
			<< nParticles << " particles * " << nRounds << " rounds\n"
			<< "attenuation / density: "
				<< attenuationPerDensity << " m^2/kg\n\n";
		
		output.resize(nParticles*nRounds);
		
		if (raytracerLike)
		{
			std::cerr << "--- \"Raytracer-like particle alignment\"\n\n";
			generateInputRaytracerLike();
		}
		else
		{
			generateInput();
		}
	}
	
	// Image output
	void outputData( const char *fn ) const
	{
	#ifdef PRINT_OUTPUT
		std::ofstream outfile( fn );
		
		for ( int i=0; i<nParticles*nRounds; ++i )
		{
			outfile << output[i];
			
		#ifdef PRINT_INPUT
			outfile << "\t" <<
				input[i].dir.x << "\t" <<
				input[i].dir.y << "\t" <<
				input[i].dir.z << "\t" <<
				input[i].t;
		#endif
				
			outfile << std::endl;
		}
	#else
		(void)fn;
	#endif
	}
	
private:
	// The random number generator, boosts' Mersenne Twister 19937
	boost::mt19937 rng; 

	/** Generate uniform random number in [0,1) */ 
	G4double generateUniform01()
	{
	// TODO: what is the correct version threshold?
	#if BOOST_VERSION <= 103301
		
		// for older versions of boost
		// (a quick, bad and dangrous solution)
		static boost::uniform_01<boost::mt19937> zeroone(rng);
		return zeroone();
		
	#else
		// Uniform real number from interval [0,1)
		boost::uniform_real<> interval;
		
		// Bind the interval to the RNG to get a Uniform(0,1)-generator
		boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
			uniform01generator(rng, interval);
			
		// Generate
		return uniform01generator();
	#endif
	}
	
	/** Generate and Exp(1) distributed random number */
	G4double generateExp()
	{
		return -std::log( 1.0 - generateUniform01() );
	}

	/** Initialize random number generator with seed */
	void randInit( int seed )
	{
		// Prevents nagging with older versions of Boost
		boost::mt19937::result_type s = seed;
		rng.seed( s );
	}
	
	/**
	 * Generate a random unit vector (uniformly distributed
	 * on the surface of the unit sphere
	 */
	G4ThreeVector generateRandomDirection()
	{
		double x, y, z;
		double r2 = 0;
		
		do
		{
			x = generateUniform01()*2.0 - 1.0;
			y = generateUniform01()*2.0 - 1.0;
			z = generateUniform01()*2.0 - 1.0;
			r2 = x*x + y*y + z*z;
		}
		while ( r2 > 1.0 || r2 <= 0 );
		
		r2 = 1.0/sqrt(r2);
		x *= r2;
		y *= r2;
		z *= r2;
		
		return G4ThreeVector_create(x,y,z);
	}
	
	void generateInputRaytracerLike()
	{
		EventOrigin event = geom->getEvent();
		
		CameraParameters camera;
		camera.dist = 1;
		camera.target_x = event.x;
		camera.target_y = event.y+1;
		camera.target_z = event.z;
		
		const int N = output.size();
		
		int xres = 400;
		int yres = N/xres;
		if (xres*yres != N)
			throw std::runtime_error("Invalid number of particles (not divisible by 400) for \"raytracer-like\" physics setup");
			
		std::vector<StubParticle> rays( N );
		createCameraJet( &(rays[0]), xres, yres, camera );
		
		randInit( randomSeed );
		
		for (int i=0; i<N; ++i)
		{
			StubParticle s = rays[i];
			ParticleWithLifetime p;
			p.pos = s.pos;
			p.dir = s.dir;
			p.t = 1 / attenuationPerDensity;
			input.push_back(p);
		}
	}
	
	void generateInput()
	{
		// generate particles
		randInit( randomSeed );
		
		EventOrigin event = geom->getEvent();
		
		const G4ThreeVector orig =
			G4ThreeVector_create(event.x, event.y, event.z);
			
		const int N = output.size();
			
		for (int i=0; i<N; ++i)
		{
			ParticleWithLifetime p;
			p.pos = orig;
			p.dir = generateRandomDirection();
			p.t = generateExp() / attenuationPerDensity;
			
			input.push_back(p);
		}
	}
	
};

// Select test case to use
#ifdef PHYSICS
typedef SimplePhysicsTestCase TestCase;
#else
typedef RaytracingTestCase TestCase;
#endif

#endif
