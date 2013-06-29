
/**
 * Entry point main host file for CUDA versions. Some CUDA host code
 * is also located in cuda.cpp
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <fstream>

#include "myutils.h"

#include "everything.h"
#include "stubParticle.h"
#include "G4Navigator.h"

#include "hostcommons.hpp"

//extern "C" {
typedef struct { const char *err, *fn; int line, errc; } my_cuda_err;

my_cuda_err cudainit( Geometry *geom, int N );
my_cuda_err cudaexec( G4double phys_step, int total, Particle *input, G4double *output );
my_cuda_err cudafinish();
//}

static void checkerr( my_cuda_err ret )
{
	if ( ret.err != NULL )
	{
		std::cerr << "CUDA error \"" << ret.err << "\" ("
			<< ret.errc << ") at " << ret.fn << ":"
			<< ret.line << std::endl;
			
		throw std::runtime_error("CUDA execution failed");
	}
}

extern "C"
{
	void myprint( const char *str )
	{
		/*std::fprintf(stderr, str);*/
	}
	
	void myprint1( const char *str, int n )
	{
		std::fprintf(stderr, str, n);
	}

	void mysleep(int n)
	{
		usleep(n);
	}

	typedef struct { int secs; int usecs; } mytimet;
	
	mytimet mytimer()
	{
		my_clock_t t = my_clock();
		mytimet s = { t.tv_sec, t.tv_usec };
		return s;
	}
	
	void myprinttdiff( mytimet a, mytimet b )
	{
		my_clock_t at = { a.secs, a.usecs };
		my_clock_t bt = { b.secs, b.usecs };
		/*std::fprintf(stderr, "%g ms\n", tdiffms(at,bt));*/
	}
}

int main(int argc, char *argv[])
{
#ifndef NOCATCH
	try
	{
#endif
		std::cout << "LABEL0" << std::endl;
		TestCase testCase( argc, argv );
		const int numInput = testCase.getSize();
		const int totalInput = numInput*testCase.getRounds();
		std::cout << "LABEL1" << std::endl;
		checkerr( cudainit( &(*testCase.geom), numInput ) );
		BEGIN_ENERGY_MEASUREMENT;
		const my_clock_t t1 = my_clock();
		std::cout << "LABEL2" << std::endl;
		checkerr( cudaexec( testCase.phys_step, totalInput, &(testCase.input[0]), &(testCase.output[0]) ) );
		const my_clock_t t2 = my_clock();
		END_ENERGY_MEASUREMENT;
		std::cout << "LABEL3" << std::endl;
		checkerr( cudafinish() );

		std::cerr << "Elapsed: " << tdiffms(t1,t2) << " ms\n\n";
		std::cout << "LABEL4" << std::endl;
		testCase.outputData( "imgcuda.txt" );

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
