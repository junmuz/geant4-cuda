
#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <stdlib.h>
/*
typedef struct timeval my_clock_t;
inline static my_clock_t my_clock()
{
	struct timeval t;
	gettimeofday( &t, NULL );
	return t;
}

/** Time difference 
inline static double tdiff( my_clock_t a, my_clock_t b )
{
	int sd = (signed)b.tv_sec - a.tv_sec;
	int ud = (signed)b.tv_usec - a.tv_usec;
	return ud * 1e-6 + sd;
}*/

/** Time difference in milliseconds */
#define tdiffms(a,b) (tdiff(a,b)*1e3)

/** Million kernel instances per second */
#define ops( N, a,b ) ((N/tdiff(a,b))*1e-6)

#define sizeMB( x ) ((x)>>20)

#ifdef __cplusplus

#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

std::string fileToString( std::string fn )
{
	std::ostringstream ss;
	std::ifstream f( fn.c_str() );
	if (!f) throw std::runtime_error( "could not open "+fn );
	ss << f.rdbuf();
	return ss.str(); 
}

/** Compare output buffers */
template <class ConA, class ConB>
static bool compareBuffers( const ConA& a, const ConB& b )
{
	const double Linf_tol = 1e-5, L1_tol = 1e-6 * a.size();
	
	typename ConA::const_iterator i = a.begin();
	typename ConB::const_iterator j = b.begin();
	double L1 = 0, Linf = 0;
	
	for ( ; i != a.end(); ++i, ++j )
	{
		const double d = std::fabs( *i - *j );
		L1 += d;
		if ( d > Linf ) Linf = d;
		
		// std::cout << "\t" << *i << " " << *j << "\n";
	}
	
	if ( L1 > L1_tol || Linf > Linf_tol )
	{
		std::cout << "Deviation:\n"
			<< "\tL1:\t" << L1 << "\n"
			<< "\tLinf:\t" << Linf << "\n";
		
		std::cout << "Deviations outside allowed tolerance !!!\n";
		return false;
	}
	else return true;
}

#endif

#endif
