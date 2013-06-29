
/**
 * Three-vector operations. Based on the Hep3Vector class in CLHEP
 */

INLINEFUNC
G4ThreeVector G4ThreeVector_create( G4double x, G4double y, G4double z )
{
	G4ThreeVector v =
	#if defined(OPENCL_HOST) && defined(NATIVE_GPU_VECTORS)
		{{x,y,z,0}};
	#else
		#ifdef THREEVECTOR_EXTRAMEMBER
			{x,y,z,0};
		#else
			{x,y,z};
		#endif
	#endif
	return v;
}

INLINEFUNC
G4ThreeVector G4ThreeVector_saxpy( G4double a, G4ThreeVector x, G4ThreeVector y )
{
	return G4ThreeVector_create(
		a*x.x + y.x,
		a*x.y + y.y,
		a*x.z + y.z );
}

INLINEFUNC
G4ThreeVector G4ThreeVector_sum( G4ThreeVector a, G4ThreeVector b )
{
	return G4ThreeVector_create( a.x+b.x, a.y+b.y, a.z+b.z );
}

INLINEFUNC
G4ThreeVector G4ThreeVector_subtract( G4ThreeVector a, G4ThreeVector b )
{
	return G4ThreeVector_create( a.x-b.x, a.y-b.y, a.z-b.z );
}

INLINEFUNC
G4ThreeVector G4ThreeVector_sum_assign( G4ThreeVector *This, G4ThreeVector b )
{
	(*This).x += b.x;
	(*This).y += b.y;
	(*This).z += b.z;
	return *This;
}

INLINEFUNC
G4ThreeVector G4ThreeVector_subtract_assign( G4ThreeVector *This, G4ThreeVector b )
{
	(*This).x -= b.x;
	(*This).y -= b.y;
	(*This).z -= b.z;
	return *This;
}

INLINEFUNC
G4ThreeVector G4ThreeVector_mult_assign( G4ThreeVector *This, G4double m )
{
	(*This).x *= m;
	(*This).y *= m;
	(*This).z *= m;
	return *This;
}

INLINEFUNC
G4ThreeVector G4ThreeVector_negation( G4ThreeVector a )
{
	return G4ThreeVector_create( -a.x, -a.y, -a.z );
}

INLINEFUNC 
G4double G4ThreeVector_mag2( G4ThreeVector v )
{
	return v.x*v.x + v.y*v.y + v.z*v.z;
}

INLINEFUNC 
G4double G4ThreeVector_mag( G4ThreeVector v )
{
	return sqrt(G4ThreeVector_mag2(v));
}

INLINEFUNC 
G4double G4ThreeVector_dot( G4ThreeVector a, G4ThreeVector b )
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

INLINEFUNC 
G4ThreeVector G4ThreeVector_cross( G4ThreeVector a, G4ThreeVector p )
{
	return G4ThreeVector_create( 
		a.y*p.z-p.y*a.z,
		a.z*p.x-p.z*a.x,
		a.x*p.y-p.x*a.y );
}

INLINEFUNC 
G4ThreeVector G4ThreeVector_mult( G4ThreeVector a, G4double m )
{
	return G4ThreeVector_create( a.x*m, a.y*m, a.z*m );
}

INLINEFUNC 
G4ThreeVector G4ThreeVector_unit( G4ThreeVector v )
{
	G4double l = G4ThreeVector_mag(v);
	if ( l > 0 )
		return G4ThreeVector_mult( v, 1.0/l );
	return v;
}

INLINEFUNC
G4bool G4ThreeVector_equal( G4ThreeVector a, G4ThreeVector b )
{	
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

INLINEFUNC
G4double G4ThreeVector_diff2( G4ThreeVector a, G4ThreeVector b )
{
	return G4ThreeVector_mag2( G4ThreeVector_subtract(a,b) );
}

#ifdef HOST_CODE
INLINEFUNC
G4ThreeVector G4ThreeVector_rotate(
				G4ThreeVector *This,
				G4double phi, 
				G4double theta, 
				G4double psi)
{
	G4double rx;
	G4double ry;
	G4double rz;

	G4double sinPhi   = sin( phi   ), cosPhi   = cos( phi   );
	G4double sinTheta = sin( theta ), cosTheta = cos( theta );
	G4double sinPsi   = sin( psi   ), cosPsi   = cos( psi   );

	rx = 	(cosPsi * cosPhi   - cosTheta * sinPsi * sinPhi)   * (*This).x  +
		(cosPsi * sinPhi   + cosTheta * sinPsi * cosPhi)   * (*This).y  +
		(sinPsi * sinTheta)				   * (*This).z  ;

	ry = 	(- sinPsi * cosPhi - cosTheta * cosPsi * sinPhi)   * (*This).x  +
		(- sinPsi * sinPhi + cosTheta * cosPsi * cosPhi)   * (*This).y  +
		(cosPsi * sinTheta)				   * (*This).z  ;

	rz = 	(sinTheta * sinPhi)				   * (*This).x  +
		(- sinTheta * cosPhi)				   * (*This).y  +
		(cosTheta)					   * (*This).z  ;

	(*This).x = rx;
	(*This).y = ry;
	(*This).z = rz;

	return *This;
}
#endif

INLINEFUNC
G4double G4ThreeVector_coord( G4ThreeVector v, EAxis axis )
{
	switch( axis )
	{
	case kXAxis: return v.x;
	case kYAxis: return v.y;
	case kZAxis: return v.z;
	default:
		myAssert(false);
		return 0;
	}
}

INLINEFUNC
void G4ThreeVector_set_coord( G4ThreeVector *v, EAxis axis, G4double val )
{
	switch( axis )
	{
	case kXAxis: v->x = val; break;
	case kYAxis: v->y = val; break;
	case kZAxis: v->z = val; break;
	default:
		myAssert(false);
		break;
	}
}
