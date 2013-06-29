
/**
 * Rotation matrix operations. Based on a C++ class in CLHEP
 */

#ifdef HOST_CODE

INLINEFUNC
G4RotationMatrix G4RotationMatrix_create_id( void )
{
	G4RotationMatrix r =
		{ 1,0,0, 0,1,0, 0,0,1
		#ifndef DOUBLE_PRECISION
		, 0
		#endif
		};
	return r;
}

#endif

INLINEFUNC
G4RotationMatrix G4RotationMatrix_create_elements
			(G4double mxx, G4double mxy, G4double mxz,
			 G4double myx, G4double myy, G4double myz,
			 G4double mzx, G4double mzy, G4double mzz)
{
	G4RotationMatrix r =
		{ mxx,mxy,mxz, myx,myy,myz, mzx,mzy,mzz
		#ifndef DOUBLE_PRECISION
		, 0
		#endif
		};
	return r;
}

INLINEFUNC
G4ThreeVector G4RotationMatrix_apply (const G4RotationMatrix *This, G4ThreeVector p)
{
  return G4ThreeVector_create(
					This->rxx*p.x + This->rxy*p.y + This->rxz*p.z,
                    This->ryx*p.x + This->ryy*p.y + This->ryz*p.z,
                    This->rzx*p.x + This->rzy*p.y + This->rzz*p.z);
}

INLINEFUNC
G4RotationMatrix G4RotationMatrix_mult (const G4RotationMatrix *This, const G4RotationMatrix *other)
{
	#define r (*other)
	return G4RotationMatrix_create_elements(
		This->rxx*r.rxx + This->rxy*r.ryx + This->rxz*r.rzx,
		This->rxx*r.rxy + This->rxy*r.ryy + This->rxz*r.rzy,
		This->rxx*r.rxz + This->rxy*r.ryz + This->rxz*r.rzz,
		This->ryx*r.rxx + This->ryy*r.ryx + This->ryz*r.rzx,
		This->ryx*r.rxy + This->ryy*r.ryy + This->ryz*r.rzy,
		This->ryx*r.rxz + This->ryy*r.ryz + This->ryz*r.rzz,
		This->rzx*r.rxx + This->rzy*r.ryx + This->rzz*r.rzx,
		This->rzx*r.rxy + This->rzy*r.ryy + This->rzz*r.rzy,
		This->rzx*r.rxz + This->rzy*r.ryz + This->rzz*r.rzz );
	#undef r
}

/*INLINEFUNC
G4RotationMatrix G4RotationMatrix_mult_assign (G4RotationMatrix *This, const G4RotationMatrix *other)
{
	*This = G4RotationMatrix_mult(This,other);
	return *This;
}*/

INLINEFUNC
G4RotationMatrix G4RotationMatrix_transform(G4RotationMatrix *This, const G4RotationMatrix *other)
{
	*This = G4RotationMatrix_mult(other,This);
	return *This;
}

INLINEFUNC
G4RotationMatrix G4RotationMatrix_inverse(const G4RotationMatrix *This)
{
	return G4RotationMatrix_create_elements(
		This->rxx, This->ryx, This->rzx, 
		This->rxy, This->ryy, This->rzy, 
		This->rxz, This->ryz, This->rzz );
}

INLINEFUNC
G4RotationMatrix G4RotationMatrix_invert(G4RotationMatrix *This)
{
	return *This = G4RotationMatrix_inverse(This);
}

#ifdef HOST_CODE

INLINEFUNC
G4RotationMatrix G4RotationMatrix_rotateAxes(
		G4RotationMatrix *This,
		G4ThreeVector newX,
		G4ThreeVector newY,
		G4ThreeVector newZ)
{
	G4RotationMatrix m = G4RotationMatrix_create_elements(
			newX.x, newY.x, newZ.x,
			newX.y, newY.y, newZ.y,
			newX.z, newY.z, newZ.z);
	return G4RotationMatrix_transform(This,&m);
}

INLINEFUNC
void G4RotationMatrix_rotateX(G4RotationMatrix *This, G4double a)
{
  G4double c = std::cos(a);
  G4double s = std::sin(a);
  G4double x = This->ryx, y = This->ryy, z = This->ryz; 
  This->ryx = c*x - s*This->rzx;
  This->ryy = c*y - s*This->rzy;
  This->ryz = c*z - s*This->rzz;
  This->rzx = s*x + c*This->rzx;
  This->rzy = s*y + c*This->rzy;
  This->rzz = s*z + c*This->rzz;
}

INLINEFUNC
void G4RotationMatrix_rotateY(G4RotationMatrix *This, G4double a)
{
  G4double c = std::cos(a);
  G4double s = std::sin(a);
  G4double x = This->rzx, y = This->rzy, z = This->rzz; 
  This->rzx = c*x - s*This->rxx;
  This->rzy = c*y - s*This->rxy;
  This->rzz = c*z - s*This->rxz;
  This->rxx = s*x + c*This->rxx;
  This->rxy = s*y + c*This->rxy;
  This->rxz = s*z + c*This->rxz;
}

INLINEFUNC
void G4RotationMatrix_rotateZ(G4RotationMatrix *This, G4double a)
{
  G4double c = std::cos(a);
  G4double s = std::sin(a);
  G4double x = This->rxx, y = This->rxy, z = This->rxz; 
  This->rxx = c*x - s*This->ryx;
  This->rxy = c*y - s*This->ryy;
  This->rxz = c*z - s*This->ryz;
  This->ryx = s*x + c*This->ryx;
  This->ryy = s*y + c*This->ryy;
  This->ryz = s*z + c*This->ryz;
}

#endif
