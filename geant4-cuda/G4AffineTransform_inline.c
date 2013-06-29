
/**
 * Affine trasform operations. Based on a C++ class in CLHEP
 */

INLINEFUNC
void G4AffineTransform_ctor_id( G4AffineTransform *This )
{
	This->rxx = 1;
	This->ryy = 1;
	This->rzz = 1;
	This->rxy = 0;
	This->rxz = 0;
	This->ryx = 0;
	This->ryz = 0;
	This->rzx = 0;
	This->rzy = 0;
	This->tx = 0;
	This->ty = 0;
	This->tz = 0;
}

INLINEFUNC 
void G4AffineTransform_ctor_vector( G4AffineTransform *This, G4ThreeVector tlate)
{
	G4AffineTransform_ctor_id( This );
	This->tx = tlate.x;
	This->ty = tlate.y;
	This->tz = tlate.z;
}

INLINEFUNC 
void G4AffineTransform_ctor_matrix( G4AffineTransform *This, G4RotationMatrix rot)
{
	G4AffineTransform_ctor_id( This );
	This->rxx = rot.rxx;
	This->ryy = rot.ryy;
	This->rzz = rot.rzz;
	This->rxy = rot.rxy;
	This->rxz = rot.rxz;
	This->ryx = rot.ryx;
	This->ryz = rot.ryz;
	This->rzx = rot.rzx;
	This->rzy = rot.rzy;
}

INLINEFUNC 
void G4AffineTransform_ctor_full(
	G4AffineTransform *This, G4RotationMatrix rot, G4ThreeVector tlate )
{
	This->rxx = rot.rxx;
	This->ryy = rot.ryy;
	This->rzz = rot.rzz;
	This->rxy = rot.rxy;
	This->rxz = rot.rxz;
	This->ryx = rot.ryx;
	This->ryz = rot.ryz;
	This->rzx = rot.rzx;
	This->rzy = rot.rzy;
	This->tx = tlate.x;
	This->ty = tlate.y;
	This->tz = tlate.z;
}

INLINEFUNC 
void G4AffineTransform_ctor_ptr(
	G4AffineTransform *This, const G4RotationMatrix *rot, G4ThreeVector tlate )
{
	if (rot) G4AffineTransform_ctor_full( This, *rot, tlate );
	else G4AffineTransform_ctor_vector( This, tlate );
}

INLINEFUNC
void G4AffineTransform_ctor_elements(
		G4AffineTransform *This,
		const G4double prxx,const G4double prxy,const G4double prxz,
		const G4double pryx,const G4double pryy,const G4double pryz,
		const G4double przx,const G4double przy,const G4double przz,
		const G4double ptx,const G4double pty,const G4double ptz)
{
	This->rxx = prxx;
	This->ryy = pryy;
	This->rzz = przz;
	This->rxy = prxy;
	This->rxz = prxz;
	This->ryx = pryx;
	This->ryz = pryz;
	This->rzx = przx;
	This->rzy = przy;
	This->tx = ptx;
	This->ty = pty;
	This->tz = ptz;
}

INLINEFUNC
G4AffineTransform G4AffineTransform_create_id(void)
{
	G4AffineTransform t;
	G4AffineTransform_ctor_id(&t);
	return t;
}

INLINEFUNC 
G4AffineTransform G4AffineTransform_create_vector(G4ThreeVector tlate)
{
	G4AffineTransform t;
	G4AffineTransform_ctor_vector(&t,tlate);
	return t;
}

INLINEFUNC 
G4AffineTransform G4AffineTransform_create_matrix( G4RotationMatrix rot )
{
	G4AffineTransform t;
	G4AffineTransform_ctor_matrix(&t,rot);
	return t;
}

INLINEFUNC 
G4AffineTransform G4AffineTransform_create_full(
	G4RotationMatrix rot, G4ThreeVector tlate )
{
	G4AffineTransform t;
	G4AffineTransform_ctor_full(&t,rot,tlate);
	return t;
}

INLINEFUNC 
G4AffineTransform G4AffineTransform_create_ptr(
	const G4RotationMatrix *rot, G4ThreeVector tlate )
{
	G4AffineTransform t;
	G4AffineTransform_ctor_ptr(&t,rot,tlate);
	return t;
}

INLINEFUNC
G4AffineTransform G4AffineTransform_create_elements(
		const G4double prxx,const G4double prxy,const G4double prxz,
		const G4double pryx,const G4double pryy,const G4double pryz,
		const G4double przx,const G4double przy,const G4double przz,
		const G4double ptx,const G4double pty,const G4double ptz)
{
	G4AffineTransform t;
	G4AffineTransform_ctor_elements(&t,
		prxx,prxy,prxz,
		pryx,pryy,pryz,
		przx,przy,przz,
		ptx,pty,ptz);
	return t;
}


#define tf (*ptrtf)
#define tf1 (*ptrtf1)
#define tf2 (*ptrtf2)

/*INLINEFUNC G4AffineTransform
G4AffineTransform_operator_mult (const G4AffineTransform *This, const G4AffineTransform* ptrtf)
{
        return G4AffineTransform_create_elements(
        This->rxx*tf.rxx+This->rxy*tf.ryx+This->rxz*tf.rzx,
        This->rxx*tf.rxy+This->rxy*tf.ryy+This->rxz*tf.rzy,
        This->rxx*tf.rxz+This->rxy*tf.ryz+This->rxz*tf.rzz,

        This->ryx*tf.rxx+This->ryy*tf.ryx+This->ryz*tf.rzx,
        This->ryx*tf.rxy+This->ryy*tf.ryy+This->ryz*tf.rzy,
        This->ryx*tf.rxz+This->ryy*tf.ryz+This->ryz*tf.rzz,

        This->rzx*tf.rxx+This->rzy*tf.ryx+This->rzz*tf.rzx,
        This->rzx*tf.rxy+This->rzy*tf.ryy+This->rzz*tf.rzy,
        This->rzx*tf.rxz+This->rzy*tf.ryz+This->rzz*tf.rzz,
        
        This->tx*tf.rxx+This->ty*tf.ryx+This->tz*tf.rzx+tf.tx,
        This->tx*tf.rxy+This->ty*tf.ryy+This->tz*tf.rzy+tf.ty,
        This->tx*tf.rxz+This->ty*tf.ryz+This->tz*tf.rzz+tf.tz);
}

INLINEFUNC G4AffineTransform
G4AffineTransform_operator_mult_assign (
		G4AffineTransform *This, const G4AffineTransform *ptrtf)
{
         // Use temporaries for `in place' compound transform computation

        G4double nrxx=This->rxx*tf.rxx+This->rxy*tf.ryx+This->rxz*tf.rzx;
        G4double nrxy=This->rxx*tf.rxy+This->rxy*tf.ryy+This->rxz*tf.rzy;
        G4double nrxz=This->rxx*tf.rxz+This->rxy*tf.ryz+This->rxz*tf.rzz;

        G4double nryx=This->ryx*tf.rxx+This->ryy*tf.ryx+This->ryz*tf.rzx;
        G4double nryy=This->ryx*tf.rxy+This->ryy*tf.ryy+This->ryz*tf.rzy;
        G4double nryz=This->ryx*tf.rxz+This->ryy*tf.ryz+This->ryz*tf.rzz;

        G4double nrzx=This->rzx*tf.rxx+This->rzy*tf.ryx+This->rzz*tf.rzx;
        G4double nrzy=This->rzx*tf.rxy+This->rzy*tf.ryy+This->rzz*tf.rzy;
        G4double nrzz=This->rzx*tf.rxz+This->rzy*tf.ryz+This->rzz*tf.rzz;
        
        G4double ntx=This->tx*tf.rxx+This->ty*tf.ryx+This->tz*tf.rzx+tf.tx;
        G4double nty=This->tx*tf.rxy+This->ty*tf.ryy+This->tz*tf.rzy+tf.ty;
        G4double ntz=This->tx*tf.rxz+This->ty*tf.ryz+This->tz*tf.rzz+tf.tz;

        This->tx=ntx; This->ty=nty; This->tz=ntz;
        This->rxx=nrxx; This->rxy=nrxy; This->rxz=nrxz;
        This->ryx=nryx; This->ryy=nryy; This->ryz=nryz;
        This->rzx=nrzx; This->rzy=nrzy; This->rzz=nrzz;

        return *This;
}*/

/*INLINEFUNC G4AffineTransform
G4AffineTransform_Product(
		G4AffineTransform *This,
		const G4AffineTransform* ptrtf1,
		const G4AffineTransform* ptrtf2)
{
        This->rxx=tf1.rxx*tf2.rxx + tf1.rxy*tf2.ryx + tf1.rxz*tf2.rzx;
        This->rxy=tf1.rxx*tf2.rxy + tf1.rxy*tf2.ryy + tf1.rxz*tf2.rzy;
        This->rxz=tf1.rxx*tf2.rxz + tf1.rxy*tf2.ryz + tf1.rxz*tf2.rzz;

        This->ryx=tf1.ryx*tf2.rxx + tf1.ryy*tf2.ryx + tf1.ryz*tf2.rzx;
        This->ryy=tf1.ryx*tf2.rxy + tf1.ryy*tf2.ryy + tf1.ryz*tf2.rzy;
        This->ryz=tf1.ryx*tf2.rxz + tf1.ryy*tf2.ryz + tf1.ryz*tf2.rzz;

        This->rzx=tf1.rzx*tf2.rxx + tf1.rzy*tf2.ryx + tf1.rzz*tf2.rzx;
        This->rzy=tf1.rzx*tf2.rxy + tf1.rzy*tf2.ryy + tf1.rzz*tf2.rzy;
        This->rzz=tf1.rzx*tf2.rxz + tf1.rzy*tf2.ryz + tf1.rzz*tf2.rzz;
        
        This->tx=tf1.tx*tf2.rxx + tf1.ty*tf2.ryx + tf1.tz*tf2.rzx   + tf2.tx;
        This->ty=tf1.tx*tf2.rxy + tf1.ty*tf2.ryy + tf1.tz*tf2.rzy   + tf2.ty;
        This->tz=tf1.tx*tf2.rxz + tf1.ty*tf2.ryz + tf1.tz*tf2.rzz   + tf2.tz; 
        
        return *This;
}*/

INLINEFUNC G4AffineTransform
G4AffineTransform_InverseProduct(
	G4AffineTransform *This,
	const G4AffineTransform* ptrtf1,
	const G4AffineTransform* ptrtf2)
{
        G4double itf2tx = - tf2.tx*tf2.rxx - tf2.ty*tf2.rxy - tf2.tz*tf2.rxz;
        G4double itf2ty = - tf2.tx*tf2.ryx - tf2.ty*tf2.ryy - tf2.tz*tf2.ryz;
        G4double itf2tz = - tf2.tx*tf2.rzx - tf2.ty*tf2.rzy - tf2.tz*tf2.rzz;

        This->rxx=tf1.rxx*tf2.rxx+tf1.rxy*tf2.rxy+tf1.rxz*tf2.rxz;
        This->rxy=tf1.rxx*tf2.ryx+tf1.rxy*tf2.ryy+tf1.rxz*tf2.ryz;
        This->rxz=tf1.rxx*tf2.rzx+tf1.rxy*tf2.rzy+tf1.rxz*tf2.rzz;

        This->ryx=tf1.ryx*tf2.rxx+tf1.ryy*tf2.rxy+tf1.ryz*tf2.rxz;
        This->ryy=tf1.ryx*tf2.ryx+tf1.ryy*tf2.ryy+tf1.ryz*tf2.ryz;
        This->ryz=tf1.ryx*tf2.rzx+tf1.ryy*tf2.rzy+tf1.ryz*tf2.rzz;

        This->rzx=tf1.rzx*tf2.rxx+tf1.rzy*tf2.rxy+tf1.rzz*tf2.rxz;
        This->rzy=tf1.rzx*tf2.ryx+tf1.rzy*tf2.ryy+tf1.rzz*tf2.ryz;
        This->rzz=tf1.rzx*tf2.rzx+tf1.rzy*tf2.rzy+tf1.rzz*tf2.rzz;
        
        This->tx=tf1.tx*tf2.rxx+tf1.ty*tf2.rxy+tf1.tz*tf2.rxz+itf2tx;
        This->ty=tf1.tx*tf2.ryx+tf1.ty*tf2.ryy+tf1.tz*tf2.ryz+itf2ty;
        This->tz=tf1.tx*tf2.rzx+tf1.ty*tf2.rzy+tf1.tz*tf2.rzz+itf2tz;

        return *This;
}

INLINEFUNC
G4ThreeVector G4AffineTransform_TransformPoint(const G4AffineTransform *This, G4ThreeVector vec)
{
        return G4ThreeVector_create(
			vec.x*This->rxx + vec.y*This->ryx + vec.z*This->rzx + This->tx,
			vec.x*This->rxy + vec.y*This->ryy + vec.z*This->rzy + This->ty,
			vec.x*This->rxz + vec.y*This->ryz + vec.z*This->rzz + This->tz  );
}

INLINEFUNC
G4ThreeVector G4AffineTransform_TransformAxis(const G4AffineTransform *This, G4ThreeVector axis)
{
		return G4ThreeVector_create(
			axis.x*This->rxx + axis.y*This->ryx + axis.z*This->rzx,
			axis.x*This->rxy + axis.y*This->ryy + axis.z*This->rzy,
			axis.x*This->rxz + axis.y*This->ryz + axis.z*This->rzz );
}

INLINEFUNC
G4AffineTransform G4AffineTransform_Inverse(const G4AffineTransform *This)
{
        return G4AffineTransform_create_elements(
				This->rxx, This->ryx, This->rzx,
				This->rxy, This->ryy, This->rzy,
				This->rxz, This->ryz, This->rzz,

				-This->tx*This->rxx - This->ty*This->rxy - This->tz*This->rxz,
				-This->tx*This->ryx - This->ty*This->ryy - This->tz*This->ryz,
				-This->tx*This->rzx - This->ty*This->rzy - This->tz*This->rzz  );
}

INLINEFUNC
G4AffineTransform G4AffineTransform_Invert(G4AffineTransform *This)
{
        G4double v1 = -This->tx*This->rxx - This->ty*This->rxy - This->tz*This->rxz;
        G4double v2 = -This->tx*This->ryx - This->ty*This->ryy - This->tz*This->ryz;
        G4double v3 = -This->tx*This->rzx - This->ty*This->rzy - This->tz*This->rzz;

        This->tx=v1; This->ty=v2; This->tz=v3;

        G4double tmp1=This->ryx; This->ryx=This->rxy; This->rxy=tmp1;
        G4double tmp2=This->rzx; This->rzx=This->rxz; This->rxz=tmp2;
        G4double tmp3=This->rzy; This->rzy=This->ryz; This->ryz=tmp3;

        return *This;
}

INLINEFUNC
G4ThreeVector G4AffineTransform_NetTranslation(const G4AffineTransform *This)
{
        return G4ThreeVector_create(This->tx,This->ty,This->tz);
}

INLINEFUNC
G4bool G4AffineTransform_IsRotated(const G4AffineTransform *This)
{
        return (This->rxx==1.0 && This->ryy==1.0 && This->rzz==1.0) ? false : true;
}


/*INLINEFUNC
G4AffineTransform G4AffineTransform_operator_sum_vector_assign(G4AffineTransform *This, G4ThreeVector tlate)
{
        This->tx += tlate.x;
        This->ty += tlate.y;
        This->tz += tlate.z;

        return *This;
}

INLINEFUNC
G4AffineTransform G4AffineTransform_operator_subtract_vector_assign(G4AffineTransform *This, G4ThreeVector tlate)
{
        This->tx -= tlate.x;
        This->ty -= tlate.y;
        This->tz -= tlate.z;

        return *This;
}*/

/*INLINEFUNC
G4double G4AffineTransform_operator_lookup (const G4AffineTransform *This, const G4int n)
{
        G4double v = 0.0;
        switch(n)
                {
                case 0:
                        v=This->rxx;
                        break;
                case 1:
                        v=This->rxy;
                        break;
                case 2:
                        v=This->rxz;
                        break;
                case 4:
                        v=This->ryx;
                        break;
                case 5:
                        v=This->ryy;
                        break;
                case 6:
                        v=This->ryz;
                        break;
                case 8:
                        v=This->rzx;
                        break;
                case 9:
                        v=This->rzy;
                        break;
                case 10:
                        v=This->rzz;
                        break;
                case 12:
                        v=This->tx;
                        break;
                case 13:
                        v=This->ty;
                        break;
                case 14:
                        v=This->tz;
                        break;
                case 3:
                case 7:
                case 11:
                        break;
                case 15:
                        v=1.0;
                        break;
                }
        return v;
}

INLINEFUNC
G4bool G4AffineTransform_IsRotated(const G4AffineTransform *This)
{
        return (This->rxx==1.0 && This->ryy==1.0 && This->rzz==1.0) ? false : true;
}

INLINEFUNC 
G4bool G4AffineTransform_IsTranslated(const G4AffineTransform *This)
{
        return (This->tx || This->ty || This->tz) ? true:false;
}

INLINEFUNC
G4RotationMatrix G4AffineTransform_NetRotation(const G4AffineTransform *This)
{
  G4RotationMatrix m = G4RotationMatrix_create_id();
  return G4RotationMatrix_rotateAxes(
			&m,
			G4ThreeVector_create(This->rxx,This->ryx,This->rzx),
            G4ThreeVector_create(This->rxy,This->ryy,This->rzy),
            G4ThreeVector_create(This->rxz,This->ryz,This->rzz));
}
*/

#undef tf
#undef tf1
#undef tf2
