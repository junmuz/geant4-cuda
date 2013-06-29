 typedef float G4double;
typedef float G4float;
typedef int G4int;
typedef int G4bool;
typedef long G4long;
const G4double kInfinity = 1.0E37;
const int BlockSize = 32;
const int Multiplier = 4;
const G4double twopi = 2.0*3.14159265358979323846264338327;
const G4double kMinExitingNormalCosine = 1E-3;
typedef enum {kOutside,kSurface,kInside} EInside;
typedef enum {kNormal,kReplica,kParameterised} EVolume;
typedef enum {kXAxis,kYAxis,kZAxis,kRho,kRadial3D,kPhi,kUndefined} EAxis;
typedef enum { kBox = 0 , kOrb, kTubs, kCons, kPolyCone, Solidcount } ESolid;
typedef struct
{
 G4double x,y,z;

 G4double w;

}
G4ThreeVector;
__device__
G4ThreeVector G4ThreeVector_create( G4double x, G4double y, G4double z )
{
 G4ThreeVector v =
   {x,y,z,0};
 return v;
}
__device__
G4ThreeVector G4ThreeVector_saxpy( G4double a, G4ThreeVector x, G4ThreeVector y )
{
 return G4ThreeVector_create(
  a*x.x + y.x,
  a*x.y + y.y,
  a*x.z + y.z );
}
__device__
G4ThreeVector G4ThreeVector_sum( G4ThreeVector a, G4ThreeVector b )
{
 return G4ThreeVector_create( a.x+b.x, a.y+b.y, a.z+b.z );
}
__device__
G4ThreeVector G4ThreeVector_subtract( G4ThreeVector a, G4ThreeVector b )
{
 return G4ThreeVector_create( a.x-b.x, a.y-b.y, a.z-b.z );
}
__device__
G4ThreeVector G4ThreeVector_sum_assign( G4ThreeVector *This, G4ThreeVector b )
{
 (*This).x += b.x;
 (*This).y += b.y;
 (*This).z += b.z;
 return *This;
}
__device__
G4ThreeVector G4ThreeVector_subtract_assign( G4ThreeVector *This, G4ThreeVector b )
{
 (*This).x -= b.x;
 (*This).y -= b.y;
 (*This).z -= b.z;
 return *This;
}
__device__
G4ThreeVector G4ThreeVector_mult_assign( G4ThreeVector *This, G4double m )
{
 (*This).x *= m;
 (*This).y *= m;
 (*This).z *= m;
 return *This;
}
__device__
G4ThreeVector G4ThreeVector_negation( G4ThreeVector a )
{
 return G4ThreeVector_create( -a.x, -a.y, -a.z );
}
__device__
G4double G4ThreeVector_mag2( G4ThreeVector v )
{
 return v.x*v.x + v.y*v.y + v.z*v.z;
}
__device__
G4double G4ThreeVector_mag( G4ThreeVector v )
{
 return sqrt(G4ThreeVector_mag2(v));
}

__device__
G4double G4ThreeVector_dot( G4ThreeVector a, G4ThreeVector b )
{
 return a.x*b.x + a.y*b.y + a.z*b.z;
}

__device__
G4ThreeVector G4ThreeVector_cross( G4ThreeVector a, G4ThreeVector p )
{
 return G4ThreeVector_create(
  a.y*p.z-p.y*a.z,
  a.z*p.x-p.z*a.x,
  a.x*p.y-p.x*a.y );
}

__device__
G4ThreeVector G4ThreeVector_mult( G4ThreeVector a, G4double m )
{
 return G4ThreeVector_create( a.x*m, a.y*m, a.z*m );
}

__device__
G4ThreeVector G4ThreeVector_unit( G4ThreeVector v )
{
 G4double l = G4ThreeVector_mag(v);
 if ( l > 0 )
  return G4ThreeVector_mult( v, 1.0/l );
 return v;
}

__device__
G4bool G4ThreeVector_equal( G4ThreeVector a, G4ThreeVector b )
{
 return a.x == b.x && a.y == b.y && a.z == b.z;
}

__device__
G4double G4ThreeVector_diff2( G4ThreeVector a, G4ThreeVector b )
{
 return G4ThreeVector_mag2( G4ThreeVector_subtract(a,b) );
}
__device__
G4double G4ThreeVector_coord( G4ThreeVector v, EAxis axis )
{
 switch( axis )
 {
 case kXAxis: return v.x;
 case kYAxis: return v.y;
 case kZAxis: return v.z;
 default:
  (void)0;
  return 0;
 }
}
__device__
void G4ThreeVector_set_coord( G4ThreeVector *v, EAxis axis, G4double val )
{
 switch( axis )
 {
 case kXAxis: v->x = val; break;
 case kYAxis: v->y = val; break;
 case kZAxis: v->z = val; break;
 default:
  (void)0;
  break;
 }
}
typedef struct
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
typedef StubParticle Particle;
__device__ void Prefix_Sum ( int * input, int * output, int length)
{
 int tid = (blockIdx.x * blockDim.x + threadIdx.x);
 int offset = 1;
 if ( tid< length)
  output[tid] = input[ tid ];
 for(int d = length>>1; d > 0; d >>=1)
 {
__syncthreads();
  if(tid<d)
  {
   int ai = offset*(2*tid + 1) - 1;
   int bi = offset*(2*tid + 2) - 1;
   output[bi] += output[ai];
  }
  offset *= 2;
 }
 if(tid == 0)
 {
  output[length - 1] = 0;
 }
 for(int d = 1; d < length ; d *= 2)
 {
  offset >>=1;
  __syncthreads();
  if(tid < d)
  {
   int ai = offset*(2*tid + 1) - 1;
   int bi = offset*(2*tid + 2) - 1;
   float t = output[ai];
   output[ai] = output[bi];
   output[bi] += t;
  }
 }
__syncthreads();
}
__device__
G4bool NoStepReduction( G4bool * noStepArray, int length )
{
 int tid = (blockIdx.x * blockDim.x + threadIdx.x);
 int offset = 1;
 for(int d = length>>1; d > 0; d >>=1)
 {
     __syncthreads();
  if(tid<d)
  {
   int ai = offset*(2*tid + 1) - 1;
   int bi = offset*(2*tid + 2) - 1;
   noStepArray[bi] = (noStepArray[ai] || noStepArray[bi]);
  }
  offset *= 2;
 }
 G4bool result = noStepArray[ length - 1 ];
 __syncthreads();
 return result;
}
typedef struct
{
 G4double
  rxx, rxy, rxz,
  ryx, ryy, ryz,
  rzx, rzy, rzz;
 G4double align;
}
G4RotationMatrix;
__device__
G4RotationMatrix G4RotationMatrix_create_elements
   (G4double mxx, G4double mxy, G4double mxz,
    G4double myx, G4double myy, G4double myz,
    G4double mzx, G4double mzy, G4double mzz)
{
 G4RotationMatrix r =
  { mxx,mxy,mxz, myx,myy,myz, mzx,mzy,mzz
  , 0
  };
 return r;
}
__device__
G4ThreeVector G4RotationMatrix_apply (const G4RotationMatrix *This, G4ThreeVector p)
{
  return G4ThreeVector_create(
     This->rxx*p.x + This->rxy*p.y + This->rxz*p.z,
                    This->ryx*p.x + This->ryy*p.y + This->ryz*p.z,
                    This->rzx*p.x + This->rzy*p.y + This->rzz*p.z);
}
__device__
G4RotationMatrix G4RotationMatrix_mult (const G4RotationMatrix *This, const G4RotationMatrix *other)
{
 return G4RotationMatrix_create_elements(
  This->rxx*(*other).rxx + This->rxy*(*other).ryx + This->rxz*(*other).rzx,
  This->rxx*(*other).rxy + This->rxy*(*other).ryy + This->rxz*(*other).rzy,
  This->rxx*(*other).rxz + This->rxy*(*other).ryz + This->rxz*(*other).rzz,
  This->ryx*(*other).rxx + This->ryy*(*other).ryx + This->ryz*(*other).rzx,
  This->ryx*(*other).rxy + This->ryy*(*other).ryy + This->ryz*(*other).rzy,
  This->ryx*(*other).rxz + This->ryy*(*other).ryz + This->ryz*(*other).rzz,
  This->rzx*(*other).rxx + This->rzy*(*other).ryx + This->rzz*(*other).rzx,
  This->rzx*(*other).rxy + This->rzy*(*other).ryy + This->rzz*(*other).rzy,
  This->rzx*(*other).rxz + This->rzy*(*other).ryz + This->rzz*(*other).rzz );
}
__device__
G4RotationMatrix G4RotationMatrix_transform(G4RotationMatrix *This, const G4RotationMatrix *other)
{
 *This = G4RotationMatrix_mult(other,This);
 return *This;
}
__device__
G4RotationMatrix G4RotationMatrix_inverse(const G4RotationMatrix *This)
{
 return G4RotationMatrix_create_elements(
  This->rxx, This->ryx, This->rzx,
  This->rxy, This->ryy, This->rzy,
  This->rxz, This->ryz, This->rzz );
}
__device__
G4RotationMatrix G4RotationMatrix_invert(G4RotationMatrix *This)
{
 return *This = G4RotationMatrix_inverse(This);
}
typedef struct
{
  G4double rxx,rxy,rxz;
  G4double ryx,ryy,ryz;
  G4double rzx,rzy,rzz;
  G4double tx,ty,tz;
}
G4AffineTransform;
__device__
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
__device__
void G4AffineTransform_ctor_vector( G4AffineTransform *This, G4ThreeVector tlate)
{
 G4AffineTransform_ctor_id( This );
 This->tx = tlate.x;
 This->ty = tlate.y;
 This->tz = tlate.z;
}
__device__
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
__device__
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
__device__
void G4AffineTransform_ctor_ptr(
 G4AffineTransform *This, const G4RotationMatrix *rot, G4ThreeVector tlate )
{
 if (rot) G4AffineTransform_ctor_full( This, *rot, tlate );
 else G4AffineTransform_ctor_vector( This, tlate );
}
__device__
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
__device__
G4AffineTransform G4AffineTransform_create_id(void)
{
 G4AffineTransform t;
 G4AffineTransform_ctor_id(&t);
 return t;
}
__device__
G4AffineTransform G4AffineTransform_create_vector(G4ThreeVector tlate)
{
 G4AffineTransform t;
 G4AffineTransform_ctor_vector(&t,tlate);
 return t;
}
__device__
G4AffineTransform G4AffineTransform_create_matrix( G4RotationMatrix rot )
{
 G4AffineTransform t;
 G4AffineTransform_ctor_matrix(&t,rot);
 return t;
}
__device__
G4AffineTransform G4AffineTransform_create_full(
 G4RotationMatrix rot, G4ThreeVector tlate )
{
 G4AffineTransform t;
 G4AffineTransform_ctor_full(&t,rot,tlate);
 return t;
}
__device__
G4AffineTransform G4AffineTransform_create_ptr(
 const G4RotationMatrix *rot, G4ThreeVector tlate )
{
 G4AffineTransform t;
 G4AffineTransform_ctor_ptr(&t,rot,tlate);
 return t;
}
__device__
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
__device__ G4AffineTransform
G4AffineTransform_InverseProduct(
 G4AffineTransform *This,
 const G4AffineTransform* ptrtf1,
 const G4AffineTransform* ptrtf2)
{
        G4double itf2tx = - (*ptrtf2).tx*(*ptrtf2).rxx - (*ptrtf2).ty*(*ptrtf2).rxy - (*ptrtf2).tz*(*ptrtf2).rxz;
        G4double itf2ty = - (*ptrtf2).tx*(*ptrtf2).ryx - (*ptrtf2).ty*(*ptrtf2).ryy - (*ptrtf2).tz*(*ptrtf2).ryz;
        G4double itf2tz = - (*ptrtf2).tx*(*ptrtf2).rzx - (*ptrtf2).ty*(*ptrtf2).rzy - (*ptrtf2).tz*(*ptrtf2).rzz;
        This->rxx=(*ptrtf1).rxx*(*ptrtf2).rxx+(*ptrtf1).rxy*(*ptrtf2).rxy+(*ptrtf1).rxz*(*ptrtf2).rxz;
        This->rxy=(*ptrtf1).rxx*(*ptrtf2).ryx+(*ptrtf1).rxy*(*ptrtf2).ryy+(*ptrtf1).rxz*(*ptrtf2).ryz;
        This->rxz=(*ptrtf1).rxx*(*ptrtf2).rzx+(*ptrtf1).rxy*(*ptrtf2).rzy+(*ptrtf1).rxz*(*ptrtf2).rzz;
        This->ryx=(*ptrtf1).ryx*(*ptrtf2).rxx+(*ptrtf1).ryy*(*ptrtf2).rxy+(*ptrtf1).ryz*(*ptrtf2).rxz;
        This->ryy=(*ptrtf1).ryx*(*ptrtf2).ryx+(*ptrtf1).ryy*(*ptrtf2).ryy+(*ptrtf1).ryz*(*ptrtf2).ryz;
        This->ryz=(*ptrtf1).ryx*(*ptrtf2).rzx+(*ptrtf1).ryy*(*ptrtf2).rzy+(*ptrtf1).ryz*(*ptrtf2).rzz;
        This->rzx=(*ptrtf1).rzx*(*ptrtf2).rxx+(*ptrtf1).rzy*(*ptrtf2).rxy+(*ptrtf1).rzz*(*ptrtf2).rxz;
        This->rzy=(*ptrtf1).rzx*(*ptrtf2).ryx+(*ptrtf1).rzy*(*ptrtf2).ryy+(*ptrtf1).rzz*(*ptrtf2).ryz;
        This->rzz=(*ptrtf1).rzx*(*ptrtf2).rzx+(*ptrtf1).rzy*(*ptrtf2).rzy+(*ptrtf1).rzz*(*ptrtf2).rzz;
        This->tx=(*ptrtf1).tx*(*ptrtf2).rxx+(*ptrtf1).ty*(*ptrtf2).rxy+(*ptrtf1).tz*(*ptrtf2).rxz+itf2tx;
        This->ty=(*ptrtf1).tx*(*ptrtf2).ryx+(*ptrtf1).ty*(*ptrtf2).ryy+(*ptrtf1).tz*(*ptrtf2).ryz+itf2ty;
        This->tz=(*ptrtf1).tx*(*ptrtf2).rzx+(*ptrtf1).ty*(*ptrtf2).rzy+(*ptrtf1).tz*(*ptrtf2).rzz+itf2tz;
        return *This;
}
__device__
G4ThreeVector G4AffineTransform_TransformPoint(const G4AffineTransform *This, G4ThreeVector vec)
{
        return G4ThreeVector_create(
   vec.x*This->rxx + vec.y*This->ryx + vec.z*This->rzx + This->tx,
   vec.x*This->rxy + vec.y*This->ryy + vec.z*This->rzy + This->ty,
   vec.x*This->rxz + vec.y*This->ryz + vec.z*This->rzz + This->tz );
}
__device__
G4ThreeVector G4AffineTransform_TransformAxis(const G4AffineTransform *This, G4ThreeVector axis)
{
  return G4ThreeVector_create(
   axis.x*This->rxx + axis.y*This->ryx + axis.z*This->rzx,
   axis.x*This->rxy + axis.y*This->ryy + axis.z*This->rzy,
   axis.x*This->rxz + axis.y*This->ryz + axis.z*This->rzz );
}
__device__
G4AffineTransform G4AffineTransform_Inverse(const G4AffineTransform *This)
{
        return G4AffineTransform_create_elements(
    This->rxx, This->ryx, This->rzx,
    This->rxy, This->ryy, This->rzy,
    This->rxz, This->ryz, This->rzz,
    -This->tx*This->rxx - This->ty*This->rxy - This->tz*This->rxz,
    -This->tx*This->ryx - This->ty*This->ryy - This->tz*This->ryz,
    -This->tx*This->rzx - This->ty*This->rzy - This->tz*This->rzz );
}
__device__
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
__device__
G4ThreeVector G4AffineTransform_NetTranslation(const G4AffineTransform *This)
{
        return G4ThreeVector_create(This->tx,This->ty,This->tz);
}
__device__
G4bool G4AffineTransform_IsRotated(const G4AffineTransform *This)
{
        return (This->rxx==1.0 && This->ryy==1.0 && This->rzz==1.0) ? false : true;
}
typedef struct
{
 G4double property;
}
StubMaterial;
struct G4SmartVoxelProxy;
typedef struct
{
 G4double fmaxExtent;
 G4double fminExtent;
 struct G4SmartVoxelProxy* * fslices;
 G4int fNumSlices;
 G4int fminEquivalent;
 G4int fmaxEquivalent;
 EAxis faxis;
 EAxis fparamAxis;
}
G4SmartVoxelHeader;
typedef struct
{
 G4int *fcontents;
 G4int fminEquivalent;
 G4int fmaxEquivalent;
 G4int fNumContents;
}
G4SmartVoxelNode;
typedef struct G4SmartVoxelProxy
{
 G4SmartVoxelHeader* fHeader;
    G4SmartVoxelNode* fNode;
}
G4SmartVoxelProxy;
__device__
void G4VoxelNode_ctor( G4SmartVoxelNode *This, G4int no )
{
 This->fmaxEquivalent = no;
 This->fminEquivalent = no;
 This->fcontents = 0;
 This->fNumContents = 0;
}
__device__ G4int
G4VoxelNode_GetNoContained( const G4SmartVoxelNode *This)
{
 return This->fNumContents;
}
__device__ G4int
G4VoxelNode_GetVolume(
 const G4SmartVoxelNode *This, G4int contentNo)
{
 (void)0;
 return This->fcontents[contentNo];
}
__device__ G4int
G4VoxelNode_GetMaxEquivalentSliceNo(
 const G4SmartVoxelNode *This )
{
 return This->fmaxEquivalent;
}
__device__ G4int
G4VoxelNode_GetMinEquivalentSliceNo(
 const G4SmartVoxelNode *This )
{
 return This->fminEquivalent;
}
__device__ G4int
G4VoxelHeader_GetMaxEquivalentSliceNo(
 const G4SmartVoxelHeader *This )
{
 return This->fmaxEquivalent;
}
__device__ G4int
G4VoxelHeader_GetMinEquivalentSliceNo(
 const G4SmartVoxelHeader *This )
{
 return This->fminEquivalent;
}
__device__ EAxis
G4VoxelHeader_GetAxis( const G4SmartVoxelHeader *This )
{
 return This->faxis;
}
__device__ G4int
G4VoxelHeader_GetNoSlices( const G4SmartVoxelHeader *This )
{
 return This->fNumSlices;
}
__device__ G4double
G4VoxelHeader_GetMinExtent( const G4SmartVoxelHeader *This )
{
 return This->fminExtent;
}
__device__ G4double
G4VoxelHeader_GetMaxExtent( const G4SmartVoxelHeader *This )
{
 return This->fmaxExtent;
}
__device__ G4SmartVoxelProxy*
G4VoxelHeader_GetSlice( const G4SmartVoxelHeader *This, G4int n )
{
 (void)0;
 return This->fslices[n];
}
__device__ G4bool
G4VoxelProxy_IsNode( const G4SmartVoxelProxy *This )
{
 return This->fNode != 0;
}
__device__ G4bool
G4VoxelProxy_IsHeader( const G4SmartVoxelProxy *This )
{
 return This->fHeader != 0;
}
__device__ G4SmartVoxelNode*
G4VoxelProxy_GetNode( const G4SmartVoxelProxy *This )
{
 return This->fNode;
}
__device__ G4SmartVoxelHeader*
G4VoxelProxy_GetHeader( const G4SmartVoxelProxy *This )
{
 return This->fHeader;
}
struct G4VPhysicalVolume;
struct G4VSolid;
typedef struct
{
 G4int fNoDaughters;
 struct G4VPhysicalVolume * *fDaughters;
   int check;
   StubMaterial* fMaterial;
 struct G4VSolid* fSolid;
 G4SmartVoxelHeader *fVoxel;
 int align;
}
G4LogicalVolume;
typedef struct G4VSolid
{
 ESolid type;
}
G4VSolid;
__device__
EInside G4VSolid_Inside( const G4VSolid *This, G4ThreeVector p);
__device__
G4ThreeVector G4VSolid_SurfaceNormal( const G4VSolid *This, G4ThreeVector p);
__device__
G4double G4VSolid_DistanceToIn_full(
    const G4VSolid *This,
    G4ThreeVector p,
    G4ThreeVector v);
__device__
G4double G4VSolid_DistanceToIn( const G4VSolid *This, G4ThreeVector p);
__device__
G4double G4VSolid_DistanceToOut_full(
      const G4VSolid *This,
      G4ThreeVector p,
      G4ThreeVector v,
      const G4bool calcNorm,
      G4bool *validNorm,
      G4ThreeVector *n);
__device__
G4double G4VSolid_DistanceToOut( const G4VSolid *This, G4ThreeVector p);
typedef struct
{
 G4VSolid solid;
    G4double fDx,fDy,fDz;
}
G4Box;
extern "C" {
__device__ EInside G4Box_Inside( const G4Box *This, G4ThreeVector p);
__device__ G4ThreeVector G4Box_SurfaceNormal( const G4Box *This, G4ThreeVector p);
__device__ G4double G4Box_DistanceToIn_full(
    const G4Box *This,
    G4ThreeVector p,
    G4ThreeVector v);
__device__ G4double G4Box_DistanceToIn( const G4Box *This, G4ThreeVector p);
__device__ G4double G4Box_DistanceToOut_full(
      const G4Box *This,
      G4ThreeVector p,
      G4ThreeVector v,
      const G4bool calcNorm,
      G4bool *validNorm,
      G4ThreeVector *n);
__device__ G4double G4Box_DistanceToOut( const G4Box *This, G4ThreeVector p);
__device__
G4ThreeVector G4Box_ApproxSurfaceNormal( const G4Box *This, G4ThreeVector p )
{
  G4double distx, disty, distz ;
  G4ThreeVector norm ;
  distx = fabs(fabs(p.x) - This->fDx) ;
  disty = fabs(fabs(p.y) - This->fDy) ;
  distz = fabs(fabs(p.z) - This->fDz) ;
  if ( distx <= disty )
  {
    if ( distx <= distz )
    {
      if ( p.x < 0 ) norm = G4ThreeVector_create(-1.0,0,0) ;
      else norm = G4ThreeVector_create( 1.0,0,0) ;
    }
    else
    {
      if ( p.z < 0 ) norm = G4ThreeVector_create(0,0,-1.0) ;
      else norm = G4ThreeVector_create(0,0, 1.0) ;
    }
  }
  else
  {
    if ( disty <= distz )
    {
      if ( p.y < 0 ) norm = G4ThreeVector_create(0,-1.0,0) ;
      else norm = G4ThreeVector_create(0, 1.0,0) ;
    }
    else
    {
      if ( p.z < 0 ) norm = G4ThreeVector_create(0,0,-1.0) ;
      else norm = G4ThreeVector_create(0,0, 1.0) ;
    }
  }
  return norm;
}
__device__
G4ThreeVector G4Box_SurfaceNormal( const G4Box *This, G4ThreeVector p)
{
  G4double distx, disty, distz ;
  G4ThreeVector norm ;
  const G4double kCarTolerance = 1E-3;
  distx = fabs(fabs(p.x) - This->fDx) ;
  disty = fabs(fabs(p.y) - This->fDy) ;
  distz = fabs(fabs(p.z) - This->fDz) ;
  const G4double delta = 0.5*kCarTolerance;
  const G4ThreeVector nX = G4ThreeVector_create( 1.0, 0,0 );
  const G4ThreeVector nmX = G4ThreeVector_create(-1.0, 0,0 );
  const G4ThreeVector nY = G4ThreeVector_create( 0, 1.0,0 );
  const G4ThreeVector nmY = G4ThreeVector_create( 0,-1.0,0 );
  const G4ThreeVector nZ = G4ThreeVector_create( 0, 0, 1.0);
  const G4ThreeVector nmZ = G4ThreeVector_create( 0, 0,- 1.0);
  G4ThreeVector
 normX = G4ThreeVector_create(0.,0.,0.),
 normY = G4ThreeVector_create(0.,0.,0.),
 normZ = G4ThreeVector_create(0.,0.,0.);
  G4ThreeVector sumnorm = G4ThreeVector_create(0., 0., 0.);
  G4int noSurfaces=0;
  if (distx <= delta)
  {
    noSurfaces ++;
    if ( p.x >= 0.){
      normX= nX ;
    }else{
      normX= nmX;
    }
    sumnorm= normX;
  }
  if (disty <= delta)
  {
    noSurfaces ++;
    if ( p.y >= 0.){
      normY= nY;
    }else{
      normY = nmY;
    }
    G4ThreeVector_sum_assign( &sumnorm, normY );
  }
  if (distz <= delta)
  {
    noSurfaces ++;
    if ( p.z >= 0.){
      normZ= nZ;
    }else{
      normZ = nmZ;
    }
    G4ThreeVector_sum_assign( &sumnorm, normZ );
  }
  const G4double invSqrt2 = 1.0 / sqrt( 2.0);
  const G4double invSqrt3 = 1.0 / sqrt( 3.0);
  norm= G4ThreeVector_create( 0., 0., 0.);
  if( noSurfaces > 0 )
  {
    if( noSurfaces == 1 ){
      norm= sumnorm;
    }else{
      if( noSurfaces == 2 ) {
        norm = G4ThreeVector_mult(sumnorm, invSqrt2);
      } else {
        norm = G4ThreeVector_mult(sumnorm, invSqrt3);
      }
    }
  }else{
     norm = G4Box_ApproxSurfaceNormal(This, p);
  }
  return norm;
}
__device__
G4double G4Box_DistanceToIn_full( const G4Box *This, G4ThreeVector p,G4ThreeVector v)
{
  G4double safx, safy, safz ;
  G4double smin=0.0, sminy, sminz ;
  G4double smax=kInfinity, smaxy, smaxz ;
  G4double stmp ;
  G4double sOut=kInfinity, sOuty=kInfinity, sOutz=kInfinity ;
  const G4double kCarTolerance = 1E-3;
  safx = fabs(p.x) - This->fDx ;
  safy = fabs(p.y) - This->fDy ;
  safz = fabs(p.z) - This->fDz ;
  if ( ((p.x*v.x >= 0.0) && safx > -kCarTolerance*0.5)
       || ((p.y*v.y >= 0.0) && safy > -kCarTolerance*0.5)
       || ((p.z*v.z >= 0.0) && safz > -kCarTolerance*0.5) )
  {
    return kInfinity ;
  }
  if ( v.x)
  {
    stmp = 1.0/fabs(v.x) ;
    if (safx >= 0.0)
    {
      smin = safx*stmp ;
      smax = (This->fDx+fabs(p.x))*stmp ;
    }
    else
    {
      if (v.x > 0) sOut = (This->fDx - p.x)*stmp ;
      if (v.x < 0) sOut = (This->fDx + p.x)*stmp ;
    }
  }
  if ( v.y)
  {
    stmp = 1.0/fabs(v.y) ;
    if (safy >= 0.0)
    {
      sminy = safy*stmp ;
      smaxy = (This->fDy+fabs(p.y))*stmp ;
      if (sminy > smin) smin=sminy ;
      if (smaxy < smax) smax=smaxy ;
      if (smin >= smax-kCarTolerance*0.5)
      {
        return kInfinity ;
      }
    }
    else
    {
      if (v.y > 0) sOuty = (This->fDy - p.y)*stmp ;
      if (v.y < 0) sOuty = (This->fDy + p.y)*stmp ;
      if( sOuty < sOut ) sOut = sOuty ;
    }
  }
  if ( v.z )
  {
    stmp = 1.0/fabs(v.z) ;
    if ( safz >= 0.0)
    {
      sminz = safz*stmp ;
      smaxz = (This->fDz+fabs(p.z))*stmp ;
      if (sminz > smin) smin = sminz ;
      if (smaxz < smax) smax = smaxz ;
      if (smin >= smax-kCarTolerance*0.5)
      {
        return kInfinity ;
      }
    }
    else
    {
      if (v.z > 0) sOutz = (This->fDz - p.z)*stmp ;
      if (v.z < 0) sOutz = (This->fDz + p.z)*stmp ;
      if( sOutz < sOut ) sOut = sOutz ;
    }
  }
  if ( sOut <= smin + 0.5*kCarTolerance)
  {
    return kInfinity ;
  }
  if (smin < 0.5*kCarTolerance) smin = 0.0 ;
  return smin ;
}
__device__
G4double G4Box_DistanceToIn( const G4Box *This, G4ThreeVector p)
{
  G4double safex, safey, safez, safe = 0.0 ;
  safex = fabs(p.x) - This->fDx ;
  safey = fabs(p.y) - This->fDy ;
  safez = fabs(p.z) - This->fDz ;
  if (safex > safe) safe = safex ;
  if (safey > safe) safe = safey ;
  if (safez > safe) safe = safez ;
  return safe ;
}
__device__
G4double G4Box_DistanceToOut_full( const G4Box *This, G4ThreeVector p,G4ThreeVector v,
                               const G4bool calcNorm,
                                G4bool *validNorm,G4ThreeVector *n)
{
  const G4double kCarTolerance = 1E-3;
  enum {kBoxUndefined,kPX,kMX,kPY,kMY,kPZ,kMZ} side = kBoxUndefined ;
  G4double pdist,stmp,snxt;
  if (calcNorm) *validNorm = true ;
  if (v.x > 0)
  {
    pdist = This->fDx - p.x ;
    if (pdist > kCarTolerance*0.5)
    {
      snxt = pdist/v.x ;
      side = kPX ;
    }
    else
    {
      if (calcNorm) *n = G4ThreeVector_create(1,0,0) ;
      return snxt = 0 ;
    }
  }
  else if (v.x < 0)
  {
    pdist = This->fDx + p.x ;
    if (pdist > kCarTolerance*0.5)
    {
      snxt = -pdist/v.x ;
      side = kMX ;
    }
    else
    {
      if (calcNorm) *n = G4ThreeVector_create(-1,0,0) ;
      return snxt = 0 ;
    }
  }
  else snxt = kInfinity ;
  if ( v.y > 0 )
  {
    pdist=This->fDy-p.y;
    if (pdist>kCarTolerance*0.5)
    {
      stmp=pdist/v.y;
      if (stmp<snxt)
      {
        snxt=stmp;
        side=kPY;
      }
    }
    else
    {
      if (calcNorm) *n = G4ThreeVector_create(0,1,0) ;
      return snxt = 0 ;
    }
  }
  else if ( v.y < 0 )
  {
    pdist = This->fDy + p.y ;
    if (pdist > kCarTolerance*0.5)
    {
      stmp=-pdist/v.y;
      if (stmp<snxt)
      {
        snxt=stmp;
        side=kMY;
      }
    }
    else
    {
      if (calcNorm) *n = G4ThreeVector_create(0,-1,0) ;
      return snxt = 0 ;
    }
  }
  if (v.z>0)
  {
    pdist=This->fDz-p.z;
    if (pdist > kCarTolerance*0.5)
    {
      stmp=pdist/v.z;
      if (stmp < snxt)
      {
        snxt=stmp;
        side=kPZ;
      }
    }
    else
    {
      if (calcNorm) *n = G4ThreeVector_create(0,0,1) ;
      return snxt = 0 ;
    }
  }
  else if (v.z<0)
  {
    pdist = This->fDz + p.z ;
    if (pdist > kCarTolerance*0.5)
    {
      stmp=-pdist/v.z;
      if (stmp < snxt)
      {
        snxt=stmp;
        side=kMZ;
      }
    }
    else
    {
      if (calcNorm) *n = G4ThreeVector_create(0,0,-1) ;
      return snxt = 0 ;
    }
  }
  if (calcNorm)
  {
    switch (side)
    {
      case kPX:
        *n=G4ThreeVector_create(1,0,0);
        break;
      case kMX:
        *n=G4ThreeVector_create(-1,0,0);
        break;
      case kPY:
        *n=G4ThreeVector_create(0,1,0);
        break;
      case kMY:
        *n=G4ThreeVector_create(0,-1,0);
        break;
      case kPZ:
        *n=G4ThreeVector_create(0,0,1);
        break;
      case kMZ:
        *n=G4ThreeVector_create(0,0,-1);
        break;
      default:
        break;
    }
  }
  return snxt;
}
__device__
G4double G4Box_DistanceToOut( const G4Box *This, G4ThreeVector p )
{
  G4double safx1,safx2,safy1,safy2,safz1,safz2,safe=0.0;
  safx1 = This->fDx - p.x ;
  safx2 = This->fDx + p.x ;
  safy1 = This->fDy - p.y ;
  safy2 = This->fDy + p.y ;
  safz1 = This->fDz - p.z ;
  safz2 = This->fDz + p.z ;
  if (safx2 < safx1) safe = safx2 ;
  else safe = safx1 ;
  if (safy1 < safe) safe = safy1 ;
  if (safy2 < safe) safe = safy2 ;
  if (safz1 < safe) safe = safz1 ;
  if (safz2 < safe) safe = safz2 ;
  if (safe < 0) safe = 0 ;
  return safe ;
}
__device__
EInside G4Box_Inside( const G4Box *This, G4ThreeVector p)
{
  const G4double kCarTolerance = 1E-3;
  EInside in = kOutside ;
  if ( fabs(p.x) <= This->fDx - kCarTolerance*0.5 )
  {
    if (fabs(p.y) <= This->fDy - kCarTolerance*0.5 )
    {
      if (fabs(p.z) <= This->fDz - kCarTolerance*0.5 ) in = kInside ;
      else if (fabs(p.z) <= This->fDz + kCarTolerance*0.5 ) in = kSurface ;
    }
    else if (fabs(p.y) <= This->fDy + kCarTolerance*0.5 )
    {
      if (fabs(p.z) <= This->fDz + kCarTolerance*0.5 ) in = kSurface ;
    }
  }
  else if (fabs(p.x) <= This->fDx + kCarTolerance*0.5 )
  {
    if (fabs(p.y) <= This->fDy + kCarTolerance*0.5 )
    {
      if (fabs(p.z) <= This->fDz + kCarTolerance*0.5) in = kSurface ;
    }
  }
  return in ;
}
}
typedef struct
{
 G4VSolid solid;
    G4double fRmax;
    G4double fRmaxTolerance;
 G4double align;
}
G4Orb;
extern "C" {
__device__ EInside G4Orb_Inside( const G4Orb *This, G4ThreeVector p);
__device__ G4ThreeVector G4Orb_SurfaceNormal( const G4Orb *This, G4ThreeVector p);
__device__ G4double G4Orb_DistanceToIn_full(
    const G4Orb *This,
    G4ThreeVector p,
    G4ThreeVector v);
__device__ G4double G4Orb_DistanceToIn( const G4Orb *This, G4ThreeVector p);
__device__ G4double G4Orb_DistanceToOut_full(
      const G4Orb *This,
      G4ThreeVector p,
      G4ThreeVector v,
      const G4bool calcNorm,
      G4bool *validNorm,
      G4ThreeVector *n);
__device__ G4double G4Orb_DistanceToOut( const G4Orb *This, G4ThreeVector p);
__device__
EInside G4Orb_Inside( const G4Orb *This, G4ThreeVector p)
{
  G4double rad2,tolRMax;
  EInside in;
  rad2 = G4ThreeVector_mag2(p);
  G4double rad = sqrt(rad2);
  tolRMax = This->fRmax - This->fRmaxTolerance*0.5 ;
  if ( rad <= tolRMax ) { in = kInside ; }
  else
  {
    tolRMax = This->fRmax + This->fRmaxTolerance*0.5 ;
    if ( rad <= tolRMax ) { in = kSurface ; }
    else { in = kOutside ; }
  }
  return in;
}
__device__
G4ThreeVector G4Orb_SurfaceNormal( const G4Orb *This, G4ThreeVector p)
{
  (void)This;
  return G4ThreeVector_unit(p);
}
__device__
G4double G4Orb_DistanceToIn_full( const G4Orb *This, G4ThreeVector p,G4ThreeVector v)
{
  G4double snxt = kInfinity ;
  G4double rad2, pDotV3d;
  G4double c, d2, s = kInfinity ;
  rad2 = G4ThreeVector_mag2(p);
  pDotV3d = G4ThreeVector_dot(p,v);
  G4double rad = sqrt(rad2);
  c = (rad - This->fRmax)*(rad + This->fRmax);
  if ( c > This->fRmaxTolerance*This->fRmax )
  {
    d2 = pDotV3d*pDotV3d - c ;
    if ( d2 >= 0 )
    {
      s = -pDotV3d - sqrt(d2) ;
      if ( s >= 0 )
      {
        return snxt = s;
      }
    }
    else
    {
      return snxt = kInfinity;
    }
  }
  else
  {
    if ( c > -This->fRmaxTolerance*This->fRmax )
    {
      d2 = pDotV3d*pDotV3d - c ;
      if ( (d2 < This->fRmaxTolerance*This->fRmax) || (pDotV3d >= 0) )
      {
        return snxt = kInfinity;
      }
      else
      {
        return snxt = 0.;
      }
    }
  }
  return snxt;
}
__device__
G4double G4Orb_DistanceToIn( const G4Orb *This, G4ThreeVector p)
{
  G4double safe = 0.0,
           rad = G4ThreeVector_mag(p);
  safe = rad - This->fRmax;
  if( safe < 0 ) { safe = 0.; }
  return safe;
}
__device__
G4double G4Orb_DistanceToOut_full( const G4Orb *This, G4ThreeVector p,G4ThreeVector v,
                               const G4bool calcNorm,
                                G4bool *validNorm,G4ThreeVector *n)
{
 G4double snxt = kInfinity;
  enum {kNull,kRMax} side = kNull;
  G4double rad2,pDotV3d;
  G4ThreeVector ipoint;
  G4double c,d2;
  rad2 = G4ThreeVector_mag2(p);
  pDotV3d = G4ThreeVector_dot(p,v);
  const G4double Rmax_plus = This->fRmax + This->fRmaxTolerance*0.5;
  G4double rad = sqrt(rad2);
  if ( rad <= Rmax_plus )
  {
    c = (rad - This->fRmax)*(rad + This->fRmax);
    if ( c < This->fRmaxTolerance*This->fRmax )
    {
      d2 = pDotV3d*pDotV3d - c;
      if( ( c > -This->fRmaxTolerance*This->fRmax) &&
          ( ( pDotV3d >= 0 ) || ( d2 < 0 )) )
      {
        if(calcNorm)
        {
          *validNorm = true ;
          *n = G4ThreeVector_create(p.x/This->fRmax,p.y/This->fRmax,p.z/This->fRmax) ;
        }
        return snxt = 0;
      }
      else
      {
        snxt = -pDotV3d + sqrt(d2);
        side = kRMax ;
      }
    }
  }
  else
  {
  }
  if (calcNorm)
  {
    switch( side )
    {
      case kRMax:
  ipoint = G4ThreeVector_saxpy(snxt,v,p);
  *n=G4ThreeVector_mult(ipoint,1.0/This->fRmax);
        *validNorm=true;
        break;
      default:
        break;
    }
  }
  return snxt;
}
__device__
G4double G4Orb_DistanceToOut( const G4Orb *This, G4ThreeVector p )
{
   G4double safe=0.0,rad = G4ThreeVector_mag(p);
  safe = This->fRmax - rad;
  if ( safe < 0. ) safe = 0.;
  return safe;
}
}
__device__
EInside G4VSolid_Inside( const G4VSolid *This, G4ThreeVector p)
{
 switch(This->type)
 {
  case kBox:
   return G4Box_Inside(( const G4Box*)This,p);
  case kOrb:
   return G4Orb_Inside(( const G4Orb*)This,p);
  default:
   (void)0;
   return kOutside;
 }
}
__device__
G4ThreeVector G4VSolid_SurfaceNormal( const G4VSolid *This, G4ThreeVector p)
{
 switch(This->type)
 {
  case kBox:
   return G4Box_SurfaceNormal(( const G4Box*)This,p);
  case kOrb:
   return G4Orb_SurfaceNormal(( const G4Orb*)This,p);
  default:
   (void)0;
   return G4ThreeVector_create(0,0,0);
 }
}
__device__
G4double G4VSolid_DistanceToIn_full(
    const G4VSolid *This,
    G4ThreeVector p,
    G4ThreeVector v)
{
 switch(This->type)
 {
  case kBox:
   return G4Box_DistanceToIn_full(( const G4Box*)This,p,v);
  case kOrb:
   return G4Orb_DistanceToIn_full(( const G4Orb*)This,p,v);
  default:
   (void)0;
   return 0;
 }
}
__device__
G4double G4VSolid_DistanceToIn( const G4VSolid *This, G4ThreeVector p)
{
 switch(This->type)
 {
  case kBox:
   return G4Box_DistanceToIn(( const G4Box*)This,p);
  case kOrb:
   return G4Orb_DistanceToIn(( const G4Orb*)This,p);
  default:
   (void)0;
   return 0;
 }
}
__device__
G4double G4VSolid_DistanceToOut_full(
      const G4VSolid *This,
      G4ThreeVector p,
      G4ThreeVector v,
      const G4bool calcNorm,
      G4bool *validNorm,
      G4ThreeVector *n)
{
 switch(This->type)
 {
  case kBox:
   return G4Box_DistanceToOut_full(( const G4Box*)This,p,v,calcNorm,validNorm,n);
  case kOrb:
   return G4Orb_DistanceToOut_full(( const G4Orb*)This,p,v,calcNorm,validNorm,n);
  default:
   (void)0;
   return 0;
 }
}
__device__
G4double G4VSolid_DistanceToOut( const G4VSolid *This, G4ThreeVector p)
{
 switch(This->type)
 {
  case kBox:
   return G4Box_DistanceToOut(( const G4Box*)This,p);
  case kOrb:
   return G4Orb_DistanceToOut(( const G4Orb*)This,p);
  default:
   (void)0;
   return 0;
 }
}
__device__
 G4SmartVoxelHeader * G4LogicalVolume_GetVoxelHeader( const G4LogicalVolume* This)
{
 return This->fVoxel;
}
__device__
G4int G4LogicalVolume_GetNoDaughters( const G4LogicalVolume* This)
{
  return This->fNoDaughters;
}
__device__
 struct G4VPhysicalVolume* G4LogicalVolume_GetDaughter( const G4LogicalVolume* This, const G4int i)
{
  return This->fDaughters[i];
}
__device__
 struct G4VSolid* G4LogicalVolume_GetSolid( const G4LogicalVolume* This)
{
  return This->fSolid;
}
__device__
 StubMaterial* G4LogicalVolume_GetMaterial( const G4LogicalVolume* This)
{
  return This->fMaterial;
}
typedef struct G4VPhysicalVolume
{
    G4RotationMatrix frot;
    G4ThreeVector ftrans;
 int guard1;
    G4LogicalVolume *flogical;
 int guard2;
 G4LogicalVolume *flmother;
 int guard3;
 int count;
 int counter_shadow;
}
G4VPhysicalVolume;
__device__
G4ThreeVector G4VPhysicalVolume_GetTranslation( const G4VPhysicalVolume *This)
{
  return This->ftrans;
}
__device__
 G4LogicalVolume* G4VPhysicalVolume_GetLogicalVolume( const G4VPhysicalVolume *This)
{
  return This->flogical;
}
__device__
 G4LogicalVolume* G4VPhysicalVolume_GetMotherLogical( const G4VPhysicalVolume *This)
{
  return This->flmother;
}
__device__
G4RotationMatrix G4VPhysicalVolume_GetObjectRotationValue( const G4VPhysicalVolume *This)
{
  return This->frot;
}
__device__
G4ThreeVector G4VPhysicalVolume_GetObjectTranslation( const G4VPhysicalVolume *This)
{
 return This->ftrans;
}
typedef struct
{
   G4AffineTransform fTransform;
   G4VPhysicalVolume* fPhysicalVolumePtr;
   EVolume fVolumeType;
}
G4NavigationLevel;
typedef struct
{
 G4NavigationLevel fNavHistory[16];
 G4int fStackDepth;
 int align;
}
G4NavigationHistory;
__device__
void G4NavigationLevel_ctor(
   G4NavigationLevel *This,
   G4VPhysicalVolume* pPhysVol,
   G4AffineTransform afTransform,
   EVolume volTp )
{
 This->fTransform = afTransform;
 This->fPhysicalVolumePtr = pPhysVol;
 This->fVolumeType = volTp;
}
__device__
void G4NavigationLevel_ctor_relative(
   G4NavigationLevel *This,
   G4VPhysicalVolume* pPhysVol,
   G4AffineTransform levelAbove,
   G4AffineTransform relativeCurrent,
   EVolume volTp )
{
 This->fPhysicalVolumePtr = pPhysVol;
 This->fVolumeType = volTp;
 G4AffineTransform_InverseProduct(&(This->fTransform), &levelAbove, &relativeCurrent );
}
__device__
G4NavigationLevel G4NavigationLevel_create(
   G4VPhysicalVolume* pPhysVol,
   G4AffineTransform afTransform,
   EVolume volTp )
{
 G4NavigationLevel lev;
 G4NavigationLevel_ctor( &lev, pPhysVol, afTransform, volTp );
 return lev;
}
__device__
G4NavigationLevel G4NavigationLevel_create_relative(
 G4VPhysicalVolume* pPhysVol,
 G4AffineTransform levelAbove,
 G4AffineTransform relativeCurrent,
 EVolume volTp)
{
 G4NavigationLevel lev;
 G4NavigationLevel_ctor_relative( &lev, pPhysVol, levelAbove, relativeCurrent, volTp );
 return lev;
}
__device__
 G4VPhysicalVolume* G4NavigationLevel_GetPhysicalVolume(
 const G4NavigationLevel *This )
{
  return This->fPhysicalVolumePtr;
}
__device__
G4AffineTransform G4NavigationLevel_GetTransform(
 const G4NavigationLevel *This )
{
  return This->fTransform;
}
__device__
const G4AffineTransform* G4NavigationLevel_GetPtrTransform(
 const G4NavigationLevel *This )
{
  return &(This->fTransform);
}
__device__
EVolume G4NavigationLevel_GetVolumeType(
 const G4NavigationLevel *This )
{
  return This->fVolumeType;
}
__device__
void G4NavigationHistory_Reset( G4NavigationHistory *This )
{
 This->fStackDepth = 0;
}
__device__
void G4NavigationHistory_Clear( G4NavigationHistory *This )
{
  G4AffineTransform origin = G4AffineTransform_create_vector(G4ThreeVector_create(0.,0.,0.));
  G4NavigationLevel tmpNavLevel = G4NavigationLevel_create(0, origin, kNormal) ;
  G4NavigationHistory_Reset( This );
  for (G4int ilev=16 -1; ilev>=0; ilev--)
  {
     This->fNavHistory[ilev] = tmpNavLevel;
  }
}
__device__
void G4NavigationHistory_ctor( G4NavigationHistory *This )
{
 This->fStackDepth = 0;
 G4NavigationHistory_Clear( This );
}
__device__
void G4NavigationHistory_dtor( G4NavigationHistory *This )
{
 (void)This;
}
__device__
void G4NavigationHistory_SetFirstEntry(
 G4NavigationHistory *This, G4VPhysicalVolume* pVol)
{
  G4ThreeVector translation = G4ThreeVector_create(0.,0.,0.);
  if( pVol!=0 )
  {
    translation = G4VPhysicalVolume_GetTranslation( pVol );
  }
  This->fNavHistory[0] =
    G4NavigationLevel_create( pVol, G4AffineTransform_create_vector(translation), kNormal );
}
__device__
const G4AffineTransform* G4NavigationHistory_GetPtrTopTransform(
 const G4NavigationHistory *This )
{
  return G4NavigationLevel_GetPtrTransform( &(This->fNavHistory[This->fStackDepth]) );
}
__device__
G4AffineTransform G4NavigationHistory_GetTopTransform(
 const G4NavigationHistory *This )
{
  return G4NavigationLevel_GetTransform( &(This->fNavHistory[This->fStackDepth]) );
}
__device__
EVolume G4NavigationHistory_GetTopVolumeType(
 const G4NavigationHistory *This )
{
  return G4NavigationLevel_GetVolumeType( &(This->fNavHistory[This->fStackDepth]) );
}
__device__
 G4VPhysicalVolume* G4NavigationHistory_GetTopVolume(
 const G4NavigationHistory *This )
{
  return G4NavigationLevel_GetPhysicalVolume( &(This->fNavHistory[This->fStackDepth]) );
}
__device__
G4int G4NavigationHistory_GetDepth(
 const G4NavigationHistory *This )
{
  return This->fStackDepth;
}
__device__
G4AffineTransform
G4NavigationHistory_GetTransform(
 const G4NavigationHistory *This, G4int n )
{
  return G4NavigationLevel_GetTransform( &(This->fNavHistory[n]) );
}
__device__
EVolume G4NavigationHistory_GetVolumeType(
 const G4NavigationHistory *This, G4int n )
{
  return G4NavigationLevel_GetVolumeType( &(This->fNavHistory[n]) );
}
__device__
 G4VPhysicalVolume* G4NavigationHistory_GetVolume(
 const G4NavigationHistory *This, G4int n )
{
  return G4NavigationLevel_GetPhysicalVolume( &(This->fNavHistory[n]) );
}
__device__
G4int G4NavigationHistory_GetMaxDepth(
 const G4NavigationHistory *This )
{
 (void)This;
 return 16;
}
__device__
void G4NavigationHistory_BackLevel( G4NavigationHistory *This )
{
  (void)0;
  This->fStackDepth--;
}
__device__
void G4NavigationHistory_NewLevel(
  G4NavigationHistory *This,
  G4VPhysicalVolume *pNewMother,
  EVolume vType )
{
  This->fStackDepth++;
  (void)0;
  This->fNavHistory[This->fStackDepth] =
    G4NavigationLevel_create_relative(
   pNewMother,
   G4NavigationLevel_GetTransform( &(This->fNavHistory[This->fStackDepth-1]) ),
   G4AffineTransform_create_full(
    G4VPhysicalVolume_GetObjectRotationValue( pNewMother ),
    G4VPhysicalVolume_GetTranslation( pNewMother )),
   vType );
}
typedef struct{
  G4VPhysicalVolume * PVolume;
  G4int trackId;
  }SolidInfo;
   typedef struct{
   float safety;
   float step;
   int trackId;
   G4VPhysicalVolume * PVolume;
   }ResultInfo;
   typedef struct{
   float safety;
   float step;
   G4VPhysicalVolume * PVolume;
   }
   FinalResult;
   typedef struct{
   G4ThreeVector Point;
   G4ThreeVector Direction;
   }PointInformation;
   __device__ void Find_minimum ( ResultInfo * Result_For_Current_Solid, FinalResult * Compacter_Result, int PrevSum, int size)
{
 int locationId = (blockIdx.x * blockDim.x + threadIdx.x);
 int i, loc ;
 float result_step, result_safety, Current_result_step, Current_result_safety;
 float Initial_result_step = (Compacter_Result [ locationId ]).step;
 float Initial_result_safety = (Compacter_Result [ locationId ]).safety;
 Current_result_step = Initial_result_step;
 Current_result_safety = Initial_result_safety;
 for( i = 0; i < size ; i++)
 {
  result_step = Result_For_Current_Solid[ PrevSum + i].step;
  result_safety = Result_For_Current_Solid[ PrevSum + i].safety;
  if ( result_step < Current_result_step)
  {
   loc = PrevSum + i;
   Current_result_step = result_step;
  }
  if ( result_safety < Current_result_safety)
  {
   Current_result_safety = result_safety;
  }
 }
 if( Current_result_step != Initial_result_step)
 {
  FinalResult final = { Current_result_safety, Current_result_step, (Result_For_Current_Solid[ loc ].PVolume)};
  Compacter_Result[ locationId ] = final;
 }
__syncthreads();
}
typedef struct
{
 G4double fVoxelSliceWidthStack[4];
 G4SmartVoxelHeader* fVoxelHeaderStack[4];
 G4int fVoxelNodeNoStack[4];
 G4int fVoxelNoSlicesStack[4];
 EAxis fVoxelAxisStack[4];
 G4int fVoxelDepth;
 G4SmartVoxelNode *fVoxelNode;
}
G4VoxelNavigation;
__device__ void G4VoxelNavigation_ctor( G4VoxelNavigation *This );
__device__ G4bool G4VoxelNavigation_LevelLocate(
 G4VoxelNavigation *This,
 G4NavigationHistory *history,
 const G4VPhysicalVolume *blockedVol,
 G4ThreeVector globalPoint,
 const G4ThreeVector* globalDirection,
 const G4bool pLocatedOnEdge,
 G4ThreeVector *localPoint );
__device__ G4SmartVoxelNode* G4VoxelNavigation_VoxelLocate(
 G4VoxelNavigation *This,
 G4SmartVoxelHeader *voxelHeader,
 G4ThreeVector point);
__device__
G4double
G4VoxelNavigation_ComputeStep(
   G4VoxelNavigation *This,
   G4ThreeVector localPoint,
   G4ThreeVector localDirection,
   const G4double currentProposedStepLength,
   G4double *newSafety,
   G4NavigationHistory *history,
   G4bool *validExitNormal,
   G4ThreeVector *exitNormal,
   G4bool *exiting,
   G4bool *entering,
   G4VPhysicalVolume *(*pBlockedPhysical)
   , G4double * Result
   );
__device__ G4double G4VoxelNavigation_ComputeSafety(
 G4VoxelNavigation *This,
 G4ThreeVector localPoint,
 const G4NavigationHistory *history);
typedef struct
{
 G4NavigationHistory fHistory;
 G4VoxelNavigation fVoxelNav;
 G4ThreeVector fStepEndPoint;
 G4ThreeVector fLastLocatedPointLocal;
 G4ThreeVector fExitNormal;
 G4ThreeVector fGrandMotherExitNormal;
 G4bool fEnteredDaughter;
 G4bool fExitedMother;
 G4bool fWasLimitedByGeometry;
 G4bool fEntering;
 G4bool fExiting;
 G4bool fLastStepWasZero;
 G4bool fLocatedOnEdge;
 G4bool fLocatedOutsideWorld;
 G4bool fValidExitNormal;
 G4bool fPushed;
 G4int fNumberZeroSteps;
 int align1;
 G4double fPreviousSafety;
 G4VPhysicalVolume *fBlockedPhysicalVolume;
 G4VPhysicalVolume *fTopPhysical;
}
G4Navigator;
__device__ void G4Navigator_ctor( G4Navigator *This );
__device__ void G4Navigator_SetWorldVolume(
 G4Navigator *This,
 G4VPhysicalVolume* pWorld );
__device__ G4VPhysicalVolume* G4Navigator_LocateGlobalPointAndSetup(
  G4Navigator *This,
  G4ThreeVector globalPoint,
  const G4ThreeVector* pGlobalDirection,
  G4bool relativeSearch,
  G4bool ignoreDirection,
  float * Result);
__device__
G4double G4Navigator_ComputeStep(
  G4Navigator *This,
  G4ThreeVector pGlobalpoint,
  G4ThreeVector pDirection,
  const G4double pCurrentProposedStepLength,
  G4double *pNewSafety
   , G4bool cur_vol_local
   , G4double * Result
  );
__device__ void G4Navigator_SetGeometricallyLimitedStep( G4Navigator *This );
__device__ G4double G4NormalNavigation_ComputeStep(
 G4ThreeVector localPoint,
 G4ThreeVector localDirection,
 const G4double currentProposedStepLength,
 G4double *newSafety,
 G4NavigationHistory *history,
 G4bool *validExitNormal,
 G4ThreeVector *exitNormal,
 G4bool *exiting,
 G4bool *entering,
 G4VPhysicalVolume *(*pBlockedPhysical));
__device__ G4double G4NormalNavigation_ComputeSafety(
 G4ThreeVector localPoint,
 const G4NavigationHistory *history );
__device__ G4bool G4NormalNavigation_LevelLocate(
 G4NavigationHistory *history,
 const G4VPhysicalVolume *blockedVol,
 G4ThreeVector* globalPoint,
 const G4ThreeVector* globalDirection,
 G4bool pLocatedOnEdge,
 G4ThreeVector* localPoint );
  __device__ void G4VoxelNavigation_ctor( G4VoxelNavigation *This );
__device__ void G4Navigator_ResetState( G4Navigator *This )
{
  This->fWasLimitedByGeometry = false;
  This->fEntering = false;
  This->fExiting = false;
  This->fLocatedOnEdge = false;
  This->fLastStepWasZero = false;
  This->fEnteredDaughter = false;
  This->fExitedMother = false;
  This->fPushed = false;
  This->fValidExitNormal = false;
  This->fExitNormal = G4ThreeVector_create(0,0,0);
  This->fPreviousSafety = 0.0;
  This->fNumberZeroSteps = 0;
  This->fBlockedPhysicalVolume = 0;
  This->fLastLocatedPointLocal = G4ThreeVector_create( 1e37, -1e37, 0.0 );
  This->fLocatedOutsideWorld = false;
}
__device__
G4ThreeVector G4Navigator_ComputeLocalAxis( const G4Navigator *This, G4ThreeVector pVec)
{
 G4AffineTransform t =
  G4NavigationHistory_GetTopTransform( &(This->fHistory) );
 return G4AffineTransform_TransformAxis(&t, pVec);
}
__device__ G4ThreeVector
G4Navigator_ComputeLocalPoint( const G4Navigator *This, G4ThreeVector pGlobalPoint)
{
 G4AffineTransform t =
  G4NavigationHistory_GetTopTransform( &(This->fHistory) );
 return G4AffineTransform_TransformPoint(&t, pGlobalPoint);
}
__device__ void G4Navigator_SetWorldVolume( G4Navigator *This, G4VPhysicalVolume* pWorld )
{
 This->fTopPhysical = pWorld;
 G4NavigationHistory_SetFirstEntry( &(This->fHistory), pWorld );
}
__device__ void G4Navigator_SetGeometricallyLimitedStep( G4Navigator *This )
{
 This->fWasLimitedByGeometry = true;
}
__device__
void G4Navigator_ResetStackAndState( G4Navigator *This )
{
 G4NavigationHistory_Reset( &(This->fHistory) );
 G4Navigator_ResetState( This );
}
__device__
EVolume G4Navigator_VolumeType( const G4Navigator *This, const G4VPhysicalVolume *pVol )
{
 (void)This;
 (void)pVol;
 return kNormal;
}
__device__ void G4Navigator_ctor( G4Navigator *This )
{
 G4NavigationHistory_ctor( &(This->fHistory) );
 G4VoxelNavigation_ctor( &(This->fVoxelNav ) );
 G4Navigator_ResetStackAndState( This );
 This->fWasLimitedByGeometry = false;
 This->fTopPhysical = 0;
 This->fPushed = false;
 This->fStepEndPoint = G4ThreeVector_create( kInfinity, kInfinity, kInfinity );
}
__device__
 G4VPhysicalVolume*
G4Navigator_LocateGlobalPointAndSetup(
  G4Navigator *This,
  G4ThreeVector globalPoint,
  const G4ThreeVector* pGlobalDirection,
  G4bool relativeSearch,
  G4bool ignoreDirection,
  float * Result
  )
{
  G4bool notKnownContained=true, noResult;
  G4VPhysicalVolume *targetPhysical;
  G4VSolid *targetSolid = 0;
  G4ThreeVector localPoint = G4ThreeVector_create(0,0,0);
  G4ThreeVector globalDirection = G4ThreeVector_create(0,0,0);
  EInside insideCode;
  G4bool considerDirection = (!ignoreDirection) || This->fLocatedOnEdge;
  if( considerDirection && pGlobalDirection != 0 )
  {
    globalDirection=*pGlobalDirection;
  }
  if ( 1 )
  {
     G4Navigator_ResetStackAndState( This );
  }
  else
  {
    if ( This->fWasLimitedByGeometry )
    {
      This->fWasLimitedByGeometry = false;
      This->fEnteredDaughter = This->fEntering;
      This->fExitedMother = This->fExiting;
      if ( This->fExiting )
      {
        if ( G4NavigationHistory_GetDepth( &(This->fHistory) ) )
        {
          This->fBlockedPhysicalVolume = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
          G4NavigationHistory_BackLevel( &(This->fHistory) );
        }
        else
        {
          This->fLastLocatedPointLocal = localPoint;
          This->fLocatedOutsideWorld = true;
          return 0;
        }
        if ( This->fLocatedOnEdge )
        {
          This->fExiting= false;
        }
      }
      else
        if ( This->fEntering )
        {
    G4NavigationHistory_NewLevel( &(This->fHistory), This->fBlockedPhysicalVolume, kNormal);
          This->fEntering = false;
          This->fBlockedPhysicalVolume = 0;
          G4AffineTransform t = G4NavigationHistory_GetTopTransform( &(This->fHistory) );
          localPoint = G4AffineTransform_TransformPoint(&t,globalPoint);
          notKnownContained = false;
        }
    }
    else
    {
      This->fBlockedPhysicalVolume = 0;
      This->fEntering = false;
      This->fEnteredDaughter = false;
      This->fExiting = false;
      This->fExitedMother = false;
    }
  }
  while (notKnownContained)
  {
 targetSolid =
   G4LogicalVolume_GetSolid(
    G4VPhysicalVolume_GetLogicalVolume(
    G4NavigationHistory_GetTopVolume(&(This->fHistory))));
 G4AffineTransform t = G4NavigationHistory_GetTopTransform( &(This->fHistory) );
 localPoint = G4AffineTransform_TransformPoint(&t,globalPoint);
 insideCode = G4VSolid_Inside(targetSolid,localPoint);
    if ( insideCode==kOutside )
    {
      if ( G4NavigationHistory_GetDepth( &(This->fHistory) ) )
      {
        This->fBlockedPhysicalVolume = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
        G4NavigationHistory_BackLevel( &(This->fHistory) );
        This->fExiting = false;
      }
      else
      {
        This->fLastLocatedPointLocal = localPoint;
        This->fLocatedOutsideWorld = true;
        return 0;
      }
    }
    else
      if ( insideCode==kSurface )
      {
        G4bool isExiting = This->fExiting;
        if( (!This->fExiting)&&considerDirection )
        {
   G4bool directionExiting = false;
   G4AffineTransform t = G4NavigationHistory_GetTopTransform( &(This->fHistory) );
   G4ThreeVector localDirection =G4AffineTransform_TransformAxis(&t,globalDirection);
   G4ThreeVector normal = G4VSolid_SurfaceNormal(targetSolid, localPoint);
   directionExiting = G4ThreeVector_dot(normal,localDirection) > 0.0;
   isExiting = isExiting || directionExiting;
        }
        if( isExiting )
        {
          if ( G4NavigationHistory_GetDepth( &(This->fHistory) ) )
          {
            This->fBlockedPhysicalVolume = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
            G4NavigationHistory_BackLevel( &(This->fHistory) );
            This->fValidExitNormal = false;
          }
          else
          {
            This->fLastLocatedPointLocal = localPoint;
            This->fLocatedOutsideWorld = true;
            return 0;
          }
        }
        else
        {
          notKnownContained=false;
        }
      }
      else
      {
        notKnownContained=false;
      }
  }
  noResult = true;
  do
  {
    targetPhysical = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
  G4LogicalVolume *targetLogical = G4VPhysicalVolume_GetLogicalVolume(targetPhysical);
    if ( G4LogicalVolume_GetVoxelHeader( targetLogical ) != 0 )
    {
  noResult =
   G4VoxelNavigation_LevelLocate(
    &(This->fVoxelNav),
    &(This->fHistory),
    This->fBlockedPhysicalVolume,
    globalPoint,
    pGlobalDirection,
    considerDirection,
    &localPoint);
 }
 else
 {
  noResult = G4NormalNavigation_LevelLocate(
    &(This->fHistory),
    This->fBlockedPhysicalVolume,
    &globalPoint,
    pGlobalDirection,
    considerDirection,
    &localPoint);
 }
    if ( noResult )
    {
      This->fBlockedPhysicalVolume = 0;
      This->fEntering = false;
      This->fEnteredDaughter = true;
    }
  } while (noResult);
  This->fLastLocatedPointLocal = localPoint;
  This->fLocatedOutsideWorld= false;
  return targetPhysical;
}
__device__ void
G4Navigator_LocateGlobalPointWithinVolume( G4Navigator *This, G4ThreeVector pGlobalpoint)
{
 This->fLastLocatedPointLocal = G4Navigator_ComputeLocalPoint( This, pGlobalpoint );
 G4VPhysicalVolume* motherPhysical = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
 G4LogicalVolume* motherLogical = G4VPhysicalVolume_GetLogicalVolume( motherPhysical );
 G4SmartVoxelHeader* pVoxelHeader = G4LogicalVolume_GetVoxelHeader( motherLogical );
 if ( pVoxelHeader )
 {
  G4VoxelNavigation_VoxelLocate( &(This->fVoxelNav), pVoxelHeader, This->fLastLocatedPointLocal );
 }
 This->fBlockedPhysicalVolume = 0;
 This->fEntering = false;
 This->fEnteredDaughter = false;
 This->fExiting = false;
 This->fExitedMother = false;
}
__device__
G4double G4Navigator_ComputeStep(
  G4Navigator *This,
  G4ThreeVector pGlobalpoint,
  G4ThreeVector pDirection,
  const G4double pCurrentProposedStepLength,
  G4double *pNewSafety
   , G4bool cur_vol_local
   , G4double * Result
  )
{
  G4ThreeVector localDirection = G4Navigator_ComputeLocalAxis(This,pDirection);
  G4double Step = 1e37;
  G4VPhysicalVolume *motherPhysical = G4NavigationHistory_GetTopVolume( &(This->fHistory) );
  const G4double kCarTolerance = 1E-3;
  G4LogicalVolume *motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  G4ThreeVector newLocalPoint = G4Navigator_ComputeLocalPoint( This, pGlobalpoint);
  if( !G4ThreeVector_equal(newLocalPoint, This->fLastLocatedPointLocal) )
  {
    G4ThreeVector oldLocalPoint = This->fLastLocatedPointLocal;
    G4double moveLenSq = G4ThreeVector_diff2(newLocalPoint,oldLocalPoint);
    if ( moveLenSq >= kCarTolerance*kCarTolerance )
    {
      G4Navigator_LocateGlobalPointWithinVolume( This, pGlobalpoint );
    }
  }
  if ( G4LogicalVolume_GetVoxelHeader(motherLogical) != 0 )
  {
 if( cur_vol_local )
 Step = G4VoxelNavigation_ComputeStep(
   &(This->fVoxelNav),
   This->fLastLocatedPointLocal,
   localDirection,
   pCurrentProposedStepLength,
   pNewSafety,
   &(This->fHistory),
   &(This->fValidExitNormal),
   &(This->fExitNormal),
   &(This->fExiting),
   &(This->fEntering),
   &(This->fBlockedPhysicalVolume)
   , Result
   );
 else
     return 0;
  }
  else
  {
 Step = G4NormalNavigation_ComputeStep(
   This->fLastLocatedPointLocal,
   localDirection,
   pCurrentProposedStepLength,
   pNewSafety,
   &(This->fHistory),
   &(This->fValidExitNormal),
   &(This->fExitNormal),
   &(This->fExiting),
   &(This->fEntering),
   &(This->fBlockedPhysicalVolume));
  }
  This->fPreviousSafety = *pNewSafety;
  This->fLocatedOnEdge = This->fLastStepWasZero && (Step==0.0);
  This->fLastStepWasZero = (Step==0.0);
  if (This->fPushed) This->fPushed = This->fLastStepWasZero;
  if ( This->fLastStepWasZero )
  {
    This->fNumberZeroSteps++;
    if( This->fNumberZeroSteps > 10 -1 )
    {
       Step += 0.9*kCarTolerance;
       This->fPushed = true;
    }
    if( This->fNumberZeroSteps > 25 -1 )
    {
   (void)0;
    }
  }
  else
  {
    if (!This->fPushed) This->fNumberZeroSteps = 0;
  }
  This->fEnteredDaughter = This->fEntering;
  This->fExitedMother = This->fExiting;
  if( This->fExiting )
  {
    if(This->fValidExitNormal)
    {
      This->fGrandMotherExitNormal= This->fExitNormal;
    }
    else
    {
      G4ThreeVector finalLocalPoint =
  G4ThreeVector_saxpy( Step, localDirection, This->fLastLocatedPointLocal );
      This->fGrandMotherExitNormal =
  G4VSolid_SurfaceNormal(
   G4LogicalVolume_GetSolid(motherLogical),finalLocalPoint);
      G4RotationMatrix mRot = G4VPhysicalVolume_GetObjectRotationValue(motherPhysical);
      G4RotationMatrix inv = G4RotationMatrix_inverse(&mRot);
      This->fGrandMotherExitNormal
       = G4RotationMatrix_apply(&inv,This->fGrandMotherExitNormal);
    }
  }
  This->fStepEndPoint =
 G4ThreeVector_saxpy(Step, pDirection, pGlobalpoint );
  if( (Step == pCurrentProposedStepLength) && (!This->fExiting) && (!This->fEntering) )
  {
    Step = kInfinity;
  }
  return Step;
}
__device__ G4bool
G4AuxiliaryNavServices_CheckPointOnSurface(
         const G4VSolid* sampleSolid,
                     G4ThreeVector localPoint,
                     const G4ThreeVector* globalDirection,
                     G4AffineTransform sampleTransform,
                     const G4bool locatedOnEdge)
{
  G4ThreeVector localDirection, sampleNormal;
  G4bool enter = false;
  EInside insideSolid =
 G4VSolid_Inside(sampleSolid, localPoint);
  if ( insideSolid!=kOutside )
  {
    G4bool checkDirection= locatedOnEdge && (globalDirection!=0);
    if( (insideSolid==kSurface) && checkDirection)
    {
      localDirection= G4AffineTransform_TransformAxis(&sampleTransform,*globalDirection);
      sampleNormal = G4VSolid_SurfaceNormal(sampleSolid,localPoint);
      if ( G4ThreeVector_dot(sampleNormal,localDirection) <= 0 )
      {
        if( G4ThreeVector_dot(sampleNormal,localDirection) == 0 )
        {
          G4double distanceToIn =
   G4VSolid_DistanceToIn_full( sampleSolid, localPoint, localDirection );
          if( distanceToIn != kInfinity )
          {
            enter = true;
          }
        }
        else
        {
          enter = true;
        }
      }
    }
    else
    {
      enter = true;
    }
  }
  return enter;
}
__device__ G4bool
G4NormalNavigation_LevelLocate(
 G4NavigationHistory *history,
 const G4VPhysicalVolume *blockedVol,
 G4ThreeVector* globalPoint,
 const G4ThreeVector* globalDirection,
 G4bool pLocatedOnEdge,
 G4ThreeVector* localPoint )
{
  G4VPhysicalVolume *targetPhysical, *samplePhysical;
  G4LogicalVolume *targetLogical;
  G4VSolid *sampleSolid;
  G4ThreeVector samplePoint;
  G4int targetNoDaughters;
  targetPhysical = G4NavigationHistory_GetTopVolume(history);
  targetLogical = G4VPhysicalVolume_GetLogicalVolume(targetPhysical);
  targetNoDaughters = G4LogicalVolume_GetNoDaughters(targetLogical);
  if (targetNoDaughters == 0) return false;
  for ( int sampleNo=targetNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
   samplePhysical =
  G4LogicalVolume_GetDaughter(targetLogical,sampleNo);
   if ( samplePhysical!=blockedVol )
   {
  G4NavigationHistory_NewLevel(history, samplePhysical, kNormal );
  sampleSolid =
   G4LogicalVolume_GetSolid(
    G4VPhysicalVolume_GetLogicalVolume(samplePhysical));
  G4AffineTransform tf =
   G4NavigationHistory_GetTopTransform(history);
  samplePoint =
   G4AffineTransform_TransformPoint( &tf, *globalPoint );
  if( G4AuxiliaryNavServices_CheckPointOnSurface(
   sampleSolid, samplePoint, globalDirection,
   tf, pLocatedOnEdge) )
  {
    *localPoint = samplePoint;
    return true;
  }
  else
  {
   G4NavigationHistory_BackLevel(history);
  }
   }
  }
  return false;
}
__device__
G4double
G4NormalNavigation_ComputeStep(
 G4ThreeVector localPoint,
 G4ThreeVector localDirection,
 const G4double currentProposedStepLength,
 G4double *newSafety,
 G4NavigationHistory *history,
 G4bool *validExitNormal,
 G4ThreeVector *exitNormal,
 G4bool *exiting,
 G4bool *entering,
 G4VPhysicalVolume *(*pBlockedPhysical))
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical, *blockedExitedVol=0;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int localNoDaughters, sampleNo;
  motherPhysical = G4NavigationHistory_GetTopVolume(history);
  motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = G4LogicalVolume_GetSolid(motherLogical);
  motherSafety = G4VSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety;
  if ( *exiting && *validExitNormal )
  {
    if ( G4ThreeVector_dot(localDirection,*exitNormal)>=kMinExitingNormalCosine )
    {
      blockedExitedVol =* pBlockedPhysical;
      ourSafety = 0;
    }
  }
  *exiting = false;
  *entering = false;
  localNoDaughters = G4LogicalVolume_GetNoDaughters(motherLogical);
  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo--)
  {
    samplePhysical = G4LogicalVolume_GetDaughter(motherLogical,sampleNo);
    if ( samplePhysical!=blockedExitedVol )
    {
      G4AffineTransform sampleTf =
       G4AffineTransform_create_full(
  G4VPhysicalVolume_GetObjectRotationValue(samplePhysical),
  G4VPhysicalVolume_GetTranslation(samplePhysical));
   G4AffineTransform_Invert(&sampleTf);
      const G4ThreeVector samplePoint =
   G4AffineTransform_TransformPoint(&sampleTf, localPoint);
      const G4VSolid *sampleSolid =
  G4LogicalVolume_GetSolid(
   G4VPhysicalVolume_GetLogicalVolume( samplePhysical ));
      const G4double sampleSafety =
  G4VSolid_DistanceToIn(sampleSolid,samplePoint);
      if ( sampleSafety<ourSafety )
      {
        ourSafety=sampleSafety;
      }
      if ( sampleSafety<=ourStep )
      {
        sampleDirection = G4AffineTransform_TransformAxis(&sampleTf, localDirection);
        const G4double sampleStep =
   G4VSolid_DistanceToIn_full(sampleSolid,samplePoint,sampleDirection);
        if ( sampleStep<=ourStep )
        {
          ourStep = sampleStep;
          *entering = true;
          *exiting = false;
          *pBlockedPhysical = samplePhysical;
        }
      }
    }
  }
  if ( currentProposedStepLength<ourSafety )
  {
    *entering = false;
    *exiting = false;
    *pBlockedPhysical = 0;
    ourStep = kInfinity;
  }
  else
  {
    if ( motherSafety<=ourStep )
    {
      G4double motherStep =
  G4VSolid_DistanceToOut_full(
   motherSolid,
   localPoint,
   localDirection,
   true,
   validExitNormal,
   exitNormal);
      if ( motherStep<=ourStep )
      {
        ourStep = motherStep;
        *exiting = true;
        *entering = false;
        if ( *validExitNormal )
        {
          G4RotationMatrix rot = G4VPhysicalVolume_GetObjectRotationValue(motherPhysical);
    G4RotationMatrix inv = G4RotationMatrix_inverse(&rot);
          *exitNormal = G4RotationMatrix_apply(&inv, *exitNormal);
        }
      }
      else
      {
        *validExitNormal = false;
      }
    }
  }
  *newSafety = ourSafety;
  return ourStep;
}
__device__ G4bool
G4AuxiliaryNavServices_CheckPointOnSurface(
         const G4VSolid* sampleSolid,
                     G4ThreeVector localPoint,
                     const G4ThreeVector* globalDirection,
                     G4AffineTransform sampleTransform,
                     const G4bool locatedOnEdge);
__device__ G4bool
G4AuxiliaryNavServices_CheckPointExiting(
       const G4VSolid* sampleSolid,
                   G4ThreeVector localPoint,
                   const G4ThreeVector* globalDirection,
                   G4AffineTransform sampleTransform );
__device__
 G4SmartVoxelNode*
G4VoxelNavigation_VoxelLocate(
   G4VoxelNavigation *This,
   G4SmartVoxelHeader* pHead,
   G4ThreeVector localPoint )
{
  G4SmartVoxelHeader *targetVoxelHeader=pHead;
  G4SmartVoxelNode *targetVoxelNode = 0;
  const G4SmartVoxelProxy *sampleProxy;
  EAxis targetHeaderAxis;
  G4double targetHeaderMin, targetHeaderNodeWidth;
  G4int targetHeaderNoSlices, targetNodeNo;
  This->fVoxelDepth = 0;
  while ( targetVoxelNode == 0 )
  {
    targetHeaderAxis = G4VoxelHeader_GetAxis(targetVoxelHeader);
    targetHeaderNoSlices = G4VoxelHeader_GetNoSlices(targetVoxelHeader);
    targetHeaderMin = G4VoxelHeader_GetMinExtent(targetVoxelHeader);
    targetHeaderNodeWidth =
  (G4VoxelHeader_GetMaxExtent(targetVoxelHeader)-targetHeaderMin)
                          / targetHeaderNoSlices;
    targetNodeNo = (G4int)(
  (G4ThreeVector_coord(localPoint,targetHeaderAxis)-targetHeaderMin)
                          / targetHeaderNodeWidth);
    if ( targetNodeNo<0 )
    {
  targetNodeNo = 0;
    }
    else if ( targetNodeNo>=targetHeaderNoSlices )
 {
  targetNodeNo = targetHeaderNoSlices-1;
 }
    This->fVoxelAxisStack[This->fVoxelDepth] = targetHeaderAxis;
    This->fVoxelNoSlicesStack[This->fVoxelDepth] = targetHeaderNoSlices;
    This->fVoxelSliceWidthStack[This->fVoxelDepth] = targetHeaderNodeWidth;
    This->fVoxelNodeNoStack[This->fVoxelDepth] = targetNodeNo;
    This->fVoxelHeaderStack[This->fVoxelDepth] = targetVoxelHeader;
    sampleProxy = G4VoxelHeader_GetSlice(targetVoxelHeader, targetNodeNo);
    if ( G4VoxelProxy_IsNode(sampleProxy) )
    {
      targetVoxelNode = G4VoxelProxy_GetNode(sampleProxy);
    }
    else
    {
      targetVoxelHeader = G4VoxelProxy_GetHeader(sampleProxy);
      This->fVoxelDepth++;
      (void)0;
    }
  }
  This->fVoxelNode = targetVoxelNode;
  return targetVoxelNode;
}
__device__
G4bool
G4VoxelNavigation_LocateNextVoxel(
   G4VoxelNavigation *This,
   G4ThreeVector localPoint,
   G4ThreeVector localDirection,
   const G4double currentStep )
{
  G4SmartVoxelHeader *workHeader=0, *newHeader=0;
  G4SmartVoxelProxy *newProxy=0;
  G4SmartVoxelNode *newVoxelNode= 0;
  G4ThreeVector targetPoint, voxelPoint;
  G4double workNodeWidth, workMinExtent, workCoord;
  G4double minVal, maxVal, newDistance=0.;
  G4double newHeaderMin, newHeaderNodeWidth;
  G4int depth=0, newDepth=0, workNodeNo=0, newNodeNo=0, newHeaderNoSlices=0;
  EAxis workHeaderAxis, newHeaderAxis;
  G4bool isNewVoxel=false;
  G4double currentDistance = currentStep;
  for (depth=0; depth<This->fVoxelDepth; depth++)
  {
    targetPoint =
  G4ThreeVector_saxpy(currentDistance,localDirection,localPoint);
    newDistance = currentDistance;
    workHeader = This->fVoxelHeaderStack[depth];
    workHeaderAxis = This->fVoxelAxisStack[depth];
    workNodeNo = This->fVoxelNodeNoStack[depth];
    workNodeWidth = This->fVoxelSliceWidthStack[depth];
    workMinExtent = G4VoxelHeader_GetMinExtent(workHeader);
    workCoord = G4ThreeVector_coord(targetPoint,workHeaderAxis);
    minVal = workMinExtent+workNodeNo*workNodeWidth;
    if ( minVal<=workCoord+1E-3*0.5 )
    {
      maxVal = minVal+workNodeWidth;
      if ( maxVal<=workCoord-1E-3*0.5 )
      {
        newNodeNo = workNodeNo+1;
        newHeader = workHeader;
        newDistance = (maxVal-G4ThreeVector_coord(localPoint,workHeaderAxis))
                    / G4ThreeVector_coord(localDirection,workHeaderAxis);
        isNewVoxel = true;
        newDepth = depth;
      }
    }
    else
    {
      newNodeNo = workNodeNo-1;
      newHeader = workHeader;
      newDistance = (minVal-G4ThreeVector_coord(localPoint,workHeaderAxis))
                  / G4ThreeVector_coord(localDirection,workHeaderAxis);
      isNewVoxel = true;
      newDepth = depth;
    }
    currentDistance = newDistance;
  }
  targetPoint =
 G4ThreeVector_saxpy(currentDistance,localDirection,localPoint);
  depth = This->fVoxelDepth;
  {
    workHeader = This->fVoxelHeaderStack[depth];
    workHeaderAxis = This->fVoxelAxisStack[depth];
    workNodeNo = This->fVoxelNodeNoStack[depth];
    workNodeWidth = This->fVoxelSliceWidthStack[depth];
    workMinExtent = G4VoxelHeader_GetMinExtent(workHeader);
    workCoord = G4ThreeVector_coord(targetPoint,workHeaderAxis);
    minVal = workMinExtent+G4VoxelNode_GetMinEquivalentSliceNo(This->fVoxelNode)*workNodeWidth;
    if ( minVal<=workCoord+1E-3*0.5 )
    {
      maxVal = workMinExtent+(G4VoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)+1)
                            *workNodeWidth;
      if ( maxVal<=workCoord-1E-3*0.5 )
      {
        newNodeNo = G4VoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)+1;
        newHeader = workHeader;
        newDistance = (maxVal-G4ThreeVector_coord(localPoint,workHeaderAxis))
                    / G4ThreeVector_coord(localDirection,workHeaderAxis);
        isNewVoxel = true;
        newDepth = depth;
      }
    }
    else
    {
      newNodeNo = G4VoxelNode_GetMinEquivalentSliceNo(This->fVoxelNode)-1;
      newHeader = workHeader;
      newDistance = (minVal-G4ThreeVector_coord(localPoint,workHeaderAxis))
                  / G4ThreeVector_coord(localDirection,workHeaderAxis);
      isNewVoxel = true;
      newDepth = depth;
    }
    currentDistance = newDistance;
  }
  if (isNewVoxel)
  {
    if ( (newNodeNo<0) || (newNodeNo>=G4VoxelHeader_GetNoSlices(newHeader)))
    {
      isNewVoxel = false;
    }
    else
    {
      voxelPoint = G4ThreeVector_saxpy(newDistance,localDirection,localPoint);
      (void)0;
      This->fVoxelNodeNoStack[newDepth] = newNodeNo;
      This->fVoxelDepth = newDepth;
      newVoxelNode = 0;
      while ( newVoxelNode == 0 )
      {
        newProxy = G4VoxelHeader_GetSlice(newHeader,newNodeNo);
        if ( G4VoxelProxy_IsNode(newProxy) )
        {
          newVoxelNode = G4VoxelProxy_GetNode(newProxy);
        }
        else
        {
          This->fVoxelDepth++;
          (void)0;
          newHeader = G4VoxelProxy_GetHeader(newProxy);
          newHeaderAxis = G4VoxelHeader_GetAxis(newHeader);
          newHeaderNoSlices = G4VoxelHeader_GetNoSlices(newHeader);
          newHeaderMin = G4VoxelHeader_GetMinExtent(newHeader);
          newHeaderNodeWidth =
   (G4VoxelHeader_GetMaxExtent(newHeader)-newHeaderMin)
                             / newHeaderNoSlices;
          newNodeNo = (G4int)(
   (G4ThreeVector_coord(voxelPoint,newHeaderAxis)-newHeaderMin)
                             / newHeaderNodeWidth );
          if ( newNodeNo<0 )
          {
            newNodeNo=0;
          }
          else if ( newNodeNo>=newHeaderNoSlices )
               {
                 newNodeNo = newHeaderNoSlices-1;
               }
          This->fVoxelAxisStack[This->fVoxelDepth] = newHeaderAxis;
          This->fVoxelNoSlicesStack[This->fVoxelDepth] = newHeaderNoSlices;
          This->fVoxelSliceWidthStack[This->fVoxelDepth] = newHeaderNodeWidth;
          This->fVoxelNodeNoStack[This->fVoxelDepth] = newNodeNo;
          This->fVoxelHeaderStack[This->fVoxelDepth] = newHeader;
        }
      }
      This->fVoxelNode = newVoxelNode;
    }
  }
  return isNewVoxel;
}
__device__
G4double
G4VoxelNavigation_ComputeVoxelSafety(
   const G4VoxelNavigation *This,
   G4ThreeVector localPoint)
{
  G4SmartVoxelHeader *curHeader;
  G4double voxelSafety, curNodeWidth;
  G4double curNodeOffset, minCurCommonDelta, maxCurCommonDelta;
  G4int minCurNodeNoDelta, maxCurNodeNoDelta;
  G4int localVoxelDepth, curNodeNo;
  EAxis curHeaderAxis;
  localVoxelDepth = This->fVoxelDepth;
  curHeader = This->fVoxelHeaderStack[localVoxelDepth];
  curHeaderAxis = This->fVoxelAxisStack[localVoxelDepth];
  curNodeNo = This->fVoxelNodeNoStack[localVoxelDepth];
  curNodeWidth = This->fVoxelSliceWidthStack[localVoxelDepth];
  curNodeOffset = curNodeNo*curNodeWidth;
  maxCurNodeNoDelta = G4VoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)-curNodeNo;
  minCurNodeNoDelta = curNodeNo-G4VoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode);
  minCurCommonDelta = G4ThreeVector_coord(localPoint,curHeaderAxis)
 - G4VoxelHeader_GetMinExtent(curHeader) - curNodeOffset;
  maxCurCommonDelta = curNodeWidth-minCurCommonDelta;
  if ( minCurNodeNoDelta<maxCurNodeNoDelta )
  {
    voxelSafety = minCurNodeNoDelta*curNodeWidth;
    voxelSafety += minCurCommonDelta;
  }
  else if (maxCurNodeNoDelta < minCurNodeNoDelta)
       {
         voxelSafety = maxCurNodeNoDelta*curNodeWidth;
         voxelSafety += maxCurCommonDelta;
        }
        else
        {
          voxelSafety = minCurNodeNoDelta*curNodeWidth;
          voxelSafety += (((minCurCommonDelta)<(maxCurCommonDelta))?(minCurCommonDelta):(maxCurCommonDelta));
        }
  while ( (localVoxelDepth>0) && (voxelSafety>0) )
  {
    localVoxelDepth--;
    curHeader = This->fVoxelHeaderStack[localVoxelDepth];
    curHeaderAxis = This->fVoxelAxisStack[localVoxelDepth];
    curNodeNo = This->fVoxelNodeNoStack[localVoxelDepth];
    curNodeWidth = This->fVoxelSliceWidthStack[localVoxelDepth];
    curNodeOffset = curNodeNo*curNodeWidth;
    minCurCommonDelta = G4ThreeVector_coord(localPoint,curHeaderAxis)
                        - G4VoxelHeader_GetMinExtent(curHeader) - curNodeOffset;
    maxCurCommonDelta = curNodeWidth-minCurCommonDelta;
    if ( minCurCommonDelta<voxelSafety )
    {
      voxelSafety = minCurCommonDelta;
    }
    if ( maxCurCommonDelta<voxelSafety )
    {
      voxelSafety = maxCurCommonDelta;
    }
  }
  if ( voxelSafety<0 )
  {
    voxelSafety = 0;
  }
  return voxelSafety;
}
__device__
void G4VoxelNavigation_ctor( G4VoxelNavigation *This )
{
 This->fVoxelDepth = -1;
 This->fVoxelNode = 0;
}
__device__
G4bool
G4VoxelNavigation_LevelLocate(
   G4VoxelNavigation *This,
   G4NavigationHistory* history,
   const G4VPhysicalVolume* blockedVol,
   G4ThreeVector globalPoint,
   const G4ThreeVector* globalDirection,
   const G4bool pLocatedOnEdge,
   G4ThreeVector *localPoint )
{
  G4SmartVoxelHeader *targetVoxelHeader;
  G4SmartVoxelNode *targetVoxelNode;
  G4VPhysicalVolume *targetPhysical, *samplePhysical;
  G4LogicalVolume *targetLogical;
  G4VSolid *sampleSolid;
  G4ThreeVector samplePoint;
  G4int targetNoDaughters;
  targetPhysical = G4NavigationHistory_GetTopVolume(history);
  targetLogical = G4VPhysicalVolume_GetLogicalVolume(targetPhysical);
  targetVoxelHeader = G4LogicalVolume_GetVoxelHeader(targetLogical);
  targetVoxelNode =
 G4VoxelNavigation_VoxelLocate(This,targetVoxelHeader,*localPoint);
  targetNoDaughters=G4VoxelNode_GetNoContained(targetVoxelNode);
  if ( targetNoDaughters==0 ) return false;
  for ( int sampleNo=targetNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    samplePhysical =
  G4LogicalVolume_GetDaughter( targetLogical,
   G4VoxelNode_GetVolume(targetVoxelNode,sampleNo));
    if ( samplePhysical!=blockedVol )
    {
      G4NavigationHistory_NewLevel(history, samplePhysical, kNormal);
      sampleSolid =
  G4LogicalVolume_GetSolid(
   G4VPhysicalVolume_GetLogicalVolume( samplePhysical ));
   G4AffineTransform tf = G4NavigationHistory_GetTopTransform( history );
      samplePoint =
  G4AffineTransform_TransformPoint( &tf, globalPoint );
      if( G4AuxiliaryNavServices_CheckPointOnSurface(
   sampleSolid, samplePoint, globalDirection,
   tf, pLocatedOnEdge) )
      {
        *localPoint = samplePoint;
        return true;
      }
      else
      {
    G4NavigationHistory_BackLevel( history );
      }
    }
  }
  return false;
}
__device__
G4double
G4VoxelNavigation_ComputeStep(
   G4VoxelNavigation *This,
   G4ThreeVector localPoint,
   G4ThreeVector localDirection,
   const G4double currentProposedStepLength,
   G4double *newSafety,
   G4NavigationHistory *history,
   G4bool *validExitNormal,
   G4ThreeVector *exitNormal,
   G4bool *exiting,
   G4bool *entering,
   G4VPhysicalVolume *(*pBlockedPhysical)
    , G4double * Result
   )
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical,
 *blockedExitedVol = 0;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int sampleNo;
  G4bool initialNode, noStep;
  const G4SmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;
  motherPhysical = G4NavigationHistory_GetTopVolume( history );
  motherLogical = G4VPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = G4LogicalVolume_GetSolid(motherLogical);
  motherSafety = G4VSolid_DistanceToOut(motherSolid, localPoint);
  ourSafety = motherSafety;
  if ( *exiting && *validExitNormal )
  {
    if ( G4ThreeVector_dot(localDirection,*exitNormal)>=kMinExitingNormalCosine )
    {
      blockedExitedVol = *pBlockedPhysical;
      ourSafety = 0;
    }
  }
  *exiting = false;
  *entering = false;
  initialNode = true;
  noStep = true;
  while ( noStep )
  {
    curVoxelNode = This->fVoxelNode;
    curNoVolumes = G4VoxelNode_GetNoContained(curVoxelNode);
    for (contentNo=curNoVolumes-1; contentNo>=0; contentNo--)
    {
      sampleNo = G4VoxelNode_GetVolume( curVoxelNode, contentNo);
        samplePhysical = G4LogicalVolume_GetDaughter(motherLogical,sampleNo);
        if ( samplePhysical!=blockedExitedVol )
        {
    G4AffineTransform sampleTf =
   G4AffineTransform_create_full(
    G4VPhysicalVolume_GetObjectRotationValue(samplePhysical),
    G4VPhysicalVolume_GetTranslation(samplePhysical));
          G4AffineTransform_Invert(&sampleTf);
          const G4ThreeVector samplePoint =
    G4AffineTransform_TransformPoint(&sampleTf,localPoint);
          const G4VSolid *sampleSolid =
    G4LogicalVolume_GetSolid(
     G4VPhysicalVolume_GetLogicalVolume(
      samplePhysical ));
          const G4double sampleSafety =
   G4VSolid_DistanceToIn(sampleSolid,samplePoint);
          if ( sampleSafety<ourSafety )
          {
            ourSafety = sampleSafety;
          }
          if ( sampleSafety<=ourStep )
          {
            sampleDirection =
    G4AffineTransform_TransformAxis( &sampleTf, localDirection );
            G4double sampleStep =
    G4VSolid_DistanceToIn_full(sampleSolid, samplePoint, sampleDirection);
            if ( sampleStep<=ourStep )
            {
              ourStep = sampleStep;
              *entering = true;
              *exiting = false;
              *pBlockedPhysical = samplePhysical;
            }
        }
      }
    }
    if (initialNode)
    {
      initialNode = false;
      voxelSafety = G4VoxelNavigation_ComputeVoxelSafety(This,localPoint);
      if ( voxelSafety<ourSafety )
      {
        ourSafety = voxelSafety;
      }
      if ( currentProposedStepLength<ourSafety )
      {
        noStep = false;
        *entering = false;
        *exiting = false;
        *pBlockedPhysical = 0;
        ourStep = kInfinity;
      }
      else
      {
        if ( motherSafety<=ourStep )
        {
          G4double motherStep =
   G4VSolid_DistanceToOut_full( motherSolid, localPoint, localDirection,
                                         true, validExitNormal, exitNormal);
          if ( motherStep<=ourStep )
          {
            ourStep = motherStep;
            *exiting = true;
            *entering = false;
            if ( *validExitNormal )
            {
    G4RotationMatrix rot = G4VPhysicalVolume_GetObjectRotationValue(motherPhysical);
    G4RotationMatrix inv = G4RotationMatrix_inverse(&rot);
    *exitNormal = G4RotationMatrix_apply( &inv, *exitNormal );
            }
          }
          else
          {
            *validExitNormal = false;
          }
        }
      }
      *newSafety = ourSafety;
    }
    if (noStep)
    {
      noStep = G4VoxelNavigation_LocateNextVoxel(This, localPoint, localDirection, ourStep);
    }
  }
  int locationId = (blockIdx.x * blockDim.x + threadIdx.x);
  return ourStep;
}
__global__ void trace(
 Particle *input,
 G4double *output,
 G4VPhysicalVolume *worldVolumeAndGeomBuffer,
 G4double phys_step,
 int totalSize
 , G4double * Result
 , SolidInfo * Solids
 , ResultInfo * Result_For_Current_Solid
 , FinalResult * Compacter_Result,
 G4SmartVoxelNode * nullVNode
 )
{
 const unsigned globalIdx = (blockIdx.x * blockDim.x + threadIdx.x);
 const unsigned localIdx = threadIdx.x;
 const unsigned locationId = globalIdx;
 if (globalIdx >= totalSize ) return;
 __shared__ int Numbers_Of_Solid[ BlockSize ];
 __shared__ int Sum_Of_Solids[ BlockSize ];
 __shared__ bool noStepArray [ BlockSize ];
 __shared__ PointInformation LocationArray[ BlockSize ];
 __shared__ G4VPhysicalVolume * info[ BlockSize ];
 G4VoxelNode_ctor( nullVNode ,1 );
 G4bool Cur_Vol_Store [ BlockSize ];
 G4Navigator navi;
 G4Navigator *nav = &navi;
 G4Navigator_ctor(nav);
 G4Navigator_SetWorldVolume( nav, worldVolumeAndGeomBuffer );
 Particle p = input[globalIdx];
 if( globalIdx == 0)
 {
 }
    const G4VPhysicalVolume * cur_vol =
  G4Navigator_LocateGlobalPointAndSetup(
   nav, p.pos, 0, false, true, Result );
 G4bool cur_vol_local = true, cur_vol_all = true;
 G4double step, safety = 0.1;
 G4double integratedDensity = 0;
 int temp = 0;
 while ( cur_vol_all )
 {
  {
  const G4double curDensity =
    G4LogicalVolume_GetMaterial( G4VPhysicalVolume_GetLogicalVolume( cur_vol ))->property;
  PointInformation NewPoint = { p.pos, p.dir };
  LocationArray[ locationId ] = NewPoint;
  if( temp == 1)
  {
   Result[ locationId ] = step;
  }
  step = G4Navigator_ComputeStep( nav, p.pos, p.dir, phys_step, &safety
       , cur_vol_local
       , Result
       );
  if ( step == kInfinity ) step = phys_step;
  const G4double nextStepIntegratedD = curDensity * step;
  int locationId = (blockIdx.x * blockDim.x + threadIdx.x);
  integratedDensity += nextStepIntegratedD;
  G4ThreeVector_sum_assign( &(p.pos), G4ThreeVector_mult( p.dir, step ) );
  G4Navigator_SetGeometricallyLimitedStep( nav );
  if( globalIdx == 0 ){
  }
  cur_vol =
   G4Navigator_LocateGlobalPointAndSetup(
    nav, p.pos, &(p.dir), true, false, Result );
  if ( !cur_vol )
   cur_vol_local = false;
  }
  Cur_Vol_Store[ locationId ] = cur_vol_local;
  __syncthreads();
  cur_vol_all = NoStepReduction( Cur_Vol_Store, BlockSize );
  __syncthreads();
  temp++;
}
 output[globalIdx] = integratedDensity;
}
__global__ void relocate ( int * ptr, void * buf, int size )
{
  typedef unsigned char byte;
  const unsigned globalidx = (blockIdx.x * blockDim.x + threadIdx.x);
  if(globalidx>=size) return;
  int destoffs, targoffs;
  destoffs = *(ptr + 2*globalidx);
  targoffs = *(ptr + 2*globalidx+ 1);
  *((byte*)buf+destoffs) = (byte)((byte*)buf + targoffs);
}
__global__ void check( G4VPhysicalVolume *worldVolumeAndGeomBuffer, unsigned long * result)
{
 unsigned int hope = ( unsigned int )worldVolumeAndGeomBuffer;
 *result = hope;
}
__global__ void test ( bool * output
 , bool * input
 )
{
 int tid = (blockIdx.x * blockDim.x + threadIdx.x);
 int offset = 1;
 G4bool result;
 if( tid == 0)
 {
 input[ 0] = true;
 input[ 1] = true;
 input[ 2] = true;
 input[ 3] = true;
 input[ 4] = true;
 input[ 5] = true;
    input[ 6] = false;
    input[ 7] = true;
 }
 __syncthreads();
}
__global__ void checkgeom( G4VPhysicalVolume *worldVolumeAndGeomBuffer, int * result, int number_of_increments)
{
 const unsigned globalid = (blockIdx.x * blockDim.x + threadIdx.x);
 if(globalid>=1) return;
 int i=0;
 G4Navigator navi;
 G4Navigator *nav = &navi;
 G4Navigator_ctor(nav);
 G4Navigator_SetWorldVolume( nav, worldVolumeAndGeomBuffer );
 G4ThreeVector pos = G4ThreeVector_create( 0.0, 0.0, 0.0);
 const G4VPhysicalVolume * cur_vol;
 unsigned int geom_start = ( unsigned int )worldVolumeAndGeomBuffer;
 pos = G4ThreeVector_create( 0.7, 1.0, 0.7);
 float x_increment = 0.2, y_increment = 0.2, z_increment = 0.2;
 for( i=0; i < number_of_increments*3 ; i+=3)
 {
  result[i] = ( int ) cur_vol->count;
  result[i + 1] = (( unsigned int )(cur_vol->flogical) - geom_start);
  result[i + 2] = ( int ) G4LogicalVolume_GetMaterial( G4VPhysicalVolume_GetLogicalVolume( cur_vol ))->property;
  pos.x+=x_increment;
  pos.y+=y_increment;
  pos.z+=z_increment;
 }
}
struct CameraParameters
{
 double
  heading,
  pitch,
  roll,
  dist,
  yfov,
  target_x,
  target_y,
  target_z;
 CameraParameters()
 :
  heading(0), pitch(0), roll(0), dist(1),
  yfov(90), target_x(0), target_y(0), target_z(0)
 {}
};
struct EventOrigin
{
 double x,y,z;
};
class Geometry
{
public:
 typedef unsigned char byte;
 virtual ~Geometry() {}
 virtual void create() = 0;
 virtual void relocate( void *newbegin ) = 0;
 virtual int size() const = 0;
 virtual int ptrs_size() const=0;
 virtual void *getBuffer() = 0;
 virtual double getScale() const = 0;
 virtual CameraParameters getCamera() const = 0;
 virtual EventOrigin getEvent() const
 {
  EventOrigin e = { 0,0,0 };
  return e;
 }
 virtual int getNumVoxelNodes() const { return 0; }
};
typedef struct { const char *err, *fn; int line, errcode; } my_cuda_err;
typedef struct { int secs; int usecs; } mytimet;
extern "C"
{
 void myprint( const char *chr );
 void myprint1( const char *chr, int n );
 mytimet mytimer();
 void myprinttdiff(mytimet a, mytimet b);
 void mysleep(int n);
}
static inline int ceilDiv( int a, int d )
{
 return a/d + ((a%d)?1:0);
}
Particle *gpuInput;
G4double *gpuOutput;
Geometry::byte *gpuGeom;
int numInput, numOutput, numInputPerRound;
const int WARP_SIZE = 32;
void createGrid( int numInput, dim3* grid, dim3* block )
{
 const int MAXSIZE = 10000000;
 const int NUMCORES = 448;
 const int NUMMULTIPROC = 14;
 const int BLOCKS_PER_MULTIPROC = 8;
 const int MAX_WARPS_PER_MULTIPROC = 48;
 const int MAX_DATA_PER_MULTIPROC = MAX_WARPS_PER_MULTIPROC*WARP_SIZE;
 int size = numInput;
 if (size > MAXSIZE) size = MAXSIZE;
 int dataPerMultiproc = ceilDiv(size,NUMMULTIPROC);
 if ( dataPerMultiproc > MAX_DATA_PER_MULTIPROC )
  dataPerMultiproc = MAX_DATA_PER_MULTIPROC;
 int blockSize = ceilDiv(dataPerMultiproc,BLOCKS_PER_MULTIPROC);
 const int MAX_BLOCK_SIZE = 1024;
 if (blockSize > MAX_BLOCK_SIZE) blockSize = MAX_BLOCK_SIZE;
 int numBlocks = ceilDiv(size,blockSize);
 int numWarps = ceilDiv(blockSize,WARP_SIZE) * numBlocks;
 if (numWarps > NUMCORES)
 {
  blockSize = ceilDiv(blockSize,WARP_SIZE)*WARP_SIZE;
  dataPerMultiproc = blockSize * BLOCKS_PER_MULTIPROC;
  if ( dataPerMultiproc > MAX_DATA_PER_MULTIPROC )
   blockSize -= WARP_SIZE;
 }
 size = blockSize*ceilDiv(size,blockSize);
 if (size > MAXSIZE) size = MAXSIZE;
 block->x = blockSize;
 block->y = block->z = 1;
 grid->x = size/blockSize;
 grid->y = 1;
 grid->z = 1;
}
my_cuda_err cudainit( Geometry *geom, int N )
{
 const mytimet t0 = mytimer();
 numOutput = numInput = numInputPerRound = N;
 do { cudaError_t errc = cudaSetDeviceFlags(0); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 157, errc }; return r; } } while(0);
 do { cudaError_t errc = cudaMalloc( (void**)&gpuInput, sizeof(Particle)*numInput ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 159, errc }; return r; } } while(0);
 do { cudaError_t errc = cudaMalloc( (void**)&gpuOutput, sizeof(G4double)*numOutput ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 160, errc }; return r; } } while(0);
 do { cudaError_t errc = cudaMalloc( (void**)&gpuGeom, geom->size() ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 161, errc }; return r; } } while(0);
 geom->relocate( gpuGeom );
 cudaFuncSetCacheConfig(trace, cudaFuncCachePreferL1);
 do { cudaError_t errc = cudaMemcpy( gpuGeom, geom->getBuffer(), geom->size(), cudaMemcpyHostToDevice ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 197, errc }; return r; } } while(0);
 const mytimet t1 = mytimer();
 myprint("Initialization: ");
 myprinttdiff(t0, t1);
 my_cuda_err ok = { 0, 0, 0, cudaSuccess }; return ok;
}
my_cuda_err cudaexec( G4double phys_step, int totalInput, Particle *input, G4double *output )
{
   for ( int i = 0; i < totalInput; i += numInput )
   {
 if ( i + numInput > totalInput ) numInput = totalInput-i;
 do { cudaError_t errc = cudaMemcpy( gpuInput, input+i, sizeof(Particle)*numInput, cudaMemcpyHostToDevice ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 333, errc }; return r; } } while(0);
 dim3 grid, block;
 createGrid( numInput, &grid, &block );
 trace <<< grid, block >>>( gpuInput, gpuOutput, (G4VPhysicalVolume*)gpuGeom, phys_step, numInput, 0, 0, 0, 0, 0 );
 do { cudaError_t errc = cudaGetLastError(); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 340, errc }; return r; } } while(0);
 do { cudaError_t errc = cudaMemcpy( output+i, gpuOutput, sizeof(G4double)*numOutput, cudaMemcpyDeviceToHost ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 343, errc }; return r; } } while(0);
   }
   my_cuda_err ok = { 0, 0, 0, cudaSuccess }; return ok;
}
my_cuda_err cudafinish()
{
 const mytimet t0 = mytimer();
 do { cudaError_t errc = cudaFree( gpuInput ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 356, errc }; return r; } } while(0);
 do { cudaError_t errc = cudaFree( gpuOutput ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 357, errc }; return r; } } while(0);
 do { cudaError_t errc = cudaFree( gpuGeom ); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 358, errc }; return r; } } while(0);
 do { cudaError_t errc = cudaThreadExit(); if (errc != cudaSuccess) { my_cuda_err r = { cudaGetErrorString(errc), "cuda.cpp", 375, errc }; return r; } } while(0);
 const mytimet t1 = mytimer();
 myprint("Finalization: ");
 myprinttdiff(t0, t1);
 my_cuda_err ok = { 0, 0, 0, cudaSuccess }; return ok;
}
