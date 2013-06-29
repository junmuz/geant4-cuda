
/**
 * Implementation of common geometry functions such as pointer
 * relocation and memory allocation from a flat array.
 */

#ifndef GEOMETRYCOMMON_HPP
#define GEOMETRYCOMMON_HPP
#define CHECK 0
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "G4LogicalVolume.h"
#include "G4RotationMatrix.h"
#include "G4VPhysicalVolume.h"



#include "geometry.hpp"

#ifdef ENABLE_VOXEL_NAVIGATION
#include "G4Voxels.c"
#include "G4VoxelLimits.c"
#include "G4VoxelHeader.cpp"
#endif

	typedef struct
	{	
		G4VPhysicalVolume * Phys_vector_pointer;
		unsigned long offset;
		int property;
	}CheckPointer;
/**
 * Geometry implementation base class, contains the
 * common methods used by the non-abstract subclasses
 */
class BasicGeometry : public Geometry
{
private:
	/*typedef struct
	{
		int destoffs;
		int targoffs;
	}
	StrayPointer;
	*/
	// The original StrayPointer struct stored the offsets as integers. For 64-32 (CPU-GPU) and 64-64 bit compatibility, it makes more sense to store unsigned long instead
	typedef struct
	{
		unsigned long destoffs;
		unsigned long targoffs;
	}
	StrayPointer;
	


	const bool userbuf;
	bool create_called;
	
	void *buf;
	int bufsz, maxbufsz;
	int default_alignment;
	
	int nvoxelnodes;
	//EDIT: Making ptrs publically acccessible. Perhaps a better implementation would be to just create a method that returns it.
		//std::vector< StrayPointer > ptrs;
	
	template <class T> void addPointer( T *& ptr )
	{
		if ( ptr != NULL )
		{
			byte **p = (byte**)(&ptr);
			StrayPointer s = {
				(byte*)p - (byte*)buf,
				*p - (byte*)buf
			};
			//REMOVE
			//std::cout<<"Before hostptr destoffs = "<<s.destoffs<<" and targoffs= "<<s.targoffs<<std::endl;
			ptrs.push_back(s);
		}
	}
	
	void *reserveThingUnaligned( int thingsz )
	{
		assert(thingsz > 0);
		
		const int needed = bufsz + thingsz;
		if (needed > maxbufsz)
		{
			if (userbuf)
				throw std::runtime_error( "Buffer size exceeded" );
				
			maxbufsz = needed;
			void *newbuf = std::realloc( buf, maxbufsz );
			if ( buf != NULL && newbuf != buf )
			{
				buf = newbuf;
				throw std::runtime_error( "Buffer reallocated to different location" );
			}
			
			buf = newbuf;
		}
		void * const bufpos = (void*)((byte*)buf+bufsz);
		bufsz = needed;
		return bufpos;
	}
	
	void *reserveThingAligned( int thingsz, int alignment )
	{
		assert(thingsz>0 && alignment>0);
		
		int padding = 0;
		if ( bufsz % alignment > 0 )
		{
			padding = alignment - (bufsz%alignment);
		}
		
		byte *m = (byte*)reserveThingUnaligned( thingsz + padding );
		
		const byte PADDING_BYTE = 0x0;
		
		// padding contents
		while ( padding-- > 0 ) *m++ = PADDING_BYTE;
		return (void*)m;
	}
	
	void *reserveThing( int thingsz )
	{
		return reserveThingAligned( thingsz, default_alignment );
	}
	
	void *addThingAligned( const void *thing, int thingsz, int alignment )
	{
		void *ptr = reserveThingAligned(thingsz, alignment);
		std::memcpy( ptr, thing, thingsz );
		return ptr;
	}
	
	void *addThingUnaligned( const void *thing, int thingsz )
	{
		void *ptr = reserveThingUnaligned(thingsz);
		std::memcpy( ptr, thing, thingsz );
		return ptr;
	}
	
	void *addThing( const void *thing, int thingsz )
	{
		void *ptr = reserveThing(thingsz);
		std::memcpy( ptr, thing, thingsz );
		return ptr;
	}
	
protected:
	
	template <class T> T* reserveThing()
	{
		return static_cast<T*>(reserveThing(sizeof(T)));
	}
	
	template <class T> T* reserveThingAligned( int alignment = sizeof(T) )
	{
		return static_cast<T*>(reserveThingAligned(sizeof(T), alignment));
	}
	
	template <class T> T* reserveThingUnaligned()
	{
		return static_cast<T*>(reserveThingUnaligned(sizeof(T)));
	}
	
	template <class T> T* reserveNThingsUnaligned( int N )
	{
		return static_cast<T*>(reserveThingUnaligned(sizeof(T)*N));
	}
	
	template <class T> T* reserveNThingsAligned( int N, int alignment = sizeof(T) )
	{
		if ( sizeof(T) % alignment )
			throw std::runtime_error( "Incorrect alignment for array" );
			
		return static_cast<T*>(reserveThingAligned(sizeof(T)*N, alignment));
	}
	
	template <class T> T* reserveNThingsSelfAligned( int N )
	{
		return reserveNThingsAligned<T>( N );
	}
	
	template <class T> T* reserveNThings( int N )
	{
		return reserveNThingsAligned<T>( N, default_alignment );
	}
	
	template <class T> T* addThing( const T &thing )
	{
		return static_cast<T*>( addThing( &thing, sizeof(T) ) );
	}
	
	template <class T> T* addThingAligned( const T &thing, int alignment = sizeof(T) )
	{
		return static_cast<T*>( addThingAligned( &thing, sizeof(T), alignment ) );
	}
	
	template <class T> T* addThingUnaligned( const T &thing )
	{
		return static_cast<T*>( addThingUnaligned( &thing, sizeof(T) ) );
	}
	
	void addLogicalVolumePointers( G4LogicalVolume *vol )
	{   //MODIFY: May need if here.
		assert( vol->fMaterial != NULL );
		assert( vol->fSolid != NULL );
		
		addPointer( vol->fMaterial );
		addPointer( vol->fSolid );
	
		addPointer( vol->fDaughters );
		for (int i=0; i<vol->fNoDaughters; ++i )
			addPointer( vol->fDaughters[i] );
		
		//vol->check = 5454; Written ony during the testing phase.


		#ifdef ENABLE_VOXEL_NAVIGATION
			addPointer( vol->fVoxel );
		#endif
	}
	
	void addPhysicalVolumePointers( G4VPhysicalVolume *vol )
	{
		assert( vol->flogical != NULL );

	
		unsigned long offset1 = (byte *)(vol->flogical)-( byte *)buf;
		unsigned long offset2 = (byte *)(vol->flmother)-( byte *)buf;
		// These offsets are part of the Geometry check, they are ccompared to offsets calculated from the GPU.


		//EDIT
			if(CHECK == 2 || CHECK==4)
			{	
				CheckPointer Pointer = { vol, offset1, vol->flogical->fMaterial->property};
				VolumeStore.push_back(Pointer);
			}
		vol->counter++;
		vol->count = vol->counter;
		// Count is an integer that each PhysicalVolume stores which is a chronological count of when it was added.
		// The counter is a static int that increments whenever a new Physical Volume is added. Unfortunately in the current OpenCL C standard
		// a struct cannot contain a static member. Therefore the definition of counter in Physical VOlume may seem odd.

		
		/*REMOVE
		if(vol->counter<200)
		{
			std::cout<<"Physical Volume Count = "<<vol->count<<" and flogical Address = "<<vol->flogical<<" and Flmother address = "<<vol->flmother<<" and buf is at "<<buf<<std::endl;
			std::cout<<" Also, flogical offset = "<<offset1<<" and Flmother offset = "<<offset2<<std::endl;
		}
		*/
		//EDIT2
		//assert(offset1<2E9);
		//assert(offset2<2E9);
		
		// The assery statement above was initially to check if the offsets were within the permissible limit. But now it seems easier to just check that the maxbufsz is within
		// the limit.


		addPointer( vol->flogical );
		addPointer( vol->flmother );
	}
	
	#ifdef ENABLE_G4POLYCONE
	void addPolyConePointers( G4PolyCone *polycone )
	{
		assert( polycone->fNumPlanes >= 2 );
		assert( polycone->fz != NULL );
		assert( polycone->frmin != NULL );
		assert( polycone->frmax != NULL );
		
		addPointer( polycone->fz );
		addPointer( polycone->frmin );
		addPointer( polycone->frmax );
	}
	#endif
	
	/**
	 * Add daughter to logical volume data structure
	 * the daughters for one volume must be added consecutively,
	 * no other "add" or "reserve" functions may be called inbetween
	 */
	void addLogicalVolumeDaughter( G4LogicalVolume *vol, G4VPhysicalVolume *d )
	{
		// Constructing a C-array --> no padding / alignment!
		G4VPhysicalVolume **dptr = addThingUnaligned( d );
		
		if (vol->fDaughters == NULL)
		{
			vol->fDaughters = dptr;
		}
		vol->fNoDaughters++;
	}
	
	#ifdef ENABLE_VOXEL_NAVIGATION
	
	G4SmartVoxelNode *storeVoxelNode( const G4SmartVoxelNode &original )
	{
		nvoxelnodes++;
		//REMOVE
		//std::cout<< "From voxel node the count apparent is "<<original.SolidType[0]<<std::endl;
		G4SmartVoxelNode *node = addThing(original);
		if ( original.fcontents != NULL )
		{
			node->fcontents = static_cast<G4int*>(
				addThing( original.fcontents, sizeof(G4int)*original.fNumContents ));
			addPointer( node->fcontents );
		}
		return node;
	}
	
	G4SmartVoxelHeader *storeVoxelHeader( const G4SmartVoxelHeader &original )
	{
		G4SmartVoxelHeader *v = addThing(original);
		assert( v->fNumSlices > 0 && v->fslices != NULL );
		
		v->fslices = reserveNThingsSelfAligned<G4SmartVoxelProxy*>( v->fNumSlices );
		addPointer( v->fslices );
		
		for ( int i=0; i<v->fNumSlices; ++i )
		{
			assert( original.fslices[i] != NULL );
			
			v->fslices[i] = addThing( *original.fslices[i] );
			addPointer( v->fslices[i] );
			
			if ( G4VoxelProxy_IsNode( original.fslices[i] ) )
			{
				v->fslices[i]->fNode = storeVoxelNode( *original.fslices[i]->fNode );
				addPointer( v->fslices[i]->fNode );
			}
			else
			{
				assert( original.fslices[i]->fHeader != NULL );
				
				v->fslices[i]->fHeader =
					storeVoxelHeader( *original.fslices[i]->fHeader );
					
				addPointer( v->fslices[i]->fHeader );
			}
		}
		return v;
	}
	
	void createVoxels( G4LogicalVolume *mother, G4double smartless = kInfinity )
	{
		G4SmartVoxelHeader head;
		G4VoxelHeader_ctor( &head, mother, 0, smartless );
		
		if ( head.fNumSlices > 0 )
			G4LogicalVolume_SetVoxelHeader(
				mother, storeVoxelHeader( head ) );
		
		G4VoxelHeader_dtor( &head );
	}
	#endif
	
	G4double uniformRand()
	{
		// TODO: better rand?
		//MODIFY: Better rand
		const G4int MAGN = 10000;
		return (std::rand() % MAGN) / (G4double)MAGN;
	}
	
	G4RotationMatrix randomRot()
	{
		G4RotationMatrix rot = G4RotationMatrix_create_id();
		
		// I know, not uniform in SO(3)
		// does not matter...
		
		const G4double angX = uniformRand() * 2.0 * M_PI;
		const G4double angY = uniformRand() * 2.0 * M_PI;
		const G4double angZ = uniformRand() * 2.0 * M_PI;
		
		G4RotationMatrix_rotateX( &rot, angX );
		G4RotationMatrix_rotateY( &rot, angY );
		G4RotationMatrix_rotateZ( &rot, angZ );
		
		return rot;
	}
	
	/**
	 * Actual "create" implementation. No guards against
	 * this begin called multiple times.
	 */
	virtual void create_impl() = 0;
	
public:

	/**
	 * Initialize geometry object, use the given buffer for the
	 * geometry. Ownership of the buffer is not transferred.
	 * 
	 * @param buffer buffer to use for the geometry
	 * @param maxsz maximum size of the buffer
	 * @param def_align default memory alignment for geometry objects
	 */
	std::vector< StrayPointer > ptrs;

	std::vector <CheckPointer>  VolumeStore;
	
			// Vector to store all Physical Volume Pointers to confirm the geometry transfer.
	
	BasicGeometry( void *buffer, int maxsz, int def_align = sizeof(G4double) )
	:
		userbuf(true),
		create_called(false),
		buf(buffer),
		bufsz(0),
		maxbufsz(maxsz),
		default_alignment(def_align),
		nvoxelnodes(0)
	{
		assert(default_alignment > 0);

		// maxsz that can be allocated should be less than the size of the gpu that we are using.
		// The size has actually to be limited by the size of GPU global memory and size that unsigned long can hold.
	}
	
	/**
	 * Initialize geometry object and allocate a geometry buffer of
	 * given size.
	 * 
	 * @param maxsz maximum size of the buffer
	 * @param def_align default memory alignment for geometry objects
	 */
	BasicGeometry( int maxsz = 0, int def_align = sizeof(G4double) )
	:
		userbuf(false),
		create_called(false),
		bufsz(0),
		maxbufsz(maxsz),
		default_alignment(def_align),
		nvoxelnodes(0)
	{
		if ( maxsz == 0 )
			buf = NULL;
		else
			buf = std::malloc( maxsz );
			
		assert(default_alignment > 0);
	}
	
	~BasicGeometry()
	{
		if (!userbuf) std::free( buf );
	}
	
	
	
	int size() const { return bufsz; }
	int ptrs_size() const {return ptrs.size();}

	std::vector<StrayPointer> pointer() { return ptrs;}
	void *getBuffer() { return buf; }
	
	/** Get default memory alignment for the geometry objects */
	int getAlignment() const { return default_alignment; }
	

	/** @return the total number of voxel nodes in the buffer */
	int getNumVoxelNodes() const { return nvoxelnodes; }
	
	void relocate( void *newbegin )
	{
		for ( unsigned i=0; i<ptrs.size(); ++i )
		{
			*(byte**)((byte*)buf+ptrs[i].destoffs) =
				(byte*)newbegin + ptrs[i].targoffs;
		}
	}
	

	void create()
	{
		printf ("Entered Create Method\n");
		if ( create_called )
			throw std::runtime_error("Geometry::create called twice");
		else create_impl();
		printf ("Leaving Create Method\n");
	}
};

#endif
