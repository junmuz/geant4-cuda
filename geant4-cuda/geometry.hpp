
/** Basic Geometry interfaces */

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

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

/**
 * Geometry interface class
 */
class Geometry
{
public:
	
	typedef unsigned char byte;
	
	virtual ~Geometry() {}

	/**
	 * Create the geometry. Must be called exactly once.
	 */
	virtual void create() = 0;
	
	/**
	 * Change the geometry buffer base address.
	 * 
	 * Does not move the array but changes the internal pointers
	 * in the geometry data structure so that the buffer may
	 * be copied to a memory segment beginning from the given address.
	 */
	virtual void relocate( void *newbegin ) = 0;
	
	
	
	/** Get geometry buffer size in bytes */
	virtual int size() const = 0;
	
	/** Get size of ptrs in bytes**/
	virtual int ptrs_size() const=0;


	
	/** Get a pointer to the geometry buffer */
	virtual void *getBuffer() = 0;
	
	/** Get the scale (some dimension) of the geometry */
	virtual double getScale() const = 0;
	
	/** Get camera parameters */
	virtual CameraParameters getCamera() const = 0;
	
	
	virtual EventOrigin getEvent() const
	{
		EventOrigin e = { 0,0,0 };
		return e;
	}

	/** Get the total number of voxel nodes in the buffer */
	virtual int getNumVoxelNodes() const { return 0; }
};

#endif
