
#ifndef __CL_HPP__
#define __CL_HPP__

#include <CL/opencl.h>

#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <list>
#include <cassert>

#include "utils.h"

namespace CL
{

// Helper function to get error string, (c) NVIDIA
// TODO: rewrite (?)
static const char* oclErrorString(cl_int error)
{
    static const char* errorString[] = {
        "CL_SUCCESS",
        "CL_DEVICE_NOT_FOUND",
        "CL_DEVICE_NOT_AVAILABLE",
        "CL_COMPILER_NOT_AVAILABLE",
        "CL_MEM_OBJECT_ALLOCATION_FAILURE",
        "CL_OUT_OF_RESOURCES",
        "CL_OUT_OF_HOST_MEMORY",
        "CL_PROFILING_INFO_NOT_AVAILABLE",
        "CL_MEM_COPY_OVERLAP",
        "CL_IMAGE_FORMAT_MISMATCH",
        "CL_IMAGE_FORMAT_NOT_SUPPORTED",
        "CL_BUILD_PROGRAM_FAILURE",
        "CL_MAP_FAILURE",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "CL_INVALID_VALUE",
        "CL_INVALID_DEVICE_TYPE",
        "CL_INVALID_PLATFORM",
        "CL_INVALID_DEVICE",
        "CL_INVALID_CONTEXT",
        "CL_INVALID_QUEUE_PROPERTIES",
        "CL_INVALID_COMMAND_QUEUE",
        "CL_INVALID_HOST_PTR",
        "CL_INVALID_MEM_OBJECT",
        "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
        "CL_INVALID_IMAGE_SIZE",
        "CL_INVALID_SAMPLER",
        "CL_INVALID_BINARY",
        "CL_INVALID_BUILD_OPTIONS",
        "CL_INVALID_PROGRAM",
        "CL_INVALID_PROGRAM_EXECUTABLE",
        "CL_INVALID_KERNEL_NAME",
        "CL_INVALID_KERNEL_DEFINITION",
        "CL_INVALID_KERNEL",
        "CL_INVALID_ARG_INDEX",
        "CL_INVALID_ARG_VALUE",
        "CL_INVALID_ARG_SIZE",
        "CL_INVALID_KERNEL_ARGS",
        "CL_INVALID_WORK_DIMENSION",
        "CL_INVALID_WORK_GROUP_SIZE",
        "CL_INVALID_WORK_ITEM_SIZE",
        "CL_INVALID_GLOBAL_OFFSET",
        "CL_INVALID_EVENT_WAIT_LIST",
        "CL_INVALID_EVENT",
        "CL_INVALID_OPERATION",
        "CL_INVALID_GL_OBJECT",
        "CL_INVALID_BUFFER_SIZE",
        "CL_INVALID_MIP_LEVEL",
        "CL_INVALID_GLOBAL_WORK_SIZE",
    };

    const int errorCount = sizeof(errorString) / sizeof(errorString[0]);

    const int index = -error;

    return (index >= 0 && index < errorCount) ? errorString[index] : "";
}

#define CHECKERR(e) Exception::checkerr(e, __FILE__, __LINE__)
#define CHECKERRD(e, d) Exception::checkerr(e, __FILE__, __LINE__, d)

class Exception : public std::runtime_error
{
private:
	cl_int code;
	const char *desc;
	const char *file;
	int line;

public:
	Exception( cl_int code_, const char *file_, int line_, const char *desc_ )
	: std::runtime_error("CL::Exception"), code(code_), desc(desc_), file(file_), line(line_)
	{}

	~Exception() throw() {}
	
	virtual const char *what() const throw()
	{
		static std::string msg;
		std::ostringstream oss;
		oss << file << ":" << line << ": " << oclErrorString(code);
		if (desc != NULL && desc[0]!='\0') oss << " (" << desc << ")";
		msg = oss.str();
		return msg.c_str();
	}

	static void checkerr( cl_int errcode, const char *file, int line, const char *desc = NULL )
	{
		if (errcode != CL_SUCCESS)
			throw Exception( errcode, file, line, desc );
	}
};

std::string getPlatformName( cl_platform_id platform )
{
	const unsigned MAXLEN = 1024;
	char buf[MAXLEN+1]; buf[MAXLEN] = '\0';
	cl_int errcode = clGetPlatformInfo(platform, CL_PLATFORM_NAME, MAXLEN, &buf, NULL);
	CHECKERR(errcode);
	return std::string(buf);
}

std::string getDeviceName( cl_device_id device )
{
	const unsigned MAXLEN = 1024;
	char buf[MAXLEN+1]; buf[MAXLEN] = '\0';
	cl_int errcode = clGetDeviceInfo(device, CL_DEVICE_NAME, MAXLEN, &buf, NULL);
	CHECKERR(errcode);
	return std::string(buf);
}

std::list<cl_platform_id> getPlatformIDs()
{
	std::list<cl_platform_id> l;

	// get number of platforms
	cl_uint numPlatforms;
	cl_int errcode = clGetPlatformIDs(0,NULL,&numPlatforms);
	CHECKERR(errcode);

	if ( numPlatforms > 0 )
	{
		// get platforms
		cl_platform_id * platforms = new cl_platform_id[numPlatforms];
		try
		{
			errcode = clGetPlatformIDs(numPlatforms, platforms, NULL);
			CHECKERR(errcode);

			for ( cl_uint i=0; i<numPlatforms; ++i ) l.push_back(platforms[i]);
		}
		catch(...)
		{
			delete [] platforms;
			throw;
		}
		delete [] platforms;
	}

	return l;
}

void listAndPrintPlatforms()
{
	std::list<cl_platform_id> platforms = getPlatformIDs();
	if ( platforms.size() == 0 ) std::cout << "No OpenCL platforms found" << std::endl;
	else
	{
		std::cout << "Available OpenCL platforms (selecting the first one)" << std::endl;
		for ( std::list<cl_platform_id>::const_iterator itr = platforms.begin();
			itr != platforms.end(); ++itr )
		{
			std::cout << "\t" << getPlatformName( *itr ) << std::endl;
		}
	}
}

cl_platform_id selectPlatform()
{
	std::list<cl_platform_id> platforms = getPlatformIDs();
	if ( platforms.size() == 0 ) throw std::runtime_error("No OpenCL platforms found");
	return *platforms.begin();
}


class DeviceList
{
private:
	cl_uint numDevices;
	cl_device_id *devices;
	
	DeviceList(const DeviceList&);
	DeviceList& operator=(const DeviceList&);
	
public:
	DeviceList( cl_platform_id platform_id, cl_device_type deviceTypes )
	: numDevices(0), devices(NULL)
	{
		// get number of device IDs
		cl_int errcode = clGetDeviceIDs(platform_id, deviceTypes, 0, NULL, &numDevices);
		CHECKERRD(errcode, "failed to get device list size");

		// get device ID list		
		devices = new cl_device_id[ numDevices ];
    	errcode = clGetDeviceIDs(platform_id, deviceTypes, numDevices, devices, NULL);
    	try { CHECKERRD(errcode, "failed to get device list"); }
    	catch( ... ) { delete [] devices; throw; }
	}
	
	// Device list associated with a program
	DeviceList( cl_program program )
	: numDevices(0), devices(NULL)
	{
		// get number of device IDs
		cl_int errcode = clGetProgramInfo( program, CL_PROGRAM_NUM_DEVICES, sizeof(numDevices), &numDevices, NULL );
		CHECKERRD(errcode, "failed to get device list size");
		
		// get device ID list
		devices = new cl_device_id[ numDevices ];
    	errcode = clGetProgramInfo(program, CL_PROGRAM_DEVICES, numDevices*sizeof(cl_device_id), devices, NULL);
		try { CHECKERRD(errcode, "failed to get device list"); }
    	catch( ... ) { delete [] devices; throw; }
	}
	
	virtual ~DeviceList()
	{
		delete [] devices;
	}

	void printInfo()
	{
		std::cout << "Listed devices" << ":" << std::endl;
		for ( cl_uint i=0; i<numDevices; ++i )
		{
			std::cout << "\t" << getDeviceName( devices[i] ) << std::endl;
		}
	}

	cl_device_id getFirstDevice() const
	{
		if ( numDevices == 0 )
			throw std::out_of_range("No devices found");

		return devices[0];
	}

	cl_uint getNumDevices() const { return numDevices; }
	cl_device_id *getDevices() { return devices; }
};

class Context
{
private:
	cl_context handle;

	Context(const Context&);
	Context& operator=(const Context&);

public:
	cl_context getHandle() { return handle; }

	Context( DeviceList &devices )
	{
		cl_int errcode;
		// No properties or error callbacks
		handle = clCreateContext(0, devices.getNumDevices(), devices.getDevices(), NULL, NULL, &errcode);
		CHECKERR(errcode);
	}

	virtual ~Context()
	{
    		clReleaseContext(handle);
	}
};

class Buffer
{
private:
	
	const size_t sz;
	
	cl_mem handle;

public:
	// EDIT: Perhaps the error is caused by this. If it is changed it should go away.
	Buffer( Context &context, cl_mem_flags flags, size_t size )
	: sz(size)
	{
		cl_int errcode;
		// do not realloc (?)
		handle = clCreateBuffer( context.getHandle(), flags, sz, NULL, &errcode );
		CHECKERRD(errcode, "failed to allocate buffer");
	}

	virtual ~Buffer()
	{
		clReleaseMemObject(handle);
	}

	inline size_t size() const { return sz; }
	inline cl_mem getHandle() { return handle; }
	inline const cl_mem* getHandlePtr() const { return &handle; }
};

class Program
{
private:
	cl_program handle;
	bool built;

	Program(const Program&);
	Program& operator=(const Program&);

protected:

	void buildProgram( const char *compilerOptions )
	{
		// Build for all devices in the context synchronously
		std::cout<<"Works till here\n";
		//EDIT: Adding -g option
		cl_int errcode = clBuildProgram(handle, 0, NULL, "-cl-single-precision-constant", NULL, NULL);
		
		if ( errcode == CL_BUILD_PROGRAM_FAILURE )
		{
			// Compiler error, query build log
			std::string logs = getBuildLog();
			
			throw std::runtime_error("Failed to compile CL source.\n"+logs);
		}
		else
			CHECKERRD(errcode, "program compilation failed");
			
		built = true;
	}
	
	void buildFromSource(
		Context &context,
		const std::string &source,
		const char *compilerOptions )
	{
		cl_int errcode;
		const size_t len = source.length();
		const char *s = source.c_str();

		// Single source string
		handle = clCreateProgramWithSource(context.getHandle(), 1, &s, &len, &errcode);
		CHECKERRD(errcode, "program creation failed");
		
		buildProgram( compilerOptions );
	}
	
	void buildFromBinaries( 
		Context &context,
		DeviceList& devices,
		const std::vector<std::string> &binaries,
		const char *compilerOptions )
	{
		assert( devices.getNumDevices() == binaries.size() );
		
		std::vector<size_t> lens;
		std::vector<const unsigned char*> bins;
		std::vector<cl_int> stats( binaries.size() );
		
		for ( std::size_t i=0; i<binaries.size(); ++i )
		{
			lens.push_back( binaries[i].size() );
			bins.push_back( (const unsigned char*)(binaries[i].c_str()) );
		}
		
		cl_int errcode;
		handle = clCreateProgramWithBinary( context.getHandle(),
				devices.getNumDevices(), devices.getDevices(),
				&(lens[0]), &(bins[0]), &(stats[0]), &errcode );
				
		CHECKERRD(errcode, "program creation failed");
		
		buildProgram( compilerOptions );
	}
	
	Program() : built(false) {}

public:
	/** Build program from source */
	Program( Context &context, const std::string &source, const char *compilerOptions = NULL )
	:
		built(false)
	{
		buildFromSource( context, source, compilerOptions );
	}
	
	/** Create program from binaries or bytecode */
	Program(
		Context &context,
		DeviceList& devices,
		const std::vector<std::string> &binaries,
		const char *compilerOptions = NULL )
	:
		built(false)
	{
		buildFromBinaries( context, devices, binaries, compilerOptions );
	}

	virtual ~Program()
	{
		if(built) clReleaseProgram(handle);
	}
	
	std::string getBuildLog()
	{
		std::ostringstream oss;
		DeviceList devlist( handle );
		const size_t numDevs = devlist.getNumDevices();
		cl_device_id *devs = devlist.getDevices();
		
		for ( size_t i=0; i < numDevs; ++i )
		{
			const cl_device_id dev = devs[i];
			
			char *log;
			size_t logsz;
			
			// query size
			cl_int errcode = clGetProgramBuildInfo( handle, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &logsz );
			CHECKERRD( errcode, "failed to get program build info size" );
			
			// get log
			log = new char[ logsz+1 ];
			errcode = clGetProgramBuildInfo( handle, dev, CL_PROGRAM_BUILD_LOG, logsz, log, NULL );
			try
			{
				CHECKERRD( errcode, "failed to get program build info" );
				oss << "Build log for " << getDeviceName( dev ) << ":\n";
				log[logsz] = '\0';
				oss << log << "\n";
			}
			catch( ... ) { delete [] log; throw; }
			delete [] log;
		}
		
		return oss.str();
	}
	
	cl_program getHandle() { return handle; }
	
	void writeBinary( std::ostream &s, size_t devno = 0 )
	{
		DeviceList devlist( handle );
		const size_t numDevs = devlist.getNumDevices();
		
		if ( devno >= numDevs )
			throw std::runtime_error("Invalid device index");
			
		std::vector<size_t> binSizes( numDevs );
		
		cl_int errcode = clGetProgramInfo( handle,
			CL_PROGRAM_BINARY_SIZES, sizeof(size_t)*numDevs,
			&(binSizes[0]), NULL);
		CHECKERRD( errcode, "failed to get program info sizes" );
		
		unsigned totalsize = 0;
		for ( unsigned i=0; i<binSizes.size(); ++i )
			totalsize += ++binSizes[i];
			
		std::vector<unsigned char> buf( totalsize );
		std::vector<unsigned char*> binaries( binSizes.size() );
		
		unsigned memh = 0;
		for ( unsigned i=0; i<binSizes.size(); ++i )
		{
			binaries[i] = &(buf[memh]);
			memh += binSizes[i];
			buf[memh-1] = '\0';
		}
		
		clGetProgramInfo( handle, CL_PROGRAM_BINARIES,
			numDevs*sizeof(unsigned char*), &(binaries[0]), NULL );
		
		std::ostringstream oss;
		s << binaries[devno];
	}
	
	std::string getBinary( size_t devno = 0 )
	{
		std::ostringstream oss;
		writeBinary( oss, devno );
		return oss.str();
	}
};

class Kernel
{
private:
	cl_kernel handle;

	Kernel(const Kernel&);
	Kernel& operator=(const Kernel&);
	
public:
	Kernel( Program &p, std::string funcname )
	{
		cl_int errcode;
		handle = clCreateKernel(p.getHandle(), funcname.c_str(), &errcode);
		CHECKERRD( errcode, "kernel selection failed" );
	}

	virtual ~Kernel()
	{
		clReleaseKernel(handle);
	}

	cl_kernel getHandle() { return handle; }

	void setArg( cl_uint index, size_t arg_size, const void *val )
	{
		const cl_int errcode = clSetKernelArg( handle, index, arg_size, val );
		CHECKERRD(errcode, "failed to set kernel argument");
	}
	
	void setArg( cl_uint index, const Buffer& buf )
	{
		setArg( index, sizeof(cl_mem), buf.getHandlePtr() );
	}
};

class Event
{
private:
	cl_event handle;
	
public:
	Event( cl_event returned ) : handle(returned) {}
	
	Event( const Event& ev ) : handle( ev.handle )
	{
		clRetainEvent( handle );
	}
	
	~Event() { clReleaseEvent( handle ); }
	
	Event& operator=( const Event& ev )
	{
		if (&ev == this) return *this;
		
		clReleaseEvent( handle );
		handle = ev.handle;
		clRetainEvent( handle );
		
		return *this;
	}
	
	cl_event getHandle() { return handle; }
};

class CommandQueue
{
private:
	cl_command_queue handle;
	Context &context;

public:
	CommandQueue( Context &con, cl_device_id device, cl_command_queue_properties properties = 0 )
	:
		context(con)
	{
		cl_int errcode;
		handle = clCreateCommandQueue(context.getHandle(), device, properties, &errcode);
		CHECKERRD(errcode, "failed to create command queue");
	}

	virtual ~CommandQueue()
	{
		clReleaseCommandQueue(handle);
	}
	
	Context& getContext() { return context; }

	/** if bytes == 0, then transfer size is set to buffer maximum */
	Event enqueueWriteBuffer( Buffer& buffer, const void *host_buffer, cl_bool blocking = CL_TRUE,
		size_t offset = 0, size_t bytes = 0, cl_uint waitListSz = 0, const cl_event * waitList = NULL )
	{
		cl_event event;
		if ( bytes == 0 ) bytes = buffer.size();
		const cl_int errcode = clEnqueueWriteBuffer( handle, buffer.getHandle(), blocking, offset,
			bytes, host_buffer, waitListSz, waitList, &event );

		CHECKERRD(errcode, "failed to enqueue buffer write");
		return Event(event);
	}

	/** if bytes == 0, then transfer size is set to buffer maximum */
	Event enqueueReadBuffer( Buffer& buffer, void *host_buffer, cl_bool blocking = CL_TRUE,
		size_t offset = 0, size_t bytes = 0, cl_uint waitListSz = 0, const cl_event * waitList = NULL )
	{
		cl_event event;
		if ( bytes == 0 ) bytes = buffer.size();
		const cl_int errcode = clEnqueueReadBuffer( handle, buffer.getHandle(), blocking, offset,
			bytes, host_buffer, waitListSz, waitList, &event );

		CHECKERRD(errcode, "failed to enqueue buffer read");
		return Event(event);
	}

	void *enqueueMapBuffer( Buffer& buffer, cl_map_flags flags, cl_bool blocking = CL_TRUE, cl_event *event = NULL,
		size_t offset = 0, size_t bytes = 0, cl_uint waitListSz = 0, const cl_event * waitList = NULL )
	{
		if ( bytes == 0 ) bytes = buffer.size();
		cl_int errcode;
		void *ptr = clEnqueueMapBuffer( handle, buffer.getHandle(), blocking, flags,
			offset, bytes, waitListSz, waitList, event, &errcode );
		CHECKERRD(errcode, "failed to enqueue buffer map");
		return ptr;
	}
	
	Event enqueueUnmapBuffer( Buffer& buffer, void* hostptr,
		cl_uint waitListSz = 0, const cl_event * waitList = NULL )
	{
		cl_event event;
		const cl_int errcode = clEnqueueUnmapMemObject( handle, buffer.getHandle(), hostptr,
			waitListSz, waitList, &event );
		CHECKERRD(errcode, "failed to enqueue buffer unmap");
		return Event(event);
	}

	void enqueueBarrier()
	{
		const cl_int errcode = clEnqueueBarrier( handle );
		CHECKERRD(errcode, "failed to enqueue barrier");
	}

	Event enqueueMarker()
	{
		cl_event event;
		const cl_int errcode = clEnqueueMarker( handle, &event );
		CHECKERRD(errcode, "failed to enqueue barrier");
		return Event(event);
	}

	void enqueueWaitForEvents( cl_uint numEvents, const cl_event * eventList )
	{
		const cl_int errcode = clEnqueueWaitForEvents( handle, numEvents, eventList );
		CHECKERRD(errcode, "failed to enqueue barrier");
	}

	void flush()
	{
		const cl_int errcode = clFlush( handle );
		CHECKERRD(errcode, "failed to flush queue");
	}

	void finish()
	{
		const cl_int errcode = clFinish( handle );
		CHECKERRD(errcode, "failed to finish queue");
	}

	Event enqueueNDRangeKernel( Kernel& kernel, cl_uint work_dim,
		const size_t *global_work_size, const size_t *local_work_size = NULL,
		cl_uint waitListSz = 0, const cl_event *waitList = NULL )
	{
		cl_event event;
		// According to OpenCL 1.0 spec, global_work_offset must currently be NULL
		
		// EDIT: REMOVE
		unsigned int * retur;
		retur = (unsigned int *)malloc(sizeof(int));
		/*clGetKernelWorkGroupInfo (	kernel.getHandle(),
			NULL ,
 			CL_KERNEL_WORK_GROUP_SIZE,
 			sizeof(size_t),
 			retur,
 			NULL);
		std::cout<< "The max work group size for kernel -> " << *retur<< "\n";
		clGetKernelWorkGroupInfo (	kernel.getHandle(),
			NULL ,
 			CL_KERNEL_LOCAL_MEM_SIZE,
 			sizeof(size_t),
 			retur,
 			NULL);
		std::cout<< "The max local mem size is -> " << *retur<< "\n";*/
		const cl_int errcode = clEnqueueNDRangeKernel( handle, kernel.getHandle(),
			work_dim, NULL, global_work_size, local_work_size,
			waitListSz, waitList, &event );
		CHECKERRD(errcode, "failed to enqueue NDRangeKernel");
		return Event(event);
	}

	/** if block_size == 0, then local_work_size = NULL (autoselect) */
	Event enqueueKernel( Kernel& kernel, size_t work_size, size_t block_size = 0,
		cl_uint waitListSz = 0, const cl_event *waitList = NULL )
	{
		const size_t * bsptr = NULL;
		if ( block_size > 0 ) bsptr = &block_size;
	
		return enqueueNDRangeKernel( kernel, 1, &work_size, bsptr, waitListSz, waitList );
	}

	Event enqueueTask( Kernel& kernel, cl_uint waitListSz = 0, const cl_event *waitList = NULL )
	{
		cl_event event;
		const cl_int errcode = clEnqueueTask( handle, kernel.getHandle(), waitListSz, waitList, &event);
		CHECKERRD(errcode, "failed to enqueue task");
		return Event(event);
	}
};

class Pin
{
protected:
	CommandQueue& com;
	Buffer& buffer;
	void *hostPtr;
	
public:
	Pin( CommandQueue& queue, Buffer& buf, cl_map_flags flags )
	:
		com(queue), buffer(buf)
	{
		hostPtr = com.enqueueMapBuffer( buffer, flags );
	}
	
	virtual ~Pin()
	{
		com.enqueueUnmapBuffer( buffer, hostPtr );
	}
	
	inline void *getHostPtr() { return hostPtr; }
};

class PinnedBuffer : protected Buffer, public Pin
{
public:
	PinnedBuffer( CommandQueue& queue, size_t size,
		cl_mem_flags mem_flags = CL_MEM_READ_WRITE,
		cl_map_flags map_flags = CL_MAP_READ | CL_MAP_WRITE )
	:
		Buffer( queue.getContext(), mem_flags | CL_MEM_ALLOC_HOST_PTR, size ),
		Pin( queue, *this, map_flags )
	{}
};

class PinnedBufferPair : public PinnedBuffer
{
private:
	Buffer deviceBuffer;
	
public:
	PinnedBufferPair( CommandQueue& queue, size_t size,
		cl_mem_flags mem_flags = CL_MEM_READ_WRITE,
		cl_map_flags map_flags = CL_MAP_READ | CL_MAP_WRITE )
	:
		PinnedBuffer( queue, size, mem_flags, map_flags ),
		deviceBuffer( queue.getContext(), mem_flags, size )
	{}
	
	Buffer &getDeviceBuffer() { return deviceBuffer; }

	
	Event transferToDevice(
		cl_bool blocking = CL_TRUE,
		size_t offset = 0, size_t bytes = 0,
		cl_uint waitListSz = 0, const cl_event * waitList = NULL )
	{
		return com.enqueueWriteBuffer( deviceBuffer, hostPtr,
			blocking, offset, bytes, waitListSz, waitList );
	}
	
	Event transferFromDevice(
		cl_bool blocking = CL_TRUE,
		size_t offset = 0, size_t bytes = 0,
		cl_uint waitListSz = 0, const cl_event * waitList = NULL)
	{
		return com.enqueueReadBuffer( deviceBuffer, hostPtr,
			blocking, offset, bytes, waitListSz, waitList );
	}
	
	using Buffer::size;
};

class SingleGPUSetup
	: protected DeviceList, public Context, public CommandQueue
{
public:
	SingleGPUSetup()
	:
		DeviceList( selectPlatform(), CL_DEVICE_TYPE_CPU ),
		Context( *static_cast<CL::DeviceList*>(this) ),
		CommandQueue( *this, getFirstDevice() )
	{}
};

class SingleFileSingleGPUSetup
	: public SingleGPUSetup, public Program 
{
public:
	SingleFileSingleGPUSetup(
		const char *filename, bool binary, const char *compilerFlags )
	{
		if ( binary )
		{
			if ( DeviceList::getNumDevices() != 1 )
				throw std::runtime_error( "Single binary for multiple devices" );
			
			std::vector<std::string> bins;
			bins.push_back( fileToString(filename) );
			
			Program::buildFromBinaries( *this, *this, bins, compilerFlags );
		}
		else
		{
			Program::buildFromSource( *this, fileToString(filename), compilerFlags );
		}
		CHECKERR( clUnloadCompiler() ); // unload compiler hint
	}
};

class SingleSourceSingleGPUSetup : public SingleFileSingleGPUSetup
{
public:
	SingleSourceSingleGPUSetup(
		const char *srcFilename, const char *compilerFlags = NULL )
	:
	SingleFileSingleGPUSetup( srcFilename, false, compilerFlags )
	{}
};

class SingleBinarySingleGPUSetup : public SingleFileSingleGPUSetup
{
public:
	SingleBinarySingleGPUSetup(
		const char *filename, const char *compilerFlags = NULL )
	:
	SingleFileSingleGPUSetup( filename, true, compilerFlags )
	{}
};

}
#endif
