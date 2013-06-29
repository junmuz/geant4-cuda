
/** Parses some geometry of the CMS detector from cms.txt */

#ifndef MYPARSER_HPP
#define MYPARSER_HPP

#include "../../geometry_common.hpp"

#include <fstream>
#include <map>
#include <list>

class CMSModel : public BasicGeometry
{
public:
	CMSModel( int maxsz = 1024*1000000 ) // ~1GB
	 :
		BasicGeometry( maxsz ),
		scalefactor( 1E-3 ) // Scale in meters!
	{}
	
	double getScale() const
	{
		return 28600*scalefactor;
	}

	CameraParameters getCamera() const
	{
		CameraParameters p;
		p.heading = 20;
		p.pitch = 60;
		p.dist = 0.22*50000*scalefactor;
		p.yfov = 80;
		p.target_z = 0.04*50000*scalefactor;
		return p;
	}

private:

	const G4double scalefactor;

	G4VPhysicalVolume *worldVolume;
	G4LogicalVolume *logWorld;
	std::map<std::string, StubMaterial*> materials;
	std::map<std::string, G4VSolid*> solids;
	std::map<std::string, G4LogicalVolume*> logvols;
	std::map<G4VPhysicalVolume*, std::string> phystologref;
	std::list<G4VPhysicalVolume*> physvols;

	void create_impl()
	{
		printf ("Entered in CMSModel.hpp\n");
		const char *filename = "benchmarks/CMS/cms.txt";
		std::ifstream stream(filename);
		
		worldVolume = reserveThing<G4VPhysicalVolume>();
		
		readMaterials( stream );
		#ifdef VERBOSE
		std::cerr << "Imported " << materials.size() << " materials\n";
		#endif
		readSolids( stream );
		#ifdef VERBOSE
		std::cerr << "Imported " << solids.size() << " solids\n";
		#endif
		readStructure( stream );
		#ifdef VERBOSE
		std::cerr << "Imported " << logvols.size() << " logical and "
			<< physvols.size()+1 << " physical volumes\n";
		#endif
		
		for ( std::list<G4VPhysicalVolume*>::const_iterator itr = physvols.begin();
			itr != physvols.end(); ++itr )
		{
			addPhysicalVolumePointers( *itr );
		}

		for ( std::map<std::string,G4LogicalVolume*>::const_iterator itr = logvols.begin();
			itr != logvols.end(); ++itr )
		{
			#ifdef ENABLE_VOXEL_NAVIGATION
			const int VOX_THRESHOLD = G4VoxelHeader_GetMinVoxelVolumesLevel1();
			const int nd = G4LogicalVolume_GetNoDaughters(itr->second);
			if ( nd >= VOX_THRESHOLD )
			{
				#ifdef VERBOSE
				//std::cerr << itr->first << ": " << nd << " daughters, voxelizing...\n";
				#endif
				createVoxels( itr->second );
			}
			#endif
			addLogicalVolumePointers( itr->second );
		}
	}
	
	template <class T> static T get( std::istream& in )
	{
		T obj;
		in >> obj;
		if (!in) throw std::runtime_error("Failed to read input");
		return obj;
	}
	
	static G4double degToRad( G4double d )
	{
		return d / 360.0 * twopi;
	}

	G4ThreeVector readVec( std::istream& in )
	{
		const G4double x = get<G4double>(in)*scalefactor;
		const G4double y = get<G4double>(in)*scalefactor;
		const G4double z = get<G4double>(in)*scalefactor;
		return G4ThreeVector_create(x,y,z);
	}
	
	G4RotationMatrix readRot( std::istream& in )
	{
		G4double elem[3][3];
		for (int i=0; i<3; ++i)
			for (int j=0; j<3; ++j)
				elem[i][j] = get<G4double>(in);
				
		return G4RotationMatrix_create_elements(
			elem[0][0], elem[0][1], elem[0][2],
			elem[1][0], elem[1][1], elem[1][2],
			elem[2][0], elem[2][1], elem[2][2] );
	}
	
	void readMaterials( std::istream& in )
	{
		const int numberOfMaterials = get<int>( in );
		
		// Possibly override material settings

		//#define IDENTICAL_MATERIALS

		const enum { OVERRIDE, INCLUDE, EXCLUDE, NONE }
			overrideMaterials =
		#if defined(PHYSICS) && defined(IDENTICAL_MATERIALS)
				INCLUDE;
		#else
				NONE;
		#endif
			
		const G4double default_property
		#if defined(PHYSICS) && defined(IDENTICAL_MATERIALS)
			 = 1.0;
		#else
			 = 0.0;
		#endif

		// The densities are given as g/cm^3, lengths are already scaled to meters
		// multiply by 1000 to get the densities in kg/m^3
		const G4double densityMultiplier = 1000.0;

		std::map<std::string, G4double> materialProperties;
		// Leave empty, no overrides
		
		StubMaterial *matarray = reserveNThings<StubMaterial>(numberOfMaterials);
		
		for ( int i=0; i<numberOfMaterials; ++i )
		{
			const std::string materialName = get<std::string>(in);
			const G4double D = get<G4double>(in);
			const G4double atom = get<G4double>(in);
			G4double property = D + 0*atom;	// Take "density" value from "D" value of material
			//property = pow(property,1.0/4);
			
			if ( materialProperties.count( materialName ) )
			{
				switch ( overrideMaterials )
				{
				case EXCLUDE: property = default_property; break;
				case OVERRIDE: property = materialProperties[materialName]; break;
				default: break;
				}
			}
			else
			{
				switch ( overrideMaterials )
				{
				case INCLUDE: property = default_property; break;
				default: break;
				}
			}
				
			matarray[i].property = property*densityMultiplier;
			
			if (materials.count(materialName))
				throw std::runtime_error("Duplicate material key "+materialName);
			
			materials[materialName] = matarray+i;
		}
	}
	
	void readSolids( std::istream& in )
	{
		const int numberOfSolids = get<int>( in );
		
		for ( int i=0; i<numberOfSolids; ++i )
		{
			const std::string type = get<std::string>(in);
			const std::string name = get<std::string>(in);
			
			// std::cerr << "Importing " << type << " " << name << "\n";
			
			G4VSolid *s = NULL;
			
			if ( type == "Box" )
			{
				const G4ThreeVector xyz = readVec(in);
				G4Box *box = reserveThing<G4Box>();
				G4Box_ctor( box, xyz.x, xyz.y, xyz.z );
				s = (G4VSolid*)box;
			}
			else if ( type == "Tube" )
			{
				const G4double rmin = get<G4double>(in)*scalefactor;
				const G4double rmax = get<G4double>(in)*scalefactor;
				const G4double z = get<G4double>(in)*scalefactor;
				const G4double sphi = degToRad(get<G4double>(in));
				const G4double dphi = degToRad(get<G4double>(in));
				
				#ifndef ENABLE_SLICED_TUBS
				if (std::fabs(dphi-twopi) > 0.1)
				{
					G4Cons *tubs = reserveThing<G4Cons>();
					s = (G4VSolid*)tubs;
					G4Cons_ctor( tubs, rmin, rmax, rmin, rmax, z, sphi, dphi );
				}
				else
				{
				#endif
					G4Tubs *tubs = reserveThing<G4Tubs>();
					s = (G4VSolid*)tubs;
					G4Tubs_ctor( tubs, rmin, rmax, z, sphi, dphi );
				#ifndef ENABLE_SLICED_TUBS
				}
				#endif
				
			}
			else if ( type == "Cone" )
			{
				G4Cons *cons = reserveThing<G4Cons>();
				s = (G4VSolid*)cons;
				
				const G4double rmin1 = get<G4double>(in)*scalefactor;
				const G4double rmax1 = get<G4double>(in)*scalefactor;
				const G4double rmin2 = get<G4double>(in)*scalefactor;
				const G4double rmax2 = get<G4double>(in)*scalefactor;
				const G4double z = get<G4double>(in)*scalefactor;
				const G4double sphi = degToRad(get<G4double>(in));
				const G4double dphi = degToRad(get<G4double>(in));
					
				G4Cons_ctor( cons, rmin1, rmax1, rmin2, rmax2, z, sphi, dphi );
			}
#ifdef ENABLE_G4POLYCONE
			else if ( type == "Polycone" )
			{
				const G4double sphi = degToRad(get<G4double>(in));
				const G4double dphi = degToRad(get<G4double>(in));
				const int nplanes = get<int>( in );
				
				G4double *zarr = reserveNThings<G4double>(nplanes);
				G4double *rminarr = reserveNThings<G4double>(nplanes);
				G4double *rmaxarr = reserveNThings<G4double>(nplanes);
				
				for ( int j=0; j<nplanes; ++j )
				{
					if ( get<std::string>(in) != "ZPlane" )
						throw std::runtime_error( "Invalid input" );
						
					const G4double r2 = get<G4double>(in)*scalefactor;
					const G4double r1 = get<G4double>(in)*scalefactor;
					const G4double z = get<G4double>(in)*scalefactor;
					
					zarr[j] = z;
					rminarr[j] = r1;
					rmaxarr[j] = r2;
				}
				
				G4PolyCone *cons = reserveThing<G4PolyCone>();
				s = (G4VSolid*)cons;
				
				G4PolyCone_ctor( cons, nplanes, zarr, rminarr, rmaxarr, sphi, dphi );
				addPolyConePointers( cons );
			}
#endif
			else
				throw std::runtime_error("Unsupported object type: "+type);
				
			assert( s != NULL );
			
			if (solids.count(name))
				throw std::runtime_error("Duplicate solid key "+name);
				
			solids[ name ] = s;
		}
	}
	
	
	void readStructure( std::istream &in )
	{
		const int numberOfLogicalVolumes = get<int>(in);
		
		for ( int i=0; i<numberOfLogicalVolumes; ++i )
		{
			const std::string volName = get<std::string>(in);
			const std::string materialKey = get<std::string>(in);
			const std::string solidKey = get<std::string>(in);
			
			G4LogicalVolume *logVol = reserveThing<G4LogicalVolume>();		
			G4LogicalVolume_ctor( logVol,
				solids.find(solidKey)->second,
				materials.find(materialKey)->second );

			if (logvols.count(volName))
				throw std::runtime_error("Duplicate volume key "+volName);
				
			logvols[volName] = logVol;
			
			if ( i == 0 )
			{
				logWorld = logVol;
				const G4RotationMatrix idRot = G4RotationMatrix_create_id();
				const G4ThreeVector zeroTrans = G4ThreeVector_create(0,0,0);
				G4VPhysicalVolume_ctor( worldVolume, idRot, zeroTrans, logVol );
				addPhysicalVolumePointers( worldVolume );
			}
			
			const int numDaughters = get<int>(in);
			if ( numDaughters > 0 )
			{
				G4VPhysicalVolume *pdarr =
					reserveNThings<G4VPhysicalVolume>(numDaughters);
				
				for ( int j=0; j<numDaughters; ++j )
				{
					if ( get<std::string>(in) != "PhysicalVolume" )
						throw std::runtime_error("Invalid input");
						// ignore 'PhysicalVolume'
						
					G4VPhysicalVolume *physvol = pdarr+j;
					phystologref[physvol] = get<std::string>(in);
					
					const G4RotationMatrix rot = readRot(in);
					const G4ThreeVector trans = readVec(in);
					
					G4VPhysicalVolume_ctor( physvol, rot, trans, NULL );
					G4VPhysicalVolume_SetMotherLogical( physvol, logVol );
					physvols.push_back(physvol);
					addLogicalVolumeDaughter( logVol, physvol );
				}
			}
		}
		
		for ( std::list<G4VPhysicalVolume*>::const_iterator itr = physvols.begin();
			itr != physvols.end(); ++itr )
		{
			G4VPhysicalVolume_SetLogicalVolume( *itr,
				logvols.find(phystologref.find(*itr)->second)->second );
		}
	}
};

#endif
