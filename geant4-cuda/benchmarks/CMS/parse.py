
#
# Strips down & modifies cms.gdml to include only supported types of
# G4 volumes. Writes results to cms.txt
#

import math
from xml.sax.handler import ContentHandler
import xml.sax
from affinetransform import AffineTransform, Rotation, Translation

INPUTFILE = "cms.gdml"
OUTPUTFILE = "cms.txt"

def readxyz( attrs ):
	x,y,z = [0,0,0]
	if "x" in attrs.keys(): x = float(attrs.get("x"))
	if "y" in attrs.keys(): y = float(attrs.get("y"))
	if "z" in attrs.keys(): z = float(attrs.get("z"))
	return [ x, y, z ]

class ThingWithUnits:
	def __init__( self, attrs ):
		self.aunit = attrs.get("aunit")
		self.lunit = attrs.get("lunit")
		self.name = attrs.get("name")
		if "lunit" in attrs.keys() and self.lunit != "mm":
			raise RuntimeError("Unsupported units")
		
	def __str__( self ):
		return self.__class__.__name__ + " " + self.name
		
	def isValid( self ):
		return True
		
	def fitsInside( self, other ):
		if self.outerRadius() > other.outerRadius(): return False
		if self.innerRadius() < other.innerRadius():
			if self.outerRadius()*2.0 > other.outerRadius()-other.innerRadius():
				return False
		return True
		
	def innerRadius( self ):
		return 0.0
		
class Tube(ThingWithUnits):
	def __init__(self, attrs):
		ThingWithUnits.__init__(self,attrs)
		self.deltaphi = attrs.get("deltaphi")
		self.startphi = attrs.get("startphi")
		self.rmax = attrs.get("rmax")
		self.rmin = attrs.get("rmin") 
		self.z = attrs.get("z")
		
	def write( self, stream ):
		stream.write(self.__class__.__name__+"\n")
		stream.write(self.name+"\n")
		stream.write(self.rmin+"\n")
		stream.write(self.rmax+"\n")
		stream.write(self.z+"\n")
		stream.write(self.startphi+"\n")
		stream.write(self.deltaphi+"\n")
		
	def outerRadius( self ):
		return float(self.rmax)
		
	def innerRadius( self ):
		return float(self.rmin)
		
	def isValid( self ):
		#if self.deltaphi != '360' or self.startphi != '0': return False
		#if float(self.rmin)/float(self.rmax) < 0.1: return False
		#if float(self.rmin)/float(self.rmax) > 0.5: return False
		#if float(self.z)/float(self.rmax) < 0.1: return False
		return True
		
class Box(ThingWithUnits):
	def __init__(self, attrs):
		ThingWithUnits.__init__(self,attrs)
		self.x, self.y, self.z = readxyz( attrs )
		
	def write( self, stream ):
		stream.write(self.__class__.__name__+"\n")
		stream.write(self.name+"\n")
		stream.write(str(self.x)+"\n")
		stream.write(str(self.y)+"\n")
		stream.write(str(self.z)+"\n")
		
	def outerRadius( self ):
		return max([self.x, self.y, self.z])*math.sqrt(2)

class Cone(ThingWithUnits):
	def __init__(self, attrs):
		ThingWithUnits.__init__(self,attrs)
		self.deltaphi = attrs.get("deltaphi")
		self.startphi = attrs.get("startphi")
		self.rmax1 = attrs.get("rmax1")
		self.rmin1 = attrs.get("rmin1")
		self.rmax2 = attrs.get("rmax2")
		self.rmin2 = attrs.get("rmin2")
		self.z = attrs.get("z")
		
	def outerRadius( self ):
		return max(float(self.rmax1),float(self.rmax2))
	
	def innerRadius( self ):
		return min(float(self.rmin1),float(self.rmin2))
		
	def write( self, stream ):
		stream.write(self.__class__.__name__+"\n")
		stream.write(self.name+"\n")
		stream.write(self.rmin1+"\n")
		stream.write(self.rmax1+"\n")
		stream.write(self.rmin2+"\n")
		stream.write(self.rmax2+"\n")
		stream.write(self.z+"\n")
		stream.write(self.startphi+"\n")
		stream.write(self.deltaphi+"\n")

class Polycone(ThingWithUnits):
	
	CONVERT = False
	
	def __init__(self,attrs,polyconeZ):
		ThingWithUnits.__init__(self,attrs)
		self.zplanes = polyconeZ
		self.deltaphi = attrs.get("deltaphi")
		self.startphi = attrs.get("startphi")
		self.fix()
		
	def write( self, stream ):
		if self.CONVERT:
			self.convert2Cone().write(stream)
		else:
			stream.write(self.__class__.__name__+"\n")
			stream.write(self.name+"\n")
			stream.write(self.startphi+"\n")
			stream.write(self.deltaphi+"\n")
			stream.write(str(len(self.zplanes))+"\n")
			for p in self.zplanes:
				p.write(stream)
			
	def convert2Cone( self ):
		if len(self.zplanes)<2: raise RuntimeError("Invalid polycone")
		attr = {}
		attr["z"] = str(self.zplanes[1].z - self.zplanes[0].z)
		attr["rmin1"] = self.zplanes[0].rmin
		attr["rmax1"] = self.zplanes[0].rmax
		attr["rmin2"] = self.zplanes[1].rmin
		attr["rmax2"] = self.zplanes[1].rmax
		attr["name"] = self.name;
		attr["deltaphi"] = self.deltaphi;
		attr["startphi"] = self.startphi;
		return Cone( attr )
		
	def fix( self ):
		newplanes = [self.zplanes[0]]
		prevz = self.zplanes[0].z
		for p in self.zplanes[1:]:
			if p.z > prevz:
				newplanes.append(p)
			prevz = p.z
		
	def isValid( self ):
		#prevz = self.zplanes[0].z
		#for p in self.zplanes[1:]:
		#	if p.z <= prevz:
		#		return False
		#	prevz = p.z
		return True

class ZPlane:
	def __init__(self,attrs):
		self.z = float(attrs.get("z"))
		self.rmax = attrs.get("rmax") 
		self.rmin = attrs.get("rmin")
		
	def write( self, stream ):
		stream.write(self.__class__.__name__+"\n")
		stream.write(self.rmax+"\n")
		stream.write(self.rmin+"\n")
		stream.write(str(self.z)+"\n")
		
class LogicalVolume:
	
	def __init__(self,name):
		self.daughters = []
		self.material = 0
		self.solid = 0
		self.cumulativesolids = 0
		self.name = name
		
	def __str__(self):
		return "LogicalVolume "+self.name + " : " + str(self.solid)
		
	def write( self, stream ):
		stream.write(self.name+"\n")
		stream.write(self.material+"\n")
		stream.write(self.solid.name+"\n")
		stream.write(str(len(self.daughters))+"\n")
		for d in self.daughters:
			d.write(stream)
	
class PhysicalVolume:
	def __init__(self):
		self.position = Translation()
		self.rotation = Rotation()
		self.logvolref = 0
		self.logvol = 0
	
	def getTransform( self ):
		return AffineTransform( self.rotation, self.position )
	
	def __str__(self):
		return "PhysicalVolume - "+str(self.logvol)
		
	def write( self, stream ):
		stream.write(self.__class__.__name__+"\n")
		stream.write(self.logvolref+"\n")
		self.rotation.write(stream)
		self.position.write(stream)

class Material:
	def __init__(self,name):
		self.name = name
		self.D = "0"
		self.atom = "0"
		
	def write( self, stream ):
		#LIMIT = 1
		#if float(self.D) < LIMIT: self.D = "0" #str(LIMIT)
		#if self.name != "AISI1018Steel0xb6a86c50": self.D = "0"
		stream.write(self.name+"\n")
		stream.write(self.D+"\n")
		stream.write(self.atom+"\n")
		

class GDMLHandler(ContentHandler):
	
	positions = {}
	rotations = {}
	materials = {}
	missedsolids = {}
	missedsolidkeys = set()
	solids = {}
	logvols = {}
	
	inPolycone = False
	inMaterial = False
	inPhysVol = False
	polyconeZ = []
	polyconeAttrs = 0
	
	curphysvol = 0
	curlogvol = 0
	worldvol = 0
	
	droppedsolids = 0
	totalsolids = 1
	
	knownSolids = set( ["box", "tube", "cone", "polycone", "reflectedSolid", "trap", "subtraction", "polyhedra"] )
	unrecognized = set()
	
	newlogvols = {}
	newsolids = {}
	newmaterials = {}

	def startElement(self, name, attrs):
	
		if name == "position":
			self.positions[ attrs.get("name") ] = Translation.Create(*readxyz( attrs ))
		elif name == "rotation":
			self.rotations[ attrs.get("name") ] = Rotation.CreateEuler(*readxyz( attrs ))
			
		elif name == "material" or name == "element":
			self.inMaterial = True
			n = attrs.get("name")
			self.curmaterial = Material(n)
			self.materials[ n ] = self.curmaterial
		elif name == "D" and self.inMaterial:
			self.curmaterial.D = attrs.get("value")
		elif name == "atom" and self.inMaterial:
			self.curmaterial.atom = attrs.get("value")
			
		elif name == "tube":
			self.solids[ attrs.get("name") ] = Tube( attrs )
		elif name == "box":
			self.solids[ attrs.get("name") ] = Box( attrs )
		elif name == "cone":
			self.solids[ attrs.get("name") ] = Cone( attrs )
			
		elif name == "polycone":
			self.inPolycone = True
			self.polyconeZ = []
			self.polyconeAttrs = attrs
		elif name == "zplane" and self.inPolycone:
			self.polyconeZ.append( ZPlane( attrs ) )
			
		elif name in self.knownSolids:
			if name in self.missedsolids: self.missedsolids[name] += 1
			else: self.missedsolids[name] = 1
			self.missedsolidkeys.add( attrs.get("name") )
			
		elif name == "volume":
			n = attrs.get("name")
			self.curlogvol = LogicalVolume( n )
			self.logvols[ n ] = self.curlogvol
			
		elif name == "physvol":
			self.inPhysVol = True
			self.curphysvol = PhysicalVolume()
			self.curlogvol.daughters.append( self.curphysvol )
			
		elif name == "materialref":
			matref = attrs.get("ref")
			if matref not in self.materials:
				raise RuntimeError("Invalid material "+matref)
			self.curlogvol.material = matref
		elif name == "solidref":
			solidref = attrs.get("ref")
			if solidref in self.solids:
				self.curlogvol.solid = self.solids[solidref]
			elif solidref not in self.missedsolidkeys:
				raise RuntimeError("Invalid solid key "+solidref)
				
		elif self.inPhysVol:
			if name == "positionref":
				self.curphysvol.position = self.positions[ attrs.get("ref") ]
			if name == "rotationref":
				self.curphysvol.rotation = self.rotations[ attrs.get("ref") ]
			if name == "volumeref":
				self.curphysvol.logvolref = attrs.get("ref")
			
		elif name == "world":
			self.worldvol = self.logvols[ attrs.get("ref") ]
		else:
			self.unrecognized.add( name )
				
	def endElement(self, name):
		if name == "physvol":
			self.inPhysVol = False
			if self.curphysvol.logvolref == 0:
				raise RuntimeError("No logical volume in physvol")
		elif name == "material" or name == "element":
			self.inMaterial = False
		elif name == "polycone":
			self.inPolycone = False
			self.solids[ self.polyconeAttrs.get("name") ] = Polycone( self.polyconeAttrs, self.polyconeZ )	


handler = GDMLHandler()
saxparser = xml.sax.make_parser()
saxparser.setContentHandler( handler )

datasource = open(INPUTFILE, "r")
saxparser.parse( datasource )

print "Total imported solids:", len(handler.solids)
print "Unsupported solids:", handler.missedsolids
print "Unhandled XML-elements:", handler.unrecognized


def isValidVolume( root ):
	
	enabledMaterials = set()
	#enabledMaterials.add("AISI1018Steel0xb6a86c50")
	
	disabledSolids = set()
	disabledSolids.add("Polycone")
	#disabledSolids.add("Cone")
	#disabledSolids.add("Tube")
	#disabledSolids.add("Box")

	disabledMaterials = set()
	#disabledMaterials.add("AISI1018Steel0xb6a86c50")
	
	validVolume = True
	
	if root.solid == 0:
		validVolume = False
	else:
		validVolume = root.solid.isValid()
	
	if root.solid.__class__.__name__ in disabledSolids: validVolume = False
	if len(enabledMaterials) and root.material not in enabledMaterials: validVolume = False
	if root.material in disabledMaterials: validVolume = False

	return validVolume
	

def traverseAndRefactor( handler, root, validroot, transform ):
	
	valid = isValidVolume(root)
	if validroot == 0: valid = True
	#elif valid:
	#	valid = root.solid.fitsInside(validroot.solid)
	#	if not valid: handler.droppednofit += 1
	
	if valid:
		transform = AffineTransform()
		validroot = root
	
	daugh = root.daughters
	root.daughters = []
	
	for d in daugh:
		if d.logvol == 0:
			d.logvol = handler.logvols[ d.logvolref ]
			
		oldtform = AffineTransform.Create(d.rotation, d.position)
		newtform = transform.mult( oldtform )
		if traverseAndRefactor( handler, d.logvol, validroot, newtform ):
			d.rotation = newtform.rot
			d.position = newtform.trans
			validroot.daughters.append( d )
			handler.totalsolids += 1
		else: handler.droppedsolids += 1
		
	return valid

def traversephys( root, handler, recdepth ):
	
	root.logvol = handler.logvols[ root.logvolref ]
		
	n = traverselog( root.logvol, handler, recdepth+1 )
	if n == 0: root.logvol = 0
	return n

def traverselog( root, validroot, handler, recdepth = 0 ):

	#print " "*recdepth, "Entering", root, "("+str(len(root.daughters))+" daughters)"

	if root.cumulativesolids > 0:
		if root.solid == 0: return 0
		return root.cumulativesolids

	newdaughters = []
	
	for d in root.daughters:
		n = traversephys( d, handler, recdepth )
		if n > 0: newdaughters.append(d)		
		root.cumulativesolids += n
	
	root.daughters = newdaughters
	
	if isValidVolume(root):
		root.cumulativesolids += 1
		return root.cumulativesolids
	else:
		handler.droppedsolids += root.cumulativesolids+1
		return 0

def traverseagain( root, handler ):
	
	handler.newlogvols[ root.name ] = root
	handler.newsolids[ root.solid.name ] = root.solid;
	handler.newmaterials[ root.material ] = handler.materials[ root.material ]
	
	for d in root.daughters:
		traverseagain( d.logvol, handler )

#if traverselog( handler.worldvol, handler ) == 0:
#	raise RuntimeError("No world volume left")

traverseAndRefactor( handler, handler.worldvol, 0, 0 )	
traverseagain( handler.worldvol, handler )

print "Dropped volumes:", handler.droppedsolids
print "Total volumes:", handler.totalsolids
#print "Did not fit:", handler.droppednofit

output = open(OUTPUTFILE, 'w')

output.write(str(len(handler.newmaterials))+"\n")
for m in handler.newmaterials.values():
	m.write(output)
	
output.write(str(len(handler.newsolids))+"\n")
for s in handler.newsolids.values():
	s.write(output)
	
output.write(str(len(handler.newlogvols))+"\n")
handler.worldvol.write(output)
for v in handler.newlogvols.values():
	if v != handler.worldvol:
		v.write(output)
