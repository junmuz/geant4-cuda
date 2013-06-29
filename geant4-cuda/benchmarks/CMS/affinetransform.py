
#
# An affine transform class for combining transformations to refactor
# a GDML/Geant4 geometry hierarchy
#

import math

class Rotation:
	
	def __init__(self):
		self.elem = [[1,0,0],[0,1,0],[0,0,1]]
		
	@staticmethod
	def CreateEuler( rx, ry, rz ):
		m = Rotation()
		m.rotateX(math.radians(rx))
		m.rotateY(math.radians(ry))
		m.rotateZ(math.radians(rz))
		return m
		
	def apply(self, vec):
		new = Translation.Create(0,0,0)
		for i in range(3):
			for j in range(3):
				new.elem[i] += self.elem[i][j]*vec.elem[j]
		return new
		
	def mult(self, other):
		new = Rotation()
		new.elem = [[0,0,0],[0,0,0],[0,0,0]]
		for k in range(3):
			for i in range(3):
				for j in range(3):
					new.elem[i][k] += self.elem[i][j]*other.elem[j][k]
		return new
		
	def isUnity( self ):
		return all([self.elem[i][i] == 1 for i in range(3)])
		
	def rotateX( self, a ):
		c = math.cos(a)
		s = math.sin(a)
		x = self.elem[1][0]
		y = self.elem[1][1]
		z = self.elem[1][2]
		
		self.elem[1][0] = c*x - s*self.elem[2][0]
		self.elem[1][1] = c*y - s*self.elem[2][1]
		self.elem[1][2] = c*z - s*self.elem[2][2]
		self.elem[2][0] = s*x + c*self.elem[2][0]
		self.elem[2][1] = s*y + c*self.elem[2][1]
		self.elem[2][2] = s*z + c*self.elem[2][2]
		
	def rotateY( self, a ):
		c = math.cos(a)
		s = math.sin(a)
		x = self.elem[2][0]
		y = self.elem[2][1]
		z = self.elem[2][2]
		
		self.elem[2][0] = c*x - s*self.elem[0][0]
		self.elem[2][1] = c*y - s*self.elem[0][1]
		self.elem[2][2] = c*z - s*self.elem[0][2]
		self.elem[0][0] = s*x + c*self.elem[0][0]
		self.elem[0][1] = s*y + c*self.elem[0][1]
		self.elem[0][2] = s*z + c*self.elem[0][2]
		
	def rotateZ( self, a ):
		c = math.cos(a)
		s = math.sin(a)
		x = self.elem[0][0]
		y = self.elem[0][1]
		z = self.elem[0][2]
		
		self.elem[0][0] = c*x - s*self.elem[1][0]
		self.elem[0][1] = c*y - s*self.elem[1][1]
		self.elem[0][2] = c*z - s*self.elem[1][2]
		self.elem[1][0] = s*x + c*self.elem[1][0]
		self.elem[1][1] = s*y + c*self.elem[1][1]
		self.elem[1][2] = s*z + c*self.elem[1][2]
	
	def write( self, s ):
		for r in self.elem:
			for e in r:
				s.write(str(e)+"\n")
		
class Translation:
	def __init__(self):
		self.elem = [0,0,0]
		
	def apply(self, vec):
		return Translation.Create(*[vec.elem[i]+self.elem[i] for i in range(3)])
		
	def mag2(self):
		return sum([x*x for x in self.elem])
		
	@staticmethod
	def Create( x, y, z ):
		t = Translation()
		t.elem = [x,y,z]
		return t
		
	def write( self, s ):
		for e in self.elem:
			s.write(str(e)+"\n")
		

class AffineTransform:
	
	def __init__(self):
		self.rot = Rotation()
		self.trans = Translation()
		
	def isUnity(self):
		return self.trans.mag2() == 0 and self.rot.isUnity()
		
	@staticmethod
	def Create( rot, trans ):
		t = AffineTransform()
		t.rot = rot
		t.trans = trans
		return t
		
	def apply(self, vec):
		return self.trans.apply(self.rot.apply(vec))
		
	def mult(self, other):
		return AffineTransform.Create(
			self.rot.mult( other.rot ),
			self.trans.apply( self.rot.apply( other.trans )))
