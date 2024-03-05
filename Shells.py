from globaldata import GlobalData
import numpy as np
from numpy import sqrt, array
from numpy.linalg import norm
from FOD import FOD

class FODShell:
    def __init__(self, atom, shape, fods):
        self.mAtom = atom
        self.mShape = shape 
        self.mfods = fods

    def __str__(self):
        return self.mShape
    
    def GetLastAtomRadius(self):
        return norm(self.mfods[-1].mPos - self.mAtom.mPos)

class Point(FODShell):
    def __init__(self, atom):
        super().__init__(atom,'point', [FOD(atom.mPos)])

class Tetra(FODShell):
    """
    Tetrahedron Class: FOD Geometry corresponding to sp3 "hybridization' geometry.
    Roadmap: There will  be different functions to create compound transformations of FODs (e.g. the base, or peak
    of the tetrahedron), and to rotate them in the proper direction as well.
    """
    def __init__(self, atom, core_amount: int) -> None:
        # The geometry of a tetrahedron in a unit circle
        super().__init__(atom, 'tetra', [
            FOD(array([0.0,0.0,1.0])),
            FOD(array([sqrt(8/9), 0.0, -1/3])),
            FOD(array([-sqrt(2/9),sqrt(2/3), -1/3])),
            FOD(array([-sqrt(2/9),-sqrt(2/3), -1/3]))
            ])
        # Scale the FODs according to radius
        scale = GlobalData.mRadii[core_amount][atom.mZ]
        [fods*scale for fods in self.mfods]
        [fods+atom.mPos for fods in self.mfods]

    #Class Methods
    def RotateTetra(self):
        pass