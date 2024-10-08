from FODLego.globaldata import GlobalData
import numpy as np
from numpy import sqrt, array
from numpy.linalg import norm
from FODLego.FOD import CFOD

class FODShell:
    """
    This is the parent class to a variety of core FOD structures.

    """
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
        super().__init__(atom,'point', [CFOD(atom, atom.mPos)])

class Tetra(FODShell):
    """
    Tetrahedron Class: FOD Geometry corresponding to sp3 "hybridization' geometry.
    Roadmap: There will  be different functions to create compound transformations of FODs (e.g. the base, or peak
    of the tetrahedron), and to rotate them in the proper direction as well.
    """
    def __init__(self, atom, core_amount: int) -> None:
        # Scale the FODs according to radius
        s = GlobalData.mRadii[core_amount][atom.mZ]
        # The geometry of a tetrahedron in a unit circle
        super().__init__(atom, 'tetra', [
            CFOD(atom, atom.mPos + s*array([0.0,0.0,1.0])),
            CFOD(atom, atom.mPos + s*array([sqrt(8/9), 0.0, -1/3])),
            CFOD(atom, atom.mPos + s*array([-sqrt(2/9),sqrt(2/3), -1/3])),
            CFOD(atom, atom.mPos + s*array([-sqrt(2/9),-sqrt(2/3), -1/3]))
            ])

        # Descriptive Stats. of Predicted FODs
        self.mPred_u_R = 0.0
        self.mPred_s2_R = 0.0
        # Descriptive Stats. of Target FODs
        self.mTarget_u_R = 0.0
        self.mTarget_s2_R = 0.0

    #Class Methods
    def RotateTetra(self):
        pass

class Triaug(FODShell):
    """
    This shape is created for SPD hybridization. It is not part of the scope of Angel's 2024 MS Thesis,
    However, its implementation could go here.
    """
    pass
