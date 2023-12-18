from globaldata import GlobalData
import numpy as np

class FODShell:
    def __init__(self, shape, fods):
        self.mShape = shape 
        self.mfods = np.array(fods)
    def __str__(self):
        return self.mShape

class Point(FODShell):
    def __init__(self):
        super().__init__('point', [0.0,0.0,0.0])

class Tetra(FODShell):
    """
    Tetrahedron Class: FOD Geometry corresponding to sp3 "hybridization' geometry.
    Roadmap: There will  be different functions to create compound transformations of FODs (e.g. the base, or peak
    of the tetrahedron), and to rotate them in the proper direction as well.
    """
    def __init__(self, core_amount: int , z: int) -> None:
        super().__init__('tetra', GlobalData.mTetraGeo)
        self.mfods *= GlobalData.mRadii[core_amount][z]
    #Class Methods
    def CreateTetra(self):
        """
        This method arranges 4 FOD points in a tetrahedral form. Depending on the bonding, and free electrons,
        the direction vector will be different. 
        """
        pass

    def RotateTetra(self):
        pass