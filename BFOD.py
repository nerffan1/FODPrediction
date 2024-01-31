from Funcs import *
from FOD import *
from ElementaryClasses import *
from globaldata import *

class BFOD(FOD):
    def __init__(self, boldAt: Atom, meekAt: Atom):
        super().__init__()
        self.mBold = boldAt
        self.mMeek = meekAt
        self.mMeekR = -1.0
        self.mBoldR = -1.0
        self.mBoldPortion = -1.0
        self.mBondDir = meekAt.mPos - boldAt.mPos  # Always in direction away fromBold atom
        self.mBondDist = np.linalg.norm(self.mBondDir)
        self.mBondDir = normalize(self.mBondDir)  

    def Calc_AxisBoldPortion(self, Zbold:int, Zmeek:int) -> float:
            """
            Finds the portion (from 0 to 1) of the bonding distance that the Bold atom covers.
            For example, the FODs might be placed at the 0.33 mark of the line starting
            from the bold atom, ending on the meek atom.   
 
            Parameters
            ----------
            - Zbold : The atomic number of the dominant (bold) atom.
            - Zmeek : The atomic number of the weaker (meek) atom.      
            
            Returns
            -------
            Float
            The portion of the bonding axis covered by the dominant atom. 
            """
            if Zbold == Zmeek:
                return 0.5 #Equal atoms, then return the midpoint
            else:
                Zbold = 0.4 if Zbold == 1 else Zbold
                r = sqrt(Zmeek/Zbold)
                g = r/(1+r)
            # The portion should always be under 0.5 since it is the closer to the
            # dominant atom.
            if g < 0.5:
                return g
            else:
                return (1-g)

class SBFOD(BFOD):
    """
    This is the Single Bonding FOD (SBFOD) class
    """
    def __init__(self, bold: Atom, meek: Atom):
        super().__init__(bold,meek)
        self.mBoldPortion = self.Calc_AxisBoldPortion(bold.mZ, meek.mZ)
        self.DetermineParamenters()

    def DetermineParamenters(self):
        """
        Add the single BFOD along the axis. 
        """
        bfod = self.mBondDir*self.mBoldPortion*self.mBondDist
        self.mPos = self.mBold.mPos + bfod
        self.mBoldR = self.mBondDist*self.mBoldPortion
        self.mMeekR = self.mBondDist*(1-self.mBoldPortion)

class DBFOD(BFOD):
    def __init__(self, bold: Atom, meek: Atom, heightdir: np.ndarray):
        super().__init__(bold,meek)
        self.mHeight = heightdir
        self.mBoldAngle = np.deg2rad(54.735) 
        self.mMeekAngle = 0.0
        self.DetermineParamenters()

    def GetHeight(self) -> float:
        """
        Return radial distance from interatomic (bonding) axis. 
        For heterogenous atoms....
        """
        if self.mBold.mZ == self.mMeek.mZ:
            return self.mBold.GetMonoCovalEdge()/2
        else:
            rad = self.mBold.GetMonoCovalRad()
            return np.sin(self.mBoldAngle)*rad

    def GetBondAxProj(self) -> float:
        """
        This function determines the ration between the projection of the FOD-ATOM 
        and the Bonding Axis.
        """
        if self.mBold.mZ == self.mMeek.mZ:
            return 0.5
        else:
            proj = np.cos(self.mBoldAngle)*self.mBold.GetMonoCovalRad()
            return proj/self.mBondDist

    def DetermineParamenters(self):
        """
        Determine the Atom-Atom-FOD angles.
        """
        # Set BondProjection
        self.mBoldPortion = self.GetBondAxProj()
        # Measure Meek Angle
        toFOD = self.mPos - self.mMeek.mPos
        self.mMeekAngle = AngleBetween(self.mBondDir,toFOD)
        # Set FOD Position 
        delta_bond = self.mBondDir*self.mBondDist*self.mBoldPortion
        delta_height = self.mHeight*self.GetHeight()
        self.mPos = self.mBold.mPos + delta_bond + delta_height

    def Duplicate(self):
        """
        Create the pair FOD. It will contain the same parameters (resulting from the naive
        solution) as the current (self) FOD, but will change the direction of the FOD.
        """    

class TBFOD(BFOD):
    def __init__(self,bold: Atom, meek: Atom):
        super().__init__(bold,meek)