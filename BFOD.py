from Funcs import *
from FOD import *
from ElementaryClasses import *
from globaldata import *

class BFOD(FOD):
    """
    The BFOD class is in charge of encapsulating the reparametrization of Bonding FODs. This class describes the FODs using the angles, distances, heights, and ratios
    along the bonding axis in order to characterize the FODs. Some of these attributes of this class will tend toward zero in the SBFOD, but we know from some examples
    that they are not always zero (e.g. they might lie slightly off the bonding axis, so there might be an angle).
    """
    def __init__(self, boldAt: Atom, meekAt: Atom, target=None):
        super().__init__()
        # Atoms
        self.mBold = boldAt
        self.mMeek = meekAt
        # Angles
        self.mBoldAngle = 0.0
        self.mMeekAngle = 0.0
        # Vectors
        self.mBondDir = meekAt.mPos - boldAt.mPos  # Always in direction away fromBold atom
        self.mHeight = np.zeros(3)
        # Distances
        self.mBondDist = np.linalg.norm(self.mBondDir)
        self.mMeekR = 0.0
        self.mBoldR = 0.0
        # Misc 
        self.mBoldPortion = 0.0
        if target != None:
            self.mPos = target
            self.ReverseDetermination(target)

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

    def IsMonoatomic(self) -> bool:
        """
        This function returns whether or not we are dealing with a monoatomic bond
        """
        if self.mMeek.mZ == self.mBold.mZ:
            return True
        else:
            return False

    def ReverseDetermination(self, targetFOD: np.ndarray):
        """"
        This function uses the target FOD (e.g. an optimized FOD from FLOSIC) and returns the parameters
        """
        # Helper data
        bold2fod = targetFOD -  self.mBold.mPos
        meek2fod = targetFOD - self.mMeek.mPos
        # Angles
        self.mBoldAngle = AngleBetween(bold2fod, self.mBondDir)
        self.mMeekAngle = AngleBetween(meek2fod, -self.mBondDir)
        # Distances
        self.mMeekR = norm(meek2fod)
        self.mBoldR = norm(bold2fod)
        # Misc 
        self.mBoldPortion = np.cos(self.mBoldAngle)*self.mBoldR/self.mBondDist

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
        bfod = self.mBondDir*self.mBoldPortion
        self.mPos = self.mBold.mPos + bfod
        self.mBoldR = self.mBondDist*self.mBoldPortion
        self.mMeekR = self.mBondDist*(1-self.mBoldPortion)

class DBFOD(BFOD):
    def __init__(self, bold: Atom, meek: Atom, heightdir: np.ndarray):
        super().__init__(bold,meek)
        self.mHeight = heightdir
        self.mBoldAngle = np.deg2rad(54.735) 
        self.DetermineParameters()
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

    def DetermineParameters(self):
        """
        Determine the Atom-Atom-FOD angles.
        """
        # Set BondProjection
        self.mBoldPortion = self.GetBondAxProj()
        # Measure Meek Angle
        toFOD = self.mPos - self.mMeek.mPos
        self.mMeekAngle = AngleBetween(self.mBondDir,toFOD)
        # Set FOD Position 
        delta_bond = self.mBondDir*self.mBoldPortion
        delta_height = self.mHeight*self.GetHeight()
        # Finalize parameters
        self.mPos = self.mBold.mPos + delta_bond + delta_height

    def Duplicate(self):
        """
        Create the pair FOD. It will contain the same parameters (resulting from the naive
        solution) as the current (self) FOD, but will change the direction of the FOD.
        """    

class TBFOD(BFOD):
    def __init__(self, bold: Atom, meek: Atom, heightdir: np.ndarray):
        super().__init__(bold,meek)
        self.mHeight = heightdir 
        self.mBoldAngle = np.deg2rad(54.735) 
        self.DetermineParameters()

    def DetermineParameters(self):
        """
        This function determines the distance of the TBFOD away from the bonding
        axis. 
        """
        # Determine the location of the FOD
        rad = self.mBold.GetMonoCovalRad()
        a = self.mBold.GetMonoCovalEdge()
        c = self.mBondDist 
        # This section modifies the height depending on the nature of the monoatomic bond
        # 2nd period elements tend to tighten BFODs.
        if self.IsMonoatomic():
            if self.mBold.mPeriod == 2 & self.mMeek.mPeriod == 2:
                theta = np.arccos((c/2)/rad)
            elif self.mBold.mPeriod > 2 & self.mMeek.mPeriod > 2:
                theta = np.arctan(rad/(c/2)) 
        else:
            #theta = np.arcsin((a/2)/rad) 
            theta = np.deg2rad(54.735)
        dr = normalize(self.mBondDir)*rad*np.cos(theta) 
        dr += self.mHeight*rad*np.sin(theta)
        self.mPos = self.mBold.mPos + dr
        
        # Determine the Meek Atom angle
        self.DetermineMeek()

    def DetermineMeek(self):
        """
        Determines the Bond-Meek-FOD angle. I.e., the angle between the bonding axis
        and the Meek-FOD vector.
        """
        self.mMeekAngle = AngleBetween(-self.mBondDir, self.mPos - self.mMeek.mPos)

    def DetermineHeight(self):
        pass