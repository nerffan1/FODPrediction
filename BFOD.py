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
        self.mBondDir = meekAt.mPos - boldAt.mPos #Always in the direction leaving the Bold atom
        self.mDistance = np.linalg.norm(self.mBondDir)
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
    def __init__(self, bold: Atom, meek: Atom):
        super().__init__(bold,meek)
        self.mBoldPortion = self.Calc_AxisBoldPortion(bold.mZ, meek.mZ)
        self.DetermineParamenters()

    def DetermineParamenters(self):
            """
            Add the single BFOD along the axis. 
            """
            bfod = AxialPoint_Simple(self.mBold, self.mMeek, self.mBondDir)
            self.mPos = self.mBold.mPos + bfod
            self.mBoldR = self.mDistance*self.mBoldPortion
            self.mMeekR = self.mDistance*(1-self.mBoldPortion)

class DBFOD(BFOD):
    def __init__(self, bold: Atom, meek: Atom, heightdir: np.ndarray):
        super().__init__(bold,meek)
        self.mBoldAngle = -1.0
        self.mMeekAngle = -1.0
        self.mHeight = np.array([0.0,0.0,0.0])


    def CalcHeight(bold: Atom, axisproj: float):
        """
        Return radial distance from interatomic (bonding) axis. 
        For heterogenous atoms....
        """
        # Radius Logic
        elecs = GlobalData.GetFullElecCount(bold.mGroup, bold.mPeriod)
        rad = GlobalData.mRadii[elecs][bold.mZ]
        return sqrt(rad**2 - ((axisproj)**2))

    def DetermineParamenters(self):
        pass

    def Duplicate(self):
        """
        Create the 
        """    

class TBFOD(BFOD):
    def __init__(self,bold: Atom, meek: Atom):
        super().__init__(bold,meek)