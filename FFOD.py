from Funcs import *
from FOD import *
from ElementaryClasses import *
from globaldata import *
from BFOD import BFOD

class FFOD(FOD):
    def __init__(self, atom: Atom):
        super().__init__()
        self.mPos = atom.mPos
        self.mAtom = atom
        self.mR = -1.0
        self.mFreeDir = self.CalcFreeDirection()
    
    def CalcFreeDirection(self) -> np.ndarray:
        """
        Calculate the direction of the Free FOD.
        TODO: Needs to be cleaned
        """
        bfods = self.mAtom.mFODStruct.mBFODs
        if bfods[0].mBold == self.mAtom:
            resultant = self.mAtom.mPos - np.sum([x.mBondDir for x in bfods], axis=0)
        else:
            resultant = self.mAtom.mPos + np.sum([x.mBondDir for x in bfods], axis=0)
        return resultant
    
class SFFOD(FFOD):
    def __init__(self, atom: Atom):
        super().__init__(atom)
        self.DetermineParameters()

    def DetermineParameters(self) -> None:
        """
        This function determines the distance away from parent atom of the SFFOD. 
        """
        if len(self.mAtom.mBonds) == 3 :
            # The vectors have been added already 
            dr = self.mFreeDir/3 
        elif len(self.mAtom.mBonds) == 2 :
            #Perhaps add a restriction here later
            dr = normalize(self.mFreeDir)*self.mAtom.GetMonoCovalRad()
        elif len(self.mAtom.mBonds) == 1:
            # The following section has 2 options. In atoms with just a 1s core shell 
            # tend to be closer to the nucleus when their valence is complete.
            if (self.mAtom.mPeriod < 3):
                dr = normalize(self.mFreeDir)*self.mAtom.GetMonoCovalEdge()/sqrt(8)  #Place at midsphere distance?
            else:
                dr = self.mFreeDir/3
        self.mPos = self.mAtom.mPos + dr

        def SB_FreeDirection():
            bfods = self.mAtom.mFODStruct.mBFODs
            if bfods[0].mBold == self.mAtom:
                resultant = self.mAtom.mPos - np.sum([x.mBondDir for x in bfods], axis=0)
            else:
                resultant = self.mAtom.mPos + np.sum([x.mBondDir for x in bfods], axis=0)
            return resultant
            

class DFFOD(FFOD):
    def __init__(self, atom: Atom):
        super.__init__(atom)

    def DetermineParameters():
        pass

class TFFOD(FFOD):
    def __init__(self, atom: Atom):
        super.__init__(atom)

    def DetermineParameters():
        pass