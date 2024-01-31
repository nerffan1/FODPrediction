from Funcs import *
from FOD import *
from ElementaryClasses import *
from globaldata import *
from BFOD import BFOD

class FFOD(FOD):
    def __init__(self, atom: Atom):
        super().__init__()
        self.mAtom = atom
        self.mR = -1.0
        self.mFreeDir = self.CalcFreeDirection()
    
    def CalcFreeDirection(self) -> np.ndarray:
       """
       Calculate the direction of the Free FOD. 
       """ 
       for fod in self.mAtom.mFODStruct.mBFODs:
           resultant += fod.mBondDir
       return resultant
    
class SFFOD(FFOD):
    pass

class DFFOD(FFOD):
    pass

class TFFOD(FFOD):
    pass