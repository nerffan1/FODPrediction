from FODLego.Funcs import *
from FODLego.FOD import FOD
from scipy.spatial import distance

class Bond:
    def __init__(self,start,end,order):
        self.mAtoms = (start,end)
        self.mOrder = order
        self.mFODs: List[FOD] = []

    # Methods
    def SetFODs(self, fods):
        """
        Set the FODs belonging to this bond.
        """
        self.mFODs = fods

    def GetDist(self):
        """
        Returns the distance between FODs. If 2 BFODs, then it simply returns their distance.
        """
        if self.mOrder == 2:
            a = self.mFODs[0].mAssocFOD
            b = self.mFODs[1].mAssocFOD
            d = dist(a.mPos, b.mPos)
        elif self.mOrder == 3:
            a = self.mFODs[0].mAssocFOD.mPos
            b = self.mFODs[1].mAssocFOD.mPos
            c = self.mFODs[2].mAssocFOD.mPos
            # Get their pairwise distance and average
            arr = np.vstack((a,b,c))
            dd = distance.cdist(arr,arr)
            dd = np.tril(dd)
            dd[dd==0] = np.nan
            d = np.nanmean(dd)

        else:
            d = 0
        return d

    def GetMeekR(self) -> tuple[float, float]:
        meeks = []
        for fod in self.mFODs:
            meeks.append(fod.mAssocFOD.mMeekR)
        return meeks

    def GetBoldR(self) -> tuple[float, float]:
        bolds = []
        for fod in self.mFODs:
            bolds.append(fod.mAssocFOD.mBoldR)
        return bolds

    def GetBoldAng(self) -> tuple[float, float]:
        bolds = []
        for fod in self.mFODs:
            bolds.append(fod.mAssocFOD.mBoldAngle)
        return bolds

    def GetMeekAng(self) -> tuple[float, float]:
        meeks = []
        for fod in self.mFODs:
            meeks.append(fod.mAssocFOD.mMeekAngle)
        return meeks
    
    def GetPort(self):
        ports = []
        for fod in self.mFODs:
            ports.append(fod.mAssocFOD.mBoldPortion)
        return ports


    ### Magic Methods ###
    def __str__(self) -> str:
        return f"{self.mAtoms[0].mName}--{self.mAtoms[1].mName}. Order: {self.mOrder}"
