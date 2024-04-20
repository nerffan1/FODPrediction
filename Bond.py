from globaldata import GlobalData

class Bond:
    def __init__(self,start,end,order):
        self.mAtoms = (start,end)
        self.mOrder = order
        self.mAtoms_p = (GlobalData.mAtoms[start], GlobalData.mAtoms[end])

    
    def __str__(self) -> str:
        return f"From {self.mAtoms[0]} to {self.mAtoms[1]}. Order: {self.mOrder}"
 