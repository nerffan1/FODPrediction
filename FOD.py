import numpy as np
from typing import List

class FOD:
    def __init__(self, pos: np.ndarray = np.zeros(3)) -> None:
        self.mPos = pos
        self.mSiblings = ()
        self.mAssocFOD = None

    def AddSibling(self, *siblings):
        """
        Add associated FODs, siblings, to the FOD. For example, the 3 FODs in a TBFOD would be the siblings, where each points to the other 2.
        """
        self.mSiblings = siblings

    #################### Operator Overloading #################### 
    def __str__(self) -> str:
        return str(self.mPos)

    def __mul__(self, factor: float):
        self.mPos *= factor
        return self

    def __add__(self, shift: float):
        """
        Overloading addition operation to add shift the mPos value. Mutates instance.
        """
        self.mPos += shift
        return self

    def __sub__(self, shift: float):
        """
        Overloading addition operation to add a negative shift the mPos value. Mutates instance.
        """
        self.mPos -= shift
        return self