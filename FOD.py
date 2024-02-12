import numpy as np

class FOD:
    def __init__(self, pos: np.ndarray = np.zeros(3)) -> None:
        self.mPos = pos

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