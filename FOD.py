import numpy as np

class FOD:
    def __init__(self) -> None:
        self.mPos = np.array([0.0,0.0,0.0])

    def __str__(self) -> str:
        return str(self.mPos)