#Author: Angel-Emilio Villegas S.
import numpy as np

def CreateTranslation(dir: np.array):
    """
    Create a translational matrix
    Returns Translation matrix
    """
    # Affine transformation
    T = np.array([1,0,0,dir[0]],
                 [0,1,0,dir[1]],
                 [0,0,1,dir[2]],
                 [0,0,0,1],)

    return T

def TranslatePoints(dir: np.array, points: [np.array]):
    """
    Translates all points equally in the same direction
    """
    #Add a last one to create proper transformation
    l = len(points[0])
    lastrow = np.full(l,1)
    points.append(lastrow)
    T = CreateTranslation(dir)
    output = np.matmul(T,points.T)
    return output[:-1]
