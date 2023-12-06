#Author: Angel-Emilio Villegas S.
import numpy as np
from typing import List

def Translation(dir: np.array):
    """
    Create a translational matrix from vector of translation
    Returns Translation matrix
    """
    # Affine transformation
    T = np.array([1,0,0,dir[0]],
                 [0,1,0,dir[1]],
                 [0,0,1,dir[2]],
                 [0,0,0,1],)

    return T

def TransformPoints(T: np.array, points: List[np.array]):
    """
    Translates all points equally in the same direction
    T: Transformation matrix
    points: A list of points  
    """
    #Add 1's to last 'row' 
    l = len(points[0])
    lastrow = np.full(l,1)
    points.append(lastrow)
    
    return (np.matmul(T,points.T)[:-1]) #Return all except last
