from  typing import List
import numpy as np
from ElementaryClasses import Atom

def AddNormals(vectors: list[np.array]) -> np.ndarray:
    """This function normalizes all the given vectors, adds their normal,
    and normalizes the resultant vector.
    
    This function is motivated by the observation that FFODs are closer to the
    bisecting angle in double FFODs, rather than a weighted """
    free_dir_norm = [1/np.linalg.norm(x) for x in vectors]
    free_dir = vectors * np.reshape(free_dir_norm,(len(vectors),1))
    free_dir = free_dir.sum(0)
    free_dir /= np.linalg.norm(free_dir)
    return free_dir

def DominantAtom(at1: Atom, at2: Atom, all=True):
    """
    This Function determines the dominant atom in the bonding and its distance to the weaker atom.
    If the 
    at1: An atom
    at2: An atom bonding to at2
    all: Boolean to return dominant and weak atom. Default only return dominant atom.
    """
    if at1.mPeriod < at2.mPeriod:
        fugal = at2.mPos - at1.mPos
        dom = at1
        weak = at2
    elif at1.mPeriod > at2.mPeriod:
        fugal = at1.mPos - at2.mPos
        dom = at2
        weak = at1
    elif at1.mZ > at2.mZ:
            fugal = at2.mPos - at1.mPos
            dom = at1
            weak = at2
    elif at1.mZ <= at2.mZ:
            fugal = at1.mPos - at2.mPos
            dom = at2
            weak = at1
    # Either return dom and sub, or just the dominant atom. 
    if all:
        return dom, weak, fugal
    else:
        return dom, fugal

def AxialPoint_Simple(dom:Atom, sub:Atom, dir:np.ndarray) -> np.ndarray:
            """
            Return FOD location for a Single FOD representing a Single Bond.
            If the bonding atoms are the same type, then the midpoint is chosen.
            Atom1 is assumed to be bigger than Atom2.  
            TODO: Should just return a length since we can calculate dominant one
            """
            Z1 = dom.mZ
            Z2 = sub.mZ
            if dom.mZ == sub.mZ:
                return dir*0.5
            else:
                Z1 = 0.4 if Z1 == 1 else Z1
                Z2 = 0.4 if Z2 == 1 else Z2
                r = sqrt(Z1/Z2)
                g = r/(1+r)
            if g < 0.5:
                return dir*g
            else:
                return dir*(1-g)