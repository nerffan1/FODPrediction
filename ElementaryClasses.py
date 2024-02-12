#Description: This file contains as set of classes that will implement the FOD Heuristic Solution 
#  following the paradigm of Object-Oriented Programming (OOP). Encapsulation for FOD_Structure, FODs,
#  FOD_Orbital, among other things will be included in this file.
# Roadmap: Use polymorphism for Closed Hybrid Shells (e.g. sp3).
#  -  Add Tetrahedra class and several attributes/methods to manipulate them
#  - In far future, somehow implement the triaugmented triangular prism that corresponds to sp3d5 ( 9 FODs, 18 electrons) 
#Author: Angel-Emilio Villegas S.
from ast import Global
from  globaldata import GlobalData
import Shells
import numpy as np
from Funcs import *
from numpy import sqrt 
from numpy.linalg import norm
from typing import List
from scipy.spatial.transform import Rotation as R
import csv
from Bond import *
from FOD import *

################# FOD STRUCTURE #################

class Atom:
    def __init__(self, index: int, Name: str, Pos):
        #Known Attributes
        self.mName = Name
        self.mPos = np.array(Pos) 
        self.mI = index
        self.mZ = GlobalData.GetElementAtt(self.mName, "AtomicNumber")
        self.mPeriod = int(GlobalData.GetZAtt(self.mZ, "Period" ))
        self.mGroup = int(GlobalData.GetZAtt(self.mZ, "Group" ))
        self.mValCount = self._FindValence()

        #Undetermined Attributes
        self.mSteric = 0
        self.mFreePairs = 0
        self.mCharge = 0  # In the future can be changed 
        self.mBonds = []
        self.mFODStruct = FODStructure(self)
        self.mCompleteVal = False
        
    #Parameters from GlobalData
    def GetMonoCovalRad(self): 
        elecs = GlobalData.GetFullElecCount(self.mGroup, self.mPeriod)
        return GlobalData.mRadii[elecs][self.mZ]

    def GetMonoCovalEdge(self):
        """
        Get the FOD edge distance of a monoatomic calculation
        """
        elecs = GlobalData.GetFullElecCount(self.mGroup, self.mPeriod)
        return GlobalData.mVert[elecs][self.mZ]

    def GetLastAtomRadius(self):
        dist = np.linalg.norm(self.mFODStruct.mCore[-1].mPos - self.mPos)
        return dist

    def GetValenceFODs(self):
        return self.mFODStruct.mValence

    def GetBonds(self):
        return self.mBonds

    def GetBFODs(self):
        return self.mFODStruct.mBFODs

    def AddBond(self, atom2: int, order: int):
        self.mBonds.append(Bond(self.mI, atom2,order))

    def AddBFOD(self, fod):
        self.mFODStruct.mBFODs.append(fod)
    
    def GetVec2BFODs(self):
        return [self.mPos - x.mPos for x in self.mFODStruct.mBFODs]

    def GetVectoNeighbors(self):
        """
        This function returns the sum of the vectors that start on
        neighboring (bonded) atoms and end on current atom.
        """
        freedir = []
        for bond in self.mBonds:
            at2 = GlobalData.mAtoms[bond.mAtoms[1]]
            freedir.append(at2.mPos - self.mPos)
        return freedir

    def AverageBFODDir(self):
        """
        Returns the sum of all the BondDir of all owned BFODs. 
        It is of particular use for finding the direction of FFODs.
        TODO: Might require reimplementation for cases where we just want the ATOM-ATOM vector, instead of the FOD vectors?
        """
        resultant = np.zeros(3)
        bfods = self.mFODStruct.mBFODs
        for bfod in bfods:
            # If the ffod atom is Meek, then the vector points in its direction! So, no issue here!
            if self == bfod.mMeek:
                resultant += bfod.mBondDir
            else: 
                # When the atom is Bold, the BondDir points away, so you must get the negative
                resultant -= bfod.mBondDir
        # Return the average
        return resultant/len(self.mFODStruct.mBFODs)

    def CalcSteric_test(self) -> None:
        """
        TODO: Need to account for systems where the valence electrons + bonding FODs
        """
        self.mFreePairs = int((self.mSteric - np.sum([2*x.mOrder for x in self.mBonds ]))/2)
        self.mSteric = self.mFreePairs + len(self.mBonds)
    
    def CalcSteric(self) -> None:
        """
        This function assumes that the system at hand only contains Closed Shell calculations.
        TODO: Need to account for systems where the valence electrons + bonding FODs
        """
        #Electrons involved in the Bond
        bondelec = np.sum([2*bond.mOrder for bond in self.mBonds])
        # The difference between the total electrons and the number of electrons that fill the shell
        self.mCharge = (self.mZ + bondelec) - GlobalData.GetFullElecCount(self.mGroup, self.mPeriod)
        self.mFreePairs = int(GlobalData.mShellCount[self.mPeriod] - bondelec)/2
        self.mSteric = self.mFreePairs + len(self.mBonds)
    
    def _FindValence(self):
        """
        This method finds the number of electrons in the valence shell by finding the 
        difference between the current Group and the last ClosedShell Group. Only works up
        to 5th period.
        TODO: This can be saved in GlobalData
        """
        if self.mGroup < 4:
            return self.mGroup
        else:
            if self.mPeriod < 4:
                return (2 + (self.mGroup - 12))
            elif self.mPeriod < 6:
                return (self.mGroup)
                    
    def _CheckFullShell(self):
        """
        Check that the atom has a closed shell.
        Future: Add a variable that saves the info so that looping every time
        this is called (if called more than once) is unnecesary
        TODO: Alternatively, just check the amount of electrons
        """
        #Determine How many electrons are needed to fill shell
        for ClGrp in GlobalData.mClosedGroups:
            if self.mGroup < ClGrp:
                checkshell = ClGrp - self.mGroup + self.mCharge
                break 
        for bond in self.mBonds:
            checkshell -= bond.mOrder
        if checkshell == 0:
            return True
        else:
            return False
    
    # Additional Functions
    def __str__(self):
        pass    

######################## FOD Structure Class  ########################
#HOW TO concatenate FODs, easily
class FODStructure:
    def __init__(self, parent: Atom):
        self.mAtom = parent
        self.mCore = [] #A list of FODShells
        self.mCoreShells = []
        self.mValence = [] #A list of FODs
        self.mBFODs = []
        self.mFFODs: List[FOD] = []
        self.mfods = [] #All Finalized FODs
        self.mLastBond = 0
    
    # Setter Functions
    def AddValFOD(self, fods: List[np.ndarray], at2=False,finalize=False):
        """
        This function adds the FODs to the valence of the current atom's FOD Structure.
        It also accepts a secondary atom, at2, in order to add the FOD to its valence. The
        Finalize parameter is True only for the first atom that is writing the FODs; by finalizing
        we mean that the FOD is added to the list of overall FODs.
        TODO: Currently only takes positions. Do we want to also add it here?
        """
        # Add to current valence
        for fod in fods:
            #Assertions
            assert isinstance(fod, np.ndarray), "The variable is not a NumPy ndarray."

            self.mValence.append(fod)
            
            #Add to bonded atom valence, if passed
            if at2 != False:
                at2.mFODStruct.AddValFOD([fod])

            #Add to the final list of atoms, without duplicating
            if finalize:
                if len(self.mfods) == 0:
                    self.mfods = fod  
                else:
                    self.mfods = np.vstack((self.mfods,fod))

    def _AddCoreShell(self, shell):
        self.mCoreShells.append(shell)
        # Add individual FODs to the electronic structure
        for fod in shell.mfods:
            GlobalData.mFODs.append(fod)
            self.mCore.append(fod)

    def PrepareShells(self, atoms: List[Atom]):
        """
        This function will determine the creation of Core FODs, those that are not 
        related to bonding. 
        Comment: The scheme is easy in the first 3 periods of the Periodic
        Table, but it will become trickier ahead if Hybridization heuristics don't work.
        Currently it only works for closed shell calculations (V 0.1).
        TODO: This will assume that we are doing up to the n=3 shell, with sp3 hybridization
        TODO: Take into account free pairs of electrons
        TODO: For mOrder=2, there are many schemes to determine direction
        TODO: Find an elegant solution to do exceptions for H bonds
        TODO: Account previous bonds formed, by talling previous FODs, or looking back at mBonds 
        TODO: GlobalData.GetFullElecCount() Can be precalculated ahead of time and placed as a member vatrable
        """
        #Lazy loading in order to 
        from FOD import FOD
        from BFOD import SBFOD, DBFOD, TBFOD  
       
        at1 = self.mAtom

        def SingleBond(at2: Atom):
           """
           Set parameters from Single Bond.
           -Accounted for in SBFOD 
           """
           dom,weak,dir = BoldMeekDir(at1,at2,True)
           bfod = AxialPoint_Simple(dom, weak, dir)
           self.AddValFOD([dom.mPos + bfod], at2, True)
           #Add the BFODs, new way
           boldMeek = BoldMeek(at1,at2)
           newFOD = SBFOD(*boldMeek)
           _AddFOD(at1,at2,newFOD)

        def _AddFOD(at1: Atom, at2: Atom, fod: FOD):
            """
            This function adds a new FOD to the individual atoms, to the list in GlobalData, and to the FODStructure
            TODO: Add the FOD to the Valence Structure?
            """ 
            # The main purpose of this function is not to not duplicate the FOD in GlobalData
            # by adding the FODs in each individual atom. This makes this class a type of 
            # FOD manager in addition to constructing the structure.
            # The reason why we don't load FODs directly
            # is because we don't know whether at1 or at2
            # is the self.mAtom
            at1.AddBFOD(fod)
            at2.AddBFOD(fod)
            GlobalData.mFODs.append(fod)          
        
        def _AddFFOD(ffod: FOD):
            """
            TODO: Put this on a bigger scope
            """
            self.mFFODs.append(ffod)
            GlobalData.mFODs.append(ffod)          

        def DoubleBond(at2: Atom, bond: Bond):
            """
            Create the FODs representing the double bond. Currently the FOD filling is unidirectional (sequential)
            and  does not account for the next atom in the iteration to see if we can further accomodate the bonding FODs 
            """
            def D_Bond_Height(dom: Atom, axisproj: float):
                """
                Return radial distance from interatomic (bonding) axis. 
                For heterogenous atoms...
                NOTE: This functionality has been moved to class DBFOD
                """
                # Radius Logic
                elecs = GlobalData.GetFullElecCount(dom.mGroup, dom.mPeriod)
                rad = GlobalData.mRadii[elecs][dom.mZ]
                return sqrt(rad**2 - ((axisproj)**2))
                            
            def HeightDir_fromNeighborBFODs():
                """
                This returns the unit-vector of the direction that a DBFOD is displaced away from 
                the bonding axis (called 'Height' throughout this code). It is obtained by a series of steps:
                1) The cross-product of the FOD-Atom-FOD is obtained
                2) Measure the angle between the Bonding Axis (BA) and the vector found in (1) 
                3) If beyond a certain threshold, then rotate the direction 
                
                """
                #(1)Obtain the FOD-Atom-FOD Cross Product. The height in a planar structure
                vector_for_cross = []
                for otherb in self.mAtom.mBonds:
                    if otherb != bond:
                        vector_for_cross.append(atoms[otherb.mAtoms[1]].mPos)
                vector_for_cross -= self.mAtom.mPos
                D = np.cross(*vector_for_cross)
                D = normalize(D) 
                #(2)Measure Angle between D and the Bonding Axis (BA) 
                BA = at2.mPos - at1.mPos
                angle = AngleBetween(BA,D)
                
                #(3) Check for orthogonality
                thres_max = 1.01*(np.pi/2) # +1% Deviation from 99%
                thres_min = .99*(np.pi/2) # -1% Deviation from 99%
                if angle > thres_max or angle < thres_min:
                    #Get axis of rotation
                    axis = np.cross(BA,D)
                    axis = normalize(axis)
                    #Get rotation angle
                    if angle > thres_max: 
                        angle_diff =  np.pi/2 - angle
                    elif angle < thres_min:
                        angle_diff = (np.pi/2) - angle
                    #Rotate
                    D = RotateVec(D, axis*angle_diff)
 
                return D

            def D_FFOD_Direction():
                if self.mAtom.mFreePairs == 0:
                    #Determine DFFOD Direction
                    if len(self.mAtom.mBonds) == 3:
                        return  HeightDir_fromNeighborBFODs()
                    elif len(self.mAtom.mBonds) == 2:
                        if self.mAtom.mBonds.index(bond) == 0:
                            return  RandomPerpDir(dir)
                        else:
                            # Cross product between atom and already-placed FODs
                            vector_for_cross = []
                            for fods in self.mValence:
                                vector_for_cross.append(fods)
                            vector_for_cross -= self.mAtom.mPos
                            return np.cross(*vector_for_cross)   
                else:
                    return RandomPerpDir(dir)

            #Information
            axis2fod = np.ndarray(3)
            dom,sub,fugal = BoldMeekDir(at1,at2, True)

            if GlobalData.GetFullElecCount(self.mAtom.mGroup,self.mAtom.mPeriod) <= 18:
                #Find perpendicular unit vector
                if self.mAtom.mFreePairs == 0:
                    axis2fod = D_FFOD_Direction()
                    _AddFOD(dom,sub, DBFOD(dom,sub,axis2fod))
                    _AddFOD(dom,sub, DBFOD(dom,sub,-axis2fod))

                #Add both FODs of the Double Bond
                midpoint = AxialPoint_Simple(dom,sub,fugal)
                # Determine Vertical Projection of DFFOD
                axis2fod *= D_Bond_Height(dom, np.linalg.norm(midpoint))

                #Add FODs
                self.AddValFOD([dom.mPos + midpoint + axis2fod],at2,True)
                self.AddValFOD([dom.mPos + midpoint - axis2fod],at2,True)

        def AddFreeElectron(free: int):
            freedir = []
            """
            TODO: Change conditionals as to create more concise code 
            TODO: Create a series of variables for chosen constants, NO magic numbers
            """
            #Variables
            vector_for_cross = []
            fulle = GlobalData.GetFullElecCount(self.mAtom.mGroup, self.mAtom.mPeriod)
            vector_for_cross = self.mAtom.mPos - [fod for fod in self.mValence]
                
            if len(self.mAtom.mBonds) == 1:
                #Useful Information
                at2 = atoms[self.mAtom.mBonds[0].mAtoms[1]]
                free_dir = self.mAtom.mPos - at2.mPos
                free_dir = normalize(free_dir)

            if free == 1:
                from FFOD import SFFOD
                # For a single atom
                if len(self.mAtom.mBonds) == 3 :
                    dr = vector_for_cross.sum(0)
                elif len(self.mAtom.mBonds) == 2 :
                    free_dir = AddNormals(vector_for_cross)
                    dr = free_dir*GlobalData.mRadii[fulle][self.mAtom.mZ]
                elif len(self.mAtom.mBonds) == 1:
                    l = GlobalData.mVert[fulle][self.mAtom.mZ]
                    if (self.mAtom.mPeriod < 3):
                        dr = free_dir*l/sqrt(8) #Place at midsphere distance
                    else:
                        dr = free_dir*GlobalData.mRadii[fulle][self.mAtom.mZ]
                self.AddValFOD([self.mAtom.mPos+dr],False,True)

            elif free == 2:
                # #direction away from atom, on plane of other bonds, or neighboring atoms
                # axis2fod = np.cross(*self.mAtom.GetVectoNeighbors())
                # axis2fod = normalize(axis2fod)
                # #get angle to atoms of ffods (ffod-atom-ffod)
                # theta = np.deg2rad(220) - anglebetween(*vector_for_cross)
                # theta /= 2 

                # if len(self.mvalence) == 2: #might be different for fods or for number of bonds
                #     #determine the free direction
                #     dxy  = addnormals(vector_for_cross)
                #     if len(at1.mbonds) == 2:
                #         dxy = at1.mpos - [atoms[bonds.matoms[1]].mpos for bonds in at1.mbonds]
                #         dxy = dxy.sum(0)
                #         dxy = normalize(dxy)

                #     #todo: must determine logic to choose distance, based of bonding as well 
                #     rad = np.linalg.norm(vector_for_cross[0]) # remember these are the atom-bfod distances
                #     dxy *= rad*np.cos(theta)
                #     axis2fod *= rad*np.sin(theta)
                # elif len(self.mvalence) == 1:
                #     #add both fods of the double bond
                #     axis2fod *= globaldata.mvert[fulle][self.matom.mz]/2
                #     dxy = self.matom.mpos + free_dir*np.linalg.norm(axis2fod)/np.tan(np.deg2rad(54.735))
                
                # self.AddValFOD([self.mAtom.mPos + dxy + axis2fod],False,True)
                # self.AddValFOD([self.mAtom.mPos + dxy - axis2fod],False,True)
                heightdir =  np.cross(*[fod.mPos for fod in self.mBFODs],axis=0)
                heightdir = normalize(heightdir)
                from FFOD import DFFOD
                _AddFFOD(DFFOD(at1,heightdir))
                _AddFFOD(DFFOD(at1,-heightdir))

            elif free == 3:
                if self.mAtom.mZ < 10:
                    #Create the starting FOD. Begin with the vertical component
                    R_f = sqrt(3/8)*np.linalg.norm(at1.mPos - at1.mFODStruct.mValence[0])
                elif at1.mPeriod > 2:                                                            
                    R_f = np.linalg.norm(at1.mPos - at1.mFODStruct.mValence[0])
                # Begin First FOD 
                axis2fod = RandomPerpDir(free_dir)
                axis2fod *= np.sin(np.deg2rad(70.5288))*R_f 
                #Create the horizontal component, parallel to the bonding axis
                horizontal = np.cos(np.deg2rad(70.5288))*R_f*free_dir
                #Add vectors, and rotate to create equilateral triangle
                dr = horizontal + axis2fod
                #Create rotations and rotated FODs
                rot_fods = RotatePoints(3, dr, free_dir)

                #Translate to the atom of interest, and add to Valence 
                self.AddValFOD(at1.mPos + rot_fods,False,True)
                # --------------------------------------------
                # TODO: Remove stuff above
                from FFOD import TFFOD
                bonddir = tofrom(at2.mPos,at1.mPos)
                dir0 = RandomPerpDir(bonddir)
                norms = RotateNormals(3, dir0, normalize(bonddir)) 
                _AddFFOD(TFFOD(at1, norms[0]))
                _AddFFOD(TFFOD(at1, norms[1]))
                _AddFFOD(TFFOD(at1, norms[2]))

        def PlaceFODs_Triple(fugal: np.ndarray, a: float, rad: float, at2: Atom, c: float):
            """ This function create an FOD at a certain distance, based of a and rad,
            which are quantities assumed to dominate the interaction.
            NOTE: Has been implemented in TBFOD Class. Deprecated
            """
            equilat_r =a/sqrt(3)
            if at1.mZ == at2.mZ:
                if at1.mPeriod ==2 & at2.mPeriod == 2:
                    theta = np.arccos((c/2)/rad)
                elif at1.mPeriod > 2 & at2.mPeriod > 2:
                    theta = np.arctan(rad/(c/2)) 
            else:
                theta = np.arcsin((a/2)/rad) 
            dr = fugal*rad*np.cos(theta) + RandomPerpDir(fugal)*rad*np.sin(theta)
            return RotatePoints(3, dr, fugal)

        def TripleBond(at2: Atom):
            """
            #TODO: Create a helper funtion for conditional statements
            """
            #Determine dominant atom and direction
            atom, fugal = BoldMeekDir(at1,at2,False)
            # Variables
            c = np.linalg.norm(fugal) # Distance
            fugal = normalize(fugal)
            elecs = GlobalData.GetFullElecCount(atom.mGroup, atom.mPeriod)
            rad = GlobalData.mRadii[elecs][atom.mZ] #Radius
            a = GlobalData.mVert[elecs][atom.mZ] #Edge
            #Place FODs
            rot_fods = PlaceFODs_Triple(fugal, a, rad, at2, c)
            self.AddValFOD(atom.mPos + rot_fods,at2,True)
            # Place BFODs, for new implementation
            bonddir = tofrom(at2.mPos,at1.mPos)
            boldmeek = BoldMeek(at1,at2)
            dir0 = RandomPerpDir(bonddir)
            norms = RotateNormals(3, dir0, normalize(bonddir)) 
            fod1 = TBFOD(*boldmeek, norms[0])
            fod2 = TBFOD(*boldmeek, norms[1])
            fod3 = TBFOD(*boldmeek, norms[2])
            _AddFOD(at1,at2,fod1)
            _AddFOD(at1,at2,fod2)
            _AddFOD(at1,at2,fod3)

        def AddBFODs():
            """
            Finish determining BFODs. This is done after initializing all initial BFODs.
            """
            for bond in self.mAtom.mBonds:
                bonded_at = atoms[bond.mAtoms[1]] 
                if bonded_at.mI  > self.mAtom.mI:
                    if bond.mOrder == 1:
                        SingleBond(bonded_at)
                    elif bond.mOrder == 2:
                        DoubleBond(bonded_at, bond)
                    elif bond.mOrder == 3:
                        TripleBond(bonded_at)

        def AddFFODs():
            from FFOD import SFFOD, DFFOD, TFFOD 
            if self.mAtom.mFreePairs == 2:
                if self.mAtom.mSteric >= 3:
                    AddFreeElectron(2)
            elif self.mAtom.mFreePairs == 1:
                if self.mAtom.mSteric >= 2:
                  #  AddFreeElectron(1)
                    _AddFFOD(SFFOD(self.mAtom))
            elif self.mAtom.mFreePairs == 3:
                AddFreeElectron(3)

        def AddCoreElectrons():
            #Count core electrons and
            core_elec = self.mAtom.mZ - self.mAtom.mValCount
            if core_elec != 0:
                for shell in GlobalData.mGeo_Ladder[core_elec]:
                    if shell == 'point':
                        self._AddCoreShell(Shells.Point(self.mAtom))
                    elif shell == 'tetra':
                        self._AddCoreShell(Shells.Tetra(self.mAtom, 10))
        
        #Prepare the valence shell first, since it will help determine the
        # orientation of the inner shells
        AddBFODs()
        AddCoreElectrons()
        AddFFODs()
            
    def FinalizeFODs(self):
        """
        Add all FODs in the FODStructure, to the atom, so that they can be
        printed when creating the XYZ output file
        """
        #Add the core FODs
        for shell in self.mCore:
            shell.mfods += self.mAtom.mPos   
            if len(self.mfods) == 0:
                self.mfods = [shell.mfods]
            else:
                self.mfods = np.vstack((self.mfods,shell.mfods))  ###HOW TO concatenate FODs, easily
    
################# ADDITIONAL FUNCTIONS #################
def BoldMeekDir(at1: Atom, at2: Atom, all=True):
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
        sub = at2
    elif at1.mPeriod > at2.mPeriod:
        fugal = at1.mPos - at2.mPos
        dom = at2
        sub = at1
    elif at1.mZ > at2.mZ:
            fugal = at2.mPos - at1.mPos
            dom = at1
            sub = at2
    elif at1.mZ <= at2.mZ:
            fugal = at1.mPos - at2.mPos
            dom = at2
            sub = at1
    # Either return dom and sub, or just the dominant atom. 
    if all:
        return dom, sub, fugal
    else:
        return dom, fugal

def BoldMeek(at1: Atom, at2: Atom):
    """
    This Function determines the dominant and meek atom in the bonding.
    at1: An atom
    at2: An atom bonding to at2
    """
    if at1.mPeriod < at2.mPeriod:
        dom = at1
        sub = at2
    elif at1.mPeriod > at2.mPeriod:
        dom = at2
        sub = at1
    elif at1.mZ > at2.mZ:
            dom = at1
            sub = at2
    elif at1.mZ <= at2.mZ:
            dom = at2
            sub = at1
    # Either return dom and sub, or just the dominant atom. 
    return dom, sub 

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
