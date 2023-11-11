#Description: This file contains as set of classes that will implement the FOD Heuristic Solution 
#  following the paradigm of Object-Oriented Programming (OOP). Encapsulation for FOD_Structure, FODs,
#  FOD_Orbital, among other things will be included in this file.
# Roadmap: Use polymorphism for Closed Hybrid Shells (e.g. sp3).
#  -  Add Tetrahedra class and several attributes/methods to manipulate them
#  - In far future, somehow implement the triaugmented triangular prism that corresponds to sp3d5 ( 9 FODs, 18 electrons) 
#Author: Angel-Emilio Villegas S.
from  globaldata import GlobalData
import math3d
import numpy as np
from numpy import sqrt 
from typing import List
from numpy.linalg import inv, det
from scipy.spatial.transform import rotation as rot
import csv

#FOD STRUCTURE

class Bond:
    def __init__(self,start,end,order):
        self.mAtoms = (start,end)
        self.mOrder = order
    def __str__(self) -> str:
        return f"From {self.mAtoms[0]} to {self.mAtoms[1]}. Order: {self.mOrder}"
    
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
        self.mCharge = 0 #Temporarily zero for all atoms. Ionic atoms can be dealt with in the future
        self.mBonds = []
        self.mFODStruct = FODStructure(self)
        self.mCompleteVal = False
        
    def _AddBond(self, atom2: int, order: int):
        self.mBonds.append(Bond(self.mI, atom2,order))

    def CalcSteric(self) -> None:
        self.mSteric = self.mValCount
        for bond in self.mBonds:
            self.mSteric += bond.mOrder 
        self.mFreePairs = int((self.mSteric - np.sum([2*x.mOrder for x in self.mBonds ]))/2)
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
        """
        #Determine How many electrons are needed to fill shell
        for ClGrp in GlobalData.mClosedGroups:
            if self.mGroup < ClGrp:
                checkshell = ClGrp - self.mGroup
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
class FODStructure:
    def __init__(self, parent: Atom):
        self.mAtom = parent
        self.mCore = [] #A list of FODShells
        self.mValence = [] #A list of FODs
        self.mfods = [] #All Finalized FODs
        self.mLastBond = 0
    
    # Setter Functions
    def AddValFOD(self, fod: np.ndarray, at2=False,finalize=False):
        """
        This function adds the FODs to the valence of the current atom's FOD Structure.
        It also accepts a secondary atom, at2, in order to add the FOD to its valence. The
        Finalize parameter is True only for the first atom that is writing the FODs; by finalizing
        we mean that the FOD is added to the list of overall FODs.
        """
        # Add to current valence
        self.mValence.append(fod)
        
        #Add to bonded atom valence, if passed
        if at2 != False:
             at2.mFODStruct.AddValFOD(fod)

        #Add to the final list of atoms, without duplicating
        if finalize:
            if self.mfods == []:
                self.mfods = fod  
            else:
                self.mfods = np.vstack((self.mfods,fod))

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
        at1 = self.mAtom
        def AxialPoint_Simple(At1: Atom, At2: Atom):
            """
            Return FOD location for a Single FOD representing a Single Bond.
            If the bonding atoms are the same type, then the midpoint is chosen.
            Atom1 is assumed to be bigger than Atom2.  
            """
            Z1 = At1.mZ
            Z2 = At2.mZ
            #TODO: To Slater, or not to Slater
            #if Z2 > 14 and Z2 < 18 : Z2 -= (2 + 2*(.35) + 8*(.85))
            #if Z1 > 14 and Z2 < 18 : Z1 -= 14
            
            if Z2 == Z1:
                #Midpoint across atoms
                g = .5
            else:
                if At1.mPeriod == At2.mPeriod:
                    #Find a weighted point in space.
                    if Z1>Z2 and Z2 != 1:
                        r = sqrt(Z1/Z2)
                        g = 1/(1+r)
                    #Z=1 Is an exception to the last case
                    elif Z2 == 1:
                        r = sqrt(Z1/.4) # A temporary factor to
                        g = r/(1+r)
                    else:
                        r = sqrt(Z2/Z1)
                        g = r/(1+r)
                elif At1.mPeriod != At2.mPeriod:
                    #Find a weighted point in space.
                    if Z1>Z2 and Z2 != 1:
                        r = sqrt(Z1/Z2)
                        g = r/(1+r)
                    #Z=1 Is an exception to the last case
                    elif Z2 == 1:
                        r = sqrt(Z1/.4) # A temporary factor to
                        g = r/(1+r)
                    else:
                        r = sqrt(Z2/Z1)
                        g = 1/(1+r)

            #Return final value with offset
            dx = (At2.mPos - At1.mPos)*g
            return (At1.mPos + dx)
        
        def SingleBond(at2: Atom):
            bfod = AxialPoint_Simple(at1, at2)
            self.AddValFOD(bfod, at2, True)

        def RandomPerpDir(dir):
            #TODO: Remove the other RandPerp
            if dir[0] == 0: return np.array([1,0,0])
            elif dir[1] == 0: return np.array([0,1,0])
            elif dir[2] == 0: return np.array([0,0,1])
            else:
                b_z = -(10*dir[0] + 2*dir[1])/dir[2]
                randperp = np.array([10,2,b_z])
                randperp /= np.linalg.norm(randperp)
                return randperp        
        
        def DoubleBond(at2: Atom):
            """ 
            Create the FODs representing the double bond. Currently the FOD filling is unidirectional (sequential)
            and  does not account for the next atom in the iteration to see if we can further accomodate the bonding FODs 
            """
            def DoubleBond_Diatomic(at2: Atom):
                """
                Return radial distance from interatomic (bonding) axis. 
                """
                # Radius Logic
                if at1.mZ == at2.mZ:
                    c = np.linalg.norm(at2.mPos - at1.mPos)
                    elecs = GlobalData.GetFullElecCount(at1.mGroup, at1.mPeriod)
                    rad = GlobalData.mRadii[elecs][at1.mZ]
                    return sqrt(rad**2 - ((c/2)**2))
                else:
                    l = np.linalg.norm(AxialPoint_Simple(at1, at2))
                    atoom = at1 if at1.mZ < at2.mZ else at2
                    elecs = GlobalData.GetFullElecCount(atoom.mGroup, atoom.mPeriod)
                    rad = GlobalData.mRadii[elecs][atoom.mZ]
                    #Do an inscrbed triangle within the Atom-BondFOD-Atom triangle
                    return sqrt(rad**2 - ((l)**2))
                            
            #Information
            axis2fod = np.ndarray(3)
            dir = at1.mPos - at2.mPos

            def CrossPerpDir():
                vector_for_cross = []
                for otherb in self.mAtom.mBonds:
                    if otherb != bond:
                        vector_for_cross.append(atoms[otherb.mAtoms[1]].mPos)
                vector_for_cross -= self.mAtom.mPos
                vec = np.cross(*vector_for_cross)
                return vec/np.linalg.norm(vec)
                
            if GlobalData.GetFullElecCount(self.mAtom.mGroup,self.mAtom.mPeriod) <= 18:
                if self.mAtom.mFreePairs == 0:
                    if len(self.mAtom.mBonds) == 3: # Planar
                        axis2fod = CrossPerpDir()
                    elif len(self.mAtom.mBonds) == 2:
                        if self.mAtom.mBonds.index(bond) == 0:
                            axis2fod = RandomPerpDir(dir)
                        else:
                            # Cross product between atom and already-placed FODs
                            vector_for_cross = []
                            for fods in self.mValence:
                                vector_for_cross.append(fods)
                            vector_for_cross -= self.mAtom.mPos
                            axis2fod = np.cross(*vector_for_cross)   
                else:
                    axis2fod = RandomPerpDir(dir)
                    
                #Find perpendicular unit vector and multiply by chosen radius
                axis2fod *= DoubleBond_Diatomic(at2)
                
                #Add both FODs of the Double Bond
                midpoint = AxialPoint_Simple(at1, at2)

                #Add FODs
                self.AddValFOD(midpoint + axis2fod,at2,True)
                self.AddValFOD(midpoint - axis2fod,at2,True)
    
        def AddFreeElectron(free: int):
            """
            TODO: Change conditionals as to create more concise code 
            TODO: Create a series of variables for chosen constants, NO magic numbers
            """
            #Variables
            vector_for_cross = []
            fulle = GlobalData.GetFullElecCount(self.mAtom.mGroup, self.mAtom.mPeriod)
            vector_for_cross = [fod for fod in self.mValence] - self.mAtom.mPos
                
            if len(self.mAtom.mBonds) == 1:
                #Useful Information
                at2 = atoms[self.mAtom.mBonds[0].mAtoms[1]]
                free_dir = self.mAtom.mPos - at2.mPos
                free_dir /= np.linalg.norm(free_dir)

            if free == 1:
                if len(self.mAtom.mBonds) > 1 :        
                    free_dir = np.sum(vector_for_cross, axis=0)
                    free_dir /= np.linalg.norm(free_dir)
                    dr = - free_dir*GlobalData.mRadii[fulle][self.mAtom.mZ]
                elif len(self.mAtom.mBonds) == 1:
                    l = GlobalData.mVert[fulle][self.mAtom.mZ]
                    dr = free_dir*l/sqrt(8) #Place at midsphere distance    
                self.AddValFOD(self.mAtom.mPos+dr,False,True)

            elif free == 2:
                #Direction away from atom, on plane of other bonds, or neighboring atoms
                axis2fod = np.cross(*vector_for_cross)

                if len(self.mAtom.mBonds) > 1 :
                    axis2fod /= np.linalg.norm(axis2fod)
                    dr =  - np.sum(vector_for_cross, axis=0)*.2
                elif len(self.mAtom.mBonds) == 1:
                    #Add both FODs of the Double Bond
                    axis2fod *= GlobalData.mVert[fulle][self.mAtom.mZ]/2
                    dr = self.mAtom.mPos + free_dir*np.linalg.norm(axis2fod)/np.tan(np.deg2rad(54.735))
                
                self.AddValFOD(self.mAtom.mPos + dr + axis2fod,False,True)
                self.AddValFOD(self.mAtom.mPos + dr - axis2fod,False,True)

            elif free == 3:
                if self.mAtom.mZ < 10:
                    #Create the starting FOD. Begin with the vertical component
                    R_f = sqrt(3/8)*np.linalg.norm(at1.mPos - at1.mFODStruct.mValence[0])
                elif self.mAtom.mZ <=18:                                                            
                    R_f = np.linalg.norm(at1.mPos - at1.mFODStruct.mValence[0])
                # Begin First FOD 
                axis2fod = RandomPerpDir(free_dir)
                axis2fod *= np.sin(np.deg2rad(70.5288))*R_f 
                #Create the horizontal component, parallel to the bonding axis
                horizontal = np.cos(np.deg2rad(70.5288))*R_f*free_dir
                #Add vectors, and rotate to create equilateral triangle
                dr = horizontal + axis2fod
                #Create rotations and rotated FODs
                rot1 = rot.Rotation.from_rotvec((2*np.pi/3)*free_dir)
                fod1 = np.matmul(rot1.as_matrix(),dr)
                rot2 = rot.Rotation.from_rotvec(-(2*np.pi/3)*free_dir)
                fod2 = np.matmul(rot2.as_matrix(),dr)

                #Translate to the atom of interest, and add to Valence 
                self.mValence.append(at1.mPos + dr)
                self.mValence.append(at1.mPos + fod1)
                self.mValence.append(at1.mPos + fod2)

        def PlaceFODs_Triple(fugal: np.ndarray, a: float, rad: float, at2: Atom, c: float):
            """ This function create an FOD at a certain distance, based of a and rad,
            which are quantities assumed to dominate the interaction.
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
            #Create rotations and rotated FODs
            rot1 = rot.Rotation.from_rotvec((2*np.pi/3)*fugal)
            fod1 = np.matmul(rot1.as_matrix(),dr)
            rot2 = rot.Rotation.from_rotvec(-(2*np.pi/3)*fugal)
            fod2 = np.matmul(rot2.as_matrix(),dr)
            return (dr,fod1,fod2)

        def TripleBond(at2: Atom):
            """
            #TODO: Create a helper funtion for conditional statements
            """
            #Data
            if at1.mPeriod < at2.mPeriod:
                fugal = at2.mPos - at1.mPos
                atom = at1
            elif at1.mPeriod > at2.mPeriod:
                fugal = at1.mPos - at2.mPos
                atom = at2
            elif at1.mZ > at2.mZ:
                    fugal = at2.mPos - at1.mPos
                    atom = at1
            elif at1.mZ <= at2.mZ:
                    fugal = at1.mPos - at2.mPos
                    atom = at2
            # Variables
            c = np.linalg.norm(fugal) # Distance
            fugal /= np.linalg.norm(fugal)  
            elecs = GlobalData.GetFullElecCount(atom.mGroup, atom.mPeriod)
            rad = GlobalData.mRadii[elecs][atom.mZ] #Radius
            a = GlobalData.mVert[elecs][atom.mZ] #Edge
            #Place FODs
            tfodPos = PlaceFODs_Triple(fugal, a, rad, at2, c)
            self.AddValFOD(atom.mPos + tfodPos[0],at2,True)
            self.AddValFOD(atom.mPos + tfodPos[1],at2,True)
            self.AddValFOD(atom.mPos + tfodPos[2],at2,True)

        def AddCoreElectrons():
            #Count core electrons and
            core_elec = self.mAtom.mZ - self.mAtom.mValCount
            if core_elec != 0:
                for shell in GlobalData.mGeo_Ladder[core_elec]:
                    if shell == 'point':
                        self.mCore.append(self.Point(self.mAtom))
                    elif shell == 'tetra':
                        self.mCore.append(self.Tetrahedron(10,self.mAtom))

        #Prepare the valence shell first, since it will help determine the
        # orientation of the inner shells
        if GlobalData.GetFullElecCount(self.mAtom.mGroup,self.mAtom.mPeriod) <= 18:
            for bond in self.mAtom.mBonds:
                if bond.mAtoms[1] > bond.mAtoms[0]:
                    if bond.mOrder == 1:
                        SingleBond(atoms[bond.mAtoms[1]])
                    elif bond.mOrder == 2:
                        DoubleBond(atoms[bond.mAtoms[1]])
                    elif bond.mOrder == 3:
                        TripleBond(atoms[bond.mAtoms[1]])
            #Add Free-Electrons
            if self.mAtom.mFreePairs == 2:
                if self.mAtom.mSteric >= 3:
                    AddFreeElectron(2)
            elif self.mAtom.mFreePairs == 1:
                if self.mAtom.mSteric >= 2:
                    AddFreeElectron(1)
            elif self.mAtom.mFreePairs == 3:
                AddFreeElectron(3)
        
        AddCoreElectrons()
             
    def FinalizeFODs(self):
        """
        Add all FODs in the FODStructure, to the atom, so that they can be
        printed when creating the XYZ output file
        """
        #Add the core FODs
        for shell in self.mCore:
            shell.mfods += self.mAtom.mPos   
            if self.mfods == []:
                self.mfods = [shell.mfods]
            else:
                self.mfods = np.vstack((self.mfods,shell.mfods))  ###HOW TO concatenate FODs, easily
    
    class FODShell:
        def __init__(self, shape, fods, owner: Atom):
            self.mOwner = owner
            self.mShape = shape 
            self.mfods = np.array(fods)
        def __str__(self):
            return self.mShape

    class Point(FODShell):
        def __init__(self, owner: Atom):
            super().__init__('point', [0.0,0.0,0.0], owner)

    class Tetrahedron(FODShell):
        """
        Tetrahedron Class: FOD Geometry corresponding to sp3 "hybridization' geometry.
        Roadmap: There will  be different functions to create compound transformations of FODs (e.g. the base, or peak
        of the tetrahedron), and to rotate them in the proper direction as well.
        """
        def __init__(self, core_amount, owner: Atom):
            super().__init__('tetra', GlobalData.mTetraGeo, owner)
            self.mfods *= GlobalData.mRadii[core_amount][self.mOwner.mZ]
        #Class Methods
        def CreateTetra(self):
            """
            This method arranges 4 FOD points in a tetrahedral form. Depending on the bonding, and free electrons,
            the direction vector will be different. 
            """
            pass

        def RotateTetra(self):
            pass
