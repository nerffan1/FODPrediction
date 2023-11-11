vector_for_cross = []
fulle = GlobalData.GetFullElecCount(self.mAtom.mGroup, self.mAtom.mPeriod)

if len(self.mAtom.mBonds) > 1 :
    for otherb in self.mAtom.mBonds:
        vector_for_cross.append(atoms[otherb.mAtoms[1]].mPos)
        #Find the cross term, to find the perpendicular vector
    vector_for_cross -= self.mAtom.mPos
elif len(self.mAtom.mBonds) == 1:
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
    if len(self.mAtom.mBonds) > 1 :
        axis2fod = np.cross(*vector_for_cross)
        axis2fod /= np.linalg.norm(axis2fod)
        #Add both FODs of the Double Bond
        dr = self.mAtom.mPos - np.sum(vector_for_cross, axis=0)*.2
                   
    elif len(self.mAtom.mBonds) == 1:
        #Similar to DoubleBond placement
        #TODO: If theres no other FODs, then you must relegate to later
        for fods in self.mValence:
            vector_for_cross.append(fods)
        vector_for_cross -= self.mAtom.mPos
        axis2fod = np.cross(*vector_for_cross)
        
        #Add both FODs of the Double Bond
        axis2fod *= GlobalData.mVert[fulle][self.mAtom.mZ]/2
        dr = self.mAtom.mPos + free_dir*np.linalg.norm(axis2fod)/np.tan(np.deg2rad(54.735))
    
    self.AddValFOD(dr + axis2fod,False,True)
    self.AddValFOD(dr - axis2fod,False,True)

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