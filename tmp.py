def SingleBond(at2: Atom):
            dom,weak,dir = DominantAtom(at1,at2,True)
            bfod = AxialPoint_Simple(dom, weak, dir)
            self.AddValFOD([dom.mPos + bfod], at2, True)
        
        def DoubleBond(at2: Atom, bond: Bond):
            """ 
            Create the FODs representing the double bond. Currently the FOD filling is unidirectional (sequential)
            and  does not account for the next atom in the iteration to see if we can further accomodate the bonding FODs 
            """
                                    
            def CrossPerpDir():
                vector_for_cross = []
                for otherb in self.mAtom.mBonds: 
                    if otherb != bond:
                        vector_for_cross.append(atoms[otherb.mAtoms[1]].mPos)
                vector_for_cross -= self.mAtom.mPos
                vec = np.cross(*vector_for_cross)
                return vec/np.linalg.norm(vec)
            
            def D_FFOD_Direction():
                if self.mAtom.mFreePairs == 0:
                    #Determine DFFOD Direction
                    if len(self.mAtom.mBonds) == 3: # Planar
                        return  CrossPerpDir()
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
            dom,sub,fugal = DominantAtom(at1,at2, True)

            if GlobalData.GetFullElecCount(self.mAtom.mGroup,self.mAtom.mPeriod) <= 18:
                #Find perpendicular unit vector
                if self.mAtom.mFreePairs == 0:
                    axis2fod = D_FFOD_Direction()
                  
                #Add both FODs of the Double Bond
                midpoint = AxialPoint_Simple(dom,sub,fugal)
                # Determine Vertical Projection of DFFOD
                axis2fod *= D_Bond_Height(dom, np.linalg.norm(midpoint))

                #Add FODs
                self.AddValFOD([dom.mPos + midpoint + axis2fod],at2,True)
                self.AddValFOD([dom.mPos + midpoint - axis2fod],at2,True)
    