#Description: The GlobalData class is dedicated to holding any static information that will be used
#  in the heuristic method of FODs. Some of the data it includes are:
#  - Number of electrons in closed shells
#  - Access to elements information found in the elements2.0 file
#RoadMap: Will eventually require access to average FOD-Tetrahedra radii (along those lines...), so that 
#  we can have a good guess at their positions in space. If there must be hybridization of FODs and Tetrahedra
#  combine, then we can also combine the shared FODs (i.e. the bonding FODs), and use the monoatomic data for reference.
#  Currently more time to look into heuristics is necessary.
#Author: Angel-Emilio Villegas S.
from numpy import where, genfromtxt, sqrt

class GlobalData:
    
    def __init__(self):
        GlobalData.mElementInfo = self.LoadElements()
        GlobalData.mElemNames = self.LoadNames()
    
    #Public Functions
    def LoadElements(self):
        with open('elements2.0', mode ='r') as file:
            return genfromtxt(file, delimiter=',', encoding=None, dtype=None)

    def LoadNames(self):
        """
        Names are loaded in their abbreviated form, since XYZ files receive them this way
        """
        return GlobalData.mElementInfo[1:,2]

    def GetElementAtt(name: str, attr: str):
        """
        Gets the attribute found in the elements2.0 file

        Parameters:
        name (str): The name of the atom, in abbreviated form.
        attr (str): The attribute you want of chosen atom, found in the first row of elements2.0

        """
        if attr == "AtomicNumber":
            attrib_i = where(GlobalData.mElemNames == name)  
            # We must offset index by +1 because the array starts at zero instead of 1
            return attrib_i[0].item() + 1      
        else:
            element_i = where(GlobalData.mElemNames == name)[0]
            attrib_i = where(GlobalData.mElementInfo[0] == attr)[0]
            #We offset by 1 because the Info list has attribute names on row 1
            return GlobalData.mElementInfo[element_i + 1,attrib_i]
        
    def GetZAtt(Z: int, attr: str):
        attrib_i = where(GlobalData.mElementInfo[0] == attr)
        #We offset by 1 because the Info list has attribute names on row 1
        return GlobalData.mElementInfo[Z,attrib_i[0].item()]

    #Degugging tests
    def _debug_samplenames():
        for att in ["AtomicNumber", "Group", "Period", "Element","Metal", "NumberofShells","NumberofValence"]:
            print(f'{att}: {GlobalData.GetElementAtt("Ga", att)}')
            print(GlobalData.GetZAtt(31, att))

    ############Class Variables############    
    mElementInfo = []
    mElemNames = []
    mClosedGroups = [2,12,18]

    # The following ladder is based of various monoatomic calculations.
    mGeo_Ladder = { 2: ['point'], 
                        4: ['point','point'],
                        10: ['point','tetra'], 
                        12: ['point','tetra', 'point'],
                        18: ['point', 'tetra', 'tetra'],
                        20: ['point', 'tetra', 'tetra', 'point'], 
                        30: ['point', 'tetra', 'tetra', 'triaug', 'point'],
                        36: ['point', 'tetra', 'tetra', 'triaug', 'tetra'],
                        54: ['point', 'tetra', 'triaug', 'triaug', 'tetra'] }
    mShellShapes = {1: ['point'], 4: ['tetra'], 9: ['triaugmented']}
    #This ladder is based of the 
    mElecConfLadder= [2,2,6,2,6,2,10,6,2,10,6,2,10,6]
    #Geometries for known shell structures
    mTetraGeo = [[0.0,0.0,1.0],
                [sqrt(8/9), 0.0, -1/3],
                [-sqrt(2/9),sqrt(2/3), -1/3],
                [-sqrt(2/9),-sqrt(2/3), -1/3]]
    
    #Molecules
    mAtoms = []