#Description: The GlobalData class is dedicated to holding any static information that will be used
#  in the heuristic method of FODs. Some of the data it includes are:
#  - Number of electrons in closed shells
#  - Access to elements information found in the elements2.0 file
#RoadMap: Will eventually require access to average FOD-Tetrahedra radii (along those lines...), so that 
#  we can have a good guess at their positions in space. If there must be hybridization of FODs and Tetrahedra
#  combine, then we can also combine the shared FODs (i.e. the bonding FODs), and use the monoatomic data for reference.
#  Currently more time to look into heuristics is necessary.
#Author: Angel-Emilio Villegas S.
from numpy import where, genfromtxt

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

    #Class Variables    
    mElementInfo = []
    mElemNames = []
    mClosedGroups = [2,12,18]
    mLadder_3p = [4,10,12,18]
    #This scheme below assumes that a 1s Core FOD will already have been placed in the FODStruct 
    mLadder_3p_geometry = { 4: ['offcenter'], 10: ['tetra'], 12: ['tetra', 'offcenter'], 18: ['tetra', 'tetra'] }
    mShellShapes = {1: ['point'], 4: ['tetra'], 9: ['triaugmented']}

    #Geometries
    
