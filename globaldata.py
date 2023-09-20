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

    def GetFullElecCount(group: int, period: int):
        """
        Receive Group number and return the amount of electrons that the atoms contains
        when its valence is fully completed. E.g., Z = 6, then we will return 10 since 10
        closes the subshell in that row of the periodic table.
        TODO: Implement the metals     
        TODO: Implement 1s shell logic   
        """
        if group > 12 and group < 18: 
            if period == 2: return 10
            elif period == 3 : return 18
            elif period == 4: return 36    

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
    mTriPlane = [[0,sqrt(3)/2,0],
                 [-sqrt(3)/2,-sqrt(3)/2,0],
                 [-sqrt(3)/2,-sqrt(3)/2,0]]
    
    #Molecules
    #Radii of Tetrahedra obtained by closing the shells of 
    #several atoms to close the sp3 shell. E.g., for Z=5, 5 
    # electrons were added. 
    mRadii = {
        10: {
            5: 1.0724622240214408,
            6: 0.9326307459817831,
            7: 0.8245934058016624,
            8: 0.6807571999707002,
            9: 0.6127101178811892, 
            10: 0.6127101178811892,
            11: 0.5005890731191966,
            13: 0.3336349313415564,
            14: 0.29211200538406207,
            15: 0.2619356755640106,
            16: 0.24136512980895986,
            17: 0.2037282349524298,
            18: 0.2037282349524298,
            31: 0.09504603503116242,
            32: 0.09131627498496026,
            33: 0.08749593856607064,
            34: 0.0782336953726361
        },
        18: {
            13: 0.7103844384774737,
            14: 1.571379042838154,
            15: 1.179809249943448,
            16: 0.956759529738896,
            17: 0.7103844384774737, 
            18: 0.7103844384774737
        }
    }
    #Average Edge length of FOD-FOD distances in the outmost shell, of this amount of 
    mVert = {
        10: {
            5: 1.751301962745536,
            6: 1.5229774321557852,
            7: 1.3465512547238012,
            8: 1.1116706527406452,
            9: 1.0005509960170376,
            10: 1.0005509960170376,
            11: 0.8174598158123096
        },
        18: {
            13: 4.173993393062307,
            14: 2.5660512319934154,
            15: 1.9266204374514642,
            16: 1.5623817696036548,
            17: 1.1600529303222396,
            18: 1.1600529303222396
        }
    }
    mAtoms = []