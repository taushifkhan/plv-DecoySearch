# -*- coding: utf-8 -*-
# @Author: taushif
# @Date:   2017-03-28 08:28:47
# @Last Modified by:   taushif
# @Last Modified time: 2017-03-28 16:57:45
#!/usr/bin/python
"""
decoy filtering and classification. Given all the decoys ina directory with separate PDB coordinate files.
MakeContactString module analyzes all coordinates and generate contact string for each conformer.

The present module takes a contact string and searched in the current availabe database of pLv_topoStatus table. 
If the exact match found, the out put shows the topological status with the option of component search. 
If an exact match could not be established, then contact string is subjected to component search. 
In this search all the smaller structural modules or legos of contact string has been extracted and then 
topology status has been performed.
"""
import pickle

class csSearch():
    """docstring for csSearch

    """
    def __init__(self,topoDict):
        self.pklTable = topoDict

    def __checkCS__(self,cs):
        cs_alt = cs.split("-")
        # import ipdb; ipdb.set_trace()
        if len(cs_alt) >= 2:
            return 1
        else:
            return 0

    def getFromTopoStat(self, cs):
        """
        @args:
        cs (str) -> contact string
        """
        outRes = {}

        if self.__checkCS__(cs):
            if cs in self.pklTable.keys():
                tres   = self.pklTable[cs]
                tres.update({'cs':cs})
                # import ipdb; ipdb.set_trace();
                print (" !!! %s  Found in TopoStatus"%cs)
                return (1,tres)
            else:
                # print ("%s Not Found in TopoStatus"%cs)
                outRes['cs'] = cs
                outRes['topoStatus'] = 'NotClassified'
                return (0,outRes)
        else:
            print ("cs %s doses not satisfy minimum requirement" %cs)
            return (0, outRes)


    def getComponenet(self,cs):
        """
        can be called from the csClass incident to identify all components
        in a topology and further classififcation
        """
        import subModuleSearch as cP
        cPp = cP.plvModule()
        if self.__checkCS__(cs):
            try:
                lego_blocks = {}
                cPp.feedMatrixOrientation(cs)
                cPp.getSSCombPatt()
                legos = cPp.subPatterns
                
                print ("COmponents Generated of Now searching TopoStatus DB")
                
                for k in legos.keys():
                    print ("Searching component of %d ..."%k)
                    for cs_module in legos[k]:
                        com_result = self.getFromTopoStat(cs_module)
                        if com_result[0]:
                            if k in lego_blocks.keys():
                                lego_blocks[k].append(com_result[1])
                            else:
                                lego_blocks[k] = []
                                lego_blocks[k].append(com_result[1])
                        else:
                            # lego_blocks[k].append(com_result[1])
                            print ("Compoent not found in TopoStatus")
                return (1, lego_blocks)

            except:
                print ("Error !!! in component identification")
                return (0,0)
        else:
            print ("cs %s doses not satisfy minimum requirement" %cs)
            return (0,0)

    def pLv_compare(self,cs):
        """
        @args:
        cs (string) -> contact string
        @return:
        topology status details

        """
        moduleSearch = {'topology':{},'modules':{}}

        for i in moduleSearch.keys():
            if i == 'topology':
                result = self.getFromTopoStat(cs)

                if result[0]:
                    print ("Found from TopoStatus database")
                    moduleSearch['topology'] = result[1]
                else:
                    moduleSearch['topology'] = result[1]

            elif i == 'modules':
                comResult = self.getComponenet(cs)
                if comResult[0]:
                    print ("components have been analyzed and listed")           
                    moduleSearch['modules'] = comResult[1]
                else:
                    print ("Components could not found\n")       

        return moduleSearch




def test():
    import sys
    cs = sys.argv[1]
    # try:
    topoDict = pickle.load(open('../DB/plVdb_topoStatus.pkl','rb'))
    # import ipdb; ipdb.set_trace();
    c = csSearch(topoDict)
    tStatus = c.pLv_compare(cs)
    import ipdb; ipdb.set_trace();
    # except:
    #     print cs
        


if __name__ == '__main__':
    test()
    # test_onlyComponent()
    # test_decoyOnlyTopoStat()


