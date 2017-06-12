# -*- coding: utf-8 -*-
# @Author: taushif
# @Date:   2017-03-28 08:49:06
# @Last Modified by:   taushif
# @Last Modified time: 2017-03-28 17:02:09

import sys
import os
import pickle
import parameterSetting as pS

sys.path.append(pS.csDir)
from decoySearch import csSearch 

topoDict = pickle.load(open(pS.DBPathTop,'rb'))

import argparse


def saveAsJson(jsonDir,cs_result,pid):
    """
    @param:
    jsonDir (str, PATH) = resultant directory as given in appSetting.py
    resDict (dict) = dictionary tobe wtitten in the json file
    pid  (str)    = pid for output file
    chain (str)  = chain of the protein
    @returns:
    1 if successfully written or 0
    """
    import json
    resFile = jsonDir+'%s_pLv.json'%(pid)
    try:
        json.dump(cs_result, open(resFile,'w'),indent=4,encoding='utf-8')
        print ("Written successfully in File: %s\n"%resFile)
        return 1
    except:
        print ("Error: 4.1: could not write json file")
        return 0

class pLv():
    """pLv doc string proLego Version
    Initiate pLv() instance by importing proLegoV3_main

    example:
    pL = pLv()
    pL.main(pdbFl, chain, pid) : writes conatca string in json format.

    """
    def __init__(self):
        self.pdbFl = ''
        self.chain = ''
        self.pid = ''

    def __cs4Chain__(self, chain, pD,cs):
        self.pLvresult = {'pdbId': '','chain': '','pdbinfo':{},'patern':{},'csModules':{}}
        k = chain
        print "Analyzing %s "%k
        cs.csGen(pS.logDir, self.pid, k, pD)
        self.pLvresult['patern'] = cs.resultCS

        cs_string = cs.resultCS['csNorm']
        # self.pLvresult['residues'] = len(cs.resultCS['resDict'].keys())

        self.pLvresult['pdbinfo']['pdbid'] = self.pid
        self.pLvresult['pdbinfo']['chain'] = self.chain
        return cs_string


    def main(self,  pdbFl, chain, pid='prot'):
        """pLv.main(pdbFl,chain,pid)
        @args:
        pdbFl (str) : full path of pdb file tobe analuzed
        chain (str) : chain of pdb file tobe analyzed
        pid (str) : pdbId_C 4letter pdb id and one letter chain (standard)
        [default value is prot]

        @returns:
        pLv.pLvresult : dictionary of all the information can be accessed.
        writes in a json file with designated file path [see saveAsJson].
        """

        self.pdbFl = pdbFl
        self.chain = chain
        self.pid   = pid

        import makeContactString_v3 as mCS
        import structureParser as sP

        pD = sP.PDB(pdbFl)
        cs = mCS.CSgenEngine(pdbFl)
        chainDetails = cs.scanChain(pD)

        finalResult = {}

        
        c = csSearch(topoDict)

        if self.chain in chainDetails.keys():
            pidC = "%s_%s"%(pid,self.chain)
            csString = self.__cs4Chain__(self.chain, pD, cs)
            tStatus = c.pLv_compare(csString)
            self.pLvresult['csModules'] = tStatus
            finalResult[pidC] = self.pLvresult

        else:
            for k in chainDetails.keys():
                pidC = "%s_%s"%(pid,k)
                csString = self.__cs4Chain__(k, pD, cs)
                tStatus = c.pLv_compare(csString)
                self.pLvresult['csModules'] = tStatus

                finalResult[pidC] =  self.pLvresult


        jval = saveAsJson(pS.jsonDir,finalResult,self.pid)

        if jval:
            print ("pLv3 for %s has been complete:"%(self.pid))
        else:
            print ("json file cuold not built")


def allChain():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p',action='store',dest='pdbFl',help='pdbFile coordinates')
    parser.add_argument('-c',action='store',dest='chain',help='protein chain [Optional]', required=False)
    parser.add_argument('-id',action='store',dest='idname',help='pdbId to save results')

    presult = parser.parse_args()

    pdbFl   =  presult.pdbFl
    pid     = presult.idname

    if presult.chain:
        chain = presult.chain
    else:
        chain = ''

    print ("Inputs\nPDB:%s\nChain:%s\n"%(pdbFl,chain))
    pL = pLv()

    if os.path.isfile(pdbFl):
        # import ipdb; ipdb.set_trace();
        pL.main(pdbFl, chain, pid)
        print pL.pLvresult.keys()
    else:
        usage()
        print ("Error 0.1: Check PDB Files, not found")


if __name__ == '__main__':
    allChain()
