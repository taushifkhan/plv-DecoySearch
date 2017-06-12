#!/usr/bin/python
"""
Modified in 17/3/2017
quicker calculation of distance matrix
in clusion of methods:
1.get_residuesContactMatrix
2.contactPresence

3. __
"""
import sys
import os
import numpy as Np
import scipy
from scipy.spatial.distance import cdist
import timeit

import parameterSetting as aS

sys.path.append(aS.csDir)

import structureParser as sP
import angleCalc as aC
# import topoGraph as tG


residueNameMap_3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','ASX':'B',\
                        'CYS':'C','SCY':'C','GLU':'E','GLN':'Q','GLX':'Z','GLY':'G',\
                        'HIS':'H','ILE':'I','LEU':'L','LYS':'K',\
                        'MET':'M','PHE':'F','PRO':'P','SER':'S',\
                        'THR':'T','TRP':'W','TYR':'Y','VAL':'V',\
                        'UNK':'X'}
# residueNameMap_3to1 : ASX : ASparagine or aspartic acid; GLX: Glutamine or Glutamic acid
# 'SCY': S-Acetyl-Cysteine 'C', charge 0

def usage():
    print '+'*80
    print"""
    generates contact string from given protein file and chain.
    usage: makeContactString_v3.py [-h] [-p PDBFL] [-c CHAIN] [-id IDNAME]

    optional arguments:
      -h, --help  show this help message and exit
      -p PDBFL    pdbFile coordinates
      -c CHAIN    protein chain
      -id IDNAME  pdbId to save results

    """
    print '+'*80

def usageAllChain():
    print '+'*80
    print"""
    Itirating over all chains in a PDB

    generates contact string from given protein file and chain.
    usage: makeContactString_v3.py [-h] [-p PDBFL] [-c CHAIN] [-id IDNAME]

    optional arguments:
      -h, --help  show this help message and exit
      -p PDBFL    pdbFile coordinates
      -id IDNAME  pdbId to save results

    """
    print '+'*80

atomVdwRad = {'C':1.70,'N':1.52,'O':1.55,'S':1.80,'H':1}
Threshold =0.6


def contactPresence(res1C,res2C):
    """
    minimum distance between two list of coordinate sets
    Atomic distance betwween two residues.
    @Args:
    res1C -> [coordinate list (x,y,z)]
    res2C -> [coordinate list (x,y,z)]
    @return
    dxy -> Minimum distance between two coordinate sets
    """
    x = Np.array(res1C.values())
    y = Np.array(res2C.values())
    dxy = Np.min(scipy.spatial.distance.cdist(x,y))
    return dxy

def get_residuesContactMatrix(PC):
    """
    N*N matrix formation for covering distances.
    """
    totalResidues = len(PC.resDict.keys())
    distMat = Np.empty((totalResidues,totalResidues),dtype=Np.float16)
    srted_resId = sorted(PC.resDict.keys())
    mapResId_matId = {}
    mapMatId_Res = {}

    for k in range(totalResidues):
        res1 = srted_resId[k]
        for l in range(k+1,totalResidues):
            res2 = srted_resId[l]
            distMat[k,l] = contactPresence(PC.resCord[res1],PC.resCord[res2])
        mapResId_matId[srted_resId[k]] = k
        mapMatId_Res[k] = srted_resId[k]

    return distMat, mapResId_matId, mapMatId_Res
 

def get_alphaBetaCount(totalSSE):
    # alp = re.compile(r'H')
    # beta = re.compile(r'E')
    sseType = {'H': 0, 'E': 0, 'seq': ''}
    for i in totalSSE:
        # if re.search(alp,i):
        if i[0] == 'H':
            sseType['H'] += 1
            sseType['seq'] += 'H'
        else:
            sseType['E'] += 1
            sseType['seq'] += 'E'
    return (sseType)



def getOrientation(sse1,sse2,dP):
    #print sse1
    cAlpS1 = dP.get_cacoord(sse1.values()[0])
    cAlpS2 = dP.get_cacoord(sse2.values()[0])
    aCC = aC.angleCalc()
    # angle,distance = aCC.getAngles(cAlpS1,cAlpS2)  # uses angleCal.calc_angles between sses
    angle = aCC.getAngles(cAlpS1,cAlpS2)  
    # import ipdb;ipdb.set_trace()

    if abs(angle) < 45:
        return ('p',angle)
    elif abs(angle) > 135:
        return ('a',angle)
    else:
        return('r',angle)


class CSgenEngine():
    def __init__(self,pdbFl):
        self.pdbFl = pdbFl
        self.resultCS = {'sseString':'','csNorm':'','resDict':{},'contactLog':[],'csFormattedFile':'','sseDetail':[]} # output return
        self.chainFl = ''
        self.dsspFl  = ''
        self.resCord = {}

    def scanChain(self,pD):
        '''
        @param:
        pD -> Instance of sP.PDB
        @return:
        Parsed instance of the chain
        '''
        if pD.readPDB():
            if pD.parserPDB():
                return pD.chainDetail
            else:
                return 0
        else:
            print "Could not parse [structureParser --> PDB --> readPDB]"
            return 0


    def prepareFiles(self,pD,dD, chain):
        '''
        @param:
        pD: Instance of sP.PDB
        dD: Instance of sP.DSSP
        @return: NULL

        build the pdb chain file in PDB format and runs dssp on that chain
        '''
        self.chain = chain
        print "Preparing Files for CS gen calculation"

        if pD.chainParser(self.chain):
            print len(pD.resDict.keys())
            chain_wv = pD.write_PDBChain(self.chain)
            if chain_wv[0]:
                self.chainFl = chain_wv[1]
                print self.chainFl

                dStat = dD.runDssp(self.chainFl)
                if dStat[0]:
                    self.dsspFl = dStat[1]
                    print self.dsspFl
                else:
                    print "Fatal Error 1.5: No dssp stop engine"
                    sys.exit(0)
            else:
                print "Fatal Error 1.4: chain has not be written "
                sys.exit(0)
        else:
            print "Error 1.3: chain parsing error\n"

        print "Chain: %s\nPDBChain Fl: %s\nDSSP File: %s\n"%(self.chain, self.chainFl, self.dsspFl)

    def gear1_SSEInfo(self,dD,pD):
        '''
        extract from the secondary structure information
        '''
        try:
            dD.readDSSP(self.dsspFl)
            dD.dsspParseChain()
            res = dD.alInfo
            sseI = dD.getMajorSSE()
            

            self.sseHE = []
            totalSSE = []

            for i in dD.majSSE:
                if i.keys()[0] != "C":
                    self.sseHE.append(i)
                    totalSSE.append(i.keys()[0])

            self.alpBeta = get_alphaBetaCount(totalSSE)
            self.resultCS['sseString']=self.alpBeta['seq']
            # import ipdb; ipdb.set_trace();
            
            for i in self.sseHE:
                fres = i.values()[0][0]
                lres = i.values()[0][-1]
                self.resultCS['sseDetail'].append({'label':i.keys()[0],'residues':\
                [(fres,pD.resDict[fres]),(lres,pD.resDict[lres])]}) # taushif
        except:
            print ("Error: 2.1: DSSP read error\n")
            sys.exit(0)

    def __checkSSEContact__(self, sse1, sse2, pD):
        contactArray = []
        s1_name = sse1.keys()[0]
        s2_name = sse2.keys()[0]
        
        s1_s = self.mapRestoMat[sse1.values()[0][0]]
        s1_e  = self.mapRestoMat[sse1.values()[0][-1]]

        s2_s = self.mapRestoMat[sse2.values()[0][0]]
        s2_e  = self.mapRestoMat[sse2.values()[0][-1]]

        contres = Np.where(self.distMat[s1_s:s1_e+1, s2_s:s2_e+1] <= 5)


        for i in range(len(contres[0])):
            #Reconstruct matrix Id
            mat_res1 = s1_s + contres[0][i]
            mat_res2 = s2_s + contres[1][i]
            # extract residue Ids
            pdb_res1 = self.mapMatId_Res[mat_res1]
            pdb_res2 = self.mapMatId_Res[mat_res2]

            res1Id = residueNameMap_3to1[pD.resDict[pdb_res1]]
            res2Id = residueNameMap_3to1[pD.resDict[pdb_res2]]

            contactArray.append((s1_name, pdb_res1, res1Id, s2_name, pdb_res2, res2Id, round(self.distMat[mat_res1,mat_res2],3)))
        # import ipdb; ipdb.set_trace();

        return contactArray


    def gear2_getCS(self,dD,pD,datalogDir,pdbId):
        """
        HH : min contact residues 2 for rest all minimum contact residues are 2
        @param:
        dD: Instance of dssp parser
        pD: Instance of PDB parser
        datalogDir : directory to write log file
        pid: protein Id to write in the file

        @return:
        complete dictionary of contact string

        """
        contType = {'HH':'H','EE':'E','HE':'C','EH':'C'}
        contres = {'HH': 3, 'EH': 2, 'HE': 2, 'EE': 2}
        contactStrngList = []
        sseHE = self.sseHE
        alpBeta = self.alpBeta

        step = 0
        #print "SSE pairs in contact"
        flog = open('%s%s_out.log'%(datalogDir,pdbId),'w')
        fresult = open('%s%s_contactResult.txt'%(datalogDir,pdbId),'w')

        # pdbid, chain,total SSE | H,E|seq of SSE in chain
        csFF_l1 = ">%s:%s:%d|%d:%d|%s;\n"%(pdbId, self.chain,len(sseHE),alpBeta['H'],\
            alpBeta['E'],alpBeta['seq'])

        flog.write(csFF_l1)
        fresult.write(csFF_l1)
        logArray = []

        startT = timeit.timeit()
        self.distMat, self.mapRestoMat, self.mapMatId_Res = get_residuesContactMatrix(pD)
        self.resultCS['resDict'] = pD.resDict
        elapseT = startT - timeit.timeit()
        print "Time for distance Calculations: (sec)", elapseT
        allContactArray = []

        for x in range(len(sseHE)):
            x = 0
            step +=1
            for y in range(x+step,len(sseHE)):
                # cont = checkContact(sseHE[x],sseHE[y],pD)
                cont = self.__checkSSEContact__(sseHE[x],sseHE[y],pD)
                conSSEs = sseHE[x].keys()[0][0]+sseHE[y].keys()[0][0]
                if len(cont) >= contres[conSSEs]:
                    typ = sseHE[x].keys()[0][0].strip()+sseHE[y].keys()[0][0].strip()
                    orient,angle = getOrientation(sseHE[x],sseHE[y],dD)
                    logArray.append([sseHE[x].keys()[0],sseHE[y].keys()[0],round(angle,3),"*Y"])
                    contactStrngList.append(contType[typ]+orient+'.')
                    allContactArray.extend(cont)

                else:
                    typ = sseHE[x].keys()[0][0].strip()+sseHE[y].keys()[0][0].strip()
                    orient,angle = getOrientation(sseHE[x],sseHE[y],dD)
                    logArray.append([sseHE[x].keys()[0],sseHE[y].keys()[0],round(angle,3),"N"])
                    contactStrngList.append('0'+'.')
                    if len(cont):
                        allContactArray.extend(cont)
                
                x += 1
            contactStrngList.append('-')

        self.resultCS['contactLog'] = allContactArray

        for i in range(1,len(contactStrngList)-1):
            if contactStrngList[i+1] == '-':
                contactStrngList[i] = contactStrngList[i][:-1]

        wCont = ''.join(contactStrngList)
        self.resultCS['csNorm'] = wCont.rstrip('-')

        csFF_l2 = "%s:%s|%s"%(pdbId,self.chain,wCont.rstrip('-'))
        fresult.write(csFF_l2)
        self.resultCS['csFormattedFile'] = csFF_l1+csFF_l2

        fresult.close()
        flog.close()
        
        return self.resultCS

    def cleanFiles(self):
        print "removing created chain files \n%s\n%s"%(self.chainFl,self.dsspFl)
        os.system('rm %s %s'%(self.chainFl,self.dsspFl))


    def csGen(self,datalogDir, pid, chain, pD):
        dD = sP.DSSP()
        self.prepareFiles(pD,dD,chain)
        self.gear1_SSEInfo(dD,pD)

        self.resultCS['Seq'] = pD.get_PDBFasta(outType=0)

        cs = self.gear2_getCS(dD, pD, datalogDir, pid)
        self.cleanFiles()

    def saveAsJson(self,resDir,pid,chain):
        import json
        resFil = resDir+'%s_%s_pLv.json'%(pid,chain)
        json.dump(self.res)




def main_csE():
    import argparse
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-p',action='store',dest='pdbFl',help='pdbFile coordinates')
        parser.add_argument('-c',action='store',dest='chain',help='protein chain')
        parser.add_argument('-id',action='store',dest='idname',help='pdbId to save results')

        presult = parser.parse_args()

        cwd     = os.getcwd()+'/'
        pdbFl   = cwd+presult.pdbFl
        chain   = presult.chain
        pid     = presult.idname
        jsonDir = cwd+'cs_Result'



        print ("Inputs\nPDB:%s\nChain:%s\n"%(pdbFl,chain))
        datalogDir = cwd+'log/'

        if not os.path.exists(datalogDir):
            print "creating log directory",datalogDir
            os.makedirs(datalogDir)

        # if not os.path.exists(jsonDir):
        #     print "Creating Dir to save results",jsonDir
        #     os.makedirs(jsonDir)

        if os.path.isfile(pdbFl):
            pD = sP.PDB(pdbFl)
            cs = CSgenEngine(pdbFl)
            chainDetails = cs.scanChain(pD)
            if chain in chainDetails.keys():
                cs.csGen(datalogDir,pid, chain, pD)
            else:
                print "Give Chain %s not found in %s"%(chain,pdbFl)
                print "Error 0.1: Provide correct chain"
                usage()

        else:
            usage()
            print ("Error 0.1: Check PDB Files, not found")
    except:
        print ("Error : give proper argument")
        usage()

def main_allChain():
    print "+"*80+'\n\n'+'='*25+'ProLego For All Chains'+'='*25
    print '\n'+"+"*80+'\n'

    import argparse

    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-p',action='store',dest='pdbFl',help='pdbFile coordinates')
        parser.add_argument('-id',action='store',dest='idname',help='pdbId to save results')

        presult = parser.parse_args()

        cwd     = os.getcwd()+'/'

        pdbFl   = cwd+presult.pdbFl
        pid     = presult.idname
        jsonDir = cwd+'cs_Result'

        print ("Inputs\nPDB:%s\n"%pdbFl)
        datalogDir = cwd+'log/'

        if not os.path.exists(datalogDir):
            print "creating log directory",datalogDir
            os.makedirs(datalogDir)

        # if not os.path.exists(jsonDir):
        #     print "Creating Dir to save results",jsonDir
        #     os.makedirs(jsonDir)

        if os.path.isfile(pdbFl):
            pD = sP.PDB(pdbFl)
            cs = CSgenEngine(pdbFl)
            chainDetails = cs.scanChain(pD)
            for k in chainDetails.keys():
                print "Analyzing %s "%k
                cs.csGen(datalogDir,pid, k, pD)
        else:
            usage()
            print ("Error 0.1: Check PDB Files, not found")
    except:
        print ("Error : give proper argument")
        usageAllChain()



if __name__ == '__main__':
    main_csE()
    # main_allChain()