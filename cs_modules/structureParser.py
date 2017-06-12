#!/usr/bin/python

import sys
import re
residueNameMap_3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','ASX':'B',\
                        'CYS':'C','SCY':'C','GLU':'E','GLN':'Q','GLX':'Z','GLY':'G',\
                        'HIS':'H','ILE':'I','LEU':'L','LYS':'K',\
                        'MET':'M','PHE':'F','PRO':'P','SER':'S',\
                        'THR':'T','TRP':'W','TYR':'Y','VAL':'V',\
                        'UNK':'X'}
# residueNameMap_3to1 : ASX : ASparagine or aspartic acid; GLX: Glutamine or Glutamic acid
# 'SCY': S-Acetyl-Cysteine 'C', charge 0

def usage():
    print """
    structrueParser classes can parse PDB and DSSP
    files important for Contact String Analysis. This is a
    part of ProLego_CS_v3 suit.

    start date: 12th July, 2016
    author    : taushif

    1. class PDB:
    get instance of the class
    read pdbFile (readPDB)
    (parserPDB) parse pdb file to get chain information of first model (incase of NMR)
    (chainParser) get atomlines and residue information from a chain
    (getatom)    coordinate information of each residue

    2. class DSSP:
    runDssp
    parseChain


    """

class PDB():
    """PDB structure parser for the code
    start the instance with file path PDB(<file path>)
    call following with instances:
        1. readPDB() --> return status of file read
        2. parsePDB() --> Parse into models [populates chain details ins dictionary chainDetails]
        3. chainParse(Chain) --> populates self.resDict   = {} self.atomLines = [] self.resCord   = {}

    Functions:
        1. getatom(self, resId) : doc in respective def
    """
    def __init__(self,fileName):
        """
        Parsing a PDB formatted file call PDB(<filename>)
        call readPDB and parsePDB
        """
        self.fileName = fileName
        self.pdb_file = []
        self.chainDetail = {}


    def readPDB(self):
        """
        read pdb file given as input
        """
        try:
            with open(self.fileName,'r') as f:
                self.pdb_file = [i.strip() for i in  f.readlines()]
            print ("PDB file has read successfully [Okay]")
            return 1
        except:
            print ("Error in File Read: Check Uploaded PDB")
            return 0

    def parserPDB(self):
        """
        parse PDB FIle
        """
        print "Follows a strict PDB format for Chain and Model Recognition.\n"
        tmp = []
        if self.readPDB():
            for li in range(len(self.pdb_file)):
                l = self.pdb_file[li]
                if l[0:6].strip() == 'ATOM':
                    tmp.append(l.strip())
                if l[0:6].strip() == 'TER':
                    chain = self.pdb_file[li-1][21]
                    self.chainDetail[chain] = tmp
                    tmp = []
                if l[0:6].strip() == 'END' or l[0:6].strip() == 'ENDMDL':
                    break

            print "Total Chain = %d "%(len(self.chainDetail.keys()))
            return 1

        else:
            print ("Could not follow next")
            return 0

    def chainParser(self,chain):
        """Chain Parser : details of a particular chain
        Args:
            chain (str): chain of a pdb chain can be found from chainDetails
        Populate the attributes
        """
        self.resDict   = {}
        self.atomLines = []
        self.resCord   = {}

        if chain in self.chainDetail.keys():
            for atom_line in self.chainDetail[chain]:
                resId = int(atom_line[22:26].strip())
                atmId = atom_line[12:16].strip()

                self.atomLines.append(atom_line)
                self.resDict[resId] = atom_line[17:20].strip()
                if resId in self.resCord.keys():
                    self.resCord[resId][atmId] = (float(atom_line[30:38]),float(atom_line[38:46]),float(atom_line[46:54]))
                else:
                    self.resCord[resId] = {}
                    self.resCord[resId][atmId] = (float(atom_line[30:38]),float(atom_line[38:46]),float(atom_line[46:54]))
            return self.resDict

        else:
            print "Error 1 !! Chain Not found."
            return 0

    def getatom(self,resId):
        """ getatom : given a residue Id returns (x,y,z) coordinates
        @Args:
        resId (int) --> residue number in the PDB
        @return:
        atomCoord (dict) --> atom wise x,y,z coordinates

        """
        atomCoord = {}
        for i in self.atomLines:
            if int(i[22:26].strip()) == resId:
                atomCoord[i[12:16].strip()] = (float(i[30:38]),float(i[38:46]),float(i[46:54]))
        if atomCoord.keys():
            return atomCoord
        else:
            print resId," is not a valid residue Number [Error 1]"
            return 0

    def write_PDBChain(self,chain):
        """writing coordinates of a chain
        @Args:
        chain (str) --> chain of a PDB file
        @return:
        status (0/1) --> fail/Success
        pdbCh_name (str) --> detail file name of the chain written
        """
        pdbCh_name = '%s_%s.pdb'%(self.fileName,chain)

        pFl = open(pdbCh_name,'w')
        if chain in self.chainDetail.keys():
            for l in self.chainDetail[chain]:
                if re.match(r'%s'%l[12:16].strip(),'H'):
                    pass
                else:
                    pFl.write("%s\n"%l)
            pFl.close()
            print ("Chain Written in %s"%pdbCh_name)
            return 1,pdbCh_name
        else:
            print "Error 1: Chain Not Present, could not write"
            return 0,'0'

    def get_PDBFasta(self,pdbId='prot_name',chain='chain_',comment='PDB Fasta',outType=1):
        """sequence in fasta format from PDB information
        Note: from the residues of a PDB file writes the fasta
        sequence.
        @Args:
        pdbId (str)-> Id of the pdb File 4 letter (default: prot_name)
        chain (str) -> chain Id (default: chain_)
        comment (str) -> comment about the structre

        @return:
        fastaString : string character '>pdbId:Chain|residues\nSEQUENCE'
        """
        self.fastaSeq = ''
        fastaHeader = ">%s_%s|%d|%s"%(pdbId,chain,len(self.resDict.keys()),comment)
        if self.resDict.keys():
            for k in self.resDict.keys():
                try:
                    self.fastaSeq += residueNameMap_3to1[self.resDict[k]]
                except 'KeyError':
                    print k,self.resDict[k],"Not found in List, \n Aborting Updatelist"

        else:
            print "residue dictionary is empty. see the PDB docstring"
        print "Got sequence of length %s"%len(self.fastaSeq)
        if outType:
            fastaString = "%s\n%s"%(fastaHeader,self.fastaSeq)
            return fastaString
        else:
            return self.fastaSeq


    def _getSSE_(self):
        """get secondary structure in formation from structur file
        Note: parse sse information from REMARKS

        """
        self.sse = {'HELIX':[],'SHEET':[],'TURN':[]}
        x = 0
        for l in self.pdb_file:
            if l[0:5].strip() in self.sse.keys():
                self.sse[l[0:6].strip()].append(l.strip())
                x = 1
            else:
                continue
        if x:
            return self.sse
        else:
            print ("Remarks seems tobe absent in th PDB file")
            return 0

    def __sseHelixGroup__(self,hlx):
        hlxGroup = {1:'rAlp',2:'rOmg',3:'rPi',4:'rGam',5:'r3t',6:'lAlp',7:'lOmg',8:'lGam',9:'ribb',10:'ppII'}
        hlxChainWise = {}
        for h in hlx:
            chain = h[19]
            helId = 'H'+h[7:10].strip()
            n_res = int(h[21:25].strip())
            c_res = int(h[33:37].strip())
            hgrp = int(h[38:40].strip())

            if chain in hlxChainWise.keys():
                hlxChainWise[chain].update({n_res:[helId,n_res,c_res,hlxGroup[hgrp]]})
            else:
                hlxChainWise[chain] = {}
                hlxChainWise[chain].update({n_res:[helId,n_res,c_res,hlxGroup[hgrp]]})
                # hlxChainWise[chain].append([])
        return hlxChainWise

    def __sseSheetGroup__(self,sheet):
        shtChainWise = {}
        for strnd in sheet:
            chain = strnd[21]
            strnId = 'E'+strnd[7:10].strip()
            n_res  = int(strnd[22:26].strip())
            c_res = int(strnd[34:37].strip())
            sheetName = strnd[11:14].strip()

            if chain in shtChainWise.keys():
                shtChainWise[chain].update({n_res:[strnId,n_res,c_res,sheetName]})
            else:
                shtChainWise[chain] = {}
                shtChainWise[chain].update({n_res:[strnId,n_res,c_res,sheetName]})
        return shtChainWise



    def pdbSSE(self):
        """
        P = PDB(pfile)
        sseChain = P.pdbSSE()

        @paaram: path of pdb file from RCSB [all detail]
        @return : Secondary structue information chain wise

        """
        self.readPDB()
        try:
            sse = self._getSSE_()
            hlxChain = self.__sseHelixGroup__(sse['HELIX'])
            shtChain = self.__sseSheetGroup__(sse['SHEET'])
            chains = list(set(hlxChain.keys()+shtChain.keys()))
            chainSSE = {}
            for c in chains:
                if c in hlxChain.keys():
                    t = hlxChain[c].copy()
                    if c in shtChain.keys():
                        t.update(shtChain[c])
                    else:
                        pass
                else:
                    t = shtChain[c].copy()

                chainSSE[c] = t
            return chainSSE
        except:
            print ("Error in PDB File ")
            return 0

class DSSP:
    def __init__(self):
        self.pdbChain = ''
        self.dsspLines = []
        self.alInfo = []
        self.majSSE= [] # list of SSE H/E in the sequenctial order
        self.pdbid = ''
        self.sse_index = {'H':[],'E':[],'I':[],'G':[],'B':[],'S':[],'T':[],'C':[]}
        self.resCalpha = {}
        self.dsspToPdbMap = {}
        self.dihedrals = {}

    def runDssp(self, pdbChain):
        self.pdbChain = pdbChain
        import subprocess
        outFile = '%s.dssp'%(self.pdbChain)
        dssp = 'dssp -i %s -o %s'%(self.pdbChain,outFile)
        try:
            print dssp
            subprocess.Popen(dssp, shell= True).wait()
            print "DSSP successfully created"

            return 1, outFile
        except:
            print ("Error 1.2:  in running dssp")
            return 0, 0

    def readDSSP(self, dsspFl):
        """
        resding dssp file return 1 or 0
        """
        try:
            with open(dsspFl,'r') as f:
                self.dsspLines = [i for i in  f.readlines()]

            print ("DSSP file has read successfully [Okay]")
            return 1
        except:
            print ("Error 2.1: DSSP Error parsing unsuccesful\n")
            return 0


    def dsspParseChain(self):
        """
        @param
        f: filehandler of file in a list
        """
        start =0

        for l in self.dsspLines:
            chk=l[10:12].strip()
            if re.search('HEADER',l):
                self.pdbid=l[62:66].strip()
                continue
            if re.search('#',l):
                start=1
                continue

            if start:   #and chk == chain:
                tmpres = {}
                try:
                    #pdb residue number
                    residue = int(l[5:10].strip())
                    self.dihedrals[residue] = {'kappa':0.0,'alpha':0.0,'phi':0.0,'psi':0.0}
                    #building empty dictionary for that residue
                    tmpres[residue] = {'aa': '', 'sse': '', 'acc': '', 'catom': ()}
                    # amino acid name
                    tmpres[residue]['aa'] = l[12:14].strip()
                    # sse type
                    sse_tmp = l[14:17].strip()
                    if re.match('[GHITEBS]', sse_tmp):
                        tmpres[residue]['sse'] = sse_tmp
                        self.sse_index[sse_tmp].append(residue)
                    else:
                        tmpres[residue]['sse'] = 'C'
                        self.sse_index['C'].append(residue)

                    #accessible solvent area    
                    tmpres[residue]['acc'] = int(l[34:38].strip())
                    # calpha x, y and z
                    x = float(l[117:122].strip())
                    y = float(l[124:129].strip())
                    z = float(l[131:].strip())
                    tmpres[residue]['catom'] = (x, y, z)

                    self.dihedrals[residue]['kappa'] = float(l[91:97].strip())
                    self.dihedrals[residue]['alpha'] = float(l[97:103].strip())
                    self.dihedrals[residue]['phi'] = float(l[103:109].strip())
                    self.dihedrals[residue]['psi'] = float(l[109:115].strip())
                    
                    # putting all in resdie information
                    self.resCalpha[int(l[5:10].strip())] = (x,y,z)
                    self.alInfo.append(tmpres)
                    self.dsspToPdbMap[int(l[5:10].strip())] = int(l[:6].strip())

                except:
                    continue

            elif chk != chk:
                start = 0

    def get_cacoord(self,resid_list):
        calpha_residues = []
        for res in resid_list:
            calpha_residues.append(self.resCalpha[res])
        return calpha_residues


    def getOnlyAlpBeta(self,sse,res,c):
        tmp = {}
        if re.match(r'[H]',sse):    # corrected Helix definition only to alpha helix (H) to NOT HGI
            c += 1
            strI = 'H'+str(c)
            tmp[strI]=res
            self.majSSE.append(tmp)
            return (c)
        elif re.match(r'E',sse):
            c += 1
            strI = 'E'+str(c)
            tmp[strI]=res
            self.majSSE.append(tmp)
            return (c)
        else:
            tmp['C']=res
            self.majSSE.append(tmp)
            return(c)

    def getMajorSSE(self):
        '''This module Stores different SSE information in separate lists
        to a key value. Input to this method is SSE list and output will be
        a dictionary.
        li is the list of sse information from and res is the residue info'''

        # declaring has and list

        sseInfo={'H':[],'G':[],'I':[],'T':[],'E':[],'B':[],'S':[],'C':[]}
        listOfsse = []
        res = []
        li = [] # sse information
        for r in self.alInfo:
            resNo = r.keys()[0]
            res.append(resNo)
            li.append(r[resNo]['sse'])

        hd=[]
        counter = 0
        # import ipdb; ipdb.set_trace();
        for i in range(len(li)-1):              # iterating over end-1 elemets of list
            if li[i] != li[i+1]:
                hd.append(res[i])
                # print hd
                sseInfo[li[i]].append(hd)
                counter = self.getOnlyAlpBeta(li[i],hd,counter)
                hd=[]
            if i+1 == len(li)-1:
                hd.append(res[i+1])
                sseInfo[li[i+1]].append(hd)
                counter = self.getOnlyAlpBeta(li[i+1],hd,counter)
            else:
                if i+1 == len(li)-1:
                    hd.append(res[i+1])
                    sseInfo[li[i]].append(hd)
                    counter = self.getOnlyAlpBeta(li[i],hd,counter)
                else:
                    hd.append(res[i])
        return sseInfo


    def get_SSEchunks(self):
        """
        get chunk of sse from the dictionary
        sse:{'H':[list of residue List],
        'E':[],
        'C':[],'h_310':[],'h_pi':,'B',[],'S':[]}
        """
        from operator import itemgetter
        from itertools import groupby
        sseGroup = {'H':[],'E':[],'I':[],'G':[],'B':[],'S':[],'T':[],'C':[]}
        # get the mapping
        for sse in self.sse_index.keys():
            sseList = self.sse_index[sse]
            #The key to the solution is differencing with a range so that consecutive 
            # numbers all appear in same group.
            for k, g in groupby(enumerate(sseList), lambda (i,x):i-x):
                sseGroup[sse].append(map(itemgetter(1), g))
        return sseGroup

def __sortArray__(x):
    sseSec = []
    tmp = []

    for i in range(len(x)-1):
        if x[i]-x[i+1] == -1:
            tmp.append(x[i])
            if i == len(x)-2:
                tmp.append(x[i+1])
                sseSec.append(tmp)
        else:
            tmp.append(x[i])
            sseSec.append(tmp)
            tmp = []
    return sseSec


def sseSeqToChunks(sseString):
    """
    @param = "HHHHHHHHHHHHTTTTTTEECCEEETTTEEEESSEETTTEECCCCCCBCCCCCSSCCCCTTEEEE"
    @retrn = sseChunk{'H':[[0,1,2,..],[...]],'E':[[x.x],[y,yx]]}

    """
    sse = {'H':[],'E':[]}
    sseChunks = {'H':[],'E':[]}
    for l in enumerate(sseString):
        if l[1] in sse.keys():
            sse[l[1]].append(l[0])
    for k in sse.keys():
        sseChunks[k] = __sortArray__(sse[k])
    return sseChunks

def getMapDict(sseChunkDict):
    sc = {}
    scMap = {}
    for k in sseChunkDict.keys():
        for l in sseChunkDict[k]:
            scMap[l[0]] = k
            sc[l[0]] = l
    return scMap,sc


def main(pfile):
    D = DSSP()
    D.readDSSP(pfile)
    D.dsspParseChain()
    
    print D.alInfo

def mainPDB(pfile):
    P = PDB(pfile)
    sse =  P.pdbSSE()


if __name__ == '__main__':
    import sys
    import os
    cwd = os.getcwd() + "/"
    fl  = cwd+sys.argv[1] 
    print fl
    mainPDB(fl)
    # main(fl)