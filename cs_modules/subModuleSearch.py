#!/usr/bin/python
import numpy as Np

class plvModule():
    '''
    subPatterns in this contactPack module extracts number of sub structures from evaluated contact string.

    Call in following order to get sub patterns of a contact string.

    q = plvModule()   # create instance of class subPatterns
    pat = 'Ea.Er.0-Ea.0-Ep' # contact string in a variable or you can directly pass to method
    q.feedMatrixOrientation(pat)    # decompose to contact matrix  [sse * sse]
    print q.sseCount                # number of secondary structures in the given contact string
    q.getSSCombPatt()               # get all possible substructure combination patterns from N to C
    #print q.subPatterns

    '''
    def __init__(self):
        self.contString = ''
        self.sseCount = 0
        self.mat = Np.array([])
        self.subPatterns = {}
    
    def feedMatrixOrientation(self,patternType):
        '''
        Build matrix from contact pattern. Build the object of contact and call feedMatrix 
        with contact pattern as argument.
        '''
        self.contString = patternType   
        tmp = patternType.split("-")
        #print tmp  
        self.sseCount = len(tmp)+1
#       mat = Np.zeros(shape = (len(tmp)+1,len(tmp)+1)) # Build matrix of helix * helix dimension
        mat = Np.chararray(shape = (len(tmp)+1,len(tmp)+1), itemsize=2) # definatino of item size in important
        # Enter Contact detail for matrix elements#######
        for d in range(len(tmp)):
            k = 0
            ltmp = tmp[d].split(".")
            #print ltmp
            for m in range(Np.shape(mat)[0]-1):
                dist = d+1
                if m+dist <= (Np.shape(mat)[0]-1):
                    #mat[m,m+dist] = tmp[k][d]
                    mat[m,m+dist] = ltmp[k]
                    #print mat[m,m+dist],ltmp[k]
                    k += 1
                else:
                    break
        self.mat = mat

    def __moduleValidity__(self,cs):
        csSplit = cs.split("-")
        step = 0
        exclusiveSets = []
        for i in range(len(csSplit)):
            csStep = csSplit[i].split(".")
            step += 1
            x = 0
            for j in range(len(csStep)):
                if csStep[j] != '0':
                    exclusiveSets.extend([j,j+step])
                else:
                    pass
        exclusiveSets = set(exclusiveSets)

        if len(exclusiveSets) == len(csSplit)+1:
            return 1
        else:
            return 0


    def getSSCombPatt(self):
        '''
        call after passing contact string to feedMatrixOrientation('<pattString>')
        parse the contact matrix into different subpatterns
        '''
        cmbHlx = {} # dictionary to store all possibilities of SSE combination
        for h in range(3,self.sseCount):
            cmbHlx[h] = []
            self.subPatterns[h] = []
            initialV = 0
            for i in range(self.sseCount):
                lastVal = initialV+h
                if lastVal <= self.sseCount:
                    cmbHlx[h].append(range(initialV,lastVal)) # Possible sequencial helix combinations
                    initialV += 1
                else:
                    break
        
        subPatterns = {}
        # Build patterns of helix combinations
        for k1 in cmbHlx.values():
            for k in k1:
                contStrng = []
                for step in range(1,len(k)):
                    i = 0
                    if len(contStrng):
                        contStrng[-1] = contStrng[-1].rstrip(".")
                        contStrng.append("-")

                    for j in range(i+step,len(k)):
                        contStrng.append(str(self.mat[k[i],k[j]])+'.')
                        i += 1
                stng = ''.join(contStrng)
                stng = stng.rstrip(".")
                if self.__moduleValidity__(stng):
                    self.subPatterns[len(k)].append(stng)   # append contact string to corresponding HGs
                else:
                    pass
        # import ipdb; ipdb.set_trace();
        #print self.subPatterns
    def help(self):
        print """
        q = plvModule()
        q.feedMatrixOrientation('pattern')
        q.getSSCombPatt()
        print q.subPatterns # A list of list subpatterns wirh different sse numbers sequenntially 
        """            
#------------------------------------------------------------------------------------------------------------####  

def main():
    import sys
    print "#############  Welcom to ProLego: SubModule Search  #############\n"

    import argparse
    parser = argparse.ArgumentParser()
    try:
        parser.add_argument('-pat',action='store',dest='contStrng',help="contact Strng as input")
        presult = parser.parse_args()
        q = plvModule()
        pat = presult.contStrng
        q.feedMatrixOrientation(pat)
        print "SSE group  :%d"%q.sseCount
        q.getSSCombPatt()
        print q.subPatterns
    except:
        print ("Provide argumrnt. Try:")
        print (" python subModuleSearch.py -pat 'Ea.Er.0-Ea.0-Ep' ")


if __name__ == "__main__":
    main()
   