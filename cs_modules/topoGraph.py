#!/usr/bin/python
'''
PreRequisite : networkx module of python
date : 1/6/15
change in topoGraph module : 
get_Nodes_Edges : assigns relative orienatation to sheets
        if in spatial contact append usin if condition either using
        elif condition.
        
'''
import numpy as Np

class contact:
    '''
    Generating distance matrix from conatcat string.
    '''
    def __init__(self):
        self.contString  = ''
        self.sseCount    = 0
        self.mat         = Np.array([])
        self.subPatterns = {}
    
    def feedMatrixOrientation(self,patternType):
        '''
        Build matrix from contact pattern. Build the object of contact and call feedMatrix 
        with contact pattern as argument.
        '''
        self.contString = patternType   
        tmp = patternType.split("-")
        self.sseCount = len(tmp)+1
        mat = Np.chararray(shape = (len(tmp)+1,len(tmp)+1), itemsize=2)
        for d in range(len(tmp)):
            k = 0
            ltmp = tmp[d].split(".")
            for m in range(Np.shape(mat)[0]-1):
                dist = d+1
                if m+dist <= (Np.shape(mat)[0]-1):
                    mat[m,m+dist] = ltmp[k]
                    k += 1
                else:
                    break
        self.mat = mat

    def getContactDict(self):
        print self.sseCount
        self.contactDict = {}
        for k in range(self.sseCount):
            #print k,k+1
            #for p in mat[k]:
            self.contactDict[k+1] = list(self.mat[k,k+1:self.sseCount])
        print self.contactDict

    def getFormattedJson(self):
        defPosition = {0:'L',1:'R',2:'+1',3:'-1'}
        self.allContact = {}
        for k in sorted(self.contactDict.keys()):
            tmp_dir = {}
            print self.contactDict[k]
            for pos in range(len(self.contactDict[k])):
                if self.contactDict[k][pos] != '0':
                    tmp_dir[defPosition[pos]]= k+1+pos

            self.allContact[k] = tmp_dir

        print self.allContact

    def get_Nodes_Edges(self):
        '''
        creates array of nodes and coresponding edges for a topology
        graph.
        assigns relative orienatation to sheets
        if in spatial contact append usin if condition either using
        elif condition.
        '''
        self.nodes = Np.arange(1,self.sseCount+1,1)
        self.edges = []
        self.orient = {}
        #print self.mat
        orientCol = {'a':'red','p':'green','r':'blue','na':'black'}
        orientStatus = 0
        lol = []

        for i in self.nodes:
            self.orient[i-1] = 2

        for r in range(Np.shape(self.mat)[0]):
            for c in range(r+1,Np.shape(self.mat)[1]):

                tmpProp = {'weight':1,'type':'sequential','color':orientCol['na']}
                tmpPropSp = {'weight':1,'type':'spatial','color':''}

                dist = abs(r-c)

                if self.mat[r,c] != '0' and dist > 1:
                    tmpPropSp['color']  = orientCol[self.mat[r,c][1]]         # Ha,Ea matrix element have sse contact and orientation.
                    tmpPropSp['weight'] = round(dist * 1.2,2)
                    
                    self.edges.append((r+1,c+1,tmpPropSp))
                                      
                    if self.mat[r,c][0] == 'E':
                        self.orient[r] = orientStatus
                        if self.mat[r,c][1] == 'a':
                            orientStatus += 1
                            orientStatus %= 2
                        self.orient[c] = orientStatus
                    #print r,c

                if self.mat[r,c][0] == 'E':
                    self.orient[r] = orientStatus

                    if self.mat[r,c][1] == 'a':
                        orientStatus += 1
                        orientStatus %= 2

                    self.orient[c] = orientStatus

                    tmpPropSp['color']  = orientCol[self.mat[r,c][1]]         # Ha,Ea matrix element have sse contact and orientation.
                    tmpPropSp['weight'] = round(dist * 1.2,2)

                    self.edges.append((r+1,c+1,tmpPropSp))
                    
                if dist == 1:
                    self.edges.append((r+1,c+1,tmpProp))
                     
        # import ipdb; ipdb.set_trace()
        #return lol
                   

    def get_networkx(self,sseDetail, pdbInfo):
        '''
        Created node and edges for a contact string and dumps into json
        object.
        sseDetail : detail of secondary structure information to be appended in nodees []
        '''
        import json 
        import networkx as nx
        from networkx.readwrite import json_graph
        G = nx.MultiGraph(pinfo = pdbInfo)

        G.add_nodes_from(self.nodes)
        G.add_edges_from(self.edges)
        #import ipdb; ipdb.set_trace();

        for n in G:
            G.node[n]["names"] = sseDetail[n-1]['label']
            G.node[n]["residues"] = sseDetail[n-1]['residues']
            if sseDetail[n-1]['label'][0] == "E":
                G.node[n]["orient"] = self.orient[n-1]
            #G.node[n]["orient"] = self.orient[n]

        # write json formatted data
        d = json_graph.node_link_data(G) # node-link format to serialize
        topoString = json.dumps(d,encoding='UTF-8',default=str)
        # import ipdb; ipdb.set_trace();
        return d

    def get_networkx_noSSE(self,sseString):
        '''
        Created node and edges for a contact string and dumps into json
        object.
        sseDetail : detail of secondary structure information to be appended in nodees []
        '''
        import json 
        import networkx as nx
        from networkx.readwrite import json_graph
        G = nx.MultiGraph(pinfo = '%s'%sseString)

        G.add_nodes_from(self.nodes)
        G.add_edges_from(self.edges)
        #import ipdb; ipdb.set_trace();

        for n in G:
            G.node[n]["names"] = sseString[n-1]
            if sseString[n-1][0] == "E":
                G.node[n]["orient"] = self.orient[n-1]
            #G.node[n]["orient"] = self.orient[n]

        # write json formatted data
        d = json_graph.node_link_data(G) # node-link format to serialize
        topoString = json.dumps(d,encoding='UTF-8',default=str)
        return d



        
if __name__ == "__main__":
    pat = 'Ea.0.Ea-0.Ea-0'
    q = contact()
    q.feedMatrixOrientation(pat)
    elol = q.get_Nodes_Edges()
    q.get_networkx(elol)
