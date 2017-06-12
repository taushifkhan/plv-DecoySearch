#!/usr/bin/python
helpDoc = """
PDB header parser. get formatted information from a PDB 
file of annotation related to a PDB Id.
"""
import urllib
import os

class pdbHeader():
	def __init__(self,pdbHeader):
		self.header={'HEADER':'','TITLE':'','COMPND':'','SOURCE':'','EXPDTA':'','NUMMDL':''}
		self.headInfo = pdbHeader

	def _decompose_cmpd__(self,compnd):
		molIds = compnd.split("MOL_ID")
		cmpndDict = {}
		
		for m in range(1,len(molIds),2):
			mids = molIds[m].split(";")
			mol  = mids[0].split(":")[-1].strip()
			cmpndDict[mol] = {}
			for k in mids[1:]:
				kmol = k.split(":")
				cmpndDict[mol][kmol[0].strip()] = kmol[1].strip()
		return cmpndDict

	def getData(self):
		for l in self.headInfo:
			nameId = l[:6].strip()
			if nameId in self.header.keys():
				self.header[nameId] += l[10:80].strip()
			else:
				pass

		cmpndDict = self._decompose_cmpd__(self.header['COMPND'])
		sourceDict =  self._decompose_cmpd__(self.header['SOURCE'])

		if cmpndDict.keys():
			self.header['COMPND'] = {}
			self.header['COMPND'] = cmpndDict
			self.header['SOURCE'] = {}
			self.header['SOURCE'] = sourceDict

		headerTitle = self.header['HEADER'][:40].strip().lower()
		headPdbId = self.header['HEADER'][52:56]
		self.header['HEADER'] = {'Name':headerTitle,'Id':headPdbId}

		return self.header

def getPDBheaderFile(pid):
	"""
	@param: 
	pid -> 4 letter pdb Id - Alpha numberic
	@return
	headerFile -> open file in an array and file path
	__docString__ = downloads only header information from rcsb
	"""
	cwd = os.getcwd()
	headFile = urllib.URLopener()
	urlHeader = "https://files.rcsb.org/header/%s.pdb"%(pid.upper())
	headsave = cwd+'/%s_header.pdb'%pid
	try:
		headFile.retrieve(urlHeader, headsave)
		headRead = open(headsave,'r').readlines()
		print "Header file:", headsave
	except Exception as e:
		raise e

	return headRead

def main(pdbId):
	headRead = getPDBheaderFile(pdbId)
	pH = pdbHeader(headRead)
	hdict = pH.getData()
	import ipdb; ipdb.set_trace();


if __name__ == '__main__':
	import sys
	main(sys.argv[1])