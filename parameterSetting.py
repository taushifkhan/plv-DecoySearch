# -*- coding: utf-8 -*-
# @Author: taushif
# @Date:   2017-03-28 08:32:18
# @Last Modified by:   taushif
# @Last Modified time: 2017-03-28 16:47:42
#!/usr/local/bin/
"""
storing global parameters
"""
import os
appRoot = '.'
csDir = appRoot + '/cs_modules/'
#DB path
DBPathTop = appRoot+ '/DB/plVdb_topoStatus.pkl'

#InputDirectory
pdbDir = ''

# Log directories
dir_pdb         = appRoot+ '/PDB_files/'
jsonDir         = appRoot+ 'cs_result/'
logDir          = appRoot+ '/log/'

# creating One time folders

if not os.path.exists(logDir):
    print "creating log directory",logDir
    os.makedirs(logDir)

if not os.path.exists(jsonDir):
    print "Creating Dir to save results",jsonDir
    os.makedirs(jsonDir)

if not os.path.exists(dir_pdb):
    print "Creating Dir to save results",dir_pdb
    os.makedirs(dir_pdb)
