# plv-DecoySearch
extension of ProLego for decoy search from a set of PDB ensemble

## ProLego Topology Evaluation 
Topology Evaluation is mainly divided into two parts. 
 1. Topology assignment and
 2. Search topology and modules in the Exclusive dataset.

Both (1 & 2) objective is done by the code: 

```
ProLegoDecoySearch$ python pLv_search.py -h

usage: pLv_search.py [-h] [-p PDBFL] [-c CHAIN] [-id IDNAME]
 optional arguments:
    -h, --help  show this help message and exit
    -p PDBFL    pdbFile coordinates
    -c CHAIN    protein chain [Optional]
    -id IDNAME  pdbId to save results
```
  
The output are in a JSON format (pid_plv.json) example 3DP8_plv.json

## Installing and Running

from your shell do following
 * [Install PIP] sudo apt-get install pip | sudo apt-get install python-pip
 * [Get Into Directory] cd ./plv-DecoySearch
 * [Get Requisistes] pip install -r requisites.txt
 * [Get the current working directory] pwd
 * [Copy the Path]
 * [Modify the Parameter Setting] Change variable appRoot <pwd>
 * [All set]
 * [run the python code]

## Example 
    []python pLv_search.py -p 3DP8.pdb -id 3DP8
    [output] ./cs_result/3DP8_plv.json
 
## About Output
The json file has chain wise information of all protein chains in the PDB file. The architecture of resulant Json file is as following,
* pdbId_Chain1
  * pdbInfo {pdbId, chain}
  * csModules {modules:{sseGroups:toplogyClass}, topology: {contactString, topology status}}
  * pattern {Seq,sseDetail, ssestring, contactLog, residueDictionary, cs, csFormattedString}
* pdbId_Chain2 ...
## Help and Bug report 
mail to : taushifkhan@gmail.com

## Licence
The MIT License (MIT) 
Copyright (c) 2016 Taushif Khan


Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: 

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
