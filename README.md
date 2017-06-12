# plv-DecoySearch
extension of ProLego for decoy search from a set of PDB ensemble

## ProLego Topology Evaluation 
Topology Evaluation is mainly divided into two parts. 
1.Topology assignment and
2.Search topology and modules in the Exclusive dataset.

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
  
##Installing and Running
 
 from your shell do following
	 *[Install PIP] sudo apt-get install pip | sudo apt-get install python-pip
	 *[Get Into Directory] cd ./ProLegoDecoySearch
	 *[Get Requisistes] pip install -r requisites.txt
	 *[Get the current working directory] pwd
	 *[Copy the Path]
	 *[Modify the Parameter Setting] Change variable appRoot <pwd>
	 *[All set]
  	*[run the python code]
 
## Example 
    []python pLv_search.py -p 3DP8.pdb -id 3DP8
    [output] ./cs_result/3DP8_plv.json
  
## Help and Bug report 
	ent mail to : taushifkhan@gmail.com

