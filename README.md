# tqZ2015Code
Code I use for tqZ analysis for Run 2

Currently 2015RunC and 2015RunD and MC is set up for 25ns.

Before running any scripts or executables, execute "source setup.sh" and run all from the root directory.

To setup file lists, execute: "bash scripts/setupFileLists.sh"

To create MC pileup, execute: "root -l scripts/createPileUpMC.C"

To create data pileup files, in cmssw release execute:
	pileupCalc.py -i MyAnalysisJSON.txt --inputLumiJSON pileup_latest.txt  --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50  MyDataPileupHistogram.root
where MyAnalysisJSON.txt is the JSON used in creating nTuples of data, and pileup_latest.txt is the latest pileup JSON ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt").
For scale up/down, just vary the inelastic cross-section by the perscribed uncertanity (currently +/- 5%).
See https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#2015_Pileup_JSON_Files for further info.


Instructions as how to update scripts/createPileUpMC.C is contained within the comments of said file.

Systematics Status:

Trigger: To-do
JEC and JER: 2015 Spring 15
bTagging: To-do
Pileup: Up to date (MC & data 2015 76X)
PDF: Up to date.
