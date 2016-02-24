# tqZ2015Code
Code I use for tqZ analysis for Run 2

Currently 2015RunC and 2015RunD and MC is set up for 25ns.

Before running any scripts or executables, execute "source setup.sh" and run all from the root directory.

To setup file lists, execute: "bash scripts/setupFileLists.sh"

To create MC pileup, execute: "root -l scripts/createPileUpMC.C"
Instructions as how to update scripts/createPileUpMC.C is contained within the comments of said file.

Systematics Status:

Trigger: To-do
JEC and JER: 2015 Spring 15
bTagging: To-do
Pileup: MC 2015, data 2012
PDF: Up to date.
