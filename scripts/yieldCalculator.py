#from ROOT import *
import ROOT

import subprocess
import sys


systs = ["","__trig__plus","__trig__minus","__jer__plus","__jer__minus","__jes__plus","__jes__minus","__pileup__plus","__pileup__minus","__bTag__plus","__bTag__minus","__met__plus","__met__minus","__pdf__plus","__pdf__minus","__ME_PS__plus","__ME_PS__minus"]

channel =  sys.argv[1]
sample =  sys.argv[2]

channelIndex = -1
if ( channel == "ee" ) : channelIndex = 1
if ( channel == "mumu" ) : channelIndex = 0

infile = ROOT.TFile.Open("mvaDirs/inputs/2015/all/mz10mw20/histofile_"+sample+".root")

nom_yield  = 0
syst_yield = 0

for syst in systs:
    tree = infile.Get("Ttree_"+sample+syst)
    for event in tree : 
        if ( channelIndex == event.Channel ) : 
            if (syst == "") : nom_yield += event.EvtWeight
            else : syst_yield += event.EvtWeight

print "nominal yield: ", nom_yield
print "syst yield: ", syst_yield


