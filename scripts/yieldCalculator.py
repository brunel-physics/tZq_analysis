#from ROOT import *
import ROOT

import subprocess
import sys


systs = ["__trig__plus","__trig__minus","__jer__plus","__jer__minus","__jes__plus","__jes__minus","__pileup__plus","__pileup__minus","__bTag__plus","__bTag__minus","__met__plus","__met__minus","__pdf__plus","__pdf__minus","__ME_PS__plus","__ME_PS__minus"]

channel =  sys.argv[1]
sample =  sys.argv[2]

channelIndex = -1
if ( channel == "ee" ) : channelIndex = 1
if ( channel == "mumu" ) : channelIndex = 0

infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaDirs/inputs/2016/all/mz5mw50/histofile_"+sample+".root")
#infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaDirs/inputs/2015/emu/mz106/histofile_"+sample+".root")

nom_yield  = 0

tree = infile.Get("Ttree_"+sample)
for event in tree :
   if ( channelIndex == event.Channel ) : nom_yield += event.EvtWeight

print "nominal yield: ", nom_yield

for syst in systs:
    syst_yield = 0
    tree = infile.Get("Ttree_"+sample+syst)
    for event in tree : 
        if ( channelIndex == event.Channel ) : syst_yield += event.EvtWeight
    print "syst yield for ", syst, " : ", syst_yield, " / abs diff : ", syst_yield-nom_yield, ", rel diff : ", (syst_yield-nom_yield)/nom_yield*100.0, "%"


