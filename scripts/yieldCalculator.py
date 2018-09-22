#from ROOT import *
import ROOT

import math
import subprocess
import sys

systs = ["__trig__plus","__trig__minus","__jer__plus","__jer__minus","__jes__plus","__jes__minus","__pileup__plus","__pileup__minus","__bTag__plus","__bTag__minus","__met__plus","__met__minus","__pdf__plus","__pdf__minus","__ME__plus","__ME__minus"]

channel =  sys.argv[1]
sample =  sys.argv[2]

channelIndex = -1
if ( channel == "emu" ) : channelIndex = 2
if ( channel == "ee" ) : channelIndex = 1
if ( channel == "mumu" ) : channelIndex = 0

#infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaDirs/inputs/2016/all/mz20mw20/histofile_"+sample+".root")
infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaDirs/inputs/2016/all/mz20mw20_zPlusJets/histofile_"+sample+".root")
#infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaDirs/inputs/2016/all/mz20mw20_zPlusJets_oldCR/histofile_"+sample+".root")
#infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaDirs/inputs/2016/emu/histofile_"+sample+".root")

doSysts = 0
nom_yield  = 0
nEvents = 0


tree = infile.Get("Ttree_"+sample)
for event in tree :
   if ( channelIndex == event.Channel ) : nom_yield += event.EvtWeight
   if ( channelIndex == event.Channel ) : nEvents += 1

print "nominal yield: ", nom_yield
print "math.sqrt(nEvents): " , math.sqrt(nEvents)
print "stat error % : " , math.sqrt(nEvents)/nEvents *100
print "stat error: " , math.sqrt(nEvents)/nEvents * nom_yield

if ( doSysts == 1 ) :
   for syst in systs:
      syst_yield = 0
      tree = infile.Get("Ttree_"+sample+syst)
      for event in tree :
         if ( channelIndex == event.Channel ) : syst_yield += event.EvtWeight
      print "syst yield for ", syst, " : ", syst_yield, " / abs diff : ", syst_yield-nom_yield, ", rel diff : ", (syst_yield-nom_yield)/(nom_yield)*100.0, "%"


