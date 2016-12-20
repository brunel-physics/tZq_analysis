#A tool to pull plots from the mva inputs and plot gaussians

from ROOT import *
import math
import os 
import sys
import subprocess 

def sortOutLeptons(tree,channel):
    ###Returns two LorentzVectors containing the two z leptons. This will be VERY useful for making all of the plots.
    #Reads the position of the z leptons from variables stored at mvaTree making time, because I'm great and finally got around to doing it.
    zLep1,zLep2 = 0,0
    #Let's try commenting this out and see if everything breaks? Hopefully it won't do...
    #if tree.numElePF2PAT < 3:
    if channel == "ee":
        zLep1 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep1Index],tree.elePF2PATGsfPy[tree.zLep1Index],tree.elePF2PATGsfPz[tree.zLep1Index],tree.elePF2PATGsfE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep2Index],tree.elePF2PATGsfPy[tree.zLep2Index],tree.elePF2PATGsfPz[tree.zLep2Index],tree.elePF2PATGsfE[tree.zLep2Index])
    if channel == "mumu":
        zLep1 = TLorentzVector(tree.muonPF2PATPx[tree.zLep1Index],tree.muonPF2PATPy[tree.zLep1Index],tree.muonPF2PATPz[tree.zLep1Index],tree.muonPF2PATE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.muonPF2PATPx[tree.zLep2Index],tree.muonPF2PATPy[tree.zLep2Index],tree.muonPF2PATPz[tree.zLep2Index],tree.muonPF2PATE[tree.zLep2Index])
    return (zLep1,zLep2)

def main():

  era = "2016"
  channel = "ee"
#  channel = "mumu"

  weighted = False

### Number of Same Sign no Fakes stuff

  zRefMass = 91.1
  zWindow = 5.0
  sameSignDY = 0
  oppSignDY = 0

  infile_DY = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz100mw50/DYJetsToLL_M-50"+channel+"mvaOut.root")
  tree_DY = infile_DY.Get("tree")

  infile_DY_SS = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz100mw50/DYJetsToLL_M-50"+channel+"invLepmvaOut.root")
  tree_DY_SS = infile_DY_SS.Get("tree")

  ## DY Histos

  DY_zMassSameOppHisto = TH1D("DY_zMassOppSignHisto","Z Mass Histo from Opposite Sign events", 300, 0.0, 300.0)
  DY_zMassSameSignHisto = TH1D("DY_zMassSameSignHisto","Z Mass Histo from Same Sign events", 300, 0.0, 300.0)

  for event in range ( tree_DY_SS.GetEntries() ) :
    tree_DY_SS.GetEntry(event)

    weight = 1
    if (weighted) : weight = tree_DY_SS.eventWeight

    (zLep1,zLep2) = sortOutLeptons(tree_DY_SS,channel)
    zMass = (zLep1+zLep2).M()

    if ( zMass < (zRefMass + zWindow) and zMass > (zRefMass - zWindow) ) : sameSignDY += 1*weight

    DY_zMassSameSignHisto.Fill(zMass,weight)

  for event in range ( tree_DY.GetEntries() ) :
    tree_DY.GetEntry(event)

    weight = 1.0
    if (weighted) : weight = tree_DYS.eventWeight

    (zLep1,zLep2) = sortOutLeptons(tree_DY,channel)
    zMass = (zLep1+zLep2).M()

    if ( zMass < (zRefMass + zWindow) and zMass > (zRefMass - zWindow) ) : oppSignDY += 1.0*weight

    DY_zMassOppSignHisto.Fill(zMass,weight)


  subprocess.call("mkdir plots/fakeLeptons/",shell=True)
  subprocess.call("mkdir plots/fakeLeptons/DY/",shell=True)

#  DY_zMassHisto.Fit("gaus")

  DY_zMassSameSignHisto.SaveAs("plots/fakeLeptons/DY/zMass_"+channel+"_SameSign.root")
  DY_zMassOppSignHisto.SaveAs("plots/fakeLeptons/DY/zMass_"+channel+"_OppSign.root")

  infile_DY_SS.Close()

##############

  listOfMCs = {"WW1l1nu2q" : "WW", "WW2l2nu":"WW","ZZ4l":"ZZ","ZZ2l2nu":"ZZ","ZZ2l2q":"ZZ","WZjets":"WZ","WZ2l2q":"WZ","WZ1l1nu2q":"WZ","sChannel":"TsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","tZq":"tZq","tHq":"THQ","ttWlnu":"TTW","ttW2q":"TTW","ttZ2l2nu":"TTZ","ttZ2q":"TTZ","ttbarInclusivePowerheg":"TT","wPlusJets":"Wjets","DYJetsToLL_M-50":"DYToLL_M50","DYJetsToLL_M-10To50":"DYToLL_M10To50"}

### Number of SS events expected from data

  sameSignMC_gen = 0
  sameSignMC_fake =0
  sameSignMC = 0
  sameSignData = 0

  #Loop over all MC samples
  for sample in listOfMCs.keys():
    print "Doing " + sample + ": ",
    sys.stdout.flush()
    infile_SS_MC = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz5mw50/"+sample+channel+"invLepmvaOut.root")
    tree_SS_MC = infile_SS_MC.Get("tree")
    try:
      print str(tree_SS_MC.GetEntriesFast())
      sys.stdout.flush()
    except AttributeError:
      print "\nAttribute Error \n"
      sys.stdout.flush()

    for event in range ( tree_SS_MC.GetEntries() ) :
      tree_SS_MC.GetEntry(event)
      weight = 1.0
      if (weighted) : weight = tree_SS_MC.eventWeight
      sameSignMC += 1.0*weight

      lep1ID, lep2ID, = False,False

      if (channel == "ee") :
        if ( abs(tree_SS_MC.genElePF2PATMotherId[tree_SS_MC.zLep1Index]) == 23 or abs(tree_SS_MC.genElePF2PATMotherId[tree_SS_MC.zLep1Index]) == 24 ) : lep1ID = True
        if ( abs(tree_SS_MC.genElePF2PATMotherId[tree_SS_MC.zLep2Index]) == 23 or abs(tree_SS_MC.genElePF2PATMotherId[tree_SS_MC.zLep2Index]) == 24 ) : lep2ID = True

        if ( lep1ID and lep2ID ) : sameSignMC_gen += 1.0*weight
        else : sameSignMC_fake += 1.0*weight
      elif (channel == "mumu") :
        if( abs(tree_SS_MC.genMuonPF2PATMotherId[tree_SS_MC.zLep1Index]) == 23 or abs(tree_SS_MC.genMuonPF2PATMotherId[tree_SS_MC.zLep1Index]) == 24 ) : lep1ID = True
        if( abs(tree_SS_MC.genMuonPF2PATMotherId[tree_SS_MC.zLep2Index]) == 23 or abs(tree_SS_MC.genMuonPF2PATMotherId[tree_SS_MC.zLep2Index]) == 24 ) : lep2ID = True

        if( lep1ID and lep2ID ): sameSignMC_gen += 1.0*weight
        else : sameSignMC_fake += 1.0*weight

    infile_SS_MC.Close()

  #Loop over data samples
  print "Doing " + channel + ": ",
  sys.stdout.flush()
  infile_SS_data = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz5mw50/"+channel+"Run2016"+channel+"invLepmvaOut.root")
  tree_SS_data = infile_SS_data.Get("tree")
  try:
    print str(tree_SS_data.GetEntriesFast())
    sys.stdout.flush()
  except AttributeError:
    print "\nAttribute Error \n"
    sys.stdout.flush()

  for event in range ( tree_SS_data.GetEntries() ) :
    tree_SS_data.GetEntry(event)
    weight = 1.0
    if (weighted) : weight = tree_SS_data.eventWeight

    sameSignData += 1.0*weight

  infile_SS_data.Close()

### Number of OS events expected from data

  oppSignMC_gen = 0
  oppSignMC_fake = 0
  oppSignMC = 0
  oppSignData = 0

  #Loop over all MC samples
  for sample in listOfMCs.keys():
    print "Doing " + sample + ": ",
    sys.stdout.flush()

    infile_OS_MC = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz5mw50/"+sample+channel+"mvaOut.root")
    tree_OS_MC = infile_OS_MC.Get("tree")
    try:
      print str(tree_OS_MC.GetEntriesFast())
      sys.stdout.flush()
    except AttributeError:
      print "\nAttribute Error \n"
      sys.stdout.flush()

    for event in range ( tree_OS_MC.GetEntries() ) :
      tree_OS_MC.GetEntry(event)
      weight = 1.0
      if (weighted) : weight = tree_OS_MC.eventWeight
      oppSignMC += 1.0*weight

      lep1ID, lep2ID, = False,False

      if (channel == "ee") :
        if( abs(tree_OS_MC.genElePF2PATMotherId[tree_OS_MC.zLep1Index]) == 23 or abs(tree_OS_MC.genElePF2PATMotherId[tree_OS_MC.zLep1Index]) == 24 ): lep1ID = True
        if( abs(tree_OS_MC.genElePF2PATMotherId[tree_OS_MC.zLep2Index]) == 23 or abs(tree_OS_MC.genElePF2PATMotherId[tree_OS_MC.zLep2Index]) == 24 ): lep2ID = True

        if( lep1ID and lep2ID ) : oppSignMC_gen += 1.0*weight
        else : oppSignMC_fake += 1.0*weight

      elif (channel == "mumu") :
        if( abs(tree_OS_MC.genMuonPF2PATMotherId[tree_OS_MC.zLep1Index]) == 23 or abs(tree_OS_MC.genMuonPF2PATMotherId[tree_OS_MC.zLep1Index]) == 24): lep1ID = True
        if( abs(tree_OS_MC.genMuonPF2PATMotherId[tree_OS_MC.zLep2Index]) == 23 or abs(tree_OS_MC.genMuonPF2PATMotherId[tree_OS_MC.zLep2Index]) == 24): lep2ID = True

        if( lep1ID and lep2ID ) : oppSignMC_gen += 1.0*weight
        else : oppSignMC_fake += 1.0*weight

    infile_OS_MC.Close()

  #Loop over data samples
  sys.stdout.flush()
  infile_OS_data = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz5mw50/"+channel+"Run2016"+channel+"mvaOut.root")
  tree_OS_data = infile_OS_data.Get("tree")
  try:
    print str(tree_OS_data.GetEntriesFast())
    sys.stdout.flush()
  except AttributeError:
    print "\nAttribute Error \n"
    sys.stdout.flush()

  for event in range ( tree_OS_data.GetEntries() ) :
    tree_OS_data.GetEntry(event)
    weight = 1.0
    if (weighted) : weight = tree_OS_data.eventWeight
    oppSignData += 1.0*weight

  infile_OS_data.Close()

##############
#Ratio between OS and SS events

  print "gen OS/SS: ", oppSignMC_gen, "/", sameSignMC_gen, " : " , oppSignMC_gen/(sameSignMC_gen + 1.0e-6)
  print "fake OS/SS: ", oppSignMC_fake, "/", sameSignMC_fake, " : " , oppSignMC_fake/(sameSignMC_fake + 1.0e-6)

##############
#Number of expected Same sign events with no fakes - DY mis-id stuff
  eff = sameSignDY/(sameSignDY + oppSignDY)

  print "sameSignDY:oppSignDY = ", sameSignDY, " : " , oppSignDY
  print "Efficiency coefficient for calculating the number of expected same sign events with no fakes = ", eff

if __name__ == "__main__":
    main()
