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


  weighted = True

  era = sys.argv[1]

  mzCut = sys.argv[2]
  mzStr = mzCut.split(".")[0]

  mwCut = sys.argv[3]
  mwStr = mwCut.split(".")[0]

### Number of Same Sign no Fakes stuff

  zRefMass = 91.1
  zWindow = 5.0

  sameSignDY_ee = 0
  oppSignDY_ee = 0

  sameSignDY_mumu = 0
  oppSignDY_mumu = 0


  infile_DY_ee = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/DYJetsToLL_M-50eemvaOut.root")
  infile_DY_mumu = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/DYJetsToLL_M-50mumumvaOut.root")
  tree_DY_ee = infile_DY_ee.Get("tree")
  tree_DY_mumu = infile_DY_mumu.Get("tree")

  infile_DY_SS_ee = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/DYJetsToLL_M-50eeinvLepmvaOut.root")
  infile_DY_SS_mumu = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/DYJetsToLL_M-50mumuinvLepmvaOut.root")
  tree_DY_SS_ee = infile_DY_SS_ee.Get("tree")
  tree_DY_SS_mumu = infile_DY_SS_mumu.Get("tree")

  ## DY Histos

  DY_zMassSameOppHisto_ee = TH1D("DY_zMassOppSignHisto_ee","Z Mass Histo (ee) from Opposite Sign events", 300, 0.0, 300.0)
  DY_zMassSameSignHisto_ee = TH1D("DY_zMassSameSignHisto_ee","Z Mass Histo from (ee) Same Sign events", 300, 0.0, 300.0)
  DY_zMassSameOppHisto_mumu = TH1D("DY_zMassOppSignHisto_mumu","Z Mass Histo from (mumu) Opposite Sign events", 300, 0.0, 300.0)
  DY_zMassSameSignHisto_mumu = TH1D("DY_zMassSameSignHisto_mumu","Z Mass Histo from (mumu) Same Sign events", 300, 0.0, 300.0)

  for event in range ( tree_DY_SS_ee.GetEntries() ) :
    tree_DY_SS_ee.GetEntry(event)

    weight = 1
    if (weighted) : weight = tree_DY_SS_ee.eventWeight

    (zLep1,zLep2) = sortOutLeptons(tree_DY_SS_ee,"ee")
    zMass = (zLep1+zLep2).M()

    if ( zMass < (zRefMass + zWindow) and zMass > (zRefMass - zWindow) ) : sameSignDY_ee += 1*weight

    DY_zMassSameSignHisto_ee.Fill(zMass,weight)

  for event in range ( tree_DY_SS_mumu.GetEntries() ) :
    tree_DY_SS_mumu.GetEntry(event)

    weight = 1
    if (weighted) : weight = tree_DY_SS_mumu.eventWeight

    (zLep1,zLep2) = sortOutLeptons(tree_DY_SS_mumu,"mumu")
    zMass = (zLep1+zLep2).M()

    if ( zMass < (zRefMass + zWindow) and zMass > (zRefMass - zWindow) ) : sameSignDY_mumu += 1*weight

    DY_zMassSameSignHisto_mumu.Fill(zMass,weight)

  for event in range ( tree_DY_ee.GetEntries() ) :
    tree_DY_ee.GetEntry(event)

    weight = 1.0
    if (weighted) : weight = tree_DY_ee.eventWeight

    (zLep1,zLep2) = sortOutLeptons(tree_DY_ee,"ee")
    zMass = (zLep1+zLep2).M()

    if ( zMass < (zRefMass + zWindow) and zMass > (zRefMass - zWindow) ) : oppSignDY_ee += 1.0*weight

    DY_zMassOppSignHisto_ee.Fill(zMass,weight)

  for event in range ( tree_DY_mumu.GetEntries() ) :
    tree_DY_mumu.GetEntry(event)

    weight = 1.0
    if (weighted) : weight = tree_DY_mumu.eventWeight

    (zLep1,zLep2) = sortOutLeptons(tree_DY_mumu,"mumu")
    zMass = (zLep1+zLep2).M()

    if ( zMass < (zRefMass + zWindow) and zMass > (zRefMass - zWindow) ) : oppSignDY_mumu += 1.0*weight

    DY_zMassOppSignHisto_mumu.Fill(zMass,weight)

  subprocess.call("mkdir plots/fakeLeptons/",shell=True)
  subprocess.call("mkdir plots/fakeLeptons/DY/",shell=True)

#  DY_zMassHisto.Fit("gaus")

  DY_zMassSameSignHisto_ee.SaveAs("plots/fakeLeptons/DY/zMass_ee_SameSign.root")
  DY_zMassOppSignHisto_ee.SaveAs("plots/fakeLeptons/DY/zMass_ee_OppSign.root")
  DY_zMassSameSignHisto_mumu.SaveAs("plots/fakeLeptons/DY/zMass_mumu_SameSign.root")
  DY_zMassOppSignHisto_mumu.SaveAs("plots/fakeLeptons/DY/zMass_mumu_OppSign.root")

  infile_DY_SS_ee.Close()
  infile_DY_SS_mumu.Close()
  infile_DY_ee.Close()
  infile_DY_mumu.Close()

##############
## With QCD MC
#  listOfMCs = {"WW1l1nu2q":"WW", "WW2l2nu":"WW","ZZ4l":"ZZ","ZZ2l2nu":"ZZ","ZZ2l2q":"ZZ","WZjets":"WZ","WZ2l2q":"WZ","WZ1l1nu2q":"WZ","sChannel":"TsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","tZq":"tZq","tHq":"THQ","ttWlnu":"TTW","ttW2q":"TTW","ttZ2l2nu":"TTZ","ttZ2q":"TTZ","ttbarDilepton_aMCatNLO":"TT","tWZ":"TWZ","wPlusJets":"Wjets","DYJetsToLL_M-50":"DYToLL_M50","DYJetsToLL_M-10To50":"DYToLL_M10To50"} #,"QCD_EMEnriched_Pt-20to30":"QCD_EM_20to30","QCD_EMEnriched_Pt-30to50":"QCD_EM_30to50","QCD_EMEnriched_Pt-50to80":"QCD_EM_80to120","QCD_EMEnriched_Pt-120to170":"QCD_EM_120to170","QCD_EMEnriched_Pt-170to300":"QCD_EM_170to300","QCD_EMEnriched_Pt-300toInf":"QCD_EM_300toInf","QCD_MuEnriched_Pt-15to20":"QCD_Mu_15to20","QCD_MuEnriched_Pt-20to30":"QCD_Mu_20to30","QCD_MuEnriched_Pt-30to50":"QCD_Mu_30to50","QCD_MuEnriched_Pt-50to80":"QCD_Mu_50to80","QCD_MuEnriched_Pt-120to170":"QCD_Mu_120to170","QCD_MuEnriched_Pt-170to300":"QCD_Mu_170to300","QCD_MuEnriched_Pt-300to470":"QCD_Mu_300to470","QCD_MuEnriched_Pt-470to600":"QCD_Mu_470to600","QCD_MuEnriched_Pt-600to800":"QCD_Mu_600to800","QCD_MuEnriched_Pt-800to1000":"QCD_Mu_800to1000","QCD_MuEnriched_Pt-1000toInf":"QCD_Mu_1000toInf"}
## Without QCD MC
#  listOfMCs = {"WW1l1nu2q":"WW", "WW2l2nu":"WW","ZZ4l":"ZZ","ZZ2l2nu":"ZZ","ZZ2l2q":"ZZ","WZjets":"WZ","WZ2l2q":"WZ","WZ1l1nu2q":"WZ","sChannel":"TsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","tZq":"tZq","tHq":"THQ","ttWlnu":"TTW","ttW2q":"TTW","ttZ2l2nu":"TTZ","ttZ2q":"TTZ","ttbarDilepton_aMCatNLO":"TT","tWZ":"TWZ","wPlusJets":"Wjets","DYJetsToLL_M-50":"DYToLL_M50","DYJetsToLL_M-10To50":"DYToLL_M10To50"}

#  listOfMCs = {"ttbarDilepton_aMCatNLO":"TT"}
  if era == "2016":
    listOfMCs = {"tZq":"tZq","tHq":"THQ","ttWlnu":"TTW","ttZ2l2nu":"TTZ","wPlusJets":"Wjets","WZjets":"WZ"}
  else:
    listOfMCs = {"tZq":"tZq","tHq":"THQ","ttWTolnu":"TTW","ttZToll":"TTZ","wPlusJets":"Wjets","WZ_3lnu":"WZ"}

### Number of SS events expected from data

  genGen_ee = 0.
  fakeGen_ee = 0.
  genGen_mumu = 0.
  fakeGen_mumu = 0.
  genGen_emu = 0.
  fakeGen_emu = 0.

  genGen_inv_ee = 0.
  fakeGen_inv_ee = 0.
  genGen_inv_mumu = 0.
  fakeGen_inv_mumu = 0.
  genGen_inv_emu = 0.
  fakeGen_inv_emu = 0.

  #Loop over all MC samples
  for sample in listOfMCs.keys():
    print "Doing " + sample + ": ",
    sys.stdout.flush()
    infile_MC_ee = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/"+sample+"eemvaOut.root")
    infile_MC_mumu = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/"+sample+"mumumvaOut.root")
    infile_MC_inv_ee = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/"+sample+"eeinvLepmvaOut.root")
    infile_MC_inv_mumu = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/"+sample+"mumuinvLepmvaOut.root")
    infile_MC_emu = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/"+sample+"emumvaOut.root")
    infile_MC_inv_emu = TFile.Open("/scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz"+mzStr+"mw"+mwStr+"/"+sample+"emuinvLepmvaOut.root")

    tree_MC_ee = infile_MC_ee.Get("tree")
    try:
      print str(tree_MC_ee.GetEntriesFast())
      sys.stdout.flush()
    except AttributeError:
      print "\nAttribute Error \n"
      sys.stdout.flush()

    for event in range ( tree_MC_ee.GetEntries() ) :
      tree_MC_ee.GetEntry(event)
      weight = 1
      if (weighted) : weight = tree_MC_ee.eventWeight

      lep1_ee = 0
      lep2_ee = 0
      if (tree_MC_ee.genElePF2PATPromptDecayed[tree_MC_ee.zLep1Index] == 1 or tree_MC_ee.genElePF2PATPromptFinalState[tree_MC_ee.zLep1Index] == 1 ) : lep1_ee = 1
      if (tree_MC_ee.genElePF2PATPromptDecayed[tree_MC_ee.zLep2Index] == 1 or tree_MC_ee.genElePF2PATPromptFinalState[tree_MC_ee.zLep2Index] == 1 ) : lep2_ee = 1

      if ( lep1_ee == 1 and lep2_ee == 1 ) : genGen_ee += 1.0 * weight
      else: fakeGen_ee += 1.0 * weight

    infile_MC_ee.Close()

    tree_MC_mumu = infile_MC_mumu.Get("tree")
    try:
      print str(tree_MC_mumu.GetEntriesFast())
      sys.stdout.flush()
    except AttributeError:
      print "\nAttribute Error \n"
      sys.stdout.flush()

    for event in range ( tree_MC_mumu.GetEntries() ) :
      tree_MC_mumu.GetEntry(event)
      weight = 1
      if (weighted) : weight = tree_MC_mumu.eventWeight

      lep1_mumu = 0
      lep2_mumu = 0

      if (tree_MC_mumu.genMuonPF2PATPromptDecayed[tree_MC_mumu.zLep1Index] == 1 or tree_MC_mumu.genMuonPF2PATPromptFinalState[tree_MC_mumu.zLep1Index] == 1 ) : lep1_mumu = 1
      if (tree_MC_mumu.genMuonPF2PATPromptDecayed[tree_MC_mumu.zLep2Index] == 1 or tree_MC_mumu.genMuonPF2PATPromptFinalState[tree_MC_mumu.zLep2Index] == 1 ) : lep2_mumu = 1

      if ( lep1_mumu == 1 and lep2_mumu == 1 ) : genGen_mumu += 1.0 * weight
      else: fakeGen_mumu += 1.0 * weight

    infile_MC_mumu.Close()

    tree_MC_emu = infile_MC_emu.Get("tree")
    try:
      print str(tree_MC_emu.GetEntriesFast())
      sys.stdout.flush()
    except AttributeError:
      print "\nAttribute Error \n"
      sys.stdout.flush()

    for event in range ( tree_MC_emu.GetEntries() ) :
      tree_MC_emu.GetEntry(event)
      weight = 1
      if (weighted) : weight = tree_MC_emu.eventWeight

      lep1_emu = 0
      lep2_emu = 0

      if (tree_MC_emu.genElePF2PATPromptDecayed[tree_MC_emu.zLep1Index] == 1 or tree_MC_emu.genElePF2PATPromptFinalState[tree_MC_emu.zLep1Index] == 1 ) : lep1_emu = 1
      if (tree_MC_emu.genMuonPF2PATPromptDecayed[tree_MC_emu.zLep2Index] == 1 or tree_MC_emu.genMuonPF2PATPromptFinalState[tree_MC_emu.zLep2Index] == 1 ) : lep2_emu = 1

      if ( lep1_emu == 1 and lep2_emu == 1 ) : genGen_emu += 1.0 * weight
      else: fakeGen_emu += 1.0 * weight

    infile_MC_emu.Close()

    tree_MC_inv_ee = infile_MC_inv_ee.Get("tree")
    try:
      print str(tree_MC_inv_ee.GetEntriesFast())
      sys.stdout.flush()
    except AttributeError:
      print "\nAttribute Error \n"
      sys.stdout.flush()

    for event in range ( tree_MC_inv_ee.GetEntries() ) :
      tree_MC_inv_ee.GetEntry(event)
      weight = 1
      if (weighted) : weight = tree_MC_inv_ee.eventWeight

      lep1_inv_ee = 0
      lep2_inv_ee = 0

      if (tree_MC_inv_ee.genElePF2PATPromptDecayed[tree_MC_inv_ee.zLep1Index] == 1 or tree_MC_inv_ee.genElePF2PATPromptFinalState[tree_MC_inv_ee.zLep1Index] == 1 ) : lep1_inv_ee = 1
      if (tree_MC_inv_ee.genElePF2PATPromptDecayed[tree_MC_inv_ee.zLep2Index] == 1 or tree_MC_inv_ee.genElePF2PATPromptFinalState[tree_MC_inv_ee.zLep2Index] == 1 ) : lep2_inv_ee = 1

      if ( lep1_inv_ee == 1 and lep2_inv_ee == 1 ) : genGen_inv_ee += 1.0 * weight
      else: fakeGen_inv_ee += 1.0 * weight

    infile_MC_inv_ee.Close()

    tree_MC_inv_mumu = infile_MC_inv_mumu.Get("tree")
    try:
      print str(tree_MC_inv_mumu.GetEntriesFast())
      sys.stdout.flush()
    except AttributeError:
      print "\nAttribute Error \n"
      sys.stdout.flush()

    for event in range ( tree_MC_inv_mumu.GetEntries() ) :
      tree_MC_inv_mumu.GetEntry(event)
      weight = 1
      if (weighted) : weight = tree_MC_inv_mumu.eventWeight

      lep1_inv_mumu = 0
      lep2_inv_mumu = 0

      if (tree_MC_inv_mumu.genMuonPF2PATPromptDecayed[tree_MC_inv_mumu.zLep1Index] == 1 or tree_MC_inv_mumu.genMuonPF2PATPromptFinalState[tree_MC_inv_mumu.zLep1Index] == 1 ) : lep1_inv_mumu = 1
      if (tree_MC_inv_mumu.genMuonPF2PATPromptDecayed[tree_MC_inv_mumu.zLep2Index] == 1 or tree_MC_inv_mumu.genMuonPF2PATPromptFinalState[tree_MC_inv_mumu.zLep2Index] == 1 ) : lep2_inv_mumu = 1

      if ( lep1_inv_mumu == 1 and lep2_inv_mumu == 1 ) : genGen_inv_mumu += 1.0 * weight
      else: fakeGen_inv_mumu += 1.0 * weight

    infile_MC_inv_mumu.Close()

    tree_MC_inv_emu = infile_MC_inv_emu.Get("tree")
    try:
      print str(tree_MC_inv_emu.GetEntriesFast())
      sys.stdout.flush()
    except AttributeError:
      print "\nAttribute Error \n"
      sys.stdout.flush()

    for event in range ( tree_MC_inv_emu.GetEntries() ) :
      tree_MC_inv_emu.GetEntry(event)
      weight = 1
      if (weighted) : weight = tree_MC_inv_emu.eventWeight

      lep1_inv_emu = 0
      lep2_inv_emu = 0

      if (tree_MC_inv_emu.genElePF2PATPromptDecayed[tree_MC_inv_emu.zLep1Index] == 1 or tree_MC_inv_emu.genElePF2PATPromptFinalState[tree_MC_inv_emu.zLep1Index] == 1 ) : lep1_inv_emu = 1
      if (tree_MC_inv_emu.genMuonPF2PATPromptDecayed[tree_MC_inv_emu.zLep2Index] == 1 or tree_MC_inv_emu.genMuonPF2PATPromptFinalState[tree_MC_inv_emu.zLep2Index] == 1 ) : lep2_inv_emu = 1

      if ( lep1_inv_emu == 1 and lep2_inv_emu == 1 ) : genGen_inv_emu += 1.0 * weight
      else: fakeGen_inv_emu += 1.0 * weight

    infile_MC_inv_emu.Close()

  print "number of prompt (final state or decayed) OS electron gen particles: ", genGen_ee
  print "number of non-prompt (final state or decayed) OS electron gen particles: ", fakeGen_ee

  print "number of prompt (final state or decayed) SS electron gen particles: ", genGen_inv_ee
  print "number of non-prompt (final state or decayed) SS electron gen particles: ", fakeGen_inv_ee

  print "number of prompt (final state or decayed) OS muon gen particles: ", genGen_mumu
  print "number of non-prompt (final state or decayed) OS muon gen particles: ", fakeGen_mumu

  print "number of prompt (final state or decayed) SS muon gen particles: ", genGen_inv_mumu
  print "number of non-prompt (final state or decayed) SS muon gen particles: ", fakeGen_inv_mumu

  print "number of prompt (final state or decayed) OS emu gen particles: ", genGen_emu
  print "number of non-prompt (final state or decayed) OS emu gen particles: ", fakeGen_emu

  print "number of prompt (final state or decayed) SS emu gen particles: ", genGen_inv_emu
  print "number of non-prompt (final state or decayed) SS emu gen particles: ", fakeGen_inv_emu

  print "number of OS/SS prompt electron gen-matched particles: ", abs( genGen_ee/(genGen_inv_ee + 1.0e-06) )
  print "number of OS/SS prompt muon gen-matched particles: ", abs( genGen_mumu/(genGen_inv_mumu + 1.0e-06) )
  print "number of OS/SS prompt emu gen-matched particles: ", abs( genGen_emu/(genGen_inv_emu + 1.0e-06) )

  print "number of OS/SS non-prompt electron gen-matched particles: ", abs( fakeGen_ee/(fakeGen_inv_ee + 1.0e-06) )
  print "number of OS/SS non-prompt muon gen-matched particles: ", abs( fakeGen_mumu/(fakeGen_inv_mumu + 1.0e-06) )
  print "number of OS/SS non-prompt emu gen-matched particles: ", abs( fakeGen_emu/(fakeGen_inv_emu + 1.0e-06) )

##############
#Number of expected Same sign events with no fakes - DY mis-id stuff
  eff_ee = sameSignDY_ee/(sameSignDY_ee + oppSignDY_ee)
  eff_mumu = sameSignDY_mumu/(sameSignDY_mumu + oppSignDY_mumu)

  print "ee sameSignDY:oppSignDY = ", sameSignDY_ee, " : " , oppSignDY_ee
  print "mumu sameSignDY:oppSignDY = ", sameSignDY_mumu, " : " , oppSignDY_mumu
  print "Efficiency coefficient for calculating the number of expected same sign events with no fakes ee/mumu = ", eff_ee, "/", eff_mumu

if __name__ == "__main__":
    main()
