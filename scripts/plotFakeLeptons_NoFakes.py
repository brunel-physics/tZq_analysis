#A tool to pull plots from the mva inputs and plot gaussians

import ROOT
import subprocess
import sys

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

  weighted = True

  infile_DY = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz100mw50/DYJetsToLL_M-50"+channel+"+invLepmvaOut.root")
  tree_DY = infile_DY.Get("tree")

  infile_DY2 = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/skims/"+era+"/mz100mw50/histofile_DYToLL_M10To50"+channel+"invLepmvaOut.root")
  tree_DY2 = infile_DY2.Get("tree")

  ## DY Histos

  DY_zMassHisto = ROOT.TH1D("DY_zMassHisto","Z Mass Histo", 300, 0.0, 300.0)

  for event in range ( tree_DY.GetEntries() ) :
    tree_DY.GetEntry(event)

    weight = 1
    if (weighted) : weight = tree_DY.eventWeight

    (zLep1,zLep2) = sortOutLeptons(tree,channel)
    zMass = (zLep1+zLep2).M()
    DY_zMassHisto.Fill(zMass,weight)

  for event in range ( tree_DY2.GetEntries() ) :
    tree_DY2.GetEntry(event)

    weight = 1
    if (weighted) : weight = tree_DY2.eventWeight

    (zLep1,zLep2) = sortOutLeptons(tree,channel)
    zMass = (zLep1+zLep2).M()    
    DY_zMassHisto.Fill(zMass,weight)

##############

  subprocess.call("mkdir plots/fakeLeptons/",shell=True)
  subprocess.call("mkdir plots/fakeLeptons/DY/",shell=True)

  DY_zMassHisto.Fit("gaus")

  DY_zMassHisto.SaveAs("plots/fakeLeptons/DY/zMass.root")

if __name__ == "__main__":
    main()
