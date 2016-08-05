#A tool to pull plots from the mva inputs and plot gaussians

import ROOT
import subprocess
import sys

def main():

  era = "2015"

  massCut = True

  infile_tZq = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_tZq.root")
  infile_TT = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_TT.root")
  infile_DY = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_DYToLL_M50.root")

  infile_DY2 = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_DYToLL_M10To50.root")
  infile_TbartChan = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_TbartChan.root")
  infile_TbartW = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_TbartW.root")
  infile_THQ = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_THQ.root")
  infile_TsChan = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_TsChan.root")
  infile_TtChan = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_TtChan.root")
  infile_TtW = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_TtW.root")
  infile_TTW = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_TTW.root")
  infile_TTZ = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_TTZ.root")
  infile_Wjets = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_Wjets.root")
  infile_WW = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_WW.root")
  infile_WZ = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_WZ.root")
  infile_ZZ = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz5mw20/histofile_ZZ.root")

  vars = ["totVecM","fourthJetPt","thirdJetPt","jjdelR","secJetPt","topMass","leadJetPt","wQuarkHt","leadJetbTag","wPaisMass","met","topPt","wPairPt","wwDelR","zlb2DelR","wTopDelR"]

  ## Histos

  totVecM_histo = ROOT.TH1D("totVecM","totVecM", 500, 0.0, 4000.0)
  fourthJetPt_histo = ROOT.TH1D("fourthJetPt","fourthJetPt", 50, 0.0, 150.0)
  thirdJetPt_histo = ROOT.TH1D("thirdJetPt","thirdJetPt", 50, 0.0, 150.0)
  jjdelR_histo = ROOT.TH1D("jjdelR","jjdelR", 50, 0.0, 8.0)
  secJetPt_histo = ROOT.TH1D("secJetPt","secJetPt", 50, 0.0, 150.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)
  _histo = ROOT.TH1D("","", 500, 0.0, 4000.0)

### ALL MC ###
  for event in infile_tZq.Ttree_tZq :

    All_chi2Histo.Fill(chi2,weight)

  for event in infile_DY.Ttree_DYToLL_M50 :


    All_chi2Histo.Fill(chi2,event.EvtWeight)

  for event in infile_TT.Ttree_TT :


  for event in infile_DY2.Ttree_DYToLL_M10To50 :

  for event in infile_TbartChan.Ttree_TbartChan :

  for event in infile_TbartW.Ttree_TbartW :

  for event in infile_THQ.Ttree_THQ :

  for event in infile_TsChan.Ttree_TsChan :

  for event in infile_TtChan.Ttree_TtChan :

  for event in infile_TtW.Ttree_TtW :

  for event in infile_TTW.Ttree_TTW :

  for event in infile_TTZ.Ttree_TTZ :

  for event in infile_Wjets.Ttree_Wjets :

  for event in infile_WW.Ttree_WW :

  for event in infile_WZ.Ttree_WZ :

  for event in infile_ZZ.Ttree_ZZ :

##############

  subprocess.call("mkdir plots/bdtVariables/",shell=True)

  All_topVsWcontrolHisto.SaveAs("plots/chiSquared/All/controlMassPlot.root")

if __name__ == "__main__":
    main()
