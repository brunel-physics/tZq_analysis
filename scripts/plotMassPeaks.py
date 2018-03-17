#A tool to pull plots from the mva inputs and plot gaussians

import ROOT
import subprocess
import sys

def main():

  era = "2016"

  massCut = False

  weighted = True
  sigCut = 50.0
  bkgCut = 290.0

  wSigma = 8.0
  topSigma = 30.0

  mzCut = "20"
  mwCut = "9999"

  infile_tZq = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_tZq.root")
  infile_TT = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_TT.root")
  infile_DY = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_DYToLL_M50.root")
  infile_DY_amcatnlo = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_DYToLL_M50_aMCatNLO.root")

  infile_DY2 = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_DYToLL_M10To50.root")
#  infile_DY2_amcatnlo = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_DYToLL_M10To50_aMCatNLO.root")
  infile_TbartChan = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_TbartChan.root")
  infile_TbartW = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_TbartW.root")
  infile_THQ = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_THQ.root")
  infile_TsChan = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_TsChan.root")
  infile_TtChan = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_TtChan.root")
  infile_TtW = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_TtW.root")
  infile_TTW = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_TTW.root")
  infile_TTZ = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_TTZ.root")
  infile_Wjets = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_Wjets.root")
  infile_WW = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_WW.root")
  infile_WZ = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_WZ.root")
  infile_ZZ = ROOT.TFile.Open(" /scratch/data/TopPhysics/mvaDirs/inputs/"+era+"/all/mz"+mzCut+"mw"+mwCut+"/histofile_ZZ.root")


  # tZq Histos

  tZq_topMassHisto = ROOT.TH1D("tZq_topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  tZq_wMassHisto = ROOT.TH1D("tZq_wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  tZq_chi2Histo = ROOT.TH1D("tZq_chi2Histo","chi2 Histo", 300, 0.0, 300.0)

  tZq_topVsWmassHisto = ROOT.TH2D("tZq_topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  tZq_topVsWsignalHisto = ROOT.TH2D("tZq_topVsWsignalHisto", "Top mass vs W Mass for signal region; m_{W}; m_{Top};", 150, 0., 150., 300, 0., 300.)
  tZq_topVsWcontrolHisto = ROOT.TH2D("tZq_topVsWcontrolHisto", "Top mass vs W Mass for control region; m_{W}; m_{Top};", 150, 0., 150., 300, 0., 300.)

  tZq_topVsChi2Histo  = ROOT.TH2D("tZq_topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  tZq_wVsChi2Histo    = ROOT.TH2D("tZq_wVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 150, 0., 150., 300, 0., 300.)

  tZq_wVsTopvsChi2Histo = ROOT.TH3D("tZq_wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

  tZq_metVsWmassHisto = ROOT.TH2D("tZq_metVsWmassHisto", "MET vs W Mass; m_{W}; MET", 150, 0., 150., 300, 0., 300.)

  ## TT Histos

  TT_topMassHisto = ROOT.TH1D("TT_topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  TT_wMassHisto = ROOT.TH1D("TT_wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  TT_chi2Histo = ROOT.TH1D("TT_chi2Histo","chi2 Histo", 300, 0.0, 300.0)

  TT_topVsWmassHisto = ROOT.TH2D("TT_topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  TT_topVsWsignalHisto = ROOT.TH2D("TT_topVsWsignalHisto", "Top mass vs W Mass for signal region, m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  TT_topVsWcontrolHisto = ROOT.TH2D("TT_topVsWcontrolHisto", "Top mass vs W Mass for control region; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)

  TT_topVsChi2Histo  = ROOT.TH2D("TT_topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  TT_wVsChi2Histo    = ROOT.TH2D("TT_wVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 300, 0., 300., 300, 0., 300.)

  TT_wVsTopvsChi2Histo = ROOT.TH3D("TT_wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

  TT_metVsWmassHisto = ROOT.TH2D("TT_metVsWmassHisto", "MET vs W Mass; m_{W}; MET", 150, 0., 150., 300, 0., 300.)

  ## DY Histos

  DY_topMassHisto = ROOT.TH1D("DY_topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  DY_wMassHisto = ROOT.TH1D("DY_wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  DY_chi2Histo = ROOT.TH1D("DY_chi2Histo","chi2 Histo", 300, 0.0, 300.0)

  DY_topVsWmassHisto = ROOT.TH2D("DY_topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  DY_topVsWsignalHisto = ROOT.TH2D("DY_topVsWsignalHisto", "Top mass vs W Mass for signal region; m_{W}; m_{Top};", 150, 0., 150., 300, 0., 300.)
  DY_topVsWcontrolHisto = ROOT.TH2D("DY_topVsWcontrolHisto", "Top mass vs W Mass for control region; m_{W}; m_{Top};", 150, 0., 150., 300, 0., 300.)

  DY_topVsChi2Histo  = ROOT.TH2D("DY_topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  DY_wVsChi2Histo    = ROOT.TH2D("DY_wVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 300, 0., 300., 300, 0., 300.)

  DY_wVsTopvsChi2Histo = ROOT.TH3D("DY_wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

  DY_metVsWmassHisto = ROOT.TH2D("DY_metVsWmassHisto", "MET vs W Mass; m_{W}; MET", 150, 0., 150., 300, 0., 300.)

  ## DY aMCatNLO Histos

  DY_amcatnlo_topMassHisto = ROOT.TH1D("DY_amcatnlo_topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  DY_amcatnlo_wMassHisto = ROOT.TH1D("DY_amcatnlo_wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  DY_amcatnlo_chi2Histo = ROOT.TH1D("DY_amcatnlo_chi2Histo","chi2 Histo", 300, 0.0, 300.0)

  DY_amcatnlo_topVsWmassHisto = ROOT.TH2D("DY_amcatnlo_topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  DY_amcatnlo_topVsWsignalHisto = ROOT.TH2D("DY_amcatnlo_topVsWsignalHisto", "Top mass vs W Mass for signal region; m_{W}; m_{Top};", 150, 0., 150., 300, 0., 300.)
  DY_amcatnlo_topVsWcontrolHisto = ROOT.TH2D("DY_amcatnlo_topVsWcontrolHisto", "Top mass vs W Mass for control region; m_{W}; m_{Top};", 150, 0., 150., 300, 0., 300.)

  DY_amcatnlo_topVsChi2Histo  = ROOT.TH2D("DY_amcatnlo_topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  DY_amcatnlo_wVsChi2Histo    = ROOT.TH2D("DY_amcatnlo_wVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 300, 0., 300., 300, 0., 300.)

  DY_amcatnlo_wVsTopvsChi2Histo = ROOT.TH3D("DY_amcatnlo_wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

  DY_amcatnlo_metVsWmassHisto = ROOT.TH2D("DY_amcatnlo_metVsWmassHisto", "MET vs W Mass; m_{W}; MET", 150, 0., 150., 300, 0., 300.)

  ## All MC Histos

  All_topMassHisto = ROOT.TH1D("All_topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  All_wMassHisto = ROOT.TH1D("All_wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  All_chi2Histo = ROOT.TH1D("All_chi2Histo","chi2 Histo", 300, 0.0, 300.0)

  All_topVsWmassHisto = ROOT.TH2D("All_topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  All_topVsWsignalHisto = ROOT.TH2D("All_topVsWsignalHisto", "Top mass vs W Mass for signal region, m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  All_topVsWcontrolHisto = ROOT.TH2D("All_topVsWcontrolHisto", "Top mass vs W Mass for control region; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)

  All_topVsChi2Histo  = ROOT.TH2D("All_topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  All_wVsChi2Histo    = ROOT.TH2D("All_wVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 300, 0., 300., 300, 0., 300.)

  All_wVsTopvsChi2Histo = ROOT.TH3D("All_wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

  All_metVsWmassHisto = ROOT.TH2D("All_metVsWmassHisto", "MET vs W Mass; m_{W}; MET", 150, 0., 150., 300, 0., 300.)

  for event in infile_tZq.Ttree_tZq :
#    if (event.Channel != 0) : continue
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 220 and event.topMass > 110 ) : tZq_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : tZq_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : tZq_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : tZq_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    tZq_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    tZq_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    tZq_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    tZq_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    tZq_chi2Histo.Fill(chi2,weight)

    tZq_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : tZq_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : tZq_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TT.Ttree_TT :
    weight = 1
#    if (event.Channel != 0) : continue
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 220 and event.topMass > 110 ) : TT_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : TT_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : TT_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : TT_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    TT_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    TT_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    TT_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    TT_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    TT_chi2Histo.Fill(chi2,weight)

    TT_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : TT_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : TT_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_DY.Ttree_DYToLL_M50 :
#    if (event.Channel != 0) : continue
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 220 and event.topMass > 110 ) : DY_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : DY_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : DY_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : DY_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    DY_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    DY_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    DY_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    DY_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    DY_chi2Histo.Fill(chi2,weight)

    DY_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : DY_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : DY_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_DY_amcatnlo.Ttree_DYToLL_M50_aMCatNLO :
#    if (event.Channel != 0) : continue
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 220 and event.topMass > 110 ) : DY_amcatnlo_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : DY_amcatnlo_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : DY_amcatnlo_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : DY_amcatnlo_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    DY_amcatnlo_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    DY_amcatnlo_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    DY_amcatnlo_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    DY_amcatnlo_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    DY_amcatnlo_chi2Histo.Fill(chi2,weight)

    DY_amcatnlo_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : DY_amcatnlo_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : DY_amcatnlo_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)


### ALL MC ###
  for event in infile_tZq.Ttree_tZq :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_DY.Ttree_DYToLL_M50 :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TT.Ttree_TT :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_DY2.Ttree_DYToLL_M10To50 :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TbartChan.Ttree_TbartChan :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TbartW.Ttree_TbartW :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_THQ.Ttree_THQ :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)


    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TsChan.Ttree_TsChan :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TtChan.Ttree_TtChan :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TtW.Ttree_TtW :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TTW.Ttree_TTW :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_TTZ.Ttree_TTZ :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_Wjets.Ttree_Wjets :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_WW.Ttree_WW :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_WZ.Ttree_WZ :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

  for event in infile_ZZ.Ttree_ZZ :
    weight = 1
    if (weighted) : weight = event.EvtWeight

    if ( massCut and event.topMass < 200 and event.topMass > 110 ) : All_topMassHisto.Fill(event.topMass,weight)
    if ( massCut and event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : All_wMassHisto.Fill(event.wPairMass,weight)
    if not ( massCut ) : All_topMassHisto.Fill(event.topMass,weight)
    if not ( massCut ) : All_wMassHisto.Fill(event.wPairMass,weight)

    wChi2Term   = (event.wPairMass - 80.3585)/wSigma
    topChi2Term = (event.topMass - 173.21)/topSigma
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    All_topVsWmassHisto.Fill(event.wPairMass,event.topMass,weight)
    All_wVsChi2Histo.Fill(event.wPairMass,chi2,weight)
    All_topVsChi2Histo.Fill(event.topMass,chi2,weight)
    All_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2,weight)
    All_chi2Histo.Fill(chi2,weight)

    All_metVsWmassHisto.Fill(event.wPairMass,event.met,weight)

    if ( chi2 < sigCut ) : All_topVsWsignalHisto.Fill(event.wPairMass,event.topMass,weight)
    if ( chi2 >= sigCut and chi2 < bkgCut ) : All_topVsWcontrolHisto.Fill(event.wPairMass,event.topMass,weight)

##############

  subprocess.call("mkdir plots/chiSquared/",shell=True)
  subprocess.call("mkdir plots/chiSquared/tZq/",shell=True)
  subprocess.call("mkdir plots/chiSquared/TT/",shell=True)
  subprocess.call("mkdir plots/chiSquared/DY/",shell=True)
  subprocess.call("mkdir plots/chiSquared/DY_amcatnlo/",shell=True)
  subprocess.call("mkdir plots/chiSquared/All/",shell=True)

  tZq_topMassHisto.Fit("gaus")
  tZq_wMassHisto.Fit("gaus")

  tZq_topVsWsignalHisto.GetZaxis().SetRangeUser(0.0, 0.1)
  tZq_topVsWcontrolHisto.GetZaxis().SetRangeUser(0.0, 0.1)

  tZq_topMassHisto.SaveAs("plots/chiSquared/tZq/topMass.root")
  tZq_wMassHisto.SaveAs("plots/chiSquared/tZq/wMass.root")
  tZq_topVsWmassHisto.SaveAs("plots/chiSquared/tZq/topVsWmass.root")
  tZq_wVsChi2Histo.SaveAs("plots/chiSquared/tZq/wVsChi2.root")
  tZq_topVsChi2Histo.SaveAs("plots/chiSquared/tZq/topVsChi2.root")
  tZq_wVsTopvsChi2Histo.SaveAs("plots/chiSquared/tZq/wVsTopvsChi2.root")
  tZq_chi2Histo.SaveAs("plots/chiSquared/tZq/chi2.root")
  tZq_topVsWsignalHisto.SaveAs("plots/chiSquared/tZq/signalMassPlot.root")
  tZq_topVsWcontrolHisto.SaveAs("plots/chiSquared/tZq/controlMassPlot.root")
  tZq_metVsWmassHisto.SaveAs("plots/chiSquared/tZq/metVsWmass.root")

  TT_topMassHisto.Fit("gaus")
  TT_wMassHisto.Fit("gaus")

  TT_topMassHisto.SaveAs("plots/chiSquared/TT/topMass.root")
  TT_wMassHisto.SaveAs("plots/chiSquared/TT/wMass.root")
  TT_topVsWmassHisto.SaveAs("plots/chiSquared/TT/topVsWmass.root")
  TT_wVsChi2Histo.SaveAs("plots/chiSquared/TT/wVsChi2.root")
  TT_topVsChi2Histo.SaveAs("plots/chiSquared/TT/topVsChi2.root")
  TT_wVsTopvsChi2Histo.SaveAs("plots/chiSquared/TT/wVsTopvsChi2.root")
  TT_chi2Histo.SaveAs("plots/chiSquared/TT/chi2.root")
  TT_topVsWsignalHisto.SaveAs("plots/chiSquared/TT/signalMassPlot.root")
  TT_topVsWcontrolHisto.SaveAs("plots/chiSquared/TT/controlMassPlot.root")
  TT_metVsWmassHisto.SaveAs("plots/chiSquared/TT/metVsWmass.root")

  DY_topMassHisto.Fit("gaus")
  DY_wMassHisto.Fit("gaus")

  DY_topMassHisto.SaveAs("plots/chiSquared/DY/topMass.root")
  DY_wMassHisto.SaveAs("plots/chiSquared/DY/wMass.root")
  DY_topVsWmassHisto.SaveAs("plots/chiSquared/DY/topVsWmass.root")
  DY_wVsChi2Histo.SaveAs("plots/chiSquared/DY/wVsChi2.root")
  DY_topVsChi2Histo.SaveAs("plots/chiSquared/DY/topVsChi2.root")
  DY_wVsTopvsChi2Histo.SaveAs("plots/chiSquared/DY/wVsTopvsChi2.root")
  DY_chi2Histo.SaveAs("plots/chiSquared/DY/chi2.root")
  DY_topVsWsignalHisto.SaveAs("plots/chiSquared/DY/signalMassPlot.root")
  DY_topVsWcontrolHisto.SaveAs("plots/chiSquared/DY/controlMassPlot.root")
  DY_metVsWmassHisto.SaveAs("plots/chiSquared/DY/metVsWmass.root")

  DY_amcatnlo_topMassHisto.SaveAs("plots/chiSquared/DY_amcatnlo/topMass.root")
  DY_amcatnlo_wMassHisto.SaveAs("plots/chiSquared/DY_amcatnlo/wMass.root")
  DY_amcatnlo_topVsWmassHisto.SaveAs("plots/chiSquared/DY_amcatnlo/topVsWmass.root")
  DY_amcatnlo_wVsChi2Histo.SaveAs("plots/chiSquared/DY_amcatnlo/wVsChi2.root")
  DY_amcatnlo_topVsChi2Histo.SaveAs("plots/chiSquared/DY_amcatnlo/topVsChi2.root")
  DY_amcatnlo_wVsTopvsChi2Histo.SaveAs("plots/chiSquared/DY_amcatnlo/wVsTopvsChi2.root")
  DY_amcatnlo_chi2Histo.SaveAs("plots/chiSquared/DY_amcatnlo/chi2.root")
  DY_amcatnlo_topVsWsignalHisto.SaveAs("plots/chiSquared/DY_amcatnlo/signalMassPlot.root")
  DY_amcatnlo_topVsWcontrolHisto.SaveAs("plots/chiSquared/DY_amcatnlo/controlMassPlot.root")
  DY_amcatnlo_metVsWmassHisto.SaveAs("plots/chiSquared/DY_amcatnlo/metVsWmass.root")

  All_topMassHisto.Fit("gaus")
  All_wMassHisto.Fit("gaus")

  All_topMassHisto.SaveAs("plots/chiSquared/All/topMass.root")
  All_wMassHisto.SaveAs("plots/chiSquared/All/wMass.root")
  All_topVsWmassHisto.SaveAs("plots/chiSquared/All/topVsWmass.root")
  All_wVsChi2Histo.SaveAs("plots/chiSquared/All/wVsChi2.root")
  All_topVsChi2Histo.SaveAs("plots/chiSquared/All/topVsChi2.root")
  All_wVsTopvsChi2Histo.SaveAs("plots/chiSquared/All/wVsTopvsChi2.root")
  All_chi2Histo.SaveAs("plots/chiSquared/All/chi2.root")
  All_topVsWsignalHisto.SaveAs("plots/chiSquared/All/signalMassPlot.root")
  All_topVsWcontrolHisto.SaveAs("plots/chiSquared/All/controlMassPlot.root")
  All_metVsWmassHisto.SaveAs("plots/chiSquared/All/metVsWmass.root")

##############

  print "##############"

  totalEvents_tZq = tZq_chi2Histo.Integral()
  sumEvents_tZq = 0.0
  tZq_flag1 = 0
  tZq_flag1a = 0 
  tZq_flag1b = 0 
  tZq_flag2 = 0
  tZq_flag3 = 0
  for x in range(0,300):
    sumEvents_tZq += tZq_chi2Histo.GetBinContent(x)
    if sumEvents_tZq >= 0.6827*totalEvents_tZq and tZq_flag1 == 0:
      print "one sigma bound for tZq chi2 at: ", x, " : " , sumEvents_tZq
      tZq_flag1 = 1
    if sumEvents_tZq >= 0.80*totalEvents_tZq and tZq_flag1a == 0:
      print "80% bound for tZq chi2 at: ", x, " : " , sumEvents_tZq
      tZq_flag1a = 1
    if sumEvents_tZq >= 0.90*totalEvents_tZq and tZq_flag1b == 0:
      print "90% bound for tZq chi2 at: ", x, " : " , sumEvents_tZq
      tZq_flag1b = 1
    if sumEvents_tZq >= 0.9545*totalEvents_tZq and tZq_flag2 == 0:
      print "two sigma bound for tZq chi2 at: ", x, " : " , sumEvents_tZq
      tZq_flag2 = 1
    if sumEvents_tZq >= 0.9973*totalEvents_tZq and tZq_flag3 == 0:
      print "three sigma bound for tZq chi2 at: ", x, " : " , sumEvents_tZq
      tZq_flag3 = 1

  print "##############"

##############

  totalEvents_DY = DY_chi2Histo.Integral()
  sumEvents_DY = 0.0
  DY_flag1 = 0
  DY_flag1a = 0
  DY_flag1b = 0
  DY_flag2 = 0
  DY_flag3 = 0
  for x in range(0,300):
    sumEvents_DY += DY_chi2Histo.GetBinContent(x)
    if sumEvents_DY >= 0.6827*totalEvents_DY and DY_flag1 == 0:
      print "one sigma bound for DY chi2 at: ", x, " : " , sumEvents_DY
      DY_flag1 = 1
    if sumEvents_DY >= 0.80*totalEvents_DY and DY_flag1a == 0:
      print "80% bound for DY chi2 at: ", x, " : " , sumEvents_DY
      DY_flag1a = 1
    if sumEvents_DY >= 0.90*totalEvents_DY and DY_flag1b == 0:
      print "90% bound for DY chi2 at: ", x, " : " , sumEvents_DY
      DY_flag1b = 1
    if sumEvents_DY >= 0.9545*totalEvents_DY and DY_flag2 == 0:
      print "two sigma bound for DY chi2 at: ", x, " : " , sumEvents_DY
      DY_flag2 = 1
    if sumEvents_DY >= 0.9973*totalEvents_DY and DY_flag3 == 0:
      print "three sigma bound for DY chi2 at: ", x, " : " , sumEvents_DY
      DY_flag3 = 1

  print "##############"

##############

  totalEvents_TT = TT_chi2Histo.Integral()
  sumEvents_TT = 0.0
  TT_flag1 = 0
  TT_flag1a = 0
  TT_flag1b = 0
  TT_flag2 = 0
  TT_flag3 = 0
  for x in range(0,300):
    sumEvents_TT += TT_chi2Histo.GetBinContent(x)
    if sumEvents_TT >= 0.6827*totalEvents_TT and TT_flag1 == 0:
      print "one sigma bound for TT chi2 at: ", x, " : " , sumEvents_TT
      TT_flag1 = 1
    if sumEvents_TT >= 0.80*totalEvents_TT and TT_flag1a == 0:
      print "80% bound for TT chi2 at: ", x, " : " , sumEvents_TT
      TT_flag1a = 1
    if sumEvents_TT >= 0.90*totalEvents_TT and TT_flag1b == 0:
      print "90% bound for TT chi2 at: ", x, " : " , sumEvents_TT
      TT_flag1b = 1
    if sumEvents_TT >= 0.9545*totalEvents_TT and TT_flag2 == 0:
      print "two sigma bound for TT chi2 at: ", x, " : " , sumEvents_TT
      TT_flag2 = 1
    if sumEvents_TT >= 0.9973*totalEvents_TT and TT_flag3 == 0:
      print "three sigma bound for TT chi2 at: ", x, " : " , sumEvents_TT
      TT_flag3 = 1

  print "##############"

##############

if __name__ == "__main__":
    main()
