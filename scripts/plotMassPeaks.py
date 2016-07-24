#A tool to pull plots from the mva inputs and plot gaussians

import ROOT

def main():

  infile_tZq = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaInput/2015/met0mtw20/histofile_tZq.root")
  infile_TT = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaInput/2015/met0mtw20/histofile_TT.root")
  infile_DY = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaInput/2015/met0mtw20/histofile_DYToLL_M50.root")

  # tZq Histos

  tZq_topMassHisto = ROOT.TH1D("tZq_topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  tZq_wMassHisto = ROOT.TH1D("tZq_wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  tZq_chi2Histo = ROOT.TH1D("tZq_chi2Histo","chi2 Histo", 300, 0.0, 300.0)
  tZq_topVsWmassHisto = ROOT.TH2D("tZq_topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  tZq_topVsChi2Histo  = ROOT.TH2D("tZq_topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  tZq_wVsChi2Histo    = ROOT.TH2D("tZq_topVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  tZq_wVsTopvsChi2Histo = ROOT.TH3D("tZq_wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

  ## TT Histos

  TT_topMassHisto = ROOT.TH1D("TT_topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  TT_wMassHisto = ROOT.TH1D("TT_wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  TT_chi2Histo = ROOT.TH1D("TT_chi2Histo","chi2 Histo", 300, 0.0, 300.0)
  TT_topVsWmassHisto = ROOT.TH2D("TT_topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  TT_topVsChi2Histo  = ROOT.TH2D("TT_topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  TT_wVsChi2Histo    = ROOT.TH2D("TT_topVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  TT_wVsTopvsChi2Histo = ROOT.TH3D("TT_wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

  ## DY Histos 

  DY_topMassHisto = ROOT.TH1D("DY_topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  DY_wMassHisto = ROOT.TH1D("DY_wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  DY_chi2Histo = ROOT.TH1D("DY_chi2Histo","chi2 Histo", 300, 0.0, 300.0)
  DY_topVsWmassHisto = ROOT.TH2D("DY_topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  DY_topVsChi2Histo  = ROOT.TH2D("DY_topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  DY_wVsChi2Histo    = ROOT.TH2D("DY_topVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  DY_wVsTopvsChi2Histo = ROOT.TH3D("DY_wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

  for event in infile_tZq.Ttree_tZq :
    if ( event.topMass < 200 and event.topMass > 120 ) : topMassHisto.Fill(event.topMass)
    if ( event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : wMassHisto.Fill(event.wPairMass)

    wChi2Term   = (event.wPairMass - 80.3585)/6.0
    topChi2Term = (event.topMass - 173.21)/10.0
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    tZq_topVsWmassHisto.Fill(event.wPairMass,event.topMass)
    tZq_wVsChi2Histo.Fill(event.wPairMass,chi2)
    tZq_topVsChi2Histo.Fill(event.topMass,chi2)
    tZq_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2)
    tZq_chi2Histo.Fill(chi2)

  for event in infile_TT.Ttree_TT :
    if ( event.topMass < 200 and event.topMass > 120 ) : topMassHisto.Fill(event.topMass)
    if ( event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : wMassHisto.Fill(event.wPairMass)

    wChi2Term   = (event.wPairMass - 80.3585)/6.0
    topChi2Term = (event.topMass - 173.21)/10.0
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    TT_topVsWmassHisto.Fill(event.wPairMass,event.topMass)
    TT_wVsChi2Histo.Fill(event.wPairMass,chi2)
    TT_topVsChi2Histo.Fill(event.topMass,chi2)
    TT_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2)
    TT_chi2Histo.Fill(chi2)

  for event in infile_DY.Ttree_DYToLL_M50 :
    if ( event.topMass < 200 and event.topMass > 120 ) : topMassHisto.Fill(event.topMass)
    if ( event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : wMassHisto.Fill(event.wPairMass)

    wChi2Term   = (event.wPairMass - 80.3585)/6.0
    topChi2Term = (event.topMass - 173.21)/10.0
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    DY_topVsWmassHisto.Fill(event.wPairMass,event.topMass)
    DY_wVsChi2Histo.Fill(event.wPairMass,chi2)
    DY_topVsChi2Histo.Fill(event.topMass,chi2)
    DY_wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2)
    DY_chi2Histo.Fill(chi2)

  subprocess.call("mkdir plots/tZq/ ,shell=True)
  subprocess.call("mkdir plots/DY/ ,shell=True)
  subprocess.call("mkdir plots/TT/ ,shell=True)

  tZq_topMassHisto.Fit("gaus")
  tZq_wMassHisto.Fit("gaus")

  tZq_topMassHisto.SaveAs("plots/tZq/topMass.root")
  tZq_wMassHisto.SaveAs("plots/tZq/wMass.root")
  tZq_topVsWmassHisto.SaveAs("plots/tZq/topVsWmass.root")
  tZq_wVsChi2Histo.SaveAs("plots/tZq/wVsChi2.root")
  tZq_topVsChi2Histo.SaveAs("plots/tZq/topVsChi2.root")
  tZq_wVsTopvsChi2Histo.SaveAs("plots/tZq/wVsTopvsChi2.root")
  tZq_chi2Histo.SaveAs("plots/tZq/chi2.root")

  TT_topMassHisto.Fit("gaus")
  TT_wMassHisto.Fit("gaus")

  TT_topMassHisto.SaveAs("plots/TT/topMass.root")
  TT_wMassHisto.SaveAs("plots/TT/wMass.root")
  TT_topVsWmassHisto.SaveAs("plots/TT/topVsWmass.root")
  TT_wVsChi2Histo.SaveAs("plots/TT/wVsChi2.root")
  TT_topVsChi2Histo.SaveAs("plots/TT/topVsChi2.root")
  TT_wVsTopvsChi2Histo.SaveAs("plots/TT/wVsTopvsChi2.root")
  TT_chi2Histo.SaveAs("plots/TT/chi2.root")

  DY_topMassHisto.Fit("gaus")
  DY_wMassHisto.Fit("gaus")

  DY_topMassHisto.SaveAs("plots_DY/topMass.root")
  DY_wMassHisto.SaveAs("plots_DY/wMass.root")
  DY_topVsWmassHisto.SaveAs("plots_DY/topVsWmass.root")
  DY_wVsChi2Histo.SaveAs("plots_DY/wVsChi2.root")
  DY_topVsChi2Histo.SaveAs("plots_DY/topVsChi2.root")
  DY_wVsTopvsChi2Histo.SaveAs("plots_DY/wVsTopvsChi2.root")
  DY_chi2Histo.SaveAs("plots_DY/chi2.root")

if __name__ == "__main__":
    main()
