#A tool to pull plots from the mva inputs and plot gaussians

import ROOT

def main():

#  infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaInput/met0mtw50/histofile_tZq.root")
#  infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaInput/met0mtw50/histofile_TT.root")
  infile = ROOT.TFile.Open("/scratch/data/TopPhysics/mvaInput/met0mtw50/histofile_DYToLL_M50.root")

  topMassHisto = ROOT.TH1D("topMassHisto","Top Mass Histo", 300, 0.0, 300.0)
  wMassHisto = ROOT.TH1D("wMassHisto","W Mass Histo", 300, 0.0, 300.0)
  chi2Histo = ROOT.TH1D("chi2Histo","chi2 Histo", 300, 0.0, 300.0)

  topVsWmassHisto = ROOT.TH2D("topVsWmassHisto", "Top mass vs W Mass; m_{W}; m_{Top};", 300, 0., 300., 300, 0., 300.)
  topVsChi2Histo  = ROOT.TH2D("topVsChi2Histo", "Chi2 vs Top mass; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300.)
  wVsChi2Histo    = ROOT.TH2D("topVsChi2Histo", "Chi2 vs W Mass; m_{W}; #chi^{2}", 300, 0., 300., 300, 0., 300.)

  wVsTopvsChi2Histo = ROOT.TH3D("wVsTopvsChi2Histo", "Chi2 vs Top vs W Masses; m_{W}; m_{Top}; #chi^{2}", 300, 0., 300., 300, 0., 300., 300, 0., 300.)

#  for event in infile.Ttree_tZq :
#  for event in infile.Ttree_TT :
  for event in string infile.Ttree_DYToLL_M50 :
    if ( event.topMass < 200 and event.topMass > 120 ) : topMassHisto.Fill(event.topMass)
    if ( event.wPairMass < 90.385 and event.wPairMass > 70.385 ) : wMassHisto.Fill(event.wPairMass)

    wChi2Term   = (event.wPairMass - 80.3585)/6.0
    topChi2Term = (event.topMass - 173.21)/10.0
    chi2 = wChi2Term*wChi2Term + topChi2Term*topChi2Term

    topVsWmassHisto.Fill(event.wPairMass,event.topMass)
    wVsChi2Histo.Fill(event.wPairMass,chi2)
    topVsChi2Histo.Fill(event.topMass,chi2)

    wVsTopvsChi2Histo.Fill(event.wPairMass,event.topMass,chi2)

    chi2Histo.Fill(chi2)

  topMassHisto.Fit("gaus")
  wMassHisto.Fit("gaus")

  c1 = ROOT.TCanvas("c1","",800,600)
  topMassHisto.Draw()
  c1.Print("topMass.root")

  c2 = ROOT.TCanvas("c2","",800,600)
  wMassHisto.Draw()
  c2.Print("wMass.root")

  c3 = ROOT.TCanvas("c3","",800,600)
  topVsWmassHisto.Draw("colz")
  c3.Print("topVsWmass.root")

  c4 = ROOT.TCanvas("c4","",800,600)
  wVsChi2Histo.Draw("colz")
  c4.Print("wVsChi2.root")

  c5 = ROOT.TCanvas("c5","",800,600)
  topVsChi2Histo.Draw("colz")
  c5.Print("topVsChi2.root")

  c6 = ROOT.TCanvas("c6","",800,600)
  wVsTopvsChi2Histo.Draw("colz")
  c6.Print("wVsTopvsChi2.root")

  c7 = ROOT.TCanvas("c7","",800,600)
  chi2Histo.Draw("colz")
  c7.Print("chi2.root")

if __name__ == "__main__":
    main()
