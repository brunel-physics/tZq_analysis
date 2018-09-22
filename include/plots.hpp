#ifndef _plots_hpp_
#define _plots_hpp_

#include <string>
#include "AnalysisEvent.hpp"
#include <map>
#include <vector>

typedef struct plot plot;

class TH1D;

class Plots{

 private:
  std::vector<plot> plotPoint;
  //Fill expressions here.
  double fillLepton1Pt(AnalysisEvent*);
  double fillLepton1Eta(AnalysisEvent*);
  double fillLepton2Pt(AnalysisEvent*);
  double fillLepton2Eta(AnalysisEvent*);
  double fillLepton3Pt(AnalysisEvent*);
  double fillLepton3Eta(AnalysisEvent*);
  double fillWbosonQuark1Pt(AnalysisEvent*);
  double fillWbosonQuark1Eta(AnalysisEvent*);
  double fillWbosonQuark2Pt(AnalysisEvent*);
  double fillWbosonQuark2Eta(AnalysisEvent*);
  double fillLepton1Phi(AnalysisEvent*);
  double fillLepton2Phi(AnalysisEvent*);
  double fillLepton3Phi(AnalysisEvent*);
  double fillWbosonQuark1Phi(AnalysisEvent*);
  double fillWbosonQuark2Phi(AnalysisEvent*);
  double fillLepton1RelIso(AnalysisEvent*);
  double fillLepton2RelIso(AnalysisEvent*);
  double fillLepton3RelIso(AnalysisEvent*);
  double fillLepton1MVA(AnalysisEvent*);
  double fillLepton2MVA(AnalysisEvent*);
  double fillLepton3MVA(AnalysisEvent*);
  double getNumberOfJets(AnalysisEvent*);
  double fillTotalJetMass(AnalysisEvent*);
  double fillTotalJetPt(AnalysisEvent*);
  double fillTotalJetEta(AnalysisEvent*);
  double fillTotalJetPhi(AnalysisEvent*);
  double fillLeadingJetPt(AnalysisEvent*);
  double fillLeadingJetEta(AnalysisEvent*);
  double fillLeadingJetPhi(AnalysisEvent*);
  double fillLeadingJetBDisc(AnalysisEvent*);
  double fillSecondJetPt(AnalysisEvent*);
  double fillSecondJetEta(AnalysisEvent*);
  double fillSecondJetPhi(AnalysisEvent*);
  double fillSecondJetBDisc(AnalysisEvent*);
  double fillThirdJetPt(AnalysisEvent*);
  double fillThirdJetEta(AnalysisEvent*);
  double fillThirdJetPhi(AnalysisEvent*);
  double fillThirdJetBDisc(AnalysisEvent*);
  double fillFourthJetPt(AnalysisEvent*);
  double fillFourthJetEta(AnalysisEvent*);
  double fillFourthJetPhi(AnalysisEvent*);
  double fillFourthJetBDisc(AnalysisEvent*);
  double fillMetPlot(AnalysisEvent*);
  double numbBJets(AnalysisEvent*);
  double fillBtagDisc(AnalysisEvent*);
  double fillZLep1Pt(AnalysisEvent*);
  double fillZLep1Eta(AnalysisEvent*);
  double fillZLep2Pt(AnalysisEvent*);
  double fillZLep2Eta(AnalysisEvent*);
  double fillWLepPt(AnalysisEvent*);
  double fillWLepEta(AnalysisEvent*);
  double fillZLep1Phi(AnalysisEvent*);
  double fillZLep2Phi(AnalysisEvent*);
  double fillWLepPhi(AnalysisEvent*);
  double fillZLep1RelIso(AnalysisEvent*);
  double fillZLep2RelIso(AnalysisEvent*);
  double fillWLepRelIso(AnalysisEvent*);
  double fillZLep1MVA(AnalysisEvent*);
  double fillZLep2MVA(AnalysisEvent*);
  double fillWLepMVA(AnalysisEvent*);
  double fillZPairMass(AnalysisEvent*);
  double fillZPairPt(AnalysisEvent*);
  double fillZPairEta(AnalysisEvent*);
  double fillZPairPhi(AnalysisEvent*);
  double fillWPair1Mass(AnalysisEvent*);
  double fillWPair2Mass(AnalysisEvent*);
  double fillWPairMass(AnalysisEvent*);
  double fillWPairPt(AnalysisEvent*);
  double fillWPairEta(AnalysisEvent*);
  double fillWPairPhi(AnalysisEvent*);
  double fillLeptonMass(AnalysisEvent*);
  double fillTopMass(AnalysisEvent*);
  double fillTopPt(AnalysisEvent*);
  double fillTopEta(AnalysisEvent*);
  double fillTopPhi(AnalysisEvent*);
  double fillLepton1D0(AnalysisEvent*);
  double fillLepton2D0(AnalysisEvent*);
  double fillLepton3D0(AnalysisEvent*);
  double fillLepton1DBD0(AnalysisEvent*);
  double fillLepton2DBD0(AnalysisEvent*);
  double fillLepton3DBD0(AnalysisEvent*);
  double fillLepton1BeamSpotCorrectedD0(AnalysisEvent*);
  double fillLepton2BeamSpotCorrectedD0(AnalysisEvent*);
  double fillLepton3BeamSpotCorrectedD0(AnalysisEvent*);
  double fillLepton1InnerTrackD0(AnalysisEvent*);
  double fillLepton2InnerTrackD0(AnalysisEvent*);
  double fillLepton3InnerTrackD0(AnalysisEvent*);
  double fillwTransverseMass(AnalysisEvent*);
  double filljjDelR(AnalysisEvent*);
  double filljjDelPhi(AnalysisEvent*);
  double fillwwDelR(AnalysisEvent*);
  double fillwwDelPhi(AnalysisEvent*);
  double fillZLepDelR(AnalysisEvent*);
  double fillZLepDelPhi(AnalysisEvent*);
  double fillZLep1Quark1DelR(AnalysisEvent*);
  double fillZLep1Quark1DelPhi(AnalysisEvent*);
  double fillZLep1Quark2DelR(AnalysisEvent*);
  double fillZLep1Quark2DelPhi(AnalysisEvent*);
  double fillZLep2Quark1DelR(AnalysisEvent*);
  double fillZLep2Quark1DelPhi(AnalysisEvent*);
  double fillZLep2Quark2DelR(AnalysisEvent*);
  double fillZLep2Quark2DelPhi(AnalysisEvent*);
  double fillZLep1BjetDelR(AnalysisEvent*);
  double fillZLep1BjetDelPhi(AnalysisEvent*);
  double fillZLep2BjetDelR(AnalysisEvent*);
  double fillZLep2BjetDelPhi(AnalysisEvent*);
  double filllbDelR(AnalysisEvent*);
  double filllbDelPhi(AnalysisEvent*);
  double fillLepHt(AnalysisEvent*);
  double fillWquarkHt(AnalysisEvent*);
  double fillJetHt(AnalysisEvent*);
  double fillTotHt(AnalysisEvent*);
  double fillTotHtOverPt(AnalysisEvent*);
  double fillTotPt(AnalysisEvent*);
  double fillTotEta(AnalysisEvent*);
  double fillTotM(AnalysisEvent*);
  double fillWZdelR(AnalysisEvent*);
  double fillWZdelPhi(AnalysisEvent*);
  double fillZquark1DelR(AnalysisEvent*);
  double fillZquark1DelPhi(AnalysisEvent*);
  double fillZquark2DelR(AnalysisEvent*);
  double fillZquark2DelPhi(AnalysisEvent*);
  double fillZtopDelR(AnalysisEvent*);
  double fillZtopDelPhi(AnalysisEvent*);
  double fillZLep1TopDelR(AnalysisEvent*);
  double fillZLep1TopDelPhi(AnalysisEvent*);
  double fillZLep2TopDelR(AnalysisEvent*);
  double fillZLep2TopDelPhi(AnalysisEvent*);

  double fillzjminR(AnalysisEvent*);
  double fillzTopDelR(AnalysisEvent*);

  const bool trileptonChannel_;

 public:
  Plots(std::vector<std::string>,std::vector<std::string>,std::vector<float>,std::vector<float>,std::vector<int>,std::vector<std::string>, std::vector<std::string>, std::vector<int>, unsigned, std::string postfixName="", const bool trileptonChannel = true);
  ~Plots();
  void fillAllPlots(AnalysisEvent*, double);
  void saveAllPlots();
  void fillOnePlot(std::string, AnalysisEvent*,double);
  void saveOnePlots(int);
  std::vector<plot> getPlotPoint(){return plotPoint;}
  std::map<std::string, double (Plots::*)(AnalysisEvent*)> getFncPtrMap();


};

struct plot{
  std::string name;
  std::string title;
  TH1D* plotHist;
  double (Plots::*fillExp)(AnalysisEvent*);
  std::string xAxisLabel;
  bool fillPlot;
};

#endif // _plots_hpp_ endif

