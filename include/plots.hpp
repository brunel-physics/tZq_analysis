#ifndef _plots_hpp_
#define _plots_hpp_

#include <string>
#include "AnalysisEvent.hpp"
#include <map>
#include <vector>

typedef struct plot plot;

class TH1F;

class Plots{

 private:
  std::vector<plot> plotPoint;
  //Fill expressions here.
  float fillLepton1Pt(AnalysisEvent*);
  float fillLepton1Eta(AnalysisEvent*);
  float fillLepton2Pt(AnalysisEvent*);
  float fillLepton2Eta(AnalysisEvent*);
  float fillLepton3Pt(AnalysisEvent*);
  float fillLepton3Eta(AnalysisEvent*);
  float fillWbosonQuark1Pt(AnalysisEvent*);
  float fillWbosonQuark1Eta(AnalysisEvent*);
  float fillWbosonQuark2Pt(AnalysisEvent*);
  float fillWbosonQuark2Eta(AnalysisEvent*);
  float fillLepton1Phi(AnalysisEvent*);
  float fillLepton2Phi(AnalysisEvent*);
  float fillLepton3Phi(AnalysisEvent*);
  float fillWbosonQuark1Phi(AnalysisEvent*);
  float fillWbosonQuark2Phi(AnalysisEvent*);
  float fillLepton1RelIso(AnalysisEvent*);
  float fillLepton2RelIso(AnalysisEvent*);
  float fillLepton3RelIso(AnalysisEvent*);
  float fillLepton1MVA(AnalysisEvent*);
  float fillLepton2MVA(AnalysisEvent*);
  float fillLepton3MVA(AnalysisEvent*);
  float getNumberOfJets(AnalysisEvent*);
  float fillLeadingJetPt(AnalysisEvent*);
  float fillLeadingJetEta(AnalysisEvent*);
  float fillLeadingJetPhi(AnalysisEvent*);
  float fillLeadingJetBDisc(AnalysisEvent*);
  float fillSecondJetPt(AnalysisEvent*);
  float fillSecondJetEta(AnalysisEvent*);
  float fillSecondJetPhi(AnalysisEvent*);
  float fillSecondJetBDisc(AnalysisEvent*);
  float fillThirdJetPt(AnalysisEvent*);
  float fillThirdJetEta(AnalysisEvent*);
  float fillThirdJetPhi(AnalysisEvent*);
  float fillThirdJetBDisc(AnalysisEvent*);
  float fillFourthJetPt(AnalysisEvent*);
  float fillFourthJetEta(AnalysisEvent*);
  float fillFourthJetPhi(AnalysisEvent*);
  float fillFourthJetBDisc(AnalysisEvent*);
  float fillMetPlot(AnalysisEvent*);
  float numbBJets(AnalysisEvent*);
  float fillBtagDisc(AnalysisEvent*);
  float fillZLep1Pt(AnalysisEvent*);
  float fillZLep1Eta(AnalysisEvent*);
  float fillZLep2Pt(AnalysisEvent*);
  float fillZLep2Eta(AnalysisEvent*);
  float fillWLepPt(AnalysisEvent*);
  float fillWLepEta(AnalysisEvent*);
  float fillZLep1Phi(AnalysisEvent*);
  float fillZLep2Phi(AnalysisEvent*);
  float fillWLepPhi(AnalysisEvent*);
  float fillZLep1RelIso(AnalysisEvent*);
  float fillZLep2RelIso(AnalysisEvent*);
  float fillWLepRelIso(AnalysisEvent*);
  float fillZLep1MVA(AnalysisEvent*);
  float fillZLep2MVA(AnalysisEvent*);
  float fillWLepMVA(AnalysisEvent*);
  float fillZPairMass(AnalysisEvent*);
  float fillZPairPt(AnalysisEvent*);
  float fillZPairEta(AnalysisEvent*);
  float fillZPairPhi(AnalysisEvent*);
  float fillWPair1Mass(AnalysisEvent*);
  float fillWPair2Mass(AnalysisEvent*);
  float fillWPairMass(AnalysisEvent*);
  float fillWPairPt(AnalysisEvent*);
  float fillWPairEta(AnalysisEvent*);
  float fillWPairPhi(AnalysisEvent*);
  float fillLeptonMass(AnalysisEvent*);
  float fillTopMass(AnalysisEvent*);
  float fillTopPt(AnalysisEvent*);
  float fillTopEta(AnalysisEvent*);
  float fillTopPhi(AnalysisEvent*);
  float fillLepton1D0(AnalysisEvent*);
  float fillLepton2D0(AnalysisEvent*);
  float fillLepton3D0(AnalysisEvent*);
  float fillLepton1DBD0(AnalysisEvent*);
  float fillLepton2DBD0(AnalysisEvent*);
  float fillLepton3DBD0(AnalysisEvent*);
  float fillLepton1BeamSpotCorrectedD0(AnalysisEvent*);
  float fillLepton2BeamSpotCorrectedD0(AnalysisEvent*);
  float fillLepton3BeamSpotCorrectedD0(AnalysisEvent*);
  float fillLepton1InnerTrackD0(AnalysisEvent*);
  float fillLepton2InnerTrackD0(AnalysisEvent*);
  float fillLepton3InnerTrackD0(AnalysisEvent*);
  float fillwTransverseMass(AnalysisEvent*);
  float filljjDelR(AnalysisEvent*);
  float filljjDelPhi(AnalysisEvent*);
  float fillwwDelR(AnalysisEvent*);
  float fillwwDelPhi(AnalysisEvent*);
  float fillZLepDelR(AnalysisEvent*);
  float fillZLepDelPhi(AnalysisEvent*);
  float fillZLep1Quark1DelR(AnalysisEvent*);
  float fillZLep1Quark1DelPhi(AnalysisEvent*);
  float fillZLep1Quark2DelR(AnalysisEvent*);
  float fillZLep1Quark2DelPhi(AnalysisEvent*);
  float fillZLep2Quark1DelR(AnalysisEvent*);
  float fillZLep2Quark1DelPhi(AnalysisEvent*);
  float fillZLep2Quark2DelR(AnalysisEvent*);
  float fillZLep2Quark2DelPhi(AnalysisEvent*);
  float fillZLep1BjetDelR(AnalysisEvent*);
  float fillZLep1BjetDelPhi(AnalysisEvent*);
  float fillZLep2BjetDelR(AnalysisEvent*);
  float fillZLep2BjetDelPhi(AnalysisEvent*);
  float filllbDelR(AnalysisEvent*);
  float filllbDelPhi(AnalysisEvent*);
  float fillLepHt(AnalysisEvent*);
  float fillWquarkHt(AnalysisEvent*);
  float fillJetHt(AnalysisEvent*);
  float fillTotHt(AnalysisEvent*);
  float fillTotHtOverPt(AnalysisEvent*);
  float fillTotPt(AnalysisEvent*);
  float fillTotEta(AnalysisEvent*);
  float fillTotM(AnalysisEvent*);
  float fillWZdelR(AnalysisEvent*);
  float fillWZdelPhi(AnalysisEvent*);
  float fillZquark1DelR(AnalysisEvent*);
  float fillZquark1DelPhi(AnalysisEvent*);
  float fillZquark2DelR(AnalysisEvent*);
  float fillZquark2DelPhi(AnalysisEvent*);
  float fillZtopDelR(AnalysisEvent*);
  float fillZtopDelPhi(AnalysisEvent*);
  float fillZLep1TopDelR(AnalysisEvent*);
  float fillZLep1TopDelPhi(AnalysisEvent*);
  float fillZLep2TopDelR(AnalysisEvent*);
  float fillZLep2TopDelPhi(AnalysisEvent*);

  const bool trileptonChannel_;

 public:
  Plots(std::vector<std::string>,std::vector<float>,std::vector<float>,std::vector<int>,std::vector<std::string>, std::vector<std::string>, std::vector<int>, unsigned, std::string postfixName="", const bool trileptonChannel = true);
  ~Plots();
  void fillAllPlots(AnalysisEvent*, float);
  void saveAllPlots();
  void fillOnePlot(std::string, AnalysisEvent*,float);
  void saveOnePlots(int);
  std::vector<plot> getPlotPoint(){return plotPoint;}
  std::map<std::string, float (Plots::*)(AnalysisEvent*)> getFncPtrMap();


};

struct plot{
  std::string name;
  TH1F* plotHist;
  float (Plots::*fillExp)(AnalysisEvent*);
  std::string xAxisLabel;
  bool fillPlot;
};

#endif // _plots_hpp_ endif

