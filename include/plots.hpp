#ifndef _plots_hpp_
#define _plots_hpp_

#include <string>
#include "TH1.h"
#include "AnalysisEvent.hpp"
#include <map>
#include <vector>

typedef struct plot plot;

class Plots{
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
  float fillSecondJetPt(AnalysisEvent*);
  float fillSecondJetEta(AnalysisEvent*);
  float fillSecondJetPhi(AnalysisEvent*);
  float fillMetPlot(AnalysisEvent*);
  float numbBJets(AnalysisEvent*);
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
  float fillWPair1Mass(AnalysisEvent*);
  float fillWPair2Mass(AnalysisEvent*);
  float fillLeptonMass(AnalysisEvent*);
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
  float fillzLepDelR(AnalysisEvent*);
  float fillzLepDelPhi(AnalysisEvent*);
  float filllbDelR(AnalysisEvent*);
  float filllbDelPhi(AnalysisEvent*);
  const bool trileptonChannel_;

 public:
  Plots(std::vector<std::string>,std::vector<float>,std::vector<float>,std::vector<int>,std::vector<std::string>, std::vector<std::string>, std::vector<int>, unsigned int, std::string postfixName="", const bool trileptonChannel = true);
  ~Plots();
  void fillAllPlots(AnalysisEvent*, float);
  void saveAllPlots();
  void fillOnePlot(std::string, AnalysisEvent*,float);
  void saveOnePlots(int);
  std::vector<plot> getPlotPoint(){return plotPoint;};
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

