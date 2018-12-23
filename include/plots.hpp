#ifndef _plots_hpp_
#define _plots_hpp_

#include "AnalysisEvent.hpp"

#include <functional>
#include <map>
#include <string>
#include <vector>

typedef struct plot plot;

class TH1D;

class Plots
{
    private:
    std::vector<plot> plotPoint;
    // Fill expressions here.
    static float fillLepton1Pt(AnalysisEvent*);
    static float fillLepton1Eta(AnalysisEvent*);
    static float fillLepton2Pt(AnalysisEvent*);
    static float fillLepton2Eta(AnalysisEvent*);
    static float fillLepton3Pt(AnalysisEvent*);
    static float fillLepton3Eta(AnalysisEvent*);
    static float fillWbosonQuark1Pt(AnalysisEvent*);
    static float fillWbosonQuark1Eta(AnalysisEvent*);
    static float fillWbosonQuark2Pt(AnalysisEvent*);
    static float fillWbosonQuark2Eta(AnalysisEvent*);
    static float fillLepton1Phi(AnalysisEvent*);
    static float fillLepton2Phi(AnalysisEvent*);
    static float fillLepton3Phi(AnalysisEvent*);
    static float fillWbosonQuark1Phi(AnalysisEvent*);
    static float fillWbosonQuark2Phi(AnalysisEvent*);
    static float fillLepton1RelIso(AnalysisEvent*);
    static float fillLepton2RelIso(AnalysisEvent*);
    static float fillLepton3RelIso(AnalysisEvent*);
    static float fillLepton1MVA(AnalysisEvent*);
    static float fillLepton2MVA(AnalysisEvent*);
    static float fillLepton3MVA(AnalysisEvent*);
    static float getNumberOfJets(AnalysisEvent*);
    static float fillTotalJetMass(AnalysisEvent*);
    static float fillTotalJetPt(AnalysisEvent*);
    static float fillTotalJetEta(AnalysisEvent*);
    static float fillTotalJetPhi(AnalysisEvent*);
    static float fillLeadingJetPt(AnalysisEvent*);
    static float fillLeadingJetEta(AnalysisEvent*);
    static float fillLeadingJetPhi(AnalysisEvent*);
    static float fillLeadingJetBDisc(AnalysisEvent*);
    static float fillSecondJetPt(AnalysisEvent*);
    static float fillSecondJetEta(AnalysisEvent*);
    static float fillSecondJetPhi(AnalysisEvent*);
    static float fillSecondJetBDisc(AnalysisEvent*);
    static float fillThirdJetPt(AnalysisEvent*);
    static float fillThirdJetEta(AnalysisEvent*);
    static float fillThirdJetPhi(AnalysisEvent*);
    static float fillThirdJetBDisc(AnalysisEvent*);
    static float fillFourthJetPt(AnalysisEvent*);
    static float fillFourthJetEta(AnalysisEvent*);
    static float fillFourthJetPhi(AnalysisEvent*);
    static float fillFourthJetBDisc(AnalysisEvent*);
    static float fillMetPlot(AnalysisEvent*);
    static float numbBJets(AnalysisEvent*);
    static float fillBtagDisc(AnalysisEvent*);
    static float fillZLep1Pt(AnalysisEvent*);
    static float fillZLep1Eta(AnalysisEvent*);
    static float fillZLep2Pt(AnalysisEvent*);
    static float fillZLep2Eta(AnalysisEvent*);
    static float fillWLepPt(AnalysisEvent*);
    static float fillWLepEta(AnalysisEvent*);
    static float fillZLep1Phi(AnalysisEvent*);
    static float fillZLep2Phi(AnalysisEvent*);
    static float fillWLepPhi(AnalysisEvent*);
    static float fillZLep1RelIso(AnalysisEvent*);
    static float fillZLep2RelIso(AnalysisEvent*);
    static float fillWLepRelIso(AnalysisEvent*);
    static float fillZLep1MVA(AnalysisEvent*);
    static float fillZLep2MVA(AnalysisEvent*);
    static float fillWLepMVA(AnalysisEvent*);
    static float fillZPairMass(AnalysisEvent*);
    static float fillZPairPt(AnalysisEvent*);
    static float fillZPairEta(AnalysisEvent*);
    static float fillZPairPhi(AnalysisEvent*);
    static float fillWPair1Mass(AnalysisEvent*);
    static float fillWPair2Mass(AnalysisEvent*);
    static float fillWPairMass(AnalysisEvent*);
    static float fillWPairPt(AnalysisEvent*);
    static float fillWPairEta(AnalysisEvent*);
    static float fillWPairPhi(AnalysisEvent*);
    static float fillLeptonMass(AnalysisEvent*);
    static float fillTopMass(AnalysisEvent*);
    static float fillTopPt(AnalysisEvent*);
    static float fillTopEta(AnalysisEvent*);
    static float fillTopPhi(AnalysisEvent*);
    static float fillLepton1D0(AnalysisEvent*);
    static float fillLepton2D0(AnalysisEvent*);
    static float fillLepton3D0(AnalysisEvent*);
    static float fillLepton1DBD0(AnalysisEvent*);
    static float fillLepton2DBD0(AnalysisEvent*);
    static float fillLepton3DBD0(AnalysisEvent*);
    static float fillLepton1BeamSpotCorrectedD0(AnalysisEvent*);
    static float fillLepton2BeamSpotCorrectedD0(AnalysisEvent*);
    static float fillLepton3BeamSpotCorrectedD0(AnalysisEvent*);
    static float fillLepton1InnerTrackD0(AnalysisEvent*);
    static float fillLepton2InnerTrackD0(AnalysisEvent*);
    static float fillLepton3InnerTrackD0(AnalysisEvent*);
    static float fillwTransverseMass(AnalysisEvent*);
    static float filljjDelR(AnalysisEvent*);
    static float filljjDelPhi(AnalysisEvent*);
    static float fillwwDelR(AnalysisEvent*);
    static float fillwwDelPhi(AnalysisEvent*);
    static float fillZLepDelR(AnalysisEvent*);
    static float fillZLepDelPhi(AnalysisEvent*);
    static float fillZLep1Quark1DelR(AnalysisEvent*);
    static float fillZLep1Quark1DelPhi(AnalysisEvent*);
    static float fillZLep1Quark2DelR(AnalysisEvent*);
    static float fillZLep1Quark2DelPhi(AnalysisEvent*);
    static float fillZLep2Quark1DelR(AnalysisEvent*);
    static float fillZLep2Quark1DelPhi(AnalysisEvent*);
    static float fillZLep2Quark2DelR(AnalysisEvent*);
    static float fillZLep2Quark2DelPhi(AnalysisEvent*);
    static float fillZLep1BjetDelR(AnalysisEvent*);
    static float fillZLep1BjetDelPhi(AnalysisEvent*);
    static float fillZLep2BjetDelR(AnalysisEvent*);
    static float fillZLep2BjetDelPhi(AnalysisEvent*);
    static float filllbDelR(AnalysisEvent*);
    static float filllbDelPhi(AnalysisEvent*);
    static float fillLepHt(AnalysisEvent*);
    static float fillWquarkHt(AnalysisEvent*);
    static float fillJetHt(AnalysisEvent*);
    static float fillTotHt(AnalysisEvent*);
    static float fillTotHtOverPt(AnalysisEvent*);
    static float fillTotPt(AnalysisEvent*);
    static float fillTotEta(AnalysisEvent*);
    static float fillTotM(AnalysisEvent*);
    static float fillWZdelR(AnalysisEvent*);
    static float fillWZdelPhi(AnalysisEvent*);
    static float fillZquark1DelR(AnalysisEvent*);
    static float fillZquark1DelPhi(AnalysisEvent*);
    static float fillZquark2DelR(AnalysisEvent*);
    static float fillZquark2DelPhi(AnalysisEvent*);
    static float fillZtopDelR(AnalysisEvent*);
    static float fillZtopDelPhi(AnalysisEvent*);
    static float fillZLep1TopDelR(AnalysisEvent*);
    static float fillZLep1TopDelPhi(AnalysisEvent*);
    static float fillZLep2TopDelR(AnalysisEvent*);
    static float fillZLep2TopDelPhi(AnalysisEvent*);

    public:
    Plots(std::vector<std::string>,
          std::vector<std::string>,
          std::vector<float>,
          std::vector<float>,
          std::vector<int>,
          std::vector<std::string>,
          std::vector<std::string>,
          std::vector<int>,
          unsigned,
          std::string postfixName = "");
    ~Plots();
    void fillAllPlots(AnalysisEvent*, float);
    void saveAllPlots();
    void fillOnePlot(std::string, AnalysisEvent*, float);
    void saveOnePlots(int);
    std::vector<plot> getPlotPoint()
    {
        return plotPoint;
    }
    std::map<std::string, std::function<float(AnalysisEvent*)>> getFncMap();
};

struct plot
{
    std::string name;
    std::string title;
    TH1D* plotHist;
    std::function<float(AnalysisEvent*)> fillExp;
    std::string xAxisLabel;
    bool fillPlot;
};

#endif // _plots_hpp_ endif
