#ifndef _triggerScaleFactorsAlgo_hpp_
#define _triggerScaleFactorsAlgo_hpp_

#include "TCanvas.h"
#include "TPad.h"
#include "dataset.hpp"

#include <LHAPDF/LHAPDF.h>
#include <TH1D.h>
#include <map>
#include <vector>

class AnalysisEvent;
class TTree;
class TFile;
class TH1F;
class TH2F;
class TProfile;
class TProfile2D;
class TLorentzVector;

class TriggerScaleFactors
{
    public:
    // Constructor
    TriggerScaleFactors();
    ~TriggerScaleFactors();

    void setBranchStatusAll(TTree* chain, bool isMC, std::string triggerFlag);
    void show_usage(std::string name);

    void parseCommandLineArguements(int argc, char* argv[]);
    void runMainAnalysis();
    void savePlots();

    private:
    std::string config;
    double usePreLumi;
    long nEvents;
    std::string outFolder;
    std::string postfix;
    int numFiles;

    std::vector<Dataset> datasets;
    double totalLumi;
    double* lumiPtr;

    bool is2016_;
    bool zCuts_;
    bool jetCuts_;
    bool bCuts_;
    bool applyHltSf_;

    // Command line arguement related variables, replacing hard coded logic
    bool HIP_ERA; // legacy hard coded variable, name kept for consistency/ease
    bool DO_HIPS; // legacy hard coded variable, name kept for consistency/ease
    bool isPart1_; // New command line related variable which helps set HIP_ERA
                   // + DO_HIPS
    bool isPart2_; // New command line related variable which helps set HIP_ERA
                   // + DO_HIPS

    bool customElectronCuts_;
    bool customMuonCuts_;
    std::vector<double> electronCutsVars;
    std::vector<double> muonCutsVars;

    bool makeJetCuts(AnalysisEvent& event, const bool isMC) const;
    TLorentzVector getJetLVec(AnalysisEvent& event,
                              const int index,
                              const bool isMC) const;
    double deltaR(const float eta1,
                  const float phi1,
                  const float eta2,
                  const float phi2) const;

    // PU reweighting
    TFile* dataPileupFile;
    TH1D* dataPU;
    TFile* mcPileupFile;
    TH1D* mcPU;
    TFile* systUpFile;
    TH1D* pileupUpHist;
    TFile* systDownFile;
    TH1D* pileupDownHist;
    TH1D* puReweight;
    TH1D* puSystUp;
    TH1D* puSystDown;

    // lepton selection
    std::vector<int> getTightElectrons(const AnalysisEvent& event) const;
    std::vector<int> getTightMuons(const AnalysisEvent& event) const;
    bool passDileptonSelection(AnalysisEvent& event,
                               const int nElectrons) const;

    // trigger cuts
    bool metTriggerCut(const AnalysisEvent& event) const;
    bool metFilters(const AnalysisEvent& event, const bool isMC) const;

    // Efficiencies
    double numberPassedElectrons[2];
    double numberTriggeredDoubleElectrons[2];

    double numberPassedMuons[2];
    double numberTriggeredDoubleMuons[2];

    double numberPassedMuonElectrons[2];
    double numberTriggeredMuonElectrons[2];

    // Systematic variables
    double numberSelectedElectrons[2];
    double numberSelectedMuons[2];
    double numberSelectedMuonElectrons[2];

    double numberSelectedDoubleElectronsTriggered[2];
    double numberSelectedDoubleMuonsTriggered[2];
    double numberSelectedMuonElectronsTriggered[2]; // Double MuonEG

    // Plots for turn on curve studies
    TProfile* p_electron1_pT_MC;
    TProfile* p_electron1_eta_MC;
    TProfile* p_electron2_pT_MC;
    TProfile* p_electron2_eta_MC;
    TProfile* p_muon1_pT_MC;
    TProfile* p_muon1_eta_MC;
    TProfile* p_muon2_pT_MC;
    TProfile* p_muon2_eta_MC;
    TProfile* p_muonElectron1_pT_MC;
    TProfile* p_muonElectron1_eta_MC;
    TProfile* p_muonElectron2_pT_MC;
    TProfile* p_muonElectron2_eta_MC;

    TProfile* p_electron1_pT_data;
    TProfile* p_electron1_eta_data;
    TProfile* p_electron2_pT_data;
    TProfile* p_electron2_eta_data;
    TProfile* p_muon1_pT_data;
    TProfile* p_muon1_eta_data;
    TProfile* p_muon2_pT_data;
    TProfile* p_muon2_eta_data;
    TProfile* p_muonElectron1_pT_data;
    TProfile* p_muonElectron1_eta_data;
    TProfile* p_muonElectron2_pT_data;
    TProfile* p_muonElectron2_eta_data;

    // Plots used to make SFs
    TProfile2D* p_electrons_pT_MC;
    TProfile2D* p_electrons_eta_MC;
    TProfile2D* p_muons_pT_MC;
    TProfile2D* p_muons_eta_MC;
    TProfile2D* p_muonElectrons_pT_MC;
    TProfile2D* p_muonElectrons_eta_MC;

    TProfile2D* p_electrons_pT_data;
    TProfile2D* p_electrons_eta_data;
    TProfile2D* p_muons_pT_data;
    TProfile2D* p_muons_eta_data;
    TProfile2D* p_muonElectrons_pT_data;
    TProfile2D* p_muonElectrons_eta_data;

    TFile* muonHltFile1;
    TFile* muonHltFile2;
    TH2F* h_muonHlt1;
    TH2F* h_muonHlt2;
};

#endif
