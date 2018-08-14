#ifndef _triggerScaleFactorsAlgo_hpp_
#define _triggerScaleFactorsAlgo_hpp_

#include <vector>
#include "TCanvas.h"
#include "TPad.h"

#include "dataset.hpp"

#include <libconfig.h++>
#include <LHAPDF/LHAPDF.h>

#include <map>

class AnalysisEvent;
class TTree;
class TFile;
class TH1F;
class TH2F;
class TProfile;
class TProfile2D;
class TLorentzVector;

class TriggerScaleFactors{

	public:
        // Constructor
        TriggerScaleFactors();
        ~TriggerScaleFactors();

	void setBranchStatusAll(TTree * chain, bool isMC, std::string triggerFlag);
	void show_usage(std::string name);

	void parseCommandLineArguements (int argc, char* argv[]);
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
        bool jetCuts_;
        bool bCuts_;
	bool applyHltSf_;

	// Command line arguement related variables, replacing hard coded logic
        bool HIP_ERA;	// legacy hard coded variable, name kept for consistency/ease
        bool DO_HIPS;	// legacy hard coded variable, name kept for consistency/ease
	bool isPart1_;  // New command line related variable which helps set HIP_ERA + DO_HIPS
	bool isPart2_;	// New command line related variable which helps set HIP_ERA + DO_HIPS

        bool customElectronCuts_;
        bool customMuonCuts_;
     	std::vector<double> electronCutsVars;
     	std::vector<double> muonCutsVars;

        bool makeJetCuts(AnalysisEvent *event, bool isMC);
        TLorentzVector getJetLVec(AnalysisEvent* event, int index, bool isMC);
        double deltaR(float eta1, float phi1, float eta2, float phi2);
 
	// PU reweighting
	TFile * dataPileupFile;
	TH1F* dataPU;
	TFile * mcPileupFile;
	TH1F* mcPU;
	TFile * systUpFile;
	TH1F* pileupUpHist;
	TFile * systDownFile;
	TH1F* pileupDownHist;
	TH1F* puReweight;
	TH1F* puSystUp;
	TH1F* puSystDown;

	// lepton selection
	std::vector<int> getTightElectrons(AnalysisEvent*);
	std::vector<int> getTightMuons(AnalysisEvent*);
	bool passDileptonSelection(AnalysisEvent*,int);

	// trigger cuts
	bool doubleElectronTriggerCut(AnalysisEvent*, bool);
	bool muonElectronTriggerCut(AnalysisEvent*, bool);
	bool doubleMuonTriggerCut(AnalysisEvent*, bool);
	bool metTriggerCut(AnalysisEvent*);
	bool metFilters(AnalysisEvent*,bool);

	//Efficiencies
	double numberPassedElectrons[2];
	double numberTriggeredDoubleElectrons[2];

	double numberPassedMuons[2];
	double numberTriggeredDoubleMuons[2];

	double numberPassedMuonElectrons[2];
	double numberTriggeredMuonElectrons[2];

	//Systematic variables
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

