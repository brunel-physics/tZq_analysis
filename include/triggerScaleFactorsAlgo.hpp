#ifndef _triggerScaleFactorsAlgo_hpp_
#define _triggerScaleFactorsAlgo_hpp_

#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"

#include "dataset.hpp"

#include <libconfig.h++>
#include <LHAPDF/LHAPDF.h>

#include <map>

class AnalysisEvent;

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
 	bool plots;
	double usePreLumi;
	long nEvents;
	std::string outFolder;
	std::string postfix;
	int numFiles;
	std::vector<int> emptyVector;

	std::vector<Dataset> datasets;
	double totalLumi;
	double* lumiPtr;

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
	bool passDileptonSelection(AnalysisEvent*,std::vector<int>,std::vector<int>);

	// trigger cuts
	bool doubleElectronTriggerCut(AnalysisEvent*);
	bool muonElectronTriggerCut(AnalysisEvent*);
	bool doubleMuonTriggerCut(AnalysisEvent*);

	//
	unsigned numberPassedElectrons[2];
	unsigned numberTriggeredElectrons[2];
	unsigned numberPassedMuons[2];
	unsigned numberTriggeredMuons[2];
};

#endif

