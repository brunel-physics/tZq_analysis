#ifndef _analysisAlgo_hpp_
#define _analysisAlgo_hpp_

#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"

#include "cutClass.hpp"
#include "histogramPlotter.hpp"
#include "dataset.hpp"

#include <libconfig.h++>
#include <LHAPDF/LHAPDF.h>

#include <map>

class AnalysisAlgo{

	public:
        // Constructor
        AnalysisAlgo();
        ~AnalysisAlgo();

	double zptSF(std::string channel, float zpt);
	void setBranchStatusAll(TTree * chain, bool isMC, std::string triggerFlag);
	void show_usage(std::string name);

	void parseCommandLineArguements (int argc, char* argv[]);
	void setupSystematics();
	void setupCuts();
	void setupPlots();
	void runMainAnalysis();
	void savePlots();

	private:

	std::string config;
 	bool plots;
	double usePreLumi;
	long nEvents;
	std::string outFolder;
	std::string postfix;
	std::string channel;
	bool infoDump;
	bool invertIsoCut; //For z+jets background estimation
	bool synchCutFlow; // For synch
	bool skipData; //utility stuff. True if flags are set and will skip either data or mc.
	bool skipMC;
	std::string* cutConfName;
	std::string* plotConfName;
	int numFiles;
	bool readEventList;
	bool dumpEventNumbers;
	bool makePostLepTree;
	bool makeMVATree;
	bool usePostLepTree;
	bool usebTagWeight;
	int systToRun;
	bool makeBTagEffPlots;
	int channelsToRun;
	bool skipTrig;
	std::string mvaDir;
	bool customJetRegion;
	float metCut;
	float mtwCut;
	bool trileptonChannel_;
	bool isFCNC_;

	std::vector<Dataset> datasets;
	double totalLumi;
	double* lumiPtr;

        // Cuts stuff
        Cuts * cutObj;

	// Plotting stuff
	std::map<std::string, std::map<std::string, std::map<std::string, Plots*> > > plotsMap;
	std::map<std::string, TH1F*> cutFlowMap;

	std::vector<std::string> stageNames;

	//A couple of things for plotting. These will soon be set in a config file.
	std::vector<std::string> legOrder;
	std::vector<std::string > plotOrder;
	std::map<std::string, datasetInfo> datasetInfos;

	std::vector<std::string> plotsVec;

	// variables for plotting. 
	std::vector<std::string> plotNames;
	std::vector<float> xMin;
	std::vector<float> xMax;
	std::vector<int> nBins,cutStage;
	std::vector<std::string> fillExp;
	std::vector<std::string> xAxisLabels;
	std::vector<int> eventNumbers;
	std::vector<unsigned int> jetRegVars;

	// Systematic Stuff
	//Making a vector of strings that will give systematics name.
	std::vector<std::string> systNames;
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
	
};

#endif

