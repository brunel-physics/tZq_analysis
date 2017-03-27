#ifndef _debugInfo_hpp_
#define _debugInfo_hpp_

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
class TTree;
class TFile;

class DebugInfo{

	public:
        // Constructor
        DebugInfo();
        ~DebugInfo();

	void setBranchStatusAll(TTree * chain, bool isMC, std::string triggerFlag);
	void show_usage(std::string name);

	void parseCommandLineArguements (int argc, char* argv[]);
	void runMainAnalysis();
	void savePlots();

	private:

        bool is2016_;

	std::string config;
	double usePreLumi;
	long nEvents;
	std::string outFolder;
	std::string postfix;
	bool makePostLepTree;
	bool usePostLepTree;
	int numFiles;

	std::vector<Dataset> datasets;
	double totalLumi;
	double* lumiPtr;
};

#endif

