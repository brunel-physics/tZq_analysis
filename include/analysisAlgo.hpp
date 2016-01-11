#ifndef _analysisAlgo_hpp_
#define _analysisAlgo_hpp_

class AnalysisAlgo{

	public:
	// Constructor
	AnalysisAlgo();
	~AnalysisAlgo();

	double zptSF(TString channel, float zpt);
	void setBranchStatusAll(TTree * chain, bool isMC, std::string triggerFlag);
	static void show_usage(std::string name);

	private:

	std::string config;
 	bool plots;
	double usePreLumi;
	long nEvents;
	std::string outFolder";
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
	bool trileptonChannel;

	// variables for plotting. 
	std::vector<std::string> plotNames;
	std::vector<float> xMin;
	std::vector<float> xMax;
	std::vector<int> nBins,cutStage;
	std::vector<std::string> fillExp;
	std::vector<std::string> xAxisLabels;
	std::vector<int> eventNumbers;
	std::vector<unsigned int> jetRegVars;


};

#endif

