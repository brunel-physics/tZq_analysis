#ifndef _analysisAlgo_hpp_
#define _analysisAlgo_hpp_

#include "cutClass.hpp"
#include "dataset.hpp"
#include "histogramPlotter.hpp"

#include <map>
#include <memory>
#include <vector>

class TH1D;
class TFile;
class TChain;
class TTree;

class AnalysisAlgo
{
    public:
    // Constructor
    AnalysisAlgo();
    ~AnalysisAlgo();

    void parseCommandLineArguements(int argc, char* argv[]);
    void setupSystematics();
    void setupCuts();
    void setupPlots();
    void runMainAnalysis();
    void savePlots();

    private:
    // functions
    std::string channelSetup(unsigned);

    // variables?
    std::string config;
    bool plots;
    bool makeHistos;
    bool useHistos;
    double usePreLumi;
    long nEvents;
    std::string outFolder;
    std::string histoDir;
    std::string postfix;
    std::string channel;
    bool invertLepCut; // For z+jets background estimation
    bool skipData; // utility stuff. True if flags are set and will skip either
                   // data or mc.
    bool skipMC;
    std::string cutConfName;
    std::string plotConfName;
    int numFiles;
    bool makePostLepTree;
    bool makeMVATree;
    bool usePostLepTree;
    bool usebTagWeight;
    int systToRun;
    int channelsToRun;
    bool skipTrig;
    std::string mvaDir;
    bool customJetRegion;
    float metCut;
    float mzCut;
    float mwCut;
    bool is2016_;
    bool doNPLs_;
    bool doZplusCR_;

    std::vector<Dataset> datasets;
    double totalLumi;

    // Cuts stuff
    Cuts* cutObj;

    // Plotting stuff
    std::map<
        std::string,
        std::map<std::string, std::map<std::string, std::shared_ptr<Plots>>>>
        plotsMap;
    std::map<std::string, TH1D*> cutFlowMap;

    std::vector<std::pair<std::string, std::string>> stageNames;

    // A couple of things for plotting. These will soon be set in a config file.
    std::vector<std::string> legOrder;
    std::vector<std::string> plotOrder;
    std::map<std::string, datasetInfo> datasetInfos;

    std::vector<std::string> plotsVec;

    // variables for plotting.
    std::vector<std::string> plotTitles;
    std::vector<std::string> plotNames;
    std::vector<float> xMin;
    std::vector<float> xMax;
    std::vector<int> nBins, cutStage;
    std::vector<std::string> fillExp;
    std::vector<std::string> xAxisLabels;
    std::vector<int> eventNumbers;
    std::vector<unsigned> jetRegVars;

    // Systematic Stuff
    // Making a vector of strings that will give systematics name.
    std::vector<std::string> systNames;
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

    // MC weight stuff
    double sumPositiveWeights_;
    double sumNegativeWeights_;
    double sumNegativeWeightsScaleUp_;
    double sumNegativeWeightsScaleDown_;
};

#endif
