#ifndef _makeMVAinputAlgo_hpp
#define _makeMVAinputAlgo_hpp_

#include "jetCorrectionUncertainty.hpp"

#include <map>
#include <vector>

class TTree;
class MvaEvent;
class TLorentzVector;

class MakeMvaInputs
{
    public:
    // Constructor
    MakeMvaInputs();
    ~MakeMvaInputs();

    void parseCommandLineArguements(int argc, char* argv[]);
    void runMainAnalysis();

    private:
    // functions
    // Simple deltaR function, because the reco namespace doesn't work or
    // something
    double deltaR(float, float, float, float);
    std::pair<TLorentzVector, TLorentzVector>
        sortOutLeptons(MvaEvent* tree, std::string channel);
    std::pair<TLorentzVector, TLorentzVector> sortOutHadronicW(MvaEvent* tree);
    std::pair<std::vector<int>, std::vector<TLorentzVector>>
        getJets(MvaEvent* tree, int syst, TLorentzVector met);
    std::pair<std::vector<int>, std::vector<TLorentzVector>> getBjets(
        MvaEvent* tree, int syst, TLorentzVector met, std::vector<int> jets);
    TLorentzVector getJetVec(MvaEvent* tree,
                             int index,
                             float smearValue,
                             TLorentzVector metVec,
                             int syst,
                             bool doMetSmear);
    TLorentzVector doUncMet(TLorentzVector met,
                            TLorentzVector zLep1,
                            TLorentzVector zLep2,
                            std::vector<TLorentzVector> jetVecs,
                            uint syst);

    std::map<std::string, float> setupInputVars();
    void setupBranches(TTree*, std::map<std::string, float>);
    void fillTree(TTree* outTreeSig,
                  TTree* outTreeSdBnd,
//                  std::map<std::string, float>&,
                  MvaEvent* tree,
                  std::string label,
                  std::string channel,
                  bool SameSignMC = false);

    // variables?

    JetCorrectionUncertainty jetUnc;

    std::map<std::string, float> inputVars;
    bool oldMetFlag;
    bool ttbarControlRegion;
    bool useSidebandRegion;
    std::string inputDir;
    std::string outputDir;
};

#endif

