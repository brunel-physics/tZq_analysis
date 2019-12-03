#ifndef _makeMVAinputAlgo_hpp
#define _makeMVAinputAlgo_hpp_

#include "jetCorrectionUncertainty.hpp"

#include <map>
#include <unordered_map>
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
    void standardAnalysis(const std::map<std::string, std::string>& listOfMCs,
                          const std::vector<std::string>& systs,
                          const std::vector<std::string>& channels,
                          const bool useSidebandRegion);
    void dataAnalysis(const std::vector<std::string>& channels,
                      const bool useSidebandRegion);
    void sameSignAnalysis(const std::map<std::string, std::string>& listOfMCs,
                          const std::vector<std::string>& channels,
                          const bool useSidebandRegion);
    std::pair<TLorentzVector, TLorentzVector>
        sortOutLeptons(const MvaEvent* tree, const std::string& channel) const;
    std::pair<TLorentzVector, TLorentzVector>
        sortOutHadronicW(const MvaEvent* tree,
                         const int syst,
                         TLorentzVector met,
                         const std::vector<int>& jets) const;
        std::pair<std::vector<int>, std::vector<TLorentzVector>> getJets(
            const MvaEvent* tree, const int syst, TLorentzVector met) const;
    std::pair<std::vector<int>, std::vector<TLorentzVector>>
        getBjets(const MvaEvent* tree,
                 const int syst,
                 TLorentzVector met,
                 const std::vector<int>& jets) const;
    TLorentzVector getJetVec(const MvaEvent* tree,
                             const int index,
                             const float smearValue,
                             TLorentzVector& metVec,
                             const int syst,
                             const bool doMetSmear) const;
    TLorentzVector doUncMet(TLorentzVector met,
                            const TLorentzVector& zLep1,
                            const TLorentzVector& zLep2,
                            const std::vector<TLorentzVector>& jetVecs,
                            const unsigned syst) const;
    void setupBranches(TTree* tree);
    void fillTree(TTree* outTreeSig,
                  TTree* outTreeSdBnd,
                  MvaEvent* tree,
                  const std::string& label,
                  const std::string& channel,
                  const bool SameSignMC = false);

    // variables?

    std::unordered_map<std::string, float> inputVars;
    bool oldMetFlag;
    bool ttbarControlRegion;
    bool useSidebandRegion;
    bool doMC;
    bool doSysts;
    bool doData;
    bool doFakes;
    bool is2016;
    std::string inputDir;
    std::string outputDir;
    std::string era;
};

#endif
