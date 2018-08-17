#ifndef _cutClass_hpp_
#define _cutClass_hpp_

#include "AnalysisEvent.hpp"
#include "BTagCalibrationStandalone.hpp"
#include "RoccoR.h"
#include "TLorentzVector.h"
#include "plots.hpp"

#include <fstream>
#include <map>
#include <vector>

class TH1F;
class TH2F;
class TH2D;
class TH3D;
class TGraphAsymmErrors;

class Cuts
{
    private:
    bool makeLeptonCuts(AnalysisEvent*,
                        float*,
                        std::map<std::string, Plots*>,
                        TH1F*,
                        int syst = 0,
                        bool isControl = false);
    bool invertIsoCut(AnalysisEvent*,
                      float*,
                      std::map<std::string, Plots*>,
                      TH1F*);
    std::pair<std::vector<int>, std::vector<float>>
        makeJetCuts(AnalysisEvent*, int, float*, bool isProper = true);
    std::vector<int> makeMetCuts(AnalysisEvent*);
    std::vector<int> makeBCuts(AnalysisEvent*, std::vector<int>, int syst = 0);
    std::vector<int>
        makeLooseBCuts(AnalysisEvent*, std::vector<int>, int syst = 0);
    std::vector<int> makeCCuts(AnalysisEvent*, std::vector<int>);

    std::vector<int> getTightEles(AnalysisEvent* event);
    std::vector<int> getInvIsoEles(AnalysisEvent* event);
    std::vector<int> getLooseEles(AnalysisEvent* event);
    std::vector<int> getTightMuons(AnalysisEvent* event);
    std::vector<int> getInvIsoMuons(AnalysisEvent* event);
    std::vector<int> getLooseMuons(AnalysisEvent* event);
    float getTrileptonZCand(AnalysisEvent*, std::vector<int>, std::vector<int>);
    bool getDileptonZCand(AnalysisEvent*, std::vector<int>, std::vector<int>);
    float getWbosonQuarksCand(AnalysisEvent*, std::vector<int>, int syst = 0);

    std::vector<std::pair<int, int>> getSynchDileptonCandidates(
        AnalysisEvent*, std::vector<int>, std::vector<int>);

    float getTopMass(AnalysisEvent*);
    bool triggerCuts(AnalysisEvent*, float*, int syst = 0);
    bool metFilters(AnalysisEvent*);

    double getChiSquared(double wMass = 0.0, double topMass = 0.0);

    // Method for running the synchronisation with Jeremy.
    bool synchCuts(AnalysisEvent* event, float* eventWeight);
    int getLooseLepsNum(AnalysisEvent* event); // Mimic preselection skims
    int getLooseElecs(AnalysisEvent* event);
    int getLooseMus(AnalysisEvent* event);
    // Methods for running tW synch
    std::vector<int> getSynchEles(AnalysisEvent* event);
    std::vector<int> getSynchMus(AnalysisEvent* event);

    // Method to do ttbar cuts for the dilepton background estimation
    bool ttbarCuts(AnalysisEvent* event,
                   float*,
                   std::map<std::string, Plots*>,
                   TH1F*,
                   int);

    // Simple deltaR function, because the reco namespace doesn't work or
    // something
    double deltaR(float, float, float, float);
    void dumpToFile(AnalysisEvent* event, int);

    // Function to get lepton SF
    float getLeptonWeight(AnalysisEvent*, int syst = 0);
    float eleSF(double, double, int syst = 0);
    float muonSF(double, double, int syst = 0);

    float singleElectronTriggerSF(double, double, int syst = 0);
    float singleMuonTriggerSF(double, double, int syst = 0);
    float muonTriggerSF(double, double, double, double, int syst = 0);

    // set to true to fill in histograms/spit out other info
    bool doPlots_;
    bool fillCutFlow_; // Fill cut flows
    bool invertLepCut_; // For background estimation
    bool synchCutFlow_; // For synch
    bool singleEventInfoDump_; // For dropping info on event for synching.
    bool makeEventDump_;
    const bool isFCNC_;
    const bool isCtag_;
    const bool is2016_;

    // Tight electron cuts
    unsigned numTightEle_;
    float tightElePt_;
    float tightElePtLeading_;
    float tightEleEta_;
    float tightEled0_;
    int tightEleMissLayers_;
    bool tightEleCheckPhotonVeto_;
    float tightEleMVA0_;
    float tightEleMVA1_;
    float tightEleMVA2_;
    float tightEleRelIso_;

    // Loose electron cuts
    unsigned numLooseEle_;
    float looseElePt_;
    float looseElePtLeading_;
    float looseEleEta_;
    float looseEleMVA0_;
    float looseEleMVA1_;
    float looseEleMVA2_;
    float looseEleRelIso_;

    // Tight muon cuts
    unsigned numTightMu_;
    float tightMuonPt_;
    float tightMuonPtLeading_;
    float tightMuonEta_;
    float tightMuonRelIso_;

    // Loose muon cuts
    unsigned numLooseMu_;
    float looseMuonPt_;
    float looseMuonPtLeading_;
    float looseMuonEta_;
    float looseMuonRelIso_;

    // z and w inv cuts
    float invZMassCut_;
    float invWMassCut_;

    // Tight jet cuts
    unsigned numJets_;
    unsigned maxJets_;
    float jetPt_;
    float jetEta_;
    int jetNConsts_;
    bool jetIDDo_;

    // B-Disc cut
    unsigned numbJets_;
    unsigned maxbJets_;
    float bDiscCut_;
    float bLooseDiscCut_;
    float bDiscSynchCut_;

    // C-Disc cut
    unsigned numcJets_;
    unsigned maxcJets_;
    float cVsLDiscCut_;
    float cVsBDiscCut_;

    // Rochester Corrections
    RoccoR rc_;

    // Temporary jet smearing prop varaible until a more elegant solution is
    // done.
    float tempSmearValue_;

    // lumi for pre-hip and post-hip era
    float lumiRunsBCDEF_;
    float lumiRunsGH_;

    // Some things that will be used for JEC uncertainties.
    std::vector<float> ptMinJEC_;
    std::vector<float> ptMaxJEC_;
    std::vector<float> etaMinJEC_;
    std::vector<float> etaMaxJEC_;
    std::vector<std::vector<float>> jecSFUp_;
    std::vector<std::vector<float>> jecSFDown_;
    void initialiseJECCors();
    float getJECUncertainty(float, float, int);
    TLorentzVector
        getJetLVec(AnalysisEvent*, int, int, bool initialRun = false);
    std::pair<float, float> jet2016SFs(float);
    std::pair<double, double> jet2017SFs(const double eta) const;

    // Histogram to be used in synchronisation.
    TH1F* synchCutFlowHist_;
    TH1I* synchNumEles_;
    TH1I* synchNumMus_;
    TH1I* synchMuonCutFlow_;
    TH1F* synchCutTopMassHist_;

    std::ofstream topMassEventDump_;
    std::ofstream step0EventDump_;
    std::ofstream step2EventDump_;
    std::ofstream step4EventDump_;
    std::ofstream step6EventDump_;
    std::ofstream step9EventDump_;
    // Sets whether to do MC or data cuts. Set every time a new dataset is
    // processed in the main loop.
    bool isMC_;
    std::string triggerFlag_;
    std::string postfixName_;
    // Set the flag used to reject non-prompt leptons when making the NPL shapes
    // for plotting purposes
    bool isNPL_;
    // Set the flag used to use the Z+jets CR
    bool isZplusCR_;
    // Set flag and vars for gen level cuts
    bool doGenMassCuts_;
    bool doGenPtCuts_;
    float minGenMassCut_;
    float maxGenMassCut_;
    float minGenPtCut_;
    float maxGenPtCut_;

    // For producing post-lepsel skims
    TTree* postLepSelTree_;

    // For removing trigger cuts. Will be set to false by default
    bool skipTrigger_;

    // For making b-tagging efficiencies. Needed for reweighting and
    // systematics.
    bool makeBTagEffPlots_;
    // And the efficiency plots.
    std::vector<TH2D*> bTagEffPlots_;
    bool getBTagWeight_;

    // bTag callibration for SFs
    BTagCalibration calib2016;
    BTagCalibration calib2017;
    BTagCalibrationReader lightReader;
    BTagCalibrationReader charmReader;
    BTagCalibrationReader beautyReader;

    double getBweight_backup(const int flavour,
                             const int type,
                             const double pt) const;

    void getBWeight(AnalysisEvent*,
                    TLorentzVector,
                    int,
                    float*,
                    float*,
                    float*,
                    float*,
                    float*,
                    float*,
                    float*,
                    float*);

    // met and mtw cut values
    float metCut_;
    float metDileptonCut_;
    float mTWCutSynch_;

    // top mass cut values
    float TopMassCutLower_;
    float TopMassCutUpper_;

    // Sets trigger from config file
    std::string cutConfTrigLabel_;

    TFile* electronHltFile;
    TFile* electronSFsFile;
    TFile* electronRecoFile;
    TH2F* h_eleHlt;
    TH2F* h_eleSFs;
    TH2F* h_eleReco;

    TFile* muonHltFile1;
    TFile* muonHltFile2;
    TFile* muonIDsFile1;
    TFile* muonIsoFile1;
    TFile* muonIDsFile2;
    TFile* muonIsoFile2;
    TH2F* h_muonHlt1;
    TH2F* h_muonHlt2;
    TH2F* h_muonIDs1;
    TH2F* h_muonIDs2;
    TH2F* h_muonPFiso1;
    TH2F* h_muonPFiso2;

    public:
    Cuts(bool,
         bool,
         bool,
         bool,
         const bool,
         const bool,
         const bool,
         const bool);
    ~Cuts();
    bool makeCuts(
        AnalysisEvent*, float*, std::map<std::string, Plots*>, TH1F*, int);
    void setTightEle(float pt = 20, float eta = 2.5, float d0 = 0.04);
    void setMC(bool isMC)
    {
        isMC_ = isMC;
    }
    void setCloneTree(TTree* tree)
    {
        postLepSelTree_ = tree;
    }
    void setNumLeps(int tightMu, int looseMu, int tightEle, int looseEle)
    {
        numTightEle_ = tightEle;
        numLooseEle_ = looseEle;
        numTightMu_ = tightMu;
        numLooseMu_ = looseMu;
    }
    void setCutConfTrigLabel(std::string newLabel)
    {
        cutConfTrigLabel_ = newLabel;
    }
    void setInvLepCut(bool invLep)
    {
        invertLepCut_ = invLep;
    }
    void setTriggerFlag(std::string triggerFlag)
    {
        triggerFlag_ = triggerFlag;
    }
    void setBTagPlots(std::vector<TH2D*> vec, bool makePlotsOrRead)
    {
        makeBTagEffPlots_ = makePlotsOrRead;
        bTagEffPlots_ = vec;
        getBTagWeight_ = !makePlotsOrRead;
    }
    void setSkipTrig(bool skip)
    {
        skipTrigger_ = skip;
    }
    void setMetCut(float cut)
    {
        metCut_ = cut;
        metDileptonCut_ = cut;
    }
    void setMWCut(float cut)
    {
        invWMassCut_ = cut;
    }
    void setMZCut(float cut)
    {
        invZMassCut_ = cut;
    }
    void setJetRegion(int nJets, int nBets, int maxJets, int maxBJets)
    {
        numJets_ = nJets;
        numbJets_ = nBets;
        maxJets_ = maxJets;
        maxbJets_ = maxBJets;
    }
    bool parse_config(std::string);
    void dumpLeptonInfo(AnalysisEvent*);
    void dumpLooseLepInfo(AnalysisEvent*);
    TH1F* getSynchCutFlow();
    int numFound()
    {
        return synchCutFlowHist_->GetBinContent(4);
    }
    void setEventInfoFlag(bool flag)
    {
        singleEventInfoDump_ = flag;
    }
    void setNplFlag(bool isNPL)
    {
        isNPL_ = isNPL;
    }
    void setZplusControlRegionFlag(bool isZplusCR)
    {
        isZplusCR_ = isZplusCR;
    }

    void setGenMassCuts(float minCut)
    {
        doGenMassCuts_ = true;
        minGenMassCut_ = minCut;
    }
    void setGenMassCuts(float minCut, float maxCut)
    {
        doGenMassCuts_ = true;
        minGenMassCut_ = minCut;
        maxGenMassCut_ = maxCut;
    }
    void setGenPtCuts(float minCut)
    {
        doGenPtCuts_ = true;
        minGenPtCut_ = minCut;
    }
    void setGenPtCuts(float minCut, float maxCut)
    {
        doGenMassCuts_ = true;
        minGenPtCut_ = minCut;
        maxGenPtCut_ = maxCut;
    }
};

#endif
