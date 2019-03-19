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

class TH1D;
class TH2F;
class TH2D;
class TH3D;
class TGraphAsymmErrors;

class Cuts
{
    private:
    bool makeLeptonCuts(AnalysisEvent& event,
                        float& eventWeight,
                        std::map<std::string, std::shared_ptr<Plots>>& plotMap,
                        TH1D& cutFlow,
                        const int syst,
                        const bool isControl = false);
    bool invertIsoCut(AnalysisEvent&,
                      float&,
                      std::map<std::string, std::shared_ptr<Plots>>&,
                      TH1D&);
    std::pair<std::vector<int>, std::vector<float>>
        makeJetCuts(AnalysisEvent& event,
                    const int syst,
                    float& eventWeight,
                    const bool isProper = true);
    std::vector<int> makeMetCuts(AnalysisEvent&);
    std::vector<int> makeBCuts(AnalysisEvent& event,
                               const std::vector<int> jets,
                               const int syst = 0);
    std::vector<int> makeLooseBCuts(AnalysisEvent& event,
                                    const std::vector<int> jets,
                                    const int syst);
    std::vector<int> makeCCuts(AnalysisEvent&, std::vector<int>);

    std::vector<int> getTightEles(const AnalysisEvent& event) const;
    std::vector<int> getInvIsoEles(const AnalysisEvent& event) const;
    std::vector<int> getLooseEles(const AnalysisEvent& event) const;
    std::vector<int> getTightMuons(const AnalysisEvent& event) const;
    std::vector<int> getInvIsoMuons(const AnalysisEvent& event) const;
    std::vector<int> getLooseMuons(const AnalysisEvent& event) const;
    float getTrileptonZCand(AnalysisEvent& event,
                            const std::vector<int> electrons,
                            const std::vector<int> muons) const;
    bool getDileptonZCand(AnalysisEvent& event,
                          const std::vector<int> electrons,
                          const std::vector<int> muons) const;
    float getWbosonQuarksCand(AnalysisEvent& event,
                              const std::vector<int> jets,
                              const int syst);

    std::vector<std::pair<int, int>>
        getSynchDileptonCandidates(const AnalysisEvent& event,
                                   const std::vector<int> eles,
                                   const std::vector<int> mus) const;

    float getTopMass(const AnalysisEvent& event) const;
    bool triggerCuts(const AnalysisEvent& event,
                     float& eventWeight,
                     const int syst = 0) const;
    bool metFilters(const AnalysisEvent& event) const;

    // Method for running the synchronisation with Jeremy.
    bool synchCuts(AnalysisEvent& event, float& eventWeight);
    int getLooseLepsNum(
        const AnalysisEvent& event) const; // Mimic preselection skims
    int getLooseElecs(const AnalysisEvent& event) const;
    int getLooseMus(const AnalysisEvent& event) const;
    // Methods for running tW synch
    std::vector<int> getSynchEles(const AnalysisEvent& event) const;
    std::vector<int> getSynchMus(const AnalysisEvent& event) const;

    // Method to do ttbar cuts for the dilepton background estimation
    bool ttbarCuts(AnalysisEvent& event,
                   float& eventWeight,
                   std::map<std::string, std::shared_ptr<Plots>>& plotMap,
                   TH1D& cutFlow,
                   const int systToRun);

    // Simple deltaR function, because the reco namespace doesn't work or
    // something
    double deltaR(const float eta1,
                  const float phi1,
                  const float eta2,
                  const float phi2) const;
    void dumpToFile(AnalysisEvent& event, const int step);

    // Function to get lepton SF
    float getLeptonWeight(const AnalysisEvent& event, const int syst) const;
    float eleSF(const double pt, const double eta, const int syst) const;
    float muonSF(const double pt, const double eta, const int syst) const;

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
    double tightElePt_;
    double tightElePtLeading_;
    double tightEleEta_;
    double tightEled0_;
    int tightEleMissLayers_;
    bool tightEleCheckPhotonVeto_;
    float tightEleMVA0_;
    float tightEleMVA1_;
    float tightEleMVA2_;
    double tightEleRelIso_;

    // Loose electron cuts
    unsigned numLooseEle_;
    double looseElePt_;
    double looseElePtLeading_;
    double looseEleEta_;
    float looseEleMVA0_;
    float looseEleMVA1_;
    float looseEleMVA2_;
    double looseEleRelIso_;

    // Tight muon cuts
    unsigned numTightMu_;
    double tightMuonPt_;
    double tightMuonPtLeading_;
    double tightMuonEta_;
    double tightMuonRelIso_;

    // Loose muon cuts
    unsigned numLooseMu_;
    double looseMuonPt_;
    double looseMuonPtLeading_;
    double looseMuonEta_;
    double looseMuonRelIso_;

    // z and w inv cuts
    float invZMassCut_;
    float invWMassCut_;

    // Tight jet cuts
    unsigned numJets_;
    unsigned maxJets_;
    double jetPt_;
    double jetEta_;
    int jetNConsts_;
    bool jetIDDo_;

    // B-Disc cut
    unsigned numbJets_;
    unsigned maxbJets_;
    float maxbJetEta_;
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
    float getJECUncertainty(const float pt,
                            const float eta,
                            const int syst) const;
    TLorentzVector getJetLVec(AnalysisEvent& event,
                              const int index,
                              const int syst,
                              const bool initialRun);
    std::pair<float, float> jet2016SFs(const float eta) const;
    std::pair<double, double> jet2017SFs(const double eta) const;

    // Histogram to be used in synchronisation.
    TH1D* synchCutFlowHist_;
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

    void getBWeight(const AnalysisEvent& event,
                    const TLorentzVector jet,
                    const int index,
                    float& mcTag,
                    float& mcNoTag,
                    float& dataTag,
                    float& dataNoTag,
                    float& err1,
                    float& err2,
                    float& err3,
                    float& err4) const;

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
    TH2D* h_muonIDs1;
    TH2F* h_muonIDs2;
    TH2D* h_muonPFiso1;
    TH2F* h_muonPFiso2;

    public:
    Cuts(const bool doPlots,
         const bool fillCutFlows,
         const bool invertLepCut,
         const bool lepCutFlow,
         const bool dumpEventNumber,
         const bool is2016,
         const bool isFCNC,
         const bool isCtag);
    ~Cuts();
    bool makeCuts(AnalysisEvent& event,
                  float& eventWeight,
                  std::map<std::string, std::shared_ptr<Plots>>& plotMap,
                  TH1D& cutFlow,
                  const int systToRun);
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
    void parse_config(const std::string confName);
    void dumpLeptonInfo(AnalysisEvent& event);
    void dumpLooseLepInfo(const AnalysisEvent& event) const;
    TH1D* getSynchCutFlow() const;
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
