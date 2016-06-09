#ifndef _cutClass_hpp_
#define _cutClass_hpp_

#include "AnalysisEvent.hpp"
#include <vector>
#include <map>
#include "plots.hpp"
#include "TH2.h"
#include <fstream>
#include "TLorentzVector.h"

#include "BTagCalibrationStandalone.hpp"

class Cuts{
  bool makeLeptonCuts(AnalysisEvent*,float*,std::map<std::string,Plots*>, TH1F*, int);
  bool invertIsoCut(AnalysisEvent*,float*,std::map<std::string,Plots*>, TH1F*);
  std::vector<int> makeJetCuts(AnalysisEvent*,int,float*);
  std::vector<int> makeMetCuts(AnalysisEvent*);
  std::vector<int> makeBCuts(AnalysisEvent*, std::vector<int>);
  std::vector<int> makeCCuts(AnalysisEvent*, std::vector<int>);
  
  std::vector<int> getTightEles(AnalysisEvent* event);
  std::vector<int> getInvIsoEles(AnalysisEvent* event);
  std::vector<int> getLooseEles(AnalysisEvent* event);
  std::vector<int> getTightMuons(AnalysisEvent* event);
  std::vector<int> getInvIsoMuons(AnalysisEvent* event);
  std::vector<int> getLooseMuons(AnalysisEvent* event);
  float getZCand(AnalysisEvent*, std::vector<int>, std::vector<int>);
  float getDileptonZCand(AnalysisEvent*, std::vector<int>, std::vector<int>);
  float getWbosonQuarksCand(AnalysisEvent*, std::vector<int>);
  float getTopMass(AnalysisEvent*);
  int getLeadingJet(AnalysisEvent *event);
  int getLeadingBjet(AnalysisEvent*);
  bool triggerCuts(AnalysisEvent*);
  
  //Method for running the synchronisation with Jeremy.
  bool synchCuts(AnalysisEvent * event);
  int getLooseLepsNum(AnalysisEvent * event); //Mimic preselection skims
  int getLooseElecs(AnalysisEvent* event);
  int getLooseMus(AnalysisEvent* event);

  //Simple deltaR function, because the reco namespace doesn't work or something
  double deltaR(float,float,float,float);
  void dumpToFile(AnalysisEvent * event, int);

  //Function to get trigger SF
  float getTriggerSF(int syst = 0, double eta1 = 999., double eta2 = 999.);

  //Function to get lepton SF
  float getLeptonWeight(AnalysisEvent*);
  float eleSF(double, double);
  float muonSF(double, double);

  //set to true to fill in histograms/spit out other info
  bool doPlots_;
  bool fillCutFlow_; // Fill cut flows
  bool invertIsoCut_; // For background estimation
  bool synchCutFlow_; //For synch
  bool singleEventInfoDump_; //For dropping info on event for synching.
  bool makeEventDump_;
  const bool trileptonChannel_;
  const bool isFCNC_;
  const bool isCtag_;

  // Tight electron cuts
  unsigned numTightEle_;
  float tightElePt_;
  float tightEleEta_;
  float tightEled0_;
  int tightEleMissLayers_;
  bool tightEleCheckPhotonVeto_;
  float tightEleMVA0_;
  float tightEleMVA1_;
  float tightEleMVA2_;
  float tightEleRelIso_;
  
  //Loose electron cuts
  unsigned numLooseEle_;
  float looseElePt_;
  float looseEleEta_;
  float looseEleMVA0_;
  float looseEleMVA1_;
  float looseEleMVA2_;
  float looseEleRelIso_;

  //Tight muon cuts
  unsigned numTightMu_;
  float tightMuonPt_;
  float tightMuonEta_;
  float tightMuonRelIso_;

  //Loose muon cuts
  unsigned numLooseMu_;
  float looseMuonPt_;
  float looseMuonEta_;
  float looseMuonRelIso_;
  
  //z and w inv cuts
  float invZMassCut_;
  float invWMassCut_;

  //Tight jet cuts
  unsigned numJets_;
  unsigned maxJets_;
  float jetPt_;
  float jetEta_;
  int jetNConsts_;
  bool jetIDDo_;

  //B-Disc cut
  unsigned numbJets_;
  unsigned maxbJets_;
  float bDiscCut_;
  float bDiscSynchCut_;

  //C-Disc cut
  unsigned numcJets_;
  unsigned maxcJets_;
  float cVsLDiscCut_;
  float cVsBDiscCut_;
  
  //Some things that will be used for JEC uncertainties.
  std::vector<float> ptMinJEC_;
  std::vector<float> ptMaxJEC_;
  std::vector<float> etaMinJEC_;
  std::vector<float> etaMaxJEC_;
  std::vector<std::vector <float> > jecSFUp_;
  std::vector<std::vector <float> > jecSFDown_;
  void initialiseJECCors();
  float getJECUncertainty(float,float,int);
  TLorentzVector getJetLVec(AnalysisEvent*,int,int);

  //Histogram to be used in synchronisation.
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
  //Sets whether to do MC or data cuts. Set every time a new dataset is processed in the main loop.
  bool isMC_;
  std::string triggerFlag_;
  std::string postfixName_;

  //For producing post-lepsel skims
  TTree* postLepSelTree_;
  TTree* postLepSelTree2_;
  TTree* postLepSelTree3_;

  //For removing trigger cuts. Will be set to false by default
  bool skipTrigger_;

  //For making b-tagging efficiencies. Needed for reweighting and systematics.
  bool makeBTagEffPlots_;
  //And the efficiency plots.
  std::vector<TH2D*> bTagEffPlots_;
  bool getBTagWeight_;
  //bTag callibration for SFs
  BTagCalibration calib;
  BTagCalibrationReader lightReader;
  BTagCalibrationReader lightReader_up;
  BTagCalibrationReader lightReader_do;
  BTagCalibrationReader charmReader;
  BTagCalibrationReader charmReader_up;
  BTagCalibrationReader charmReader_do;
  BTagCalibrationReader beautyReader;
  BTagCalibrationReader beautyReader_up;
  BTagCalibrationReader beautyReader_do;

  void getBWeight(AnalysisEvent *, TLorentzVector, int, float*, float*, float*, float*, float*, float*, float*, float*);

  //met and mtw cut values
  float metCut_;
  float mTWCut_;
  float mTWCutSynch_;

  // top mass cut values
  float TopMassCutLower_;
  float TopMassCutUpper_;

  //Sets trigger from config file
  std::string cutConfTrigLabel_;

 public:
  Cuts(bool, bool, bool, bool, bool, const bool, const bool, const bool);
  ~Cuts();
  bool makeCuts(AnalysisEvent*,float*,std::map<std::string,Plots*>, TH1F*,int);
  void setTightEle(float pt = 20, float eta = 2.5, float d0 = 0.04);
  void setMC(bool isMC) {isMC_ = isMC;}
  void setCloneTree(TTree* tree, TTree* tree2, TTree* tree3) {postLepSelTree_ = tree; postLepSelTree2_ = tree2; postLepSelTree3_ = tree3;}
  void setNumLeps(int tightMu, int looseMu, int tightEle, int looseEle){numTightEle_ = tightEle; numLooseEle_ = looseEle; numTightMu_ = tightMu; numLooseMu_ = looseMu;}
  void setCutConfTrigLabel(std::string newLabel){cutConfTrigLabel_ = newLabel;}
  void setInvIsoCut(bool invIso){invertIsoCut_ = invIso;}
  void setTriggerFlag(std::string triggerFlag) {triggerFlag_ = triggerFlag;}
  void setBTagPlots(std::vector<TH2D*> vec, bool makePlotsOrRead) {makeBTagEffPlots_ = makePlotsOrRead; bTagEffPlots_ = vec;getBTagWeight_ = !makePlotsOrRead;}
  void setSkipTrig(bool skip){skipTrigger_ = skip;}
  void setMetCut(float cut){metCut_ = cut;}
  void setMTWCut(float cut){mTWCut_ = cut;}
  void setJetRegion(int nJets, int nBets, int maxJets, int maxBJets){numJets_ = nJets; numbJets_ = nBets; maxJets_ = maxJets; maxbJets_ = maxBJets;}
  bool parse_config(std::string);
  void dumpLeptonInfo(AnalysisEvent*);
  void dumpLooseLepInfo(AnalysisEvent*);
  TH1F* getSynchCutFlow();
  int numFound(){return synchCutFlowHist_->GetBinContent(4);}
  void setEventInfoFlag(bool flag){singleEventInfoDump_ = flag;}

  private:

  TFile* electronSFsFile;
  TH2F* h_eleSFs;
  TFile* muonIDsFile;
  TFile* muonIsoFile;  
  TH2F* h_muonIDs;
  TH2F* h_muonPFiso;

};

#endif
