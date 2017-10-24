#include "TTree.h"
#include "TMVA/Timer.h"

#include "triggerScaleFactorsAlgo.hpp"
#include "config_parser.hpp"
#include "AnalysisEvent.hpp"

#include <boost/program_options.hpp>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>

#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"

const bool HIP_ERA (false);
const bool DO_HIPS (false);

Double_t ptBins[] = { 0, 10, 15, 18, 22, 24, 26, 30, 40, 50, 60, 80, 120, 500 };
Int_t numPt_bins = {13};
Double_t etaBins[] = { -2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4 };
Int_t numEta_bins = {13};

//Float_t sf_ptBins[]{ 15, 20, 25, 30, 40, 120, 200 };
//numPt_bins{6};
//Float_t sf_etaBins[]{ -2.4, -1.2, 0.0, 1.2, 2.4 };
//numEta_bins{4};

TriggerScaleFactors::TriggerScaleFactors():
  //For efficiencies
  numberPassedElectrons(),
  numberTriggeredDoubleElectrons(),

  numberPassedMuons(),
  numberTriggeredDoubleMuons(),

  numberPassedMuonElectrons(),
  numberTriggeredMuonElectrons(),

  // alpha systematics

  numberSelectedElectrons(),
  numberSelectedMuons(),
  numberSelectedMuonElectrons(),

  numberSelectedDoubleElectronsTriggered(),
  numberSelectedDoubleMuonsTriggered(),
  numberSelectedMuonElectronsTriggered()
{
  //// Plots for turn on curve studies

  // MC histos

  p_electron1_pT_MC = new TProfile("electron1_pT_MC","p_{T} turn-on curve for leading electron MC", numPt_bins, ptBins);
  p_electron1_eta_MC = new TProfile("electron1_eta_MC","#eta turn-on curve for leading MC", numEta_bins, etaBins);
  p_electron2_pT_MC = new TProfile("electron2_pT_MC","p_{T} turn-on curve for subleading electron MC", numPt_bins, ptBins);
  p_electron2_eta_MC = new TProfile("electron2_eta_MC","#eta turn-on curve for subleading MC", numEta_bins, etaBins);

  p_muon1_pT_MC = new TProfile("muon1_pT_MC","p_{T} turn-on curve for leading #mu MC", numPt_bins, ptBins);
  p_muon1_eta_MC = new TProfile("muon1_eta_MC","#eta turn-on curve for leading #mu MC", numEta_bins, etaBins);
  p_muon2_pT_MC = new TProfile("muon2_pT_MC","p_{T} turn-on curve for subleading #mu MC", numPt_bins, ptBins);
  p_muon2_eta_MC = new TProfile("muon2_eta_MC","#eta turn-on curve for subleading #mu MC", numEta_bins, etaBins);

  p_muonElectron1_pT_MC = new TProfile("muonElectron1_pT_MC","p_{T} turn-on curve for leading lep e#mu MC", numPt_bins, ptBins);
  p_muonElectron1_eta_MC = new TProfile("muonElectron1_eta_MC","#eta turn-on curve for leading lep e#mu MC", numEta_bins, etaBins);
  p_muonElectron2_pT_MC = new TProfile("muonElectron2_pT_MC","p_{T} turn-on curve for subleading lep e#mu MC", numPt_bins, ptBins);
  p_muonElectron2_eta_MC = new TProfile("muonElectron2_eta_MC","#eta turn-on curve for subleading lep e#mu MC", numEta_bins, etaBins);

  p_electrons_pT_MC = new TProfile2D("p_electrons_pT_MC","", numPt_bins, ptBins, numPt_bins, ptBins);
  p_electrons_eta_MC = new TProfile2D("p_electrons_eta_MC","", numEta_bins, etaBins, numEta_bins, etaBins);
  p_muons_pT_MC = new TProfile2D("p_muons_pT_MC","", numPt_bins, ptBins, numPt_bins, ptBins);
  p_muons_eta_MC = new TProfile2D("p_muons_eta_MC","", numEta_bins, etaBins, numEta_bins, etaBins);
  p_muonElectrons_pT_MC = new TProfile2D("p_muonElectrons_pT_MC","", numPt_bins, ptBins, numPt_bins, ptBins);
  p_muonElectrons_eta_MC = new TProfile2D("p_muonElectrons_eta_MC","", numEta_bins, etaBins, numEta_bins, etaBins);

  // Data histos

  p_electron1_pT_data = new TProfile("electron1_pT_data","p_{T} turn-on curve for leading electron data", numPt_bins, ptBins);
  p_electron1_eta_data = new TProfile("electron1_eta_data","#eta turn-on curve for leading electron data", numEta_bins, etaBins);
  p_electron2_pT_data = new TProfile("electron2_pT_data","p_{T} turn-on curve for subleading electron data", numPt_bins, ptBins);
  p_electron2_eta_data = new TProfile("electron2_eta_data","#eta turn-on curve for subleading electron data", numEta_bins, etaBins);

  p_muon1_pT_data = new TProfile("muon1_pT_data","p_{T} turn-on curve for leading #mu data", numPt_bins, ptBins);
  p_muon1_eta_data = new TProfile("muon1_eta_data","#eta turn-on curve for leading #mu data", numEta_bins, etaBins);
  p_muon2_pT_data = new TProfile("muon2_pT_data","p_{T} turn-on curve for subleading #mu data", numPt_bins, ptBins);
  p_muon2_eta_data = new TProfile("muon2_eta_data","#eta turn-on curve for subleading #mu data", numEta_bins, etaBins);

  p_muonElectron1_pT_data = new TProfile("muonElectron1_pT_data","p_{T} turn-on curve for leading lep e#mu data", numPt_bins, ptBins);
  p_muonElectron1_eta_data = new TProfile("muonElectron1_eta_data","#eta turn-on curve for leading lep e#mu data", numEta_bins, etaBins);
  p_muonElectron2_pT_data = new TProfile("muonElectron2_pT_data","p_{T} turn-on curve for subleading lep e#mu data", numPt_bins, ptBins);
  p_muonElectron2_eta_data = new TProfile("muonElectron2_eta_data","#eta turn-on curve for subleading lep e#mu data", numEta_bins, etaBins);

  p_electrons_pT_data = new TProfile2D("p_electrons_pT_data","", numPt_bins, ptBins, numPt_bins, ptBins);
  p_electrons_eta_data = new TProfile2D("p_electrons_eta_data","", numEta_bins, etaBins, numEta_bins, etaBins);
  p_muons_pT_data = new TProfile2D("p_muons_pT_data","", numPt_bins, ptBins, numPt_bins, ptBins);
  p_muons_eta_data = new TProfile2D("p_muons_eta_data","", numEta_bins, etaBins, numEta_bins, etaBins);
  p_muonElectrons_pT_data = new TProfile2D("p_muonElectrons_pT_data","", numPt_bins, ptBins, numPt_bins, ptBins);
  p_muonElectrons_eta_data = new TProfile2D("p_muonElectrons_eta_data","", numEta_bins, etaBins, numEta_bins, etaBins);

}

TriggerScaleFactors::~TriggerScaleFactors(){}

//This method is here to set up a load of branches in the TTrees that I will be analysing. Because it's vastly quicker to not load the whole damned thing.
void TriggerScaleFactors::setBranchStatusAll(TTree * chain, bool isMC, std::string triggerFlag){
  //Get electron branches
  chain->SetBranchStatus("numElePF2PAT",1);
  chain->SetBranchStatus("elePF2PATPT",1);
  chain->SetBranchStatus("elePF2PATPX",1);
  chain->SetBranchStatus("elePF2PATPY",1);
  chain->SetBranchStatus("elePF2PATPZ",1);
  chain->SetBranchStatus("elePF2PATE",1);
  chain->SetBranchStatus("elePF2PATIsGsf",1);
  chain->SetBranchStatus("elePF2PATGsfPx",1);
  chain->SetBranchStatus("elePF2PATGsfPy",1);
  chain->SetBranchStatus("elePF2PATGsfPz",1);
  chain->SetBranchStatus("elePF2PATGsfE",1);
  chain->SetBranchStatus("elePF2PATEta",1);
  chain->SetBranchStatus("elePF2PATPhi",1);
  chain->SetBranchStatus("elePF2PATBeamSpotCorrectedTrackD0",1);
  chain->SetBranchStatus("elePF2PATMissingInnerLayers",1);
  chain->SetBranchStatus("elePF2PATPhotonConversionVeto",1);
  chain->SetBranchStatus("elePF2PATMVA",1);
  chain->SetBranchStatus("elePF2PATComRelIsoRho",1);
  chain->SetBranchStatus("elePF2PATComRelIsodBeta",1);
  chain->SetBranchStatus("elePF2PATComRelIso",1);
  chain->SetBranchStatus("elePF2PATChHadIso",1);
  chain->SetBranchStatus("elePF2PATNtHadIso",1);
  chain->SetBranchStatus("elePF2PATGammaIso",1);
  chain->SetBranchStatus("elePF2PATRhoIso",1);
  chain->SetBranchStatus("elePF2PATAEff03",1);
  chain->SetBranchStatus("elePF2PATCharge",1);
  chain->SetBranchStatus("elePF2PATTrackD0",1);
  chain->SetBranchStatus("elePF2PATTrackDBD0",1);
  chain->SetBranchStatus("elePF2PATD0PV",1);
  chain->SetBranchStatus("elePF2PATBeamSpotCorrectedTrackD0",1);
  chain->SetBranchStatus("elePF2PATSCEta",1);
  //get muon branches
  chain->SetBranchStatus("muonPF2PATIsPFMuon",1);
  chain->SetBranchStatus("muonPF2PATGlobalID",1);
  chain->SetBranchStatus("muonPF2PATTrackID",1);
  chain->SetBranchStatus("numMuonPF2PAT",1);
  chain->SetBranchStatus("muonPF2PATPt",1);
  chain->SetBranchStatus("muonPF2PATPX",1);
  chain->SetBranchStatus("muonPF2PATPY",1);
  chain->SetBranchStatus("muonPF2PATPZ",1);
  chain->SetBranchStatus("muonPF2PATE",1);
  chain->SetBranchStatus("muonPF2PATEta",1);
  chain->SetBranchStatus("muonPF2PATPhi",1);
  chain->SetBranchStatus("muonPF2PATCharge",1);
  chain->SetBranchStatus("muonPF2PATComRelIsodBeta",1);
  chain->SetBranchStatus("muonPF2PATTrackDBD0",1);
  chain->SetBranchStatus("muonPF2PATD0",1);
  chain->SetBranchStatus("muonPF2PATDBInnerTrackD0",1);
  chain->SetBranchStatus("muonPF2PATTrackDBD0",1);
  chain->SetBranchStatus("muonPF2PATBeamSpotCorrectedD0",1);
  chain->SetBranchStatus("muonPF2PATD0",1);
  chain->SetBranchStatus("muonPF2PATChi2",1);
  chain->SetBranchStatus("muonPF2PATNDOF",1);
  chain->SetBranchStatus("muonPF2PATVertX",1);
  chain->SetBranchStatus("muonPF2PATVertY",1);
  chain->SetBranchStatus("muonPF2PATVertZ",1);
  chain->SetBranchStatus("muonPF2PATNChambers",1);
  chain->SetBranchStatus("muonPF2PATTrackNHits",1);
  chain->SetBranchStatus("muonPF2PATMuonNHits",1);
  chain->SetBranchStatus("muonPF2PATTkLysWithMeasurements",1);
  chain->SetBranchStatus("muonPF2PATGlbTkNormChi2",1);
  chain->SetBranchStatus("muonPF2PATDBPV",1);
  chain->SetBranchStatus("muonPF2PATDZPV",1);
  chain->SetBranchStatus("muonPF2PATVldPixHits",1);
  chain->SetBranchStatus("muonPF2PATMatchedStations",1);
  chain->SetBranchStatus("muonPF2PATGlbTkNormChi2",1);
  if (is2016_)
  {
    chain->SetBranchStatus("muonPF2PATValidFraction",1);
    chain->SetBranchStatus("muonPF2PATChi2LocalPosition",1);
    chain->SetBranchStatus("muonPF2PATTrkKick",1);
    chain->SetBranchStatus("muonPF2PATSegmentCompatibility",1);
  }
  //MET variables - for plotting (no cuts on these)
  chain->SetBranchStatus("metPF2PATEt",1);
  chain->SetBranchStatus("metPF2PATPt",1);
  //primary vertex info. For muon cut
  chain->SetBranchStatus("pvX",1);
  chain->SetBranchStatus("pvY",1);
  chain->SetBranchStatus("pvZ",1);
  //Event info
  chain->SetBranchStatus("eventNum",1);
  chain->SetBranchStatus("eventRun",1);
  chain->SetBranchStatus("eventLumiblock",1);

  if ( !is2016_ ) {
    chain->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2",1);
    chain->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3",1);
    chain->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1",1);
    chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v2",1);
    chain->SetBranchStatus("HLT_PFMET170_JetIdCleaned_v2",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v2",1);
    chain->SetBranchStatus("HLT_PFHT350_PFMET100_v1",1);
    chain->SetBranchStatus("HLT_PFHT800_v2",1);
    chain->SetBranchStatus("HLT_MET250_v1",1);
    chain->SetBranchStatus("HLT_PFHT750_4JetPt50_v3",1);
  }
  else {
    chain->SetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf_v1",1);
    chain->SetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf_v2",1);
    chain->SetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf_v3",1);
    chain->SetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf_v4",1);
    chain->SetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf_v5",1);
    chain->SetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf_v6",1);
    chain->SetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf_v7",1);
    chain->SetBranchStatus("HLT_Ele27_WPTight_Gsf_v1",1);
    chain->SetBranchStatus("HLT_Ele27_WPTight_Gsf_v2",1);
    chain->SetBranchStatus("HLT_Ele27_WPTight_Gsf_v3",1);
    chain->SetBranchStatus("HLT_Ele27_WPTight_Gsf_v4",1);
    chain->SetBranchStatus("HLT_Ele27_WPTight_Gsf_v5",1);
    chain->SetBranchStatus("HLT_Ele27_WPTight_Gsf_v6",1);
    chain->SetBranchStatus("HLT_Ele27_WPTight_Gsf_v7",1);
    chain->SetBranchStatus("HLT_Ele32_eta2p1_WPTight_Gsf_v2",1);
    chain->SetBranchStatus("HLT_Ele32_eta2p1_WPTight_Gsf_v3",1);
    chain->SetBranchStatus("HLT_Ele32_eta2p1_WPTight_Gsf_v4",1);
    chain->SetBranchStatus("HLT_Ele32_eta2p1_WPTight_Gsf_v5",1);
    chain->SetBranchStatus("HLT_Ele32_eta2p1_WPTight_Gsf_v6",1);
    chain->SetBranchStatus("HLT_Ele32_eta2p1_WPTight_Gsf_v7",1);
    chain->SetBranchStatus("HLT_Ele32_eta2p1_WPTight_Gsf_v8",1);
    chain->SetBranchStatus("HLT_IsoMu24_v1",1);
    chain->SetBranchStatus("HLT_IsoMu24_v2",1);
    chain->SetBranchStatus("HLT_IsoMu24_v3",1);
    chain->SetBranchStatus("HLT_IsoMu24_v4",1);
    chain->SetBranchStatus("HLT_IsoTkMu24_v1",1);
    chain->SetBranchStatus("HLT_IsoTkMu24_v2",1);
    chain->SetBranchStatus("HLT_IsoTkMu24_v3",1);
    chain->SetBranchStatus("HLT_IsoTkMu24_v4",1);
    chain->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3",1);
    chain->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4",1);
    chain->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5",1);
    chain->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6",1);
    chain->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7",1);
    chain->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8",1);
    chain->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3",1);
    chain->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4",1);
    chain->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3",1);
    chain->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4",1);
    chain->SetBranchStatus("HLT_MET200_v1",1);
    chain->SetBranchStatus("HLT_MET200_v2",1);
    chain->SetBranchStatus("HLT_MET200_v3",1);
    chain->SetBranchStatus("HLT_MET200_v4",1);
    chain->SetBranchStatus("HLT_MET200_v5",1);
    chain->SetBranchStatus("HLT_MET250_v1",1);
    chain->SetBranchStatus("HLT_MET250_v2",1);
    chain->SetBranchStatus("HLT_MET250_v3",1);
    chain->SetBranchStatus("HLT_MET250_v4",1);
    chain->SetBranchStatus("HLT_MET250_v5",1);
    chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v2",1);
    chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v3",1);
    chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v4",1);
    chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v5",1);
    chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v6",1);
    chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v7",1);
    chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v8",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v2",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v3",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v4",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v5",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v6",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v7",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v8",1);
    chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v9",1);
    chain->SetBranchStatus("HLT_PFHT800_v2",1);
    chain->SetBranchStatus("HLT_PFHT800_v3",1);
    chain->SetBranchStatus("HLT_PFHT800_v4",1);
    chain->SetBranchStatus("HLT_PFHT800_v5",1);
    chain->SetBranchStatus("HLT_PFHT900_v4",1);
    chain->SetBranchStatus("HLT_PFHT900_v5",1);
    chain->SetBranchStatus("HLT_PFHT900_v6",1);
    chain->SetBranchStatus("HLT_PFHT750_4JetPt50_v3",1);
    chain->SetBranchStatus("HLT_PFHT750_4JetPt50_v4",1);
    chain->SetBranchStatus("HLT_PFHT750_4JetPt50_v5",1);
    chain->SetBranchStatus("HLT_PFHT750_4JetPt50_v6",1);
    chain->SetBranchStatus("HLT_PFHT750_4JetPt70_v1",1);
    chain->SetBranchStatus("HLT_PFHT750_4JetPt70_v2",1);
    chain->SetBranchStatus("HLT_PFHT750_4JetPt80_v2",1);
    chain->SetBranchStatus("HLT_PFHT300_PFMET100_v1",1);
    chain->SetBranchStatus("HLT_PFHT300_PFMET100_v2",1);
    chain->SetBranchStatus("HLT_PFHT300_PFMET100_v3",1);
    chain->SetBranchStatus("HLT_PFHT300_PFMET100_v4",1);
    chain->SetBranchStatus("HLT_PFHT300_PFMET110_v4",1);
    chain->SetBranchStatus("HLT_PFHT300_PFMET110_v5",1);
    chain->SetBranchStatus("HLT_PFHT300_PFMET110_v6",1);
  }
}

void TriggerScaleFactors::parseCommandLineArguements(int argc, char* argv[])
{
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Print this message.")
    ("config,c", po::value<std::string>(&config)->required(),
     "The configuration file to be used.")
    (",n", po::value<long>(&nEvents)->default_value(0),
     "The number of events to be run over. All if set to 0.")
    ("outFolder,o", po::value<std::string>(&outFolder)->default_value("plots/scaleFactors/"),
     "The output directory for the plots. Overrides the config file.")
    ("postfix,s", po::value<std::string>(&postfix)->default_value("default"),
     "Set postfix for plots. Overrides the config file.")
    ("2016", po::bool_switch(&is2016_), "Use 2016 conditions (SFs, et al.)")
    ("nFiles,f", po::value<int>(&numFiles)->default_value(-1),
     "Number of files to run over. All if set to -1.");
  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
      std::cout << desc;
      std::exit(0);
    }

    po::notify(vm);
  }
  catch (const po::error& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    std::cerr << desc;
    std::exit(1);
  }

  gErrorIgnoreLevel = kInfo;

  //Set up environment a little.
  std::cout << std::setprecision(5) << std::fixed;

  //Some vectors that will be filled in the parsing.
  totalLumi = 0;
  lumiPtr = &totalLumi;
  if (!Parser::parse_config(config,&datasets,lumiPtr)){
    std::cerr << "There was an error parsing the config file.\n";
    exit(0);
  }
}

void TriggerScaleFactors::runMainAnalysis(){

  //PU reweighting

  if ( !is2016_ ) {
    //Make pileupReweighting stuff here
    dataPileupFile = new TFile("pileup/2015/truePileupTest.root","READ");
    dataPU = (TH1F*)(dataPileupFile->Get("pileup")->Clone());
    mcPileupFile = new TFile("pileup/2015/pileupMC.root","READ");
    mcPU = (TH1F*)(mcPileupFile->Get("pileup")->Clone());

    //Get systematic files too.
    systUpFile = new TFile("pileup/2015/truePileupUp.root","READ");
    pileupUpHist = (TH1F*)(systUpFile->Get("pileup")->Clone());
    systDownFile = new TFile("pileup/2015/truePileupDown.root","READ");
    pileupDownHist = (TH1F*)(systDownFile->Get("pileup")->Clone());
  }

  else {
    //Make pileupReweighting stuff here
    dataPileupFile = new TFile("pileup/2016/truePileupTest.root","READ");
    dataPU = (TH1F*)(dataPileupFile->Get("pileup")->Clone());
    mcPileupFile = new TFile("pileup/2016/pileupMC.root","READ");
    mcPU = (TH1F*)(mcPileupFile->Get("pileup")->Clone());

    //Get systematic files too.
    systUpFile = new TFile("pileup/2016/truePileupUp.root","READ");
    pileupUpHist = (TH1F*)(systUpFile->Get("pileup")->Clone());
    systDownFile = new TFile("pileup/2016/truePileupDown.root","READ");
    pileupDownHist = (TH1F*)(systDownFile->Get("pileup")->Clone());
  }

  puReweight = (TH1F*)(dataPU->Clone());
  puReweight->Scale(1.0/puReweight->Integral());
  mcPU->Scale(1.0/mcPU->Integral());
  puReweight->Divide(mcPU);
  puReweight->SetDirectory(nullptr);

  /// And do the same for systematic sampl
  puSystUp = (TH1F*)(pileupUpHist->Clone());
  puSystUp->Scale(1.0/puSystUp->Integral());
  puSystUp->Divide(mcPU);
  puSystUp->SetDirectory(nullptr);
  puSystDown = (TH1F*)(pileupDownHist->Clone());
  puSystDown->Scale(1.0/puSystDown->Integral());
  puSystDown->Divide(mcPU);
  puSystDown->SetDirectory(nullptr);

  dataPileupFile->Close();
  mcPileupFile->Close();
  systUpFile->Close();
  systDownFile->Close();

  bool datasetFilled = false;

  if (totalLumi == 0.) totalLumi = usePreLumi;
  std::cout << "Using lumi: " << totalLumi << std::endl;
  for (auto dataset = datasets.begin(); dataset!=datasets.end(); ++dataset){
    datasetFilled = false;
    TChain * datasetChain = new TChain(dataset->treeName().c_str());

    std::cerr << "Processing dataset " << dataset->name() << std::endl;
    if (!datasetFilled){
      if (!dataset->fillChain(datasetChain,numFiles)){
        std::cerr << "There was a problem constructing the chain for " << dataset->name() << ". Continuing with next dataset.\n";
	continue;
      }
      datasetFilled = true;
    }

    else{
      std::string inputPostfix{};
      inputPostfix += postfix;
      std::cout << "/scratch/data/TopPhysics/miniSkims2015/"+dataset->name()+inputPostfix + "SmallSkim.root" << std::endl;
      datasetChain->Add(("/scratch/data/TopPhysics/miniSkims2015/"+dataset->name()+inputPostfix + "SmallSkim.root").c_str());
      std::ifstream secondTree{"/scratch/data/TopPhysics/miniSkims2015/"+dataset->name()+inputPostfix + "SmallSkim1.root"};
      if (secondTree.good()) datasetChain->Add(("/scratch/data/TopPhysics/miniSkims2015/"+dataset->name()+inputPostfix + "SmallSkim1.root").c_str());
      std::ifstream thirdTree{"/scratch/data/TopPhysics/miniSkims2015/"+dataset->name()+inputPostfix + "SmallSkim2.root"};
      if (thirdTree.good()) datasetChain->Add(("/scratch/data/TopPhysics/miniSkims2015/"+dataset->name()+inputPostfix + "SmallSkim2.root").c_str());

    }

    std::cout << "Trigger flag: " << dataset->getTriggerFlag() << std::endl;

    AnalysisEvent * event = new AnalysisEvent(dataset->isMC(),dataset->getTriggerFlag(),datasetChain, is2016_, true);

    double eventWeight = 1.0;

    double pileupWeight = puReweight->GetBinContent(puReweight->GetXaxis()->FindBin(event->numVert));
    std::cout << "pileupWeight: " << pileupWeight << std::endl;
    if ( dataset->isMC() ) eventWeight *= pileupWeight;
    std::cout << "eventWeight: " << eventWeight << std::endl;

    int numberOfEvents = datasetChain->GetEntries();
    if (nEvents && nEvents < numberOfEvents) numberOfEvents = nEvents;
    auto  lEventTimer = new TMVA::Timer (numberOfEvents, "Running over dataset ...", false);
    lEventTimer->DrawProgressBar(0, "");
    for (int i = 0; i < numberOfEvents; i++) {
      lEventTimer->DrawProgressBar(i);
      event->GetEntry(i);

//      std::cout << __LINE__ << " : " << __FILE__ << std::endl;

      if ( HIP_ERA && event->eventRun >= 278820 && !(dataset->isMC()) && DO_HIPS ) continue;
      if ( !HIP_ERA && event->eventRun < 278820 && !(dataset->isMC()) && DO_HIPS ) continue;

      if ( !metFilters(event, dataset->isMC()) ) continue;

      //Does this event pass tight electron cut?
      //Create electron index
      event->electronIndexTight = getTightElectrons( event );
      bool passDoubleElectronSelection ( passDileptonSelection( event, 2 ) );
      //Does this event pass tight muon cut?
      //Create muon index
      event->muonIndexTight = getTightMuons( event );
      bool passDoubleMuonSelection ( passDileptonSelection( event, 0 ) );

      bool passMuonElectronSelection ( passDileptonSelection( event , 1 ) );

      //Triggering stuff
      int triggerDoubleEG (0), triggerDoubleMuon (0), triggerMuonElectron (0); // Passes Double Lepton Trigger
      int triggerMetDoubleEG (0), triggerMetDoubleMuon (0), triggerMetMuonElectron (0); // Passes Double Lepton and MET triggers

      //Passes event selection and MET triggers
      int triggerMetElectronSelection (0), triggerMetMuonSelection (0), triggerMetMuonElectronSelection (0); // Passes lepton selection and MET triggers

      //Does event pass Single/Double EG trigger and the electron selection?
      if ( passDoubleElectronSelection ) triggerDoubleEG 	= ( doubleElectronTriggerCut( event, dataset->isMC() ) );
      if ( passDoubleElectronSelection ) triggerMetDoubleEG 	= ( doubleElectronTriggerCut( event, dataset->isMC() )*metTriggerCut( event ) );
      //Does event pass Single/Double Muon trigger and the muon selection?
      if ( passDoubleMuonSelection )  triggerDoubleMuon = ( doubleMuonTriggerCut( event, dataset->isMC() ) );
      if ( passDoubleMuonSelection )  triggerMetDoubleMuon = ( doubleMuonTriggerCut( event, dataset->isMC() )*metTriggerCut( event ) );
      //Does event pass Single Electron/Single Muon/MuonEG trigger and the muon selection?
      if ( passMuonElectronSelection )  triggerMuonElectron = ( muonElectronTriggerCut( event, dataset->isMC() ) );
      if ( passMuonElectronSelection )  triggerMetMuonElectron = ( muonElectronTriggerCut( event, dataset->isMC() )*metTriggerCut( event ) );

      //
      //Does event pass either double lepton seletion and the MET triggers?
      if ( passDoubleElectronSelection ) triggerMetElectronSelection = ( metTriggerCut( event ) );
      if ( passDoubleMuonSelection ) triggerMetMuonSelection = ( metTriggerCut( event ) );
      if ( passMuonElectronSelection )  triggerMetMuonElectronSelection = ( metTriggerCut( event ) );

      if ( dataset->isMC() ) { // If is MC
	//SFs bit
	numberPassedElectrons[0] += triggerMetElectronSelection*eventWeight; //Number of electrons passing the cross trigger and electron selection
	numberTriggeredDoubleElectrons[0] += triggerMetDoubleEG*eventWeight; //Number of electrons passing both cross trigger+electron selection AND double EG trigger
	numberPassedMuons[0] += triggerMetMuonSelection*eventWeight; //Number of muons passing the cross trigger and muon selection
	numberTriggeredDoubleMuons[0] += triggerMetDoubleMuon*eventWeight; //Number of muons passing both cross trigger+muon selection AND double muon trigger
        numberPassedMuonElectrons[0] += triggerMetMuonElectronSelection*eventWeight; //Number of muonEGs passing cross trigger and muonEG selection
        numberTriggeredMuonElectrons[0] += triggerMetMuonElectron*eventWeight; //Number muonEGs passing both cross trigger+muonEG selection AND muonEG trigger

	// Systematic stuff
	numberSelectedElectrons[0] += passDoubleElectronSelection*eventWeight;
	numberSelectedMuons[0] += passDoubleMuonSelection*eventWeight;
	numberSelectedMuonElectrons[0] += passMuonElectronSelection*eventWeight;

	numberSelectedDoubleElectronsTriggered[0] += triggerDoubleEG*eventWeight;;
	numberSelectedDoubleMuonsTriggered[0] += triggerDoubleMuon*eventWeight;
	numberSelectedMuonElectronsTriggered[0] += triggerMuonElectron*eventWeight;

	//Histos bit
        if ( triggerMetElectronSelection > 0 ) { //If passed event selection, then will want to add to denominator
	  p_electron1_pT_MC->Fill( event->zPairLeptons.first.Pt(), triggerMetDoubleEG/triggerMetElectronSelection );
	  p_electron1_eta_MC->Fill( event->zPairLeptons.first.Eta(), triggerMetDoubleEG/triggerMetElectronSelection );
	  p_electron2_pT_MC->Fill( event->zPairLeptons.second.Pt(), triggerMetDoubleEG/triggerMetElectronSelection );
	  p_electron2_eta_MC->Fill( event->zPairLeptons.second.Eta(), triggerMetDoubleEG/triggerMetElectronSelection );

	  p_electrons_pT_MC->Fill( event->zPairLeptons.first.Pt(), event->zPairLeptons.second.Pt(), triggerMetDoubleEG/triggerMetElectronSelection );
	  p_electrons_eta_MC->Fill( event->zPairLeptons.first.Eta(), event->zPairLeptons.second.Eta(), triggerMetDoubleEG/triggerMetElectronSelection );
        }
        if ( triggerMetMuonSelection > 0 ) { //If passed event selection, then will want to add to denominator
          p_muon1_pT_MC->Fill( event->zPairLeptons.first.Pt(), triggerMetDoubleMuon/triggerMetMuonSelection );
          p_muon1_eta_MC->Fill( event->zPairLeptons.first.Eta(), triggerMetDoubleMuon/triggerMetMuonSelection );
          p_muon2_pT_MC->Fill( event->zPairLeptons.second.Pt(), triggerMetDoubleMuon/triggerMetMuonSelection );
          p_muon2_eta_MC->Fill( event->zPairLeptons.second.Eta(), triggerMetDoubleMuon/triggerMetMuonSelection );

	  p_muons_pT_MC->Fill( event->zPairLeptons.first.Pt(), event->zPairLeptons.second.Pt(), triggerMetDoubleMuon/triggerMetMuonSelection );
	  p_muons_eta_MC->Fill( event->zPairLeptons.first.Eta(), event->zPairLeptons.second.Eta(), triggerMetDoubleMuon/triggerMetMuonSelection );
        }
        if ( triggerMetMuonElectronSelection > 0 ) { //If passed event selection, then will want to add to denominator
	  p_muonElectron1_pT_MC->Fill( event->zPairLeptons.first.Pt(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
	  p_muonElectron1_eta_MC->Fill( event->zPairLeptons.first.Eta(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
	  p_muonElectron2_pT_MC->Fill( event->zPairLeptons.second.Pt(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
	  p_muonElectron2_eta_MC->Fill( event->zPairLeptons.second.Eta(), triggerMetMuonElectron/triggerMetMuonElectronSelection );

	  p_muonElectrons_pT_MC->Fill( event->zPairLeptons.first.Pt(), event->zPairLeptons.second.Pt(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
	  p_muonElectrons_eta_MC->Fill( event->zPairLeptons.first.Eta(), event->zPairLeptons.second.Eta(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
        }
      }
      else { // Else is data
	//SFs bit
	numberPassedElectrons[1] += triggerMetElectronSelection*eventWeight; //Number of electrons passing the cross trigger and electron selection
	numberTriggeredDoubleElectrons[1] += triggerMetDoubleEG*eventWeight; //Number of electrons passing both cross trigger+electron selection AND double EG trigger
	numberPassedMuons[1] += triggerMetMuonSelection*eventWeight; //Number of muons passing the cross trigger and muon selection
	numberTriggeredDoubleMuons[1] += triggerMetDoubleMuon*eventWeight; //Number of muons passing both cross trigger+muon selection AND double muon trigger
        numberPassedMuonElectrons[1] += triggerMetMuonElectronSelection*eventWeight; //Number of muonEGs passing cross trigger and muonEG selection
        numberTriggeredMuonElectrons[1] += triggerMetMuonElectron*eventWeight; //Number muonEGs passing both cross trigger+muonEG selection AND muonEG trigger

	// NB No systematic stuff required for data

	//Histos bit
        if ( triggerMetElectronSelection > 0 ) { //If passed event selection, then will want to add to denominator
	  p_electron1_pT_data->Fill( event->zPairLeptons.first.Pt(), triggerMetDoubleEG/triggerMetElectronSelection );
	  p_electron1_eta_data->Fill( event->zPairLeptons.first.Eta(), triggerMetDoubleEG/triggerMetElectronSelection );
	  p_electron2_pT_data->Fill( event->zPairLeptons.second.Pt(), triggerMetDoubleEG/triggerMetElectronSelection );
	  p_electron2_eta_data->Fill( event->zPairLeptons.second.Eta(), triggerMetDoubleEG/triggerMetElectronSelection );

	  p_electrons_pT_data->Fill( event->zPairLeptons.first.Pt(), event->zPairLeptons.second.Pt(), triggerMetDoubleEG/triggerMetElectronSelection );
	  p_electrons_eta_data->Fill( event->zPairLeptons.first.Eta(), event->zPairLeptons.second.Eta(), triggerMetDoubleEG/triggerMetElectronSelection );
        }
        if ( triggerMetMuonSelection > 0 ) { //If passed event selection, then will want to add to denominator
          p_muon1_pT_data->Fill( event->zPairLeptons.first.Pt(), triggerMetDoubleMuon/triggerMetMuonSelection);
          p_muon1_eta_data->Fill( event->zPairLeptons.first.Eta(), triggerMetDoubleMuon/triggerMetMuonSelection);
          p_muon2_pT_data->Fill( event->zPairLeptons.second.Pt(), triggerMetDoubleMuon/triggerMetMuonSelection);
          p_muon2_eta_data->Fill( event->zPairLeptons.second.Eta(), triggerMetDoubleMuon/triggerMetMuonSelection);

	  p_muons_pT_data->Fill( event->zPairLeptons.first.Pt(), event->zPairLeptons.second.Pt(), triggerMetDoubleMuon/triggerMetMuonSelection );
	  p_muons_eta_data->Fill( event->zPairLeptons.first.Eta(), event->zPairLeptons.second.Eta(), triggerMetDoubleMuon/triggerMetMuonSelection );
        }
        if ( triggerMetMuonElectronSelection > 0 ) { //If passed event selection, then will want to add to denominator
	  p_muonElectron1_pT_data->Fill( event->zPairLeptons.first.Pt(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
	  p_muonElectron1_eta_data->Fill( event->zPairLeptons.first.Eta(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
	  p_muonElectron2_pT_data->Fill( event->zPairLeptons.second.Pt(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
	  p_muonElectron2_eta_data->Fill( event->zPairLeptons.second.Eta(), triggerMetMuonElectron/triggerMetMuonElectronSelection );

	  p_muonElectrons_pT_data->Fill( event->zPairLeptons.first.Pt(), event->zPairLeptons.second.Pt(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
	  p_muonElectrons_eta_data->Fill( event->zPairLeptons.first.Eta(), event->zPairLeptons.second.Eta(), triggerMetMuonElectron/triggerMetMuonElectronSelection );
        }
      }

    }

    delete datasetChain;
  } //end dataset loop
}

std::vector<int> TriggerScaleFactors::getTightElectrons(AnalysisEvent* event) {
  std::vector<int> electrons;
  for (int i{0}; i < event->numElePF2PAT; i++){
    if (!event->elePF2PATIsGsf[i]) continue;
    TLorentzVector tempVec{event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]};

    if ( electrons.size() < 1 && tempVec.Pt() <= 20. && !is2016_ ) continue;
    else if ( electrons.size() >= 1 && tempVec.Pt() <= 15. && !is2016_ ) continue;

    if ( electrons.size() < 1 && tempVec.Pt() <= 25. && is2016_ ) continue;
    else if ( electrons.size() >= 1 && tempVec.Pt() <= 20. && is2016_ ) continue;

    if ( std::abs(event->elePF2PATSCEta[i]) > 2.50 ) continue;

    // 2015 cuts
    if ( !is2016_ ) {
      if ( (std::abs(event->elePF2PATSCEta[i]) > 1.4442 && std::abs(event->elePF2PATSCEta[i]) < 1.566) || std::abs(event->elePF2PATSCEta[i]) > 2.50 ) continue;

      // VID cut
      if ( event->elePF2PATCutIdTight[i] < 1 ) continue;
    }

    // 2016 cuts
    else {
      if ( (std::abs(event->elePF2PATSCEta[i]) > 1.4442 && std::abs(event->elePF2PATSCEta[i]) < 1.566) || std::abs(event->elePF2PATSCEta[i]) > 2.50 ) continue;

      // VID cut
      if ( event->elePF2PATCutIdTight[i] < 1 ) continue;

      // Cuts not part of the tuned ID
      if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){
        if ( std::abs(event->elePF2PATD0PV[i]) >= 0.05  ) continue;
        if ( std::abs(event->elePF2PATDZPV[i]) >= 0.10  ) continue;
      }
      else if ( std::abs(event->elePF2PATSCEta[i]) > 1.479 && std::abs(event->elePF2PATSCEta[i]) < 2.50 ){
        if ( std::abs(event->elePF2PATD0PV[i]) >= 0.10 ) continue;
        if ( std::abs(event->elePF2PATDZPV[i]) >= 0.20 ) continue;
      }
    }
    electrons.emplace_back(i);
  }

  return electrons;
}

std::vector<int> TriggerScaleFactors::getTightMuons(AnalysisEvent* event) {
  std::vector<int> muons;

  for (int i = 0; i < event->numMuonPF2PAT; i++){

    if (!event->muonPF2PATIsPFMuon[i]) continue;

    if ( muons.size() < 1 && event->muonPF2PATPt[i] <= 20. && !is2016_ ) continue;
    else if ( muons.size() >= 1 && event->muonPF2PATPt[i] <= 20. && !is2016_) continue;

    if ( muons.size() < 1 && event->muonPF2PATPt[i] <= 30. && is2016_ ) continue;
    else if ( muons.size() >= 1 && event->muonPF2PATPt[i] <= 30. && is2016_) continue;

    if (std::abs(event->muonPF2PATEta[i]) >= 2.40) continue;
    if (event->muonPF2PATComRelIsodBeta[i] >= 0.15) continue;

//    if ( !is2016_ ) {
      //Do a little test of muon id stuff here.
      if (!event->muonPF2PATTrackID[i]) continue;
      if (!event->muonPF2PATGlobalID[i]) continue;

      if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] >= 10.) continue;
      if (event->muonPF2PATTkLysWithMeasurements[i] <= 5) continue;
      if (std::abs(event->muonPF2PATDBPV[i]) >= 0.2) continue;
      if (std::abs(event->muonPF2PATDZPV[i]) >= 0.5) continue;
      if (event->muonPF2PATMuonNHits[i] < 1) continue;
      if (event->muonPF2PATVldPixHits[i] < 1) continue;
      if (event->muonPF2PATMatchedStations[i] < 2) continue;
//    }
    // 2016 cuts
    else {

      if (!event->muonPF2PATTrackID[i]) continue;
      if (!event->muonPF2PATGlobalID[i]) continue;

      if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] >= 10.) continue;
      if (event->muonPF2PATMatchedStations[i] < 2) continue;
      if (std::abs(event->muonPF2PATDBPV[i]) >= 0.2) continue;
      if (std::abs(event->muonPF2PATDZPV[i]) >= 0.5) continue;
      if (event->muonPF2PATMuonNHits[i] < 1) continue;
      if (event->muonPF2PATVldPixHits[i] < 1) continue;
      if (event->muonPF2PATTkLysWithMeasurements[i] <= 5) continue;

    }

  muons.emplace_back(i);
  }
  return muons;
}

bool TriggerScaleFactors::passDileptonSelection( AnalysisEvent *event, int nElectrons ){

  //Check if there are at least two electrons first. Otherwise use muons.

  float invMass (0.0);
  float pT (0.0);

  //DoubleEG

  if (nElectrons == 2){
    std::vector<int> leptons = event->electronIndexTight;
    for ( unsigned i = 0; i < leptons.size(); i++ ){
      for ( unsigned j = i + 1; j < leptons.size(); j++ ){
        if (event->elePF2PATCharge[leptons[i]] * event->elePF2PATCharge[leptons[j]] >= 0) continue; // check electron pair have correct charge.
        TLorentzVector lepton1 = TLorentzVector(event->elePF2PATGsfPx[leptons[i]],event->elePF2PATGsfPy[leptons[i]],event->elePF2PATGsfPz[leptons[i]],event->elePF2PATGsfE[leptons[i]]);
        TLorentzVector lepton2 = TLorentzVector(event->elePF2PATGsfPx[leptons[j]],event->elePF2PATGsfPy[leptons[j]],event->elePF2PATGsfPz[leptons[j]],event->elePF2PATGsfE[leptons[j]]);
        float candidateMass  = (lepton1 + lepton2).M();
	if ( std::abs( (lepton1 + lepton2).Pt() ) > std::abs(pT) ){
        	event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
        	event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? leptons[i]:leptons[j];
        	event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
        	event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? leptons[j]:leptons[i];
		invMass = candidateMass;
		pT = (lepton1 + lepton2).Pt();
      		}
	}
    }
  }

  // DoubleMuon
  else if (nElectrons == 0){
    std::vector<int> leptons = event->muonIndexTight;
    for ( unsigned i = 0; i < leptons.size(); i++ ){
      for ( unsigned j = i + 1; j < leptons.size(); j++ ){
	if (event->muonPF2PATCharge[leptons[i]] * event->muonPF2PATCharge[leptons[j]] >= 0) continue; 
	TLorentzVector lepton1 = TLorentzVector(event->muonPF2PATPX[leptons[i]],event->muonPF2PATPY[leptons[i]],event->muonPF2PATPZ[leptons[i]],event->muonPF2PATE[leptons[i]]);
	TLorentzVector lepton2 = TLorentzVector(event->muonPF2PATPX[leptons[j]],event->muonPF2PATPY[leptons[j]],event->muonPF2PATPZ[leptons[j]],event->muonPF2PATE[leptons[j]]);
	float candidateMass = (lepton1 + lepton2).M();
	if ( std::abs( (lepton1 + lepton2).Pt() ) > std::abs(pT) ){
		event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
		event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? leptons[i]:leptons[j];
		event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
		event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? leptons[j]:leptons[i];
		invMass = candidateMass;
		pT = (lepton1 + lepton2).Pt();
		}
      }
    }
  }

  // MuonEG
  else if (nElectrons == 1){
    std::vector<int> electrons = event->electronIndexTight;
    std::vector<int> muons = event->muonIndexTight;
    for ( unsigned i = 0; i < electrons.size(); i++ ){
      for ( unsigned j = 0; j < muons.size(); j++ ){
        if ( !(event->elePF2PATCharge[electrons[i]] * event->muonPF2PATCharge[muons[j]] >= 0) ) continue; // check muon-electron pair have correct (same) charge.
        TLorentzVector lepton1 = TLorentzVector(event->elePF2PATGsfPx[electrons[i]],event->elePF2PATGsfPy[electrons[i]],event->elePF2PATGsfPz[electrons[i]],event->elePF2PATGsfE[electrons[i]]);
	TLorentzVector lepton2 = TLorentzVector(event->muonPF2PATPX[muons[j]],event->muonPF2PATPY[muons[j]],event->muonPF2PATPZ[muons[j]],event->muonPF2PATE[muons[j]]);
	float candidateMass = (lepton1 + lepton2).M();
	if ( std::abs( (lepton1 + lepton2).Pt() ) > std::abs(pT) ){
        	event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
        	event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? electrons[i]:muons[j];
        	event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
        	event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? muons[j]:electrons[i];
		invMass = candidateMass;
		pT = (lepton1 + lepton2).Pt();
		}
      }
    }
  }

  else {
    std::cout << "Only dilepton searches currently supported. Exiting ..." << std::endl;
    exit(888);
  }

  if ( invMass > 20.0 ) return true;
  else return false;
}

bool TriggerScaleFactors::doubleElectronTriggerCut( AnalysisEvent* event, bool isMC ) {
  bool eTrig{false};

  if ( is2016_ )  { // Single lepton paths not implemented for 2015
    if ( !isMC ) {
      if ( event->HLT_Ele32_eta2p1_WPTight_Gsf_v2 > 0 ) eTrig = true;
      if ( event->HLT_Ele32_eta2p1_WPTight_Gsf_v3 > 0 ) eTrig = true;
      if ( event->HLT_Ele32_eta2p1_WPTight_Gsf_v4 > 0 ) eTrig = true;
      if ( event->HLT_Ele32_eta2p1_WPTight_Gsf_v5 > 0 ) eTrig = true;
      if ( event->HLT_Ele32_eta2p1_WPTight_Gsf_v6 > 0 ) eTrig = true;
      if ( event->HLT_Ele32_eta2p1_WPTight_Gsf_v7 > 0 ) eTrig = true;
      if ( event->HLT_Ele32_eta2p1_WPTight_Gsf_v8 > 0 ) eTrig = true;
    }
    else {
      if ( event->HLT_Ele32_eta2p1_WPTight_Gsf_v8 > 0 ) eTrig = true;
    }
  }

  bool eeTrig{false};

  if ( !is2016_ ) {
    if ( !isMC ) {
      if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) eeTrig = true;
      if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 ) eeTrig = true;
      if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) eeTrig = true;
    }
    else {
      if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) eeTrig = true;
    }
  }
  else {
    if ( !isMC ) {
      if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) eeTrig = true;
      if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) eeTrig = true;
      if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5 > 0 ) eeTrig = true;
      if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6 > 0 ) eeTrig = true;
      if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7 > 0 ) eeTrig = true;
      if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8 > 0 ) eeTrig = true;
      if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9 > 0 ) eeTrig = true;
    }
    else {
      if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9 > 0 ) eeTrig = true; 
   }
  }

  if ( eeTrig == true || eTrig == true ) return true;
  else return false;
}


bool TriggerScaleFactors::muonElectronTriggerCut( AnalysisEvent* event, bool isMC ) {
  bool muEGTrig{false};
  if ( !is2016_ ) {
    if ( !isMC ) {
      if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
    }
    else {
      if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
    }
  }
  else {
    if ( !isMC ) {
      if ( event->eventRun < 280919 ) { // Mu8 leg disabled and non-DZ versions prescaled for Run2016H 
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9 > 0 ) muEGTrig = true;
      }
      if ( event->eventRun >= 280919 ) {
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) muEGTrig = true;
	if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) muEGTrig = true;
        if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) muEGTrig = true;
      }
    }
    else { 
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) muEGTrig = true;
    }
  }
  if ( muEGTrig == true ) return true;
  else return false;
}

bool TriggerScaleFactors::doubleMuonTriggerCut( AnalysisEvent* event, bool isMC ) {
  bool muTrig{false};

  if ( is2016_ ) { // Single lepton paths not implemented for 2015
    if ( !isMC ) {

      if ( event->HLT_IsoMu24_v1 > 0 ) muTrig = true;
      if ( event->HLT_IsoMu24_v2 > 0 ) muTrig = true;
      if ( event->HLT_IsoMu24_v3 > 0 ) muTrig = true;
      if ( event->HLT_IsoMu24_v4 > 0 ) muTrig = true;
      if ( event->HLT_IsoTkMu24_v1 > 0 ) muTrig = true;
      if ( event->HLT_IsoTkMu24_v2 > 0 ) muTrig = true;
      if ( event->HLT_IsoTkMu24_v3 > 0 ) muTrig = true;
      if ( event->HLT_IsoTkMu24_v4 > 0 ) muTrig = true;

//Run B
//      if ( event->HLT_IsoMu24_v1 > 0 && event->eventRun >= 272007 && event->eventRun < 275657 ) muTrig = true; // RunB
//      if ( event->HLT_IsoMu24_v2 > 0 && event->eventRun >= 272007 && event->eventRun < 275657 ) muTrig = true; // RunB
//      if ( event->HLT_IsoTkMu24_v1 > 0 && event->eventRun >= 272007 && event->eventRun < 275657 ) muTrig = true; // RunB
//      if ( event->HLT_IsoTkMu24_v2 > 0 && event->eventRun >= 272007 && event->eventRun < 275657 ) muTrig = true; // RunB

//Run C
//      if ( event->HLT_IsoMu24_v2 > 0 && event->eventRun >= 275657 && event->eventRun < 276315 ) muTrig = true; // RunC
//      if ( event->HLT_IsoTkMu24_v3 > 0 && event->eventRun >= 275657 && event->eventRun < 276315 ) muTrig = true; // RunC

//Run D
//      if ( event->HLT_IsoMu24_v2 > 0 && event->eventRun >= 276315 && event->eventRun < 276831 ) muTrig = true; // RunD
//      if ( event->HLT_IsoTkMu24_v3 > 0 && event->eventRun >= 276315 && event->eventRun < 276831 ) muTrig = true; // RunD

//Run E
//      if ( event->HLT_IsoMu24_v2 > 0 && event->eventRun >= 276831 && event->eventRun < 277772 ) muTrig = true; // RunE
//      if ( event->HLT_IsoTkMu24_v3 > 0 && event->eventRun >= 276831 && event->eventRun < 277772 ) muTrig = true; // RunE

//Run F
//      if ( event->HLT_IsoMu24_v2 > 0 && event->eventRun >= 277772 && event->eventRun < 278820 ) muTrig = true; // RunF
//      if ( event->HLT_IsoTkMu24_v3 > 0 && event->eventRun >= 277772 && event->eventRun < 278820 ) muTrig = true; // RunF

//Run G
//      if ( event->HLT_IsoMu24_v2 > 0 && event->eventRun >= 278820 && event->eventRun < 280919 ) muTrig = true; // RunG
//      if ( event->HLT_IsoTkMu24_v3 > 0 && event->eventRun >= 278820 && event->eventRun < 280919 ) muTrig = true; // RunG

//Run H
//      if ( event->HLT_IsoMu24_v4 > 0 && event->eventRun >= 280919 ) muTrig = true; // RunH
//     if ( event->HLT_IsoTkMu24_v4 > 0 && event->eventRun >= 280919 ) muTrig = true; // RunH


    }
    else {
      if ( event->HLT_IsoMu24_v4 > 0 ) muTrig = true;
      if ( event->HLT_IsoTkMu24_v4 > 0 ) muTrig = true;
    }
  }

  bool mumuTrig{false};

  if ( !is2016_ ) {
    if ( !isMC ) {
      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1 > 0 ) mumuTrig = true;
      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1 > 0 ) mumuTrig = true;
      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
    }
    else {
      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
    }
  }
  else {
    if ( !isMC ) {

      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 && HIP_ERA ) mumuTrig = true; //pre-HIP
      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 > 0 && HIP_ERA ) mumuTrig = true; //pre-HIP
      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 ) mumuTrig = true; //pre-HIP & post-HIP
      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 > 0 && !HIP_ERA ) mumuTrig = true; //post-HIP
      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 && HIP_ERA ) mumuTrig = true; //pre-HIP
      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 ) mumuTrig = true; //pre-HIP & post-HIP
      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6 > 0 && !HIP_ERA ) mumuTrig = true; //post-HIP
    

//Run B
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 && event->eventRun >= 272007 && event->eventRun < 275657 ) mumuTrig = true; // RunB
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 > 0 && event->eventRun >= 272007 && event->eventRun < 275657 ) mumuTrig = true; // RunB
//      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 && event->eventRun >= 272007 && event->eventRun < 275657 ) mumuTrig = true; // RunB
//      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 && event->eventRun >= 272007 && event->eventRun < 275657 ) mumuTrig = true; // RunB

//Run C
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 > 0 && event->eventRun >= 275657 && event->eventRun < 276315 ) mumuTrig = true; // RunC
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 && event->eventRun >= 275657 && event->eventRun < 276315 ) mumuTrig = true; // RunC
//      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 && event->eventRun >= 275657 && event->eventRun < 276315 ) mumuTrig = true; // RunC

//Run D
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 && event->eventRun >= 276315 && event->eventRun < 276831 ) mumuTrig = true; // RunD
//      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 && event->eventRun >= 276315 && event->eventRun < 276831 ) mumuTrig = true; // RunD

//Run E
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 && event->eventRun >= 276831 && event->eventRun < 277772 ) mumuTrig = true; // RunE
//      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 && event->eventRun >= 276831 && event->eventRun < 277772 ) mumuTrig = true; // RunE

//Run F
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 && event->eventRun >= 277772 && event->eventRun < 278820 ) mumuTrig = true; // RunF
//      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 && event->eventRun >= 277772 && event->eventRun < 278820 ) mumuTrig = true; // RunF

//Run G
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 && event->eventRun >= 278820 && event->eventRun < 280919 ) mumuTrig = true; // RunG
//      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 && event->eventRun >= 278820 && event->eventRun < 280919 ) mumuTrig = true; // RunG

//Run H
//      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 > 0 && event->eventRun >= 280919 ) mumuTrig = true; // RunH
 //     if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6 > 0 && event->eventRun >= 280919 ) mumuTrig = true; // RunH

    }
    else {
      if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 > 0 ) mumuTrig = true;
      if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6 > 0 ) mumuTrig = true; //post-HIP
    }
  }

  if ( mumuTrig == true || muTrig == true ) return true;
  else return false;
}

bool TriggerScaleFactors::metTriggerCut( AnalysisEvent* event ) {
  bool metTrig {false}; 
  if ( !is2016_ ) {
    if ( event->HLT_PFMET120_PFMHT120_IDTight_v2 > 0 ) metTrig = true;
    if ( event->HLT_PFMET170_JetIdCleaned_v2 > 0 ) metTrig = true;
    if ( event->HLT_PFMET170_HBHECleaned_v2 > 0 ) metTrig = true;
    if ( event->HLT_PFHT800_v2 > 0 ) metTrig = true;
    if ( event->HLT_MET250_v1 > 0 ) metTrig = true;
  }
  else {
    if ( event->HLT_MET200_v1 > 0 ) metTrig = true;
    if ( event->HLT_MET200_v2 > 0 ) metTrig = true;
    if ( event->HLT_MET200_v3 > 0 ) metTrig = true;
    if ( event->HLT_MET200_v4 > 0 ) metTrig = true;
    if ( event->HLT_MET200_v5 > 0 ) metTrig = true;
    if ( event->HLT_MET250_v1 > 0 ) metTrig = true;
    if ( event->HLT_MET250_v2 > 0 ) metTrig = true;
    if ( event->HLT_MET250_v3 > 0 ) metTrig = true;
    if ( event->HLT_MET250_v4 > 0 ) metTrig = true;
    if ( event->HLT_MET250_v5 > 0 ) metTrig = true;
    if ( event->HLT_PFMET120_PFMHT120_IDTight_v2 > 0 ) metTrig = true;
    if ( event->HLT_PFMET120_PFMHT120_IDTight_v3 > 0 ) metTrig = true;
    if ( event->HLT_PFMET120_PFMHT120_IDTight_v4 > 0 )metTrig = true;
    if ( event->HLT_PFMET120_PFMHT120_IDTight_v5 > 0 )metTrig = true;
    if ( event->HLT_PFMET120_PFMHT120_IDTight_v6 > 0 )metTrig = true;
    if ( event->HLT_PFMET120_PFMHT120_IDTight_v7 > 0 )metTrig = true;
    if ( event->HLT_PFMET120_PFMHT120_IDTight_v8 > 0 )metTrig = true;
    if ( event->HLT_PFMET170_HBHECleaned_v2 > 0 ) metTrig = true;
    if ( event->HLT_PFMET170_HBHECleaned_v3 > 0 ) metTrig = true;
    if ( event->HLT_PFMET170_HBHECleaned_v4 > 0 ) metTrig = true; 
    if ( event->HLT_PFMET170_HBHECleaned_v5 > 0 ) metTrig = true;
    if ( event->HLT_PFMET170_HBHECleaned_v6 > 0 ) metTrig = true;
    if ( event->HLT_PFMET170_HBHECleaned_v7 > 0 ) metTrig = true;
    if ( event->HLT_PFMET170_HBHECleaned_v8 > 0 ) metTrig = true;
    if ( event->HLT_PFMET170_HBHECleaned_v9 > 0 ) metTrig = true;

//
// EXCLUDED DUE TO BEING UNABLE TO VERIFY THESE PATHS WERE NOT PRESCALED
//
/*    if ( event->HLT_PFHT800_v2 > 0 ) metTrig = true;
    if ( event->HLT_PFHT800_v3 > 0 ) metTrig = true;
    if ( event->HLT_PFHT800_v4 > 0 ) metTrig = true;
    if ( event->HLT_PFHT800_v5 > 0 ) metTrig = true;
    if ( event->HLT_PFHT900_v4 > 0 ) metTrig = true;
    if ( event->HLT_PFHT900_v5 > 0 ) metTrig = true;
    if ( event->HLT_PFHT900_v6 > 0 ) metTrig = true;
    if ( event->HLT_PFHT750_4JetPt50_v3 > 0 ) metTrig = true;
    if ( event->HLT_PFHT750_4JetPt50_v4 > 0 ) metTrig = true;
    if ( event->HLT_PFHT750_4JetPt50_v5 > 0 ) metTrig = true;
    if ( event->HLT_PFHT750_4JetPt50_v6 > 0 ) metTrig = true;
    if ( event->HLT_PFHT750_4JetPt70_v1 > 0 ) metTrig = true;
    if ( event->HLT_PFHT750_4JetPt70_v2 > 0 ) metTrig = true;
    if ( event->HLT_PFHT750_4JetPt80_v2 > 0 ) metTrig = true;
*/
//
//
//

    if ( event->HLT_PFHT300_PFMET100_v1 > 0 ) metTrig = true;
    if ( event->HLT_PFHT300_PFMET100_v2 > 0 ) metTrig = true;
    if ( event->HLT_PFHT300_PFMET100_v3 > 0 ) metTrig = true;
    if ( event->HLT_PFHT300_PFMET100_v4 > 0 ) metTrig = true;
    if ( event->HLT_PFHT300_PFMET110_v4 > 0 ) metTrig = true;
    if ( event->HLT_PFHT300_PFMET110_v5 > 0 ) metTrig = true;
    if ( event->HLT_PFHT300_PFMET110_v6 > 0 ) metTrig = true;
  }
  if ( metTrig == true ) return true;
  else return false;
}

bool TriggerScaleFactors::metFilters(AnalysisEvent* event, bool isMC) {
  if ( event->Flag_HBHENoiseFilter <= 0 ) return false;
  if ( event->Flag_HBHENoiseIsoFilter <= 0 ) return false;
  if ( event->Flag_EcalDeadCellTriggerPrimitiveFilter <= 0 ) return false;
  if ( event->Flag_goodVertices <= 0 ) return false;
  if ( event->Flag_eeBadScFilter <= 0 ) return false;
  if ( !is2016_ && event->Flag_CSCTightHalo2015Filter <= 0 ) return false;
  if ( is2016_ ) {
    if ( event->Flag_globalTightHalo2016Filter <= 0 ) return false;
    if ( event->Flag_chargedHadronTrackResolutionFilter <= 0 ) return false;
    if ( event->Flag_muonBadTrackFilter <= 0 ) return false;
    if ( event->Flag_ecalLaserCorrFilter <= 0 ) return false;
    if ( event->Flag_noBadMuons <= 0 && !isMC ) return false;
  }
  return true;
}

void TriggerScaleFactors::savePlots()
{

  double level = 0.60; //ClopperPearson interval level

  // Histos first
  TFile *outFile{new TFile{ (outFolder+"triggerPlots.root").c_str(), "RECREATE"}};
  
  // Do pT errors
  for ( Int_t bin = 1; bin != numPt_bins+1; bin++ ) {
    double errUp, errDown, error;

    // electrons MC
    //ele 1 pT MC
    errUp = ( p_electron1_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron1_pT_MC->GetBinEntries(bin), p_electron1_pT_MC->GetBinEntries(bin) * p_electron1_pT_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_electron1_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron1_pT_MC->GetBinEntries(bin), p_electron1_pT_MC->GetBinEntries(bin) * p_electron1_pT_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_electron1_pT_MC->SetBinError(bin, error);
    //ele 2 pT MC
    errUp = ( p_electron2_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron2_pT_MC->GetBinEntries(bin), p_electron2_pT_MC->GetBinEntries(bin) * p_electron2_pT_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_electron2_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron2_pT_MC->GetBinEntries(bin), p_electron2_pT_MC->GetBinEntries(bin) * p_electron2_pT_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_electron2_pT_MC->SetBinError(bin, error);

    // muons MC
    //muon 1 pT MC
    errUp = ( p_muon1_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon1_pT_MC->GetBinEntries(bin), p_muon1_pT_MC->GetBinEntries(bin) * p_muon1_pT_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_muon1_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon1_pT_MC->GetBinEntries(bin), p_muon1_pT_MC->GetBinEntries(bin) * p_muon1_pT_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muon1_pT_MC->SetBinError(bin, error);
    //muon 2 pT MC
    errUp = ( p_muon2_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon2_pT_MC->GetBinEntries(bin), p_muon2_pT_MC->GetBinEntries(bin) * p_muon2_pT_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_muon2_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon2_pT_MC->GetBinEntries(bin), p_muon2_pT_MC->GetBinEntries(bin) * p_muon2_pT_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muon2_pT_MC->SetBinError(bin, error);

    // muonEG MC
    //muonEG 1 pT MC
    errUp = ( p_muonElectron1_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron1_pT_MC->GetBinEntries(bin), p_muonElectron1_pT_MC->GetBinEntries(bin) * p_muonElectron1_pT_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_muonElectron1_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron1_pT_MC->GetBinEntries(bin), p_muonElectron1_pT_MC->GetBinEntries(bin) * p_muonElectron1_pT_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muonElectron1_pT_MC->SetBinError(bin, error);
    //muonEG 2 pT MC
    errUp = ( p_muonElectron2_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron2_pT_MC->GetBinEntries(bin), p_muonElectron2_pT_MC->GetBinEntries(bin) * p_muonElectron2_pT_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_muonElectron2_pT_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron2_pT_MC->GetBinEntries(bin), p_muonElectron2_pT_MC->GetBinEntries(bin) * p_muonElectron2_pT_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muonElectron2_pT_MC->SetBinError(bin, error);

    // electrons data
    //ele 1 pT data
    errUp = ( p_electron1_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron1_pT_data->GetBinEntries(bin), p_electron1_pT_data->GetBinEntries(bin) * p_electron1_pT_data->GetBinContent( bin ), level, true) );
    errDown = ( p_electron1_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron1_pT_data->GetBinEntries(bin), p_electron1_pT_data->GetBinEntries(bin) * p_electron1_pT_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_electron1_pT_data->SetBinError(bin, error);
    //ele 2 pT data
    errUp = ( p_electron2_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron2_pT_data->GetBinEntries(bin), p_electron2_pT_data->GetBinEntries(bin) * p_electron2_pT_data->GetBinContent( bin ), level, true) );
    errDown = ( p_electron2_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron2_pT_data->GetBinEntries(bin), p_electron2_pT_data->GetBinEntries(bin) * p_electron2_pT_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_electron2_pT_data->SetBinError(bin, error);

    // muons data
    //muon 1 pT data
    errUp = ( p_muon1_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon1_pT_data->GetBinEntries(bin), p_muon1_pT_data->GetBinEntries(bin) * p_muon1_pT_data->GetBinContent( bin ), level, true) );
    errDown = ( p_muon1_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon1_pT_data->GetBinEntries(bin), p_muon1_pT_data->GetBinEntries(bin) * p_muon1_pT_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muon1_pT_data->SetBinError(bin, error);
    //muon 2 pT data
    errUp = ( p_muon2_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon2_pT_data->GetBinEntries(bin), p_muon2_pT_data->GetBinEntries(bin) * p_muon2_pT_data->GetBinContent( bin ), level, true) );
    errDown = ( p_muon2_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon2_pT_data->GetBinEntries(bin), p_muon2_pT_data->GetBinEntries(bin) * p_muon2_pT_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muon2_pT_data->SetBinError(bin, error);

    // muonEG data
    //muonEG 1 pT data
    errUp = ( p_muonElectron1_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron1_pT_data->GetBinEntries(bin), p_muonElectron1_pT_data->GetBinEntries(bin) * p_muonElectron1_pT_data->GetBinContent( bin ), level, true) );
    errDown = ( p_muonElectron1_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron1_pT_data->GetBinEntries(bin), p_muonElectron1_pT_data->GetBinEntries(bin) * p_muonElectron1_pT_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muonElectron1_pT_data->SetBinError(bin, error);
    //muonEG 2 pT data
    errUp = ( p_muonElectron2_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron2_pT_data->GetBinEntries(bin), p_muonElectron2_pT_data->GetBinEntries(bin) * p_muonElectron2_pT_data->GetBinContent( bin ), level, true) );
    errDown = ( p_muonElectron2_pT_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron2_pT_data->GetBinEntries(bin), p_muonElectron2_pT_data->GetBinEntries(bin) * p_muonElectron2_pT_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muonElectron2_pT_data->SetBinError(bin, error);

  }

  // Do eta errors
  for ( Int_t bin = 1; bin != numEta_bins+1; bin++ ) {
    double errUp, errDown, error;

    // electrons MC
    //ele 1 eta MC
    errUp = ( p_electron1_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron1_eta_MC->GetBinEntries(bin), p_electron1_eta_MC->GetBinEntries(bin) * p_electron1_eta_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_electron1_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron1_eta_MC->GetBinEntries(bin), p_electron1_eta_MC->GetBinEntries(bin) * p_electron1_eta_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_electron1_eta_MC->SetBinError(bin, error);
    //ele 2 eta MC
    errUp = ( p_electron2_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron2_eta_MC->GetBinEntries(bin), p_electron2_eta_MC->GetBinEntries(bin) * p_electron2_eta_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_electron2_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron2_eta_MC->GetBinEntries(bin), p_electron2_eta_MC->GetBinEntries(bin) * p_electron2_eta_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_electron2_eta_MC->SetBinError(bin, error);

    // muons MC
    //muon 1 eta MC
    errUp = ( p_muon1_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon1_eta_MC->GetBinEntries(bin), p_muon1_eta_MC->GetBinEntries(bin) * p_muon1_eta_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_muon1_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon1_eta_MC->GetBinEntries(bin), p_muon1_eta_MC->GetBinEntries(bin) * p_muon1_eta_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muon1_eta_MC->SetBinError(bin, error);
    //muon 2 eta MC
    errUp = ( p_muon2_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon2_eta_MC->GetBinEntries(bin), p_muon2_eta_MC->GetBinEntries(bin) * p_muon2_eta_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_muon2_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon2_eta_MC->GetBinEntries(bin), p_muon2_eta_MC->GetBinEntries(bin) * p_muon2_eta_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muon2_eta_MC->SetBinError(bin, error);

    // muonEG MC
    //muonEG 1 eta MC
    errUp = ( p_muonElectron1_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron1_eta_MC->GetBinEntries(bin), p_muonElectron1_eta_MC->GetBinEntries(bin) * p_muonElectron1_eta_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_muonElectron1_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron1_eta_MC->GetBinEntries(bin), p_muonElectron1_eta_MC->GetBinEntries(bin) * p_muonElectron1_eta_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muonElectron1_eta_MC->SetBinError(bin, error);
    //muonEG 2 eta MC
    errUp = ( p_muonElectron2_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron2_eta_MC->GetBinEntries(bin), p_muonElectron2_eta_MC->GetBinEntries(bin) * p_muonElectron2_eta_MC->GetBinContent( bin ), level, true) );
    errDown = ( p_muonElectron2_eta_MC->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron2_eta_MC->GetBinEntries(bin), p_muonElectron2_eta_MC->GetBinEntries(bin) * p_muonElectron2_eta_MC->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muonElectron2_eta_MC->SetBinError(bin, error);

    // electrons data
    //ele 1 eta data
    errUp = ( p_electron1_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron1_eta_data->GetBinEntries(bin), p_electron1_eta_data->GetBinEntries(bin) * p_electron1_eta_data->GetBinContent( bin ), level, true) );
    errDown = ( p_electron1_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron1_eta_data->GetBinEntries(bin), p_electron1_eta_data->GetBinEntries(bin) * p_electron1_eta_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_electron1_eta_data->SetBinError(bin, error);
    //ele 2 eta data
    errUp = ( p_electron2_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron2_eta_data->GetBinEntries(bin), p_electron2_eta_data->GetBinEntries(bin) * p_electron2_eta_data->GetBinContent( bin ), level, true) );
    errDown = ( p_electron2_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_electron2_eta_data->GetBinEntries(bin), p_electron2_eta_data->GetBinEntries(bin) * p_electron2_eta_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_electron2_eta_data->SetBinError(bin, error);

    // muons data
    //muon 1 eta data
    errUp = ( p_muon1_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon1_eta_data->GetBinEntries(bin), p_muon1_eta_data->GetBinEntries(bin) * p_muon1_eta_data->GetBinContent( bin ), level, true) );
    errDown = ( p_muon1_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon1_eta_data->GetBinEntries(bin), p_muon1_eta_data->GetBinEntries(bin) * p_muon1_eta_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muon1_eta_data->SetBinError(bin, error);
    //muon 2 eta data
    errUp = ( p_muon2_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon2_eta_data->GetBinEntries(bin), p_muon2_eta_data->GetBinEntries(bin) * p_muon2_eta_data->GetBinContent( bin ), level, true) );
    errDown = ( p_muon2_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muon2_eta_data->GetBinEntries(bin), p_muon2_eta_data->GetBinEntries(bin) * p_muon2_eta_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muon2_eta_data->SetBinError(bin, error);

    // muonEG data
    //muonEG 1 eta data
    errUp = ( p_muonElectron1_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron1_eta_data->GetBinEntries(bin), p_muonElectron1_eta_data->GetBinEntries(bin) * p_muonElectron1_eta_data->GetBinContent( bin ), level, true) );
    errDown = ( p_muonElectron1_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron1_eta_data->GetBinEntries(bin), p_muonElectron1_eta_data->GetBinEntries(bin) * p_muonElectron1_eta_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muonElectron1_eta_data->SetBinError(bin, error);
    //muonEG 2 eta data
    errUp = ( p_muonElectron2_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron2_eta_data->GetBinEntries(bin), p_muonElectron2_eta_data->GetBinEntries(bin) * p_muonElectron2_eta_data->GetBinContent( bin ), level, true) );
    errDown = ( p_muonElectron2_eta_data->GetBinContent( bin ) - TEfficiency::ClopperPearson (p_muonElectron2_eta_data->GetBinEntries(bin), p_muonElectron2_eta_data->GetBinEntries(bin) * p_muonElectron2_eta_data->GetBinContent( bin ), level, false) );
    error = errUp>errDown?errUp:errDown;
    p_muonElectron2_eta_data->SetBinError(bin, error);
  }

  // SF histos

  TH2F* electronPtSF = new TH2F ("electronPtSF","ee p_{T} trigger SFs", numPt_bins, ptBins, numPt_bins, ptBins);
  TH2F* electronEtaSF = new TH2F ("electronEtaSF","ee #eta trigger SFs", numEta_bins, etaBins, numEta_bins, etaBins);

  TH2F* muonPtSF = new TH2F ("muonPtSF","#mu#mu p_{T} trigger SFs", numPt_bins, ptBins, numPt_bins, ptBins);
  TH2F* muonEtaSF = new TH2F ("muonEtaSF","#mu#mu #eta trigger SFs", numEta_bins, etaBins, numEta_bins, etaBins);

  TH2F* muonElectronPtSF = new TH2F ("muonElectronPtSF","e#mu p_{T} trigger SFs", numPt_bins, ptBins, numPt_bins, ptBins);
  TH2F* muonElectronEtaSF = new TH2F ("muonElectronEtaSF","e#mu #eta trigger SFs", numEta_bins, etaBins, numEta_bins, etaBins);

  for ( Int_t bin1 = 1; bin1 != numPt_bins+1; bin1++ ) {
    for ( Int_t bin2 = 1; bin2 != numPt_bins+1; bin2++ ) {
      double pT, pT_err;

      pT = p_electrons_pT_data->GetBinContent( bin1, bin2 )/(p_electrons_pT_MC->GetBinContent( bin1, bin2 )+1.0e-6);
      pT_err = ( p_electrons_pT_data->GetBinContent( bin1, bin2 )+p_electrons_pT_data->GetBinError( bin1, bin2 ) ) / ( p_electrons_pT_MC->GetBinContent( bin1, bin2 ) - p_electrons_pT_MC->GetBinError( bin1, bin2) );
      electronPtSF->SetBinContent( bin1, bin2, pT);
      electronPtSF->SetBinError( bin1, bin2, pT_err );

      pT = p_muons_pT_data->GetBinContent( bin1, bin2 )/(p_muons_pT_MC->GetBinContent( bin1, bin2 )+1.0e-6);
      pT_err = ( p_muons_pT_data->GetBinContent( bin1, bin2 )+p_muons_pT_data->GetBinError( bin1, bin2 ) ) / ( p_muons_pT_MC->GetBinContent( bin1, bin2 ) - p_muons_pT_MC->GetBinError( bin1, bin2) );
      muonPtSF->SetBinContent( bin1, bin2, pT);
      muonPtSF->SetBinError( bin1, bin2, pT_err );

      pT = p_muonElectrons_pT_data->GetBinContent( bin1, bin2 )/(p_muonElectrons_pT_MC->GetBinContent( bin1, bin2 )+1.0e-6);
      pT_err = ( p_muonElectrons_pT_data->GetBinContent( bin1, bin2 )+p_muonElectrons_pT_data->GetBinError( bin1, bin2 ) ) / ( p_muonElectrons_pT_MC->GetBinContent( bin1, bin2 ) - p_muonElectrons_pT_MC->GetBinError( bin1, bin2) );
      muonElectronPtSF->SetBinContent( bin1, bin2, pT);
      muonElectronPtSF->SetBinError( bin1, bin2, pT_err );
    }
  }

  for ( Int_t bin1 = 1; bin1 != numEta_bins+1; bin1++ ) {
    for ( Int_t bin2 = 1; bin2 != numEta_bins+1; bin2++ ) {
      double eta, eta_err;

      eta = p_electrons_eta_data->GetBinContent( bin1, bin2 )/(p_electrons_eta_MC->GetBinContent( bin1, bin2 )+1.0e-6);
      eta_err = ( p_electrons_eta_data->GetBinContent( bin1, bin2 )+p_electrons_eta_data->GetBinError( bin1, bin2 ) ) / ( p_electrons_eta_MC->GetBinContent( bin1, bin2 ) - p_electrons_eta_MC->GetBinError( bin1, bin2) );
      electronEtaSF->SetBinContent( bin1, bin2, eta);
      electronEtaSF->SetBinError( bin1, bin2, eta_err );

      eta = p_muons_eta_data->GetBinContent( bin1, bin2 )/(p_muons_eta_MC->GetBinContent( bin1, bin2 )+1.0e-6);
      eta_err = ( p_muons_eta_data->GetBinContent( bin1, bin2 )+p_muons_eta_data->GetBinError( bin1, bin2 ) ) / ( p_muons_eta_MC->GetBinContent( bin1, bin2 ) - p_muons_eta_MC->GetBinError( bin1, bin2) );
      muonEtaSF->SetBinContent( bin1, bin2, eta);
      muonEtaSF->SetBinError( bin1, bin2, eta_err );

      eta = p_muonElectrons_eta_data->GetBinContent( bin1, bin2 )/(p_muonElectrons_eta_MC->GetBinContent( bin1, bin2 )+1.0e-6);
      eta_err = ( p_muonElectrons_eta_data->GetBinContent( bin1, bin2 )+p_muonElectrons_eta_data->GetBinError( bin1, bin2 ) ) / ( p_muonElectrons_eta_MC->GetBinContent( bin1, bin2 ) - p_muonElectrons_eta_MC->GetBinError( bin1, bin2) );
      muonElectronEtaSF->SetBinContent( bin1, bin2, eta);
      muonElectronEtaSF->SetBinError( bin1, bin2, eta_err );
    }
  }

  // Write Histos

  p_electron1_pT_MC->Write();
  p_electron1_eta_MC->Write();
  p_electron2_pT_MC->Write();
  p_electron2_eta_MC->Write();
  p_muon1_pT_MC->Write();
  p_muon1_eta_MC->Write();
  p_muon2_pT_MC->Write();
  p_muon2_eta_MC->Write();
  p_muonElectron1_pT_MC->Write();
  p_muonElectron1_eta_MC->Write();
  p_muonElectron2_pT_MC->Write();
  p_muonElectron2_eta_MC->Write();

  p_electron1_pT_data->Write();
  p_electron1_eta_data->Write();
  p_electron2_pT_data->Write();
  p_electron2_eta_data->Write();
  p_muon1_pT_data->Write();
  p_muon1_eta_data->Write();
  p_muon2_pT_data->Write();
  p_muon2_eta_data->Write();
  p_muonElectron1_pT_data->Write();
  p_muonElectron1_eta_data->Write();
  p_muonElectron2_pT_data->Write();
  p_muonElectron2_eta_data->Write();

  p_electrons_pT_MC->Write();
  p_electrons_eta_MC->Write();
  p_muons_pT_MC->Write();
  p_muons_eta_MC->Write();
  p_muonElectrons_pT_MC->Write();
  p_muonElectrons_eta_MC->Write();

  p_electrons_pT_data->Write();
  p_electrons_eta_data->Write();
  p_muons_pT_data->Write();
  p_muons_eta_data->Write();
  p_muonElectrons_pT_data->Write();
  p_muonElectrons_eta_data->Write();

  electronPtSF->Write();
  electronEtaSF->Write();
  muonPtSF->Write();
  muonEtaSF->Write();
  muonElectronPtSF->Write();
  muonElectronEtaSF->Write();

  outFile->Close();

  // Calculate MC efficiency

  //// LeptonTriggers
  double doubleElectronEfficiencyMC 	= numberTriggeredDoubleElectrons[0]/(numberPassedElectrons[0]+1.0e-6);
  double doubleMuonEfficiencyMC 	= numberTriggeredDoubleMuons[0]/(numberPassedMuons[0]+1.0e-6);
  double muonElectronEfficiencyMC 	= numberTriggeredMuonElectrons[0]/(numberPassedMuonElectrons[0]+1.0e-6);

  // Calculate Data efficiency

  //// DoubleLeptonTriggers
  double doubleElectronEfficiencyData	= numberTriggeredDoubleElectrons[1]/(numberPassedElectrons[1]+1.0e-6);
  double doubleMuonEfficiencyData 	= numberTriggeredDoubleMuons[1]/(numberPassedMuons[1]+1.0e-6);
  double muonElectronEfficiencyData 	= numberTriggeredMuonElectrons[1]/(numberPassedMuonElectrons[1]+1.0e-6);

  // Calculate SF

  //// LeptonTriggers
  double doubleElectronSF 	= doubleElectronEfficiencyData/(doubleElectronEfficiencyMC+1.0e-6);
  double doubleMuonSF 		= doubleMuonEfficiencyData/(doubleMuonEfficiencyMC+1.0e-6);
  double muonElectronSF 	= muonElectronEfficiencyData/(muonElectronEfficiencyMC+1.0e-6);

  // Calculate alphas
  double alphaDoubleElectron    = ( (numberSelectedDoubleElectronsTriggered[0]/numberSelectedElectrons[0])*(numberPassedElectrons[0]/numberSelectedElectrons[0]) )/(numberTriggeredDoubleElectrons[0]/numberSelectedElectrons[0]+1.0e-6);
  double alphaDoubleMuon        = ( (numberSelectedDoubleMuonsTriggered[0]/numberSelectedMuons[0])*(numberPassedMuons[0]/numberSelectedMuons[0]) )/(numberTriggeredDoubleMuons[0]/numberSelectedMuons[0]+1.0e-6);
  double alphaMuonElectron 	= ( (numberSelectedMuonElectronsTriggered[0]/numberSelectedMuonElectrons[0])*(numberPassedMuonElectrons[0]/numberSelectedMuonElectrons[0]) )/(numberTriggeredMuonElectrons[0]/numberSelectedMuonElectrons[0]+1.0e-6);

  // Calculate uncertainities

  //// LeptonTriggers

  double doubleElectronDataUpperUncert = doubleElectronEfficiencyData-TEfficiency::ClopperPearson(numberPassedElectrons[1], numberTriggeredDoubleElectrons[1], level, true);
  double doubleElectronMcUpperUncert   = doubleElectronEfficiencyMC-TEfficiency::ClopperPearson(numberPassedElectrons[0], numberTriggeredDoubleElectrons[0], level, true);
  double doubleElectronDataLowerUncert = doubleElectronEfficiencyData-TEfficiency::ClopperPearson(numberPassedElectrons[1], numberTriggeredDoubleElectrons[1], level, false);
  double doubleElectronMcLowerUncert   = doubleElectronEfficiencyMC-TEfficiency::ClopperPearson(numberPassedElectrons[0], numberTriggeredDoubleElectrons[0], level, false);

  double doubleMuonDataUpperUncert = doubleMuonEfficiencyData-TEfficiency::ClopperPearson(numberPassedMuons[1], numberTriggeredDoubleMuons[1], level, true);
  double doubleMuonMcUpperUncert   = doubleMuonEfficiencyMC-TEfficiency::ClopperPearson(numberPassedMuons[0], numberTriggeredDoubleMuons[0], level, true);
  double doubleMuonDataLowerUncert = doubleMuonEfficiencyData-TEfficiency::ClopperPearson(numberPassedMuons[1], numberTriggeredDoubleMuons[1], level, false);
  double doubleMuonMcLowerUncert   = doubleMuonEfficiencyMC-TEfficiency::ClopperPearson(numberPassedMuons[0], numberTriggeredDoubleMuons[0], level, false);

  double muonElectronDataUpperUncert = muonElectronEfficiencyData-TEfficiency::ClopperPearson(numberPassedMuonElectrons[1], numberTriggeredMuonElectrons[1], level, true);
  double muonElectronMcUpperUncert   = muonElectronEfficiencyMC-TEfficiency::ClopperPearson(numberPassedMuonElectrons[0], numberTriggeredMuonElectrons[0], level, true);
  double muonElectronDataLowerUncert = muonElectronEfficiencyData-TEfficiency::ClopperPearson(numberPassedMuonElectrons[1], numberTriggeredMuonElectrons[1], level, false);
  double muonElectronMcLowerUncert   = muonElectronEfficiencyMC-TEfficiency::ClopperPearson(numberPassedMuonElectrons[0], numberTriggeredMuonElectrons[0], level, false);

  double doubleEleSfUp = (doubleElectronEfficiencyData+doubleElectronDataUpperUncert)/(doubleElectronEfficiencyMC-doubleElectronMcLowerUncert+1.0e-6) - doubleElectronSF;
  double doubleEleSfDown = (doubleElectronEfficiencyData+doubleElectronDataLowerUncert)/(doubleElectronEfficiencyMC-doubleElectronMcUpperUncert+1.0e-6) - doubleElectronSF;
  double doubleEleSfUncert = 0.0;
  if ( doubleEleSfUp > doubleEleSfDown ) doubleEleSfUncert = doubleEleSfUp;
  else doubleEleSfUncert = doubleEleSfDown;

  double doubleMuonSfUp = (doubleMuonEfficiencyData+doubleMuonDataUpperUncert)/(doubleMuonEfficiencyMC-doubleMuonMcLowerUncert+1.0e-6) - doubleMuonSF;
  double doubleMuonSfDown = (doubleMuonEfficiencyData+doubleMuonDataLowerUncert)/(doubleMuonEfficiencyMC-doubleMuonMcUpperUncert+1.0e-6) - doubleMuonSF;
  double doubleMuonSfUncert = 0.0;
  if ( doubleMuonSfUp > doubleMuonSfDown ) doubleMuonSfUncert = doubleMuonSfUp;
  else doubleMuonSfUncert = doubleMuonSfDown;

  double muonElectronSfUp = (muonElectronEfficiencyData+muonElectronDataUpperUncert)/(muonElectronEfficiencyMC-muonElectronMcLowerUncert+1.0e-6) - muonElectronSF;
  double muonElectronSfDown = (muonElectronEfficiencyData+muonElectronDataLowerUncert)/(muonElectronEfficiencyMC-muonElectronMcUpperUncert+1.0e-6) - muonElectronSF;
  double muonElectronSfUncert = 0.0;
  if ( muonElectronSfUp > muonElectronSfDown ) muonElectronSfUncert = muonElectronSfUp;
  else muonElectronSfUncert = muonElectronSfDown;

  // Print output

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Double Electron data efficiency: " << doubleElectronEfficiencyData << " +/- " << doubleElectronDataUpperUncert << "/" << doubleElectronDataLowerUncert << std::endl;
  std::cout << "Double Electron MC efficiency: " << doubleElectronEfficiencyMC << " +/- " << doubleElectronMcUpperUncert << "/" << doubleElectronMcLowerUncert << std::endl;
  std::cout << "Double Electron trigger SF: " << doubleElectronSF << " +/- " << doubleEleSfUncert << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Double Muon data efficiency: " << doubleMuonEfficiencyData << " +/- " << doubleMuonDataUpperUncert << "/" << doubleMuonDataLowerUncert << std::endl;
  std::cout << "Double Muon MC efficiency: " << doubleMuonEfficiencyMC << " +/- " << doubleMuonMcUpperUncert << "/" << doubleMuonMcLowerUncert << std::endl;
  std::cout << "Double Muon trigger SF: " << doubleMuonSF << " +/- " << doubleMuonSfUncert << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "MuonEG data efficiency: " << muonElectronEfficiencyData << " +/- " << muonElectronDataUpperUncert << "/" << muonElectronDataLowerUncert << std::endl;
  std::cout << "MuonEG MC efficiency: " << muonElectronEfficiencyMC << " +/- " << muonElectronMcUpperUncert << " / " << muonElectronMcLowerUncert << std::endl;
  std::cout << "MuonEG trigger SF: " << muonElectronSF << " +/- " << muonElectronSfUncert << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "alpha for DoubleEG/DoubleMuon/MuonEG Triggers: " << alphaDoubleElectron << "/" << alphaDoubleMuon << "/" << alphaMuonElectron << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

}

