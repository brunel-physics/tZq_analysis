#include "TTree.h"
#include "TMVA/Timer.h"

#include "debugInfo.hpp"
#include "config_parser.hpp"
#include "AnalysisEvent.hpp"

#include <boost/program_options.hpp>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>

#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TLatex.h"

int main(int argc, char* argv[]){

  DebugInfo debugInfo;

  debugInfo.parseCommandLineArguements(argc, argv);
  debugInfo.runMainAnalysis();
}

DebugInfo::DebugInfo(){}

DebugInfo::~DebugInfo(){}

void DebugInfo::parseCommandLineArguements(int argc, char* argv[])
{
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Print this message.")
    ("config,c", po::value<std::string>(&config)->required(),
     "The configuration file to be used.")
    (",n", po::value<long>(&nEvents)->default_value(0),
     "The number of events to be run over. All if set to 0.")
    ("outFolder,o", po::value<std::string>(&outFolder)->default_value("plots/debug/"),
     "The output directory for the plots. Overrides the config file.")
    ("postfix,s", po::value<std::string>(&postfix)->default_value("default"),
     "Set postfix for plots. Overrides the config file.")
    ("2016", po::bool_switch(&is2016_), "Use 2016 conditions (SFs, et al.)")
    ("MC,m", po::bool_switch(&isMC_),
     "Running over Monte Carlo.")
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
  std::cout << std::setprecision(3) << std::fixed;

  //Some vectors that will be filled in the parsing.
  totalLumi = 0;
  lumiPtr = &totalLumi;
  if (!Parser::parse_config(config,&datasets,lumiPtr)){
    std::cerr << "There was an error parsing the config file.\n";
    exit(0);
  }

  if (!is2016_) {
    std::cerr << "Not setup for debuging 2015 data.\n";
    exit(1);
  }

}

void DebugInfo::runMainAnalysis(){

  TH1F* histChannel = new TH1F("histChannel", "channel yields vs RunH/RunsB-G;channel;RunH/RunsB-G",2,0.5,2.5);
  TH1F* histTrileptonChannel = new TH1F("histTrileptonChannel", "trilepton channel yields vs RunH/RunsB-G;channel;RunH/RunsB-G",4,0.5,4.5);
  TH1F* histNumJets = new TH1F("histNumJets", "number of jets vs RunH/RunsB-G;# jets;RunH/RunsB-G",10,-0.5,9.5);
  TH1F* histNumBJets = new TH1F("histNumBJets", "number of b-jets vs RunH/RunsB-G;# b-jets;RunH/RunsB-G",5,-0.5,4.5);

  std::pair< int, int > numElectrons {0,0};
  std::pair< int, int > numMuons {0,0};
  std::pair< int, int > numEEE {0,0};
  std::pair< int, int > numEEMU {0,0};
  std::pair< int, int > numEMUMU {0,0};
  std::pair< int, int > numMUMUMU {0,0};
  std::pair< int, int > numJets[10] {{0,0}};
  std::pair< int, int > numBJets[5] {{0,0}};

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

    AnalysisEvent * event = new AnalysisEvent(dataset->isMC(),dataset->getTriggerFlag(),datasetChain, is2016_, true);

    int numberOfEvents = datasetChain->GetEntries();
    if (nEvents && nEvents < numberOfEvents) numberOfEvents = nEvents;
    auto  lEventTimer = new TMVA::Timer (numberOfEvents, "Running over dataset ...", false);
    lEventTimer->DrawProgressBar(0, "");
    for (int i = 0; i < numberOfEvents; i++) {
      lEventTimer->DrawProgressBar(i);
      event->GetEntry(i);

// MET Triggers
/*
      bool metTrig = false;
      if( event->HLT_MET250_v2 > 0 ) metTrig = true;
      if( event->HLT_MET250_v3 > 0 ) metTrig = true;
      if( event->HLT_MET250_v4 > 0 ) metTrig = true;
      if( event->HLT_MET250_v5 > 0 ) metTrig = true;

      if( event->HLT_PFHT300_PFMET100_v1 > 0 ) > 0 ) metTrig = true;
      if( event->HLT_PFHT300_PFMET100_v2 > 0 ) metTrig = true;
      if( event->HLT_PFHT300_PFMET100_v3 > 0 ) metTrig = true;
      if( event->HLT_PFHT300_PFMET100_v4 > 0 ) metTrig = true;
      if( event->HLT_PFHT300_PFMET110_v4 > 0 ) metTrig = true;
      if( event->HLT_PFHT300_PFMET110_v5 > 0 ) metTrig = true;
      if( event->HLT_PFHT300_PFMET110_v6 > 0 ) metTrig = true;

      if( event->HLT_PFMET120_PFMHT120_IDTight_v2 > 0 ) metTrig = true;
      if( event->HLT_PFMET120_PFMHT120_IDTight_v3 > 0 ) metTrig = true;
      if( event->HLT_PFMET120_PFMHT120_IDTight_v4 > 0 ) metTrig = true;
      if( event->HLT_PFMET120_PFMHT120_IDTight_v5 > 0 ) metTrig = true;
      if( event->HLT_PFMET120_PFMHT120_IDTight_v6 > 0 ) metTrig = true;
      if( event->HLT_PFMET120_PFMHT120_IDTight_v7 > 0 ) metTrig = true;
      if( event->HLT_PFMET120_PFMHT120_IDTight_v8 > 0 ) metTrig = true;
      if( event->HLT_PFMET170_HBHECleaned_v2 > 0 ) metTrig = true;
      if( event->HLT_PFMET170_HBHECleaned_v3 > 0 ) metTrig = true;
      if( event->HLT_PFMET170_HBHECleaned_v4 > 0 ) metTrig = true;
      if( event->HLT_PFMET170_HBHECleaned_v5 > 0 ) metTrig = true;
      if( event->HLT_PFMET170_HBHECleaned_v6 > 0 ) metTrig = true;
      if( event->HLT_PFMET170_HBHECleaned_v7 > 0 ) metTrig = true;
      if( event->HLT_PFMET170_HBHECleaned_v8 > 0 ) metTrig = true;
      if( event->HLT_PFMET170_HBHECleaned_v9 > 0 ) metTrig = true;

      if ( !metTrig ) continue;
*/

      if ( event->numElePF2PAT == 2 && event->numMuonPF2PAT == 0 ) {
        if ( event->elePF2PATPT[0] < 15.0 ) continue;
        if ( event->elePF2PATPT[1] < 15.0 ) continue;
        if ( event->elePF2PATCharge[0] * event->elePF2PATCharge[1] >= 0 )  continue; // check electron pair have correct charge.
        TLorentzVector lepton1{event->elePF2PATGsfPx[0],event->elePF2PATGsfPy[0],event->elePF2PATGsfPz[0],event->elePF2PATGsfE[0]};
        TLorentzVector lepton2{event->elePF2PATGsfPx[1],event->elePF2PATGsfPy[1],event->elePF2PATGsfPz[1],event->elePF2PATGsfE[1]};
        double invMass{(lepton1 + lepton2).M() -91.1};
	if (std::abs(invMass) > 30.0 ) continue;
        if (!isMC_){
          if ( event->eventRun <= 280385 ) numElectrons.first += 1; // If Runs B-G
          else numElectrons.second += 1; // else if Run H
        }
        else numElectrons.first += 1; // just MC
      }

      if ( event->numMuonPF2PAT == 2 && event->numElePF2PAT == 0 ) {
        if ( event->muonPF2PATPt[0] < 15.0 ) continue;
        if ( event->muonPF2PATPt[1] < 15.0 ) continue;
        if ( event->muonPF2PATCharge[0] * event->muonPF2PATCharge[1] >= 0 ) continue;
	TLorentzVector lepton1{event->muonPF2PATPX[0],event->muonPF2PATPY[0],event->muonPF2PATPZ[0],event->muonPF2PATE[0]};
	TLorentzVector lepton2{event->muonPF2PATPX[1],event->muonPF2PATPY[1],event->muonPF2PATPZ[1],event->muonPF2PATE[1]};
        double invMass{(lepton1 + lepton2).M() -91.1};
	if (std::abs(invMass) > 30.0 ) continue;
        if (!isMC_){
          if ( event->eventRun <= 280385 ) numMuons.first += 1; // If Runs B-G
          else numMuons.second += 1; // else if Run H
        }
        else numMuons.first += 1; // just MC
      }

      if ( event->numElePF2PAT == 3 && event->numMuonPF2PAT == 0 ) {
//        if ( event->elePF2PATPT[0] < 15.0 ) continue;
//        if ( event->elePF2PATPT[1] < 15.0 ) continue;
//        if ( event->elePF2PATPT[2] < 15.0 ) continue;
        if (!isMC_){
          if ( event->eventRun <= 280385 ) numEEE.first += 1; // If Runs B-G
          else numEEE.second += 1; // else if Run H
        }
        else numEEE.first += 1; // just MC
      }
      if ( event->numElePF2PAT == 2 && event->numMuonPF2PAT == 1 ) {
//        if ( event->elePF2PATPT[0] < 15.0 ) continue;
//        if ( event->elePF2PATPT[1] < 15.0 ) continue;
//        if ( event->muonPF2PATPt[0] < 15.0 ) continue;
//        if ( event->elePF2PATCharge[0] * event->elePF2PATCharge[1] >= 0 )  continue; // check electron pair have correct charge.
        if (!isMC_){
          if ( event->eventRun <= 280385 ) numEEMU.first += 1; // If Runs B-G
          else numEEMU.second += 1; // else if Run H
        }
        else numEEMU.first += 1; // just MC
      }
      if ( event->numElePF2PAT == 1 && event->numMuonPF2PAT == 2 ) {
//        if ( event->elePF2PATPT[0] < 15.0 ) continue;
//        if ( event->muonPF2PATPt[0] < 15.0 ) continue;
//        if ( event->muonPF2PATPt[1] < 15.0 ) continue;
//        if ( event->muonPF2PATCharge[0] * event->muonPF2PATCharge[1] >= 0 ) continue;
        if (!isMC_){
          if ( event->eventRun <= 280385 ) numEMUMU.first += 1; // If Runs B-G
          else numEMUMU.second += 1; // else if Run H
        }
        else numEMUMU.first += 1; // just MC
      }
      if ( event->numElePF2PAT == 0 && event->numMuonPF2PAT == 3 ) {
//        if ( event->muonPF2PATPt[0] < 15.0 ) continue;
//        if ( event->muonPF2PATPt[1] < 15.0 ) continue;
//        if ( event->muonPF2PATPt[2] < 15.0 ) continue;
        if (!isMC_){
          if ( event->eventRun <= 280385 ) numMUMUMU.first += 1; // If Runs B-G
          else numMUMUMU.second += 1; // else if Run H
        }
        else numMUMUMU.first += 1; // just MC
      }

      if ( event->numJetPF2PAT < 10 ) {
        if (!isMC_){
          if ( event->eventRun <= 280385) numJets[event->numJetPF2PAT].first += 1;
          else numJets[event->numJetPF2PAT].second += 1;
        }
        else numJets[event->numJetPF2PAT].first += 1; // just MC
      }
      unsigned int bJets {0};
      for (int j = 0; j< event->numJetPF2PAT; j++) {
        if ( event->jetPF2PATBDiscriminator[i] > 0.5426 ) bJets += 1;
      }
      if ( bJets < 5 ) {
        if (!isMC_){
          if ( event->eventRun <= 280385) numBJets[bJets].first += 1;
          else numBJets[bJets].second += 1;
        }
        else numBJets[bJets].first += 1; // just MC
      }
    }
    delete datasetChain;
  } //end dataset loop

//  std::cout << "numElectrons.first/numElectrons.second : " << numElectrons.first << "/" << numElectrons.second << std::endl;
//  std::cout << "numMuons.first/numMuons.second : " << numMuons.first << "/" << numMuons.second << std::endl;

  float electronFraction = float(numElectrons.second)/(float(numElectrons.first)+1.0e-06);
  float muonFraction = float(numMuons.second)/(float(numMuons.first)+1.0e-06);

  if (!isMC_) histChannel->Fill(1, electronFraction);
  else histChannel->Fill(1, numElectrons.first);
  if (!isMC_) histChannel->Fill(2, muonFraction);
  else histChannel->Fill(2, numMuons.first);

  float eeeFraction = float(numEEE.second)/(float(numEEE.first)+1.0e-06);
  float eemuFraction = float(numEEMU.second)/(float(numEEMU.first)+1.0e-06);
  float emumuFraction = float(numEMUMU.second)/(float(numEMUMU.first)+1.0e-06);
  float mumumuFraction = float(numMUMUMU.second)/(float(numMUMUMU.first)+1.0e-06);

  if (!isMC_) histTrileptonChannel->Fill(1, eeeFraction);
  else histTrileptonChannel->Fill(1, numEEE.first);
  if (!isMC_) histTrileptonChannel->Fill(2, eemuFraction);
  else histTrileptonChannel->Fill(2, numEEMU.first);
  if (!isMC_) histTrileptonChannel->Fill(3, emumuFraction);
  else histTrileptonChannel->Fill(3, numEMUMU.first);
  if (!isMC_) histTrileptonChannel->Fill(4, mumumuFraction);
  else histTrileptonChannel->Fill(4, numMUMUMU.first);

  for ( int i = 0; i < 10; i++ ) {
    if (!isMC_) histNumJets->Fill(i, numJets[i].second/(numJets[i].first+1.0e-06));
    else histNumJets->Fill(i, numJets[i].first);
  }
  for ( int i = 0; i < 5; i++ ) {
    if (!isMC_) histNumBJets->Fill(i, numBJets[i].second/(numBJets[i].first+1.0e-06));
    else histNumBJets->Fill(i, numBJets[i].first);
  }

  histChannel->GetXaxis()->SetBinLabel(1,"ee");
  histChannel->GetXaxis()->SetBinLabel(2,"#mu#mu");

  histTrileptonChannel->GetXaxis()->SetBinLabel(1,"eee");
  histTrileptonChannel->GetXaxis()->SetBinLabel(2,"ee#mu");
  histTrileptonChannel->GetXaxis()->SetBinLabel(3,"e#mu#mu");
  histTrileptonChannel->GetXaxis()->SetBinLabel(4,"#mu#mu#mu");

  histChannel->SetMinimum(0.0);
  histNumJets->SetMinimum(0.0);
  histNumBJets->SetMinimum(0.0);

  mkdir( (outFolder).c_str(),0700);
  TFile* outFile = new TFile ( (outFolder+postfix+".root").c_str(), "RECREATE" );
  histChannel->Write();
  histTrileptonChannel->Write();
  histNumJets->Write();
  histNumBJets->Write();
  outFile->Close();
}

//This method is here to set up a load of branches in the TTrees that I will be analysing. Because it's vastly quicker to not load the whole damned thing.
void DebugInfo::setBranchStatusAll(TTree * chain, bool isMC, std::string triggerFlag){
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


