#include <array>
#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <set>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include "AnalysisEvent.hpp"


int main (int argc, char* argv[])
{
  std::string dileptonFileName = argv[1];
  std::string singleLeptonFileName = argv[2];
  std::string datasetName = argv[3];
  std::string channel = argv[4];

  int singleElectron {0}, dupElectron {0};
  int singleMuon {0}, dupMuon {0};

  const std::string postTriggerSkimDir{std::string{"/scratch/data/TopPhysics/postTriggerSkims2016/"}};

  std::ifstream dileptonFileList(dileptonFileName.c_str());
  std::string dileptonLine;

  std::ifstream singleLeptonFileList(singleLeptonFileName.c_str());
  std::string singleLeptonLine;

  int fileNum{0};
  
  std::set< std::pair < Int_t, Int_t > > triggerDoubleCountCheck;
 
  while(getline(dileptonFileList,dileptonLine)){

    std::string numName;
    std::ostringstream convert;
    convert << fileNum;
    numName = convert.str();
      
    TChain datasetChain{"tree"};
    datasetChain.Add(dileptonLine.c_str());

    std::cout << "Dilepton input file: " << dileptonLine << std::endl;

    TTree * const outTree = datasetChain.CloneTree(0);
      
    std::string outFilePath{ postTriggerSkimDir + datasetName + "/triggerSkim" + numName + ".root" };
    TFile outFile{outFilePath.c_str(), "RECREATE"};

    std::cout << "Output file: " << outFilePath << std::endl;
	  
    const long long int numberOfEvents{datasetChain.GetEntries()};
    
    AnalysisEvent event{false, "", &datasetChain, true};

    for (long long int i{0}; i < numberOfEvents; i++) {
      if (i % 500 < 0.01) std::cerr << i << "/" << numberOfEvents << " (" << 100*float(i)/numberOfEvents << "%)\r";
      event.GetEntry(i);
	
      if ( channel == "ee" ) {
	bool eeTrig{false};
	if ( event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) eeTrig = true;
	if ( event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) eeTrig = true;
	if ( event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5 > 0 ) eeTrig = true;
	if ( event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6 > 0 ) eeTrig = true;
	if ( event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7 > 0 ) eeTrig = true;
	if ( event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8 > 0 ) eeTrig = true;
	if ( event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9 > 0 ) eeTrig = true;
	  
	if ( eeTrig ) {
	  triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
	  outTree->Fill();
	}
      }
	
      if ( channel == "mumu" ) {
        bool mumuTrig{false};
        if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 > 0 ) mumuTrig = true;
        if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5 > 0 ) mumuTrig = true;
	if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6 > 0 ) mumuTrig = true;
	  
	if ( mumuTrig ) {
	  triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
	  outTree->Fill();
	}
      }
	
      if ( channel == "emu" ) {
	bool muEGTrig{false};
	// Runs B-G only
	if ( event.eventRun < 280919 ) {
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9 > 0 ) muEGTrig = true;
	}
	// Run H only
	if ( event.eventRun >= 280919 ) {
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) muEGTrig = true;
	  if ( event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) muEGTrig = true;
	}
	if ( muEGTrig ) {
	  triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
	  outTree->Fill();
	}
      } // end emu	
    }
    outFile.cd();	
    outTree->Write();

    outFile.Write();
    outFile.Close();
	
    fileNum++;

    std::cerr << "" << std::endl;
  }
  
  while(getline(singleLeptonFileList,singleLeptonLine)){

    std::string numName;
    std::ostringstream convert;
    convert << fileNum;
    numName = convert.str();

    TChain datasetChain{"tree"};
    datasetChain.Add(singleLeptonLine.c_str());

    std::cout << "Single lepton input file: " << singleLeptonLine << std::endl;

    TTree * const outTree = datasetChain.CloneTree(0);
      
    std::string outFilePath{ postTriggerSkimDir + datasetName + "/triggerSkim" + numName + ".root" };
    TFile outFile{outFilePath.c_str(), "RECREATE"};

    std::cout << "Output file: " << outFilePath << std::endl;
	  
    const long long int numberOfEvents{datasetChain.GetEntries()};

    AnalysisEvent event{false, "", &datasetChain, true};

    for (long long int i{0}; i < numberOfEvents; i++) {
      event.GetEntry(i);
      if ( channel == "ee" ) {
        bool eTrig{false};
 	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v2 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v3 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v4 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v5 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v6 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v7 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v8 > 0 ) eTrig = true;

	if ( eTrig ) {
          auto it = triggerDoubleCountCheck.find( std::make_pair( event.eventRun, event.eventNum ) );
	  singleElectron++;
	  // If event has already been found ... skip event
	  if ( it != triggerDoubleCountCheck.end() ) dupElectron++;
	  // If event has not already been found, add to new skim
	  else {
            triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
	    outTree->Fill();
	  }
	}
      } // end single electron check for ee
	
      if ( channel == "mumu") {
	bool muTrig{false};
	if ( event.HLT_IsoMu24_v1 > 0 ) muTrig = true;
	if ( event.HLT_IsoMu24_v2 > 0 ) muTrig = true;
	if ( event.HLT_IsoMu24_v3 > 0 ) muTrig = true;
	if ( event.HLT_IsoMu24_v4 > 0 ) muTrig = true;
	if ( event.HLT_IsoTkMu24_v1 > 0 ) muTrig = true;
	if ( event.HLT_IsoTkMu24_v2 > 0 ) muTrig = true;
	if ( event.HLT_IsoTkMu24_v3 > 0 ) muTrig = true;
	if ( event.HLT_IsoTkMu24_v4 > 0 ) muTrig = true;

	// If single Muon triggered fired, check to see if event also fired a DoubleMuon trigger
	if ( muTrig ) {
	  singleMuon++;
          auto it = triggerDoubleCountCheck.find( std::make_pair( event.eventRun, event.eventNum ) );
	  // If event has already been found ... skip event
	  if ( it != triggerDoubleCountCheck.end() ) dupMuon++;
	  // If event has not already been found, add to new skim
          else {
            triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
	    outTree->Fill();
	  }
	}

      } // end single muon check for mumu
	
      if ( channel == "emu" ) {
        // check eTrigger for emu first
        bool eTrig{false};
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v2 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v3 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v4 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v5 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v6 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v7 > 0 ) eTrig = true;
	if ( event.HLT_Ele32_eta2p1_WPTight_Gsf_v8 > 0 ) eTrig = true;
	// then check muTrigger for emu
	bool muTrig{false};
	if ( event.HLT_IsoMu24_v1 > 0 ) muTrig = true;
	if ( event.HLT_IsoMu24_v2 > 0 ) muTrig = true;
	if ( event.HLT_IsoMu24_v3 > 0 ) muTrig = true;
	if ( event.HLT_IsoMu24_v4 > 0 ) muTrig = true;
	if ( event.HLT_IsoTkMu24_v1 > 0 ) muTrig = true;
	if ( event.HLT_IsoTkMu24_v2 > 0 ) muTrig = true;
	if ( event.HLT_IsoTkMu24_v3 > 0 ) muTrig = true;
	if ( event.HLT_IsoTkMu24_v4 > 0 ) muTrig = true;

	// If either single lepton  triggered fired, check to see if event also fired a MuonEG trigger
	if ( eTrig || muTrig ) {
          if ( eTrig ) singleElectron++;
          if ( muTrig ) singleMuon++;
          auto it = triggerDoubleCountCheck.find( std::make_pair( event.eventRun, event.eventNum ) );
	  // If event has already been found ... skip event
	  if ( it != triggerDoubleCountCheck.end() ) {
            if ( eTrig )  dupElectron++;
            if ( muTrig ) dupMuon++;
	  }
	  // If event has not already been found, add to new skim
          else {
            triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
	    outTree->Fill();
	  }
	}
      } // end single lepton check for emu
    }
    outFile.cd();	
    outTree->Write();

    outFile.Write();
    outFile.Close();

    fileNum++;
	
    std::cerr << "" << std::endl;
  }

  if ( channel == "ee" || channel == "emu" )    std::cout << "Single electron trigger fired with double lepton trigger/Total single electron triggers fired: " << dupElectron << " / " << singleElectron << std::endl;
  if ( channel == "mumu" || channel == "emu" )  std::cout << "Single muon trigger fired with double lepton trigger/Total single muon triggers fired: " << dupMuon << " / " << singleMuon << std::endl;

}
