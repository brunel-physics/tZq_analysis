#include "AnalysisEvent.hpp"
#include <libconfig.h++>
#include <dirent.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMVA/Timer.h"

#include <stdlib.h> 
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

static void show_usage(std::string name){
  std::cerr << "Usage: " << name << " <options>"
	    << "Options:\n"
	    << "\t-o  \tOUTFILE\tOutput file path for plots. If set overwrites the default.\n"
	    << "\t-i \tOUTFILE\tInput file path for input nTuples. If set overwrites the default.\n"
	    << "\t-h  --help\t\t\tShow this help message\n"
	    << std::endl;
}

int main(int argc, char* argv[]) {

  std::string inputDir = "";
  std::string outFileString = "plots/distributions/output.root";
  bool inputFlag (false);

  for (int i = 1; i < argc; ++i){
    std::string arg = argv[i];
    if ((arg=="-h") || (arg == "--help")){ // Display help stuff
      show_usage(argv[0]);
      return 0;
    }
    else if (arg=="-o"){//Set output file
      if (i + 1 < argc){
	outFileString = argv[++i];
      } 
      else{
	std::cerr << "requires a string for output folder name." << std::endl;;
      }
    }
    else if (arg=="-i"){//Set input folder
      if (i + 1 < argc){
	inputDir = argv[++i];
    if ( inputDir.back() != '/' ) inputDir += '/';
	if ( inputDir == "" ){
	  std::cerr << "requires a non-null input dir to be run over!" << std::endl;;
	  return 0;
	}
	else{
	  inputFlag = true;
	}
      }
    }
  }


  if (!inputFlag){
    std::cerr << "Program requires an input dir to be run over!" << std::endl;;
	return 0;
  }

  DIR *dp;
  struct dirent *dirp;

  if((dp  = opendir( inputDir.c_str() )) == nullptr) {
    std::cout << "Error opening Directory" << std::endl;
    std::cout << inputDir.c_str() << " is not a valid directory" << std::endl;
    return 0;
  }

  std::vector<TTree*> inputTrees;

  std::cout << "Attaching files to TTree ... " << std::endl;

  while ( (dirp = readdir(dp)) != nullptr) {
    std::string line (dirp->d_name);
    if ( line == "." || line == "..")
    continue;
    TFile *inputFile = new TFile ((inputDir+line).c_str()) ;
    TTree *lTempTree = dynamic_cast<TTree*>(inputFile->Get("tree"));
    inputTrees.push_back(lTempTree);
  }

  std::cout << "Attached all files to TTree!" << std::endl;

  auto  histElePt      = new TH1F ("histElePt"    , "Distribution of reco-electron p_{T}" , 500, 0.0  , 500.0);
  auto  histEleEta     = new TH1F ("histEleEta"   , "Distribution of reco-electron #eta"  , 500, -2.50, 2.5);
  auto  histEleGenPt   = new TH1F ("histEleGenPt" , "Distribution of gen-electron p_{T}"  , 500, 0.0  , 500.0);
  auto  histEleGenEta  = new TH1F ("histEleGenEta", "Distribution of gen-electron #eta"   , 500, -2.5 , 2.5);

  auto  histMuPt       = new TH1F ("histMuPt"     , "Distribution of reco-muon p_{T}"     , 500, 0.0  , 500.0);
  auto  histMuEta      = new TH1F ("histMuEta"    , "Distribution of reco-muon #eta"      , 500, -2.50, 2.5);
  auto  histMuGenPt    = new TH1F ("histMuGenPt"  , "Distribution of gen-muon p_{T}"      , 500, 0.0  , 500.0);
  auto  histMuGenEta   = new TH1F ("histMuGenEta" , "Distribution of gen-muon #eta"       , 500, -2.5 , 2.5);

  auto histEleGenPtEta = new TH2F{"histEleGenPtEta",
      "Distribution of reco-electron p_{T} against #eta", 500, 0, 300, 500, -3, 3};
  auto histMuGenPtEta = new TH2F{"histMuGenPtEta",
      "Distribution of reco-muon p_{T} against #eta", 500, 0, 300, 500, -3, 3};

  auto  lTimer = new TMVA::Timer ( inputTrees.size(), "Running over trees", false );

  lTimer->DrawProgressBar(0, "");

  Int_t lCounter (1);

  // pT thresholds
  const unsigned ele1pTThresh{17};
  const unsigned ele1pTThreshProposed{23};
  const unsigned ele2pTThresh{12};
  // Event counters
  unsigned totalEvents{0};
  unsigned pTCut{0};
  unsigned eleNumCut{0};
  unsigned newlyCut{0};

  for ( std::vector<TTree*>::const_iterator lIt = inputTrees.begin(); lIt != inputTrees.end(); ++lIt ){
    
    AnalysisEvent* lEvent = new AnalysisEvent(true, "null", *lIt);

    Int_t lNumEvents = (*lIt)->GetEntries();
    totalEvents += lNumEvents;
    
    for ( Int_t j = 0; j < lNumEvents; j++ ){
      (*lIt)->GetEvent(j);
      
      for ( Int_t k = 0; k < lEvent->numLooseElePF2PAT; k++){
	histElePt->Fill(lEvent->elePF2PATlooseElectronSortedPt[k]);
	histEleEta->Fill(lEvent->elePF2PATlooseElectronSortedEta[k]);
	histEleGenPt->Fill(lEvent->genLooseElePF2PATPT[k]);
	histEleGenEta->Fill(lEvent->genLooseElePF2PATEta[k]);

        histEleGenPtEta->Fill(lEvent->elePF2PATlooseElectronSortedPt[k], 
            lEvent->elePF2PATlooseElectronSortedEta[k]);

      }
      for ( Int_t k = 0; k < lEvent->numMuonPF2PAT; k++){
	histMuPt->Fill(lEvent->muonPF2PATlooseMuonSortedPt[k]);
	histMuEta->Fill(lEvent->muonPF2PATlooseMuonSortedEta[k]);
  	histMuGenPt->Fill(lEvent->genLooseMuonPF2PATPT[k]);
    	histMuGenEta->Fill(lEvent->genLooseMuonPF2PATEta[k]);

        histMuGenPtEta->Fill(lEvent->muonPF2PATlooseMuonSortedPt[k], 
            lEvent->muonPF2PATlooseMuonSortedEta[k]);
      }

      // Count the number of events which will be cut
      if (lEvent->numLooseElePF2PAT >= 2){  // electron no. cut
        std::vector<Float_t> elePts;

        for (Int_t k = 0; k < lEvent->numLooseElePF2PAT; k++){
          elePts.emplace_back(lEvent->elePF2PATlooseElectronSortedPt[k]);
        }
        std::nth_element(elePts.begin(), elePts.begin()+1, elePts.end(),
            std::greater<Float_t>());
        if (elePts.at(1) < ele2pTThresh || elePts.at(0) < ele1pTThresh){
          pTCut++;
        }
        else if (elePts.at(0) < ele1pTThreshProposed){
          pTCut++;
          newlyCut++;
        }
      }
      else {
        eleNumCut++;
      }
    }
    lTimer->DrawProgressBar(lCounter++, "");
  }

  std::cout << std::endl;
  std::cout << "Total no. of evets:\t\t" <<  totalEvents << std::endl;
  std::cout << "Cut due to pT requirements:\t" << pTCut << std::endl;
  std::cout << "Cut due to too few electrons:\t" << eleNumCut << std::endl;
  std::cout << "Total no. cut events:\t\t" << pTCut + eleNumCut << std::endl;
  std::cout << "Cut due to new proposals:\t" << newlyCut << std::endl;

  auto outFile = new TFile ( outFileString.c_str(), "RECREATE" );
  
  histElePt->Write();
  histEleEta->Write();
  histEleGenPt->Write();
  histEleGenEta->Write();

  histMuPt->Write();
  histMuEta->Write();
  histMuGenPt->Write();
  histMuGenEta->Write();

  histEleGenPtEta->Write();
  histMuGenPtEta->Write();

  outFile->Close();
  std::cout << "\n Finished." << std::endl;
}

