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
  const unsigned ele1PtThresh{17};
  const unsigned ele2PtThresh{12};
  const unsigned mu1PtThresh{17};
  const unsigned mu2PtThresh{8};
  const unsigned ele1PtThreshProposed{23};
  // Event counters
  int totalEvents{0};
  int elePtCut{0};
  int eleNumCut{0};
  int muPtCut{0};
  int muNumCut{0};
  int newlyCut{0};
  int totalCut{0};

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
      bool cut{true};
      if (lEvent->numLooseElePF2PAT >= 2)  // electron no. cut
      {
        std::vector<Float_t> elePTs{lEvent->genLooseElePF2PATPT,
          lEvent->genLooseElePF2PATPT + lEvent->numLooseElePF2PAT};

        std::nth_element(elePTs.begin(), elePTs.begin()+1, elePTs.end(),
            std::greater<Float_t>());

        if (elePTs.at(1) < ele2PtThresh || elePTs.at(0) < ele1PtThresh)  // ele pT cut
        {
          elePtCut++;
        }
        else if (elePTs.at(0) < ele1PtThreshProposed)
        {
          elePtCut++;
          newlyCut++;
        }
        else
        {
          cut = false;
        }
      }
      else
      {
        eleNumCut++;
      }

      if (lEvent->numMuonPF2PAT >=2)  // muon no. cut
      {
        std::vector<Float_t> muPTs{lEvent->genLooseMuonPF2PATPT,
          lEvent->genLooseMuonPF2PATPT + lEvent->numMuonPF2PAT};

        std::nth_element(muPTs.begin(), muPTs.begin()+1, muPTs.end(),
            std::greater<Float_t>());

        if (muPTs.at(1) < mu2PtThresh || muPTs.at(0) < mu1PtThresh)  // mu pT cut
        {
          muPtCut++;
        }
        else
        {
          cut = false;
        }
      }
      else
      {
        muNumCut++;
      }

      if (cut)
      {
        totalCut++;
      }
    }
    lTimer->DrawProgressBar(lCounter++, "");
  }

  std::cout << std::endl << std::endl;
  std::cout << "Total no. of events:\t\t\t" << totalEvents << std::endl;
  std::cout << std::endl;
  std::cout << "Fail electron pT requirements:\t\t" << elePtCut << std::endl;
  std::cout << "Fail due to too few electrons:\t\t" << eleNumCut << std::endl;
  std::cout << "Total failing electron requirements:\t" << elePtCut + eleNumCut << std::endl;
  std::cout << "Change due to new proposals:\t\t" << newlyCut << std::endl;
  std::cout << std::endl;
  std::cout << "Fail muon pT requirements:\t\t" << muPtCut << std::endl;
  std::cout << "Fail due to too few muons:\t\t" << muNumCut << std::endl;
  std::cout << "Total failing muon requirements:\t" << muPtCut + muNumCut << std::endl;
  std::cout << std::endl;
  std::cout << "Total no. of cut events\t\t\t" << totalCut << std::endl;

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
