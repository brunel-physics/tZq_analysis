#include "AnalysisEvent.hpp"
#include <libconfig.h++>
#include <dirent.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMVA/Timer.h"

#include <boost/numeric/conversion/cast.hpp>
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

  std::string inputDir{};
  std::string outFileString{"plots/distributions/output.root"};
  bool inputFlag{false};

  for (int i{1}; i < argc; ++i){
    std::string arg{argv[i]};
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
    std::cout << inputDir << " is not a valid directory" << std::endl;
    return 0;
  }

  std::vector<TTree*> inputTrees;

  std::cout << "Attaching files to TTree" << std::flush;

  while ( (dirp = readdir(dp)) != nullptr) {
    std::string line{dirp->d_name};
    if ( line == "." || line == "..")
    continue;
    TFile *inputFile{new TFile {(inputDir+line).c_str()}};
    TTree *lTempTree{dynamic_cast<TTree*>(inputFile->Get("tree"))};
    inputTrees.emplace_back(lTempTree);

    std::cout << '.' << std::flush;
  }
  std::cout << std::endl;

  std::cout << "Attached all files to TTree!" << std::endl;

  TH1F *histElePt    {new TH1F{"histElePt"    , "Distribution of reco-electron p_{T}", 500, 0.0  , 500.0}};
  TH1F *histEleEta   {new TH1F{"histEleEta"   , "Distribution of reco-electron #eta" , 500, -2.50, 2.5}};
  TH1F *histEleGenPt {new TH1F{"histEleGenPt" , "Distribution of gen-electron p_{T}" , 500, 0.0  , 500.0}};
  TH1F *histEleGenEta{new TH1F{"histEleGenEta", "Distribution of gen-electron #eta"  , 500, -2.5 , 2.5}};

  TH1F *histMuPt     {new TH1F{"histMuPt"     , "Distribution of reco-muon p_{T}"    , 500, 0.0  , 500.0}};
  TH1F *histMuEta    {new TH1F{"histMuEta"    , "Distribution of reco-muon #eta"     , 500, -2.50, 2.5}};
  TH1F *histMuGenPt  {new TH1F{"histMuGenPt"  , "Distribution of gen-muon p_{T}"     , 500, 0.0  , 500.0}};
  TH1F *histMuGenEta {new TH1F{"histMuGenEta" , "Distribution of gen-muon #eta"      , 500, -2.5 , 2.5}};

  TH2F *histEleGenPtEta{new TH2F{"histEleGenPtEta",
        "Distribution of reco-electron p_{T} against #eta", 500, 0, 300, 500, -3, 3}};
  TH2F *histMuGenPtEta{new TH2F{"histMuGenPtEta",
        "Distribution of reco-muon p_{T} against #eta", 500, 0, 300, 500, -3, 3}};

  TMVA::Timer *lTimer{new TMVA::Timer{boost::numeric_cast<int>(inputTrees.size()), "Running over trees", false}};

  lTimer->DrawProgressBar(0, "");

  Int_t lCounter{1};

  // pT thresholds
  const Float_t ele1PtThresh{17};
  const Float_t ele2PtThresh{12};
  const Float_t mu1PtThresh{17};
  const Float_t mu2PtThresh{8};
  const Float_t ele1PtThreshProposed{23};
  // eta thresholds
  const Float_t eleEtaThresh{2.5};
  const Float_t muEtaThresh{2.4};
  // Event counters
  int totalEvents{0};
  int passEleNum{0};
  int passElePt{0};
  int passElePtDiff{0};
  int passEleEta{0};
  int passMuNum{0};
  int passMuPt{0};
  int passMuEta{0};

  int newlyCut{0};
  int totalCut{0};

  for ( std::vector<TTree*>::const_iterator lIt = inputTrees.begin(); lIt != inputTrees.end(); ++lIt ){
    
    AnalysisEvent* lEvent{new AnalysisEvent{true, "null", *lIt}};

    Long64_t lNumEvents{(*lIt)->GetEntries()};
    totalEvents += lNumEvents;
    
    for ( Int_t j{0}; j < lNumEvents; j++ ){
      (*lIt)->GetEvent(j);
      
      for ( Int_t k{0}; k < lEvent->numElePF2PAT; k++){
	histElePt->Fill(lEvent->elePF2PATPT[k]);
	histEleEta->Fill(lEvent->elePF2PATEta[k]);
	histEleGenPt->Fill(lEvent->genElePF2PATPT[k]);
	histEleGenEta->Fill(lEvent->genElePF2PATEta[k]);

        histEleGenPtEta->Fill(lEvent->elePF2PATPT[k], 
            lEvent->elePF2PATEta[k]);

      }
      for ( Int_t k{0}; k < lEvent->numMuonPF2PAT; k++){
	histMuPt->Fill(lEvent->muonPF2PATPt[k]);
	histMuEta->Fill(lEvent->muonPF2PATEta[k]);
  	histMuGenPt->Fill(lEvent->genMuonPF2PATPT[k]);
    	histMuGenEta->Fill(lEvent->genMuonPF2PATEta[k]);

        histMuGenPtEta->Fill(lEvent->muonPF2PATPt[k], 
            lEvent->muonPF2PATEta[k]);
      }

      // Count the number of events which will be cut
      struct particle
      {
        particle(const Float_t a, const Float_t b, const Int_t c, const Int_t d) :
          eta(a),
          pt(b),
          Id(c),
          motherId(d)
        {}

        Float_t eta;
        Float_t pt;
        Int_t Id;
        Int_t motherId;
      };

      bool cut{true};
      bool cutDiff{false};

      std::vector<particle> eles;
      std::vector<particle> mus;
      for(int k = 0; k < lEvent->nGenPar; k++)
      {
        switch (std::abs(lEvent->genParId[k]))
        {
          case 11:  // electron
            eles.emplace_back(lEvent->genParEta[k],
                lEvent->genParPt[k],
                lEvent->genParId[k],
                lEvent->genParMotherId[k]);
            break;
          case 13:  // muon
            mus.emplace_back(lEvent->genParEta[k],
                lEvent->genParPt[k],
                lEvent->genParId[k],
                lEvent->genParMotherId[k]);
            break;
          default:
            break;
        }
      }

      if (eles.size() >= 2)  // electron no. cut
      {
        passEleNum++;

        eles.erase(std::remove_if(eles.begin(), eles.end(),
            [&](const particle &p) -> bool {return std::abs(p.eta) > eleEtaThresh;}),
          eles.end());

        if (eles.size() >= 2)  // electron eta cut
        {
          passEleEta++;

          std::nth_element(eles.begin(), eles.begin()+1, eles.end(),
              [](const particle &a, const particle &b) -> bool
              {return a.pt > b.pt;});

          if (eles.at(0).pt > ele1PtThresh && eles.at(1).pt > ele2PtThresh)  // ele pT cut
          {
            if (eles.at(0).pt > ele1PtThreshProposed)
            {
              passElePt++;
              cut = false;
            }
            else
            {  // fail proposed pT cut
              passElePtDiff++;
              cutDiff = true;
            }
          }
        }
      }

      if (mus.size() >= 2)  // muon no. cut
      {
        passMuNum++;

        mus.erase(std::remove_if(mus.begin(), mus.end(),
            [&](const particle &p) -> bool {return std::abs(p.eta) > muEtaThresh;}),
          mus.end());

        if (mus.size() >= 2)  // muon eta cut
        {
          passMuEta++;

          std::nth_element(mus.begin(), mus.begin()+1, mus.end(),
              [](const particle &a, const particle &b) -> bool
              {return a.pt > b.pt;});

          if (mus.at(0).pt > mu1PtThresh && mus.at(1).pt > mu2PtThresh)  // muon pT cut
          {
            passMuPt++;
            cut = false;
            cutDiff = false;
          }
        }
      }

      if (cut)
      {
        totalCut++;
      }
      if (cutDiff)
      {
        newlyCut++;
      }
    }
    lTimer->DrawProgressBar(lCounter++, "");
  }

  std::cout << std::endl << std::endl;
  std::cout << "Total no. of events:\t\t\t" << totalEvents << std::endl;
  std::cout << std::endl;
  std::cout << "Containing at least two electrons:\t" << passEleNum << std::endl;
  std::cout << "...of which pass eta requirements:\t" << passEleEta << std::endl;
  std::cout << "...of which pass pT requirements:\t" << passElePt << std::endl;
  std::cout << "Change due to new proposals:\t\t" << passElePtDiff << std::endl;
  std::cout << std::endl;
  std::cout << "Containing at least two muons:\t\t" << passMuNum << std::endl;
  std::cout << "...of which pass eta requirements:\t" << passMuEta << std::endl;
  std::cout << "...of which pass pT requirements:\t" << passMuPt << std::endl;
  std::cout << std::endl;
  std::cout << "Total no. of cut events\t\t\t" << totalCut << std::endl;
  std::cout << "Change due to new proposals:\t\t" << newlyCut << std::endl;

  TFile *outFile{new TFile{outFileString.c_str(), "RECREATE"}};

  histElePt->Write();
  histEleEta->Write();
  histEleGenPt->Write();
  histEleGenEta->Write();

  histMuPt->Write();
  histMuEta->Write();
  histMuGenPt->Write();
  histMuGenEta->Write();

  histEleGenPtEta->SetStats(kFALSE);
  histEleGenPtEta->SetTitle("");
  histEleGenPtEta->GetXaxis()->SetTitle("p_{T}");
  histEleGenPtEta->GetYaxis()->SetTitle("#eta");
  histEleGenPtEta->SetTitleOffset(0.5f, "Y");
  histEleGenPtEta->Write();

  histMuGenPtEta->SetStats(kFALSE);
  histMuGenPtEta->SetTitle("");
  histMuGenPtEta->GetXaxis()->SetTitle("p_{T}");
  histMuGenPtEta->GetYaxis()->SetTitle("#eta");
  histMuGenPtEta->SetTitleOffset(0.5f, "Y");
  histMuGenPtEta->Write();

  outFile->Close();
  std::cout << "\n Finished." << std::endl;
}
