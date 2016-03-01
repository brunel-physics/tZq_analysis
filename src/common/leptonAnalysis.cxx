#include "AnalysisEvent.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include <stdlib.h> 
#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char* argv[]) {

  Int_t lNumberOfFiles = argc;
  std::cout << "lNumberOfFiles: " << lNumberOfFiles-1 << std::endl;

  TH1F* histElePt      = new TH1F ("histElePt"    , "Distribution of reco-electron p_{T}" , 500, 0.0  , 500.0);
  TH1F* histEleEta     = new TH1F ("histEleEta"   , "Distribution of reco-electron #eta"  , 500, -2.50, 2.5);
  TH1F* histEleGenPt   = new TH1F ("histEleGenPt" , "Distribution of gen-electron p_{T}"  , 500, 0.0  , 500.0);
  TH1F* histEleGenEta  = new TH1F ("histEleGenEta", "Distribution of gen-electron #eta"   , 500, -2.5 , 2.5);

  TH1F* histMuPt       = new TH1F ("histMuPt"     , "Distribution of reco-muon p_{T}"     , 500, 0.0  , 500.0);
  TH1F* histMuEta      = new TH1F ("histMuEta"    , "Distribution of reco-muon #eta"      , 500, -2.50, 2.5);
  TH1F* histMuGenPt    = new TH1F ("histMuGenPt"  , "Distribution of gen-muon p_{T}"      , 500, 0.0  , 500.0);
  TH1F* histMuGenEta   = new TH1F ("histMuGenEta" , "Distribution of gen-muon #eta"       , 500, -2.5 , 2.5);


  for ( Int_t i = 1; i < lNumberOfFiles; i++ ){
    std::stringstream lSStr;
    lSStr << argv[i];
    std::cout << "Loading file: " << lSStr.str().c_str() << std::endl;

    TFile *inputFile = new TFile ( (lSStr.str()).c_str() );
    TTree *lTree = (TTree*)inputFile->Get("tree");       
    AnalysisEvent* lEvent = new AnalysisEvent(true, "null", lTree);

    Int_t lNumEvents = lTree->GetEntries();
    std::cout << "lNumEvents: " << lNumEvents << std::endl;

    for ( Int_t j = 0; j < lNumEvents; j++ ){
      lTree->GetEvent(j);

      for ( Int_t k = 0; k < lEvent->numElePF2PAT; k++){
	histElePt->Fill(lEvent->elePF2PATPT[k]);
	histEleEta->Fill(lEvent->elePF2PATEta[k]);
	//	  histEleGenPt->Fill(lEvent->genElePATPT[k]);
	histEleGenEta->Fill(lEvent->genElePF2PATEta[k]);
      }
      for ( Int_t l = 0; l < lEvent->numMuonPF2PAT; l++){
	histMuPt->Fill(lEvent->muonPF2PATPt[l]);
	histMuEta->Fill(lEvent->muonPF2PATEta[l]);
	//	  histMuGenPt->Fill(lEvent->genMuonPATPT[l]);
	histMuGenEta->Fill(lEvent->genMuonPF2PATEta[l]);
      }
    }
  }

  TFile *outFile = new TFile ( "pT_eta_distributions.root", "RECREATE" );

  histElePt->Write();
  histEleEta->Write();
  //    histEleGenPt->Write();
  histEleGenEta->Write();

  histMuPt->Write();
  histMuEta->Write();
  //    histMuGenPt->Write();
  histMuGenEta->Write();

  outFile->Close();
}

