#include "AnalysisEvent.hpp"
#include <libconfig.h++>
#include <dirent.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMVA/Timer.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

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

  // Read in root files from directory
  DIR *dp;
  struct dirent *dirp;

  if((dp  = opendir( inputDir.c_str() )) == nullptr) {
    std::cout << "Error opening Directory" << std::endl;
    std::cout << inputDir << " is not a valid directory" << std::endl;
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

  const int discriminatorIncrement (100);

  // Initialise Histograms.
  auto  h_CvsBjetPassEfficiency = new TH2F ("h_CvsBjetPassEfficiency", "c/b jet (vs b) Pass Effifiency", discriminatorIncrement, -1.0, 1.0, 100, 0.0, 1.0);
  auto  h_CvsLjetPassEfficiency = new TH2F ("h_CvsLjetPassEfficiency", "c/l jet (vs l) Pass Effifiency", discriminatorIncrement, -1.0, 1.0, 100, 0.0, 1.0);

  auto  h_CvsB_cJetPassEfficiency = new TH2F ("h_CvsBcJetPassEfficiency", "c-jet (vs b) Pass Effifiency", discriminatorIncrement, -1.0, 1.0, 100, 0.0, 1.0);
  auto  h_CvsL_cJetPassEfficiency = new TH2F ("h_CvsLcJetPassEfficiency", "c-jet (vs l) Pass Effifiency", discriminatorIncrement, -1.0, 1.0, 100, 0.0, 1.0);

  // Setup counters.
  //  int lTotalNumJets (0);
  int lTotalNumLjets (0);
  int lTotalNumBjets (0);
  int lTotalNumCjets (0);
  int lTotalNumLjetsCvsL [discriminatorIncrement] = {0};
  int lTotalNumBjetsCvsB [discriminatorIncrement] = {0};
  int lTotalNumCjetsCvsL [discriminatorIncrement] = {0};
  int lTotalNumCjetsCvsB [discriminatorIncrement] = {0};

  double lPassedCvsBjets [discriminatorIncrement] = {0.0};
  double lPassedCvsLjets [discriminatorIncrement] = {0.0};
  double lPassedCvsBjetFraction [discriminatorIncrement] = {0.0};
  double lPassedCvsLjetFraction [discriminatorIncrement] = {0.0};
  double lPassedCvsB_cJetFraction [discriminatorIncrement] = {0.0};
  double lPassedCvsL_cJetFraction [discriminatorIncrement] = {0.0};

  // Initial progress bar.
  auto  lTimer = new TMVA::Timer ( inputTrees.size(), "Running over trees", true );
  auto  lTimer2 = new TMVA::Timer ( discriminatorIncrement, "Filling plots", true );
  lTimer->DrawProgressBar(0, "");

  Int_t lCounter (1);

  // Start looping over jets.

  for ( std::vector<TTree*>::const_iterator lIt = inputTrees.begin(); lIt != inputTrees.end(); ++lIt ){
    
    AnalysisEvent* lEvent = new AnalysisEvent(true, "null", *lIt);

    Int_t lNumEvents = (*lIt)->GetEntries();
    
    for ( Int_t j = 0; j < lNumEvents; j++ ){
      (*lIt)->GetEvent(j);
      
      for (Int_t k = 0; k < lEvent->numJetPF2PAT; k++){

	Double_t lPt = (lEvent->jetPF2PATPt[k]);
	Double_t lEta = (lEvent->jetPF2PATEta[k]);
	Int_t lFlavour = (lEvent->genJetPF2PATPID[k]);
	Double_t lCvsBdisc = (lEvent->jetPF2PATCvsBDiscriminator[k]);
	Double_t lCvsLdisc = (lEvent->jetPF2PATCvsLDiscriminator[k]);

	if ( lPt < 30.0 && std::abs( lEta > 5.0 )) continue; // Check it isn't out of eta range or is too soft.
	bool jetID (false);

	if (
	    ( std::abs(lEta<=3.0) && lEvent->jetPF2PATNeutralHadronEnergyFraction[k] < 0.99 && lEvent->jetPF2PATNeutralEmEnergyFraction[k] < 0.99 &&  (lEvent->jetPF2PATNeutralMultiplicity[k] && lEvent->jetPF2PATChargedMultiplicity[k]) > 1 && 
	      ( ( std::abs(lEta<=2.40) && (lEvent->jetPF2PATChargedHadronEnergyFraction[k] > 0 && lEvent->jetPF2PATChargedMultiplicity[k] > 0 && lEvent->jetPF2PATChargedEmEnergyFraction[k] < 0.99)) && std::abs(lEta > 2.40) ) ) 
	    || (std::abs(lEta)>3.0 && lEvent->jetPF2PATNeutralEmEnergyFraction[k] < 0.90 && lEvent->jetPF2PATNeutralMultiplicity[k] > 10) 
	    )
	  jetID = true;
	
	if (!jetID) continue;

	if ( std::abs(lFlavour) <= 3 && std::abs(lFlavour) != 0 ) ++lTotalNumLjets;
	if ( std::abs(lFlavour) == 4 ) ++lTotalNumCjets;
	if ( std::abs(lFlavour) == 5 ) ++lTotalNumBjets;


	int lArrayIt(0);
	for (float lDiscIt = -1.0; lDiscIt <= 1.0; lDiscIt += (discriminatorIncrement/2.0), ++lArrayIt){

	  if (lCvsBdisc >= lArrayIt) ++lPassedCvsBjets[lArrayIt];
	  if (lCvsLdisc >= lArrayIt) ++lPassedCvsLjets[lArrayIt];

	  if (lCvsBdisc >= lArrayIt && std::abs(lFlavour)<= 3 && lFlavour != 0) ++lTotalNumLjetsCvsL[lArrayIt];
	  if (lCvsBdisc >= lArrayIt && std::abs(lFlavour) == 5) ++lTotalNumBjetsCvsB[lArrayIt];

	  if (lCvsBdisc >= lArrayIt && std::abs(lFlavour) == 4) ++lTotalNumCjetsCvsB[lArrayIt];
	  if (lCvsLdisc >= lArrayIt && std::abs(lFlavour) == 4) ++lTotalNumCjetsCvsL[lArrayIt];
	} 

      }
    }
  lTimer->DrawProgressBar(lCounter++, "");
  }

  double lTempDiscr (-1.0);
  std::cout << "\n" << std::endl;
  lTimer2->DrawProgressBar(0, "");
  Int_t lCounter2 (1);

  for (int i = 0; i != discriminatorIncrement; i++, lTempDiscr += 2/discriminatorIncrement){

    lPassedCvsBjetFraction[i] = lPassedCvsBjets[i]/(lTotalNumBjets+lTotalNumCjets);
    lPassedCvsLjetFraction[i] = lPassedCvsLjets[i]/(lTotalNumLjets+lTotalNumCjets);

    h_CvsBjetPassEfficiency->Fill(lTempDiscr, lPassedCvsBjetFraction[i]);
    h_CvsLjetPassEfficiency->Fill(lTempDiscr, lPassedCvsLjetFraction[i]);

    lPassedCvsB_cJetFraction [i] = lTotalNumCjetsCvsB[i]/lPassedCvsBjets[i];
    lPassedCvsL_cJetFraction [i] = lTotalNumCjetsCvsL[i]/lPassedCvsLjets[i];

    h_CvsB_cJetPassEfficiency->Fill( lTempDiscr, lPassedCvsB_cJetFraction[i]);
    h_CvsL_cJetPassEfficiency->Fill( lTempDiscr, lPassedCvsL_cJetFraction[i]);

    lTimer2->DrawProgressBar(lCounter2++, "");
  }

  auto outFile = new TFile ( outFileString.c_str(), "RECREATE" );
   
  h_CvsBjetPassEfficiency->Write();
  h_CvsLjetPassEfficiency->Write();

  h_CvsB_cJetPassEfficiency->Write();
  h_CvsL_cJetPassEfficiency->Write();

  outFile->Close();
  std::cout << "\n Finished." << std::endl;
}
