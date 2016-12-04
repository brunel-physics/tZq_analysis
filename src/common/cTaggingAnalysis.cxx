#include "AnalysisEvent.hpp"
#include <libconfig.h++>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMVA/Timer.h"

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {

  std::string inputDir;
  std::string outFileString{"plots/distributions/output.root"};

  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
      ("help,h", "Print this message.")
      ("indir,i", po::value<std::string>(&inputDir)->required(), "Input folder for nTuples.")
      ("outfile,o", po::value<std::string>(&outFileString)->default_value(outFileString),
       "Output file for plots.");
  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
      std::cout << desc;
      return 0;
    }

    po::notify(vm);
  }
  catch (const po::error& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }

  std::vector<TTree*> inputTrees;

  if (boost::filesystem::is_directory(inputDir))
  {
      const long int count{std::distance
          (boost::filesystem::directory_iterator{inputDir},
           boost::filesystem::directory_iterator())};

      TMVA::Timer lTimer{count, "Attaching files to TTree", false};
      Int_t lCounter{1};

      lTimer.DrawProgressBar(0, "");

      for(const auto& file: boost::make_iterator_range(
                  boost::filesystem::directory_iterator{inputDir}, {}))
      {
          TFile *inputFile{new TFile {file.path().string().c_str()}};
          TTree *lTempTree{dynamic_cast<TTree*>(inputFile->Get("tree"))};
          inputTrees.emplace_back(lTempTree);

          lTimer.DrawProgressBar(lCounter++, "");
      }
  }
  else
  {
      std::cout << "ERROR: " << inputDir << "is not a valid directory" <<
          std::endl;
      return 1;
  }

  std::cout << std::endl;

  constexpr int discriminatorIncrement{100};

  // Initialise Histograms.
  TH2F* h_CvsBjetPassEfficiency{new TH2F {"h_CvsBjetPassEfficiency", "c/b jet (vs b) Pass Effifiency", discriminatorIncrement, -1.0, 1.0, 100, 0.0, 1.0}};
  TH2F* h_CvsLjetPassEfficiency{new TH2F {"h_CvsLjetPassEfficiency", "c/l jet (vs l) Pass Effifiency", discriminatorIncrement, -1.0, 1.0, 100, 0.0, 1.0}};

  TH2F* h_CvsB_cJetPassEfficiency{new TH2F {"h_CvsBcJetPassEfficiency", "c-jet (vs b) Pass Effifiency", discriminatorIncrement, -1.0, 1.0, 100, 0.0, 1.0}};
  TH2F* h_CvsL_cJetPassEfficiency{new TH2F {"h_CvsLcJetPassEfficiency", "c-jet (vs l) Pass Effifiency", discriminatorIncrement, -1.0, 1.0, 100, 0.0, 1.0}};

  // Setup counters.
  //  int lTotalNumJets (0);
  int lTotalNumLjets{0};
  int lTotalNumBjets{0};
  int lTotalNumCjets{0};
  int lTotalNumLjetsCvsL[discriminatorIncrement] {0};
  int lTotalNumBjetsCvsB[discriminatorIncrement] {0};
  int lTotalNumCjetsCvsL[discriminatorIncrement] {0};
  int lTotalNumCjetsCvsB[discriminatorIncrement] {0};

  double lPassedCvsBjets[discriminatorIncrement] {0};
  double lPassedCvsLjets[discriminatorIncrement] {0};
  double lPassedCvsBjetFraction[discriminatorIncrement] {0};
  double lPassedCvsLjetFraction[discriminatorIncrement] {0};
  double lPassedCvsB_cJetFraction[discriminatorIncrement] {0};
  double lPassedCvsL_cJetFraction[discriminatorIncrement] {0};

  // Initial progress bar.
  TMVA::Timer *lTimer{new TMVA::Timer {boost::numeric_cast<int>(inputTrees.size()), "Running over trees", false}};
  TMVA::Timer *lTimer2{new TMVA::Timer {discriminatorIncrement, "Filling plots", false}};
  lTimer->DrawProgressBar(0, "");

  Int_t lCounter{1};

  // Start looping over jets.

  for ( std::vector<TTree*>::const_iterator lIt{inputTrees.begin()}; lIt != inputTrees.end(); ++lIt ){

    AnalysisEvent* lEvent{new AnalysisEvent{true, "null", *lIt}};

    Long64_t lNumEvents{(*lIt)->GetEntries()};

    for ( Int_t j{0}; j < lNumEvents; j++ ){
      (*lIt)->GetEvent(j);

      for (Int_t k{0}; k < lEvent->numJetPF2PAT; k++){

        Double_t lPt{lEvent->jetPF2PATPt[k]};
        Double_t lEta{lEvent->jetPF2PATEta[k]};
        Int_t lFlavour{lEvent->genJetPF2PATPID[k]};
        Double_t lCvsBdisc{lEvent->jetPF2PATCvsBDiscriminator[k]};
        Double_t lCvsLdisc{(lEvent->jetPF2PATCvsLDiscriminator[k])};

	if ( lPt < 30.0 && std::abs( lEta > 5.0 )) continue; // Check it isn't out of eta range or is too soft.
	bool jetID{false};

	if (
	    ( std::abs(lEta<=3.0) && lEvent->jetPF2PATNeutralHadronEnergyFraction[k] < 0.99 && lEvent->jetPF2PATNeutralEmEnergyFraction[k] < 0.99 &&  (lEvent->jetPF2PATNeutralMultiplicity[k] && lEvent->jetPF2PATChargedMultiplicity[k]) > 1.00 &&
	      ( ( std::abs(lEta<=2.40) && (lEvent->jetPF2PATChargedHadronEnergyFraction[k] > 0 && lEvent->jetPF2PATChargedMultiplicity[k] > 0 && lEvent->jetPF2PATChargedEmEnergyFraction[k] < 0.99)) && std::abs(lEta > 2.40) ) )
	    || (std::abs(lEta)>3.0 && lEvent->jetPF2PATNeutralEmEnergyFraction[k] < 0.90 && lEvent->jetPF2PATNeutralMultiplicity[k] > 10)
	    )
	  jetID = true;
	
	if (!jetID) continue;

	if ( std::abs(lFlavour) <= 3 && std::abs(lFlavour) != 0 ) ++lTotalNumLjets;
	if ( std::abs(lFlavour) == 4 ) ++lTotalNumCjets;
	if ( std::abs(lFlavour) == 5 ) ++lTotalNumBjets;


	int lArrayIt{0};
	for (float lDiscIt{-1}; lDiscIt <= 1.0; lDiscIt += (discriminatorIncrement/2.0), ++lArrayIt){

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

  double lTempDiscr{-1};
  std::cout << "\n" << std::endl;
  lTimer2->DrawProgressBar(0, "");
  Int_t lCounter2{1};

  for (int i{0}; i != discriminatorIncrement; i++, lTempDiscr += 2/discriminatorIncrement){

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

  TFile *outFile{new TFile{outFileString.c_str(), "RECREATE"}};

  h_CvsBjetPassEfficiency->Write();
  h_CvsLjetPassEfficiency->Write();

  h_CvsB_cJetPassEfficiency->Write();
  h_CvsL_cJetPassEfficiency->Write();

  outFile->Close();
  std::cout << "\n Finished." << std::endl;
}
