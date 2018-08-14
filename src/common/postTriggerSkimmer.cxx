#include <array>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/range/iterator_range.hpp>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include "AnalysisEvent.hpp"


namespace fs = boost::filesystem;

int main (int argc, char* argv[])
{

  std::vector<std::string> dileptonDirs;
  std::vector<std::string> singleLeptonDirs;
  std::string datasetName;
  std::string channel;

  int singleElectron {0}, dupElectron {0};
  int singleMuon {0}, dupMuon {0};

  const std::string postTriggerSkimDir{std::string{"/data0/data/TopPhysics/postTriggerSkims2016/"}};

  // Define command-line flags
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Print this message.")
    ("channel,c", po::value<std::string>(&channel)->required(),
     "Channel to operate over. Either ee, emu or mumu.")
    ("dileptonDirs,d", po::value<std::vector<std::string>>
     (&dileptonDirs)->multitoken()->required(),
     "Directories in which to look for double lepton datasets.")
    ("singleLeptonDirs,s", po::value<std::vector<std::string>>
     (&singleLeptonDirs)->multitoken()->required(),
     "Directories in which to look for single lepton datasets.")
    ("datasetName,o", po::value<std::string>(&datasetName)->required(),
     "Output dataset name.");
  po::variables_map vm;

  // Parse arguments
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

  
  const std::regex mask{".*\\.root"};
  int fileNum{0};
  
  std::set< std::pair < Int_t, Int_t > > triggerDoubleCountCheck;
 
  for (const auto& dileptonDir: dileptonDirs) {  // for each dilepton input directory
    for (const auto& file : boost::make_iterator_range(fs::directory_iterator{dileptonDir}, {}))  {  // for each file in directory
      const std::string path {file.path().string()};
      std::cout << "path: " << path << std::endl;

      if (!fs::is_regular_file(file.status()) || !std::regex_match(path, mask)) {
	continue;  // skip if not a root file
      }

      const std::string numName{std::to_string(fileNum)};
      const std::string numNamePlus{std::to_string(fileNum + 2)};
      
      if (fs::is_regular_file( postTriggerSkimDir + datasetName + "/triggerSkim" + numNamePlus + ".root" )) {
	// don't overwrite existing skim files, except for the last two
	fileNum++;
	continue;
      }
      
      TChain datasetChain{"tree"};
      datasetChain.Add(path.c_str());

      TTree * const outTree = datasetChain.CloneTree(0);
      
      std::string outFilePath{ postTriggerSkimDir + datasetName + "/triggerSkim" + numName + ".root" };
      TFile outFile{outFilePath.c_str(), "RECREATE"};
	  
      const long long int numberOfEvents{datasetChain.GetEntries()};
      
      boost::progress_display progress{numberOfEvents, std::cout, outFilePath + "\n"};
      AnalysisEvent event{false, "", &datasetChain, true};

      for (long long int i{0}; i < numberOfEvents; i++) {
	++progress;  // update progress bar (++ must be prefix)
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

          // non-DZ legs are prescaled for Run2016H
          if ( event.eventRun < 280919 ) {
            if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2 > 0 ) mumuTrig = true;
            if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3 > 0 ) mumuTrig = true;
            if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4 > 0 ) mumuTrig = true;
            if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6 > 0 ) mumuTrig = true;

            if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2  > 0 ) mumuTrig = true;
            if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3  > 0 ) mumuTrig = true;
            if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5  > 0 ) mumuTrig = true;
          }

          // DZ legs avaliable all the time but inefficient in data for Runs B-F -> hence uses of non-DZ legs

	  if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
	  if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 > 0 ) mumuTrig = true;
	  if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 ) mumuTrig = true;
	  if ( event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 > 0 ) mumuTrig = true;
	  if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
	  if ( event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 ) mumuTrig = true;
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
      std::cout << std::endl;

     delete outTree;
    }
  }
  
  for (const auto& singleLeptonDir: singleLeptonDirs) { // for each single lepton input directory
    for (const auto& file : boost::make_iterator_range(fs::directory_iterator{singleLeptonDir}, {}))  {  // for each file in directory
      const std::string path {file.path().string()};
      std::cout << "path: " << path << std::endl;

      if (!fs::is_regular_file(file.status()) || !std::regex_match(path, mask)) {	
	continue;  // skip if not a root file
      }

      const std::string numName{std::to_string(fileNum)};
      const std::string numNamePlus{std::to_string(fileNum + 2)};
      
      if (fs::is_regular_file( postTriggerSkimDir + datasetName + "/triggerSkim" + numNamePlus + ".root" )) {
	// don't overwrite existing skim files, except for the last two
	fileNum++;
	continue;
      }

      TChain datasetChain{"tree"};
      datasetChain.Add(path.c_str());
      TTree * const outTree = datasetChain.CloneTree(0);
      
      std::string outFilePath{ postTriggerSkimDir + datasetName + "/triggerSkim" + numName + ".root" };
      TFile outFile{outFilePath.c_str(), "RECREATE"};
	  
      const long long int numberOfEvents{datasetChain.GetEntries()};
      
      boost::progress_display progress{numberOfEvents, std::cout, outFilePath + "\n"};
      AnalysisEvent event{false, "", &datasetChain, true};

      for (long long int i{0}; i < numberOfEvents; i++) {
	++progress;  // update progress bar (++ must be prefix)
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
              //triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
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
              //triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
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
              //triggerDoubleCountCheck.emplace( std::make_pair(event.eventRun, event.eventNum) );
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
      std::cout << std::endl;

      delete outTree;
    }
  }

  if ( channel == "ee" || channel == "emu" )    std::cout << "Single electron trigger fired with double lepton trigger/Total single electron triggers fired: " << dupElectron << " / " << singleElectron << std::endl;
  if ( channel == "mumu" || channel == "emu" )  std::cout << "Single muon trigger fired with double lepton trigger/Total single muon triggers fired: " << dupMuon << " / " << singleMuon << std::endl;

}
