#include "TTree.h"
#include "TMVA/Timer.h"

#include "triggerScaleFactorsAlgo.hpp"
#include "config_parser.hpp"
#include "AnalysisEvent.hpp"

#include <iomanip>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>

TriggerScaleFactors::TriggerScaleFactors():
  config(""),
  plots(false),
  nEvents(0.),
  outFolder("plots/scaleFactors/"),
  postfix("default"),
  numFiles(-1),
  makePostLepTree(false),
  usePostLepTree(false),
  postLepSelTree_{nullptr},
  postLepSelTree2_{nullptr},
  postLepSelTree3_{nullptr},

  numberPassedElectrons(),
  numberTriggeredElectrons(),
  numberPassedMuons(),
  numberTriggeredMuons(),
  numberSelectedElectrons(),
  numberSelectedMuons(),
  numberSelectedElectronsTriggered(),
  numberSelectedMuonsTriggered()
{}

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

  if (!isMC){
    chain->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2",1);
    chain->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3",1);
  }
  else{
    chain->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1",1);
  }

  chain->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_v2",1);
  chain->SetBranchStatus("HLT_PFMET170_JetIdCleaned_v2",1);
  chain->SetBranchStatus("HLT_PFMET170_HBHECleaned_v2",1);
  chain->SetBranchStatus("HLT_PFHT350_PFMET100_v1",1);
  chain->SetBranchStatus("HLT_PFHT800_v2",1);
  chain->SetBranchStatus("HLT_MET250_v1",1);
  chain->SetBranchStatus("HLT_PFHT750_4JetPt50_v3",1);

}

void TriggerScaleFactors::show_usage(std::string name){
  std::cerr << "Usage: " << name << " <options>"
	    << "Options:\n"
	    << "\t-c  --config\tCONFIGURATION\tThe configuration file to be run over.\n"
            << "\t-n\t\t\t\tSet the number of events to run over. Leave blank for all.\n"
	    << "\t-o  --outFolder\tOUTFOLDER\tOutput folder for plots. If set overwrites what may be in the config file.\n"
	    << "\t-s  --postfix\tPOSTFIX\t\tPostfix for produced plots. Over-rides anything set in a configuration file.\n"
	    << "\t-d\t\t\t\tDump event info. For now this is the yield at each stage. May also include event lists later on. \n\t\t\t\t\tIf this flag is set all event weights are 1.\n"
	    << "\t-f  --nFiles \tNFILES\t\tUses a specific number of files to run over. \n\t\t\t\t\tThis is useful if testing stuff so that it doesn't have to access the T2 a lot etc.\n"
	    << std::endl;
}
		       

void TriggerScaleFactors::parseCommandLineArguements(int argc, char* argv[])
{
  gErrorIgnoreLevel = kInfo;
  //Set up environment a little.
  std::cout << std::setprecision(3) << std::fixed;
  // "This is the main function. It basically just loads a load of other stuff.";
  //Parse command line arguments - looking for config file.
  if (argc < 3){
    TriggerScaleFactors::show_usage(argv[0]);
    exit(1);
  }

  // Loop for parsing command line arguments.
  for (int i = 1; i < argc; ++i){
    std::string arg = argv[i];
    if ((arg=="-h") || (arg == "--help")){ // Display help stuff
      TriggerScaleFactors::show_usage(argv[0]);
      exit(0);
    }
    else if ((arg=="-c")||(arg=="--config")){ // Sets configuration file - Required!
      if (i + 1 < argc) {
      config = argv[++i];
      } else{
	std::cerr << "--config requires an argument!";
	exit(0);
      }
    }
    else if (arg=="-n") { // using this option sets the number of entries to run over.
      if (i+1 < argc){
	nEvents = atol(argv[++i]);
      }else{
	std::cerr << "-n requires a number of events to run over! You idiot!";
      }
    }
    else if (arg=="-p") {
      plots = true;
    }
    else if ((arg=="-o")||(arg=="--outFolder")){//Set output folder
      if (i + 1 < argc){
	outFolder = argv[++i];
      } else{
	std::cerr << "requires a string for output folder name.";
      }
    }
    else if ((arg=="-s")||(arg=="--postfix")){//Set plot postfix
      if (i + 1 < argc){
	postfix = argv[++i];
      } else{
	std::cerr << "requires a string for plot postfix.";
      }
    }
    else if (arg == "-g"){
      makePostLepTree = true;
    }
    else if (arg == "-u"){
      usePostLepTree = true;
    }
    else if (arg == "-f" || arg == "--nFiles"){
      if (i+1 < argc){
	numFiles = atoi(argv[++i]);
      }else{
	 std::cerr << "-f requires an int";
	 exit(0);
      }
    }

  } // End command line arguments loop.
  if (config == ""){
    std::cerr << "We need a configuration file! Type -h for usage. Error";
    exit(0);

  }
  //Some vectors that will be filled in the parsing.
  totalLumi = 0;
  lumiPtr = &totalLumi;
  if (!Parser::parse_config(config,&datasets,lumiPtr)){
    std::cerr << "There was an error parsing the config file.\n";
    exit(0);
  }
}

void TriggerScaleFactors::runMainAnalysis(){

  //Make pileupReweighting stuff here
  dataPileupFile = new TFile("pileup/truePileupTest.root","READ");
  dataPU = (TH1F*)(dataPileupFile->Get("pileup")->Clone());
  mcPileupFile = new TFile("pileup/pileupMC.root","READ");
  mcPU = (TH1F*)(mcPileupFile->Get("pileup")->Clone());

  //Get systematic files too.
  systUpFile = new TFile("pileup/truePileupUp.root","READ");
  pileupUpHist = (TH1F*)(systUpFile->Get("pileup")->Clone());
  systDownFile = new TFile("pileup/truePileupDown.root","READ");
  pileupDownHist = (TH1F*)(systDownFile->Get("pileup")->Clone());

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
    if (!usePostLepTree){
      if (!datasetFilled){
	if (!dataset->fillChain(datasetChain,numFiles)){
	  std::cerr << "There was a problem constructing the chain for " << dataset->name() << ". Continuing with next dataset.\n";
	  continue;
	}
	datasetFilled = true;
      }
    }

    else{
      std::string inputPostfix{};
      inputPostfix += postfix;
      std::cout << "skims/"+dataset->name()+inputPostfix + "SmallSkim.root" << std::endl;
      datasetChain->Add(("skims/"+dataset->name()+inputPostfix + "SmallSkim.root").c_str());
      std::ifstream secondTree{"skims/"+dataset->name()+inputPostfix + "SmallSkim1.root"};
      if (secondTree.good()) datasetChain->Add(("skims/"+dataset->name()+inputPostfix + "SmallSkim1.root").c_str());
      std::ifstream thirdTree{"skims/"+dataset->name()+inputPostfix + "SmallSkim2.root"};
      if (thirdTree.good()) datasetChain->Add(("skims/"+dataset->name()+inputPostfix + "SmallSkim2.root").c_str());

    }

    std::cout << "Trigger flag: " << dataset->getTriggerFlag() << std::endl;

    AnalysisEvent * event = new AnalysisEvent(dataset->isMC(),dataset->getTriggerFlag(),datasetChain, true);

    //Adding in some stuff here to make a skim file out of post lep sel stuff
    TFile * outFile1{nullptr};
    TTree * cloneTree{nullptr};

    TFile * outFile2{nullptr};
    TTree * cloneTree2{nullptr};

    TFile * outFile3{nullptr};
    TTree * cloneTree3{nullptr};

    if (makePostLepTree){
      outFile1 = new TFile{("skims/"+dataset->name() + postfix + "SmallSkim.root").c_str(),"RECREATE"};
      outFile2 = new TFile{("skims/"+dataset->name() + postfix + "SmallSkim1.root").c_str(),"RECREATE"};
      outFile3 = new TFile{("skims/"+dataset->name() + postfix + "SmallSkim2.root").c_str(),"RECREATE"};
      cloneTree = datasetChain->CloneTree(0);
      cloneTree->SetDirectory(outFile1);
      cloneTree2 = datasetChain->CloneTree(0);
      cloneTree2->SetDirectory(outFile2);
      cloneTree3 = datasetChain->CloneTree(0);
      cloneTree3->SetDirectory(outFile3);
      postLepSelTree_ = cloneTree;
      postLepSelTree2_ = cloneTree2;
      postLepSelTree3_ = cloneTree3;
    }

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
      
      //Does this event pass tight electron cut?
      //Create electron index
      event->electronIndexTight = getTightElectrons( event );
      bool passDoubleElectronSelection ( passDileptonSelection( event, event->electronIndexTight, true ) );
      //Does this event pass tight muon cut?
      //Create muon index
      event->muonIndexTight = getTightMuons( event );
      bool passDoubleMuonSelection ( passDileptonSelection( event, event->muonIndexTight, false ) );

      if ( (passDoubleElectronSelection || passDoubleMuonSelection) && makePostLepTree ) {
	if(postLepSelTree_) {
	  if (postLepSelTree_->GetEntriesFast() < 40000) postLepSelTree_->Fill();
	  else {
	    if (postLepSelTree2_->GetEntriesFast() < 40000) postLepSelTree2_->Fill();
	    else postLepSelTree3_->Fill();
	  }
	}
      }

      int triggerDoubleEG (0), triggerDoubleMuon (0), triggerMetDoubleEG (0), triggerMetDoubleMuon (0);
      int triggerMetElectronSelection (0), triggerMetMuonSelection (0);


      //Does event pass Double EG trigger and the electron selection?
      if ( passDoubleElectronSelection ) triggerDoubleEG 	= ( doubleElectronTriggerCut( event ) );
      if ( passDoubleElectronSelection ) triggerMetDoubleEG 	= ( doubleElectronTriggerCut( event )*metTriggerCut( event ) );
      //Does event pass Double Muon trigger and the muon selection?
      if ( passDoubleMuonSelection )  triggerDoubleMuon = ( doubleMuonTriggerCut( event ) );
      if ( passDoubleMuonSelection )  triggerMetDoubleMuon = ( doubleMuonTriggerCut( event )*metTriggerCut( event ) );
      //Does event pass either double lepton seletion and the MET triggers?
      if ( passDoubleElectronSelection ) triggerMetElectronSelection = ( metTriggerCut( event ) );
      if ( passDoubleMuonSelection ) triggerMetMuonSelection = ( metTriggerCut( event ) );

      if ( dataset->isMC() ) {
	numberPassedElectrons[0] += triggerMetElectronSelection*eventWeight; //Number of electrons passing the cross trigger and electron selection
	numberTriggeredElectrons[0] += triggerMetDoubleEG*eventWeight; //Number of electrons passing both cross trigger+electron selection AND double EG trigger
	numberPassedMuons[0] += triggerMetMuonSelection*eventWeight; //Number of muons passing the cross trigger and muon selection
	numberTriggeredMuons[0] += triggerMetDoubleMuon*eventWeight; //Number of muons passing both cross trigger+muon selection AND double muon trigger
	
	// Systematic stuff
	numberSelectedElectrons[0] += passDoubleElectronSelection*eventWeight;
	numberSelectedMuons[0] += passDoubleMuonSelection*eventWeight;
	numberSelectedElectronsTriggered[0] += triggerDoubleEG*eventWeight;;
	numberSelectedMuonsTriggered[0] += triggerDoubleMuon*eventWeight;;

      }
      else {
	numberPassedElectrons[1] += triggerMetElectronSelection*eventWeight; //Number of electrons passing the cross trigger and electron selection
	numberTriggeredElectrons[1] += triggerMetDoubleEG*eventWeight; //Number of electrons passing both cross trigger+electron selection AND double EG trigger
	numberPassedMuons[1] += triggerMetMuonSelection*eventWeight; //Number of muons passing the cross trigger and muon selection
	numberTriggeredMuons[1] += triggerMetDoubleMuon*eventWeight; //Number of muons passing both cross trigger+muon selection AND double muon trigger
      }

    }

    if (makePostLepTree){
      outFile1->cd();
      std::cout << "\nPrinting some info on the tree " <<dataset->name() << " " << cloneTree->GetEntries() << std::endl;
      std::cout << "But there were :" <<  datasetChain->GetEntries() << " entries in the original tree" << std::endl;
      cloneTree->Write();

      delete cloneTree;
      cloneTree = nullptr;
      outFile1->Write();
      outFile1->Close();

      //If we have any events in the second tree:
      if (cloneTree2->GetEntries() > 0){
	std::cout << "There are " << cloneTree2->GetEntries() << " entries in the second tree!" << std::endl;
	outFile2->cd();
	cloneTree2->Write();
	outFile2->Write();
      }
      if (cloneTree3->GetEntries() > 0){
	std::cout << "There are " << cloneTree3->GetEntries() << " entries in the third tree! What a lot of trees we've made." << std::endl;
	outFile3->cd();
	cloneTree3->Write();
	outFile3->Write();
      }
      delete cloneTree2;
      delete cloneTree3;
      cloneTree2 = nullptr;
      cloneTree3 = nullptr;
      outFile2->Close();
      outFile3->Close();
    }


    delete datasetChain;
  } //end dataset loop
}

std::vector<int> TriggerScaleFactors::getTightElectrons(AnalysisEvent* event) {
  std::vector<int> electrons;

  for (int i = 0; i < event->numElePF2PAT; i++){
    if (!event->elePF2PATIsGsf[i]) continue;
    TLorentzVector tempVec(event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]);

    if (tempVec.Pt() <= 20.0) continue;
    if (std::abs(tempVec.Eta()) >= 2.40) continue;

    // ECAL Gap
    if ( std::abs(event->elePF2PATSCEta[i]) > 1.4442 && std::abs(event->elePF2PATSCEta[i]) < 1.566 ) continue;

    // Spring15 Cut-based ID - Tight
    if (event->elePF2PATPhotonConversionTag[i]) continue;
      if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){
	  if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0101 ) continue;
	  if ( std::abs(event->elePF2PATDeltaEtaSC[i]) >= 0.0103 ) continue;
	  if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.0336 ) continue;
	  if ( event->elePF2PATHoverE[i] >= 0.0876 ) continue;
	  if ( event->elePF2PATComRelIsoRho[i] >= 0.0766 ) continue;
	  if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.0174 ) continue;
	  if ( std::abs(event->elePF2PATD0PV[i]) >= 0.0118 )continue;
	  if ( std::abs(event->elePF2PATDZPV[i]) >= 0.373 ) continue;
	  if ( event->elePF2PATMissingInnerLayers[i] > 2 ) continue;
	}
      else if ( std::abs(event->elePF2PATSCEta[i]) > 1.479 && std::abs(event->elePF2PATSCEta[i]) < 2.50 ){ // Endcap cut-based ID
	  if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0283 ) continue;
	  if ( std::abs(event->elePF2PATDeltaEtaSC[i]) >= 0.00733 ) continue;
	  if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.114 ) continue;
	  if ( event->elePF2PATHoverE[i] >= 0.0678 ) continue;
	  if ( event->elePF2PATComRelIsoRho[i] >= 0.0678 ) continue;
	  if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.0898 ) continue;
	  if ( std::abs(event->elePF2PATD0PV[i]) >= 0.0739 )continue;
	  if ( std::abs(event->elePF2PATDZPV[i]) >= 0.602 ) continue;
	  if ( event->elePF2PATMissingInnerLayers[i] > 1 ) continue;
	  }
      else continue;
      electrons.push_back(i);

    }
    return electrons;
}

std::vector<int> TriggerScaleFactors::getTightMuons(AnalysisEvent* event) {
  std::vector<int> muons;
  for (int i = 0; i < event->numMuonPF2PAT; i++){
    if (!event->muonPF2PATIsPFMuon[i]) continue;
    if (!event->muonPF2PATTrackID[i]) continue;
    if (!event->muonPF2PATGlobalID[i]) continue;

    if (event->muonPF2PATPt[i] <= 20.0) continue;
    if (std::abs(event->muonPF2PATEta[i]) >= 2.50) continue;
    if (event->muonPF2PATComRelIsodBeta[i] >= 0.15) continue;

    if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] >= 10.) continue;   
    if (event->muonPF2PATTkLysWithMeasurements[i] <= 5) continue;
    if (std::abs(event->muonPF2PATDBPV[i]) >= 0.2) continue;
    if (std::abs(event->muonPF2PATDZPV[i]) >= 0.5) continue;
    if (event->muonPF2PATMuonNHits[i] < 1) continue;
    if (event->muonPF2PATVldPixHits[i] < 1) continue;
    if (event->muonPF2PATMatchedStations[i] < 2) continue;

    muons.push_back(i);
  }
  return muons;
}

bool TriggerScaleFactors::passDileptonSelection( AnalysisEvent *event, std::vector<int> leptons, bool isElectron ){

  //Check if there are at least two electrons first. Otherwise use muons.
  
  float invMass (0.0);

  if (isElectron){
    for ( unsigned i = 0; i < leptons.size(); i++ ){
      for ( unsigned j = i + 1; j < leptons.size(); j++ ){
        if (event->elePF2PATCharge[leptons[i]] * event->elePF2PATCharge[leptons[j]] >= 0) continue; // check electron pair have correct charge.
        TLorentzVector lepton1 = TLorentzVector(event->elePF2PATGsfPx[leptons[i]],event->elePF2PATGsfPy[leptons[i]],event->elePF2PATGsfPz[leptons[i]],event->elePF2PATGsfE[leptons[i]]);
        TLorentzVector lepton2 = TLorentzVector(event->elePF2PATGsfPx[leptons[j]],event->elePF2PATGsfPy[leptons[j]],event->elePF2PATGsfPz[leptons[j]],event->elePF2PATGsfE[leptons[j]]);
        float candidateMass  = (lepton1 + lepton2).M();
	if (std::abs(candidateMass) > std::abs(invMass)){
        	event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
        	event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? leptons[i]:leptons[j];
        	event->zPairRelIso.first = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[leptons[i]]:event->elePF2PATComRelIsoRho[leptons[j]];
        	event->zPairRelIso.second = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[leptons[j]]:event->elePF2PATComRelIsoRho[leptons[i]];
        	event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
        	event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? leptons[j]:leptons[i];
		invMass = candidateMass;
      		}
	}
    } 
  }

  else {
    for ( unsigned i = 0; i < leptons.size(); i++ ){
      for ( unsigned j = i + 1; j < leptons.size(); j++ ){
	if (event->muonPF2PATCharge[leptons[i]] * event->muonPF2PATCharge[leptons[j]] >= 0) continue;
	TLorentzVector lepton1 = TLorentzVector(event->muonPF2PATPX[leptons[i]],event->muonPF2PATPY[leptons[i]],event->muonPF2PATPZ[leptons[i]],event->muonPF2PATE[leptons[i]]);
	TLorentzVector lepton2 = TLorentzVector(event->muonPF2PATPX[leptons[j]],event->muonPF2PATPY[leptons[j]],event->muonPF2PATPZ[leptons[j]],event->muonPF2PATE[leptons[j]]);
	float candidateMass = (lepton1 + lepton2).M() -91.1;
	if (std::abs(candidateMass) > std::abs(invMass)){
		event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
		event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? leptons[i]:leptons[j];
		event->zPairRelIso.first = event->muonPF2PATComRelIsodBeta[leptons[i]];
		event->zPairRelIso.second = event->muonPF2PATComRelIsodBeta[leptons[j]];
		event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
		event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? leptons[j]:leptons[i];
		invMass = candidateMass;
		}
      }
    }
  }
  if ( invMass > 20.0 ) return true;
  else return false;
}

bool TriggerScaleFactors::doubleElectronTriggerCut( AnalysisEvent* event ) {
  if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 || event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 || event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) return true;
  else return false;
}

bool TriggerScaleFactors::muonElectronTriggerCut( AnalysisEvent* event ) {
  if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2 > 0 || event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) return true;
  else return false;
}

bool TriggerScaleFactors::doubleMuonTriggerCut( AnalysisEvent* event ) {
  if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 || event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 ) return true;
  else return false;
}

bool TriggerScaleFactors::metTriggerCut( AnalysisEvent* event ) {
  if ( event->HLT_PFMET120_PFMHT120_IDTight_v2 > 0 || event->HLT_PFMET170_JetIdCleaned_v2 > 0 || event->HLT_PFMET170_HBHECleaned_v2 > 0 || event->HLT_PFHT800_v2 > 0 || event->HLT_MET250_v1 > 0 ) return true;
  else return false;
}

void TriggerScaleFactors::savePlots()
{
  // Calculate MC efficiency

  double electronEfficiencyMC 	= numberTriggeredElectrons[0]/(numberPassedElectrons[0]+1.0e-6);
  double muonEfficiencyMC 	= numberTriggeredMuons[0]/(numberPassedMuons[0]+1.0e-6);

  // Calculate Data efficiency

  double electronEfficiencyData	= numberTriggeredElectrons[1]/(numberPassedElectrons[1]+1.0e-6);
  double muonEfficiencyData 	= numberTriggeredMuons[1]/(numberPassedMuons[1]+1.0e-6);

  // Calculate SF

  double electronSF 	= electronEfficiencyData/(electronEfficiencyMC+1.0e-6);
  double muonSF 	= muonEfficiencyData/(muonEfficiencyMC+1.0e-6);

  // Calculate alphas
  double alphaElectron = ( (numberSelectedElectronsTriggered[0]/numberSelectedElectrons[0])*(numberPassedElectrons[0]/numberSelectedElectrons[0]) )/(numberTriggeredElectrons[0]/numberSelectedElectrons[0]+1.0e-6);
  double electronSystUncert = (1.0-alphaElectron);
  double alphaMuon = ( (numberSelectedMuonsTriggered[0]/numberSelectedMuons[0])*(numberPassedMuons[0]/numberSelectedMuons[0]) )/(numberTriggeredMuons[0]/numberSelectedMuons[0]+1.0e-6); 
  double muonSystUncert = (1.0-alphaMuon);

  // Save output

  std::cout << "electron/muon data efficiencies: " << electronEfficiencyData << " / " << muonEfficiencyData << std::endl;
  std::cout << "electron/muon MC efficiencies: " << electronEfficiencyMC << " / " << muonEfficiencyMC << std::endl;
  std::cout << "electron/muon trigger SFs: " << electronSF << " / " << muonSF << std::endl;

  TH2F* histEfficiencies = new TH2F("histEleEfficiencies", "Efficiencies and Scale Factors of dilepton triggers.; ; Channel.",3,0.0,3.0, 2,0.0,2.0);
  histEfficiencies->Fill(0.5,0.5, electronEfficiencyData);
  histEfficiencies->Fill(1.5,0.5, electronEfficiencyMC);
  histEfficiencies->Fill(2.5,0.5, electronSF);

  histEfficiencies->Fill(0.5,1.5, muonEfficiencyData);
  histEfficiencies->Fill(1.5,1.5, muonEfficiencyMC);
  histEfficiencies->Fill(2.5,1.5, muonSF);

  histEfficiencies->GetXaxis()->SetBinLabel(1,"Data #epsilon");
  histEfficiencies->GetXaxis()->SetBinLabel(2,"MC #epsilon");
  histEfficiencies->GetXaxis()->SetBinLabel(2,"SF");
  histEfficiencies->GetYaxis()->SetBinLabel(1,"ee trigger");
  histEfficiencies->GetYaxis()->SetBinLabel(2,"#mu#mu trigger");

  histEfficiencies->SetStats(kFALSE);

  mkdir( (outFolder+postfix).c_str(),0700);
  TFile* outFile = new TFile ( (outFolder+postfix+".root").c_str(), "RECREATE" );
  histEfficiencies->Write();
  outFile->Close();
}

