#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h> 
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"

#include "analysisAlgo.hpp"
#include "config_parser.hpp"
#include "AnalysisEvent.hpp"
#include "cutClass.hpp"
#include <libconfig.h++>
#include "histogramPlotter.hpp"

#include <iomanip>
#include <map>
#include <math.h>
#include <LHAPDF/LHAPDF.h>

double analysisAlgo::zptSF(TString channel, float zpt{

  double param1 = 0;
  double param2 = 0;
  double param3 = 0;

  if(channel == "mumumu"){
    //mumumu
    //1  p0           1.64891e+00   9.46744e-02   8.94955e-05 5.38436e-04
    //2  p1          -3.18363e-02   1.83020e-03   1.43202e-06 2.28732e-02
    //3  p2           1.96813e-01   1.80560e-02   2.59898e-05 7.20918e-03

    //1  p0           9.12190e-01   6.16931e-02   8.40864e-05 3.40620e-05
    //2  p1          -2.12648e-02   1.47982e-03   1.38426e-06 2.38949e-01
    //3  p2           2.32868e-01   2.61906e-02   3.56068e-05 8.38502e-03


    param1 = 9.12190e-01;
    param2 =-2.12648e-02;
    param3 = 2.32868e-01;

  }

  if(channel == "emumu"){
    //mumue
    //1  p0           1.08009e+00   2.20412e-01   1.81954e-04 -2.77029e-04
    //2  p1          -1.83319e-02   2.98128e-03   1.46500e-06 -2.16667e-02
    //3  p2          -3.79236e-03   2.77700e-02   1.94305e-05 -8.63841e-05

    //1  p0           5.88293e-01   5.43378e-02   6.56657e-05 -2.55726e-03
    //2  p1          -9.58817e-03   1.49703e-03   6.91871e-07 1.64841e-01
    //3  p2          -3.15588e-02   7.50287e-02   3.63099e-05 1.24242e-03


    param1 =  5.88293e-01;
    param2 = -9.58817e-03;
    param3 = -3.15588e-02;
  }

  if(channel == "eemu"){
    //eemu
    //1  p0           1.81997e+00   1.09691e-01   1.27075e-04 2.67625e-03
    //2  p1          -3.53330e-02   2.11348e-03   2.01050e-06 2.95414e-01
    //3  p2           2.00004e-01   1.93575e-02   3.33897e-05 1.34863e-02

    //1  p0           1.03732e+00   6.79924e-02   1.10651e-04 -4.52533e-02
    //2  p1          -2.11550e-02   1.34032e-03   1.62803e-06 -2.88549e+00
    //3  p2           1.52830e-01   2.17467e-02   4.20291e-05 -5.57304e-02


    param1 = 1.03732e+00;
    param2 =-2.11550e-02;
    param3 = 1.52830e-01;
  }


  if(channel == "eee"){
    //eee
    //1  p0           1.66655e+00   2.04856e-01   1.22417e-04 -8.87600e-06
    // 2  p1          -2.90064e-02   3.37196e-03   1.67677e-06 1.94266e-05
    //3  p2           1.12276e-01   2.87604e-02   2.89272e-05 -1.94049e-07

    //1  p0           8.23251e-01   8.60477e-02   6.95364e-05 3.23597e-03
    //2  p1          -1.74036e-02   2.04299e-03   1.02005e-06 2.12854e-01
    //3  p2           1.64031e-01   4.57851e-02   3.12269e-05 7.55832e-03
    param1 = 8.23251e-01;
    param2 = -1.74036e-02;
    param3 = 1.64031e-01;
  }

  // placeholder dilepton values
  if (channel == "mumu"){
    param1 = 9.12190e-01;
    param2 =-2.12648e-02;
    param3 = 2.32868e-01;
  }

  if (channel == "emu"){
    param1 = 1.03732e+00;
    param2 =-2.11550e-02;
    param3 = 1.52830e-01;
  }

  if (channel == "ee"){
    param1 = 8.23251e-01;
    param2 = -1.74036e-02;
    param3 = 1.64031e-01;
  }

  return  (exp(param1+param2*zpt) +param3 );
}





//This method is here to set up a load of branches in the TTrees that I will be analysing. Because it's vastly quicker to not load the whole damned thing.
void analysisAlgo::setBranchStatusAll(TTree * chain, bool isMC, std::string triggerFlag){
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
    //Jet variables
  chain->SetBranchStatus("numJetPF2PAT",1);
  chain->SetBranchStatus("jetPF2PATPx",1);
  chain->SetBranchStatus("jetPF2PATPy",1);
  chain->SetBranchStatus("jetPF2PATPz",1);
  chain->SetBranchStatus("jetPF2PATE",1);
  chain->SetBranchStatus("jetPF2PATEt",1);
  chain->SetBranchStatus("jetPF2PATPt",1);
  chain->SetBranchStatus("jetPF2PATPtRaw",1);
  chain->SetBranchStatus("jetPF2PATUnCorPt",1);
  chain->SetBranchStatus("jetPF2PATEta",1);
  chain->SetBranchStatus("jetPF2PATPhi",1);
  chain->SetBranchStatus("jetPF2PATNConstituents",1);
  chain->SetBranchStatus("jetPF2PAT*EnergyFractionCorr",1);
  chain->SetBranchStatus("jetPF2PAT*EnergyFraction",1);
  chain->SetBranchStatus("jetPF2PATChargedMultiplicity",1);
  chain->SetBranchStatus("jetPF2PATdRClosestLepton",1);
  //BTag
  chain->SetBranchStatus("jetPF2PATBDiscriminator",1);
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
    chain->SetBranchStatus("HLT_IsoMu20_v2",1);
    chain->SetBranchStatus("HLT_IsoMu20_eta2p1_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2",1);
    chain->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3",1);
    chain->SetBranchStatus("HLT_IsoMu18_v1",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3",1);
  }
  else{
    chain->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1",1);
    chain->SetBranchStatus("HLT_IsoMu20_v1",1);
    chain->SetBranchStatus("HLT_IsoMu20_eta2p1_v1",1);
    chain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1",1);
    chain->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1",1);
  }

}

static void analysisAlgo::show_usage(std::string name){
  std::cerr << "Usage: " << name << " <options>"
	    << "Options:\n"
	    << "\t-c  --config\tCONFIGURATION\tThe configuration file to be run over\n"
            << "\t--dilepton \t The dilepton channel is set to be run over instead of the default trilepton channel."
	    << "\t-p\t\t\t\tMake all plots. Currently segfaults if this isn't set, I believe.\n"
            << "\t-n\t\t\t\tSet the number of events to run over. Leave blank for all.\n"
            << "\t-l\t\t\t\tIf this option is set, scale MC plots to a fixed lumi. Default is lumi from data samples.\n"
	    << "\t-o  --outFolder\tOUTFOLDER\tOutput folder for plots. If set overwrites what may be in the config file.\n"
	    << "\t-s  --postfix\tPOSTFIX\t\tPostfix for produced plots. Over-rides anything set in a configuration file.\n"
	    << "\t-d\t\t\t\tDump event info. For now this is the yield at each stage. May also include event lists later on. \n\t\t\t\t\tIf this flag is set all event weights are 1.\n"
	    << "\t-x  --cutConf\tCUTCONF\t\tOverrides the cut configuration given in the usual configuration file.\n\t\t\t\t\tThis is mostly so that MC can be run on different cuts without having to make wqhole new confs.\n"
	    << "\t    --plotConf\tPLOTCONF\tOverrides the plot configuration file in the usual configuration file. \n\t\t\t\t\tFor various reasons I guess. Also sets -p flag automatically. If you don't want plots, DON'T USE THIS OPTION.\n"
	    << "\t-i\t\t\t\tInvert the isolation cut of the third lepton. This is for background estimation purposes. \n\t\t\t\t\tWho knows how I am supposed to use that though.\n"
	    << "\t-a  --synch\t\t\tMakes cutflows for synch exercise i.e. detailed lepSel cutflows. Doesn't do full event selection.\n"
	    << "\t-e\t\t\t\tGive a comma separated list of events to run on. This is for synch, but might be useful later?\n"
	    << "\t-f  --nFiles \tNFILES\t\tUses a specific number of files to run over. \n\t\t\t\t\tThis is useful if testing stuff so that it doesn't have to access the T2 a lot etc.\n"
	    << "\t-m\t\t\t\tMonte carlo only mode. Will not run over any data in the configuration.\n"
	    << "\t-b\t\t\t\tData only mode. Only runs over data, skips all MC.\n"
	    << "\t-t\t\t\t\tUse b-tagging efficiencies to reweight MC\n"
	    << "\t-y\t\t\t\tProduces a file of event dumps for stages of the synch.\n"
	    << "\t-g\t\t\t\tMakes post-lepSel tree\n"
	    << "\t-u\t\t\t\tUses post-lepSel trees\n"
	    << "\t-z  --makeMVATree\t\tProduce a tree after event selection for MVA purposes\n"
	    << "\t-v  --syst  \tSYST\t\tDo the desired systematic. Brief workaround here, not final yet\n"
	    << "\t-j\t\t\t\tMake b-tagging efficiency histograms. Probably doesn't need to be run too many times.\n"
	    << "\t-k          \tCHANS\t\tBit mask dealy for the channels. 1 - eee 2 - eemu 4 - emumu 8 - mumumu 16 through 128 are same but for inverted third lep iso.\n"
	    << "\t    --skipTrig\t\t\tSkip running triggers. Used for trigger studies or something.\n"
	    << "\t    --mvaDir \tDIR\t\tChange the name of the folder the mva outputs go to. mvaTest/ by default. Include the /.\n"
	    << "\t    --jetRegion \t\tSet the jet region that the analysis will look at. Takes arguments NJETS,NBJETS,MAXJETS,MAXBJETS.\n"
	    << "\t    --metCut \tCUT\t\tAlter the MET cut of the analysis. 0 by default.\n"
	    << "\t    --mtwCut \tCUT\t\tAlter the mTW cut of the analysis. 0 by default.\n"
	    << "\t-h  --help\t\t\tShow this help message\n"
	    << std::endl;
}
		       

void analysisAlgo::parseCommandLineArguements(int argc, char* argv[])
{
  // Loop for parsing command line arguments.
  for (int i = 1; i < argc; ++i){
    std::string arg = argv[i];
    if ((arg=="-h") || (arg == "--help")){ // Display help stuff
      show_usage(argv[0]);
      return 0;
    }
    else if ((arg=="-c")||(arg=="--config")){ // Sets configuration file - Required!
      if (i + 1 < argc) {
      config = argv[++i];
      } else{
	std::cerr << "--config requires an argument!";
	return 0;
      }
    }
    else if (arg=="--dilepton"){ // Sets whether trilepton or dilepton channel is to be analysed.
      trileptonChannel = false;
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
    else if (arg=="-l"){
      if (i+1 < argc){
	usePreLumi = atof(argv[++i]);
      }else{
	std::cerr << "-l requries a float!";
	return 0;
      }
    }
    else if (arg == "-d"){//Option for dumping event information about selected events.
      infoDump = true;
    }
    else if (arg == "-x" || arg == "--cutConf"){
      if (i+1 < argc){
	*cutConfName = argv[++i];
      }else{
	std::cerr <<" -x requires an argument";
	return 0;
      }
    }
    else if (arg == "--plotConf"){
      if (i+1 < argc){
	*plotConfName = argv[++i];
	plots = true;
      }else{
	std::cerr <<" --plotConf requires an argument";
	return 0;
      }
    }
    else if (arg == "-i"){ //Set up anti-isolation cut for 
      invertIsoCut = true;
    }
    else if (arg == "-a" || arg == "--synch"){ // Change to synch exercise cut flow.
      synchCutFlow = true;
    }
    else if (arg == "-m"){
      skipData = true;
    }
    else if (arg == "-b"){
      skipMC = true;
    }
    else if (arg == "-y"){
      dumpEventNumbers = true;
    }
    else if (arg == "-t"){
      usebTagWeight = true;
    }
    else if (arg == "-f" || arg == "--nFiles"){
      if (i+1 < argc){
	numFiles = atoi(argv[++i]);
      }else{
	 std::cerr << "-f requires an int";
	 return 0;
      }
    }
    else if (arg == "-e"){
      readEventList = true;
      std::stringstream ss(argv[++i]);
      std::string item;
      while (std::getline(ss,item,',')){
	jetRegVars.push_back(atoi(item.c_str()));
      }
    }
    else if (arg == "-g"){
      makePostLepTree = true;
    }
    else if (arg == "-z" || arg == "--makeMVATree"){
      makeMVATree = true;
    }
    else if (arg == "-u"){
      usePostLepTree = true;
    }
    else if (arg == "-v" || arg == "--syst"){
      if (i+1 < argc){
	systToRun = atoi(argv[++i]);
      }
      else{
	std::cerr << "-v requires an int";
	return 0;
      }
    }
    else if (arg == "-j"){
      makeBTagEffPlots = true;
    }
    else if (arg == "-k"){
      if (i+1 < argc){
	channelsToRun = atoi(argv[++i]);
      }
      else{
	std::cerr << "-k needs a config file!";
	return 0;
      }
    }
    else if (arg == "--skipTrig"){
      skipTrig = true;
      std::cout << "Note that you're skipping the trigger!" << std::endl;
    }
    else if (arg == "--mvaDir"){
      mvaDir = argv[++i];
    }
    else if (arg == "--jetRegion"){
      customJetRegion = true;
      std::stringstream ss(argv[++i]);
      std::string item;
      while (std::getline(ss,item,',')){
	jetRegVars.push_back(atoi(item.c_str()));
      }
      std::cout << "CAUTION! Using a custom jet region of "<< jetRegVars[0] << "-" << jetRegVars[2] << " jets, and " << jetRegVars[1] << "-" << jetRegVars[3] << " b-jets" <<std::endl;

    }
    else if (arg == "--metCut"){
      metCut = atof(argv[++i]);
      std::cout << "Non zero MET cut! Applied " << metCut << " cut." << std::endl;
    }
    else if (arg == "--mtwCut"){
      mtwCut = atof(argv[++i]);
      std::cout << "Non zero mTW cut! Applied " << mtwCut << " cut." << std::endl;
    }
    
  } // End command line arguments loop.
  if (config == ""){
    std::cerr << "We need a configuration file! Type -h for usage. Error";
    return 0;
  }
  if (usebTagWeight && !usePostLepTree){
    std::cerr << "At the moment only getting btag weights from post lep-sel trees is supported. Sorry everyone.";
    return 0;
  }
  if (usebTagWeight && makeBTagEffPlots){
    std::cerr << "I doubt that set of options is a good idea, so I'm going to just quietly exit here. Don't generate and use b-tag efficiency numbers at the same time...";
    return 0;
  }
  //Make some vectors that will be filled in the parsing.
  std::vector<Dataset> datasets;
  double totalLumi = 0;
  double* lumiPtr = &totalLumi;
  if (!Parser::parse_config(config,&datasets,lumiPtr,&plotNames,&xMin,&xMax,&nBins,&fillExp,&xAxisLabels,&cutStage,cutConfName,plotConfName,&outFolder,&postfix,&channel)){
    std::cerr << "There was an error parsing the config file.\n";
    return 0;
  }

  if (channelsToRun && trileptonChannel){
    std::cout << "Running over the channels: " << std::endl;
    for (unsigned int channelInd = 1; channelInd != 256; channelInd = channelInd << 1){
      if (!(channelInd & channelsToRun) && channelsToRun) continue;
      if (channelInd & 17){
	std::cout << "eee ";
      }
      if (channelInd & 34){ //eemu channels
	std::cout << "eemu ";
      }
      if (channelInd & 68){ // emumu channels
	std::cout << "emumu ";
      }
      if (channelInd & 136){ // mumumu channels
	std::cout << "mumumu ";
      }
      if (channelInd & 15){ //nominal samples
	std::cout << "nominal" << std::endl;
      }
      if (channelInd & 240){ //inv iso samples
	std::cout << "inverted" << std::endl;
      }
    }
  }

  if (channelsToRun && !trileptonChannel){
    std::cout << "Running over the channels: " << std::endl;
    for (unsigned int channelInd = 1; channelInd != 64; channelInd = channelInd << 1){
      if (!(channelInd & channelsToRun) && channelsToRun) continue;
      if (channelInd & 9){
	std::cout << "ee ";
      }
      if (channelInd & 18){ //emu channels
	std::cout << "emu ";
      }
      if (channelInd & 40){ // mumu channels
	std::cout << "mumu ";
      }
      if (channelInd & 7){ //nominal samples
	std::cout << "nominal" << std::endl;
      }
      if (channelInd & 56){ //inv iso samples
	std::cout << "inverted" << std::endl;
      }
    }
  }

}

void analysisAlgo::setupSystematics()
{
  systNames.push_back("");
  systNames.push_back("__trig__plus");
  systNames.push_back("__trig__minus");
  systNames.push_back("__jer__plus");
  systNames.push_back("__jer__minus");
  systNames.push_back("__jes__plus");
  systNames.push_back("__jes__minus");
  systNames.push_back("__pileup__plus");
  systNames.push_back("__pileup__minus");
  systNames.push_back("__bTag__plus");
  systNames.push_back("__bTag__minus");
  systNames.push_back("__pdf__plus");
  systNames.push_back("__pdf__minus");

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

  puReweight = (TH1F*)dataPU->Clone();
  puReweight->Scale(1.0/puReweight->Integral());
  mcPU->Scale(1.0/mcPU->Integral());
  puReweight->Divide(mcPU);
  puReweight->SetDirectory(0);

  /// And do the same for systematic sampl
  puSystUp = (TH1F*)pileupUpHist->Clone();
  puSystUp->Scale(1.0/puSystUp->Integral());
  puSystUp->Divide(mcPU);
  puSystUp->SetDirectory(0);
  puSystDown = (TH1F*)pileupDownHist->Clone();
  puSystDown->Scale(1.0/puSystDown->Integral());
  puSystDown->Divide(mcPU);
  puSystDown->SetDirectory(0);

  dataPileupFile->Close();
  mcPileupFile->Close();
  systUpFile->Close();
  systDownFile->Close();

  //Initialise PDFs
  if (systToRun & 1024 || systToRun & 2048){
    LHAPDF::initPDFSet(1, "CT10nnlo.LHgrid");
    //    LHAPDF::initPDFSet(1, "cteq6ll.LHpdf");
    //    LHAPDF::initPDFSet(1, "cteq6lg.LHgrid");
  }
  //    LHAPDF::initPDFSet(1, "CT10nnlo.LHgrid");


}

void analysisAlgo::setupCuts()
{
  //Make cuts object. The methods in it should perhaps just be i nthe AnalysisEvent class....
  Cuts * cutObj = new Cuts(plots,plots||infoDump,invertIsoCut,synchCutFlow,dumpEventNumbers);
  if (!cutObj->parse_config(*cutConfName)){
    std::cerr << "There was a problem with parsing the config!" << std::endl;
    return 0;
  };
  //For studying some trigger things. Default is false.
  cutObj->setSkipTrig(skipTrig);
  if (customJetRegion) cutObj->setJetRegion(jetRegVars[0],jetRegVars[1],jetRegVars[2],jetRegVars[3]);
  cutObj->setMetCut(metCut);
  cutObj->setMTWCut(mtwCut);
}

void analysisAlgo::setupPlots();
{
  //Do a little initialisation for the plots here. Will later on be done in a config file.
  
  //Initialise plot stage names.
  stageNames = {"lepSel","zMass","jetSel","bTag"};
}

bool datasetFilled = false;

if (totalLumi == 0.) totalLumi = usePreLumi;
std::cout << "Using lumi: " << totalLumi << std::endl;
for (std::vector<Dataset>::iterator dataset = datasets.begin(); dataset!=datasets.end(); ++dataset){
  datasetFilled = false;
  TChain * datasetChain = new TChain(dataset->treeName().c_str());
  uint channelIndMax = 256;
  if ( !trileptonChannel ){ channelIndMax = 64; }
  for (unsigned int channelInd = 1; channelInd != channelIndMax; channelInd = channelInd << 1){
    std::string chanName = "";
    if (!(channelInd & channelsToRun) && channelsToRun) continue;
    if (channelsToRun && trileptonChannel){
      if (channelInd & 17){ // eee channels
	cutObj->setNumLeps(0,0,3,3);
	cutObj->setCutConfTrigLabel("e");
	channel = "eee";
	postfix = "eee";
	chanName += "eee";
      }
      if (channelInd & 34){ //eemu channels
	cutObj->setNumLeps(1,1,2,2);
	cutObj->setCutConfTrigLabel("d");
	channel = "eemu";
	postfix = "eemu";
	chanName += "eemu";
      }
      if (channelInd & 68){ // emumu channels
	cutObj->setNumLeps(2,2,1,1);
	cutObj->setCutConfTrigLabel("d");
	channel = "emumu";
	postfix = "emumu";
	chanName += "emumu";
      }
      if (channelInd & 136){ // mumumu channels
	cutObj->setNumLeps(3,3,0,0);
	cutObj->setCutConfTrigLabel("m");
	channel = "mumumu";
	postfix = "mumumu";
	chanName += "mumumu";
      }
      if (channelInd & 15){ //nominal samples
	cutObj->setInvIsoCut(false);
	invertIsoCut = false;
	chanName += "nom";
      }
      if (channelInd & 240){ //inv iso samples
	cutObj->setInvIsoCut(true);
	invertIsoCut = true;
	chanName += "inv";
      }
    } 
    if (channelsToRun && !trileptonChannel){
      if (channelInd & 9){ // ee channels
	cutObj->setNumLeps(0,0,2,2);
	cutObj->setCutConfTrigLabel("e");
	channel = "ee";
	postfix = "ee";
	chanName += "ee";
      }
      if (channelInd & 18){ //emu channels
	cutObj->setNumLeps(1,1,1,1);
	cutObj->setCutConfTrigLabel("d");
	channel = "emu";
	postfix = "emu";
	chanName += "emu";
      }
      if (channelInd & 36){ // mumu channels
	cutObj->setNumLeps(2,2,0,0);
	cutObj->setCutConfTrigLabel("m");
	channel = "mumu";
	postfix = "mumu";
	chanName += "mumu";
      }
      if (channelInd & 7){ //nominal samples
	cutObj->setInvIsoCut(false);
	invertIsoCut = false;
	chanName += "nom";
      }
      if (channelInd & 56){ //inv iso samples
	cutObj->setInvIsoCut(true);
	invertIsoCut = true;
	chanName += "inv";
      }
    } 
    if (dataset->isMC() && skipMC) continue;
    if (!dataset->isMC() && skipData) continue;
    if (plots||infoDump) { // Initialise a load of stuff that's required by the plotting macro.
      int systMask = 1;
      for (unsigned int systInd = 0; systInd < systNames.size(); systInd++){
	if (systInd > 0 && !(systToRun & systMask)){
	  systMask = systMask << 1;
	  continue;
	} 
	if (cutFlowMap.find(dataset->getFillHisto()+systNames[systInd]) == cutFlowMap.end()){
	  cutFlowMap[dataset->getFillHisto()] = new TH1F((dataset->getFillHisto()+systNames[systInd]+"cutFlow").c_str(),(dataset->getFillHisto()+systNames[systInd]+"cutFlow").c_str(),4,0,4); //Hopefully make this configurable later on. Same deal as the rest of the plots I guess, work out libconfig.
	  if (systInd == 0 && datasetInfos.find(dataset->getFillHisto()) == datasetInfos.end()){
	    legOrder.push_back(dataset->getFillHisto());
	    plotOrder.push_back(dataset->getFillHisto());
	    datasetInfos[dataset->getFillHisto()] = datasetInfo();
	    datasetInfos[dataset->getFillHisto()].colour = dataset->getColour();
	    datasetInfos[dataset->getFillHisto()].legLabel = dataset->getPlotLabel();
	    datasetInfos[dataset->getFillHisto()].legType = dataset->getPlotType();
	  }
	  if (plots){ // Only make all the plots if it's entirely necessary.
	    std::cout << "Made plots under" << (dataset->getFillHisto()+systNames[systInd]+channel).c_str() << std::endl; 
	    if (plotsMap.find(channel) == plotsMap.end()){
	      plotsVec.push_back(systNames[systInd]+channel);
	    }
	    plotsMap[systNames[systInd]+channel][(dataset->getFillHisto()).c_str()] = std::map<std::string,Plots*>();
	    for (unsigned int j = 0; j < stageNames.size(); j++){
	      plotsMap[systNames[systInd]+channel][(dataset->getFillHisto()).c_str()][stageNames[j]] = new Plots(plotNames, xMin, xMax,nBins, fillExp, xAxisLabels, cutStage, j, dataset->getFillHisto()+"_"+stageNames[j]+systNames[systInd]+"_"+channel);
	    }
	  }
	}//end cutFlow find loop
	if (systInd > 0) systMask = systMask << 1;
      }//end systematic loop
      
    } //end plots if
      //If making either plots or doing the event dump, make cut flow object.
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
      std::string inputPostfix = "";
      inputPostfix += postfix;
      inputPostfix += invertIsoCut?"invIso":"";
      std::cout << "skims/"+dataset->name()+inputPostfix + "SmallSkim.root" << std::endl;
      datasetChain->Add(("skims/"+dataset->name()+inputPostfix + "SmallSkim.root").c_str());
      std::ifstream secondTree(("skims/"+dataset->name()+inputPostfix + "SmallSkim1.root").c_str());
      if (secondTree.good()) datasetChain->Add(("skims/"+dataset->name()+inputPostfix + "SmallSkim1.root").c_str());
      std::ifstream thirdTree(("skims/"+dataset->name()+inputPostfix + "SmallSkim2.root").c_str());
      if (thirdTree.good()) datasetChain->Add(("skims/"+dataset->name()+inputPostfix + "SmallSkim2.root").c_str());
    }
    cutObj->setMC(dataset->isMC());
    cutObj->setEventInfoFlag(readEventList);
    cutObj->setTriggerFlag(dataset->getTriggerFlag());
    std::cout << "Trigger flag: " << dataset->getTriggerFlag() << std::endl;

    //Here we will initialise the b-tag eff plots if we are doing b-tag efficiencies
    std::vector<TH2D*> bTagEffPlots;
    std::vector<std::string> denomNum {"Denom","Num"};
    std::vector<std::string> typesOfEff {"b","c","uds","g"};
    if (makeBTagEffPlots && dataset->isMC()){
      int ptBins = 4, etaBins = 4;
      float ptMin = 0., ptMax = 200., etaMin = 0., etaMax = 2.4;
      for (int unsigned denNum = 0; denNum < denomNum.size(); denNum++){
	for (int unsigned type = 0; type < typesOfEff.size(); type++){
	  bTagEffPlots.push_back(new TH2D(("bTagEff_"+denomNum[denNum]+"_"+typesOfEff[type]).c_str(),("bTagEff_"+denomNum[denNum]+"_"+typesOfEff[type]).c_str(),ptBins,ptMin,ptMax,etaBins,etaMin,etaMax));
	}
      }
      cutObj->setBTagPlots(bTagEffPlots,true);
    }//end btag eff plots.
    if (usePostLepTree && usebTagWeight && dataset->isMC()){
      //Get efficiency plots from the file. Will have to be from post-lep sel trees I guess.
      std::string inputPostfix = "";
      inputPostfix += postfix;
      inputPostfix += invertIsoCut?"invIso":"";
      TFile * datasetFileForHists = new TFile(("skims/"+dataset->name() + inputPostfix + "SmallSkim.root").c_str(), "READ");
      for (int unsigned denNum = 0; denNum < denomNum.size(); denNum++){
	for (int unsigned eff = 0; eff < typesOfEff.size(); eff++){
	  bTagEffPlots.push_back((TH2D*)(datasetFileForHists->Get(("bTagEff_"+denomNum[denNum]+"_"+typesOfEff[eff]).c_str())->Clone()));
	}
      }
      for (int unsigned plotIt = 0; plotIt < bTagEffPlots.size(); plotIt++){
	bTagEffPlots[plotIt]->SetDirectory(0);
      }
      cutObj->setBTagPlots(bTagEffPlots,false);
      datasetFileForHists->Close();
    }
    //extract the dataset weight.
    float datasetWeight = dataset->getDatasetWeight(totalLumi);


    //Apply trigger SF here. Also does systematic for trigger +-
    if (infoDump) datasetWeight = 1;
    std::cout << datasetChain->GetEntries() << " number of items in tree. Dataset weight: " << datasetWeight << std::endl;
    AnalysisEvent * event = new AnalysisEvent(dataset->isMC(),dataset->getTriggerFlag(),datasetChain);

    //Adding in some stuff here to make a skim file out of post lep sel stuff
    TTree * cloneTree = 0;
    TTree * cloneTree2 = 0;
    TTree * cloneTree3 = 0;
    if (makePostLepTree){
      cloneTree = datasetChain->CloneTree(0);
      cloneTree2 = datasetChain->CloneTree(0);
      cloneTree3 = datasetChain->CloneTree(0);
      cutObj->setCloneTree(cloneTree,cloneTree2,cloneTree3);
    }
    //If we're making the MVA tree, set it up here. 
    std::vector<TTree *> mvaTree;
    //Add a few variables into the MVA tree for easy access of stuff like lepton index etc
    float eventWeight = 0;
    int zLep1Index = -1; // Addresses in elePF2PATWhatever of the z lepton
    int zLep2Index = -1;
    int wLepIndex = -1;
    int jetInd[15];  // The index of the selected jets;
    int bJetInd[10]; // Index of selected b-jets;
    //Now add in the branches:
    
    if (makeMVATree){
      int systMask = 1;
      std::cout << "Making systematic trees for " << dataset->name() << ": ";
      for (unsigned int systIn = 0; systIn < systNames.size(); systIn++){
	std::cout << systNames[systIn] << " ";
	//	std::cout << "Making systs: " << systMask << " " << systToRun << " " << systIn << " " << (systMask & systToRun) << std::endl;
	/*	if (systIn > 0 && !(systMask & systToRun)){
		if (systIn > 0) systMask = systMask << 1;
		continue;
		}*/
	mvaTree.push_back(datasetChain->CloneTree(0));
	mvaTree[systIn]->SetName((mvaTree[systIn]->GetName()+systNames[systIn]).c_str());
	mvaTree[systIn]->Branch("eventWeight", &eventWeight, "eventWeight/F");
	mvaTree[systIn]->Branch("zLep1Index",&zLep1Index,"zLep1Index/I");
	mvaTree[systIn]->Branch("zLep2Index",&zLep2Index,"zLep2Index/I");
	mvaTree[systIn]->Branch("wLepIndex",&wLepIndex,"wLepIndex/I");
	mvaTree[systIn]->Branch("jetInd",jetInd,"jetInd[15]/I");
	mvaTree[systIn]->Branch("bJetInd",bJetInd,"jetInd[10]/I");

	if (systIn > 0) systMask = systMask << 1;
      }
      std::cout <<std::endl;
    }
    /*    else{
	  event->fChain->SetBranchStatus("*",0); //Should disable most branches.
	  setBranchStatusAll(event->fChain,dataset->isMC(),dataset->getTriggerFlag());
	  }*/

    int numberOfEvents = datasetChain->GetEntries();
    if (nEvents && nEvents < numberOfEvents) numberOfEvents = nEvents;
    //    datasetChain->Draw("numElePF2PAT","numMuonPF2PAT > 2");
    //    TH1F * htemp = (TH1F*)gPad->GetPrimitive("htemp");
    //    htemp->SaveAs("tempCanvas.png");
    int foundEvents = 0;
    for (int i = 0; i < numberOfEvents; i++){
      if (i % 500 < 0.01) std::cerr << i << " (" << 100*float(i)/numberOfEvents << "%) with " << event->numElePF2PAT << " electrons. Found " << (synchCutFlow?cutObj->numFound():foundEvents) << " events.\r";
      event->GetEntry(i);
      //Do the systematics indicated by the systematic flag, oooor just do data if that's your thing. Whatevs.
      int systMask = 1;
      for (unsigned int systInd = 0; systInd < systNames.size(); systInd++){
	if (!dataset->isMC() && systInd > 0) break;
	//	std::cout << systInd << " " << systMask << std::endl;
	if (systInd > 0 && !(systMask & systToRun)) {
	  if (systInd > 0) systMask = systMask << 1;
	  continue;
	}
	eventWeight = 1;
	//apply trigger weights here.
	if (dataset->isMC()){
	  float pileupWeight = puReweight->GetBinContent(puReweight->GetXaxis()->FindBin(event->numVert));
	  if (systMask == 64) pileupWeight = puSystUp->GetBinContent(puSystUp->GetXaxis()->FindBin(event->numVert));
	  if (systMask == 128) pileupWeight = puSystDown->GetBinContent(puSystDown->GetXaxis()->FindBin(event->numVert));
	  eventWeight *= pileupWeight;
	  // trilepton stuff
	  if (channel == "eee"){
	    float twgt = 0.987;
	    if (systInd > 0 && (systMask == 1)) twgt += 0.036;
	    if (systInd > 0 && (systMask == 2)) twgt -= 0.036;
	    eventWeight *= twgt;
	  }
	  else if (channel == "eemu"){
	    float twgt = 0.987;
	    if (systInd > 0 && (systMask == 1)) twgt += 0.035;
	    if (systInd > 0 && (systMask == 2)) twgt -= 0.035;
	    eventWeight *= twgt;
	  }
	  if (channel == "emumu"){
	    float twgt = 0.886;
	    if (systInd > 0 && (systMask == 1)) twgt += 0.042;
	    if (systInd > 0 && (systMask == 2)) twgt -= 0.042;
	    eventWeight *= twgt;
	  }
	  if (channel == "mumumu"){
	    float twgt = 0.9871;
	    if (systInd > 0 && (systMask == 1)) twgt += 0.0242;
	    if (systInd > 0 && (systMask == 2)) twgt -= 0.0212;
	    eventWeight *= twgt;
	  }
	  // dilepton stuff, placeholder for now
	  if (channel == "ee"){
	    float twgt = 0.987;
	    if (systInd > 0 && (systMask == 1)) twgt += 0.036;
	    if (systInd > 0 && (systMask == 2)) twgt -= 0.036;
	    eventWeight *= twgt;
	  }
	  else if (channel == "emu"){
	    float twgt = 0.987;
	    if (systInd > 0 && (systMask == 1)) twgt += 0.035;
	    if (systInd > 0 && (systMask == 2)) twgt -= 0.035;
	    eventWeight *= twgt;
	  }
	  if (channel == "mumu"){
	    float twgt = 0.886;
	    if (systInd > 0 && (systMask == 1)) twgt += 0.042;
	    if (systInd > 0 && (systMask == 2)) twgt -= 0.042;
	    eventWeight *= twgt;
	  }
	}
	if (infoDump) eventWeight = 1;
	if (readEventList) {
	  bool tempBool = false;
	  for (unsigned int i = 0; i < eventNumbers.size(); i++){
	    if (eventNumbers[i] == event->eventNum) {
	      tempBool = true;
	      break;
	    }
	  }
	  if (!tempBool) continue;
	  std::cout << event->eventNum << " " << event->eventRun << " " << event->eventLumiblock << " " << datasetChain->GetFile()->GetName() << std::endl;
	  cutObj->dumpLooseLepInfo(event);
	  cutObj->dumpLeptonInfo(event);
	
	}
	eventWeight*=datasetWeight;
	if (!cutObj->makeCuts(event,&eventWeight,plotsMap[systNames[systInd]+channel][dataset->getFillHisto()],cutFlowMap[dataset->getFillHisto()+systNames[systInd]],systInd?systMask:systInd)) {
	  if (systInd) systMask = systMask << 1;
	  continue;
	}
	//	  std::cout << std::setprecision(9) << eventWeight << " " << datasetWeight << std::endl;
	//Do PDF reweighting things here
	if (systMask == 1024 || systMask == 2048){
	  //std::cout << std::setprecision(15) << eventWeight << " ";
	  LHAPDF::usePDFMember(1,0);
	  float q = event->genPDFScale;
	  float x1 = event->genPDFx1;
	  float x2 = event->genPDFx2;
	  int id1 = event->genPDFf1;
	  int id2 = event->genPDFf2;
	  if (id2 == 21) id2 = 0;
	  if (id1 == 21) id1 = 0;
	  double xpdf1 = LHAPDF::xfx(1, x1, q, id1);
	  double xpdf2 = LHAPDF::xfx(1, x2, q, id2);
	  std::vector<float> pdf_weights;
	  //std::cout << q << " " << x1 << " " << x2 << " " << id1 << " " << id2 << " ";
	  //std::cout << xpdf1 << " " << xpdf2 << " " << xpdf1 * xpdf2 << " ";
	  float min = 1.0;
	  float max = 1.0;
	  float pdfWeightUp = 0.0;
	  float pdfWeightDown = 0.0;
	  for (int i = 1; i <=50; ++i){
	    LHAPDF::usePDFMember(1,i);
	    double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
	    double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
	    //std::cout << " " << x1 << " " << id1 << " " << x2 << " " << id2 << " " << q << " " <<xpdf1 << " " << xpdf2 << " " << xpdf1_new << " " << xpdf2_new << " ";
	    double weight = 1.;
	    if( (xpdf1 * xpdf2) > 0.00001)
	      weight = xpdf1_new * xpdf2_new / (xpdf1 * xpdf2);
	    pdf_weights.push_back(weight);
	    if (weight > 1.0) pdfWeightUp += (1-weight) * (1-weight);
	    if (weight < 1.0) pdfWeightDown += (1-weight) * (1-weight);
	    if (weight > max) max = weight;
	    if (weight < min) min = weight;
	    //	      std::cout << " " << xpdf1_new << " " << xpdf2_new << " " << weight << " ";
			
	  }
	  if (systMask == 1024) eventWeight *= max;
	  if (systMask == 2048) eventWeight *= min;
	  //std::cout << eventWeight << std::setprecision(4) << max << " " << min << " " << 1+std::sqrt(pdfWeightUp) << " " << 1-std::sqrt(pdfWeightDown) << std::endl;
	  //std::cout << std::setprecision(9) << " " << min << " " << max << " " << eventWeight << std::endl;
	}
	//      if (synchCutFlow){
	//	std::cout << event->eventNum << " " << event->eventRun << " " << event->eventLumiblock << " " << std::endl;
	//}
	//Do the Zpt reweighting here
	if (invertIsoCut){
	  float zPT = (event->zPairLeptons.first+event->zPairLeptons.second).Pt();
	  eventWeight *= zptSF(channel,zPT);
	}
	if (makeMVATree){
	  zLep1Index = event->zPairIndex.first;
	  zLep2Index = event->zPairIndex.second;
	  wLepIndex = event->wLepIndex;
	  for (unsigned int jetIndexIt = 0; jetIndexIt < 15; jetIndexIt++){
	    if (jetIndexIt < event->jetIndex.size()) jetInd[jetIndexIt] = event->jetIndex[jetIndexIt];
	    else jetInd[jetIndexIt] = -1;
	  }
	  for (unsigned int bJetIt = 0; bJetIt < 10; bJetIt++){
	    if (bJetIt < event->bTagIndex.size()) bJetInd[bJetIt] = event->bTagIndex[bJetIt];
	    else bJetInd[bJetIt] = -1;
	  }
	  mvaTree[systInd]->Fill();
	  
	}

	foundEvents++;
	if (systInd > 0) systMask = systMask << 1;
      }// End systematics loop.
    } //end event loop

      //If we're making post lepSel skims save the tree here
    if (makePostLepTree){
      TFile outFile(("skims/"+dataset->name() + postfix + (invertIsoCut?"invIso":"") + "SmallSkim.root").c_str(),"RECREATE");
      outFile.cd();
      std::cout << "\nPrinting some info on the tree " <<dataset->name() << " " << cloneTree->GetEntries() << std::endl;
      std::cout << "But there were :" <<  datasetChain->GetEntries() << " entries in the original tree" << std::endl;
      cloneTree->Write();
      //If we're doing b-tag efficiencies, let's save them here.
      if (makeBTagEffPlots){
	for (int unsigned i = 0; i < bTagEffPlots.size(); i++){
	  bTagEffPlots[i]->Write();
	}
      }
      outFile.Write();
      outFile.Close();
      delete cloneTree;

      //If we have any events in the second tree:
      if (cloneTree2->GetEntries() > 0){
	std::cout << "There are " << cloneTree2->GetEntries() << " entries in the second tree!" << std::endl;
	TFile outFile1(("skims/"+dataset->name() + postfix + (invertIsoCut?"invIso":"") + "SmallSkim1.root").c_str(),"RECREATE");
	outFile1.cd();
	cloneTree2->Write();
	outFile1.Write();
	outFile1.Close();
      }
      if (cloneTree3->GetEntries() > 0){
	std::cout << "There are " << cloneTree3->GetEntries() << " entries in the third tree! What a lot of trees we've made." << std::endl;
	TFile outFile1(("skims/"+dataset->name() + postfix + (invertIsoCut?"invIso":"") + "SmallSkim2.root").c_str(),"RECREATE");
	outFile1.cd();
	cloneTree3->Write();
	outFile1.Write();
	outFile1.Close();
      }
      delete cloneTree2;
      delete cloneTree3;
    }


    
    //Save mva outputs
    if (makeMVATree){
      std::cout << (mvaDir + dataset->name() + postfix + (invertIsoCut?"invIso":"")  +  "mvaOut.root") << std::endl;
      TFile mvaOutFile((mvaDir + dataset->name() + postfix + (invertIsoCut?"invIso":"")  +  "mvaOut.root").c_str(),"RECREATE");
      mvaOutFile.cd();
      std::cout << std::endl;
      int systMask = 1;
      std::cout << "Saving Systematics: ";
      for (unsigned int systInd = 0; systInd < systNames.size(); systInd++){
	if (systInd > 0 && !(systToRun & systMask)){
	  systMask = systMask << 1;
	  continue;
	}
	std::cout << systNames[systInd] << ": " << mvaTree[systInd]->GetEntriesFast() << " " << std::flush;
	mvaTree[systInd]->Write();
	if (systInd > 0) systMask = systMask << 1;
	if (!dataset->isMC()) break;
      }
      std::cout << std::endl;
      //Save the efficiency plots for b-tagging here if we're doing that.
      if (makeBTagEffPlots){
	for (int unsigned i = 0; i < bTagEffPlots.size(); i++){
	  bTagEffPlots[i]->Write();
	}
      }
      mvaOutFile.Write();
      mvaOutFile.Close();
      for (unsigned int i = 0; i < mvaTree.size(); i++){
	delete mvaTree[i];
      }
    }
    if (infoDump){
      std::cout << "In dataset " << dataset->getFillHisto() << " the cut flow looks like:" << std::endl;
      for (int i = 0; i < cutFlowMap[dataset->getFillHisto()]->GetNbinsX(); i++){
	std::cout << stageNames[i] << "\t" << cutFlowMap[dataset->getFillHisto()]->GetBinContent(i+1) << std::endl;
      }
    }
    std::cerr << "\nFound " << foundEvents << " in " << dataset->name() << std::endl;
    //Delete plots from out btag vector. Avoid memory leaks, kids.
    if (makeBTagEffPlots){
      for (int unsigned i = 0; i < bTagEffPlots.size(); i++){
	delete bTagEffPlots[i];
      }
    }
    //datasetChain->MakeClass("AnalysisEvent");
  } // end channel loop.
  delete datasetChain;
 } //end dataset loop

//Save all plot objects. For testing purposes.
  
void analysisAlgo::savePlots();
{
  //Now test out the histogram plotter class I just wrote.
  //Make the plotting object.
  if (plots||infoDump){
    HistogramPlotter plotObj = HistogramPlotter(legOrder,plotOrder,datasetInfos);
    plotObj.setLabelOne("CMS Preliminary");
    plotObj.setLabelTwo("Some amount of lumi");
    plotObj.setPostfix("");
    plotObj.setOutputFolder(outFolder);

    for (unsigned int i = 0; i < plotsVec.size(); i++){
      std::cout << plotsVec[i] << std::endl;
      if (plots)
	plotObj.plotHistos(plotsMap[plotsVec[i]]);
    }
    
    std::vector<std::string> cutFlowLabels (4);
    cutFlowLabels = {"lepSel","zMass","jetSel","bTag"};
    //std::cout << cutFlowMap << std::endl;
    plotObj.makePlot(cutFlowMap,"cutFlow",cutFlowLabels);
  }
  if (synchCutFlow){
    cutObj->getSynchCutFlow();
  }

  /*  for (std::vector<Dataset>::iterator dataset = datasets.begin(); dataset!=datasets.end(); ++dataset){
    for (unsigned int j = 0; j < stageNames.size(); j++){
    plotsMap[dataset->getFillHisto()][stageNames[j]]->saveAllPlots();
    }
    cutFlowMap[dataset->getFillHisto()]->SaveAs(("plots/"+dataset->name()+"_cutFlow.root").c_str());    
    cutFlowMap[dataset->getFillHisto()]->Draw();
    cutFlowMap[dataset->getFillHisto()]->SaveAs(("plots/"+dataset->name()+"_cutFlow.png").c_str());
    }*/

  //Delete all the plot objects.

  std::cerr << "Gets to the delete bit" << std::endl;
  if (plots || infoDump){
    for (std::vector<Dataset>::iterator dataset = datasets.begin(); dataset!=datasets.end(); ++dataset){
      if (cutFlowMap.find(dataset->getFillHisto()) == cutFlowMap.end()) continue;
      delete cutFlowMap[dataset->getFillHisto()];
      if (!plots) continue;
      for (unsigned int j = 0; j < stageNames.size(); j++){
	int systMask = 1;
	for (unsigned int systInd = 0; systInd < systNames.size(); systInd++){
	  if (systInd > 0 && !(systInd & systMask)) {
	    systMask = systMask << 1;
	    continue;
	  }
	  delete plotsMap[systNames[systInd]][dataset->getFillHisto()][stageNames[j]];
	  if (systInd > 0) systMask = systMask << 1;
	}
      }
    }

  }
  delete cutConfName;
  delete plotConfName;
  std::cerr  << "But not past it" << std::endl;
}

