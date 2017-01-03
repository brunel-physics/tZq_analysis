#include "cutClass.hpp"

#include "TLorentzVector.h"
#include "TRandom.h"

#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>

#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraphAsymmErrors.h"

#include <libconfig.h++>

Cuts::Cuts( bool doPlots, bool fillCutFlows,bool invertLepCut, bool lepCutFlow, bool dumpEventNumber, const bool trileptonChannel, const bool is2016, const bool isFCNC, const bool isCtag ):

  //Do plots?
  doPlots_{doPlots},
  fillCutFlow_{fillCutFlows},
  //background estimation. May not be possible
  invertLepCut_{invertLepCut},
  //Synchronisation cut flow.
  synchCutFlow_{lepCutFlow},
  //Synchronisation cut flow.
  singleEventInfoDump_{false},
  makeEventDump_{dumpEventNumber},
  trileptonChannel_{trileptonChannel},
  isFCNC_{isFCNC},
  isCtag_{isCtag},
  is2016_{is2016},

  // Set all default parameters. These will be editable later on, probably.
  numTightEle_{3},
  tightElePt_{20.},
  tightEleEta_{2.5},
  tightEled0_{0.011811},
  tightEleMissLayers_{3}, //Cut based has barrel =2; endcap: veto=3, others=1
  tightEleCheckPhotonVeto_{true},
  tightEleMVA0_{0.972153}, // Medium cut
  tightEleMVA1_{0.922126}, // Medium cut
  tightEleMVA2_{0.610764}, // Medium cut
  //tightEleMVA0_{0.988153}, // Tight cut
  //tightEleMVA1_{0.967910}, // Tight cut
  //tightEleMVA2_{0.841729}, // Tight cut
  tightEleRelIso_{0.107587},
  //Loose electron initialisation
  numLooseEle_{3},
  looseElePt_{20},
  looseEleEta_{2.5},
  looseEleMVA0_{0.972153},
  looseEleMVA1_{0.922126},
  looseEleMVA2_{0.610764},
  looseEleRelIso_{0.15},
  //Tight muon initialisation
  numTightMu_{0},
  tightMuonPt_{20.},
  tightMuonEta_{2.4},
  tightMuonRelIso_{0.15},
  //Loose muons
  numLooseMu_{0},
  looseMuonPt_{20.},
  looseMuonEta_{2.4},
  looseMuonRelIso_{0.25},
  //zMass cuts
  invZMassCut_{10.},
  invWMassCut_{20.},
  //Jet initialisation
  numJets_{2},
  maxJets_{4},
  jetPt_{30.},
  jetEta_{5.0},
  jetNConsts_{2},
  jetIDDo_{true},
  //B-discriminator cut
  numbJets_{1},
  maxbJets_{2},
  looseBjetVeto_{0},
  bDiscCut_{0.9535}, // Tight cut
  //bDiscCut_{0.8484}, // Medium level
  //bDiscCut_{0.5426}, // Loose cut
  bLooseDiscCut_{0.5426}, // Loose cut
  bDiscSynchCut_{0.5426},
  //C-discriminator cut
  numcJets_{1},
  maxcJets_{1},
  //cVsLDiscCut_{0.69}, // Tight cut
  cVsLDiscCut_{-0.1}, // Medium level
  //cVsLDiscCut_{-0.48}, // Loose cut
  cVsBDiscCut_{0.45}, // Tight cut
  //cVsBDiscCut_{0.08}, // Medium level
  //cVsBDiscCut_{-0.17}, // Loose cut

  tempSmearValue_{1.0}, // Temporary solution to smearing propagation bug fix. A more elegant solution is needed!

  //Set isMC. Default is true, but it's called everytime a new dataset is processed anyway.
  isMC_{true},
  //Same for trigger flag.
  triggerFlag_{},
  //Make cloned tree false for now
  postLepSelTree_{nullptr},
  postLepSelTree2_{nullptr},
  postLepSelTree3_{nullptr},
  //Skips running trigger stuff
  skipTrigger_{false},
  //Are we making b-tag efficiency plots?
  makeBTagEffPlots_{false},
  getBTagWeight_{false},

  //MET and mTW cuts go here.
  metCut_{0.0},
  metDileptonCut_{100.0},
  mTWCut_{20.0},
  mTWCutSynch_{20.0},
  TopMassCutLower_{95.},
  TopMassCutUpper_{200.}

{
  //Space here in case other stuff needs to be done.

  //If doing synchronisation., initialise that here.
  if (synchCutFlow_){
    synchCutFlowHist_ = new TH1F{"synchCutFlow","synchCutFlow",11,0,11};
    synchNumEles_ = new TH1I{"synchNumEles","synchNumEles",11,0,11};
    synchNumMus_ = new TH1I{"synchNumMuos","synchNumMuos",11,0,11};
    synchMuonCutFlow_ = new TH1I{"synchMuonCutFlow","synchMuonCutFlow",11,0,11};
    synchCutTopMassHist_ = new TH1F{"synchCutTopMassHist", "synchCutTopMassHist", 200, 0., 200.};
  }

  std::cout << "\nInitialises fine" << std::endl;
  initialiseJECCors();
  std::cout << "Gets past JEC Cors" << std::endl;

  if ( !is2016_ ) {
    std::cout << "\nLoad 2015 electron SFs from root file ... " << std::endl;
    electronSFsFile = new TFile("scaleFactors/2015/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root"); // Electron cut-based Tight ID
    h_eleSFs = dynamic_cast<TH2F*>(electronSFsFile->Get("EGamma_SF2D"));
    electronRecoFile = new TFile{"scaleFactors/2015/eleRECO.txt.egamma_SF2D.root"}; // Electron Reco SF
    h_eleReco = dynamic_cast<TH2F*>(electronRecoFile->Get("EGamma_SF2D"));
    std::cout << "Got 2015 electron SFs!\n" << std::endl;

    std::cout << "Load 2015 muon SFs from root file ... " << std::endl;
    muonIDsFile = new TFile{"scaleFactors/2015/MuonID_Z_RunCD_Reco76X_Feb15.root"};
    muonIsoFile = new TFile{"scaleFactors/2015/MuonIso_Z_RunCD_Reco76X_Feb15.root"};
    muonIDsFile->cd("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1"); // Tight ID
    h_muonIDs = dynamic_cast<TH2F*>(muonIDsFile->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")); // Tight ID
    muonIsoFile->cd("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"); // Tight ID
    h_muonPFiso = dynamic_cast<TH2F*>(muonIsoFile->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio")); // Tight ID
    std::cout << "Got 2015 muon SFs!\n" << std::endl;
  }
  else {
    std::cout << "\nLoad 2016 electron SFs from root file ... " << std::endl;
    electronSFsFile = new TFile("scaleFactors/2016/egammaEffi.txt_SF2D.root"); // Electron cut-based Tight ID
    h_eleSFs = dynamic_cast<TH2F*>(electronSFsFile->Get("EGamma_SF2D"));
    electronRecoFile = new TFile{"scaleFactors/2016/egammaRecoEffi.txt_SF2D.root"}; // Electron Reco SF
    h_eleReco = dynamic_cast<TH2F*>(electronRecoFile->Get("EGamma_SF2D"));
    std::cout << "Got 2016 electron SFs!\n" << std::endl;

    std::cout << "Load 2016 muon SFs from root file ... " << std::endl;
    muonIDsFile = new TFile{"scaleFactors/2016/MuonID_Z_RunBCD_prompt80X_7p65.root"};
    muonIsoFile = new TFile{"scaleFactors/2016/MuonIso_Z_RunBCD_prompt80X_7p65.root"};
    muonRecoFile = new TFile{"scaleFactors/2016/MuonReco_RunBCD_prompt80X_9p2.root"};

    muonIDsFile->cd("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1"); // Tight ID
    h_muonIDs = dynamic_cast<TH2F*>(muonIDsFile->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")); // Medium ID for ICHEP
    muonIsoFile->cd("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"); // Tight ID
    h_muonPFiso = dynamic_cast<TH2F*>(muonIsoFile->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio")); // Medium ID for ICHEP
    h_muonRecoGraph = dynamic_cast<TGraphAsymmErrors*>(muonRecoFile->Get("ratio_eta")); // HIP issue means muon tracking RECO efficiency != ~1
    std::cout << "Got 2016 muon SFs!\n" << std::endl;
  }

  // Setup bTag calibration code (2015/2016)
  // bTag calib code
  calib2015 = BTagCalibration ("CSVv2", "scaleFactors/2015/CSVv2.csv");
  calib2016 = BTagCalibration ("CSVv2", "scaleFactors/2016/CSVv2.csv");

  // udsg jets
  lightReader = BTagCalibrationReader (BTagEntry::OP_TIGHT,"central",{"up","down"});// operating point
  // c/b jets
  charmReader = BTagCalibrationReader (BTagEntry::OP_TIGHT, "central",{"up","down"});	//central
  beautyReader = BTagCalibrationReader (BTagEntry::OP_TIGHT, "central",{"up","down"});	//central


  // if doing bTag SFs, load stuff here ...
  // N.B. 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG
  if ( getBTagWeight_ ) {
    if ( !is2016_ ) {
      lightReader.load(calib2015, BTagEntry::FLAV_UDSG, "incl");
      charmReader.load(calib2015, BTagEntry::FLAV_C, "mujets");
      beautyReader.load(calib2015, BTagEntry::FLAV_B, "mujets");
    }
    else {
      lightReader.load(calib2016, BTagEntry::FLAV_UDSG, "incl");
      charmReader.load(calib2016, BTagEntry::FLAV_C, "mujets");
      beautyReader.load(calib2016, BTagEntry::FLAV_B, "mujets");
    }
  }
}

Cuts::~Cuts(){

  electronSFsFile->Close();
  electronRecoFile->Close();
  muonIDsFile->Close();
  muonIsoFile->Close();

  if (synchCutFlow_) {
    delete synchCutFlowHist_;
    delete synchNumEles_;
    delete synchNumMus_;
    delete synchMuonCutFlow_;
    delete synchCutTopMassHist_;
    if (makeEventDump_) {
      topMassEventDump_.close();
      step0EventDump_.close();
      step2EventDump_.close();
      step4EventDump_.close();
      step6EventDump_.close();
      step9EventDump_.close();
    }
  }
}

bool Cuts::parse_config(std::string confName){
  libconfig::Config cutConf;

  //Get the configuration file
  try {
    cutConf.readFile(confName.c_str());
  }
  catch (const libconfig::FileIOException &exep){
    std::cerr << "Error opening cut configuration file" << std::endl;
    return 0;
  }
  catch(const libconfig::ParseException &e){
    std::cerr << "Parse error in cut configuration at " << e.getFile() << ":" << e.getLine() << " - " << e.getError() << std::endl;
    return 0;
  }
  const libconfig::Setting& root = cutConf.getRoot();

  if (! root.exists("cuts")){
    std::cerr << "I don't know what you think you're doing, but there aren't any cuts here. Get back to the drawing board." << std::endl;
    return false;
  }
  if (root.exists("trigLabel")) root.lookupValue("trigLabel",cutConfTrigLabel_);
  if (root.exists("plotPostfix")) root.lookupValue("plotPostfix",postfixName_);
  libconfig::Setting& cuts = root["cuts"];
  if (cuts.exists("tightElectrons")){
    libconfig::Setting& eles = cuts["tightElectrons"];
    eles.lookupValue("pt",tightElePt_);
    eles.lookupValue("eta",tightEleEta_);
    eles.lookupValue("relIso",tightEleRelIso_);
    eles.lookupValue("number",numTightEle_);
  }
  if (cuts.exists("looseElectrons")){
    libconfig::Setting& eles = cuts["looseElectrons"];
    eles.lookupValue("pt",looseElePt_);
    eles.lookupValue("eta",looseEleEta_);
    eles.lookupValue("relIso",looseEleRelIso_);
    eles.lookupValue("number",numLooseEle_);
  }
  if (cuts.exists("tightMuons")){
    libconfig::Setting& mus = cuts["tightMuons"];
    mus.lookupValue("pt",tightMuonPt_);
    mus.lookupValue("eta",tightMuonEta_);
    mus.lookupValue("relIso",tightMuonRelIso_);
    mus.lookupValue("number",numTightMu_);
  }
  if (cuts.exists("looseMuons")){
    libconfig::Setting& mus = cuts["looseMuons"];
    mus.lookupValue("pt",looseMuonPt_);
    mus.lookupValue("eta",looseMuonEta_);
    mus.lookupValue("relIso",looseMuonRelIso_);
    mus.lookupValue("number",numLooseMu_);
  }
  if (cuts.exists("jets")){
    libconfig::Setting& jets = cuts["jets"];
    jets.lookupValue("pt",jetPt_);
    jets.lookupValue("eta",jetEta_);
    jets.lookupValue("numJets",numJets_);
    jets.lookupValue("maxJets",maxJets_);
    jets.lookupValue("numbJets",numbJets_);
    jets.lookupValue("maxbJets",maxbJets_);
    jets.lookupValue("numcJets",numcJets_);
    jets.lookupValue("maxcJets",maxcJets_);
    jets.lookupValue("looseBjetVeto",looseBjetVeto_);
  }

  std::cerr << "And so it's looking for " << numTightMu_ << " muons and " << numTightEle_ << " electrons" << std::endl;

  if (makeEventDump_) {
    topMassEventDump_.open("topMassEventDump"+postfixName_+".txt");
    step0EventDump_.open("step0EventDump"+postfixName_+".txt");
    step2EventDump_.open("step2EventDump"+postfixName_+".txt");
    step4EventDump_.open("step4EventDump"+postfixName_+".txt");
    step6EventDump_.open("step6EventDump"+postfixName_+".txt");
    step9EventDump_.open("step9EventDump"+postfixName_+".txt");
  }

  return true;
}

bool Cuts::makeCuts(AnalysisEvent *event, float *eventWeight, std::map<std::string,Plots*> plotMap, TH1F* cutFlow, int systToRun){

  //If we're doing synchronisation, do this function.
  if (synchCutFlow_){
    return synchCuts(event, eventWeight);
  }

  //If emu and dilepton - doing ttbar background estimation

  if ( !trileptonChannel_ && cutConfTrigLabel_.find("d") != std::string::npos ){
     return ttbarCuts(event, eventWeight, plotMap, cutFlow, systToRun);
  }

  if( !skipTrigger_ ) {
    if ( !is2016_ ) if (!triggerCuts(event)) return false; // Do trigger on MC and data for 2015
    if ( !isMC_ && is2016_ ) if (!triggerCuts(event)) return false; // Do trigger for data and 2016, exclude MC.

    if( !is2016_ && isMC_ ) *eventWeight *= get2015TriggerSF (systToRun); // Apply SFs to MC if 2015
    else if ( is2016_ && isMC_ ) *eventWeight *= get2016TriggerSF (systToRun); // Apply data efficiencies onto MC if 2016 due to incomplete trigger menu
  }

  if (!metFilters(event)) return false;

  //Make lepton cuts. Does the inverted iso cuts if necessary.
  if (invertLepCut_ && trileptonChannel_) {
     if ( !invertIsoCut(event,eventWeight, plotMap,cutFlow) ) return false;
  }
  else {
    if ( !makeLeptonCuts(event,eventWeight,plotMap,cutFlow,systToRun) ) return false;
  }

  event->jetIndex = makeJetCuts(event, systToRun, eventWeight);
  if (event->jetIndex.size() < numJets_) return false;
  if (event->jetIndex.size() > maxJets_) return false;
  if (doPlots_||fillCutFlow_) cutFlow->Fill(2.5,*eventWeight);
  if (doPlots_) plotMap["jetSel"]->fillAllPlots(event,*eventWeight);

  event->bTagIndex = makeBCuts(event,event->jetIndex, systToRun);
  if ( looseBjetVeto_ == 1 ) event->bTagLooseIndex = makeLooseBCuts(event,event->jetIndex, systToRun);
  if (event->bTagIndex.size() < numbJets_) return false;
  if (event->bTagIndex.size() > maxbJets_) return false;
  if ( looseBjetVeto_ == 1 && ( event->bTagIndex.size() != event->bTagLooseIndex.size() ) ) return false;
  if (doPlots_) plotMap["bTag"]->fillAllPlots(event,*eventWeight);
  if (doPlots_||fillCutFlow_) cutFlow->Fill(3.5,*eventWeight);

  if ( !trileptonChannel_ && !isFCNC_ ) { // Do wMass stuff
    float invWmass{0.};
    invWmass = getWbosonQuarksCand(event,event->jetIndex,systToRun);

   // Debug chi2 cut
//   float topMass = getTopMass(event);
//   float topTerm = ( topMass-173.21 )/30.0;
//   float wTerm = ( (event->wPairQuarks.first + event->wPairQuarks.second).M() - 80.3585 )/8.0;

//   float chi2 = topTerm*topTerm + wTerm*wTerm;
//   if ( chi2 < 2.0 && chi2 > 7.0 ) return false; // control region
//   if ( chi2 >= 2.0 ) return false; //signal region

    if ( std::abs(invWmass) > invWMassCut_ ) return false;
    if ( doPlots_ ) plotMap["wMass"]->fillAllPlots(event,*eventWeight);
    if ( doPlots_ || fillCutFlow_ ) cutFlow->Fill(4.5,*eventWeight);
  }

  if (!trileptonChannel_ && isFCNC_) // Do FCNC stuff
  {
      // Leading jet cannot be b-tagged
      if (event->jetPF2PATBDiscriminator[0] > bDiscCut_)
      {
          return false;
      }

      if (isCtag_) // Do cTagging
      {
          event->cTagIndex = makeCCuts(event,event->jetIndex);

          if (event->cTagIndex.size() < numcJets_
              || event->cTagIndex.size() > maxJets_)
          {
              return false;
          }
          if (doPlots_ || fillCutFlow_)
          {

              if (doPlots_)
              {
                  plotMap["cTag"]->fillAllPlots(event,*eventWeight);
              }

              cutFlow->Fill(4.5,*eventWeight);
          }
      }
  }

  //Apply met and mtw cuts here. By default these are 0, so don't do anything.
  if (trileptonChannel_ && !isFCNC_ && event->metPF2PATPt < metCut_) return false;
//  if (!trileptonChannel_ && !isFCNC_ && event->metPF2PATPt > metDileptonCut_) return false;
  TLorentzVector tempMet;
  tempMet.SetPtEtaPhiE(event->metPF2PATPt,0,event->metPF2PATPhi,event->metPF2PATEt);
  double mtw{std::sqrt(2*event->metPF2PATPt*event->wLepton.Pt()*(1-std::cos(event->metPF2PATPhi - event->wLepton.Phi())))};
  if (trileptonChannel_ && !isFCNC_ && mtw < mTWCut_) return false;
  return true;
}

//Make lepton cuts. Will become customisable in a config later on.
bool Cuts::makeLeptonCuts(AnalysisEvent* event,float * eventWeight,std::map<std::string,Plots*> plotMap, TH1F* cutFlow, int syst, bool isControl){

  //Do lepton selection.
  event->electronIndexTight = getTightEles(event);
  if (event->electronIndexTight.size() != numTightEle_) return false;
  event->electronIndexLoose = getLooseEles(event);
  if (event->electronIndexLoose.size() != numLooseEle_) return false;

  event->muonIndexTight = getTightMuons(event);
  if (event->muonIndexTight.size() != numTightMu_) return false;
  event->muonIndexLoose = getLooseMuons(event);
  if (event->muonIndexLoose.size() != numLooseMu_) return false;

  //This is to make some skims for faster running. Do lepSel and save some files.
  if(postLepSelTree_ && !synchCutFlow_) {
    if (postLepSelTree_->GetEntriesFast() < 40000) postLepSelTree_->Fill();
    else {
      if (postLepSelTree2_->GetEntriesFast() < 40000) postLepSelTree2_->Fill();
      else postLepSelTree3_->Fill();

    }
  }

  //Should I make it return which leptons are the zMass candidate? Probably.
  float invZmass{9999.};

  if (trileptonChannel_ == true){
    invZmass = getZCand(event, event->electronIndexTight, event->muonIndexTight);
  }
  else if (trileptonChannel_ == false && !isControl){ //Control doesn't need Z plots
    invZmass = getDileptonZCand(event, event->electronIndexTight, event->muonIndexTight);
  }
  else if ( isControl ){
    invZmass = getTTbarCand(event, event->electronIndexTight, event->muonIndexTight);
  }
  else if (!isControl){
    std::cout << "You've somehow managed to set the trilepton/dilepton bool flag to a value other than these two!" << std::endl;
    std::cout << "HOW?! Well done for breaking this ..." << std::endl;
    exit(0);
  }

    * eventWeight *= getLeptonWeight(event,syst);

  if(doPlots_) plotMap["lepSel"]->fillAllPlots(event,*eventWeight);
  if(doPlots_||fillCutFlow_) cutFlow->Fill(0.5,*eventWeight);

  if (std::abs(invZmass) > invZMassCut_ && !isControl) return false;
  if (std::abs(invZmass) < 106 && isControl) return false;

  if(doPlots_) plotMap["zMass"]->fillAllPlots(event,*eventWeight);
  if (doPlots_||fillCutFlow_) cutFlow->Fill(1.5,*eventWeight);

  return true;
}

std::vector<int> Cuts::getTightEles(AnalysisEvent* event) {
  std::vector<int> electrons;
  for (int i{0}; i < event->numElePF2PAT; i++){

    if (!event->elePF2PATIsGsf[i]) continue;
    TLorentzVector tempVec{event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]};

    if (tempVec.Pt() <= tightElePt_) continue;
    if (std::abs(tempVec.Eta()) >= tightEleEta_) continue;

    // 2015 cuts
    if ( !is2016_ ) {
      if (postLepSelTree_) { //If post-lep tree creation, keep more info (medium cuts)
        if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){
	  if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0101 ) continue;
	  if ( std::abs(event->elePF2PATDeltaEtaSC[i]) >= 0.0103 ) continue;
	  if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.0336 ) continue;
	  if ( event->elePF2PATHoverE[i] >= 0.0876 ) continue;
	  if ( event->elePF2PATComRelIsoRho[i] >= 0.0766 ) continue;
	  if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.0174 ) continue;
	  if ( std::abs(event->elePF2PATD0PV[i]) >= 0.0118 )continue;
	  if ( std::abs(event->elePF2PATDZPV[i]) >= 0.0373 ) continue;
	  if ( event->elePF2PATMissingInnerLayers[i] > 2 ) continue;
          if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ ) continue;
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
          if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ ) continue;
	}
	else continue;
      }

      else { // Else do tight cut-based ID
        if ( std::abs(event->elePF2PATSCEta[i]) > 1.4442 && std::abs(event->elePF2PATSCEta[i]) < 1.566 ) continue;
   	  // Barrel cut-based ID
          if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){
          if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0101 ) continue;
	  if ( std::abs(event->elePF2PATDeltaEtaSC[i]) >= 0.00926 ) continue;
	  if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.0336 ) continue;
	  if ( event->elePF2PATHoverE[i] >= 0.0597 ) continue;
	  if ( event->elePF2PATComRelIsoRho[i] >= 0.0354 ) continue;
	  if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.012 ) continue;
	  if ( std::abs(event->elePF2PATD0PV[i]) >= 0.0111 )continue;
	  if ( std::abs(event->elePF2PATDZPV[i]) >= 0.0466 ) continue;
	  if ( event->elePF2PATMissingInnerLayers[i] > 2 ) continue;
          if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ ) continue;

	}
        else if ( std::abs(event->elePF2PATSCEta[i]) > 1.479 && std::abs(event->elePF2PATSCEta[i]) < 2.50 ){ // Endcap cut-based ID
	  if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0279 ) continue;
	  if ( std::abs(event->elePF2PATDeltaEtaSC[i]) >= 0.00724 ) continue;
	  if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.0918 ) continue;
	  if ( event->elePF2PATHoverE[i] >= 0.0615 ) continue;
	  if ( event->elePF2PATComRelIsoRho[i] >= 0.0646 ) continue;
	  if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.00999 ) continue;
	  if ( std::abs(event->elePF2PATD0PV[i]) >= 0.0351 )continue;
	  if ( std::abs(event->elePF2PATDZPV[i]) >= 0.417 ) continue;
	  if ( event->elePF2PATMissingInnerLayers[i] > 1 ) continue;
          if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ ) continue;
	}
	else continue;
      }
    }

    // 2016 cuts
    else {

      // Ensure we aren't in the barrel/endcap gap and below the max safe eta range
      if ( (std::abs(event->elePF2PATSCEta[i]) > 1.4442 && std::abs(event->elePF2PATSCEta[i]) < 1.566) || std::abs(event->elePF2PATSCEta[i]) > 2.50 ) continue;

      // VID cut
      if ( event->elePF2PATCutIdTight[i] < 1 ) continue;

      // Cuts not part of the tuned ID
      if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){
        if ( std::abs(event->elePF2PATD0PV[i]) >= 0.05  ) continue;
        if ( std::abs(event->elePF2PATDZPV[i]) >= 0.10  ) continue;
      }
      else if ( std::abs(event->elePF2PATSCEta[i]) > 1.479 && std::abs(event->elePF2PATSCEta[i]) < 2.50 ){
        if ( std::abs(event->elePF2PATD0PV[i]) >= 0.10 ) continue;
        if ( std::abs(event->elePF2PATDZPV[i]) >= 0.20 ) continue;
      }

      // Old hand cuts for cut based IDs
/*
      // Barrel cut-based ID
      if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){
        if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.00998 ) continue;
        if ( std::abs(event->elePF2PATDeltaEtaSeedSC[i]) >= 0.00308 ) continue;
        if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.0816 ) continue;
        if ( event->elePF2PATHoverE[i] >= 0.0414 ) continue;
        if ( event->elePF2PATComRelIsoRho[i] >= 0.0588 ) continue;
        if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.0129 ) continue;
        if ( event->elePF2PATMissingInnerLayers[i] > 1 ) continue;
        if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ ) continue;
      }
      else if ( std::abs(event->elePF2PATSCEta[i]) > 1.479 && std::abs(event->elePF2PATSCEta[i]) < 2.50 ){ // Endcap cut-based ID
        if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0292 ) continue;
        if ( std::abs(event->elePF2PATDeltaEtaSeedSC[i]) >= 0.00605 ) continue;
        if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.0394 ) continue;
        if ( event->elePF2PATHoverE[i] >= 0.0641 ) continue;
        if ( event->elePF2PATComRelIsoRho[i] >= 0.0571 ) continue;
        if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.0129 ) continue;
        if ( event->elePF2PATMissingInnerLayers[i] > 1 ) continue;
        if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ ) continue;
        }
	else continue; // If outside barrel range, kill
*/
      }
    electrons.emplace_back(i);
  }
  return electrons;
}

std::vector<int> Cuts::getLooseEles(AnalysisEvent* event){
  std::vector<int> electrons;
  for (int i{0}; i < event->numElePF2PAT; i++){
    TLorentzVector tempVec{event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]};

    if (tempVec.Pt() < tightElePt_) continue;
    if (std::abs(tempVec.Eta()) > tightEleEta_)continue;

    // 2015 cuts
    if ( !is2016_ ) {
    // Barrel cut-based Veto ID
      if ( std::abs(event->elePF2PATSCEta[i]) > 1.4442 && std::abs(event->elePF2PATSCEta[i]) < 1.566 && postLepSelTree_ ) continue;
      if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){
	if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0114 ) continue;
	if ( std::abs(event->elePF2PATDeltaEtaSC[i]) >= 0.0152 ) continue;
	if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.216 ) continue;
	if ( event->elePF2PATHoverE[i] >= 0.181 ) continue;
	if ( synchCutFlow_ && event->elePF2PATComRelIsoRho[i] >= 0.0354 ) continue; // Use same rel iso as tight
	if ( !synchCutFlow_ && event->elePF2PATComRelIsoRho[i] >= 0.126 ) continue; // Use same rel iso as loose
	if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.207 ) continue;
	if ( std::abs(event->elePF2PATD0PV[i]) >= 0.0564)continue;
	if ( std::abs(event->elePF2PATDZPV[i]) >= 0.472 ) continue;
	if ( event->elePF2PATMissingInnerLayers[i] > 2  ) continue;
        if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ )continue;
      }
      else if ( std::abs(event->elePF2PATSCEta[i]) > 1.479 && std::abs(event->elePF2PATSCEta[i]) < 2.50 ){
	if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0352 ) continue;
	if ( std::abs(event->elePF2PATDeltaEtaSC[i]) >= 0.0113 ) continue;
	if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.237 ) continue;
	if ( event->elePF2PATHoverE[i] >= 0.116 ) continue;
	if ( synchCutFlow_ && event->elePF2PATComRelIsoRho[i] >= 0.0646 ) continue; // Use same rel iso as tight
	if ( !synchCutFlow_ && event->elePF2PATComRelIsoRho[i] >= 0.144 ) continue; // Use same rel iso as loose
	if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.174 ) continue;
	if ( std::abs(event->elePF2PATD0PV[i]) >= 0.222 )continue;
	if ( std::abs(event->elePF2PATDZPV[i]) >= 0.921 ) continue;
	if ( event->elePF2PATMissingInnerLayers[i] > 3 ) continue;
        if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ )continue;
      }
      else continue;
    }

    // 2016 cuts
    else {
      // Barrel cut-based Veto ID

      // Ensure we aren't in the barrel/endcap gap and below the max safe eta range
      if ( (std::abs(event->elePF2PATSCEta[i]) > 1.4442 && std::abs(event->elePF2PATSCEta[i]) < 1.566) || std::abs(event->elePF2PATSCEta[i]) > 2.50 ) continue;

      // VID cut
      if ( event->elePF2PATCutIdVeto[i] < 1 ) continue;

      // Cuts not part of the tuned ID
      if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){
        if ( std::abs(event->elePF2PATD0PV[i]) >= 0.05  ) continue;
        if ( std::abs(event->elePF2PATDZPV[i]) >= 0.10  ) continue;
      }
      else if ( std::abs(event->elePF2PATSCEta[i]) > 1.479 && std::abs(event->elePF2PATSCEta[i]) < 2.50 ){
        if ( std::abs(event->elePF2PATD0PV[i]) >= 0.10 ) continue;
        if ( std::abs(event->elePF2PATDZPV[i]) >= 0.20 ) continue;
      }

      // Old hand cuts for cut based IDs
/*
      if ( std::abs(event->elePF2PATSCEta[i]) <= 1.479 ){

        if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.0115 ) continue;
        if ( std::abs(event->elePF2PATDeltaEtaSeedSC[i]) >= 0.00749 ) continue;
        if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.228 ) continue;
        if ( event->elePF2PATHoverE[i] >= 0.356 ) continue;
        if ( event->elePF2PATComRelIsoRho[i] >= 0.175 ) continue;
        if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.299 ) continue;
        if ( event->elePF2PATMissingInnerLayers[i] > 2  ) continue;
        if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ )continue;
      }
      else if ( std::abs(event->elePF2PATSCEta[i]) > 1.479 && std::abs(event->elePF2PATSCEta[i]) < 2.50 ){

        if ( event->elePF2PATSCSigmaIEtaIEta5x5[i] >= 0.037 ) continue;
        if ( std::abs(event->elePF2PATDeltaEtaSeedSC[i]) >= 0.00895 ) continue;
        if ( std::abs(event->elePF2PATDeltaPhiSC[i]) >= 0.213 ) continue;
        if ( event->elePF2PATHoverE[i] >= 0.211 ) continue;
        if ( event->elePF2PATComRelIsoRho[i] >= 0.159 ) continue;
        if ( (std::abs(1.0 - event->elePF2PATSCEoverP[i])*(1.0/event->elePF2PATEcalEnergy[i])) >= 0.15 ) continue;
        if ( event->elePF2PATMissingInnerLayers[i] > 3 ) continue;
        if ( event->elePF2PATPhotonConversionTag[i] && tightEleCheckPhotonVeto_ )continue;
      }
      else continue;
*/
    }
    electrons.emplace_back(i);

  }
  return electrons;
}

std::vector<int> Cuts::getTightMuons(AnalysisEvent* event){
  std::vector<int> muons;
  for (int i{0}; i < event->numMuonPF2PAT; i++){
    //if (i == 0) std::cout << "Starts doing first..." << event->muonPF2PATIsPFMuon[i];
    //if (i > 0) std::cout << "   " << event->muonPF2PATIsPFMuon[i];
    //    std::cout << i << "\t" << event->muonPF2PATPt[i] << "\t" << event->muonPF2PATEta[i] << "\t" << event->muonPF2PATComRelIsodBeta[i] << "\t" << event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] << "\t" << event->muonPF2PATGlobalID[i] << "\t" << event->muonPF2PATIsPFMuon[i] << "\t" << event->muonPF2PATTrackID[i] << "\t" <<  event->muonPF2PATTkLysWithMeasurements[i] << "\t" << event->muonPF2PATDBPV[i] << "\t" << event->muonPF2PATTrackNHits[i] << "\t" << event->muonPF2PATMuonNHits[i] << "\t" << event->muonPF2PATVldPixHits[i] << "\t" << event->muonPF2PATMatchedStations[i] << "\t" << std::abs(event->pvZ - event->muonPF2PATVertZ[i]) << std::endl;

    if (!event->muonPF2PATIsPFMuon[i]) continue;

    if (event->muonPF2PATPt[i] <= tightMuonPt_) continue;
    if (std::abs(event->muonPF2PATEta[i]) >= tightMuonEta_) continue;
    if (event->muonPF2PATComRelIsodBeta[i] >= tightMuonRelIso_) continue;

    // 2015 cuts
//    if ( !is2016_ ) {
      //Do a little test of muon id stuff here.
      if (!event->muonPF2PATTrackID[i]) continue;
      if (!event->muonPF2PATGlobalID[i]) continue;

      if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] >= 10.) continue;
      if (event->muonPF2PATTkLysWithMeasurements[i] <= 5) continue;
      if (std::abs(event->muonPF2PATDBPV[i]) >= 0.2) continue;
      if (std::abs(event->muonPF2PATDZPV[i]) >= 0.5) continue;
      if (event->muonPF2PATMuonNHits[i] < 1) continue;
      if (event->muonPF2PATVldPixHits[i] < 1) continue;
      if (event->muonPF2PATMatchedStations[i] < 2) continue;
//    }
    // 2016 cuts
    else {

      // If not either track muon and global muon ...
      if ( !(event->muonPF2PATTrackID[i]) && !(event->muonPF2PATGlobalID[i]) ) continue; // Normal loose ID on top of ICHEP cuts
      if ( event->muonPF2PATValidFraction[i] <= 0.49 ) continue;

      bool goodGlobalMuon (true), tightSegmentCompatible (true);

      if (!event->muonPF2PATTrackID[i]) goodGlobalMuon = false;
      if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] >= 3.) goodGlobalMuon = false;
      if (event->muonPF2PATChi2LocalPosition[i] >= 12.) goodGlobalMuon = false;
      if (event->muonPF2PATTrkKick[i] >= 20.) goodGlobalMuon = false;
      if (event->muonPF2PATSegmentCompatibility[i] <= 0.303) goodGlobalMuon = false;

      if (event->muonPF2PATSegmentCompatibility[i] <= 0.451) tightSegmentCompatible = false;

      // If both good global muon and tight segment compatible are not true ...
      if ( !(goodGlobalMuon) && !(tightSegmentCompatible) ) continue;
    }
    muons.emplace_back(i);
  }
  return muons;
}

std::vector<int> Cuts::getLooseMuons(AnalysisEvent* event){
  std::vector<int> muons;
  for (int i{0}; i < event->numMuonPF2PAT; i++){
    if (!event->muonPF2PATIsPFMuon[i]) continue;
    if (event->muonPF2PATPt[i] <= tightMuonPt_) continue;
    if (std::abs(event->muonPF2PATEta[i]) >= looseMuonEta_) continue;
    if (event->muonPF2PATComRelIsodBeta[i] >= looseMuonRelIso_) continue;
    if (event->muonPF2PATGlobalID[i] || event->muonPF2PATTrackID[i])
	muons.emplace_back(i);
    else continue;
  }
  return muons;
}

float Cuts::getZCand(AnalysisEvent *event, std::vector<int> electrons, std::vector<int> muons){
  float closestMass{9999.};
  //Use electrons if there are at least 2. Otherwise use muons.
  if (electrons.size() > 1){ // electrons.size() == number of electrons for selected channel
    for (unsigned i{0}; i < electrons.size(); i++){
      for (unsigned j{i + 1}; j < electrons.size(); j++) {
	if (event->elePF2PATCharge[electrons[i]] * event->elePF2PATCharge[electrons[j]] > 0) continue;
	TLorentzVector lepton1{event->elePF2PATGsfPx[electrons[i]],event->elePF2PATGsfPy[electrons[i]],event->elePF2PATGsfPz[electrons[i]],event->elePF2PATGsfE[electrons[i]]};
    TLorentzVector lepton2{event->elePF2PATGsfPx[electrons[j]],event->elePF2PATGsfPy[electrons[j]],event->elePF2PATGsfPz[electrons[j]],event->elePF2PATGsfE[electrons[j]]};
	double invMass{(lepton1 + lepton2).M() -91.1};
	if (std::abs(invMass) < std::abs(closestMass)){
	  // set up the tlorentz vectors in the event. For plotting and jazz.
	  event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
	  event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? electrons[i]:electrons[j];
	  event->zPairRelIso.first = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[electrons[i]]:event->elePF2PATComRelIsoRho[electrons[j]];
	  event->zPairRelIso.second = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[electrons[j]]:event->elePF2PATComRelIsoRho[electrons[i]];
	  event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
	  event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? electrons[j]:electrons[i];
	  closestMass = invMass;
	  //Now set up W lepton ...
	  if (electrons.size() == 2) {
	    event->wLepton = TLorentzVector(event->muonPF2PATPX[muons[0]],event->muonPF2PATPY[muons[0]],event->muonPF2PATPZ[muons[0]],event->muonPF2PATE[muons[0]]);
	    event->wLeptonRelIso = event->muonPF2PATComRelIsodBeta[muons[0]];
	    event->wLepIndex = muons[0];
	  }
	  else {
	    for (unsigned k{0}; k < electrons.size(); k++){
	      if (k == i || k == j) continue;
	      event->wLepton = TLorentzVector(event->elePF2PATGsfPx[electrons[k]],event->elePF2PATGsfPy[electrons[k]],event->elePF2PATGsfPz[electrons[k]],event->elePF2PATGsfE[electrons[k]]);
	      event->wLeptonRelIso = event->elePF2PATComRelIsoRho[electrons[k]];
	      event->wLepIndex = electrons[k];
	    }
	  }
	}
      }
    }
  }
  else {
    for (unsigned i{0}; i < muons.size(); i++){
      for (unsigned j{i + 1}; j < muons.size(); j++) {
	if (event->muonPF2PATCharge[muons[i]] * event->muonPF2PATCharge[muons[j]] > 0) continue;
	TLorentzVector lepton1{event->muonPF2PATPX[muons[i]],event->muonPF2PATPY[muons[i]],event->muonPF2PATPZ[muons[i]],event->muonPF2PATE[muons[i]]};
	TLorentzVector lepton2{event->muonPF2PATPX[muons[j]],event->muonPF2PATPY[muons[j]],event->muonPF2PATPZ[muons[j]],event->muonPF2PATE[muons[j]]};
	double invMass{(lepton1 + lepton2).M() -91};
	if (std::abs(invMass) < std::abs(closestMass)){
	  // set up the tlorentz vectors in the event. For plotting and jazz.
	  event->zPairLeptons.first = lepton1;
	  event->zPairIndex.first = muons[i];
	  event->zPairLeptons.second = lepton2;
	  event->zPairIndex.second = muons[j];
	  event->zPairRelIso.first = event->muonPF2PATComRelIsodBeta[muons[i]];
	  event->zPairRelIso.second = event->muonPF2PATComRelIsodBeta[muons[j]];
	  closestMass = invMass;
	  //Now set up W lepton
	  if (muons.size() == 2){
	    event->wLepton = TLorentzVector(event->elePF2PATGsfPx[electrons[0]],event->elePF2PATGsfPy[electrons[0]],event->elePF2PATGsfPz[electrons[0]],event->elePF2PATGsfE[electrons[0]]);
	    event->wLeptonRelIso = event->elePF2PATComRelIsoRho[electrons[0]];
	    event->wLepIndex = electrons[0];
	  }
	  else {
	    for (unsigned k{0}; k < muons.size(); k++){
	      if (k == i || k == j) continue;
	      event->wLepton = TLorentzVector(event->muonPF2PATPX[muons[k]],event->muonPF2PATPY[muons[k]],event->muonPF2PATPZ[muons[k]],event->muonPF2PATE[muons[k]]);
	      event->wLeptonRelIso = event->muonPF2PATComRelIsodBeta[muons[k]];
	      event->wLepIndex = muons[k];
	    }
	  }
	}
      }
    }
  }
  return closestMass;
}

float Cuts::getDileptonZCand(AnalysisEvent *event, std::vector<int> electrons, std::vector<int> muons){

  float closestMass {9999.};

  //Check if there are at least two electrons first. Otherwise use muons.

  if (electrons.size() > 1){
    for ( unsigned i{0}; i < electrons.size(); i++ ){
      for ( unsigned j{i + 1}; j < electrons.size(); j++ ){
        if ( !invertLepCut_) {
          if ( event->elePF2PATCharge[electrons[i]] * event->elePF2PATCharge[electrons[j]] >= 0 )  continue; // check electron pair have correct charge.
        }
        else {
          if ( !(event->elePF2PATCharge[electrons[i]] * event->elePF2PATCharge[electrons[j]] >= 0) ) continue; // check electron pair have correct charge for same sign control region.
        }
        TLorentzVector lepton1{event->elePF2PATGsfPx[electrons[i]],event->elePF2PATGsfPy[electrons[i]],event->elePF2PATGsfPz[electrons[i]],event->elePF2PATGsfE[electrons[i]]};
        TLorentzVector lepton2{event->elePF2PATGsfPx[electrons[j]],event->elePF2PATGsfPy[electrons[j]],event->elePF2PATGsfPz[electrons[j]],event->elePF2PATGsfE[electrons[j]]};
        double invMass{(lepton1 + lepton2).M() -91.1};
	if (std::abs(invMass) < std::abs(closestMass)){
	  event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
	  event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? electrons[i]:electrons[j];
	  event->zPairRelIso.first = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[electrons[i]]:event->elePF2PATComRelIsoRho[electrons[j]];
	  event->zPairRelIso.second = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[electrons[j]]:event->elePF2PATComRelIsoRho[electrons[i]];
	  event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
	  event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? electrons[j]:electrons[i];
	  closestMass = invMass;
	}
      }
    }
  }

  else if (muons.size() > 1){
    for ( unsigned i {0}; i < muons.size(); i++ ){
      for ( unsigned j{i + 1}; j < muons.size(); j++ ){
        if ( !invertLepCut_) {
  	  if ( event->muonPF2PATCharge[muons[i]] * event->muonPF2PATCharge[muons[j]] >= 0 ) continue;
        }
        else {
  	  if ( !(event->muonPF2PATCharge[muons[i]] * event->muonPF2PATCharge[muons[j]] >= 0) ) continue;
        }
	TLorentzVector lepton1{event->muonPF2PATPX[muons[i]],event->muonPF2PATPY[muons[i]],event->muonPF2PATPZ[muons[i]],event->muonPF2PATE[muons[i]]};
	TLorentzVector lepton2{event->muonPF2PATPX[muons[j]],event->muonPF2PATPY[muons[j]],event->muonPF2PATPZ[muons[j]],event->muonPF2PATE[muons[j]]};
	double invMass{(lepton1 + lepton2).M() -91.1};
	if (std::abs(invMass) < std::abs(closestMass)){
	  event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
	  event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? muons[i]:muons[j];
	  event->zPairRelIso.first = event->muonPF2PATComRelIsodBeta[muons[i]];
	  event->zPairRelIso.second = event->muonPF2PATComRelIsodBeta[muons[j]];
	  event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
	  event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? muons[j]:muons[i];
	  closestMass = invMass;
	}
      }
    }
  }

  return closestMass;
}

float Cuts::getWbosonQuarksCand(AnalysisEvent *event, std::vector<int> jets, int syst){
  float closestWmass {9999.};
  if ( jets.size() > 2 )
    for ( unsigned k{0}; k < jets.size(); k++ ){
      for ( unsigned l{k + 1}; l < jets.size(); l++ ){
	// Now ensure that the leading b jet isn't one of these!
        if ( event->bTagIndex.size() > 0 ) {
  	  if ( event->jetPF2PATBDiscriminator[jets[k]] > bDiscCut_ ){
	    if( event->jetIndex[event->bTagIndex[0]] == jets[k] ) continue;
	  }
	  else if ( event->jetPF2PATBDiscriminator[jets[l]] > bDiscCut_ ){
	    if( event->jetIndex[event->bTagIndex[0]] == jets[l] ) continue;
	  }
        }
        TLorentzVector jetVec1{getJetLVec(event,jets[k],syst)};
        TLorentzVector jetVec2{getJetLVec(event,jets[l],syst)};

	TLorentzVector wQuark1{jetVec1.Px(),jetVec1.Py(),jetVec1.Pz(),jetVec1.E()};
	TLorentzVector wQuark2{jetVec2.Px(),jetVec2.Py(),jetVec2.Pz(),jetVec2.E()};
	double invWbosonMass{(wQuark1 + wQuark2).M() - 80.385};

	if( std::abs(invWbosonMass) < std::abs(closestWmass) ){
	  event->wPairQuarks.first = wQuark1.Pt() > wQuark2.Pt()?wQuark1:wQuark2;
	  event->wPairIndex.first = wQuark1.Pt() > wQuark2.Pt() ? jets[k]:jets[l];
	  event->wPairQuarks.second = wQuark1.Pt() > wQuark2.Pt()?wQuark2:wQuark1;
	  event->wPairIndex.second = wQuark1.Pt() > wQuark2.Pt() ? jets[l]:jets[k];
	  // 	    std::cout << "wQuarks pT: " << event->wPairQuarks.first.Pt() << "/" << event->wPairQuarks.second.Pt() << std::endl;
	  //	    std::cout << "closestMass/invWmass: " << closestWmass << "/" << invWbosonMass << std::endl;
	  closestWmass = invWbosonMass;
	}
      }
    }
  return closestWmass;
}

float Cuts::getTTbarCand(AnalysisEvent *event, std::vector<int> electrons, std::vector<int> muons){

  float closestMass {-1.0};

  for (unsigned i{0}; i < electrons.size(); i++){
    for (unsigned j{0}; j < muons.size(); j++){
      if (event->elePF2PATCharge[electrons[i]] * event->muonPF2PATCharge[muons[j]] > 0) continue;
      TLorentzVector lepton1{event->elePF2PATGsfPx[electrons[i]],event->elePF2PATGsfPy[electrons[i]],event->elePF2PATGsfPz[electrons[i]],event->elePF2PATGsfE[electrons[i]]};
      TLorentzVector lepton2{event->muonPF2PATPX[muons[j]],event->muonPF2PATPY[muons[j]],event->muonPF2PATPZ[muons[j]],event->muonPF2PATE[muons[j]]};
      float invMass{(lepton1+lepton2).M()};
      if( std::abs(invMass) > std::abs(closestMass) ){
        event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
        event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
        closestMass = invMass;
      }
    }
  }
  return closestMass;
}

float Cuts::getTopMass(AnalysisEvent *event){
  TLorentzVector metVec{event->metPF2PATPx,event->metPF2PATPy,0,event->metPF2PATEt};
  TLorentzVector bVec(event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
  float topMass{-1.0};
  if (trileptonChannel_) topMass = (event->wLepton + bVec + metVec).M();
  else topMass = ( bVec + event->wPairQuarks.first + event->wPairQuarks.second ).M();
  return topMass;
}

std::vector<int> Cuts::makeJetCuts(AnalysisEvent *event, int syst, float * eventWeight){

  std::vector<int> jets;
  float mcTag{1.};
  float mcNoTag{1.};
  float dataTag{1.};
  float dataNoTag{1.};
  // b-tagging errors
  float err1{0.};
  float err2{0.};
  float err3{0.};
  float err4{0.};

  //  std::cout << event->eventNum << std::endl << "Jets: " << std::endl;
  for (int i{0}; i < event->numJetPF2PAT; i++){
    //if (std::sqrt(event->jetPF2PATPx[i] * event->jetPF2PATPx[i] + event->jetPF2PATPy[i] * event->jetPF2PATPy[i]) < jetPt_) continue;
    TLorentzVector jetVec{getJetLVec(event,i,syst)};
    // std::cout << getJECUncertainty(sqrt(jetPx*jetPx + jetPy*jetPy), event->jetPF2PATEta[i],syst) << " " << syst << std::endl;
    if (jetVec.Pt() <= jetPt_) continue;
    if (std::abs(jetVec.Eta()) >= jetEta_) continue;

    bool jetId{true};

    // Jet ID == loose
    if ( jetIDDo_ ) {
      if ( std::abs( jetVec.Eta() ) <= 2.7 ) { // for cases where jet eta <= 2.7

        // for all jets with eta <= 2.7
        if ( event->jetPF2PATNeutralHadronEnergyFraction[i] >= 0.99 ) jetId = false;
        if ( event->jetPF2PATNeutralEmEnergyFraction[i] >= 0.99 ) jetId = false;
        if ( ( event->jetPF2PATChargedMultiplicity[i] + event->jetPF2PATNeutralMultiplicity[i] ) <= 1 ) jetId = false;

        //for jets with eta <= 2.40
        if ( std::abs(jetVec.Eta()) <= 2.40 ) {
  	  if ( event->jetPF2PATChargedHadronEnergyFraction[i] <= 0.0 ) jetId = false;
          if ( event->jetPF2PATChargedMultiplicity[i] <= 0.0 ) jetId = false;
          if ( event->jetPF2PATChargedEmEnergyFraction[i] >= 0.99 ) jetId = false;
	}
      }
      else if ( std::abs( jetVec.Eta() ) <= 3.0 && std::abs( jetVec.Eta() ) > 2.40 ) {
        if ( event->jetPF2PATNeutralEmEnergyFraction[i] >= 0.90 ) jetId = false;
        if ( event->jetPF2PATNeutralMultiplicity[i] <= 2 ) jetId = false;
      }
      else if ( std::abs(jetVec.Eta()) > 3.0 ) { // for cases where jet eta > 3.0 and less than 5.0 (or max).
        if ( event->jetPF2PATNeutralEmEnergyFraction[i] >= 0.90 ) jetId = false;
        if ( event->jetPF2PATNeutralMultiplicity[i] <= 10 ) jetId = false;
      }
    }

    if (!jetId) continue;

    double deltaLep{10000};

    if (deltaLep > deltaR(event->zPairLeptons.first.Eta(),event->zPairLeptons.first.Phi(),jetVec.Eta(),jetVec.Phi()))
      deltaLep = deltaR(event->zPairLeptons.first.Eta(),event->zPairLeptons.first.Phi(),jetVec.Eta(),jetVec.Phi());
    if (deltaLep > deltaR(event->zPairLeptons.second.Eta(),event->zPairLeptons.second.Phi(),jetVec.Eta(),jetVec.Phi()))
      deltaLep = deltaR(event->zPairLeptons.second.Eta(),event->zPairLeptons.second.Phi(),jetVec.Eta(),jetVec.Phi());

    if (trileptonChannel_ == true && deltaLep > deltaR(event->wLepton.Eta(),event->wLepton.Phi(),jetVec.Eta(),jetVec.Phi()))
      deltaLep = deltaR(event->wLepton.Eta(),event->wLepton.Phi(),jetVec.Eta(),jetVec.Phi());

    //std::cout << event->jetPF2PATPtRaw[i] << " " << deltaLep << std::endl;
  if (deltaLep < 0.4) continue;

    //    if (event->jetPF2PATdRClosestLepton[i] < 0.5) continue;
   if (isMC_ && makeBTagEffPlots_){
      //Fill eff info here if needed.
      if (std::abs(event->jetPF2PATPID[i]) == 5){ // b-jets
	bTagEffPlots_[0]->Fill(jetVec.Pt(),jetVec.Eta());
	if (event->jetPF2PATBDiscriminator[i] > bDiscCut_)
	  bTagEffPlots_[4]->Fill(jetVec.Pt(),jetVec.Eta());
      }
      if (std::abs(event->jetPF2PATPID[i]) == 4){ // charm
	bTagEffPlots_[1]->Fill(jetVec.Pt(),jetVec.Eta());
	if (event->jetPF2PATBDiscriminator[i] > bDiscCut_)
	  bTagEffPlots_[5]->Fill(jetVec.Pt(),jetVec.Eta());
      }
      if (std::abs(event->jetPF2PATPID[i]) > 0 && std::abs(event->jetPF2PATPID[i]) < 4){ // light jets
	bTagEffPlots_[2]->Fill(jetVec.Pt(),jetVec.Eta());
	if (event->jetPF2PATBDiscriminator[i] > bDiscCut_)
	  bTagEffPlots_[6]->Fill(jetVec.Pt(),jetVec.Eta());
      }
      if (std::abs(event->jetPF2PATPID[i]) == 21){ // gluons
	bTagEffPlots_[3]->Fill(jetVec.Pt(),jetVec.Eta());
	if (event->jetPF2PATBDiscriminator[i] > bDiscCut_)
	  bTagEffPlots_[7]->Fill(jetVec.Pt(),jetVec.Eta());
      }
    }

    jets.emplace_back(i);
    event->jetSmearValue.emplace_back(tempSmearValue_);
    tempSmearValue_ = 1.0; // Reset temporary smear value. Need a cleaner solution than this!

    if (getBTagWeight_ && ( !synchCutFlow_ || (synchCutFlow_ && !trileptonChannel_) )){
      getBWeight(event,jetVec,i,&mcTag,&mcNoTag,&dataTag,&dataNoTag,&err1,&err2,&err3,&err4);
    }
  }
  //Evaluate b-tag weight for event here.
  if (getBTagWeight_){
    float bWeight{(dataNoTag * dataTag)/(mcNoTag * mcTag)};
    if (mcNoTag == 0 || mcTag == 0 || dataNoTag == 0 || dataTag == 0 || mcNoTag != mcNoTag || mcTag != mcTag || dataTag != dataTag || dataNoTag != dataNoTag){
      bWeight = 1.;
    }
    float bWeightErr{std::sqrt( pow(err1+err2,2) + pow(err3 + err4, 2)) * bWeight};
    if (syst == 256)
      bWeight += bWeightErr;
    if (syst == 512)
      bWeight -= bWeightErr;

    *eventWeight *= bWeight;
  }

  return jets;
}

std::vector<int> Cuts::makeBCuts(AnalysisEvent* event, std::vector<int> jets, int syst){

  std::vector<int> bJets;
  for (auto i = 0; i != jets.size(); i++){
    TLorentzVector jetVec{getJetLVec(event,jets[i],syst)};
    if (singleEventInfoDump_) std::cout << __LINE__ << "/" << __FILE__ << ": " << event->jetPF2PATPtRaw[i] << " " << event->jetPF2PATBDiscriminator[jets[i]] << std::endl;
    if (event->jetPF2PATBDiscriminator[jets[i]] <= bDiscSynchCut_ && (synchCutFlow_ && trileptonChannel_)) continue;
    if (event->jetPF2PATBDiscriminator[jets[i]] <= bDiscCut_ ) continue;
    if (jetVec.Eta() >= 2.40) continue;
//    if (event->jetPF2PATEta[ jets[i] ] >= 2.40) continue;
    bJets.emplace_back(i);
  }
  return bJets;
}

std::vector<int> Cuts::makeLooseBCuts(AnalysisEvent* event, std::vector<int> jets, int syst){

  std::vector<int> bJets;

  for (auto i = 0; i != jets.size(); i++){

    TLorentzVector jetVec{getJetLVec(event,jets[i],syst)};
    if (singleEventInfoDump_) std::cout << __LINE__ << "/" << __FILE__ << ": " << event->jetPF2PATPtRaw[i] << " " << event->jetPF2PATBDiscriminator[jets[i]] << std::endl;
    if (event->jetPF2PATBDiscriminator[jets[i]] <= bLooseDiscCut_) continue;
    if (jetVec.Eta() >= 2.40) continue;
//    if (event->jetPF2PATEta[ jets[i] ] >= 2.40) continue;
    bJets.emplace_back(i);
  }
  return bJets;
}

std::vector<int> Cuts::makeCCuts(AnalysisEvent* event, std::vector<int> jets){

  std::vector<int> cJets;
  for (unsigned i{0}; i < jets.size(); i++){
    if (singleEventInfoDump_) std::cout << event->jetPF2PATPtRaw[jets[i]] << " " << event->jetPF2PATCvsLDiscriminator[jets[i]] << std::endl;
//      if (event->jetPF2PATJetCharge[jets[i]] <= 0) continue; // If a negatively charged jet ... I.e. if not a  u or c ...
    if (event->jetPF2PATCvsLDiscriminator[jets[i]] < cVsLDiscCut_) continue; // If doesn't pass c vs light discriminator
    if (event->jetPF2PATBDiscriminator[jets[i]] > bDiscCut_) continue; // If a b jet, continue
    cJets.emplace_back(i);
  }
  return cJets;

}

void Cuts::setTightEle(float,float,float)
{
}

bool Cuts::triggerCuts(AnalysisEvent* event){
  if (skipTrigger_) return true;

  if (synchCutFlow_){
    if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2 > 0 || event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3 > 0 || event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 || event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 || event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 || event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 || event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0)
      return true;
  }

  //MuEG triggers
  bool muEGTrig{false};
  if ( !is2016_ ) {
    if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 > 0 ) muEGTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2 > 0 ) muEGTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
    if ( event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1 > 0 ) muEGTrig = true;
    if ( event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2 > 0 ) muEGTrig = true;
    if ( event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
  }
  else {
    if ( event->eventRun < 280919 ) { // Mu8 leg disabled and non-DZ versions prescaled for Run2016H 
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9 > 0 ) muEGTrig = true;
    }
    if ( event->eventRun >= 280919 ) {
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) muEGTrig = true;
      if ( event->HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) muEGTrig = true;
    }
  }

  //double electron triggers
  bool eeTrig{false};
  if ( !is2016_ ) {
    if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) eeTrig = true;
    if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 ) eeTrig = true;
    if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) eeTrig = true;
  }
  else {
    if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) eeTrig = true;
    if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0 ) eeTrig = true;
    if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5 > 0 ) eeTrig = true;
    if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6 > 0 ) eeTrig = true;
    if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7 > 0 ) eeTrig = true;
    if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8 > 0 ) eeTrig = true;
    if ( event->HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9 > 0 ) eeTrig = true;
  }

  //double muon triggers
  bool mumuTrig{false};
  if ( !is2016_ ) {
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
  }
  else {
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5 > 0 ) mumuTrig = true;
    if ( event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6 > 0 ) mumuTrig = true;
  }

  if ( cutConfTrigLabel_.find("d1") != std::string::npos || cutConfTrigLabel_.find("d2") != std::string::npos || cutConfTrigLabel_.find("d") != std::string::npos ){if (muEGTrig) return true;}
  if ( cutConfTrigLabel_.find("e") != std::string::npos ){ if ( eeTrig && !(muEGTrig || mumuTrig) ) return true; }
  if ( cutConfTrigLabel_.find("m") != std::string::npos ){ if ( mumuTrig && !(eeTrig || muEGTrig) ) return true; }

  return false;
}

// Does event pass MET Filter
bool Cuts::metFilters(AnalysisEvent* event){
  if ( event->Flag_HBHENoiseFilter <= 0 ) return false;
  if ( event->Flag_HBHENoiseIsoFilter <= 0 ) return false;
  if ( event->Flag_EcalDeadCellTriggerPrimitiveFilter <= 0 ) return false;
  if ( event->Flag_goodVertices <= 0 ) return false;
  if ( event->Flag_eeBadScFilter <= 0 ) return false;
  if ( !is2016_ && event->Flag_CSCTightHalo2015Filter <= 0 ) return false;
  if ( is2016_ ) {
    if ( event->Flag_globalTightHalo2016Filter <= 0 ) return false;
    if ( event->Flag_chargedHadronTrackResolutionFilter <= 0 ) return false;
    if ( event->Flag_muonBadTrackFilter <= 0 ) return false;
    if ( event->Flag_ecalLaserCorrFilter <= 0 ) return false;
  }
  else return true;
}

//Does the cuts required for the synchronisation.
bool Cuts::synchCuts(AnalysisEvent* event, float *eventWeight){
  //Trigger stuff would go first, but in MC at least (what I'm starting with) I don't have that saved. Idiot.

  //This is to make some skims for faster running. Do lepSel and save some files.
  if(postLepSelTree_ && synchCutFlow_) {
    if (postLepSelTree_->GetEntriesFast() < 40000) postLepSelTree_->Fill();
    else {
      if (postLepSelTree2_->GetEntriesFast() < 40000) postLepSelTree2_->Fill();
      else postLepSelTree3_->Fill();

    }
  }

  if (singleEventInfoDump_){
//    std::cout << std::setprecision(6) << std::fixed;
  }

  if (makeEventDump_) dumpToFile(event,9);

  if ( !trileptonChannel_ ){
    synchCutFlowHist_->Fill(0.5, *eventWeight); // Total events
//    std::cout << std::setprecision(6) << std::fixed;

    if (!triggerCuts(event)) return false;
    if (!metFilters(event)) return false;

    synchCutFlowHist_->Fill(1.5, *eventWeight); // Trigger cuts - Step 0
    event->electronIndexTight = getTightEles(event);
    event->muonIndexTight = getTightMuons(event);
    synchNumEles_->Fill(event->electronIndexTight.size());
    synchNumMus_->Fill(event->muonIndexTight.size());

    if (event->electronIndexTight.size() != numTightEle_) return false;
    if (event->muonIndexTight.size() != numTightMu_) return false;
    if ( (event->electronIndexTight.size() + event->muonIndexTight.size()) != 2 ) return false;

    synchCutFlowHist_->Fill(2.5, *eventWeight); // 2 Tight Leptons - step 1
    if ( (event->electronIndexTight.size() != getLooseEles(event).size()) || (event->muonIndexTight.size() != getLooseMuons(event).size()) ) return false;
    synchCutFlowHist_->Fill(3.5, *eventWeight); // Lepton Veto - step 2

    if (std::abs(getDileptonZCand(event,event->electronIndexTight,event->muonIndexTight)) > invZMassCut_) return false;
    synchCutFlowHist_->Fill(4.5, *eventWeight); // Z veto - step 3

    *eventWeight *= getLeptonWeight(event,0);

    event->jetIndex = makeJetCuts(event, 0, eventWeight);
    if (event->jetIndex.size() < numJets_) return false;
    if (event->jetIndex.size() > maxJets_) return false;
    synchCutFlowHist_->Fill(5.5, *eventWeight); // jet selection - step 4

    event->bTagIndex = makeBCuts(event,event->jetIndex, 0);
    if (event->bTagIndex.size() < numbJets_) return false;
    if (event->bTagIndex.size() > maxbJets_) return false;
    synchCutFlowHist_->Fill(6.5, *eventWeight); // b-jet selection - step 5

    synchCutFlowHist_->Fill(7.5, *eventWeight); // no Met cut applied

    double wMass = getWbosonQuarksCand(event,event->jetIndex,0);

    if ( std::abs(wMass) > invWMassCut_ ) return false;

    TLorentzVector bVec(event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);

//    double topMass = (bVec + event->wPairQuarks.first + event->wPairQuarks.second).M();
    double topMass = getTopMass(event);
//    std::cout << topMass << " / " << getTopMass(event) << std::endl;

/*    if ( topMass > 130 ){
       if ( std::abs(event->jetPF2PATPID[event->jetIndex[event->bTagIndex[0]]]) == 5)
       std::cout << "jet PID/b mass/w Mass/leading b pT: " << event->jetPF2PATPID[event->jetIndex[event->bTagIndex[0]]] << " / " << bVec.M() << " / " << (event->wPairQuarks.first + event->wPairQuarks.second).M() << " / " << event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]] << "\t " << std::endl;
    }*/

    return true;
  }

  synchCutFlowHist_->Fill(0.5, *eventWeight); // Total events
//  std::cout << std::setprecision(6) << std::fixed;

  if (makeEventDump_) dumpToFile(event,0);

  if (!triggerCuts(event)) return false;
//  if( !skipTrigger_ && isMC_ ) {
//    if( !is2016_ ) *eventWeight *= get2015TriggerSF (systToRun);
//    else *eventWeight *= get2016TriggerSF (systToRun);
//  }

  synchCutFlowHist_->Fill(1.5, *eventWeight); // Trigger cuts - Step 0

  // Check number of leptons is correct
  if (singleEventInfoDump_) std::cout << "Correct number of leptons and loose: " << getLooseEles(event).size() << " " << getLooseMuons(event).size() << std::endl;

  event->electronIndexTight = getTightEles(event);
  event->muonIndexTight = getTightMuons(event);

  // Check at exactly three tight leptons
  synchNumEles_->Fill(event->electronIndexTight.size());
  synchNumMus_->Fill(event->muonIndexTight.size());
  if (event->electronIndexTight.size() != numTightEle_) return false;
  if (event->muonIndexTight.size() != numTightMu_) return false;
  if ( (event->electronIndexTight.size() + event->muonIndexTight.size()) != 3 ) return false;
  synchCutFlowHist_->Fill(2.5, *eventWeight); // 3 Tight Leptons - step 1

//  for ( int i = 0; i != event->numMuonPF2PAT; i++ ) {
//    topMassEventDump_ << "EvtNb="<< event->eventNum << "mu_pt=" << event->muonPF2PATPt[i] <<" mu_eta=" << event->muonPF2PATEta[i] << " mu_phi=" << event->muonPF2PATPhi[i] << " mu_iso=" << event->muonPF2PATComRelIsodBeta[i] << std::endl;
//  }

  //loose lepton veto
  //  int looseLeps = getLooseLepsNum(event);
  //  if (isMC_ && looseLeps < 2) return false;
  //  if (!isMC_ && looseLeps < 3) return false;

  if ( (event->electronIndexTight.size() != getLooseEles(event).size()) || (event->muonIndexTight.size() != getLooseMuons(event).size()) ) return false;
  if (singleEventInfoDump_) std::cout << " and passes veto too." << std::endl;
  if ( makeEventDump_ ) dumpToFile(event,2);
  synchCutFlowHist_->Fill(3.5, *eventWeight); // Lepton Veto - step 2

  // Z selection
  if (singleEventInfoDump_) std::cout << "Z mass: " << getZCand(event,event->electronIndexTight,event->muonIndexTight) << std::endl;
  if (std::abs(getZCand(event,event->electronIndexTight,event->muonIndexTight)) > invZMassCut_) return false;
  synchCutFlowHist_->Fill(4.5, *eventWeight); // Z veto - step 3

//  *eventWeight *= getLeptonWeight(event,0);

  //Add in extra steps here.

  // Jet selection
  event->jetIndex = makeJetCuts(event, 0, eventWeight);
  if (singleEventInfoDump_) std::cout << "Number of jets: " << event->jetIndex.size() << std::endl;
  if (event->jetIndex.size() < 1) return false;
  if (makeEventDump_) dumpToFile(event,4);
  //  std::cout << event->jetIndex.size() << std::endl;
  synchCutFlowHist_->Fill(5.5, *eventWeight); // jet selection - step 4

  // bTag selection
  event->bTagIndex = makeBCuts(event,event->jetIndex, 0);
  synchCutTopMassHist_->Fill(getTopMass(event)); // Plot top mass distribution for all top candidates - all sanity checks done, Z mass exists, got b jets too.
  if (singleEventInfoDump_) std::cout << "One bJet: " << event->bTagIndex.size() << std::endl;
  if (event->bTagIndex.size() != 1) return false;
  synchCutFlowHist_->Fill(6.5, *eventWeight); // b-jet selection - step 5

  // MET cut
  if (singleEventInfoDump_) std::cout << "MET: " << event->metPF2PATPt << std::endl;
  //  if (event->metPF2PATPt < metCut_) return false; // MET Cut not currently applied
  synchCutFlowHist_->Fill(7.5, *eventWeight);

  // mTW cut
  if (singleEventInfoDump_) std::cout << "mTW: " << event->metPF2PATPt << std::endl;
  double mtW{std::sqrt( 2*event->metPF2PATPt*event->wLepton.Pt()*(1-cos(event->metPF2PATPhi - event->wLepton.Phi())) )};
  if ( mtW < mTWCutSynch_ ) return false;
  synchCutFlowHist_->Fill(8.5, *eventWeight); // mWT cut - step 6

  // Top Mass cut
  if (singleEventInfoDump_) std::cout << "top mass cut: " << getTopMass(event)  << std::endl;
  double topMass (getTopMass(event));

  if (topMass > TopMassCutUpper_ || topMass < TopMassCutLower_) return false;
  synchCutFlowHist_->Fill(9.5, *eventWeight); // top mass cut - step 7

  if (!metFilters(event)) return false;
  synchCutFlowHist_->Fill(10.5, *eventWeight);

  if (singleEventInfoDump_) std::cout << "Passes all cuts! Yay!" << std::endl;
  if (makeEventDump_) dumpToFile(event,6);
//  if (singleEventInfoDump_) std::cout << std::setprecision(1) << std::fixed;
  return true;
}

TH1F* Cuts::getSynchCutFlow(){
  std::cout << "Eles: " << numTightEle_ << " Muons: " << numTightMu_ << std::endl;
  char const *names[] {"Total Events","Trigger","3 Leptons", "Lepton Veto", "zMass","1 jet","1 b-tag","MET","mTW", "topMass", "metFilters"};
  for (unsigned i{1}; i < 12; ++i){
    std::cout << names[i-1] << ": " << synchCutFlowHist_->GetBinContent(i) << std::endl;
  }
  std::cout << "The number of leptons in the passing trigger events:" << std::endl;
  std::cout << "Leps\tEle\tMuon" << std::endl;
  for (unsigned i{0}; i < 12; i++){
    std::cout << i << "\t" << synchNumEles_->GetBinContent(i) << "\t" << synchNumMus_->GetBinContent(i) << std::endl;
  }
  char const *labels[] = {"In", "ID", "PtEtaIso", "chi2","tklay","DBPV","TrackHit","MuHits","PixHits","MtStats","DZPV"};
  for (unsigned i{1}; i < 12; ++i){
    std::cout << labels[i-1] << ": \t" << synchMuonCutFlow_->GetBinContent(i) << std::endl;
  }

  return synchCutFlowHist_;
}

//Method used for the synchronisation. Mimics the IPHC preselection skims.
int Cuts::getLooseLepsNum(AnalysisEvent* event){
  return getLooseElecs(event) + getLooseMus(event);
}

int Cuts::getLooseElecs(AnalysisEvent* event){
  int looseLeps{0};
  for (int i{0}; i < event->numElePF2PAT; i++){
    if (!event->elePF2PATIsGsf[i]) continue;
    TLorentzVector tempVec{event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]};
    if (tempVec.Pt() < 20) continue;
    if (std::abs(tempVec.Eta()) > 2.5)continue;
    looseLeps++;
  }
  return looseLeps;
}

int Cuts::getLooseMus(AnalysisEvent* event){
  int looseLeps{0};
  for (int i{0}; i < event->numMuonPF2PAT; i++){
    if (!event->muonPF2PATGlobalID[i]||!event->muonPF2PATTrackID[i]) continue;
    if (event->muonPF2PATPt[i] < 20) continue;
    if (std::abs(event->muonPF2PATEta[i]) > 2.4) continue;
    looseLeps++;
  }
  return looseLeps;
}


//First tentative attempt at doing the background isolation.
bool Cuts::invertIsoCut(AnalysisEvent* event,float *eventWeight,std::map<std::string,Plots*> plotMap, TH1F* cutFlow){

  if (trileptonChannel_ == false){
    std::cout << "Invert Iso Cut is not avaliable for the dilepton channel." << std::endl;
    return false;
  }
  //Check there are exactly 2 tight leptons with the correct isolation cut.
  event->electronIndexTight = getTightEles(event);
  event->muonIndexTight = getTightMuons(event);

  if ((event->electronIndexTight.size() + event->muonIndexTight.size()) != 2) return false;
  if (event->electronIndexTight.size() == 1) return false;

  //Check they are a valid zCandidate. (Just call the zCand method here? - No, that won't actually work.)
  float invMass{100.};
  if (numTightEle_ > 1){
    if (event->electronIndexTight.size() < 2) return false;
    if (event->elePF2PATCharge[event->electronIndexTight[0]] * event->elePF2PATCharge[event->electronIndexTight[1]] > 0) return false;
    TLorentzVector lep1{event->elePF2PATGsfPx[event->electronIndexTight[0]],event->elePF2PATGsfPy[event->electronIndexTight[0]],event->elePF2PATGsfPz[event->electronIndexTight[0]],event->elePF2PATGsfE[event->electronIndexTight[0]]};
    TLorentzVector lep2{event->elePF2PATGsfPx[event->electronIndexTight[1]],event->elePF2PATGsfPy[event->electronIndexTight[1]],event->elePF2PATGsfPz[event->electronIndexTight[1]],event->elePF2PATGsfE[event->electronIndexTight[1]]};
    event->zPairLeptons.first = lep1.Pt() > lep2.Pt() ? lep1:lep2;
    event->zPairIndex.first = lep1.Pt() > lep2.Pt() ? event->electronIndexTight[0] : event->electronIndexTight[1];
    event->zPairLeptons.second = lep1.Pt() > lep2.Pt() ? lep2:lep1;
    event->zPairIndex.second = lep1.Pt() > lep2.Pt() ? event->electronIndexTight[1] : event->electronIndexTight[0];
    invMass = (event->zPairLeptons.first + event->zPairLeptons.second).M() -91.1;
  }
  else{
    if (event->muonIndexTight.size() < 2) return false;
    if (event->muonPF2PATCharge[event->muonIndexTight[0]] * event->muonPF2PATCharge[event->muonIndexTight[1]] > 0) return false;
    event->zPairLeptons.first = TLorentzVector(event->muonPF2PATPX[event->muonIndexTight[0]],event->muonPF2PATPY[event->muonIndexTight[0]],event->muonPF2PATPZ[event->muonIndexTight[0]],event->muonPF2PATE[event->muonIndexTight[0]]);
    event->zPairIndex.first = event->muonIndexTight[0];
    event->zPairLeptons.second = TLorentzVector(event->muonPF2PATPX[event->muonIndexTight[1]],event->muonPF2PATPY[event->muonIndexTight[1]],event->muonPF2PATPZ[event->muonIndexTight[1]],event->muonPF2PATE[event->muonIndexTight[1]]);
    event->zPairIndex.second = event->muonIndexTight[1];
    invMass = (event->zPairLeptons.first + event->zPairLeptons.second).M() -91.1;
  }

  //Get rev iso candidates.
  std::vector<int> invIsoEle{getInvIsoEles(event)};
  std::vector<int> invIsoMus{getInvIsoMuons(event)};

  //Check we have the right number of leptons.
  if (event->electronIndexTight.size() + invIsoEle.size() != numTightEle_) return false;
  if (event->muonIndexTight.size() + invIsoMus.size() != numTightMu_) return false;

  //Debugging
  /*  std::cout << event->numElePF2PAT << " " << event->electronIndexTight.size() << " " << invIsoEle.size();
  std::cout << " tight index: ";
  for (unsigned i = 0; i < event->electronIndexTight.size(); i++) std:: cout << " " << event->electronIndexTight[i];
  std::cout << " inv index: ";
  for (unsigned i = 0; i < invIsoEle.size(); i++) std::cout << " " << invIsoEle[i] << " " << "relIso: " << event->elePF2PATComRelIsoRho[invIsoEle[i]]/event->elePF2PATPT[invIsoEle[i]] ;
  std::cout << std::endl;*/

  //Put extra lepton into W boson thing.
  if (invIsoEle.size() == 1){
    event->wLepton = TLorentzVector(event->elePF2PATGsfPx[invIsoEle[0]],event->elePF2PATGsfPy[invIsoEle[0]],event->elePF2PATGsfPz[invIsoEle[0]],event->elePF2PATGsfE[invIsoEle[0]]);
    event->wLepIndex = invIsoEle[0];
    event->wLeptonRelIso = event->elePF2PATComRelIsoRho[invIsoEle[0]];
  }
  else{
    event->wLepton = TLorentzVector(event->muonPF2PATPX[invIsoMus[0]],event->muonPF2PATPY[invIsoMus[0]],event->muonPF2PATPZ[invIsoMus[0]],event->muonPF2PATE[invIsoMus[0]]);
    event->wLepIndex = invIsoMus[0];
    event->wLeptonRelIso = event->muonPF2PATComRelIsodBeta[invIsoMus[0]];
  }

  if(doPlots_) plotMap["lepSel"]->fillAllPlots(event,*eventWeight);
  if(doPlots_||fillCutFlow_) cutFlow->Fill(0.5,*eventWeight);

  if (std::abs(invMass) > invZMassCut_) return false;
  if (doPlots_) plotMap["zMass"]->fillAllPlots(event,*eventWeight);
  if(doPlots_||fillCutFlow_) cutFlow->Fill(1.5,*eventWeight);
  return true;

}

std::vector<int> Cuts::getInvIsoEles(AnalysisEvent* event) {
  std::vector<int> electrons;
  for (int i{0}; i < event->numElePF2PAT; i++){
    bool same{false};
    /*    for (int j = 0; j < event->electronIndexTight.size(); j++){
      if (i == event->electronIndexTight[j]) same = true;
      }*/
    if (same) continue;
    if (!event->elePF2PATIsGsf[i]) continue;
    TLorentzVector tempVec{event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]};
    if (tempVec.Pt() < tightElePt_) continue;
    if (std::abs(tempVec.Eta()) > tightEleEta_)continue;
    if (std::abs(event->elePF2PATBeamSpotCorrectedTrackD0[i]) > tightEled0_)continue;
    if (event->elePF2PATMissingInnerLayers[i] > tightEleMissLayers_) continue;
    if (!event->elePF2PATPhotonConversionVeto[i] && tightEleCheckPhotonVeto_)continue;
    if ( event->elePF2PATMVAcategory[i] == 0 && (event->elePF2PATMVA[i] < tightEleMVA0_) ) continue;
    if ( event->elePF2PATMVAcategory[i] == 1 && (event->elePF2PATMVA[i] < tightEleMVA1_) ) continue;
    if ( event->elePF2PATMVAcategory[i] == 2 && (event->elePF2PATMVA[i] < tightEleMVA2_) ) continue;
    if (event->elePF2PATComRelIsoRho[i] < tightEleRelIso_)continue;

    electrons.emplace_back(i);
  }
  return electrons;
}

std::vector<int> Cuts::getInvIsoMuons(AnalysisEvent* event){
  std::vector<int> muons;
  for (int i{0}; i < event->numMuonPF2PAT; i++){
    if (!event->muonPF2PATGlobalID[i] && !event->muonPF2PATTrackID[i]) continue;
    if (event->muonPF2PATPt[i] < tightMuonPt_) continue;
    if (std::abs(event->muonPF2PATEta[i]) > tightMuonEta_) continue;
    if (event->muonPF2PATComRelIsodBeta[i] < tightMuonRelIso_) continue;

    //Do a little test of muon id stuff here.
    if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] > 10.) continue;
    if (std::abs(event->muonPF2PATDBInnerTrackD0[i]) > 0.2) continue;
    if (event->muonPF2PATNChambers[i] < 2) continue;

    muons.emplace_back(i);
  }
  return muons;
}

double Cuts::getChiSquared( double wMass, double topMass ){

  double topMassTerm = ( (topMass - 173.21) / 50.0 );
  double wMassTerm = ( (wMass - 80.3585 ) / 10.0 );
  return std::sqrt( topMassTerm*topMassTerm + wMassTerm*wMassTerm );
}

bool Cuts::ttbarCuts(AnalysisEvent* event, float *eventWeight, std::map<std::string,Plots*> plotMap, TH1F* cutFlow, int systToRun){

  if( !skipTrigger_ ) {
    if ( !is2016_ ) if (!triggerCuts(event)) return false; // Do trigger on MC and data for 2015
    if ( !isMC_ && is2016_ ) if (!triggerCuts(event)) return false; // Do trigger for data and 2016, exclude MC.

    if( !is2016_ && isMC_ ) *eventWeight *= get2015TriggerSF (systToRun); // Apply SFs to MC if 2015
    else if ( is2016_ && isMC_ ) *eventWeight *= get2016TriggerSF (systToRun); // Apply data efficiencies onto MC if 2016 due to incomplete trigger menu
  }

  if (!metFilters(event)) return false;

  if (!(makeLeptonCuts(event,eventWeight,plotMap,cutFlow,systToRun,true))) return false;

  event->jetIndex = makeJetCuts(event, systToRun, eventWeight);
  if (event->jetIndex.size() < numJets_) return false;
  if (event->jetIndex.size() > maxJets_) return false;
  if (doPlots_||fillCutFlow_) cutFlow->Fill(2.5,*eventWeight);
  if (doPlots_) plotMap["jetSel"]->fillAllPlots(event,*eventWeight);

  event->bTagIndex = makeBCuts(event,event->jetIndex, systToRun);
  if ( looseBjetVeto_ == 1 ) event->bTagLooseIndex = makeLooseBCuts(event,event->jetIndex, systToRun);
  if (event->bTagIndex.size() < numbJets_) return false;
  if (event->bTagIndex.size() > maxbJets_) return false;
  if ( looseBjetVeto_ == 1  && (event->bTagIndex.size() != event->bTagLooseIndex.size()) ) return false;

  if (doPlots_) plotMap["bTag"]->fillAllPlots(event,*eventWeight);
  if (doPlots_||fillCutFlow_) cutFlow->Fill(3.5,*eventWeight);

  return true;
}

//For synchronisation I am dumping the lepton information here.
void Cuts::dumpLeptonInfo(AnalysisEvent* event){
//  std::cout << std::setprecision(6) << std::fixed;
//  std::cerr << std::setprecision(6) << std::fixed;
  //Dump electron info first
  event->electronIndexTight = getTightEles(event);
  event->muonIndexTight = getTightMuons(event);
  std::cout << "Electrons: found " << event->electronIndexTight.size() << std::endl;
  for (unsigned i{0}; i < event->electronIndexTight.size(); i++){
    TLorentzVector tempVec{event->elePF2PATGsfPx[event->electronIndexTight[i]],event->elePF2PATGsfPy[event->electronIndexTight[i]],event->elePF2PATGsfPz[event->electronIndexTight[i]],event->elePF2PATGsfE[event->electronIndexTight[i]]};
    //    std::cout << i << " | " << event->elePF2PATPT[event->electronIndexTight[i]];
    std::cout << " | " << tempVec.Pt();
    //std::cout << " | " << event->elePF2PATEta[event->electronIndexTight[i]];
    std::cout << " | " << tempVec.Eta();
    std::cout << " | " << tempVec.Phi();
    std::cout << " | " << event->elePF2PATPhi[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATD0PV[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATMVA[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATMVAcategory[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATMissingInnerLayers[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATComRelIsoRho[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATPhotonConversionVeto[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATCharge[event->electronIndexTight[i]];
    /*std::cout << " | " << event->elePF2PATRhoIso[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATAEff03[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATChHadIso[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATNtHadIso[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATGammaIso[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATComRelIsodBeta[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATComRelIso[event->electronIndexTight[i]]; */
    std::cout << std::endl;
  }

  std::cout << "Muons: found " << event->muonIndexTight.size() << std::endl;
  for (unsigned i{0}; i <  event->muonIndexTight.size(); i++){
    std::cout << i;
    std::cout << " | " << event->muonPF2PATPt[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATEta[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATPhi[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATChi2[event->muonIndexTight[i]]/event->muonPF2PATNDOF[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATTkLysWithMeasurements[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATMuonNHits[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATDBPV[event->muonIndexTight[i]];
    std::cout << " | " << event->pvZ - event->muonPF2PATVertZ[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATVldPixHits[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATMatchedStations[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATComRelIsodBeta[event->muonIndexTight[i]];
    std::cout << " | " << event->muonPF2PATCharge[event->muonIndexTight[i]];
    std::cout << std::endl;
  }
  int numbJets{event->numJetPF2PAT > 4?4:event->numJetPF2PAT};
  std::cout << "Jets: " << event->numJetPF2PAT << std::endl;
  for (int i{0}; i < numbJets; i++){
    std::cout << i;
    //    std::cout << " | " << event->jetPF2PATPt[i];
    //std::cout << " | " << std::sqrt(event->jetPF2PATPx[i] * event->jetPF2PATPx[i] + event->jetPF2PATPy[i] * event->jetPF2PATPy[i]);
    // std::cout << " | " << event->jetPF2PATEt[i];
    std::cout << " | " << event->jetPF2PATPtRaw[i];
    //std::cout << " | " << event->jetPF2PATUnCorPt[i];
    std::cout << " | " << event->jetPF2PATEta[i];
    std::cout << " | " << event->jetPF2PATPhi[i];
    std::cout << " | " << event->jetPF2PATNConstituents[i];
    std::cout << " | " << (event->jetPF2PATNeutralHadronEnergyFraction[i] < 0.99 && event->jetPF2PATNeutralEmEnergyFraction[i] < 0.99) && ((std::abs(event->jetPF2PATEta[i]) > 2.4) || (event->jetPF2PATChargedEmEnergyFraction[i] < 0.99 && event->jetPF2PATChargedHadronEnergyFraction[i] > 0. && event->jetPF2PATChargedMultiplicity[i] > 0.));
    std::cout << " | " << event->jetPF2PATdRClosestLepton[i];
    std::cout << std::endl;
  }
  std::cout << "MET: " << event->metPF2PATEt << " | " << event->metPF2PATPt << std::endl;

//  std::cout << std::setprecision(1) << std::fixed;
//  std::cerr << std::setprecision(1) << std::fixed;
}

void Cuts::dumpLooseLepInfo(AnalysisEvent* event){
//  std::cout << std::setprecision(6) << std::fixed;
//  std::cerr << std::setprecision(6) << std::fixed;

  std::cout << "Electrons: " << event->numElePF2PAT << std::endl;
  for (int i{0}; i < event->numElePF2PAT; i++){

    TLorentzVector tempVec{event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]};
    std::cout << i;
    std::cout << " | " << tempVec.Pt();
    std::cout << " | " << tempVec.Eta();
    std::cout << " | " << event->elePF2PATSCEta[i];
    std::cout << " | " << event->elePF2PATComRelIsoRho[i];
    std::cout << std::endl;
  }
  std::cout << "Muons: " << event->numMuonPF2PAT << std::endl;
  for (int i{0}; i < event->numMuonPF2PAT; i++){
    std::cout << i;
    std::cout << " | " << event->muonPF2PATPt[i];
    std::cout << " | " << event->muonPF2PATEta[i];
    std::cout << " | " << event->muonPF2PATComRelIsodBeta[i];
    std::cout << " | " << event->muonPF2PATGlobalID[i];
    std::cout << " | " << event->muonPF2PATTrackID[i];
    std::cout << std::endl;
  }

//  std::cout << std::setprecision(1) << std::fixed;
//  std::cerr << std::setprecision(1) << std::fixed;
}

double Cuts::deltaR(float eta1, float phi1, float eta2, float phi2){
  double dEta{eta1-eta2};
  double dPhi{phi1-phi2};
  while (std::abs(dPhi) > M_PI){
    dPhi += (dPhi > 0.? -2*M_PI:2*M_PI);
  }
  //  if(singleEventInfoDump_)  std::cout << eta1 << " " << eta2 << " phi " << phi1 << " " << phi2 << " ds: " << eta1-eta2 << " " << phi1-phi2 << " dR: " << std::sqrt((dEta*dEta)+(dPhi*dPhi)) << std::endl;
  return std::sqrt((dEta*dEta)+(dPhi*dPhi));
}

//For dumping contents of step 4. Bit more complicated than old, so doing it elsewhere.
void Cuts::dumpToFile(AnalysisEvent* event, int step){

  std::vector<TLorentzVector> tempLepVec;
  unsigned triggerFlag[3]{};
  std::string channel{"nan"};
  std::pair<int,int> leadingLeptons[3]{}; // Initalise as empty

  // lepton ID for step0?
  event->electronIndexTight = getTightEles(event);
  event->muonIndexTight = getTightMuons(event);

  if ( step == 9 ) { // Used for 2016 debug synch
    for ( int i = 0; i < event->numElePF2PAT; i++) {
//      step9EventDump_.precision(3);
      step9EventDump_ << "|" << event->eventNum << "|" << event->elePF2PATPT[i] << "|" << event->elePF2PATCutIdVeto[i] << "|" << event->elePF2PATCutIdTight[i] << "|" << std::endl;
    }
  }

  if ( step == 0 ) { // Used for 2015/2016 synch

    // Get trigger bit setup
    if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2 > 0 || event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) triggerFlag[2] = 1; // Set Z=1 if MuonEG trigger fires
    if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 || event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 || event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) triggerFlag[1] = 1; // Set Y=1 if DoubleEG trigger fires
    if ( event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1 > 0 || event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0 || event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0 ) triggerFlag[0] = 1; // Set X=1 if DoubleMuon trigger fires

  // Get leading 3 leptons pT
  // Search over electrons
    for ( auto electronIt = event->electronIndexTight.begin(); electronIt != event->electronIndexTight.end(); electronIt++) {
    //    for ( int electronIt = 0; electronIt != event->numElePF2PAT; electronIt++ ) {
    float elePt{event->elePF2PATPT[*electronIt]};
    float itPt[3]{};
    for ( unsigned j{0}; j != 3; j++ ){
      if ( leadingLeptons[j].second == 1 ) itPt[j] = event->elePF2PATPT[leadingLeptons[j].first];
      else if ( leadingLeptons[j].second == 2 ) itPt[j] = event->muonPF2PATPt[leadingLeptons[j].first];
    }

    if ( elePt > itPt[2] && elePt <= itPt[1] ){
      leadingLeptons[2] = std::make_pair(*electronIt,1);
    }
    else if ( elePt > itPt[1] && elePt <= itPt[0] ){
      leadingLeptons[2] = leadingLeptons[1];
      leadingLeptons[1] = std::make_pair(*electronIt,1);
    }
    else if ( elePt > itPt[0] ){
      leadingLeptons[2] = leadingLeptons[1];
      leadingLeptons[1] = leadingLeptons[0];
      leadingLeptons[0] = std::make_pair(*electronIt,1);
    }
  }

  // Search over muons
    for ( auto muonIt = event->muonIndexTight.begin(); muonIt != event->muonIndexTight.end(); muonIt++) {
    //    for ( int  muonIt = 0; muonIt != event->numMuonPF2PAT; muonIt++ ) {
    float muonPt{event->muonPF2PATPt[*muonIt]};
    float itPt[3]{};
    for ( unsigned j{0}; j != 3; j++ ){
      if ( leadingLeptons[j].second == 1 ) itPt[j] = event->elePF2PATPT[leadingLeptons[j].first];
      else if ( leadingLeptons[j].second == 2 ) itPt[j] = event->muonPF2PATPt[leadingLeptons[j].first];
    }

    if ( muonPt > itPt[2] && muonPt <= itPt[1] ){
      leadingLeptons[2] = std::make_pair(*muonIt,2);
    }
    else if ( muonPt > itPt[1] && muonPt <= itPt[0] ){
      leadingLeptons[2] = leadingLeptons[1];
      leadingLeptons[1] = std::make_pair(*muonIt,2);
    }
    else if ( muonPt > itPt[0] ){
      leadingLeptons[2] = leadingLeptons[1];
      leadingLeptons[1] = leadingLeptons[0];
      leadingLeptons[0] = std::make_pair(*muonIt,2);
    }
  }

    // Setup channel label
    int numEles{0};
    int numMuons{0};
    for ( unsigned i{0}; i != 3; ++i ){
    if (leadingLeptons[i].second == 1) numEles++;
    if (leadingLeptons[i].second == 2) numMuons++;
  }

  if (  numEles == 3 &&  numMuons == 0 ) channel = "eee";
    else if ( numEles == 2 &&  numMuons == 1 ) channel = "eem";
    else if ( numEles == 1 &&  numMuons == 2 ) channel = "emm";
    else if ( numEles == 0 &&  numMuons == 3 ) channel = "mmm";
  }

  if ( step > 2 && trileptonChannel_ == true){
    if (event->zPairLeptons.first.Pt() < event->wLepton.Pt()){
      tempLepVec.emplace_back(event->wLepton);
      tempLepVec.emplace_back(event->zPairLeptons.first);
      tempLepVec.emplace_back(event->zPairLeptons.second);
    }
    else if (event->wLepton.Pt() > event->zPairLeptons.second.Pt()){
      tempLepVec.emplace_back(event->zPairLeptons.first);
      tempLepVec.emplace_back(event->wLepton);
      tempLepVec.emplace_back(event->zPairLeptons.second);
    }
    else {
      tempLepVec.emplace_back(event->zPairLeptons.first);
      tempLepVec.emplace_back(event->zPairLeptons.second);
      tempLepVec.emplace_back(event->wLepton);
    }
  }
  if (step==2 && trileptonChannel_ == true){
    int muUsed{0};
    for (int i{0}; i < event->numElePF2PAT; i++){
      TLorentzVector tempVec{event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]};
      for (int j{muUsed}; j < event->numMuonPF2PAT; j++){
	if (tempVec.Pt() > event->muonPF2PATPt[i]){
	  tempLepVec.emplace_back(tempVec);
	  break;
	}
	else {
	  tempLepVec.emplace_back(event->muonPF2PATPX[j],event->muonPF2PATPY[j],event->muonPF2PATPZ[j],event->muonPF2PATE[j]);
	  muUsed++;
	  if (tempLepVec.size() > 2) break;
	}
	if (tempLepVec.size() > 2) break;
      }
    }
    if (tempLepVec.size() < 3){
      for (int j{0}; j < 3; j++){
	tempLepVec.emplace_back(event->muonPF2PATPX[j],event->muonPF2PATPY[j],event->muonPF2PATPZ[j],event->muonPF2PATE[j]);
      }
    }
  }


  switch (step) {
  case 0:
    step0EventDump_.precision(3);
    step0EventDump_ << "|" << event->eventNum << "|" << triggerFlag[0] << triggerFlag[1] << triggerFlag[2] << "|" << channel << "|";
  case 2:
    step2EventDump_ << event->eventRun << " " << event->eventNum << " ";
    break;
  case 4:
    step4EventDump_ << event->eventRun << " " << event->eventNum << " ";
    break;
  case 6:
    step6EventDump_ << event->eventRun << " " << event->eventNum << " ";
    break;
  }
  for (unsigned i{0}; i < 3; i++){
    switch (step) {
    case 0:
      step0EventDump_.precision(3);
      if ( leadingLeptons[i].second == 1 ) step0EventDump_ << event->elePF2PATPT[leadingLeptons[i].first] << "|";
      else if ( leadingLeptons[i].second == 2 )	step0EventDump_ << event->muonPF2PATPt[leadingLeptons[i].first] << "|";
      else step0EventDump_ << 0 << "|";
      break;
    case 2:
      step2EventDump_ << tempLepVec[i].Pt() << " " << tempLepVec[i].Eta() << " ";
      break;
    case 4:
      step4EventDump_ << tempLepVec[i].Pt() << " " << tempLepVec[i].Eta() << " ";
      break;
    case 6:
      step6EventDump_ << tempLepVec[i].Pt() << " " << tempLepVec[i].Eta() << " ";
      break;
    }
  }

  // Rel iso for step0
  for (unsigned i{0}; i < 3; i++){
    switch (step) {
    case 0:
      step0EventDump_.precision(3);
      if ( leadingLeptons[i].second == 1 ) step0EventDump_ << event->elePF2PATComRelIsoRho[leadingLeptons[i].first] << "|";
      else if ( leadingLeptons[i].second == 2 ) step0EventDump_ << event->muonPF2PATComRelIsodBeta[leadingLeptons[i].first] << "|";
      else step0EventDump_ << 0 << "|";
      break;
    }
  }

// Id flag superfluous as now we are doing selection on leading leptons
/*  switch (step) {
  case 0:
    for (unsigned i{0}; i < 3; i++){

      bool IdFlag{false};
      if ( leadingLeptons[i].second == 1 ) {
	for ( unsigned j{0}; j != event->electronIndexTight.size(); j++ ) {
	  if ( event->electronIndexTight[j] == leadingLeptons[i].first ) IdFlag = 1;
	}
      }
      else if ( leadingLeptons[i].second == 2 ) {
	for ( unsigned j{0}; j != event->muonIndexTight.size(); j++ ) {
	  if ( event->muonIndexTight[j] == leadingLeptons[i].first ) IdFlag = 1;
	}
      }
      step0EventDump_ << IdFlag << "|";
    }
    break;
  }*/

  float tempWeight{1.};
  event->jetIndex = makeJetCuts(event, 0, &tempWeight);

  switch (step){
  case 0:
    step0EventDump_.precision(3);
    step0EventDump_ << event->jetPF2PATPt[0] << "|" << event->jetPF2PATBDiscriminator[0] << "|";
    break;
  }

  for (unsigned i{0}; i < 4; i++){
    switch (step) {
    case 2:
      step2EventDump_ << ((i < event->jetIndex.size())?event->jetPF2PATPtRaw[event->jetIndex[i]]:-666) << " " << ((i < event->jetIndex.size())?event->jetPF2PATEta[event->jetIndex[i]]:-666) << " ";
      break;
    case 4:
      step4EventDump_ << ((i < event->jetIndex.size())?event->jetPF2PATPtRaw[event->jetIndex[i]]:-666) << " " << ((i < event->jetIndex.size())?event->jetPF2PATEta[event->jetIndex[i]]:-666) << " ";
      break;
    case 6:
      step6EventDump_ << ((i < event->jetIndex.size())?event->jetPF2PATPtRaw[event->jetIndex[i]]:-666) << " " << ((i < event->jetIndex.size())?event->jetPF2PATEta[event->jetIndex[i]]:-666) << " ";
      break;
    }
  }

  switch(step){
  case 0:
    step0EventDump_.precision(3);
    step0EventDump_ << event->metPF2PATPt << "|";
    // Synch Cut Flow stuff
    if ( triggerCuts(event) ) step0EventDump_ << "1"; // Trigger Selection - step 0
    else step0EventDump_ << "0";
    if ( (event->electronIndexTight.size() + event->muonIndexTight.size()) == 3 ) step0EventDump_ << "1"; // 3 tight lepton selection - step 1
    else step0EventDump_ << "0";
    if ( (event->electronIndexTight.size() == getLooseEles(event).size()) && (event->muonIndexTight.size() == getLooseMuons(event).size()) ) step0EventDump_ << "1"; // no additional loose leptons - step 2
    else step0EventDump_ << "0";
    if ( (event->electronIndexTight.size() + event->muonIndexTight.size()) < 3 ) step0EventDump_ << "0"; // Check to ensure there are at least three leptons - otherwise memory leak occurs.
    else if (std::abs(getZCand(event,event->electronIndexTight,event->muonIndexTight)) > invZMassCut_) step0EventDump_ << "1"; // Z selection - step 3
    else step0EventDump_ << "0";
    event->jetIndex = makeJetCuts(event, 0, &tempWeight);
    if ( event->jetIndex.size() >= 1 ) step0EventDump_ << "1"; // Jet selection, at least one jet - step 4
    else step0EventDump_ << "0";
    event->bTagIndex = makeBCuts(event,event->jetIndex);
    if ( event->bTagIndex.size() == 1 ) step0EventDump_ << "1"; // b-Tag selection, exactly one b-jet - step 5
    else step0EventDump_ << "0";
    if ( std::sqrt(2*event->metPF2PATPt*event->wLepton.Pt()*(1-std::cos(event->metPF2PATPhi - event->wLepton.Phi()))) > mTWCutSynch_ ) step0EventDump_ << "1"; // MET selection, step 6
    else step0EventDump_ << "0";
    if ( (getTopMass(event) < TopMassCutUpper_) && (getTopMass(event) > TopMassCutLower_) ) step0EventDump_ << "1"; // Top Mass cut, step 7
    else step0EventDump_ << "0";
    step0EventDump_ << std::endl;
    break;
  case 2:
    step2EventDump_ << event->metPF2PATPt;
    step2EventDump_ << std::endl;
    break;
  case 4:
    step4EventDump_ << event->metPF2PATPt;
    step4EventDump_ << std::endl;
    break;
  case 6:
    step6EventDump_ << event->metPF2PATPt;
    step6EventDump_ << std::endl;
    break;
  }

}

float Cuts::getLeptonWeight(AnalysisEvent * event, int syst){
  //If number of electrons is > 1  then both z pair are electrons, so get their weight
  if (!isMC_) return 1.;
  float leptonWeight{1.};
  if (trileptonChannel_ == true){
    if (numTightEle_ > 1){
      leptonWeight *= eleSF(event->elePF2PATPT[event->zPairIndex.first],event->elePF2PATSCEta[event->zPairIndex.first],syst);
      leptonWeight *= eleSF(event->elePF2PATPT[event->zPairIndex.second],event->elePF2PATSCEta[event->zPairIndex.second],syst);
    }
    else{
      leptonWeight *= muonSF(event->zPairLeptons.first.Pt(),event->zPairLeptons.first.Eta(),syst);
      leptonWeight *= muonSF(event->zPairLeptons.second.Pt(),event->zPairLeptons.second.Eta(),syst);
    }
    if (numTightEle_ == 3 || numTightEle_ == 1){
      leptonWeight *= eleSF(event->elePF2PATPT[event->wLepIndex],event->elePF2PATSCEta[event->wLepIndex],syst);
    }
    else{
      leptonWeight *= eleSF(event->elePF2PATPT[event->wLepIndex],event->elePF2PATSCEta[event->wLepIndex],syst);
    }
  }
  else if(trileptonChannel_ == false){
    if (numTightEle_ == 2){
      leptonWeight *= eleSF(event->elePF2PATPT[event->zPairIndex.first],event->elePF2PATSCEta[event->zPairIndex.first],syst);
      leptonWeight *= eleSF(event->elePF2PATPT[event->zPairIndex.second],event->elePF2PATSCEta[event->zPairIndex.second],syst);
    }
    else if (numTightMu_ == 2){
      leptonWeight *= muonSF(event->zPairLeptons.first.Pt(),event->zPairLeptons.first.Eta(),syst);
      leptonWeight *= muonSF(event->zPairLeptons.second.Pt(),event->zPairLeptons.second.Eta(),syst);
    }
    else if (numTightEle_ == 1 && numTightMu_ == 1){
      leptonWeight *= eleSF(event->elePF2PATPT[event->electronIndexTight[0]],event->elePF2PATSCEta[event->electronIndexTight[0]],syst);
      leptonWeight *= muonSF(event->muonPF2PATPt[event->muonIndexTight[0]],event->muonPF2PATEta[event->muonIndexTight[0]],syst);
    }
  }

  else{
    std::cout << "You've somehow managed to set the trilepton/dilepton bool flag to a value other than these two!" << std::endl;
    std::cout << "HOW?! Well done for breaking this ..." << std::endl;
    exit(0);
  }

  return leptonWeight;

}

float Cuts::get2015TriggerSF(int syst, double eta1, double eta2){

  std::string channel = "";

  //Dilepton channels
  if ( !trileptonChannel_ && cutConfTrigLabel_.find("e") != std::string::npos ) channel = "ee";
  if ( !trileptonChannel_ && cutConfTrigLabel_.find("d") != std::string::npos ) channel = "emu";
  if ( !trileptonChannel_ && cutConfTrigLabel_.find("m") != std::string::npos ) channel = "mumu";
  //Trilepton channels
  if ( trileptonChannel_ && cutConfTrigLabel_.find("e")  != std::string::npos ) channel = "eee";
  if ( trileptonChannel_ && cutConfTrigLabel_.find("d1") != std::string::npos ) channel = "eemu";
  if ( trileptonChannel_ && cutConfTrigLabel_.find("d2") != std::string::npos ) channel = "emumu";
  if ( trileptonChannel_ && cutConfTrigLabel_.find("m")  != std::string::npos ) channel = "mumumu";

  //If using SFs based on eta?
  if ( std::abs(eta1) <= 5.0 && std::abs(eta2) <= 5.0 ) {
	std::cout << "Not implemented yet. Don't pass eta values!" << std::endl;
	exit(999);
  }

  //If no pT or eta dependancy ...
  else {
    //Dilepton channels
    if (channel == "ee"){
      float twgt = 0.954; // tight=0.954; medium=0.958
      if (syst == 1) twgt = 0.963;
      if(syst == 2) twgt = 0.945;
      return twgt;
    }
    if (channel == "mumu"){
      float twgt = 0.927; // tight=0.934; medium=0.931
      if (syst == 1) twgt = 0.932;
      if (syst == 2) twgt = 0.922;
      return twgt;
    }
    if (channel == "emu"){
      float twgt = 0.969;
      if (syst == 1) twgt += 0.006;
      if (syst == 2) twgt -= 0.006;
      return twgt;
    }
    //Trilepton channels
    if (channel == "eee"){
      float twgt = 0.987;
      if (syst == 1) twgt += 0.036;
      if (syst == 2) twgt -= 0.036;
      return twgt;
    }
    if (channel == "eemu"){
      float twgt = 0.987;
      if (syst == 1) twgt += 0.035;
      if (syst == 2) twgt -= 0.035;
      return twgt;
    }
    if (channel == "emumu"){
      float twgt = 0.886;
      if (syst == 1) twgt += 0.042;
      if (syst == 2) twgt -= 0.042;
      return twgt;
    }
    if (channel == "mumumu"){
      float twgt = 0.9871;
      if (syst == 1) twgt += 0.0242;
      if (syst == 2) twgt -= 0.0212;
      return twgt;
    }
  }

 std::cout << "Trigger not found!" << std::endl;
 return 0.0; // Return 0.0 if trigger isn't found.
}

float Cuts::get2016TriggerSF(int syst, double eta1, double eta2){

  std::string channel = "";

  //Dilepton channels
  if ( !trileptonChannel_ && cutConfTrigLabel_.find("e") != std::string::npos ) channel = "ee";
  if ( !trileptonChannel_ && cutConfTrigLabel_.find("d") != std::string::npos ) channel = "emu";
  if ( !trileptonChannel_ && cutConfTrigLabel_.find("m") != std::string::npos ) channel = "mumu";
  //Trilepton channels
  if ( trileptonChannel_ && cutConfTrigLabel_.find("e")  != std::string::npos ) channel = "eee";
  if ( trileptonChannel_ && cutConfTrigLabel_.find("d1") != std::string::npos ) channel = "eemu";
  if ( trileptonChannel_ && cutConfTrigLabel_.find("d2") != std::string::npos ) channel = "emumu";
  if ( trileptonChannel_ && cutConfTrigLabel_.find("m")  != std::string::npos ) channel = "mumumu";

  //If using SFs based on eta?
  if ( std::abs(eta1) <= 5.0 && std::abs(eta2) <= 5.0 ) {
	std::cout << "Not implemented yet. Don't pass eta values!" << std::endl;
	exit(999);
  }

  //If no pT or eta dependancy ...
  else {
    //Dilepton channels
    if (channel == "ee"){
      float twgt = 0.931; // 0.931 for eff; 0.968 for SF
      if (syst == 1) twgt += 0.003; // 0.003 for eff; 0.002 for SF
      if (syst == 2) twgt -= 0.003;
      return twgt;
    }
    if (channel == "mumu"){
      float twgt = 0.772; // 0.772 for eff; 0.853 for SF
      if (syst == 1) twgt += 0.003; // 0.003 for eff; 0.002 for SF
      if (syst == 2) twgt -= 0.003;
      return twgt;
    }
    if (channel == "emu"){
      float twgt = 0.895; // 0.895 for eff; 0.986 for SF
      if (syst == 1) twgt += 0.003; // 0.003 for eff; 0.002 for SF
      if (syst == 2) twgt -= 0.003;
      return twgt;
    }
    //Trilepton channels
    if (channel == "eee"){
      float twgt = 0.894;
      if (syst == 1) twgt += 0.003;
      if (syst == 2) twgt -= 0.003;
      return twgt;
    }
    if (channel == "eemu"){
      float twgt = 0.987;
      if (syst == 1) twgt += 0.035;
      if (syst == 2) twgt -= 0.035;
      return twgt;
    }
    if (channel == "emumu"){
      float twgt = 0.886;
      if (syst == 1) twgt += 0.042;
      if (syst == 2) twgt -= 0.042;
      return twgt;
    }
    if (channel == "mumumu"){
      float twgt = 0.9871;
      if (syst == 1) twgt += 0.0242;
      if (syst == 2) twgt -= 0.0212;
      return twgt;
    }
  }

 std::cout << "Trigger not found!" << std::endl;
 return 0.0; // Return 0.0 if trigger isn't found.
}

float Cuts::eleSF(double pt, double eta, int syst){

  double maxPt{h_eleSFs->GetYaxis()->GetXmax()};
  unsigned bin1{0}, bin2{0};

  // If cut-based, std::abs eta, else just eta
  if ( pt <= maxPt ){
    bin1 = h_eleSFs->FindBin(std::abs(eta),pt);
    bin2 = h_eleReco->FindBin(eta,pt);
  }
  else {
    bin1 = h_eleSFs->FindBin(std::abs(eta),maxPt);
    bin2 = h_eleReco->FindBin(eta,maxPt);
  }

  float eleIdSF   = h_eleSFs->GetBinContent(bin1);
  float eleRecoSF = h_eleReco->GetBinContent(bin2);

  if ( syst == 1 ) {
    eleIdSF   += h_eleSFs->GetBinError(bin1);
    eleRecoSF += h_eleReco->GetBinError(bin2);
  }

  if ( syst == 2 ) {
    eleIdSF   -= h_eleSFs->GetBinError(bin1);
    eleRecoSF -= h_eleReco->GetBinError(bin2);
  }

  return eleIdSF*eleRecoSF;
}

float Cuts::muonSF(double pt, double eta, int syst){

  double maxIdPt{h_muonIDs->GetYaxis()->GetXmax()};
  double maxIsoPt{h_muonPFiso->GetYaxis()->GetXmax()};
  unsigned binId {0};
  unsigned binIso {0};

  if ( pt <= maxIdPt ) binId = h_muonIDs->FindBin(std::abs(eta),pt);
  else binId = h_muonIDs->FindBin(std::abs(eta),maxIdPt);

  if ( pt <= maxIsoPt ) binIso = h_muonPFiso->FindBin(std::abs(eta),pt);
  else binIso = h_muonPFiso->FindBin(std::abs(eta),maxIsoPt);

  float muonIdSF    = h_muonIDs->GetBinContent(binId);
  float muonPFisoSF = h_muonPFiso->GetBinContent(binIso);
  float muonRecoSF = 1.0;

  if ( is2016_ ) {
    muonRecoSF = h_muonRecoGraph->Eval(eta);
  }

  if ( syst == 1 ) {
    muonIdSF    += h_muonIDs->GetBinError(binId);
    muonPFisoSF += h_muonPFiso->GetBinError(binIso);
    if ( is2016_ ) muonRecoSF += muonRecoSF*0.01;
  }

  if ( syst == 2 ) {
    muonIdSF    -= h_muonIDs->GetBinError(binId);
    muonPFisoSF -= h_muonPFiso->GetBinError(binIso);
    if ( is2016_ ) muonRecoSF -= muonRecoSF*0.01;
  }

  return muonIdSF*muonPFisoSF*muonRecoSF;
}

void Cuts::initialiseJECCors(){
  std::ifstream jecFile;
  if ( !is2016_ ) jecFile.open( "scaleFactors/2015/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt", std::ifstream::in );
  else jecFile.open( "scaleFactors/2016/Spring16_23Sep2016V2_MC_Uncertainty_AK4PFchs.txt", std::ifstream::in );
  std::string line;
  bool first{true};

  if (!jecFile.is_open()){
    std::cout << "Unable to open jecFile." << std::endl;
    exit(0);
  }

  while(getline(jecFile,line)){
    std::vector<std::string> tempVec;
    std::stringstream lineStream{line};
    std::string item;
    while (std::getline(lineStream,item,' ')){
      tempVec.emplace_back(item);
    }
    std::vector<float> tempUp;
    std::vector<float> tempDown;

    etaMinJEC_.emplace_back(std::stof(tempVec[0]));
    etaMaxJEC_.emplace_back(std::stof(tempVec[1]));
    for (unsigned i{1}; i < tempVec.size()/3; i++){
      unsigned ind{i * 3};
      if (first){
      	ptMinJEC_.emplace_back(std::stof(tempVec[ind]));
	ptMaxJEC_.emplace_back((ind+3 >= tempVec.size()?10000.:std::stof(tempVec[ind+3])));
      }
      tempUp.emplace_back(std::stof(tempVec[ind+1]));
      tempDown.emplace_back(std::stof(tempVec[ind+2]));
    }
    jecSFUp_.emplace_back(tempUp);
    jecSFDown_.emplace_back(tempDown);
    first = false;
  }
}

float Cuts::getJECUncertainty(float pt, float eta, int syst){
  if (!(syst == 4 || syst == 8)){
    return 0.;
  }
  unsigned ptBin{0};
  unsigned etaBin{0};
  for (unsigned i{0}; i < ptMinJEC_.size(); i++){
    if (pt > ptMinJEC_[i] && pt < ptMaxJEC_[i]){
      ptBin = i;
      break;
    }
  }
  for (unsigned i{0}; i < etaMinJEC_.size(); i++){
    if (eta > etaMinJEC_[i] && eta < etaMaxJEC_[i]){
      etaBin = i;
      break;
    }
  }

  float lowFact{syst==4?jecSFUp_[etaBin][ptBin]:jecSFDown_[etaBin][ptBin]};
  float hiFact{syst==4?jecSFUp_[etaBin][ptBin+1]:jecSFDown_[etaBin][ptBin+1]};
  //Now do some interpolation
  float a{(hiFact - lowFact)/(ptMaxJEC_[ptBin]-ptMinJEC_[ptBin])};
  float b{(lowFact * (ptMaxJEC_[ptBin]) - hiFact*ptMinJEC_[ptBin])/(ptMaxJEC_[ptBin] - ptMinJEC_[ptBin])};
  return (syst == 4 ? a * pt + b : -(a * pt + b));

}

TLorentzVector Cuts::getJetLVec(AnalysisEvent* event, int index, int syst){

  TLorentzVector returnJet;
  float newSmearValue{1.0};

  if ( synchCutFlow_ && trileptonChannel_ ) {
	returnJet.SetPxPyPzE(event->jetPF2PATPx[index],event->jetPF2PATPy[index],event->jetPF2PATPz[index],event->jetPF2PATE[index]);
	return returnJet;
  }

  if ( !isMC_ ) {
        event->jetSmearValue.emplace_back(newSmearValue);
	returnJet.SetPxPyPzE(event->jetPF2PATPx[index],event->jetPF2PATPy[index],event->jetPF2PATPz[index],event->jetPF2PATE[index]);
	return returnJet;
  }

  float jerSF{1.};
  float jerSigma{0.};
  std::pair< float, float > jetSFs{};

  if ( !is2016_ ) jetSFs = jet2015SFs( std::abs(event->jetPF2PATEta[index]) );
  else jetSFs = jet2016SFs( std::abs(event->jetPF2PATEta[index]) );

  jerSF = jetSFs.first;
  jerSigma = jetSFs.second;

  if ( isMC_ && event->genJetPF2PATPT[index] > -990.){
    if ( deltaR(event->genJetPF2PATEta[index],event->genJetPF2PATPhi[index],event->jetPF2PATEta[index],event->jetPF2PATPhi[index]) < 0.4/2.0 /*&& std::abs(event->jetPF2PATPtRaw[index] - event->genJetPF2PATPT[index]) < 3.0*jerSigma*/ ) { // If matching from GEN to RECO using dR<Rcone/2 and dPt < 3*sigma, just scale, just scale
      if (syst == 16) jerSF += jerSigma;
      else if (syst == 32) jerSF -= jerSigma;
      newSmearValue = std::max(0.0, event->jetPF2PATPtRaw[index] + (event->jetPF2PATPtRaw[index] - event->genJetPF2PATPT[index]) * jerSF)/event->jetPF2PATPtRaw[index];
      returnJet.SetPxPyPzE(newSmearValue*event->jetPF2PATPx[index],newSmearValue*event->jetPF2PATPy[index],newSmearValue*event->jetPF2PATPz[index],newSmearValue*event->jetPF2PATE[index]);
    }
      else { // If not, randomly smear
      newSmearValue = 1.0+TRandom(rand()).Gaus(0.0, std::sqrt(jerSF*jerSF-1)*jerSigma);
      returnJet.SetPxPyPzE(newSmearValue*event->jetPF2PATPx[index],newSmearValue*event->jetPF2PATPy[index],newSmearValue*event->jetPF2PATPz[index],newSmearValue*event->jetPF2PATE[index]);
    }

  }

  else returnJet.SetPxPyPzE(event->jetPF2PATPx[index],event->jetPF2PATPy[index],event->jetPF2PATPz[index],event->jetPF2PATE[index]);

  if (isMC_){
    float jerUncer{getJECUncertainty(returnJet.Pt(),returnJet.Eta(),syst)};
    returnJet *= 1+jerUncer;
  }
  tempSmearValue_ = newSmearValue;
  return returnJet;
}

std::pair< float, float > Cuts::jet2015SFs( float eta ) {
  // JER Scaling Factors and uncertainities for 2015
  float jerSF{0.};
  float jerSigma{0.};

  if (eta <= 0.5) {
    jerSF = 1.095;
    jerSigma = 0.018;
  }
  else if (eta <= 0.8) {
    jerSF = 1.120;
    jerSigma = 0.028;
  }
  else if (eta <= 1.1) {
    jerSF = 1.097;
    jerSigma = 0.017;
  }
  else if (eta <= 1.3) {
    jerSF = 1.103;
    jerSigma = 0.033;
  }
  else if (eta <= 1.7) {
    jerSF = 1.118;
    jerSigma = 0.014;
  }
  else if (eta <= 1.9) {
    jerSF = 1.100;
    jerSigma = 0.033;
  }
  else if (eta <= 2.1) {
    jerSF = 1.162;
    jerSigma = 0.044;
  }
  else if (eta <= 2.3) {
    jerSF = 1.160;
    jerSigma = 0.048;
  }
  else if (eta <= 2.5) {
    jerSF = 1.161;
    jerSigma = 0.060;
  }
  else if (eta <= 2.8) {
    jerSF = 1.209;
    jerSigma = 0.059;
  }
  else if (eta <= 3.0){
    jerSF = 1.564;
    jerSigma = 0.321;
  }
  else if (eta <= 3.2){
    jerSF = 1.384;
    jerSigma = 0.033;
  }
  else {
    jerSF = 1.216;
    jerSigma = 0.050;
  }

  return std::make_pair( jerSF, jerSigma );
}

std::pair< float, float > Cuts::jet2016SFs( float eta ) {
  // JER Scaling Factors and uncertainities for 2016
  float jerSF{0.};
  float jerSigma{0.};

  if (eta <= 0.5) {
    jerSF = 1.109;
    jerSigma = 0.008;
  }
  else if (eta <= 0.8) {
    jerSF = 1.138;
    jerSigma = 0.013;
  }
  else if (eta <= 1.1) {
    jerSF = 1.114;
    jerSigma = 0.013;
  }
  else if (eta <= 1.3) {
    jerSF = 1.123;
    jerSigma = 0.024;
  }
  else if (eta <= 1.7) {
    jerSF = 1.084;
    jerSigma = 0.011;
  }
  else if (eta <= 1.9) {
    jerSF = 1.082;
    jerSigma = 0.035;
  }
  else if (eta <= 2.1) {
    jerSF = 1.140;
    jerSigma = 0.047;
  }
  else if (eta <= 2.3) {
    jerSF = 1.067;
    jerSigma = 0.053;
  }
  else if (eta <= 2.5) {
    jerSF = 1.177;
    jerSigma = 0.041;
  }
  else if (eta <= 2.8) {
    jerSF = 1.364;
    jerSigma = 0.039;
  }
  else if (eta <= 3.0){
    jerSF = 1.857;
    jerSigma = 0.071;
  }
  else if (eta <= 3.2){
    jerSF = 1.328;
    jerSigma = 0.022;
  }
  else {
    jerSF = 1.160;
    jerSigma = 0.029;
  }

  return std::make_pair( jerSF, jerSigma );
}

void Cuts::getBWeight(AnalysisEvent* event, TLorentzVector jet, int index, float * mcTag, float * mcNoTag, float * dataTag, float * dataNoTag, float * err1, float * err2, float * err3, float * err4){
  //Use b-tagging efficiencies and scale factors.
  //Firstly get efficiency for pt/eta bin here.
  float eff{1.};
  int partonFlavour{std::abs(event->jetPF2PATPID[index])};

  if (partonFlavour == 0) return;
  if (partonFlavour == 5){
    eff = bTagEffPlots_[4]->GetBinContent(bTagEffPlots_[4]->GetXaxis()->FindBin(jet.Pt()),bTagEffPlots_[4]->GetYaxis()->FindBin(std::abs(jet.Eta()))) / bTagEffPlots_[0]->GetBinContent(bTagEffPlots_[0]->GetXaxis()->FindBin(jet.Pt()),bTagEffPlots_[0]->GetYaxis()->FindBin(std::abs(jet.Eta())));
  }
  if (partonFlavour == 4){
    eff = bTagEffPlots_[5]->GetBinContent(bTagEffPlots_[5]->GetXaxis()->FindBin(jet.Pt()),bTagEffPlots_[5]->GetYaxis()->FindBin(std::abs(jet.Eta()))) / bTagEffPlots_[1]->GetBinContent(bTagEffPlots_[1]->GetXaxis()->FindBin(jet.Pt()),bTagEffPlots_[1]->GetYaxis()->FindBin(std::abs(jet.Eta())));
  }
  if (partonFlavour < 4){
    eff = bTagEffPlots_[6]->GetBinContent(bTagEffPlots_[6]->GetXaxis()->FindBin(jet.Pt()),bTagEffPlots_[6]->GetYaxis()->FindBin(std::abs(jet.Eta()))) / bTagEffPlots_[2]->GetBinContent(bTagEffPlots_[2]->GetXaxis()->FindBin(jet.Pt()),bTagEffPlots_[2]->GetYaxis()->FindBin(std::abs(jet.Eta())));
  }
  if (partonFlavour == 21){
    eff = bTagEffPlots_[7]->GetBinContent(bTagEffPlots_[7]->GetXaxis()->FindBin(jet.Pt()),bTagEffPlots_[7]->GetYaxis()->FindBin(std::abs(jet.Eta()))) / bTagEffPlots_[3]->GetBinContent(bTagEffPlots_[3]->GetXaxis()->FindBin(jet.Pt()),bTagEffPlots_[3]->GetYaxis()->FindBin(std::abs(jet.Eta())));
  }

  //Get SF
  // Initalise variables.
  double jet_scalefactor{1.};
  double jet_scalefactor_up{1.};
  double jet_scalefactor_do{1.};

  float SFerr{0.};
  float jetPt = jet.Pt();
  float maxBjetPt{670};
  float maxLjetPt{1000.0};
  bool doubleUncertainty{false};
  //Do some things if it's a b or c

  if ( partonFlavour == 5 ){
    if (jetPt > maxBjetPt){
      jetPt = maxBjetPt;
      doubleUncertainty = true;
    }
    //jet_scalefactor = beautyReader.eval_auto_bounds("central", BTagEntry::FLAV_B, jet.Eta(), jetPt);
    //jet_scalefactor_up = beautyReader.eval_auto_bounds("up", BTagEntry::FLAV_B, jet.Eta(), jetPt);
    //jet_scalefactor_do = beautyReader.eval_auto_bounds("down", BTagEntry::FLAV_B, jet.Eta(), jetPt);
    jet_scalefactor = getBweight_backup( 0, 0, jetPt );
    jet_scalefactor_up = getBweight_backup( 0, 1, jetPt );
    jet_scalefactor_do = getBweight_backup( 0, -1, jetPt );
  }

  else if ( partonFlavour == 4 ){
    if (jetPt > maxBjetPt){
      jetPt = maxBjetPt;
      doubleUncertainty = true;
    }
    //jet_scalefactor = charmReader.eval_auto_bounds("central", BTagEntry::FLAV_C, jet.Eta(), jetPt);
    //jet_scalefactor_up = charmReader.eval_auto_bounds("up", BTagEntry::FLAV_C, jet.Eta(), jetPt);
    //jet_scalefactor_do = charmReader.eval_auto_bounds("down", BTagEntry::FLAV_C, jet.Eta(), jetPt);
    jet_scalefactor = getBweight_backup( 1, 0, jetPt );
    jet_scalefactor_up = getBweight_backup( 1, 1, jetPt );
    jet_scalefactor_do = getBweight_backup( 1, -1, jetPt );
  }

  //Light jets
  else {
    if (jetPt > maxLjetPt){
      jetPt = maxLjetPt;
      doubleUncertainty = true;
    }
    //jet_scalefactor = lightReader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, jet.Eta(), jetPt);
    //jet_scalefactor_up = lightReader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG, jet.Eta(), jetPt);
    //jet_scalefactor_do = lightReader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG, jet.Eta(), jetPt);
    jet_scalefactor = getBweight_backup( 2, 0, jetPt );
    jet_scalefactor_up = getBweight_backup( 2, 1, jetPt );
    jet_scalefactor_do = getBweight_backup( 2, -1, jetPt );
  }

  if (doubleUncertainty) {
    jet_scalefactor_up = 2*(jet_scalefactor_up - jet_scalefactor) + jet_scalefactor;
    jet_scalefactor_do = 2*(jet_scalefactor_do - jet_scalefactor) + jet_scalefactor;
  }


  SFerr = std::abs(jet_scalefactor_up - jet_scalefactor)>std::abs(jet_scalefactor_do - jet_scalefactor)? std::abs(jet_scalefactor_up - jet_scalefactor):std::abs(jet_scalefactor_do - jet_scalefactor);

  //Apply the weight of the jet and set the error
  if (event->jetPF2PATBDiscriminator[index] > bDiscCut_){
    *mcTag *= eff;
    *dataTag *= eff*jet_scalefactor;

    if (partonFlavour == 5 || partonFlavour == 4) *err1 += SFerr/jet_scalefactor;
    else *err3 += SFerr/jet_scalefactor;
  }
  else{
    *mcNoTag *= (1-eff);
    *dataNoTag *= (1-eff*jet_scalefactor);

    if (partonFlavour == 5 || partonFlavour == 4) *err2 += (-eff*SFerr)/(1-eff*jet_scalefactor);
    else *err4 += (-eff*SFerr)/(1-eff*jet_scalefactor);

  }

}

// Backup temporary method to do Btag Scale Factors whilst debugging is ongoing.

float Cuts::getBweight_backup(int flavour, int type, float pt){

  float sf = 1.0;
  float x = pt;

  if (!is2016_){

    if (flavour == 0) { // B flavour
      if (type == 0) sf = 0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x)));

      if (type == 1) { // scale up
	if ( pt < 50.0 )   sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.019803794100880623;
	if ( pt < 70.0 )   sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.026958625763654709;
	if ( pt < 100.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.024285079911351204;
	if ( pt < 140.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.028512096032500267;
	if ( pt < 200.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.029808893799781799;
	if ( pt < 300.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.026503190398216248;
	else sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.042264193296432495;
      }

      if (type == -1) { // scale down
	if ( pt < 50.0 )   sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.019803794100880623;
	if ( pt < 70.0 )   sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.026958625763654709;
	if ( pt < 100.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.024285079911351204;
	if ( pt < 140.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.028512096032500267;
	if ( pt < 200.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.029808893799781799;
	if ( pt < 300.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.026503190398216248;
	else sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.042264193296432495;
      }
    }

    if (flavour == 1) { // C flavour
      if (type == 0) sf = 0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x)));

      if ( type == 1 ) {
	if ( pt < 50.0 )   sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.039607588201761246;
	if ( pt < 70.0 )   sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.053917251527309418;
	if ( pt < 100.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.048570159822702408;
	if ( pt < 140.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.057024192065000534;
	if ( pt < 200.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.059617787599563599;
	if ( pt < 300.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.053006380796432495;
	else sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))+0.08452838659286499;
     }

      if ( type == -1 ) {
	if ( pt < 50.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.039607588201761246;
	if ( pt < 70.0 )  sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.053917251527309418;
	if ( pt < 100.0 ) sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.048570159822702408;
	if ( pt < 140.0 ) sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.057024192065000534;
	if ( pt < 200.0 ) sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.057024192065000534;
	if ( pt < 300.0 ) sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.059617787599563599;
	else sf = (0.886376*((1.+(0.00250226*x))/(1.+(0.00193725*x))))-0.08452838659286499;
       }
    }
    if (flavour == 2) { // UDSG flavour
      if (type == 0)  sf = 0.992339;
      if (type == 1)  sf = 1.17457;
      if (type == -1) sf = 0.810103;
    }
  }
  else {
    if (flavour == 0) { // B flavour
      if ( type == 0 ) sf = 0.857294+(3.75846e-05*x);
      if ( type == 1 ) {
        if ( pt < 50.0 )  sf = (0.857294+(3.75846e-05*x))+0.01438605971634388;
        if ( pt < 70.0 )  sf = (0.857294+(3.75846e-05*x))+0.013547157868742943;
        if ( pt < 100.0 ) sf = (0.857294+(3.75846e-05*x))+0.01491188071668148;
        if ( pt < 140.0 ) sf = (0.857294+(3.75846e-05*x))+0.019157662987709045;
        if ( pt < 200.0 ) sf = (0.857294+(3.75846e-05*x))+0.021716726943850517;
        if ( pt < 300.0 ) sf = (0.857294+(3.75846e-05*x))+0.046845909208059311;
        else sf = (0.857294+(3.75846e-05*x))+0.042299006134271622;
      }
      if ( type == -1 ) {
        if ( pt < 50.0 )  sf = 0.857294+((3.75846e-05*x)-0.028772119432687759);
        if ( pt < 70.0 )  sf = 0.857294+((3.75846e-05*x)-0.027094315737485886);
        if ( pt < 100.0 ) sf = 0.857294+((3.75846e-05*x)-0.029823761433362961);
        if ( pt < 140.0 ) sf = 0.857294+((3.75846e-05*x)-0.038315325975418091);
        if ( pt < 200.0 ) sf = 0.857294+((3.75846e-05*x)-0.043433453887701035);
        if ( pt < 300.0 ) sf = 0.857294+((3.75846e-05*x)-0.093691818416118622);
        else sf = 0.857294+((3.75846e-05*x)-0.084598012268543243);
      }
    }
    if (flavour == 1) { // C flavour
      if ( type == 0 ) sf = 0.857294+(3.75846e-05*x);
      if ( type == 1 ) {
        if ( pt < 50.0 )  sf = (0.857294+(3.75846e-05*x))+0.028772119432687759;
        if ( pt < 70.0 )  sf = (0.857294+(3.75846e-05*x))+0.027094315737485886;
        if ( pt < 100.0 ) sf = (0.857294+(3.75846e-05*x))+0.029823761433362961;
        if ( pt < 140.0 ) sf = (0.857294+(3.75846e-05*x))+0.038315325975418091;
        if ( pt < 200.0 ) sf = (0.857294+(3.75846e-05*x))+0.043433453887701035;
        if ( pt < 300.0 ) sf = (0.857294+(3.75846e-05*x))+0.093691818416118622;
        else sf = (0.857294+(3.75846e-05*x))+0.084598012268543243;
      }
      if ( type == -1 ) {
        if ( pt < 50.0 )  sf = 0.857294+((3.75846e-05*x)-0.01438605971634388);
        if ( pt < 70.0 )  sf = 0.857294+((3.75846e-05*x)-0.013547157868742943);
        if ( pt < 100.0 ) sf = 0.857294+((3.75846e-05*x)-0.01491188071668148);
        if ( pt < 140.0 ) sf = 0.857294+((3.75846e-05*x)-0.019157662987709045);
        if ( pt < 200.0 ) sf = 0.857294+((3.75846e-05*x)-0.021716726943850517);
        if ( pt < 300.0 ) sf = 0.857294+((3.75846e-05*x)-0.046845909208059311);
        else sf = 0.857294+((3.75846e-05*x)-0.042299006134271622);
      }
    }
    if (flavour == 2) { // UDSG flavour
      if (type == 0)  sf = 0.688619+260.84/(x*x);
      if (type == 1)  sf = (0.688619+260.84/(x*x))*(1+(0.144982+0.000116685*x+-1.0021e-07*x*x));
      if (type == -1) sf = (0.688619+260.84/(x*x))*(1-(0.144982+0.000116685*x+-1.0021e-07*x*x));
    }
  }
  return sf;
}
