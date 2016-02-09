#include "cutClass.hpp"
#include <sstream>
#include <cmath>
#include <math.h>
#include "TLorentzVector.h"
#include <iostream>
#include <libconfig.h++>
#include <iomanip>
#include <fstream>

Cuts::Cuts(bool doPlots, bool fillCutFlows,bool invertIsoCut, bool lepCutFlow, bool dumpEventNumber, const bool trileptonChannel):

  //Do plots?
  doPlots_(doPlots),
  fillCutFlow_(fillCutFlows),
  //background estimation. May not be possible
  invertIsoCut_(invertIsoCut),
  //Synchronisation cut flow.
  synchCutFlow_(lepCutFlow),
  //Synchronisation cut flow.
  singleEventInfoDump_(dumpEventNumber),
  trileptonChannel_(trileptonChannel),

  // Set all default parameters. These will be editable later on, probably.
  numTightEle_(3),
  tightElePt_(20.),
  tightEleEta_(2.5),
  tightEled0_(0.011811),
  tightEleMissLayers_(0),
  tightEleCheckPhotonVeto_(true),
  tightEleMVA_(0.5),
  tightEleRelIso_(0.15),
  //Loose electron initialisation
  numLooseEle_(3),
  looseElePt_(10),
  looseEleEta_(2.5),
  looseEleMVA_(0.5),
  looseEleRelIso_(0.15),
  //Tight muon initialisation
  numTightMu_(0),
  tightMuonPt_(20.),
  tightMuonEta_(2.4),
  tightMuonRelIso_(0.15),
  //Loose muons
  numLooseMu_(0),
  looseMuonPt_(10.),
  looseMuonEta_(2.4),
  looseMuonRelIso_(0.25),
  //zMass cuts
  invZMassCut_(15.),
  invWMassCut_(5.),
  //Jet initialisation
  numJets_(2),
  maxJets_(4),
  jetPt_(30.),
  jetEta_(2.5),
  jetNConsts_(2),
  jetIDDo_(true),
  //B-discriminator cut
  numbJets_(1),
  maxbJets_(2),
  //bDiscCut_(0.935), // Tight cut
  bDiscCut_(0.80), // Medium level
  //bDiscCut_(0.460), // Loose cut
  //Set isMC. Default is true, but it's called everytime a new dataset is processed anyway.
  isMC_(true),
  //Same for trigger flag.
  triggerFlag_(""),
  //Make cloned tree false for now
  postLepSelTree_(0),
  postLepSelTree2_(0),
  postLepSelTree3_(0),
  //Skips running trigger stuff
  skipTrigger_(false),
  //Are we making b-tag efficiency plots?
  makeBTagEffPlots_(false),
  getBTagWeight_(false),
  //MET and mTW cuts go here.
  metCut_(0.),
  mTWCut_(0.),
  TopMassCutLower_(91.),
  TopMassCutUpper_(155.)
{
  //Space here in case other stuff needs to be done.

  //If doing synchronisation., initialise that here.
  if (synchCutFlow_){
    synchCutFlowHist_ = new TH1F("synchCutFlow","synchCutFlow",9,0,9);
    synchNumEles_ = new TH1I("synchNumEles","synchNumEles",10,0,10);
    synchNumMus_ = new TH1I("synchNumMuos","synchNumMuos",10,0,10);
    synchMuonCutFlow_ = new TH1I("synchMuonCutFlow","synchMuonCutFlow",11,0,11);
    synchCutTopMassHist_ = new TH1F("synchCutTopMassHist", "synchCutTopMassHist", 200, 0., 200.);
  }
  std::cout << "Initialises fine" << std::endl;
  initialiseJECCors();
  std::cout << "Gets past JEC Cors" << std::endl;
  
  //Initialise array of b-tagging errors here.
  SFb_error = { 0.033408, 0.015446, 0.0146992, 0.0183964, 0.0185363, 0.0145547, 0.0176743, 0.0203609, 0.0143342, 0.0148771, 0.0157936, 0.0176496, 0.0209156, 0.0278529, 0.0346877, 0.0350101 };
}

Cuts::~Cuts(){
  if (synchCutFlow_) {
    delete synchCutFlowHist_;
    delete synchNumEles_;
    delete synchNumMus_;
    delete synchMuonCutFlow_;
    delete synchCutTopMassHist_;
    if (makeEventDump_){ step0EventDump_.close();
      step2EventDump_.close();
      step4EventDump_.close();
      step6EventDump_.close();
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
    eles.lookupValue("MVA",tightEleMVA_);
    eles.lookupValue("number",numTightEle_);
  }
  if (cuts.exists("looseElectrons")){
    libconfig::Setting& eles = cuts["looseElectrons"];
    eles.lookupValue("pt",looseElePt_);
    eles.lookupValue("eta",looseEleEta_);
    eles.lookupValue("relIso",looseEleRelIso_);
    eles.lookupValue("MVA",looseEleMVA_);
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
  std::cerr << "And so it's looking for " << numTightMu_ << " muons and " << numTightEle_ << " electrons" << std::endl;
  
  if (makeEventDump_){ 
    step0EventDump_.open("step0EventDump"+postfixName_+".txt");
    step2EventDump_.open("step2EventDump"+postfixName_+".txt");
    step4EventDump_.open("step4EventDump"+postfixName_+".txt");
    step6EventDump_.open("step6EventDump"+postfixName_+".txt");
  }

  
  return true;
}

bool Cuts::makeCuts(AnalysisEvent *event, float *eventWeight, std::map<std::string,Plots*> plotMap, TH1F* cutFlow, int systToRun){

  //If we're doing synchronisation, do this function.
  if (synchCutFlow_){
    return synchCuts(event);
  }
  
  if (!isMC_) if (!triggerCuts(event)) return false;
  //Make lepton cuts. Does the inverted iso cuts if necessary.
  if (!(invertIsoCut_?invertIsoCut(event,eventWeight, plotMap,cutFlow):makeLeptonCuts(event,eventWeight, plotMap,cutFlow))) return false;
  //  if (!makeLeptonCuts(event,eventWeight,plotMap,cutFlow)) return false;
  //This is to make some skims for faster running. Do lepSel and save some files.
  if(postLepSelTree_) {
    if (postLepSelTree_->GetEntriesFast() < 40000) postLepSelTree_->Fill();
    else {
      if (postLepSelTree2_->GetEntriesFast() < 40000) postLepSelTree2_->Fill();
      else postLepSelTree3_->Fill();

    }
  }

  event->jetIndex = makeJetCuts(event, systToRun, eventWeight);
  if (doPlots_) plotMap["zMass"]->fillAllPlots(event,*eventWeight);
  if (event->jetIndex.size() < numJets_) return false;
  if (event->jetIndex.size() > maxJets_) return false;

  if (doPlots_||fillCutFlow_) cutFlow->Fill(2.5,*eventWeight);
  event->bTagIndex = makeBCuts(event,event->jetIndex);
  if (doPlots_) plotMap["jetSel"]->fillAllPlots(event,*eventWeight);
  if (event->bTagIndex.size() < numbJets_) return false;
  if (event->bTagIndex.size() > maxbJets_) return false;

  if (doPlots_) plotMap["bTag"]->fillAllPlots(event,*eventWeight);
  if (doPlots_||fillCutFlow_) cutFlow->Fill(3.5,*eventWeight);
  if ( !trileptonChannel_ && getWbosonQuarksCand(event,event->jetIndex) > invWMassCut_ ) return false;
  if ( doPlots_ && !trileptonChannel_ ) plotMap["wMass"]->fillAllPlots(event,*eventWeight);
  if ( doPlots_ && !trileptonChannel_ ) cutFlow->Fill(4.5,*eventWeight);

  //Apply met and mtw cuts here. By default these are 0, so don't do anything.
  if (trileptonChannel_ && event->metPF2PATPt < metCut_) return false;
  TLorentzVector tempMet;
  tempMet.SetPtEtaPhiE(event->metPF2PATPt,0,event->metPF2PATPhi,event->metPF2PATEt);
  float mtw = std::sqrt(2*event->metPF2PATPt*event->wLepton.Pt()*(1-cos(event->metPF2PATPhi - event->wLepton.Phi())));
  if (trileptonChannel_ && mtw < mTWCut_) return false;

  return true;
}

//Make lepton cuts. Will become customisable in a config later on.
bool Cuts::makeLeptonCuts(AnalysisEvent* event,float * eventWeight,std::map<std::string,Plots*> plotMap, TH1F* cutFlow){

  //Do lepton selection. 
  event->electronIndexTight = getTightEles(event);
  if (event->electronIndexTight.size() != numTightEle_) return false;
  event->electronIndexLoose = getLooseEles(event);
  if (event->electronIndexLoose.size() != numLooseEle_) return false;
  event->muonIndexTight = getTightMuons(event);
  if (event->muonIndexTight.size() != numTightMu_) return false;
  event->muonIndexLoose = getLooseMuons(event);
  if (event->muonIndexLoose.size() != numLooseMu_) return false;


  //Should I make it return which leptons are the zMass candidate? Probably.
  float invZmass (9999.);

  if (trileptonChannel_ == true){
    invZmass = getZCand(event, event->electronIndexTight, event->muonIndexTight);
  }
  else if (trileptonChannel_ == false){
    invZmass = getDileptonZCand(event, event->electronIndexTight, event->muonIndexTight);
  }
  else{
    std::cout << "You've somehow managed to set the trilepton/dilepton bool flag to a value other than these two!" << std::endl;
    std::cout << "HOW?! Well done for breaking this ..." << std::endl;
    exit(0);
  }

  * eventWeight *= getLeptonWeight(event);
  if(doPlots_) plotMap["lepSel"]->fillAllPlots(event,*eventWeight);
  if(doPlots_||fillCutFlow_){
    cutFlow->Fill(0.5,*eventWeight);
  }

  if (fabs(invZmass) > invZMassCut_) return false;

  //  plotMap["zMass"]->fillAllPlots(event,eventWeight);
  if(doPlots_||fillCutFlow_) cutFlow->Fill(1.5,*eventWeight);
  return true;
}

std::vector<int> Cuts::getTightEles(AnalysisEvent* event) {
  std::vector<int> electrons;
  for (int i = 0; i < event->numElePF2PAT; i++){
    if (!event->elePF2PATIsGsf[i]) continue;
    TLorentzVector tempVec(event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]);
    /*    std::cout << tempVec.Pt();
    std::cout << " | " << event->elePF2PATComRelIsoRho[i]/tempVec.Pt();
    std::cout << " | " << event->elePF2PATRhoIso[i];
    std::cout << " | " << event->elePF2PATAEff03[i];
    std::cout << " | " << event->elePF2PATChHadIso[i];
    std::cout << " | " << event->elePF2PATNtHadIso[i];
    std::cout << " | " << event->elePF2PATGammaIso[i] << std::endl;*/

    if (tempVec.Pt() < tightElePt_) continue;
    if (std::abs(tempVec.Eta()) > tightEleEta_)continue;
    if (std::abs(event->elePF2PATD0PV[i]) > tightEled0_)continue;
    if (event->elePF2PATMissingInnerLayers[i] > tightEleMissLayers_) continue;
    if (!event->elePF2PATPhotonConversionVeto[i] && tightEleCheckPhotonVeto_)continue;
    if (event->elePF2PATMVA[i] < tightEleMVA_)continue;
    if (event->elePF2PATComRelIsoRho[i]/tempVec.Pt() > tightEleRelIso_)continue;
    
    electrons.push_back(i);
  }
  return electrons;
}

std::vector<int> Cuts::getLooseEles(AnalysisEvent* event){
  std::vector<int> electrons;
  for (int i = 0; i < event->numElePF2PAT; i++){
    TLorentzVector tempVec(event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]);
    if (tempVec.Pt() < looseElePt_) continue;
    if (std::abs(tempVec.Eta()) > looseEleEta_)continue;
    if (event->elePF2PATMVA[i] < looseEleMVA_)continue;
    if (event->elePF2PATComRelIsoRho[i]/tempVec.Pt() > looseEleRelIso_)continue;    

    electrons.push_back(i);
  }
  return electrons;
}

std::vector<int> Cuts::getTightMuons(AnalysisEvent* event){
  std::vector<int> muons;
  for (int i = 0; i < event->numMuonPF2PAT; i++){
    //if (i == 0) std::cout << "Starts doing first..." << event->muonPF2PATIsPFMuon[i];
    //if (i > 0) std::cout << "   " << event->muonPF2PATIsPFMuon[i];
    //    std::cout << i << "\t" << event->muonPF2PATPt[i] << "\t" << event->muonPF2PATEta[i] << "\t" << event->muonPF2PATComRelIsodBeta[i] << "\t" << event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] << "\t" << event->muonPF2PATGlobalID[i] << "\t" << event->muonPF2PATIsPFMuon[i] << "\t" << event->muonPF2PATTrackID[i] << "\t" <<  event->muonPF2PATTkLysWithMeasurements[i] << "\t" << event->muonPF2PATDBPV[i] << "\t" << event->muonPF2PATTrackNHits[i] << "\t" << event->muonPF2PATMuonNHits[i] << "\t" << event->muonPF2PATVldPixHits[i] << "\t" << event->muonPF2PATMatchedStations[i] << "\t" << std::abs(event->pvZ - event->muonPF2PATVertZ[i]) << std::endl;
    if(!event->muonPF2PATIsPFMuon[i]) continue;
    if (!event->muonPF2PATTrackID[i]) continue;
    if (!event->muonPF2PATGlobalID[i]) continue;
    //    if (!event->muonPF2PATIsPFMuon[i]) continue;
    if (event->muonPF2PATPt[i] < tightMuonPt_) continue;
    if (std::abs(event->muonPF2PATEta[i]) > tightMuonEta_) continue;
    if (event->muonPF2PATComRelIsodBeta[i] > tightMuonRelIso_) continue;

    //Do a little test of muon id stuff here.
    if (event->muonPF2PATChi2[i] < 0.) continue;
    if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] >= 10.) continue;
    if (std::abs(event->muonPF2PATDBInnerTrackD0[i]) > 0.2) continue;
    //if (event->muonPF2PATNChambers[i] < 2) continue;
    //if (i == 0) std::cout << "gets to tighter ";
    //if (i == 0) std::cout << "First muon ";
    //if (i > 0) std::cout << "Checking second muon";
    
    if (event->muonPF2PATTkLysWithMeasurements[i] < 5) continue;
    if (fabs(event->muonPF2PATDBPV[i]) >= 0.02) continue;
    //      if (event->muonPF2PATTrackNHits[i] < 11) continue;
    if (event->muonPF2PATMuonNHits[i] <= 1) continue;
    if (event->muonPF2PATVldPixHits[i] < 0) continue;
    if (event->muonPF2PATMatchedStations[i] < 1) continue;
    if (std::abs(event->pvZ - event->muonPF2PATVertZ[i]) > 0.5) continue;
    //if(i == 0) std::cout << "does first ";
    //if (i > 0) std::cout << "allows second muon";
    muons.push_back(i);
  }
  //  std::cout << muons.size() << std::endl;
  return muons;
}

std::vector<int> Cuts::getLooseMuons(AnalysisEvent* event){
  std::vector<int> muons;
  for (int i = 0; i < event->numMuonPF2PAT; i++){
    if (!event->muonPF2PATGlobalID[i] && !event->muonPF2PATTrackID[i]) continue;
    if (event->muonPF2PATPt[i] < looseMuonPt_) continue;
    if (std::abs(event->muonPF2PATEta[i]) > looseMuonEta_) continue;
    if (event->muonPF2PATComRelIsodBeta[i] > looseMuonRelIso_) continue;

    muons.push_back(i);
  }
  return muons;
}

float Cuts::getZCand(AnalysisEvent *event, std::vector<int> electrons, std::vector<int> muons){
  float closestMass = 9999.;
  //Use electrons if there are at least 2. Otherwise use muons.
  if (electrons.size() > 1){ // electrons.size() == number of electrons for selected channel
    for (unsigned int i = 0; i < electrons.size(); i++){
      for (unsigned int j = i + 1; j < electrons.size(); j++) {
	if (event->elePF2PATCharge[electrons[i]] * event->elePF2PATCharge[electrons[j]] > 0) continue;
	TLorentzVector lepton1 = TLorentzVector(event->elePF2PATGsfPx[electrons[i]],event->elePF2PATGsfPy[electrons[i]],event->elePF2PATGsfPz[electrons[i]],event->elePF2PATGsfE[electrons[i]]);
	TLorentzVector lepton2 = TLorentzVector(event->elePF2PATGsfPx[electrons[j]],event->elePF2PATGsfPy[electrons[j]],event->elePF2PATGsfPz[electrons[j]],event->elePF2PATGsfE[electrons[j]]);
	float invMass = (lepton1 + lepton2).M() -91.1;
	if (fabs(invMass) < fabs(closestMass)){
	  // set up the tlorentz vectors in the event. For plotting and jazz.
	  event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
	  event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? electrons[i]:electrons[j];
	  event->zPairRelIso.first = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[electrons[i]]/lepton1.Pt():event->elePF2PATComRelIsoRho[electrons[j]]/lepton2.Pt();
	  event->zPairRelIso.second = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[electrons[j]]/lepton2.Pt():event->elePF2PATComRelIsoRho[electrons[i]]/lepton1.Pt();
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
	    for (unsigned int k = 0; k < electrons.size(); k++){
	      if (k == i || k == j) continue;
	      event->wLepton = TLorentzVector(event->elePF2PATGsfPx[electrons[k]],event->elePF2PATGsfPy[electrons[k]],event->elePF2PATGsfPz[electrons[k]],event->elePF2PATGsfE[electrons[k]]);
	      event->wLeptonRelIso = event->elePF2PATComRelIsoRho[electrons[k]]/event->wLepton.Pt();
	      event->wLepIndex = electrons[k];
	    }
	  }
	}
      }
    }
  } else {
    for (unsigned int i = 0; i < muons.size(); i++){
      for (unsigned int j = i + 1; j < muons.size(); j++) {
	if (event->muonPF2PATCharge[muons[i]] * event->muonPF2PATCharge[muons[j]] > 0) continue;
	TLorentzVector lepton1 = TLorentzVector(event->muonPF2PATPX[muons[i]],event->muonPF2PATPY[muons[i]],event->muonPF2PATPZ[muons[i]],event->muonPF2PATE[muons[i]]);
	TLorentzVector lepton2 = TLorentzVector(event->muonPF2PATPX[muons[j]],event->muonPF2PATPY[muons[j]],event->muonPF2PATPZ[muons[j]],event->muonPF2PATE[muons[j]]);
	float invMass = (lepton1 + lepton2).M() -91;
	if (fabs(invMass) < fabs(closestMass)){
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
	    event->wLeptonRelIso = event->elePF2PATComRelIsoRho[electrons[0]]/event->wLepton.Pt();
	    event->wLepIndex = electrons[0];
	  }
	  else {
	    for (unsigned int k = 0; k < muons.size(); k++){
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

  float closestMass (9999.); 

  //Check if there are at least two electrons first. Otherwise use muons.
  
  if (electrons.size() > 1){
    for ( uint i = 0; i < electrons.size(); i++ ){
      for ( uint j = i + 1; j < electrons.size(); j++ ){
        if (event->elePF2PATCharge[electrons[i]] * event->elePF2PATCharge[electrons[j]] >= 0) continue; // check electron pair have correct charge.
        TLorentzVector lepton1 = TLorentzVector(event->elePF2PATGsfPx[electrons[i]],event->elePF2PATGsfPy[electrons[i]],event->elePF2PATGsfPz[electrons[i]],event->elePF2PATGsfE[electrons[i]]);
        TLorentzVector lepton2 = TLorentzVector(event->elePF2PATGsfPx[electrons[j]],event->elePF2PATGsfPy[electrons[j]],event->elePF2PATGsfPz[electrons[j]],event->elePF2PATGsfE[electrons[j]]);
        float invMass  = (lepton1 + lepton2).M() -91.1;
	if (fabs(invMass) < fabs(closestMass)){
        	event->zPairLeptons.first = lepton1.Pt() > lepton2.Pt()?lepton1:lepton2;
        	event->zPairIndex.first = lepton1.Pt() > lepton2.Pt() ? electrons[i]:electrons[j];
        	event->zPairRelIso.first = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[electrons[i]]/lepton1.Pt():event->elePF2PATComRelIsoRho[electrons[j]]/lepton2.Pt();
        	event->zPairRelIso.second = lepton1.Pt() > lepton2.Pt()?event->elePF2PATComRelIsoRho[electrons[j]]/lepton2.Pt():event->elePF2PATComRelIsoRho[electrons[i]]/lepton1.Pt();
        	event->zPairLeptons.second = lepton1.Pt() > lepton2.Pt()?lepton2:lepton1;
        	event->zPairIndex.second = lepton1.Pt() > lepton2.Pt() ? electrons[j]:electrons[i];
		closestMass = invMass;
      		}
	}
    } 
  }

  else if (muons.size() > 1){
    for ( uint i = 0; i < muons.size(); i++ ){
      for ( uint j = i + 1; j < muons.size(); j++ ){
	if (event->muonPF2PATCharge[muons[i]] * event->muonPF2PATCharge[muons[j]] >= 0) continue;
	TLorentzVector lepton1 = TLorentzVector(event->muonPF2PATPX[muons[i]],event->muonPF2PATPY[muons[i]],event->muonPF2PATPZ[muons[i]],event->muonPF2PATE[muons[i]]);
	TLorentzVector lepton2 = TLorentzVector(event->muonPF2PATPX[muons[j]],event->muonPF2PATPY[muons[j]],event->muonPF2PATPZ[muons[j]],event->muonPF2PATE[muons[j]]);
	float invMass = (lepton1 + lepton2).M() -91.1;
	if (fabs(invMass) < fabs(closestMass)){
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

float Cuts::getWbosonQuarksCand(AnalysisEvent *event, std::vector<int> jets){
  float closestWmass (9999.);
  if ( jets.size() > 3 )
    for ( uint k = 0; k < jets.size(); k++ ){
      for ( uint l = k + 1; l < jets.size(); l++ ){
	// check that one of the two compared jets aren't neutral.
//	std::cout << "jet charge: " << event->jetPF2PATJetCharge[jets[k]] * event->jetPF2PATJetCharge[jets[l]] << std::endl;
	if ( event->jetPF2PATJetCharge[jets[k]] * event->jetPF2PATJetCharge[jets[l]] == 0 ) continue;
	// Now ensure that the leading b jet isn't one of these!
	if ( event->jetPF2PATBDiscriminator[k] > bDiscCut_ ){
	  if( getLeadingBjetPt(event,event->bTagIndex) == event->jetPF2PATPt[jets[k]]) continue;
	}
	else if ( event->jetPF2PATBDiscriminator[l] > bDiscCut_ ){
	  if( getLeadingBjetPt(event,event->bTagIndex) == event->jetPF2PATPt[jets[l]]) continue;
	}
	TLorentzVector wQuark1 = TLorentzVector(event->jetPF2PATPx[jets[k]],event->jetPF2PATPy[jets[k]],event->jetPF2PATPz[jets[k]],event->jetPF2PATE[jets[k]]);
	TLorentzVector wQuark2 = TLorentzVector(event->jetPF2PATPx[jets[l]],event->jetPF2PATPy[jets[l]],event->jetPF2PATPz[jets[l]],event->jetPF2PATE[jets[l]]);
	float invWbosonMass = (wQuark1 + wQuark2).M() - 80.;
	if( fabs(invWbosonMass) < fabs(closestWmass) ){
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

float Cuts::getTopMass(AnalysisEvent *event, std::vector<int> bTagIndex){
  
  float leadingBjetMass = getLeadingBjetMass( event, bTagIndex );
  if (leadingBjetMass == -1.0) return -1.0;
  float topMass = ((event->wLepton).M() + leadingBjetMass + event->metPF2PATPt);
  return topMass;
}

float Cuts::getLeadingBjetMass(AnalysisEvent *event, std::vector<int> bJets){
    
    float leadingBjetPt = -1.0;
    int tempIt = 9999;
    float leadingBJetMass = 9999999.0;


    for (std::vector<int>::const_iterator lIt = bJets.begin(); lIt != bJets.end(); ++lIt){

      if (bJets[*lIt] < int(numJets_)) return -1.0;
      if (bJets[*lIt] > int(maxJets_)) return -1.0;

      if ( event->jetPF2PATPtRaw[bJets[*lIt]] > leadingBjetPt ){
	leadingBjetPt = event->jetPF2PATPtRaw[bJets[*lIt]];
	tempIt = *lIt;
      }
    }
    
    //    if ( TLorentzVector(event->jetPF2PATPx[tempIt],event->jetPF2PATPy[tempIt],event->jetPF2PATPz[tempIt],event->jetPF2PATE[tempIt]).M() < leadingBJetMass ){
    leadingBJetMass = TLorentzVector(event->jetPF2PATPx[tempIt],event->jetPF2PATPy[tempIt],event->jetPF2PATPz[tempIt],event->jetPF2PATE[tempIt]).M();
      //      }
    return leadingBJetMass;
}

float Cuts::getLeadingBjetPt(AnalysisEvent *event, std::vector<int> bJets){
    
    float leadingBjetPt = -1.0;
    int tempIt = 9999;


    for (std::vector<int>::const_iterator lIt = bJets.begin(); lIt != bJets.end(); ++lIt){

      if (bJets[*lIt] < int(numJets_)) return 9999.;
      if (bJets[*lIt] > int(maxJets_)) return 9999.;

      if ( event->jetPF2PATPtRaw[bJets[*lIt]] > leadingBjetPt ){
	leadingBjetPt = event->jetPF2PATPtRaw[bJets[*lIt]];
	tempIt = *lIt;
      }
    }
    
    return leadingBjetPt;
}


std::vector<int> Cuts::makeJetCuts(AnalysisEvent *event, int syst, float * eventWeight){
  
  std::vector<int> jets;
  float mcTag = 1., mcNoTag = 1., dataTag = 1., dataNoTag = 1., errTag = 0., errNoTag = 0., err1 = 0., err2 = 0., err3 = 0., err4 = 0.;
  //  std::cout << event->eventNum << std::endl << "Jets: " << std::endl;
  for (int i = 0; i < event->numJetPF2PAT; i++){
    //if (std::sqrt(event->jetPF2PATPx[i] * event->jetPF2PATPx[i] + event->jetPF2PATPy[i] * event->jetPF2PATPy[i]) < jetPt_) continue;
    TLorentzVector jetVec = getJetLVec(event,i,syst);
    //    std::cout << getJECUncertainty(sqrt(jetPx*jetPx + jetPy*jetPy), event->jetPF2PATEta[i],syst) << " " << syst << std::endl;
    if (jetVec.Pt() < jetPt_) continue;
    if (std::abs(jetVec.Eta()) > jetEta_) continue;
    if (event->jetPF2PATNConstituents[i] < jetNConsts_) continue;
    //std::cerr << ((event->jetPF2PATNeutralHadronEnergyFractionCorr[i] < 0.99 && event->jetPF2PATNeutralEmEnergyFractionCorr[i] < 0.99)) << "\t" << std::abs(event->jetPF2PATEta[i]) << "\t" <<  (event->jetPF2PATChargedEmEnergyFraction[i] < 0.99) <<  "\t" << (event->jetPF2PATChargedHadronEnergyFraction[i] > 0.) << "\t" << (event->jetPF2PATChargedMultiplicity[i] > 0.)  << "\t" << ((event->jetPF2PATNeutralHadronEnergyFractionCorr[i] < 0.99 && event->jetPF2PATNeutralEmEnergyFractionCorr[i] < 0.99) && ((std::abs(event->jetPF2PATEta[i]) > 2.4) || (event->jetPF2PATChargedEmEnergyFraction[i] < 0.99 && event->jetPF2PATChargedHadronEnergyFraction[i] > 0. && event->jetPF2PATChargedMultiplicity[i] > 0.))) << std::endl;
    if (jetIDDo_ && !((event->jetPF2PATNeutralHadronEnergyFractionCorr[i] < 0.99 && event->jetPF2PATNeutralEmEnergyFractionCorr[i] < 0.99) && ((std::abs(event->jetPF2PATEta[i]) > 2.4) || (event->jetPF2PATChargedEmEnergyFraction[i] < 0.99 && event->jetPF2PATChargedHadronEnergyFraction[i] > 0. && event->jetPF2PATChargedMultiplicity[i] > 0.)))) continue;
    double deltaLep = 10000.;
//    double deltaQuark = 10000.;
    //Testing out if electron only cleaning works.
    /*    if (event->electronIndexTight.size() > 1){
      if (deltaLep > deltaR(event->zPairLeptons.first.Eta(),event->zPairLeptons.first.Phi(),event->jetPF2PATEta[i],event->jetPF2PATPhi[i]))
      deltaLep = deltaR(event->zPairLeptons.first.Eta(),event->zPairLeptons.first.Phi(),event->jetPF2PATEta[i],event->jetPF2PATPhi[i]);
    if (deltaLep > deltaR(event->zPairLeptons.second.Eta(),event->zPairLeptons.second.Phi(),event->jetPF2PATEta[i],event->jetPF2PATPhi[i]))
      deltaLep = deltaR(event->zPairLeptons.second.Eta(),event->zPairLeptons.second.Phi(),event->jetPF2PATEta[i],event->jetPF2PATPhi[i]);
    }

    else if (event->electronIndexTight.size() > 0){
      if (deltaLep > deltaR(event->wLepton.Eta(),event->wLepton.Phi(),event->jetPF2PATEta[i],event->jetPF2PATPhi[i]))
	deltaLep = deltaR(event->wLepton.Eta(),event->wLepton.Phi(),event->jetPF2PATEta[i],event->jetPF2PATPhi[i]);
	}*/
    //    else { deltaLep = 10000.;}

    if (deltaLep > deltaR(event->zPairLeptons.first.Eta(),event->zPairLeptons.first.Phi(),jetVec.Eta(),jetVec.Phi()))
      deltaLep = deltaR(event->zPairLeptons.first.Eta(),event->zPairLeptons.first.Phi(),jetVec.Eta(),jetVec.Phi());
    if (deltaLep > deltaR(event->zPairLeptons.second.Eta(),event->zPairLeptons.second.Phi(),jetVec.Eta(),jetVec.Phi()))
      deltaLep = deltaR(event->zPairLeptons.second.Eta(),event->zPairLeptons.second.Phi(),jetVec.Eta(),jetVec.Phi());

    if (trileptonChannel_ == true && deltaLep > deltaR(event->wLepton.Eta(),event->wLepton.Phi(),jetVec.Eta(),jetVec.Phi()))
      deltaLep = deltaR(event->wLepton.Eta(),event->wLepton.Phi(),jetVec.Eta(),jetVec.Phi());

/*
    if (!trileptonChannel_ && deltaQuark > deltaR(event->wPairQuarks.first.Eta(),event->wPairQuarks.first.Phi(),jetVec.Eta(),jetVec.Phi()))
      deltaQuark = deltaR(event->wPairQuarks.first.Eta(),event->wPairQuarks.first.Phi(),jetVec.Eta(),jetVec.Phi());
	
    if (!trileptonChannel_ && deltaQuark > deltaR(event->wPairQuarks.second.Eta(),event->wPairQuarks.second.Phi(),jetVec.Eta(),jetVec.Phi()))
      deltaQuark = deltaR(event->wPairQuarks.second.Eta(),event->wPairQuarks.second.Phi(),jetVec.Eta(),jetVec.Phi());
*/

    //std::cout << event->jetPF2PATPtRaw[i] << " " << deltaLep << std::endl;
  if (deltaLep < 0.5) continue;
//    if (deltaQuark < 1.0 && !trileptonChannel_) continue;
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
    jets.push_back(i);

    if (getBTagWeight_){
      getBWeight(event,jetVec,i,&mcTag,&mcNoTag,&dataTag,&dataNoTag,&errTag,&errNoTag,&err1,&err2,&err3,&err4);
    }
  }
  //Evaluate b-tag weight for event here.
  if (getBTagWeight_){
    float bWeight = (dataNoTag * dataTag)/(mcNoTag * mcTag);
    if (mcNoTag == 0 || mcTag == 0 || dataNoTag == 0 || dataTag == 0 || mcNoTag != mcNoTag || mcTag != mcTag || dataTag != dataTag || dataNoTag != dataNoTag)
      bWeight = 1.;
    float bWeightErr = std::sqrt( pow(err1+err2,2) + pow(err3 + err4, 2)) * bWeight;
    if (syst == 256)
      bWeight += bWeightErr;
    if (syst == 512)
      bWeight -= bWeightErr;
    *eventWeight *= bWeight;
  }
  return jets;
}

std::vector<int> Cuts::makeBCuts(AnalysisEvent* event, std::vector<int> jets){
  
  std::vector<int> bJets;

  for (unsigned int i = 0; i < jets.size(); i++){
    if (singleEventInfoDump_) std::cout << event->jetPF2PATPtRaw[jets[i]] << " " << event->jetPF2PATBDiscriminator[jets[i]] << std::endl;
    if (event->jetPF2PATBDiscriminator[jets[i]] < bDiscCut_) continue;

    bJets.push_back(i);
  }
  return bJets;
}

void Cuts::setTightEle(float,float,float)
{
}

//This is only called if the thing is data. There's also a little clause to run over certain triggers if they exist. Because I failed miserably to get them all first time through processing...
bool Cuts::triggerCuts(AnalysisEvent* event){

  if (skipTrigger_) return true;

  //MuEG triggers
  bool muEGTrig = false;
  if (!isMC_) {if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2 > 0 || event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3 > 0 ) muEGTrig = true;}
  else if ( event->HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 > 0 || event->HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1 > 0 ){
    muEGTrig = true;
    //    std::cout << "muEGTrig for MC fired." << std::endl;
  }

  //double electron triggers
  bool eeTrig = false;
  //  std::cout << "event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 : " << event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 << std::endl;
  if (!isMC_) {if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0 || event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0 ) eeTrig = true;}
  else if ( event->HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0 ) {
    eeTrig = true;
    //    std::cout << "eeTrig for MC fired." << std::endl;   
  }

  //double muon triggers
  bool mumuTrig = false;
  if (!isMC_) {if ( event->HLT_IsoMu20_v2  > 0 || event->HLT_IsoMu20_eta2p1_v2 > 0 || event->HLT_IsoMu18_v1 > 0 ) mumuTrig = true;}

  else if (event->HLT_IsoMu20_v1 > 0 || event->HLT_IsoMu20_eta2p1_v1 > 0){
    mumuTrig = true;
    //    std::cout << "mumuTrig for MC fired." << std::endl;  
  }

  if (cutConfTrigLabel_.find("d") != std::string::npos){if (muEGTrig) return true;}
  if (cutConfTrigLabel_.find("e") != std::string::npos){if (eeTrig && !(muEGTrig || mumuTrig)) return true;}
  if (cutConfTrigLabel_.find("m") != std::string::npos){if (mumuTrig && !(eeTrig || muEGTrig)) return true;}
  return false;
}

//Does the cuts required for the synchronisation.
bool Cuts::synchCuts(AnalysisEvent* event){
  //Trigger stuff would go first, but in MC at least (what I'm starting with) I don't have that saved. Idiot.

  if ( trileptonChannel_ == false ){
    std::cout << "Not setup for dilepton synch cuts yet. Exiting program!" << std::endl;
    exit(2);
  }
  
  if (singleEventInfoDump_){
    std::cout << std::setprecision(6) << std::fixed;
  }
  
  int looseLeps = getLooseLepsNum(event);
  if (isMC_ && looseLeps < 2) return false;
  if (!isMC_ && looseLeps < 3) return false;
  if (/*!isMC_ &&*/ !triggerCuts(event)) return false; // Commented out this so that it runs on both MC and data.
  synchCutFlowHist_->Fill(0.5);
  if (makeEventDump_) {step0EventDump_ << event->eventNum << std::endl;}
  
  event->electronIndexTight = getTightEles(event);
  event->muonIndexTight = getTightMuons(event);
  synchNumEles_->Fill(event->electronIndexTight.size());
  synchNumMus_->Fill(event->muonIndexTight.size());
  //If electrons are expected to be the Z, check there are an oppositely charged pair.
  bool zCand = false;
  if (numTightEle_ > 1){
    for (unsigned int i = 0; i < event->electronIndexTight.size();++i){
      for (unsigned int j = i +1; j < event->electronIndexTight.size();++j){
	if (event->elePF2PATCharge[event->electronIndexTight[i]] * event->elePF2PATCharge[event->electronIndexTight[j]] > 0) continue;
	zCand = true;
      }
    }

  }
  else if (numTightMu_ > 1){
    for (unsigned int i = 0; i < event->muonIndexTight.size();++i){
      for (unsigned int j = i +1; j < event->muonIndexTight.size();++j){
	if (event->muonPF2PATCharge[event->muonIndexTight[i]] * event->muonPF2PATCharge[event->muonIndexTight[j]] > 0) continue;
	zCand = true;
      }
    }
  }
  if (!zCand) return false;
  synchCutFlowHist_->Fill(1.5); 
  if (singleEventInfoDump_) std::cout << "Gets past z cand stuff with " << event->electronIndexTight.size() << " tight electrons and " << event->muonIndexTight.size() << " tight muons" << std::endl;
  if (event->electronIndexTight.size() != numTightEle_) return false;
  if (event->muonIndexTight.size() != numTightMu_) return false;
  //loose lepton veto
  if (singleEventInfoDump_) std::cout << "Correct number of leptons and loose: " << getLooseEles(event).size() << " " << getLooseMuons(event).size() << std::endl;
  if (event->electronIndexTight.size() != getLooseEles(event).size()) return false;
  if (event->muonIndexTight.size() != getLooseMuons(event).size()) return false;
  if (singleEventInfoDump_) std::cout << " and passes veto too." << std::endl;
  synchCutFlowHist_->Fill(2.5);
  if (makeEventDump_){dumpToFile(event,2);}
  if (singleEventInfoDump_) std::cout << "Z mass: " << getZCand(event,event->electronIndexTight,event->muonIndexTight) << std::endl;
  if (fabs(getZCand(event,event->electronIndexTight,event->muonIndexTight)) > invZMassCut_) return false;
  synchCutFlowHist_->Fill(3.5);

  //Add in extra steps here.
  float tempWeight = 1.;
  event->jetIndex = makeJetCuts(event, 0, &tempWeight);
  if (singleEventInfoDump_) std::cout << "Number of jets: " << event->jetIndex.size() << std::endl;
  if (event->jetIndex.size() < 1) return false;
  if (makeEventDump_){
    dumpToFile(event,4);
  }
  //  std::cout << event->jetIndex.size() << std::endl;
  synchCutFlowHist_->Fill(4.5);
  event->bTagIndex = makeBCuts(event,event->jetIndex);
  synchCutTopMassHist_->Fill(getTopMass(event, event->jetIndex)); // Plot top mass distribution for all top candidates - all sanity checks done, Z mass exists, got b jets too.
  if (singleEventInfoDump_) std::cout << "One bJet: " << event->bTagIndex.size() << std::endl;
  if (!event->bTagIndex.size() == 1) return false;
  synchCutFlowHist_->Fill(5.5);
  if (singleEventInfoDump_) std::cout << "MET: " << event->metPF2PATPt << std::endl;
  if (event->metPF2PATPt < metCut_) return false;
  synchCutFlowHist_->Fill(6.5);
  if (singleEventInfoDump_) std::cout << "mTW: " << event->metPF2PATPt << std::endl;
  if (std::sqrt(2*event->metPF2PATPt*event->wLepton.Pt()*(1-cos(event->metPF2PATPhi - event->wLepton.Phi()))) < mTWCut_) return false;
  synchCutFlowHist_->Fill(7.5);
  if (singleEventInfoDump_) std::cout << "top mass cut: " << getTopMass(event, event->jetIndex)  << std::endl;
  if (getTopMass(event, event->jetIndex) > TopMassCutUpper_ || getTopMass(event, event->jetIndex) < TopMassCutLower_) return false;
  synchCutFlowHist_->Fill(8.5);
  if (singleEventInfoDump_) std::cout << "Passes all cuts! Yay!" << std::endl;
  if (makeEventDump_){
    dumpToFile(event,6);
  }
  if (singleEventInfoDump_) std::cout << std::setprecision(1) << std::fixed;
  return true;
}

TH1F* Cuts::getSynchCutFlow(){
  std::cout << "Eles: " << numTightEle_ << " Muons: " << numTightMu_ << std::endl;
  char const *names[] = {"Trigger","Lep Pair", "3 Leptons", "zMass","1 jet","1 b-tag","MET","mTW", "topMass"};
  for (unsigned int i = 1; i < 10; ++i){
    std::cout << names[i-1] << ": " << synchCutFlowHist_->GetBinContent(i) << std::endl;
  }
  std::cout << "The number of leptons in the passing trigger events:" << std::endl;
  std::cout << "Leps\tEle\tMuon" << std::endl;
  for (unsigned int i = 0; i < 10; i++){
    std::cout << i << "\t" << synchNumEles_->GetBinContent(i) << "\t" << synchNumMus_->GetBinContent(i) << std::endl;
  }
  char const *labels[] = {"In", "ID", "PtEtaIso", "chi2","tklay","DBPV","TrackHit","MuHits","PixHits","MtStats","DZPV"};
  for (unsigned int i = 1; i < 12; ++i){
    std::cout << labels[i-1] << ": \t" << synchMuonCutFlow_->GetBinContent(i) << std::endl;
  }

  synchCutFlowHist_->SaveAs("plots/synch/synchCutFlowHist.root");
  synchCutTopMassHist_->SaveAs("plots/synch/synchCutTopMassHist.root");

  return synchCutFlowHist_;
}

//Method used for the synchronisation. Mimics the IPHC preselection skims.
int Cuts::getLooseLepsNum(AnalysisEvent* event){
  return getLooseElecs(event) + getLooseMus(event);
}

int Cuts::getLooseElecs(AnalysisEvent* event){
  int looseLeps = 0.;
  for (int i = 0; i < event->numElePF2PAT; i++){
    if (!event->elePF2PATIsGsf[i]) continue;
    TLorentzVector tempVec(event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]);
    if (tempVec.Pt() < 15) continue;
    if (std::abs(tempVec.Eta()) > 2.5)continue;
    looseLeps++;
  }
  return looseLeps;
}

int Cuts::getLooseMus(AnalysisEvent* event){
  int looseLeps = 0.;
  for (int i = 0; i < event->numMuonPF2PAT; i++){
    if (!event->muonPF2PATGlobalID[i]||!event->muonPF2PATTrackID[i]) continue;
    if (event->muonPF2PATPt[i] < 15) continue;
    if (std::abs(event->muonPF2PATEta[i]) > 2.5) continue;
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
  float invMass = 100.;
  if (numTightEle_ > 1){
    if (event->electronIndexTight.size() < 2) return false;
    if (event->elePF2PATCharge[event->electronIndexTight[0]] * event->elePF2PATCharge[event->electronIndexTight[1]] > 0) return false;
    TLorentzVector lep1 = TLorentzVector(event->elePF2PATGsfPx[event->electronIndexTight[0]],event->elePF2PATGsfPy[event->electronIndexTight[0]],event->elePF2PATGsfPz[event->electronIndexTight[0]],event->elePF2PATGsfE[event->electronIndexTight[0]]);
    TLorentzVector lep2 = TLorentzVector(event->elePF2PATGsfPx[event->electronIndexTight[1]],event->elePF2PATGsfPy[event->electronIndexTight[1]],event->elePF2PATGsfPz[event->electronIndexTight[1]],event->elePF2PATGsfE[event->electronIndexTight[1]]); 
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
  std::vector<int> invIsoEle = getInvIsoEles(event);
  std::vector<int> invIsoMus = getInvIsoMuons(event);
  
  //Check we have the right number of leptons.
  if (event->electronIndexTight.size() + invIsoEle.size() != numTightEle_) return false;
  if (event->muonIndexTight.size() + invIsoMus.size() != numTightMu_) return false;

  //Debugging
  /*  std::cout << event->numElePF2PAT << " " << event->electronIndexTight.size() << " " << invIsoEle.size();
  std::cout << " tight index: ";
  for (unsigned int i = 0; i < event->electronIndexTight.size(); i++) std:: cout << " " << event->electronIndexTight[i];
  std::cout << " inv index: ";
  for (unsigned int i = 0; i < invIsoEle.size(); i++) std::cout << " " << invIsoEle[i] << " " << "relIso: " << event->elePF2PATComRelIsoRho[invIsoEle[i]]/event->elePF2PATPT[invIsoEle[i]] ;
  std::cout << std::endl;*/

  //Put extra lepton into W boson thing.
  if (invIsoEle.size() == 1){
    event->wLepton = TLorentzVector(event->elePF2PATGsfPx[invIsoEle[0]],event->elePF2PATGsfPy[invIsoEle[0]],event->elePF2PATGsfPz[invIsoEle[0]],event->elePF2PATGsfE[invIsoEle[0]]);
    event->wLepIndex = invIsoEle[0];
    event->wLeptonRelIso = event->elePF2PATComRelIsoRho[invIsoEle[0]]/event->wLepton.Pt();;
  }
  else{
    event->wLepton = TLorentzVector(event->muonPF2PATPX[invIsoMus[0]],event->muonPF2PATPY[invIsoMus[0]],event->muonPF2PATPZ[invIsoMus[0]],event->muonPF2PATE[invIsoMus[0]]);
    event->wLepIndex = invIsoMus[0];
    event->wLeptonRelIso = event->muonPF2PATComRelIsodBeta[invIsoMus[0]];
  }

  if(doPlots_) plotMap["lepSel"]->fillAllPlots(event,*eventWeight);
  if(doPlots_||fillCutFlow_){
    cutFlow->Fill(0.5,*eventWeight);
  }

  if (std::abs(invMass) > invZMassCut_) return false;
  if(doPlots_||fillCutFlow_) cutFlow->Fill(1.5,*eventWeight);
  return true;

}

std::vector<int> Cuts::getInvIsoEles(AnalysisEvent* event) {                                                                                                                                                                                  
  std::vector<int> electrons;
  for (int i = 0; i < event->numElePF2PAT; i++){
    bool same = false;
    /*    for (int j = 0; j < event->electronIndexTight.size(); j++){
      if (i == event->electronIndexTight[j]) same = true;
      }*/
    if (same) continue;
    if (!event->elePF2PATIsGsf[i]) continue;
    TLorentzVector tempVec(event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]);
    if (tempVec.Pt() < tightElePt_) continue;
    if (std::abs(tempVec.Eta()) > tightEleEta_)continue;
    if (std::abs(event->elePF2PATBeamSpotCorrectedTrackD0[i]) > tightEled0_)continue;
    if (event->elePF2PATMissingInnerLayers[i] > tightEleMissLayers_) continue;
    if (!event->elePF2PATPhotonConversionVeto[i] && tightEleCheckPhotonVeto_)continue;
    if (event->elePF2PATMVA[i] < tightEleMVA_)continue;
    if (event->elePF2PATComRelIsoRho[i]/tempVec.Pt() < tightEleRelIso_)continue;

    electrons.push_back(i);
  }
  return electrons;
}

std::vector<int> Cuts::getInvIsoMuons(AnalysisEvent* event){
  std::vector<int> muons;
  for (int i = 0; i < event->numMuonPF2PAT; i++){
    if (!event->muonPF2PATGlobalID[i] && !event->muonPF2PATTrackID[i]) continue;
    if (event->muonPF2PATPt[i] < tightMuonPt_) continue;
    if (std::abs(event->muonPF2PATEta[i]) > tightMuonEta_) continue;
    if (event->muonPF2PATComRelIsodBeta[i] < tightMuonRelIso_) continue;

    //Do a little test of muon id stuff here.
    if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] > 10.) continue;
    if (std::abs(event->muonPF2PATDBInnerTrackD0[i]) > 0.2) continue;
    if (event->muonPF2PATNChambers[i] < 2) continue;

    muons.push_back(i);
  }
  return muons;
}

std::vector<int> Cuts::synchTightMuons(AnalysisEvent* event){
  std::vector<int> muons;
  for (int i = 0; i < event->numMuonPF2PAT; i++){
    synchMuonCutFlow_->Fill(0.5);
    //if (i == 0) std::cout << "Starts doing first..." << event->muonPF2PATIsPFMuon[i];
    //if (i > 0) std::cout << "   " << event->muonPF2PATIsPFMuon[i];
    if (!event->muonPF2PATGlobalID[i] && !event->muonPF2PATTrackID[i]) continue;
    synchMuonCutFlow_->Fill(1.5);
    //    if (!event->muonPF2PATIsPFMuon[i]) continue;
    if (event->muonPF2PATPt[i] < tightMuonPt_) continue;
    if (std::abs(event->muonPF2PATEta[i]) > tightMuonEta_) continue;
    if (event->muonPF2PATComRelIsodBeta[i] > tightMuonRelIso_) continue;
    synchMuonCutFlow_->Fill(2.5);
    //Do a little test of muon id stuff here.
    if (event->muonPF2PATChi2[i] < 0.) continue;
    if (event->muonPF2PATChi2[i]/event->muonPF2PATNDOF[i] >= 10.) continue;
    //    if (std::abs(event->muonPF2PATDBInnerTrackD0[i]) > 0.2) continue;
    //if (event->muonPF2PATNChambers[i] < 2) continue;
    //if (i == 0) std::cout << "gets to tighter ";
    synchMuonCutFlow_->Fill(3.5);
    //if (i == 0) std::cout << "First muon ";
    //if (i > 0) std::cout << "Checking second muon";
    //      std::cout << event->muonPF2PATTkLysWithMeasurements[i] << "\t" << event->muonPF2PATDBPV[i] << "\t" << event->muonPF2PATTrackNHits[i] << "\t" << event->muonPF2PATMuonNHits[i] << "\t" << event->muonPF2PATVldPixHits[i] << "\t" << event->muonPF2PATMatchedStations[i] << "\t" << std::abs(event->pvZ - event->muonPF2PATVertZ[i]) << std::endl;
    if (event->muonPF2PATTkLysWithMeasurements[i] < 5 ) continue;
    synchMuonCutFlow_->Fill(4.5);
    if (event->muonPF2PATDBPV[i] > 0.02) continue;
    synchMuonCutFlow_->Fill(5.5);
    if (event->muonPF2PATTrackNHits[i] < 11) continue;
    synchMuonCutFlow_->Fill(6.5);
    if (event->muonPF2PATMuonNHits[i] <= 1) continue;
    synchMuonCutFlow_->Fill(7.5);
    if (event->muonPF2PATVldPixHits[i] < 0) continue;
    synchMuonCutFlow_->Fill(8.5);
    if (event->muonPF2PATMatchedStations[i] < 1) continue;
    synchMuonCutFlow_->Fill(9.5);
    if (std::abs(event->pvZ - event->muonPF2PATVertZ[i]) > 0.5) continue;
    synchMuonCutFlow_->Fill(10.5);
    //if(i == 0) std::cout << "does first ";
    //if (i > 0) std::cout << "allows second muon";
    muons.push_back(i);
  }
  //  std::cout << muons.size() << std::endl;
  return muons;
}

//For synchronisation I am dumping the lepton information here.
void Cuts::dumpLeptonInfo(AnalysisEvent* event){
  std::cout << std::setprecision(6) << std::fixed; 
  std::cerr << std::setprecision(6) << std::fixed; 
  //Dump electron info first
  event->electronIndexTight = getTightEles(event);
  event->muonIndexTight = getTightMuons(event); 
  std::cout << "Electrons: found " << event->electronIndexTight.size() << std::endl;
  for (unsigned int i = 0; i < event->electronIndexTight.size(); i++){
    TLorentzVector tempVec(event->elePF2PATGsfPx[event->electronIndexTight[i]],event->elePF2PATGsfPy[event->electronIndexTight[i]],event->elePF2PATGsfPz[event->electronIndexTight[i]],event->elePF2PATGsfE[event->electronIndexTight[i]]);
    //    std::cout << i << " | " << event->elePF2PATPT[event->electronIndexTight[i]];
    std::cout << " | " << tempVec.Pt();
    //std::cout << " | " << event->elePF2PATEta[event->electronIndexTight[i]];
    std::cout << " | " << tempVec.Eta();
    std::cout << " | " << tempVec.Phi();
    std::cout << " | " << event->elePF2PATPhi[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATD0PV[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATMVA[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATMissingInnerLayers[event->electronIndexTight[i]];
    std::cout << " | " << event->elePF2PATComRelIsoRho[event->electronIndexTight[i]]/tempVec.Pt();
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
  for (unsigned int i = 0; i <  event->muonIndexTight.size(); i++){
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
  int numbJets = event->numJetPF2PAT > 4?4:event->numJetPF2PAT;
  std::cout << "Jets: " << event->numJetPF2PAT << std::endl;
  for (int i = 0; i < numbJets; i++){
    std::cout << i;
    //    std::cout << " | " << event->jetPF2PATPt[i];
    //std::cout << " | " << std::sqrt(event->jetPF2PATPx[i] * event->jetPF2PATPx[i] + event->jetPF2PATPy[i] * event->jetPF2PATPy[i]);
    // std::cout << " | " << event->jetPF2PATEt[i];
    std::cout << " | " << event->jetPF2PATPtRaw[i];
    //std::cout << " | " << event->jetPF2PATUnCorPt[i];
    std::cout << " | " << event->jetPF2PATEta[i];
    std::cout << " | " << event->jetPF2PATPhi[i];
    std::cout << " | " << event->jetPF2PATNConstituents[i];
    std::cout << " | " << (event->jetPF2PATNeutralHadronEnergyFractionCorr[i] < 0.99 && event->jetPF2PATNeutralEmEnergyFractionCorr[i] < 0.99) && ((std::abs(event->jetPF2PATEta[i]) > 2.4) || (event->jetPF2PATChargedEmEnergyFraction[i] < 0.99 && event->jetPF2PATChargedHadronEnergyFraction[i] > 0. && event->jetPF2PATChargedMultiplicity[i] > 0.));
    std::cout << " | " << event->jetPF2PATdRClosestLepton[i];
    std::cout << std::endl;
  }
  std::cout << "MET: " << event->metPF2PATEt << " | " << event->metPF2PATPt << std::endl;

  std::cout << std::setprecision(1) << std::fixed; 
  std::cerr << std::setprecision(1) << std::fixed; 
}

void Cuts::dumpLooseLepInfo(AnalysisEvent* event){
  std::cout << std::setprecision(6) << std::fixed; 
  std::cerr << std::setprecision(6) << std::fixed; 

  std::cout << "Electrons: " << event->numElePF2PAT << std::endl;
  for (int i = 0; i < event->numElePF2PAT; i++){

    TLorentzVector tempVec(event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]);
    std::cout << i;
    std::cout << " | " << tempVec.Pt();
    std::cout << " | " << tempVec.Eta();
    std::cout << " | " << event->elePF2PATSCEta[i];
    std::cout << " | " << event->elePF2PATComRelIsoRho[i]/tempVec.Pt();
    std::cout << std::endl;
  }
  std::cout << "Muons: " << event->numMuonPF2PAT << std::endl;
  for (int i = 0; i < event->numMuonPF2PAT; i++){
    std::cout << i;
    std::cout << " | " << event->muonPF2PATPt[i];
    std::cout << " | " << event->muonPF2PATEta[i];
    std::cout << " | " << event->muonPF2PATComRelIsodBeta[i];
    std::cout << " | " << event->muonPF2PATGlobalID[i];
    std::cout << " | " << event->muonPF2PATTrackID[i];
    std::cout << std::endl;
  }

  std::cout << std::setprecision(1) << std::fixed; 
  std::cerr << std::setprecision(1) << std::fixed; 
}

double Cuts::deltaR(float eta1, float phi1, float eta2, float phi2){
  double dEta = eta1-eta2;
  double dPhi = phi1-phi2;
  while (fabs(dPhi) > 3.14159265359){
    dPhi += (dPhi > 0.? -2*3.14159265359:2*3.14159265359);
  }
  //  if(singleEventInfoDump_)  std::cout << eta1 << " " << eta2 << " phi " << phi1 << " " << phi2 << " ds: " << eta1-eta2 << " " << phi1-phi2 << " dR: " << std::sqrt((dEta*dEta)+(dPhi*dPhi)) << std::endl;
  return std::sqrt((dEta*dEta)+(dPhi*dPhi));
}

//For dumping contents of step 4. Bit more complicated than old, so doing it elsewhere.
void Cuts::dumpToFile(AnalysisEvent* event, int step){
  std::vector<TLorentzVector> tempLepVec;
  if (step > 2 && trileptonChannel_ == true){
    if (event->zPairLeptons.first.Pt() < event->wLepton.Pt()){
      tempLepVec.push_back(event->wLepton);
      tempLepVec.push_back(event->zPairLeptons.first);
      tempLepVec.push_back(event->zPairLeptons.second);
    }
    else if (event->wLepton.Pt() > event->zPairLeptons.second.Pt()){
      tempLepVec.push_back(event->zPairLeptons.first);
      tempLepVec.push_back(event->wLepton);
      tempLepVec.push_back(event->zPairLeptons.second);
    }
    else {
      tempLepVec.push_back(event->zPairLeptons.first);
      tempLepVec.push_back(event->zPairLeptons.second);
      tempLepVec.push_back(event->wLepton);
    }
  }
  if (step==2 && trileptonChannel_ == true){
    int muUsed = 0;
    for (int i = 0; i < event->numElePF2PAT; i++){
      TLorentzVector tempVec(event->elePF2PATGsfPx[i],event->elePF2PATGsfPy[i],event->elePF2PATGsfPz[i],event->elePF2PATGsfE[i]);
      for (int j = muUsed; j < event->numMuonPF2PAT; j++){
	if (tempVec.Pt() > event->muonPF2PATPt[i]){
	  tempLepVec.push_back(tempVec);
	  break;
	}
	else {
	  tempLepVec.push_back(TLorentzVector(event->muonPF2PATPX[j],event->muonPF2PATPY[j],event->muonPF2PATPZ[j],event->muonPF2PATE[j]));
	  muUsed++;
	  if (tempLepVec.size() > 2) break;
	}
	if (tempLepVec.size() > 2) break;
      }
    }
    if (tempLepVec.size() < 3){
      for (int j = 0; j < 3; j++){
	tempLepVec.push_back(TLorentzVector(event->muonPF2PATPX[j],event->muonPF2PATPY[j],event->muonPF2PATPZ[j],event->muonPF2PATE[j]));
      }
    }
  }
  switch (step) {
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
  for (unsigned int i = 0; i < 3; i++){
    switch (step) {
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
  float tempWeight = 1.;
  event->jetIndex = makeJetCuts(event, 0, &tempWeight);
  for (unsigned int i = 0; i < 4; i++){
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

float Cuts::getLeptonWeight(AnalysisEvent * event){
  //If number of electrons is > 1  then both z pair are electrons, so get their weight
  if (!isMC_) return 1.;
  float leptonWeight = 1.;
  if (trileptonChannel_ == true){
    if (numTightEle_ > 1){
      leptonWeight *= eleSF(event->zPairLeptons.first.Pt(),fabs(event->zPairLeptons.first.Eta()));
      leptonWeight *= eleSF(event->zPairLeptons.second.Pt(),fabs(event->zPairLeptons.second.Eta()));
    }
    else{
      leptonWeight *= muonSF(event->zPairLeptons.first.Pt(),fabs(event->zPairLeptons.first.Eta()));
      leptonWeight *= muonSF(event->zPairLeptons.second.Pt(),fabs(event->zPairLeptons.second.Eta()));
    }
    if (numTightEle_ == 3 || numTightEle_ == 1){
      leptonWeight *= eleSF(event->wLepton.Pt(),fabs(event->wLepton.Eta()));
    }
    else{
      leptonWeight *= eleSF(event->wLepton.Pt(),fabs(event->wLepton.Eta()));
    }
  }

  else if(trileptonChannel_ == false){
    if (numTightEle_ == 2){
      leptonWeight *= eleSF(event->zPairLeptons.first.Pt(),fabs(event->zPairLeptons.first.Eta()));
      leptonWeight *= eleSF(event->zPairLeptons.second.Pt(),fabs(event->zPairLeptons.second.Eta()));
    }
    else if (numTightMu_ == 2){
      leptonWeight *= muonSF(event->zPairLeptons.first.Pt(),fabs(event->zPairLeptons.first.Eta()));
      leptonWeight *= muonSF(event->zPairLeptons.second.Pt(),fabs(event->zPairLeptons.second.Eta()));
    }
  }

  else{
    std::cout << "You've somehow managed to set the trilepton/dilepton bool flag to a value other than these two!" << std::endl;
    std::cout << "HOW?! Well done for breaking this ..." << std::endl;
    exit(0);
  }

  return leptonWeight;
  
}

float Cuts::eleSF(float pt, float eta){
  if ( pt < 30 ){
    if (eta < 0.8) return 0.969;
    else if (eta < 1.4442) return 0.935;
    else if (eta < 1.5660) return 1.032;
    else if (eta < 2.5) return 0.919;
  }
  else if (pt < 40){
    if (eta < 0.8) return 0.926;        
    else if (eta < 1.4442) return 0.945;
    else if (eta < 1.5660) return 0.907;
    else if (eta < 2.5) return 0.926;   
  }
  else if (pt < 50){
    if (eta < 0.8) return 0.969;        
    else if (eta < 1.4442) return 0.964;
    else if (eta < 1.5660) return 0.957;
    else if (eta < 2.5) return 0.952;   
  }
  else{
    if (eta < 0.8) return 0.975;        
    else if (eta < 1.4442) return 0.974;
    else if (eta < 1.5660) return 0.877;
    else if (eta < 2.5) return 0.950;   
  }
  return 1.;
}

float Cuts::muonSF(float pt, float eta){
  if (eta < 0.9) return 0.9925;
  else if (eta < 1.2) return 0.9928 * 1.0006;
  else if (eta < 2.1) return 0.9960 * 1.0006;
  else return 0.9952 * 1.0004;
}

void Cuts::initialiseJECCors(){
  std::ifstream jecFile("scripts/Summer15_25nsV6M3_DATA_Uncertainty_AK4PFchs.txt");
  std::string line;
  bool first = true;

  if (!jecFile.is_open()){
    std::cout << "Unable to open jecFile." << std::endl;
    exit(0);
  }

  while(getline(jecFile,line)){
    std::vector<std::string> tempVec;
    std::stringstream lineStream(line);
    std::string item;
    while (std::getline(lineStream,item,' ')){
      tempVec.push_back(item);
    }
    std::vector<float> tempUp;
    std::vector<float> tempDown;

    etaMinJEC_.push_back(atof(tempVec[0].c_str()));
    etaMaxJEC_.push_back(atof(tempVec[1].c_str()));
    for (unsigned int i = 1; i < tempVec.size()/3; i++){
      unsigned int ind = i * 3;
      if (first){
      	ptMinJEC_.push_back(atof(tempVec[ind].c_str()));
	ptMaxJEC_.push_back((ind+3 >= tempVec.size()?10000.:atof(tempVec[ind+3].c_str())));
      }
      tempUp.push_back(atof(tempVec[ind+1].c_str()));
      tempDown.push_back(atof(tempVec[ind+2].c_str()));
    }
    jecSFUp_.push_back(tempUp);
    jecSFDown_.push_back(tempDown);
    first = false;
  }
}

float Cuts::getJECUncertainty(float pt, float eta, int syst){
  if (!(syst == 4 || syst == 8)){
    return 0.;
  }
  unsigned int ptBin = 0, etaBin = 0;
  for (unsigned int i = 0; i < ptMinJEC_.size(); i++){
    if (pt > ptMinJEC_[i] && pt < ptMaxJEC_[i]){
      ptBin = i;
      break;
    }
  }
  for (unsigned int i = 0; i < etaMinJEC_.size(); i++){
    if (eta > etaMinJEC_[i] && eta < etaMaxJEC_[i]){
      etaBin = i;
      break;
    }
  }

  float lowFact = (syst==4?jecSFUp_[etaBin][ptBin]:jecSFDown_[etaBin][ptBin]);
  float hiFact = (syst==4?jecSFUp_[etaBin][ptBin+1]:jecSFDown_[etaBin][ptBin+1]);
  //Now do some interpolation
  float a = (hiFact - lowFact)/(ptMaxJEC_[ptBin]-ptMinJEC_[ptBin]);
  float b = (lowFact * (ptMaxJEC_[ptBin]) - hiFact*ptMinJEC_[ptBin])/(ptMaxJEC_[ptBin] - ptMinJEC_[ptBin]);
  return (syst == 4 ? a * pt + b : -(a * pt + b));
  
}

TLorentzVector Cuts::getJetLVec(AnalysisEvent* event, int index, int syst){
  
  float oldSmearCorr = 0.;
  if (fabs(event->jetPF2PATEta[index]) <= 1.1) oldSmearCorr = 0.066;
  else if (fabs(event->jetPF2PATEta[index]) <= 1.7) oldSmearCorr = 0.191;
  else if (fabs(event->jetPF2PATEta[index]) <= 2.3) oldSmearCorr = 0.096;
  else oldSmearCorr = 0.166;
  float oldSmearValue = 1.0;
  if (isMC_ && event->genJetPF2PATPT[index] > -990.) oldSmearValue = std::max(0.0, event->jetPF2PATPtRaw[index] + (event->jetPF2PATPtRaw[index] - event->genJetPF2PATPT[index]) * oldSmearCorr)/event->jetPF2PATPtRaw[index];
  float newJECCorr = 0.;
  if (fabs(event->jetPF2PATEta[index]) <= 0.8) {
    newJECCorr = 1.061;
    if (syst == 16) newJECCorr = 1.084;
    else if (syst == 32) newJECCorr = 1.038;
  }
  else if (fabs(event->jetPF2PATEta[index]) <= 1.3) {
    newJECCorr = 1.088;
    if (syst == 16) newJECCorr = 1.117;
    else if (syst == 32) newJECCorr = 1.059;
  }
  else if (fabs(event->jetPF2PATEta[index]) <= 1.9) {
    newJECCorr = 1.106;
    if (syst == 16) newJECCorr = 1.136;
    else if (syst == 32) newJECCorr = 1.076;
  }
  else if (fabs(event->jetPF2PATEta[index]) <= 2.5) {
    newJECCorr = 1.126;
    if (syst == 16) newJECCorr = 1.22;
    else if (syst == 32) newJECCorr = 1.032;
  }
  else if (fabs(event->jetPF2PATEta[index]) <= 3.0){
    newJECCorr = 1.343;
    if (syst == 16) newJECCorr = 1.466;
    else if (syst == 32) newJECCorr = 1.22;
  }
  else if (fabs(event->jetPF2PATEta[index]) <= 3.2){
    newJECCorr = 1.303;
    if (syst == 16) newJECCorr = 1.414;
    else if (syst == 32) newJECCorr = 1.192;
  } 
  else {
    newJECCorr = 1.320;
    if (syst == 16) newJECCorr = 1.606;
    else if (syst == 32) newJECCorr = 1.304;
  }
  float newSmearValue = 1.0;
  if (isMC_ && event->genJetPF2PATPT[index] > -990.) newSmearValue = std::max(0.0, event->jetPF2PATPtRaw[index] + (event->jetPF2PATPtRaw[index] - event->genJetPF2PATPT[index]) * newJECCorr)/event->jetPF2PATPtRaw[index];
  //  std::cout << event->jetPF2PATPtRaw[index] << std::setprecision(7) << " " << oldSmearValue << " " << newSmearValue << std::endl;
  TLorentzVector returnJet;
  if (newSmearValue < 0.01) {
    returnJet.SetPxPyPzE(0.01,0.01,0.01,0.01);
    return returnJet;
  }
  else returnJet.SetPxPyPzE(newSmearValue*event->jetPF2PATPx[index]/oldSmearValue,newSmearValue*event->jetPF2PATPy[index]/oldSmearValue,newSmearValue*event->jetPF2PATPz[index]/oldSmearValue,newSmearValue*event->jetPF2PATE[index]/oldSmearValue);
  if (isMC_){
    float jerUncer = getJECUncertainty(returnJet.Pt(),returnJet.Eta(),syst);
    returnJet *= 1+jerUncer;
  }
  //  if (returnJet.Pt() == 0.) std::cout << event->jetPF2PATPtRaw[index]  << " " << event->genJetPF2PATPT[index] << " " << event->jetPF2PATPx[index] << " " << event->jetPF2PATPy[index] << " " << event->jetPF2PATPz[index] << " " << event->jetPF2PATE[index] << " " << oldSmearValue << " " << newSmearValue <<std::endl;
  return returnJet;

}

void Cuts::getBWeight(AnalysisEvent* event, TLorentzVector jet, int index, float * mcTag, float * mcNoTag, float * dataTag, float * dataNoTag, float * errTag, float * errNoTag, float * err1, float * err2, float * err3, float * err4){
  //Use b-tagging efficiencies and scale factors.
  //Firstly get efficiency for pt/eta bin here.
  float eff = 1.;
  int partonFlavour = std::abs(event->jetPF2PATPID[index]);
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
  float SF = 1.;
  float SFerr = 0.;
  
  float x = jet.Pt();
  //Do some things if it's a b or c
  if (partonFlavour == 4 || partonFlavour == 5){
    SF = 1.00572*((1.+(0.013676*x))/(1.+(0.0143279*x)));
    int ptBin = getptbin_for_btag(jet.Pt());
    SFerr = SFb_error[ptBin];
    if (partonFlavour == 4) SFerr*=2;
  }
  //Light jets
  else {
    float min;
    float max;
    if (std::abs(jet.Eta()) < 0.5){
      SF = ((1.01177+(0.0023066*x))+(-4.56052e-06*(x*x)))+(2.57917e-09*(x*(x*x)));
      min = ((0.977761+(0.00170704*x))+(-3.2197e-06*(x*x)))+(1.78139e-09*(x*(x*x)));
      max = ((1.04582+(0.00290226*x))+(-5.89124e-06*(x*x)))+(3.37128e-09*(x*(x*x)));
    }
    else if (std::abs(jet.Eta()) < 1.0){
      SF = ((0.975966+(0.00196354*x))+(-3.83768e-06*(x*x)))+(2.17466e-09*(x*(x*x)));
      min = ((0.945135+(0.00146006*x))+(-2.70048e-06*(x*x)))+(1.4883e-09*(x*(x*x)));
      max = ((1.00683+(0.00246404*x))+(-4.96729e-06*(x*x)))+(2.85697e-09*(x*(x*x)));
    }
    else if (std::abs(jet.Eta()) < 1.5){
      SF = ((0.93821+(0.00180935*x))+(-3.86937e-06*(x*x)))+(2.43222e-09*(x*(x*x)));
      min = ((0.911657+(0.00142008*x))+(-2.87569e-06*(x*x)))+(1.76619e-09*(x*(x*x)));
      max = ((0.964787+(0.00219574*x))+(-4.85552e-06*(x*x)))+(3.09457e-09*(x*(x*x)));
    }
    else{
      SF = ((1.00022+(0.0010998*x))+(-3.10672e-06*(x*x)))+(2.35006e-09*(x*(x*x)));
      min = ((0.970045+(0.000862284*x))+(-2.31714e-06*(x*x)))+(1.68866e-09*(x*(x*x)));
      max = ((1.03039+(0.0013358*x))+(-3.89284e-06*(x*x)))+(3.01155e-09*(x*(x*x)));
    }
    SFerr = std::abs(max-SF)>fabs(min-SF)? std::abs(max-SF):std::abs(min-SF);
  }


  //Apply the weight of the jet and set the error
  if (event->jetPF2PATBDiscriminator[index] > bDiscCut_){
    *mcTag *= eff;
    *dataTag *= eff*SF;

    if (partonFlavour == 5 || partonFlavour == 4) *err1 += SFerr/SF;
    else *err3 += SFerr/SF;
  }
  else{
    *mcNoTag *= (1-eff);
    *dataNoTag *= (1-eff*SF);

    if (partonFlavour == 5 || partonFlavour == 4) *err2 += (-eff*SFerr)/(1-eff*SF);
    else *err4 += (-eff*SFerr)/(1-eff*SF);
    
  }

}

int Cuts::getptbin_for_btag(float pt){
  if(pt<30) return 0;
  else if(pt<40) return 1;
  else if(pt<50) return 2;
  else if(pt<60) return 3;
  else if(pt<70) return 4;
  else if(pt<80) return 5;
  else if(pt<100) return 6;
  else if(pt<120) return 7;
  else if(pt<160) return 8;
  else if(pt<210) return 9;
  else if(pt<260) return 10;
  else if(pt<320) return 11;
  else if(pt<400) return 12;
  else if(pt<500) return 13;
  else if(pt<600) return 14;
  else return 15;
  
}
