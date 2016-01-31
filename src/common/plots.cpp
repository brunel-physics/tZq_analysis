#include "plots.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>
#include "TLorentzVector.h"

Plots::Plots(std::vector<std::string> names, std::vector<float> xMins, std::vector<float> xMaxs, std::vector<int> nBins, std::vector<std::string> fillExps, std::vector<std::string>  xAxisLabels, std::vector<int> cutStage, unsigned int thisCutStage, std::string postfixName, const bool trileptonChannel):

  trileptonChannel_(trileptonChannel)

{
  //Get the function pointer map for later custopmisation. This is gonna be great, I promise.
  std::map<std::string, float (Plots::*)(AnalysisEvent*)> functionPointerMap = getFncPtrMap();
  plotPoint = std::vector<plot>(names.size());
  for (unsigned int i = 0; i < names.size(); i++){
    std::string plotName = names[i] + "_" + postfixName;
    plotPoint[i].name = plotName;
    plotPoint[i].plotHist = new TH1F(plotName.c_str(),plotName.c_str(),nBins[i],xMins[i],xMaxs[i]);
    plotPoint[i].fillExp = functionPointerMap[fillExps[i]];
    plotPoint[i].xAxisLabel = xAxisLabels[i];
    plotPoint[i].fillPlot = (unsigned int)cutStage[i] <= thisCutStage;
  }
  
}

Plots::~Plots(){
  for (unsigned int i = 0; i < plotPoint.size(); i++){
    delete plotPoint[i].plotHist;
  }
}

std::map<std::string, float (Plots::*)(AnalysisEvent*)> Plots::getFncPtrMap(){
  std::map<std::string, float (Plots::*)(AnalysisEvent*)> functionPointerMap;
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1Pt",&Plots::fillLepton1Pt));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1Eta",&Plots::fillLepton1Eta));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2Pt",&Plots::fillLepton2Pt));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2Eta",&Plots::fillLepton2Eta));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1RelIso",&Plots::fillLepton1RelIso));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2RelIso",&Plots::fillLepton2RelIso));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1Phi",&Plots::fillLepton1Phi));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2Phi",&Plots::fillLepton2Phi));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1MVA",&Plots::fillLepton1MVA));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2MVA",&Plots::fillLepton2MVA));

  if (trileptonChannel_){
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3Pt",&Plots::fillLepton3Pt));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3Eta",&Plots::fillLepton3Eta));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3Phi",&Plots::fillLepton3Phi));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3RelIso",&Plots::fillLepton3RelIso));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3MVA",&Plots::fillLepton3MVA));
  }
  else if (!trileptonChannel_){
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wQuark1Pt",&Plots::fillWbosonQuark1Pt));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wQuark1Eta",&Plots::fillWbosonQuark1Eta));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wQuark1Phi",&Plots::fillWbosonQuark1Phi));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wQuark2Pt",&Plots::fillWbosonQuark2Pt));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wQuark2Eta",&Plots::fillWbosonQuark2Eta));
    functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wQuark2Phi",&Plots::fillWbosonQuark2Phi));
  }

  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("met",&Plots::fillMetPlot));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("numbJets",&Plots::getNumberOfJets));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("leadingJetPt",&Plots::fillLeadingJetPt));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("leadingJetEta",&Plots::fillLeadingJetEta));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("leadingJetPhi",&Plots::fillLeadingJetPhi));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("secondJetPt",&Plots::fillSecondJetPt));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("secondJetEta",&Plots::fillSecondJetEta));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("secondJetPhi",&Plots::fillSecondJetPhi));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("numbBJets",&Plots::numbBJets));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepton1Pt",&Plots::fillZLep1Pt));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepton1Eta",&Plots::fillZLep1Eta));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepton2Pt",&Plots::fillZLep2Pt));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepton2Eta",&Plots::fillZLep2Eta));
  if (trileptonChannel_)  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wLeptonPt",&Plots::fillWLepPt));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wLeptonEta",&Plots::fillWLepEta));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepton1RelIso",&Plots::fillZLep1RelIso));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepton2RelIso",&Plots::fillZLep2RelIso));
  if (trileptonChannel_)  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wLeptonRelIso",&Plots::fillWLepRelIso));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepton1Phi",&Plots::fillZLep1Phi));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepton2Phi",&Plots::fillZLep2Phi));
  if (trileptonChannel_)  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wLeptonPhi",&Plots::fillWLepPhi));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zPairMass",&Plots::fillZPairMass));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zPairPt",&Plots::fillZPairPt));
  if (trileptonChannel_)  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wPair1Mass",&Plots::fillWPair1Mass));
  if (trileptonChannel_)  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wPair2Mass",&Plots::fillWPair2Mass));
  if (trileptonChannel_)  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lepMass",&Plots::fillLeptonMass));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1D0",&Plots::fillLepton1D0));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2D0",&Plots::fillLepton2D0));
  if (trileptonChannel_)  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3D0",&Plots::fillLepton3D0));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1DBD0",&Plots::fillLepton1DBD0));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2DBD0",&Plots::fillLepton2DBD0));
  if (trileptonChannel_)  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3DBD0",&Plots::fillLepton3DBD0));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1BeamSpotCorrectedD0",&Plots::fillLepton1BeamSpotCorrectedD0));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2BeamSpotCorrectedD0",&Plots::fillLepton2BeamSpotCorrectedD0));
  if (trileptonChannel_)functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3BeamSpotCorrectedD0",&Plots::fillLepton3BeamSpotCorrectedD0));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep1InnerTrackD0",&Plots::fillLepton1InnerTrackD0));
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep2InnerTrackD0",&Plots::fillLepton2InnerTrackD0));
  if (trileptonChannel_)functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lep3InnerTrackD0",&Plots::fillLepton3InnerTrackD0));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("wTransverseMass",&Plots::fillwTransverseMass));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("jjDelR",&Plots::filljjDelR));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepDelR",&Plots::fillzLepDelR));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("zLepDelPhi",&Plots::fillzLepDelPhi));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lbDelR",&Plots::filllbDelR));  
  functionPointerMap.insert(std::pair<std::string, float(Plots::*)(AnalysisEvent*)>("lbDelPhi",&Plots::filllbDelPhi));  

  return functionPointerMap;
}

float Plots::fillLepton1Pt(AnalysisEvent* event){
  if (event->electronIndexTight.size() >1){
    return event->elePF2PATPT[event->electronIndexTight[0]];
  }
  else{
    return event->muonPF2PATPt[event->muonIndexTight[0]];
  }
  return -10;
}

float Plots::fillLepton1Eta(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return fabs(event->elePF2PATEta[event->electronIndexTight[0]]);
  else return  fabs(event->muonPF2PATEta[event->muonIndexTight[0]]);
  return -10;
}
float Plots::fillLepton2Pt(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return event->elePF2PATPT[event->electronIndexTight[1]];
  else return  event->muonPF2PATPt[1];
  return -10;
}

float Plots::fillLepton2Eta(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return fabs(event->elePF2PATEta[event->electronIndexTight[1]]);
  else return  fabs(event->muonPF2PATEta[event->muonIndexTight[1]]);
  return -10;
}
float Plots::fillLepton3Pt(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return event->elePF2PATPT[event->electronIndexTight[2]];
  if (event->muonIndexTight.size() > 2)
    return event->muonPF2PATPt[event->muonIndexTight[2]];
  if (event->electronIndexTight.size() > 1)
    return event->muonPF2PATPt[event->muonIndexTight[0]];
  else
    return event->elePF2PATPT[event->electronIndexTight[0]];
  
  return -10;
}

float Plots::fillLepton3Eta(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return fabs(event->elePF2PATEta[event->electronIndexTight[2]]);
  if (event->muonIndexTight.size() > 2)
    return event->muonPF2PATEta[event->muonIndexTight[2]];
  if (event->electronIndexTight.size() > 1)
    return event->muonPF2PATEta[event->muonIndexTight[0]];
  else
    return event->elePF2PATEta[event->electronIndexTight[0]];
  return -10;
}

float Plots::fillWbosonQuark1Pt(AnalysisEvent* event){
  return event->wPairQuarks.first.Pt();
}

float Plots::fillWbosonQuark1Eta(AnalysisEvent* event){
  return event->wPairQuarks.first.Eta();
}

float Plots::fillWbosonQuark1Phi(AnalysisEvent* event){
  return event->wPairQuarks.first.Phi();
}

float Plots::fillWbosonQuark2Pt(AnalysisEvent* event){
  return event->wPairQuarks.second.Pt();
}

float Plots::fillWbosonQuark2Eta(AnalysisEvent* event){
  return event->wPairQuarks.second.Eta();
}

float Plots::fillWbosonQuark2Phi(AnalysisEvent* event){
  return event->wPairQuarks.second.Phi();
}

float Plots::fillLepton1RelIso(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return event->elePF2PATComRelIsoRho[event->electronIndexTight[0]]/event->elePF2PATPT[event->electronIndexTight[0]];
  else return  event->muonPF2PATComRelIsodBeta[event->muonIndexTight[0]];
  return -10;
}

float Plots::fillLepton2RelIso(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return event->elePF2PATComRelIsoRho[event->electronIndexTight[1]]/event->elePF2PATPT[event->electronIndexTight[1]];
  else return event->muonPF2PATComRelIsodBeta[event->muonIndexTight[1]];
  return -10;
}

float Plots::fillLepton3RelIso(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return event->elePF2PATComRelIsoRho[event->electronIndexTight[2]]/event->elePF2PATPT[event->electronIndexTight[2]];
  if (event->muonIndexTight.size() > 2)
    return event->muonPF2PATComRelIsodBeta[event->muonIndexTight[2]];
  if (event->electronIndexTight.size() > 1)
    return event->muonPF2PATComRelIsodBeta[event->muonIndexTight[0]];
  else
    return event->elePF2PATComRelIsoRho[event->electronIndexTight[0]]/event->elePF2PATPT[event->electronIndexTight[0]];
  return -10;
}

float Plots::fillLepton1Phi(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return (event->elePF2PATPhi[event->electronIndexTight[0]]);
 else return  (event->muonPF2PATPhi[event->muonIndexTight[0]]);
  return -10;
}
float Plots::fillLepton2Phi(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return (event->elePF2PATPhi[event->electronIndexTight[1]]);
  else return  (event->muonPF2PATPhi[event->muonIndexTight[1]]);
  return -10;
}
float Plots::fillLepton3Phi(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return (event->elePF2PATPhi[event->electronIndexTight[2]]);
  if (event->muonIndexTight.size() > 2)
    return event->muonPF2PATPhi[event->muonIndexTight[2]];
  if (event->electronIndexTight.size() > 1)
    return event->muonPF2PATPhi[event->muonIndexTight[0]];
  else
    return event->elePF2PATPhi[event->electronIndexTight[0]];
  return -10;
}

float Plots::fillLepton1MVA(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return event->elePF2PATMVA[event->electronIndexTight[0]];
  else return -5;
  return -10;
}
float Plots::fillLepton2MVA(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return event->elePF2PATMVA[event->electronIndexTight[1]];
  else return -5;
  return -10;
}
float Plots::fillLepton3MVA(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return event->elePF2PATMVA[event->electronIndexTight[2]];
  if (event->muonIndexTight.size() > 2)
    return -5;
  if (event->electronIndexTight.size() > 1)
    return -5;
  else
    return event->elePF2PATMVA[event->electronIndexTight[0]];

  return -10;
}

float Plots::fillMetPlot(AnalysisEvent* event){
  return event->metPF2PATEt;
}

float Plots::getNumberOfJets(AnalysisEvent* event){
  return event->jetIndex.size();
}

float Plots::fillLeadingJetPt(AnalysisEvent* event){
  if (event->jetIndex.size() > 0) return (std::sqrt(event->jetPF2PATPx[event->jetIndex[0]] * event->jetPF2PATPx[event->jetIndex[0]] + event->jetPF2PATPy[event->jetIndex[0]] * event->jetPF2PATPy[event->jetIndex[0]]) );
  return -10;
}

float Plots::fillSecondJetPt(AnalysisEvent* event){
  if (event->jetIndex.size() > 1) return (std::sqrt(event->jetPF2PATPx[event->jetIndex[1]] * event->jetPF2PATPx[event->jetIndex[1]] + event->jetPF2PATPy[event->jetIndex[1]] * event->jetPF2PATPy[event->jetIndex[1]]) );
  return -10;
}

float Plots::fillLeadingJetEta(AnalysisEvent* event){
  if (event->jetIndex.size() > 0) return fabs(event->jetPF2PATEta[event->jetIndex[0]]);
  return -10;
}

float Plots::fillSecondJetEta(AnalysisEvent* event){
  if (event->jetIndex.size() > 1) return fabs(event->jetPF2PATEta[event->jetIndex[1]]);
  return -10;
}

float Plots::fillLeadingJetPhi(AnalysisEvent* event){
  if (event->jetIndex.size() > 0) return event->jetPF2PATPhi[event->jetIndex[0]];
  return -10;
}

float Plots::fillSecondJetPhi(AnalysisEvent* event){
  if (event->jetIndex.size() > 1) return event->jetPF2PATPhi[event->jetIndex[1]];
  return -10;
}

float Plots::numbBJets(AnalysisEvent* event){
  return event->bTagIndex.size();
}

float Plots::fillZLep1Pt(AnalysisEvent* event){
  return event->zPairLeptons.first.Pt();
}

float Plots::fillZLep2Pt(AnalysisEvent* event){
  return event->zPairLeptons.second.Pt();
}

float Plots::fillZLep1Eta(AnalysisEvent* event){
  return fabs(event->zPairLeptons.first.Eta());
}

float Plots::fillZLep2Eta(AnalysisEvent* event){
  return fabs(event->zPairLeptons.second.Eta());
}

float Plots::fillZLep1Phi(AnalysisEvent* event){
  return event->zPairLeptons.first.Phi();
}

float Plots::fillZLep2Phi(AnalysisEvent* event){
  return event->zPairLeptons.second.Phi();
}

float Plots::fillZLep1RelIso(AnalysisEvent* event){
  return event->zPairRelIso.first;
}

float Plots::fillZLep2RelIso(AnalysisEvent* event){
  return event->zPairRelIso.second;
}

float Plots::fillWLepPt(AnalysisEvent* event){
  return event->wLepton.Pt();
}

float Plots::fillWLepEta(AnalysisEvent* event){
  return fabs(event->wLepton.Eta());
}

float Plots::fillWLepPhi(AnalysisEvent* event){
  return event->wLepton.Phi();
}

float Plots::fillWLepRelIso(AnalysisEvent* event){
  return event->wLeptonRelIso;
}

float Plots::fillZPairMass(AnalysisEvent* event){
  return (event->zPairLeptons.first + event->zPairLeptons.second).M();
}

float Plots::fillZPairPt(AnalysisEvent* event){
  return (event->zPairLeptons.first + event->zPairLeptons.second).Pt();
}

float Plots::fillWPair1Mass(AnalysisEvent* event){
  return (event->zPairLeptons.first + event->wLepton).M();
}

float Plots::fillWPair2Mass(AnalysisEvent* event){
  return (event->zPairLeptons.second + event->wLepton).M();
}

float Plots::fillLeptonMass(AnalysisEvent* event){
  return (event->zPairLeptons.second + event->zPairLeptons.first + event->wLepton).M();
}

float Plots::fillLepton1D0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return (event->elePF2PATTrackD0[event->electronIndexTight[0]]);
 else return  (event->muonPF2PATD0[event->muonIndexTight[0]]);
  return -10;
}
float Plots::fillLepton2D0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return (event->elePF2PATTrackD0[event->electronIndexTight[1]]);
  else return  (event->muonPF2PATD0[event->muonIndexTight[1]]);
  return -10;
}
float Plots::fillLepton3D0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return (event->elePF2PATTrackD0[event->electronIndexTight[2]]);
  if (event->muonIndexTight.size() > 2)
    return event->muonPF2PATD0[event->muonIndexTight[2]];
  if (event->electronIndexTight.size() > 1)
    return event->muonPF2PATD0[event->muonIndexTight[0]];
  else
    return event->elePF2PATTrackD0[event->electronIndexTight[0]];
  return -10;
}

float Plots::fillLepton1DBD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return (event->elePF2PATTrackDBD0[event->electronIndexTight[0]]);
 else return  (event->muonPF2PATTrackDBD0[event->muonIndexTight[0]]);
  return -10;
}
float Plots::fillLepton2DBD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return (event->elePF2PATTrackDBD0[event->electronIndexTight[1]]);
  else return  (event->muonPF2PATTrackDBD0[event->muonIndexTight[1]]);
  return -10;
}
float Plots::fillLepton3DBD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return (event->elePF2PATTrackDBD0[event->electronIndexTight[2]]);
  if (event->muonIndexTight.size() > 2)
    return event->muonPF2PATTrackDBD0[event->muonIndexTight[2]];
  if (event->electronIndexTight.size() > 1)
    return event->muonPF2PATTrackDBD0[event->muonIndexTight[0]];
  else
    return event->elePF2PATTrackDBD0[event->electronIndexTight[0]];
  return -10;
}

float Plots::fillLepton1BeamSpotCorrectedD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return (event->elePF2PATBeamSpotCorrectedTrackD0[event->electronIndexTight[0]]);
 else return  (event->muonPF2PATBeamSpotCorrectedD0[event->muonIndexTight[0]]);
  return -10;
}
float Plots::fillLepton2BeamSpotCorrectedD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return (event->elePF2PATBeamSpotCorrectedTrackD0[event->electronIndexTight[1]]);
  else return  (event->muonPF2PATBeamSpotCorrectedD0[event->muonIndexTight[1]]);
  return -10;
}
float Plots::fillLepton3BeamSpotCorrectedD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return (event->elePF2PATBeamSpotCorrectedTrackD0[event->electronIndexTight[2]]);
  if (event->muonIndexTight.size() > 2)
    return event->muonPF2PATBeamSpotCorrectedD0[event->muonIndexTight[2]];
  if (event->electronIndexTight.size() > 1)
    return event->muonPF2PATBeamSpotCorrectedD0[event->muonIndexTight[0]];
  else
    return event->elePF2PATBeamSpotCorrectedTrackD0[event->electronIndexTight[0]];
  return -10;
}

float Plots::fillLepton1InnerTrackD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return -10;
 else return  (event->muonPF2PATDBInnerTrackD0[event->muonIndexTight[0]]);
  return -10;
}
float Plots::fillLepton2InnerTrackD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 1)
    return -10;
  else return  (event->muonPF2PATDBInnerTrackD0[event->muonIndexTight[1]]);
  return -10;
}
float Plots::fillLepton3InnerTrackD0(AnalysisEvent* event){
  if (event->electronIndexTight.size() > 2)
    return -10;
  if (event->muonIndexTight.size() > 2)
    return event->muonPF2PATDBInnerTrackD0[event->muonIndexTight[2]];
  if (event->electronIndexTight.size() > 1)
    return event->muonPF2PATDBInnerTrackD0[event->muonIndexTight[0]];
  else
    return -10;
  return -10;
}

float Plots::fillwTransverseMass(AnalysisEvent* event){
  if ( trileptonChannel_ ) {
  TLorentzVector tempMet;
  tempMet.SetPtEtaPhiE(event->metPF2PATPt,0,event->metPF2PATPhi,event->metPF2PATEt);
  return std::sqrt(2*event->metPF2PATPt*event->wLepton.Pt()*(1-cos(event->metPF2PATPhi - event->wLepton.Phi())));
  }
  else if ( !trileptonChannel_ ) {
  return std::sqrt(2*event->jetPF2PATPt[event->wPairIndex.first]*event->jetPF2PATPt[event->wPairIndex.second]*(1-cos(event->jetPF2PATPhi[event->wPairIndex.first] - event->jetPF2PATPhi[event->wPairIndex.second])));
  }

}

float Plots::filljjDelR(AnalysisEvent* event){
  TLorentzVector tempJet1;
  TLorentzVector tempJet2;
  if (event->jetIndex.size() < 2) return -1.;
  tempJet1.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[0]],event->jetPF2PATPy[event->jetIndex[0]],event->jetPF2PATPz[event->jetIndex[0]],event->jetPF2PATE[event->jetIndex[0]]);
  tempJet2.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[1]],event->jetPF2PATPy[event->jetIndex[1]],event->jetPF2PATPz[event->jetIndex[1]],event->jetPF2PATE[event->jetIndex[1]]);
  return tempJet1.DeltaR(tempJet2);
}

float Plots::fillzLepDelR(AnalysisEvent* event){
  return (event->zPairLeptons.first.DeltaR(event->zPairLeptons.second));
}

float Plots::fillzLepDelPhi(AnalysisEvent* event){
  return (event->zPairLeptons.first.DeltaPhi(event->zPairLeptons.second));
}

float Plots::filllbDelR(AnalysisEvent* event){
  TLorentzVector tempJet1;
  if (event->bTagIndex.size() < 1) return -10.;
  tempJet1.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
  return tempJet1.DeltaR(event->wLepton);

}

float Plots::filllbDelPhi(AnalysisEvent* event){
  TLorentzVector tempJet1;
  if (event->bTagIndex.size() < 1) return -10.;
  tempJet1.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
  return tempJet1.DeltaPhi(event->wLepton);
}

void Plots::fillAllPlots(AnalysisEvent* event, float eventWeight){
  for (unsigned int i = 0; i < plotPoint.size(); i++){
    if (plotPoint[i].fillPlot){
      plotPoint[i].plotHist->Fill((this->*plotPoint[i].fillExp)(event),eventWeight);
    }
  }
}

void Plots::saveAllPlots(){
  for (unsigned int i = 0; i < plotPoint.size(); i++){
    plotPoint[i].plotHist->SaveAs(("plots/"+plotPoint[i].name + ".root").c_str());
    plotPoint[i].plotHist->SaveAs(("plots/"+plotPoint[i].name + ".png").c_str());
  }
}
