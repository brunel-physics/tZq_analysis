#include "plots.hpp"

#include "TH1D.h"
#include "TLorentzVector.h"

#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

Plots::Plots(const std::vector<std::string> titles,
             const std::vector<std::string> names,
             const std::vector<float> xMins,
             const std::vector<float> xMaxs,
             const std::vector<int> nBins,
             const std::vector<std::string> fillExps,
             const std::vector<std::string> xAxisLabels,
             const std::vector<int> cutStage,
             const unsigned thisCutStage,
             const std::string postfixName)
{
    // Get the function pointer map for later custopmisation. This is gonna be
    // great, I promise.
    const auto functionMap{getFncMap()};

    plotPoint = std::vector<plot>(names.size());
    for (unsigned i{0}; i < names.size(); i++)
    {
        std::string plotName = names[i] + "_" + postfixName;
        plotPoint[i].name = plotName;
        plotPoint[i].title = titles[i];
        plotPoint[i].fillExp = functionMap.at(fillExps[i]);
        plotPoint[i].xAxisLabel = xAxisLabels[i];
        plotPoint[i].plotHist =
            new TH1D{plotName.c_str(),
                     (plotName + ";" + plotPoint[i].xAxisLabel).c_str(),
                     nBins[i],
                     xMins[i],
                     xMaxs[i]};
        plotPoint[i].fillPlot =
            boost::numeric_cast<unsigned>(cutStage[i]) <= thisCutStage;
    }
}

Plots::~Plots()
{
    for (unsigned i{0}; i < plotPoint.size(); i++)
    {
        delete plotPoint[i].plotHist;
    }
}

std::unordered_map<std::string, std::function<float(const AnalysisEvent&)>>
    Plots::getFncMap() const
{
    return {
        {"lep1Pt",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 TLorentzVector tempVec{
                     event.elePF2PATPX[event.electronIndexTight[0]],
                     event.elePF2PATPY[event.electronIndexTight[0]],
                     event.elePF2PATPZ[event.electronIndexTight[0]],
                     event.elePF2PATE[event.electronIndexTight[0]]};
                 return tempVec.Pt();
             }
             else
             {
                 TLorentzVector tempVec{
                     event.muonPF2PATPX[event.muonIndexTight[0]],
                     event.muonPF2PATPY[event.muonIndexTight[0]],
                     event.muonPF2PATPZ[event.muonIndexTight[0]],
                     event.muonPF2PATE[event.muonIndexTight[0]]};
                 tempVec *= event.muonMomentumSF[0];
                 return tempVec.Pt();
             }
         }},
        {"lep1Eta",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return std::abs(
                     event.elePF2PATSCEta[event.electronIndexTight[0]]);
             }
             else
             {
                 TLorentzVector tempVec{
                     event.muonPF2PATPX[event.muonIndexTight[0]],
                     event.muonPF2PATPY[event.muonIndexTight[0]],
                     event.muonPF2PATPZ[event.muonIndexTight[0]],
                     event.muonPF2PATE[event.muonIndexTight[0]]};
                 tempVec *= event.muonMomentumSF[0];
                 return tempVec.Eta();
             }
         }},
        {"lep2Pt",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 TLorentzVector tempVec{
                     event.elePF2PATPX[event.electronIndexTight[1]],
                     event.elePF2PATPY[event.electronIndexTight[1]],
                     event.elePF2PATPZ[event.electronIndexTight[1]],
                     event.elePF2PATE[event.electronIndexTight[1]]};
                 return tempVec.Pt();
             }
             else
             {
                 TLorentzVector tempVec{
                     event.muonPF2PATPX[event.muonIndexTight[1]],
                     event.muonPF2PATPY[event.muonIndexTight[1]],
                     event.muonPF2PATPZ[event.muonIndexTight[1]],
                     event.muonPF2PATE[event.muonIndexTight[1]]};
                 tempVec *= event.muonMomentumSF[1];
                 return tempVec.Pt();
             }
         }},
        {"lep2Eta",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return std::abs(
                     event.elePF2PATSCEta[event.electronIndexTight[1]]);
             }
             else
             {
                 TLorentzVector tempVec{
                     event.muonPF2PATPX[event.muonIndexTight[1]],
                     event.muonPF2PATPY[event.muonIndexTight[1]],
                     event.muonPF2PATPZ[event.muonIndexTight[1]],
                     event.muonPF2PATE[event.muonIndexTight[1]]};
                 tempVec *= event.muonMomentumSF[1];
                 return tempVec.Pt();
             }
         }},
        {"lep1RelIso",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return event
                     .elePF2PATComRelIsoRho[event.electronIndexTight[0]];
             }
             else
             {
                 return event
                     .muonPF2PATComRelIsodBeta[event.muonIndexTight[0]];
             }
         }},
        {"lep2RelIso",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return event
                     .elePF2PATComRelIsoRho[event.electronIndexTight[1]];
             }
             else
             {
                 return event
                     .muonPF2PATComRelIsodBeta[event.muonIndexTight[1]];
             }
         }},
        {"lep1Phi",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return (event.elePF2PATPhi[event.electronIndexTight[0]]);
             }
             else
             {
                 TLorentzVector tempVec{
                     event.muonPF2PATPX[event.muonIndexTight[0]],
                     event.muonPF2PATPY[event.muonIndexTight[0]],
                     event.muonPF2PATPZ[event.muonIndexTight[0]],
                     event.muonPF2PATE[event.muonIndexTight[0]]};
                 tempVec *= event.muonMomentumSF[0];
                 return tempVec.Phi();
             }
         }},
        {"lep2Phi",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return (event.elePF2PATPhi[event.electronIndexTight[1]]);
             }
             else
             {
                 TLorentzVector tempVec{
                     event.muonPF2PATPX[event.muonIndexTight[1]],
                     event.muonPF2PATPY[event.muonIndexTight[1]],
                     event.muonPF2PATPZ[event.muonIndexTight[1]],
                     event.muonPF2PATE[event.muonIndexTight[1]]};
                 tempVec *= event.muonMomentumSF[1];
                 return tempVec.Phi();
             }
         }},
        {"wQuark1Pt",
         [](const AnalysisEvent& event) -> float {
             return event.wPairQuarks.first.Pt();
         }},
        {"wQuark1Eta",
         [](const AnalysisEvent& event) -> float {
             return event.wPairQuarks.first.Eta();
         }},
        {"wQuark1Phi",
         [](const AnalysisEvent& event) -> float {
             return event.wPairQuarks.first.Phi();
         }},
        {"wQuark2Pt",
         [](const AnalysisEvent& event) -> float {
             return event.wPairQuarks.second.Pt();
         }},
        {"wQuark2Eta",
         [](const AnalysisEvent& event) -> float {
             return std::abs(event.wPairQuarks.second.Eta());
         }},
        {"wQuark2Phi",
         [](const AnalysisEvent& event) -> float {
             return event.wPairQuarks.second.Phi();
         }},
        {"met",
         [](const AnalysisEvent& event) -> float { return event.metPF2PATEt; }},
        {"numbJets",
         [](const AnalysisEvent& event) -> float { return event.jetIndex.size(); }},
        {"totalJetMass",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector totalJet;
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt = event.jetIndex.begin();
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue = event.jetSmearValue[*jetIt];
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totalJet += tempJet;
                 }
                 return totalJet.M();
             }
             else
             {
                 return -10;
             }
         }},
        {"totalJetPt",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector totalJet;
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt = event.jetIndex.begin();
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue = event.jetSmearValue[*jetIt];
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totalJet += tempJet;
                 }
                 return totalJet.Pt();
             }
             else
             {
                 return -10;
             }
         }},
        {"totalJetEta",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector totalJet;
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt = event.jetIndex.begin();
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue = event.jetSmearValue[*jetIt];
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totalJet += tempJet;
                 }
                 return totalJet.Eta();
             }
             else
             {
                 return -10;
             }
         }},
        {"totalJetPhi",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector totalJet;
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt = event.jetIndex.begin();
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue = event.jetSmearValue[*jetIt];
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totalJet += tempJet;
                 }
                 return totalJet.Phi();
             }
             else
             {
                 return -10;
             }
         }},
        {"leadingJetPt",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 0)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[0]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[0]],
                                    event.jetPF2PATPy[event.jetIndex[0]],
                                    event.jetPF2PATPz[event.jetIndex[0]],
                                    event.jetPF2PATE[event.jetIndex[0]]);
                 tempJet *= smearValue;
                 return tempJet.Pt();
             }
             else
             {
                 return -10;
             }
         }},
        {"leadingJetEta",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 0)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[0]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[0]],
                                    event.jetPF2PATPy[event.jetIndex[0]],
                                    event.jetPF2PATPz[event.jetIndex[0]],
                                    event.jetPF2PATE[event.jetIndex[0]]);
                 tempJet *= smearValue;
                 return tempJet.Eta();
             }
             else
             {
                 return -10;
             }
         }},
        {"leadingJetPhi",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 0)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[0]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[0]],
                                    event.jetPF2PATPy[event.jetIndex[0]],
                                    event.jetPF2PATPz[event.jetIndex[0]],
                                    event.jetPF2PATE[event.jetIndex[0]]);
                 tempJet *= smearValue;
                 return tempJet.Phi();
             }
             else
             {
                 return -10;
             }
         }},
        {"leadingJetBDisc",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 0)
             {
                 return event
                     .jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                         [event.jetIndex[0]];
             }
             else
             {
                 return -10;
             }
         }},
        {"secondJetPt",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 1)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[1]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[1]],
                                    event.jetPF2PATPy[event.jetIndex[1]],
                                    event.jetPF2PATPz[event.jetIndex[1]],
                                    event.jetPF2PATE[event.jetIndex[1]]);
                 tempJet *= smearValue;
                 return tempJet.Pt();
             }
             else
             {
                 return -10;
             }
         }},
        {"secondJetEta",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 1)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[1]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[1]],
                                    event.jetPF2PATPy[event.jetIndex[1]],
                                    event.jetPF2PATPz[event.jetIndex[1]],
                                    event.jetPF2PATE[event.jetIndex[1]]);
                 tempJet *= smearValue;
                 return tempJet.Eta();
             }
             else
             {
                 return -10;
             }
         }},
        {"secondJetPhi",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 1)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[1]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[1]],
                                    event.jetPF2PATPy[event.jetIndex[1]],
                                    event.jetPF2PATPz[event.jetIndex[1]],
                                    event.jetPF2PATE[event.jetIndex[1]]);
                 tempJet *= smearValue;
                 return tempJet.Phi();
             }
             else
             {
                 return -10;
             }
         }},
        {"secondJetBDisc",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 1)
             {
                 return event
                     .jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                         [event.jetIndex[1]];
             }
             else
             {
                 return -10;
             }
         }},
        {"thirdJetPt",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 2)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[2]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[2]],
                                    event.jetPF2PATPy[event.jetIndex[2]],
                                    event.jetPF2PATPz[event.jetIndex[2]],
                                    event.jetPF2PATE[event.jetIndex[2]]);
                 tempJet *= smearValue;
                 return tempJet.Pt();
             }
             else
             {
                 return -10;
             }
         }},
        {"thirdJetEta",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 2)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[2]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[2]],
                                    event.jetPF2PATPy[event.jetIndex[2]],
                                    event.jetPF2PATPz[event.jetIndex[2]],
                                    event.jetPF2PATE[event.jetIndex[2]]);
                 tempJet *= smearValue;
                 return tempJet.Eta();
             }
             else
             {
                 return -10;
             }
         }},
        {"thirdJetPhi",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 2)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[2]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[2]],
                                    event.jetPF2PATPy[event.jetIndex[2]],
                                    event.jetPF2PATPz[event.jetIndex[2]],
                                    event.jetPF2PATE[event.jetIndex[2]]);
                 tempJet *= smearValue;
                 return tempJet.Phi();
             }
             else
             {
                 return -10;
             }
         }},
        {"thirdJetBDisc",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 2)
             {
                 return event
                     .jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                         [event.jetIndex[2]];
             }
             else
             {
                 return -10;
             }
         }},
        {"fourthJetPt",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 3)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[3]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[3]],
                                    event.jetPF2PATPy[event.jetIndex[3]],
                                    event.jetPF2PATPz[event.jetIndex[3]],
                                    event.jetPF2PATE[event.jetIndex[3]]);
                 tempJet *= smearValue;
                 return tempJet.Pt();
             }
             else
             {
                 return -10;
             }
         }},
        {"fourthJetEta",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 3)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[3]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[3]],
                                    event.jetPF2PATPy[event.jetIndex[3]],
                                    event.jetPF2PATPz[event.jetIndex[3]],
                                    event.jetPF2PATE[event.jetIndex[3]]);
                 tempJet *= smearValue;
                 return tempJet.Eta();
             }
             else
             {
                 return -10;
             }
         }},
        {"fourthJetPhi",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 3)
             {
                 TLorentzVector tempJet;
                 float smearValue = event.jetSmearValue[event.jetIndex[3]];
                 tempJet.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[3]],
                                    event.jetPF2PATPy[event.jetIndex[3]],
                                    event.jetPF2PATPz[event.jetIndex[3]],
                                    event.jetPF2PATE[event.jetIndex[3]]);
                 tempJet *= smearValue;
                 return tempJet.Phi();
             }
             else
             {
                 return -10;
             }
         }},
        {"fourthJetBDisc",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() > 3)
             {
                 return event
                     .jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                         [event.jetIndex[3]];
             }
             else
             {
                 return -10;
             }
         }},
        {"numbBJets",
         [](const AnalysisEvent& event) -> float { return event.bTagIndex.size(); }},
        {"bTagDisc",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() > 0)
             {
                 return event
                     .jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                         [event.jetIndex[event.bTagIndex[0]]];
             }
             return -10;
         }},
        {"zLepton1Pt",
         [](const AnalysisEvent& event) -> float {
             return event.zPairLeptons.first.Pt();
         }},
        {"zLepton1Eta",
         [](const AnalysisEvent& event) -> float {
             return std::abs(event.zPairLeptons.first.Eta());
         }},
        {"zLepton2Pt",
         [](const AnalysisEvent& event) -> float {
             return event.zPairLeptons.second.Pt();
         }},
        {"zLepton2Eta",
         [](const AnalysisEvent& event) -> float {
             return std::abs(event.zPairLeptons.second.Eta());
         }},
        {"zLepton1RelIso",
         [](const AnalysisEvent& event) -> float {
             return event.zPairRelIso.first;
         }},
        {"zLepton2RelIso",
         [](const AnalysisEvent& event) -> float {
             return event.zPairRelIso.second;
         }},
        {"zLepton1Phi",
         [](const AnalysisEvent& event) -> float {
             return event.zPairLeptons.first.Phi();
         }},
        {"zLepton2Phi",
         [](const AnalysisEvent& event) -> float {
             return event.zPairLeptons.second.Phi();
         }},
        {"zPairMass",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .M();
         }},
        {"zPairPt",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .Pt();
         }},
        {"zPairEta",
         [](const AnalysisEvent& event) -> float {
             return std::abs(
                 (event.zPairLeptons.first + event.zPairLeptons.second)
                     .Eta());
         }},
        {"zPairPhi",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .Phi();
         }},
        {"wPairMass",
         [](const AnalysisEvent& event) -> float {
             return (event.wPairQuarks.first + event.wPairQuarks.second).M();
         }},
        {"topMass",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() > 0)
             {
                 TLorentzVector tempBjet;
                 float smearValue = event.jetSmearValue[event.bTagIndex[0]];
                 tempBjet.SetPtEtaPhiE(
                     event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
                 tempBjet *= smearValue;
                 return (tempBjet + event.wPairQuarks.first
                         + event.wPairQuarks.second)
                     .M();
             }
             else
             {
                 return -10;
             }
         }},
        {"topPt",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() > 0)
             {
                 TLorentzVector tempBjet;
                 float smearValue = event.jetSmearValue[event.bTagIndex[0]];
                 tempBjet.SetPtEtaPhiE(
                     event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
                 tempBjet *= smearValue;
                 return (tempBjet + event.wPairQuarks.first
                         + event.wPairQuarks.second)
                     .Pt();
             }
             else
             {
                 return -10;
             }
         }},
        {"topEta",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() > 0)
             {
                 TLorentzVector tempBjet;
                 float smearValue = event.jetSmearValue[event.bTagIndex[0]];
                 tempBjet.SetPtEtaPhiE(
                     event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
                 tempBjet *= smearValue;
                 return std::abs((tempBjet + event.wPairQuarks.first
                                  + event.wPairQuarks.second)
                                     .Eta());
             }
             else
             {
                 return -10;
             }
         }},
        {"topPhi",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() > 0)
             {
                 TLorentzVector tempBjet;
                 float smearValue = event.jetSmearValue[event.bTagIndex[0]];
                 tempBjet.SetPtEtaPhiE(
                     event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                     event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
                 tempBjet *= smearValue;
                 return (tempBjet + event.wPairQuarks.first
                         + event.wPairQuarks.second)
                     .Phi();
             }
             else
             {
                 return -10;
             }
         }},
        {"lep1D0",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return (event.elePF2PATD0PV[event.electronIndexTight[0]]);
             }
             else
             {
                 return (event.muonPF2PATDBPV[event.muonIndexTight[0]]);
             }
         }},
        {"lep2D0",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return (event.elePF2PATD0PV[event.electronIndexTight[1]]);
             }
             else
             {
                 return (event.muonPF2PATDBPV[event.muonIndexTight[1]]);
             }
         }},
        {"lep1DBD0",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return (
                     event.elePF2PATTrackDBD0[event.electronIndexTight[0]]);
             }
             else
             {
                 return (event.muonPF2PATTrackDBD0[event.muonIndexTight[0]]);
             }
         }},
        {"lep2DBD0",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return (
                     event.elePF2PATTrackDBD0[event.electronIndexTight[1]]);
             }
             else
             {
                 return (event.muonPF2PATTrackDBD0[event.muonIndexTight[1]]);
             }
         }},
        {"lep1BeamSpotCorrectedD0",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return (event.elePF2PATBeamSpotCorrectedTrackD0
                             [event.electronIndexTight[0]]);
             }
             else
             {
                 return (event.muonPF2PATBeamSpotCorrectedD0
                             [event.muonIndexTight[0]]);
             }
         }},
        {"lep2BeamSpotCorrectedD0",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return (event.elePF2PATBeamSpotCorrectedTrackD0
                             [event.electronIndexTight[1]]);
             }
             else
             {
                 return (event.muonPF2PATBeamSpotCorrectedD0
                             [event.muonIndexTight[1]]);
             }
         }},
        {"lep1InnerTrackD0",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return -10;
             }
             else
             {
                 return (
                     event.muonPF2PATDBInnerTrackD0[event.muonIndexTight[0]]);
             }
         }},
        {"lep2InnerTrackD0",
         [](const AnalysisEvent& event) -> float {
             if (event.electronIndexTight.size() > 1)
             {
                 return -10;
             }
             else
             {
                 return (
                     event.muonPF2PATDBInnerTrackD0[event.muonIndexTight[1]]);
             }
         }},
        {"wTransverseMass",
         [](const AnalysisEvent& event) -> float {
             return std::sqrt(
                 2 * event.wPairQuarks.first.Pt()
                 * event.wPairQuarks.second.Pt()
                 * (1
                    - std::cos(event.wPairQuarks.first.Phi()
                               - event.wPairQuarks.second.Phi())));
         }},
        {"jjDelR",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempJet1;
             TLorentzVector tempJet2;
             if (event.jetIndex.size() < 2)
             {
                 return -1.;
             }
             float smearValue1 = event.jetSmearValue[event.jetIndex[0]];
             float smearValue2 = event.jetSmearValue[event.jetIndex[1]];
             tempJet1.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[0]],
                                 event.jetPF2PATPy[event.jetIndex[0]],
                                 event.jetPF2PATPz[event.jetIndex[0]],
                                 event.jetPF2PATE[event.jetIndex[0]]);
             tempJet2.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[1]],
                                 event.jetPF2PATPy[event.jetIndex[1]],
                                 event.jetPF2PATPz[event.jetIndex[1]],
                                 event.jetPF2PATE[event.jetIndex[1]]);
             tempJet1 *= smearValue1;
             tempJet2 *= smearValue2;
             return tempJet1.DeltaR(tempJet2);
         }},
        {"jjDelPhi",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempJet1;
             TLorentzVector tempJet2;
             if (event.jetIndex.size() < 2)
             {
                 return -1.;
             }
             float smearValue1 = event.jetSmearValue[event.jetIndex[0]];
             float smearValue2 = event.jetSmearValue[event.jetIndex[1]];
             tempJet1.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[0]],
                                 event.jetPF2PATPy[event.jetIndex[0]],
                                 event.jetPF2PATPz[event.jetIndex[0]],
                                 event.jetPF2PATE[event.jetIndex[0]]);
             tempJet2.SetPxPyPzE(event.jetPF2PATPx[event.jetIndex[1]],
                                 event.jetPF2PATPy[event.jetIndex[1]],
                                 event.jetPF2PATPz[event.jetIndex[1]],
                                 event.jetPF2PATE[event.jetIndex[1]]);
             tempJet1 *= smearValue1;
             tempJet2 *= smearValue2;
             return tempJet1.DeltaPhi(tempJet2);
         }},
        {"wwDelR",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() < 3)
             {
                 return -1.;
             }
             return event.wPairQuarks.first.DeltaR(event.wPairQuarks.second);
         }},
        {"wwDelPhi",
         [](const AnalysisEvent& event) -> float {
             if (event.jetIndex.size() < 3)
             {
                 return -1.;
             }
             return event.wPairQuarks.first.DeltaPhi(
                 event.wPairQuarks.second);
         }},
        {"lbDelR",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempJet1;
             if (event.bTagIndex.size() < 1)
             {
                 return -10.;
             }
             float smearValue =
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]];
             tempJet1.SetPxPyPzE(
                 event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempJet1 *= smearValue;
             return tempJet1.DeltaR(event.wLepton);
         }},
        {"lbDelPhi",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempJet1;
             if (event.bTagIndex.size() < 1)
             {
                 return -10.;
             }
             float smearValue =
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]];
             tempJet1.SetPxPyPzE(
                 event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempJet1 *= smearValue;
             return tempJet1.DeltaPhi(event.wLepton);
         }},
        {"zLepDelR",
         [](const AnalysisEvent& event) -> float {
             return (
                 event.zPairLeptons.first.DeltaR(event.zPairLeptons.second));
         }},
        {"zLepDelPhi",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first.DeltaPhi(
                 event.zPairLeptons.second));
         }},
        {"zLep1Quark1DelR",
         [](const AnalysisEvent& event) -> float {
             return (
                 event.zPairLeptons.first.DeltaR(event.wPairQuarks.first));
         }},
        {"zLep1Quark1DelPhi",
         [](const AnalysisEvent& event) -> float {
             return (
                 event.zPairLeptons.first.DeltaPhi(event.wPairQuarks.first));
         }},
        {"zLep1Quark2DelR",
         [](const AnalysisEvent& event) -> float {
             return (
                 event.zPairLeptons.first.DeltaR(event.wPairQuarks.second));
         }},
        {"zLep1Quark2DelPhi",
         [](const AnalysisEvent& event) -> float {
             return (
                 event.zPairLeptons.first.DeltaPhi(event.wPairQuarks.second));
         }},
        {"zLep2Quark1DelR",
         [](const AnalysisEvent& event) -> float {
             return (
                 event.zPairLeptons.second.DeltaR(event.wPairQuarks.first));
         }},
        {"zLep2Quark1DelPhi",
         [](const AnalysisEvent& event) -> float {
             return (
                 event.zPairLeptons.second.DeltaPhi(event.wPairQuarks.first));
         }},
        {"zLep2Quark2DelR",
         [](const AnalysisEvent& event) -> float {
             return (
                 event.zPairLeptons.second.DeltaR(event.wPairQuarks.second));
         }},
        {"zLep2Quark2DelPhi",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.second.DeltaPhi(
                 event.wPairQuarks.second));
         }},
        {"zLep1BjetDelR",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() < 1)
             {
                 return -10.;
             }
             TLorentzVector tempJet1;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempJet1.SetPxPyPzE(
                 event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempJet1 *= smearValue;
             return (event.zPairLeptons.first.DeltaR(tempJet1));
         }},
        {"zLep1BjetDelPhi",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() < 1)
             {
                 return -10.;
             }
             TLorentzVector tempJet1;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempJet1.SetPxPyPzE(
                 event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempJet1 *= smearValue;
             return (event.zPairLeptons.first.DeltaPhi(tempJet1));
         }},
        {"zLep2BjetDelR",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() < 1)
             {
                 return -10.;
             }
             TLorentzVector tempJet1;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempJet1.SetPxPyPzE(
                 event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempJet1 *= smearValue;
             return (event.zPairLeptons.second.DeltaR(tempJet1));
         }},
        {"zLep2BjetDelPhi",
         [](const AnalysisEvent& event) -> float {
             if (event.bTagIndex.size() < 1)
             {
                 return -10.;
             }
             TLorentzVector tempJet1;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempJet1.SetPxPyPzE(
                 event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempJet1 *= smearValue;
             return (event.zPairLeptons.second.DeltaPhi(tempJet1));
         }},
        {"lepHt",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .Pt();
         }},
        {"wQuarkHt",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .Pt();
         }},
        {"jetHt",
         [](const AnalysisEvent& event) -> float {
             float jetHt{0.0};
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt = event.jetIndex.begin();
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue{event.jetSmearValue[*jetIt]};
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     jetHt += tempJet.Pt();
                 }
             }
             return jetHt;
         }},
        {"totHt",
         [](const AnalysisEvent& event) -> float {
             float totHt{0.0};
             totHt +=
                 (event.zPairLeptons.first + event.zPairLeptons.second).Pt();
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt{event.jetIndex.begin()};
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue{event.jetSmearValue[*jetIt]};
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totHt += tempJet.Pt();
                 }
             }
             return totHt;
         }},
        {"totHtOverPt",
         [](const AnalysisEvent& event) -> float {
             float totHt{0.0};
             totHt +=
                 (event.zPairLeptons.first + event.zPairLeptons.second).Pt();
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt{event.jetIndex.begin()};
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue{event.jetSmearValue[*jetIt]};
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totHt += tempJet.Pt();
                 }
             }
             float totPx{0.0};
             float totPy{0.0};
             totPx +=
                 (event.zPairLeptons.first + event.zPairLeptons.second).Px();
             totPy +=
                 (event.zPairLeptons.first + event.zPairLeptons.second).Py();
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt{event.jetIndex.begin()};
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue{event.jetSmearValue[*jetIt]};
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totPx += tempJet.Px();
                     totPy += tempJet.Py();
                 }
             }

             return totHt / std::sqrt(totPx * totPx + totPy * totPy);
         }},
        {"totPt",
         [](const AnalysisEvent& event) -> float {
             float totPx{0.0};
             float totPy{0.0};
             totPx +=
                 (event.zPairLeptons.first + event.zPairLeptons.second).Px();
             totPy +=
                 (event.zPairLeptons.first + event.zPairLeptons.second).Py();
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt{event.jetIndex.begin()};
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue{event.jetSmearValue[*jetIt]};
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totPx += tempJet.Px();
                     totPy += tempJet.Py();
                 }
             }
             return std::sqrt(totPx * totPx + totPy * totPy);
         }},
        {"totEta",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector totVec;
             totVec = event.zPairLeptons.first + event.zPairLeptons.second;
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt{event.jetIndex.begin()};
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue{event.jetSmearValue[*jetIt]};
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                 }
             }
             return std::abs(totVec.Eta());
         }},
        {"totM",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector totVec;
             totVec = event.zPairLeptons.first + event.zPairLeptons.second;
             if (event.jetIndex.size() > 0)
             {
                 for (auto jetIt{event.jetIndex.begin()};
                      jetIt != event.jetIndex.end();
                      ++jetIt)
                 {
                     TLorentzVector tempJet;
                     float smearValue{event.jetSmearValue[*jetIt]};
                     tempJet.SetPxPyPzE(event.jetPF2PATPx[*jetIt],
                                        event.jetPF2PATPy[*jetIt],
                                        event.jetPF2PATPz[*jetIt],
                                        event.jetPF2PATE[*jetIt]);
                     tempJet *= smearValue;
                     totVec += tempJet;
                 }
             }
             return totVec.M();
         }},
        {"wzDelR",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .DeltaR(event.wPairQuarks.first + event.wPairQuarks.second);
         }},
        {"wzDelPhi",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .DeltaPhi(event.wPairQuarks.first
                           + event.wPairQuarks.second);
         }},
        {"zQuark1DelR",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .DeltaR(event.wPairQuarks.first);
         }},
        {"zQuark1DelPhi",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .DeltaPhi(event.wPairQuarks.first);
         }},
        {"zQuark2DelR",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .DeltaR(event.wPairQuarks.second);
         }},
        {"zQuark2DelPhi",
         [](const AnalysisEvent& event) -> float {
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .DeltaPhi(event.wPairQuarks.second);
         }},
        {"zTopDelR",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempBjet;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempBjet.SetPtEtaPhiE(
                 event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempBjet *= smearValue;
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .DeltaR(event.wPairQuarks.first + event.wPairQuarks.second
                         + tempBjet);
         }},
        {"zTopDelPhi",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempBjet;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempBjet.SetPtEtaPhiE(
                 event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempBjet *= smearValue;
             return (event.zPairLeptons.first + event.zPairLeptons.second)
                 .DeltaPhi(event.wPairQuarks.first + event.wPairQuarks.second
                           + tempBjet);
         }},
        {"zl1TopDelR",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempBjet;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempBjet.SetPtEtaPhiE(
                 event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempBjet *= smearValue;
             return (event.zPairLeptons.first)
                 .DeltaR(event.wPairQuarks.first + event.wPairQuarks.second
                         + tempBjet);
         }},
        {"zl1TopDelPhi",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempBjet;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempBjet.SetPtEtaPhiE(
                 event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempBjet *= smearValue;
             return (event.zPairLeptons.first)
                 .DeltaPhi(event.wPairQuarks.first + event.wPairQuarks.second
                           + tempBjet);
         }},
        {"zl2TopDelR",
         [](const AnalysisEvent& event) -> float {
             TLorentzVector tempBjet;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempBjet.SetPtEtaPhiE(
                 event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempBjet *= smearValue;
             return (event.zPairLeptons.second)
                 .DeltaR(event.wPairQuarks.first + event.wPairQuarks.second
                         + tempBjet);
         }},
        {"zl2TopDelPhi", [](const AnalysisEvent& event) -> float {
             TLorentzVector tempBjet;
             float smearValue{
                 event.jetSmearValue[event.jetIndex[event.bTagIndex[0]]]};
             tempBjet.SetPtEtaPhiE(
                 event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATEta[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATPhi[event.jetIndex[event.bTagIndex[0]]],
                 event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
             tempBjet *= smearValue;
             return (event.zPairLeptons.second)
                 .DeltaPhi(event.wPairQuarks.first + event.wPairQuarks.second
                           + tempBjet);
         }}};
}

void Plots::fillAllPlots(const AnalysisEvent& event, const float eventWeight)
{
    for (unsigned i{0}; i < plotPoint.size(); i++)
    {
        if (plotPoint[i].fillPlot)
        {
            plotPoint[i].plotHist->Fill((this->plotPoint[i].fillExp)(event),
                                        eventWeight);
        }
    }
}

void Plots::saveAllPlots()
{
    for (unsigned i{0}; i < plotPoint.size(); i++)
    {
        plotPoint[i].plotHist->SaveAs(
            ("plots/" + plotPoint[i].name + ".root").c_str());
        //    plotPoint[i].plotHist->SaveAs(("plots/"+plotPoint[i].name +
        //    ".pdf").c_str());
    }
}
