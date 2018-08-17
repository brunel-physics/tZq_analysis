#include "plots.hpp"

#include "TH1F.h"
#include "TLorentzVector.h"

#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

Plots::Plots(std::vector<std::string> titles,
             std::vector<std::string> names,
             std::vector<float> xMins,
             std::vector<float> xMaxs,
             std::vector<int> nBins,
             std::vector<std::string> fillExps,
             std::vector<std::string> xAxisLabels,
             std::vector<int> cutStage,
             unsigned thisCutStage,
             std::string postfixName)
{
    // Get the function pointer map for later custopmisation. This is gonna be
    // great, I promise.
    std::map<std::string, float (Plots::*)(AnalysisEvent*)> functionPointerMap{
        getFncPtrMap()};

    plotPoint = std::vector<plot>(names.size());
    for (unsigned i{0}; i < names.size(); i++)
    {
        std::string plotName = names[i] + "_" + postfixName;
        plotPoint[i].name = plotName;
        plotPoint[i].title = titles[i];
        plotPoint[i].fillExp = functionPointerMap[fillExps[i]];
        plotPoint[i].xAxisLabel = xAxisLabels[i];
        plotPoint[i].plotHist =
            new TH1F{plotName.c_str(),
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

// clang-format off
std::map<std::string, float (Plots::*)(AnalysisEvent*)> Plots::getFncPtrMap()
{
    std::map<std::string, float (Plots::*)(AnalysisEvent*)> functionPointerMap;
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep1Pt", &Plots::fillLepton1Pt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep1Eta", &Plots::fillLepton1Eta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep2Pt", &Plots::fillLepton2Pt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep2Eta", &Plots::fillLepton2Eta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep1RelIso", &Plots::fillLepton1RelIso));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep2RelIso", &Plots::fillLepton2RelIso));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep1Phi", &Plots::fillLepton1Phi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep2Phi", &Plots::fillLepton2Phi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wQuark1Pt", &Plots::fillWbosonQuark1Pt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wQuark1Eta", &Plots::fillWbosonQuark1Eta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wQuark1Phi", &Plots::fillWbosonQuark1Phi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wQuark2Pt", &Plots::fillWbosonQuark2Pt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wQuark2Eta", &Plots::fillWbosonQuark2Eta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wQuark2Phi", &Plots::fillWbosonQuark2Phi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("met", &Plots::fillMetPlot));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("numbJets", &Plots::getNumberOfJets));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totalJetMass", &Plots::fillTotalJetMass));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totalJetPt", &Plots::fillTotalJetPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totalJetEta", &Plots::fillTotalJetEta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totalJetPhi", &Plots::fillTotalJetPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("leadingJetPt", &Plots::fillLeadingJetPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("leadingJetEta", &Plots::fillLeadingJetEta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("leadingJetPhi", &Plots::fillLeadingJetPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("leadingJetBDisc", &Plots::fillLeadingJetBDisc));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("secondJetPt", &Plots::fillSecondJetPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("secondJetEta", &Plots::fillSecondJetEta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("secondJetPhi", &Plots::fillSecondJetPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("secondJetBDisc", &Plots::fillSecondJetBDisc));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("thirdJetPt", &Plots::fillThirdJetPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("thirdJetEta", &Plots::fillThirdJetEta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("thirdJetPhi", &Plots::fillThirdJetPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("thirdJetBDisc", &Plots::fillThirdJetBDisc));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("fourthJetPt", &Plots::fillFourthJetPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("fourthJetEta", &Plots::fillFourthJetEta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("fourthJetPhi", &Plots::fillFourthJetPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("fourthJetBDisc", &Plots::fillFourthJetBDisc));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("numbBJets", &Plots::numbBJets));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("bTagDisc", &Plots::fillBtagDisc));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepton1Pt", &Plots::fillZLep1Pt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepton1Eta", &Plots::fillZLep1Eta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepton2Pt", &Plots::fillZLep2Pt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepton2Eta", &Plots::fillZLep2Eta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepton1RelIso", &Plots::fillZLep1RelIso));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepton2RelIso", &Plots::fillZLep2RelIso));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepton1Phi", &Plots::fillZLep1Phi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepton2Phi", &Plots::fillZLep2Phi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zPairMass", &Plots::fillZPairMass));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zPairPt", &Plots::fillZPairPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zPairEta", &Plots::fillZPairEta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zPairPhi", &Plots::fillZPairPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("topMass", &Plots::fillTopMass));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("topPt", &Plots::fillTopPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("topEta", &Plots::fillTopEta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("topPhi", &Plots::fillTopPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep1D0", &Plots::fillLepton1D0));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep2D0", &Plots::fillLepton2D0));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep1DBD0", &Plots::fillLepton1DBD0));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep2DBD0", &Plots::fillLepton2DBD0));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep1BeamSpotCorrectedD0", &Plots::fillLepton1BeamSpotCorrectedD0));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep2BeamSpotCorrectedD0", &Plots::fillLepton2BeamSpotCorrectedD0));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep1InnerTrackD0", &Plots::fillLepton1InnerTrackD0));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lep2InnerTrackD0", &Plots::fillLepton2InnerTrackD0));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wTransverseMass", &Plots::fillwTransverseMass));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("jjDelR", &Plots::filljjDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("jjDelPhi", &Plots::filljjDelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wwDelR", &Plots::fillwwDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wwDelPhi", &Plots::fillwwDelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lbDelR", &Plots::filllbDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lbDelPhi", &Plots::filllbDelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepDelR", &Plots::fillZLepDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLepDelPhi", &Plots::fillZLepDelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep1Quark1DelR", &Plots::fillZLep1Quark1DelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep1Quark1DelPhi", &Plots::fillZLep1Quark1DelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep1Quark2DelR", &Plots::fillZLep1Quark2DelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep1Quark2DelPhi", &Plots::fillZLep1Quark2DelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep2Quark1DelR", &Plots::fillZLep2Quark1DelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep2Quark1DelPhi", &Plots::fillZLep2Quark1DelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep2Quark2DelR", &Plots::fillZLep2Quark2DelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep2Quark2DelPhi", &Plots::fillZLep2Quark2DelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep1BjetDelR", &Plots::fillZLep1BjetDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep1BjetDelPhi", &Plots::fillZLep1BjetDelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep2BjetDelR", &Plots::fillZLep2BjetDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zLep2BjetDelPhi", &Plots::fillZLep2BjetDelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("lepHt", &Plots::fillLepHt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wQuarkHt", &Plots::fillWquarkHt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("jetHt", &Plots::fillJetHt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totHt", &Plots::fillTotHt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totHtOverPt", &Plots::fillTotHtOverPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totPt", &Plots::fillTotPt));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totEta", &Plots::fillTotEta));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("totM", &Plots::fillTotM));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wzDelR", &Plots::fillWZdelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("wzDelPhi", &Plots::fillWZdelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zQuark1DelR", &Plots::fillZquark1DelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zQuark1DelPhi", &Plots::fillZquark1DelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zQuark2DelR", &Plots::fillZquark2DelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zQuark2DelPhi", &Plots::fillZquark2DelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zTopDelR", &Plots::fillZtopDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zTopDelPhi", &Plots::fillZtopDelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zl1TopDelR", &Plots::fillZLep1TopDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zl1TopDelPhi", &Plots::fillZLep1TopDelPhi));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zl2TopDelR", &Plots::fillZLep2TopDelR));
    functionPointerMap.insert(std::pair<std::string, float (Plots::*)(AnalysisEvent*)>("zl2TopDelPhi", &Plots::fillZLep2TopDelPhi));

    return functionPointerMap;
}
// clang-format on

float Plots::fillLepton1Pt(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        TLorentzVector tempVec{event->elePF2PATPX[event->electronIndexTight[0]],
                               event->elePF2PATPY[event->electronIndexTight[0]],
                               event->elePF2PATPZ[event->electronIndexTight[0]],
                               event->elePF2PATE[event->electronIndexTight[0]]};
        return tempVec.Pt();
    }
    else
    {
        TLorentzVector tempVec{event->muonPF2PATPX[event->muonIndexTight[0]],
                               event->muonPF2PATPY[event->muonIndexTight[0]],
                               event->muonPF2PATPZ[event->muonIndexTight[0]],
                               event->muonPF2PATE[event->muonIndexTight[0]]};
        tempVec *= event->muonMomentumSF[0];
        return tempVec.Pt();
    }
    return -10;
}

float Plots::fillLepton1Eta(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return std::abs(event->elePF2PATSCEta[event->electronIndexTight[0]]);
    }
    else
    {
        TLorentzVector tempVec{event->muonPF2PATPX[event->muonIndexTight[0]],
                               event->muonPF2PATPY[event->muonIndexTight[0]],
                               event->muonPF2PATPZ[event->muonIndexTight[0]],
                               event->muonPF2PATE[event->muonIndexTight[0]]};
        tempVec *= event->muonMomentumSF[0];
        return tempVec.Eta();
    }
    return -10;
}
float Plots::fillLepton2Pt(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        TLorentzVector tempVec{event->elePF2PATPX[event->electronIndexTight[1]],
                               event->elePF2PATPY[event->electronIndexTight[1]],
                               event->elePF2PATPZ[event->electronIndexTight[1]],
                               event->elePF2PATE[event->electronIndexTight[1]]};
        return tempVec.Pt();
    }
    else
    {
        TLorentzVector tempVec{event->muonPF2PATPX[event->muonIndexTight[1]],
                               event->muonPF2PATPY[event->muonIndexTight[1]],
                               event->muonPF2PATPZ[event->muonIndexTight[1]],
                               event->muonPF2PATE[event->muonIndexTight[1]]};
        tempVec *= event->muonMomentumSF[1];
        return tempVec.Pt();
    }
    return -10;
}

float Plots::fillLepton2Eta(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return std::abs(event->elePF2PATSCEta[event->electronIndexTight[1]]);
    }
    else
    {
        TLorentzVector tempVec{event->muonPF2PATPX[event->muonIndexTight[1]],
                               event->muonPF2PATPY[event->muonIndexTight[1]],
                               event->muonPF2PATPZ[event->muonIndexTight[1]],
                               event->muonPF2PATE[event->muonIndexTight[1]]};
        tempVec *= event->muonMomentumSF[1];
        return tempVec.Pt();
    }
    return -10;
}
float Plots::fillLepton3Pt(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 2)
    {
        return event->elePF2PATPT[event->electronIndexTight[2]];
    }
    if (event->muonIndexTight.size() > 2)
    {
        return event->muonPF2PATPt[event->muonIndexTight[2]];
    }
    if (event->electronIndexTight.size() > 1)
    {
        return event->muonPF2PATPt[event->muonIndexTight[0]];
    }
    else
    {
        return event->elePF2PATPT[event->electronIndexTight[0]];
    }

    return -10;
}

float Plots::fillLepton3Eta(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 2)
    {
        return std::abs(event->elePF2PATEta[event->electronIndexTight[2]]);
    }
    if (event->muonIndexTight.size() > 2)
    {
        return std::abs(event->muonPF2PATEta[event->muonIndexTight[2]]);
    }
    if (event->electronIndexTight.size() > 1)
    {
        return std::abs(event->muonPF2PATEta[event->muonIndexTight[0]]);
    }
    else
    {
        return std::abs(event->elePF2PATEta[event->electronIndexTight[0]]);
    }
    return -10;
}

float Plots::fillWbosonQuark1Pt(AnalysisEvent* event)
{
    return event->wPairQuarks.first.Pt();
}

float Plots::fillWbosonQuark1Eta(AnalysisEvent* event)
{
    return std::abs(event->wPairQuarks.first.Eta());
}

float Plots::fillWbosonQuark1Phi(AnalysisEvent* event)
{
    return event->wPairQuarks.first.Phi();
}

float Plots::fillWbosonQuark2Pt(AnalysisEvent* event)
{
    return event->wPairQuarks.second.Pt();
}

float Plots::fillWbosonQuark2Eta(AnalysisEvent* event)
{
    return std::abs(event->wPairQuarks.second.Eta());
}

float Plots::fillWbosonQuark2Phi(AnalysisEvent* event)
{
    return event->wPairQuarks.second.Phi();
}

float Plots::fillLepton1RelIso(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return event->elePF2PATComRelIsoRho[event->electronIndexTight[0]];
    }
    else
    {
        return event->muonPF2PATComRelIsodBeta[event->muonIndexTight[0]];
    }
    return -10;
}

float Plots::fillLepton2RelIso(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return event->elePF2PATComRelIsoRho[event->electronIndexTight[1]];
    }
    else
    {
        return event->muonPF2PATComRelIsodBeta[event->muonIndexTight[1]];
    }
    return -10;
}

float Plots::fillLepton3RelIso(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 2)
    {
        return event->elePF2PATComRelIsoRho[event->electronIndexTight[2]];
    }
    if (event->muonIndexTight.size() > 2)
    {
        return event->muonPF2PATComRelIsodBeta[event->muonIndexTight[2]];
    }
    if (event->electronIndexTight.size() > 1)
    {
        return event->muonPF2PATComRelIsodBeta[event->muonIndexTight[0]];
    }
    else
    {
        return event->elePF2PATComRelIsoRho[event->electronIndexTight[0]];
    }
    return -10;
}

float Plots::fillLepton1Phi(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return (event->elePF2PATPhi[event->electronIndexTight[0]]);
    }
    else
    {
        TLorentzVector tempVec{event->muonPF2PATPX[event->muonIndexTight[0]],
                               event->muonPF2PATPY[event->muonIndexTight[0]],
                               event->muonPF2PATPZ[event->muonIndexTight[0]],
                               event->muonPF2PATE[event->muonIndexTight[0]]};
        tempVec *= event->muonMomentumSF[0];
        return tempVec.Phi();
    }
    return -10;
}
float Plots::fillLepton2Phi(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return (event->elePF2PATPhi[event->electronIndexTight[1]]);
    }
    else
    {
        TLorentzVector tempVec{event->muonPF2PATPX[event->muonIndexTight[1]],
                               event->muonPF2PATPY[event->muonIndexTight[1]],
                               event->muonPF2PATPZ[event->muonIndexTight[1]],
                               event->muonPF2PATE[event->muonIndexTight[1]]};
        tempVec *= event->muonMomentumSF[1];
        return tempVec.Phi();
    }
    return -10;
}
float Plots::fillLepton3Phi(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 2)
    {
        return (event->elePF2PATPhi[event->electronIndexTight[2]]);
    }
    if (event->muonIndexTight.size() > 2)
    {
        return event->muonPF2PATPhi[event->muonIndexTight[2]];
    }
    if (event->electronIndexTight.size() > 1)
    {
        return event->muonPF2PATPhi[event->muonIndexTight[0]];
    }
    else
    {
        return event->elePF2PATPhi[event->electronIndexTight[0]];
    }
    return -10;
}

float Plots::fillMetPlot(AnalysisEvent* event)
{
    return event->metPF2PATEt;
}

float Plots::getNumberOfJets(AnalysisEvent* event)
{
    return event->jetIndex.size();
}

float Plots::fillTotalJetMass(AnalysisEvent* event)
{
    TLorentzVector totalJet;
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totalJet += tempJet;
        }
        return totalJet.M();
    }
    else
    {
        return -10;
    }
}

float Plots::fillTotalJetPt(AnalysisEvent* event)
{
    TLorentzVector totalJet;
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totalJet += tempJet;
        }
        return totalJet.Pt();
    }
    else
    {
        return -10;
    }
}

float Plots::fillTotalJetEta(AnalysisEvent* event)
{
    TLorentzVector totalJet;
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totalJet += tempJet;
        }
        return totalJet.Eta();
    }
    else
    {
        return -10;
    }
}

float Plots::fillTotalJetPhi(AnalysisEvent* event)
{
    TLorentzVector totalJet;
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totalJet += tempJet;
        }
        return totalJet.Phi();
    }
    else
    {
        return -10;
    }
}

float Plots::fillLeadingJetPt(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 0)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[0]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[0]],
                           event->jetPF2PATPy[event->jetIndex[0]],
                           event->jetPF2PATPz[event->jetIndex[0]],
                           event->jetPF2PATE[event->jetIndex[0]]);
        tempJet *= smearValue;
        return tempJet.Pt();
    }
    else
    {
        return -10;
    }
}

float Plots::fillSecondJetPt(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 1)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[1]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[1]],
                           event->jetPF2PATPy[event->jetIndex[1]],
                           event->jetPF2PATPz[event->jetIndex[1]],
                           event->jetPF2PATE[event->jetIndex[1]]);
        tempJet *= smearValue;
        return tempJet.Pt();
    }
    else
    {
        return -10;
    }
}

float Plots::fillThirdJetPt(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 2)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[2]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[2]],
                           event->jetPF2PATPy[event->jetIndex[2]],
                           event->jetPF2PATPz[event->jetIndex[2]],
                           event->jetPF2PATE[event->jetIndex[2]]);
        tempJet *= smearValue;
        return tempJet.Pt();
    }
    else
    {
        return -10;
    }
}

float Plots::fillFourthJetPt(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 3)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[3]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[3]],
                           event->jetPF2PATPy[event->jetIndex[3]],
                           event->jetPF2PATPz[event->jetIndex[3]],
                           event->jetPF2PATE[event->jetIndex[3]]);
        tempJet *= smearValue;
        return tempJet.Pt();
    }
    else
    {
        return -10;
    }
}

float Plots::fillLeadingJetEta(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 0)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[0]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[0]],
                           event->jetPF2PATPy[event->jetIndex[0]],
                           event->jetPF2PATPz[event->jetIndex[0]],
                           event->jetPF2PATE[event->jetIndex[0]]);
        tempJet *= smearValue;
        return tempJet.Eta();
    }
    else
    {
        return -10;
    }
}

float Plots::fillSecondJetEta(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 1)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[1]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[1]],
                           event->jetPF2PATPy[event->jetIndex[1]],
                           event->jetPF2PATPz[event->jetIndex[1]],
                           event->jetPF2PATE[event->jetIndex[1]]);
        tempJet *= smearValue;
        return tempJet.Eta();
    }
    else
    {
        return -10;
    }
}

float Plots::fillThirdJetEta(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 2)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[2]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[2]],
                           event->jetPF2PATPy[event->jetIndex[2]],
                           event->jetPF2PATPz[event->jetIndex[2]],
                           event->jetPF2PATE[event->jetIndex[2]]);
        tempJet *= smearValue;
        return tempJet.Eta();
    }
    else
    {
        return -10;
    }
}

float Plots::fillFourthJetEta(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 3)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[3]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[3]],
                           event->jetPF2PATPy[event->jetIndex[3]],
                           event->jetPF2PATPz[event->jetIndex[3]],
                           event->jetPF2PATE[event->jetIndex[3]]);
        tempJet *= smearValue;
        return tempJet.Eta();
    }
    else
    {
        return -10;
    }
}

float Plots::fillLeadingJetPhi(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 0)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[0]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[0]],
                           event->jetPF2PATPy[event->jetIndex[0]],
                           event->jetPF2PATPz[event->jetIndex[0]],
                           event->jetPF2PATE[event->jetIndex[0]]);
        tempJet *= smearValue;
        return tempJet.Phi();
    }
    else
    {
        return -10;
    }
}

float Plots::fillSecondJetPhi(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 1)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[1]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[1]],
                           event->jetPF2PATPy[event->jetIndex[1]],
                           event->jetPF2PATPz[event->jetIndex[1]],
                           event->jetPF2PATE[event->jetIndex[1]]);
        tempJet *= smearValue;
        return tempJet.Phi();
    }
    else
    {
        return -10;
    }
}

float Plots::fillThirdJetPhi(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 2)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[2]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[2]],
                           event->jetPF2PATPy[event->jetIndex[2]],
                           event->jetPF2PATPz[event->jetIndex[2]],
                           event->jetPF2PATE[event->jetIndex[2]]);
        tempJet *= smearValue;
        return tempJet.Phi();
    }
    else
    {
        return -10;
    }
}

float Plots::fillFourthJetPhi(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 3)
    {
        TLorentzVector tempJet;
        float smearValue = event->jetSmearValue[event->jetIndex[3]];
        tempJet.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[3]],
                           event->jetPF2PATPy[event->jetIndex[3]],
                           event->jetPF2PATPz[event->jetIndex[3]],
                           event->jetPF2PATE[event->jetIndex[3]]);
        tempJet *= smearValue;
        return tempJet.Phi();
    }
    else
    {
        return -10;
    }
}

float Plots::fillLeadingJetBDisc(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 0)
    {
        return event->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
            [event->jetIndex[0]];
    }
    else
    {
        return -10;
    }
}

float Plots::fillSecondJetBDisc(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 1)
    {
        return event->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
            [event->jetIndex[1]];
    }
    else
    {
        return -10;
    }
}

float Plots::fillThirdJetBDisc(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 2)
    {
        return event->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
            [event->jetIndex[2]];
    }
    else
    {
        return -10;
    }
}

float Plots::fillFourthJetBDisc(AnalysisEvent* event)
{
    if (event->jetIndex.size() > 3)
    {
        return event->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
            [event->jetIndex[3]];
    }
    else
    {
        return -10;
    }
}

float Plots::numbBJets(AnalysisEvent* event)
{
    return event->bTagIndex.size();
}

float Plots::fillBtagDisc(AnalysisEvent* event)
{
    if (event->bTagIndex.size() > 0)
    {
        return event->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
            [event->jetIndex[event->bTagIndex[0]]];
    }
    return -10;
}

float Plots::fillZLep1Pt(AnalysisEvent* event)
{
    return event->zPairLeptons.first.Pt();
}

float Plots::fillZLep2Pt(AnalysisEvent* event)
{
    return event->zPairLeptons.second.Pt();
}

float Plots::fillZLep1Eta(AnalysisEvent* event)
{
    return std::abs(event->zPairLeptons.first.Eta());
}

float Plots::fillZLep2Eta(AnalysisEvent* event)
{
    return std::abs(event->zPairLeptons.second.Eta());
}

float Plots::fillZLep1Phi(AnalysisEvent* event)
{
    return event->zPairLeptons.first.Phi();
}

float Plots::fillZLep2Phi(AnalysisEvent* event)
{
    return event->zPairLeptons.second.Phi();
}

float Plots::fillZLep1RelIso(AnalysisEvent* event)
{
    return event->zPairRelIso.first;
}

float Plots::fillZLep2RelIso(AnalysisEvent* event)
{
    return event->zPairRelIso.second;
}

float Plots::fillWLepPt(AnalysisEvent* event)
{
    return event->wLepton.Pt();
}

float Plots::fillWLepEta(AnalysisEvent* event)
{
    return std::abs(event->wLepton.Eta());
}

float Plots::fillWLepPhi(AnalysisEvent* event)
{
    return event->wLepton.Phi();
}

float Plots::fillWLepRelIso(AnalysisEvent* event)
{
    return event->wLeptonRelIso;
}

float Plots::fillZPairMass(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second).M();
}

float Plots::fillZPairPt(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second).Pt();
}

float Plots::fillZPairEta(AnalysisEvent* event)
{
    return std::abs(
        (event->zPairLeptons.first + event->zPairLeptons.second).Eta());
}

float Plots::fillZPairPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second).Phi();
}

float Plots::fillWPair1Mass(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->wLepton).M();
}

float Plots::fillWPair2Mass(AnalysisEvent* event)
{
    return (event->zPairLeptons.second + event->wLepton).M();
}

float Plots::fillWPairMass(AnalysisEvent* event)
{
    return (event->wPairQuarks.first + event->wPairQuarks.second).M();
}

float Plots::fillWPairPt(AnalysisEvent* event)
{
    return (event->wPairQuarks.first + event->wPairQuarks.second).Pt();
}

float Plots::fillWPairEta(AnalysisEvent* event)
{
    return (event->wPairQuarks.first + event->wPairQuarks.second).Eta();
}

float Plots::fillWPairPhi(AnalysisEvent* event)
{
    return (event->wPairQuarks.first + event->wPairQuarks.second).Phi();
}

float Plots::fillLeptonMass(AnalysisEvent* event)
{
    return (event->zPairLeptons.second + event->zPairLeptons.first
            + event->wLepton)
        .M();
}

float Plots::fillTopMass(AnalysisEvent* event)
{
    if (event->bTagIndex.size() > 0)
    {
        TLorentzVector tempBjet;
        float smearValue = event->jetSmearValue[event->bTagIndex[0]];
        tempBjet.SetPtEtaPhiE(
            event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
        tempBjet *= smearValue;
        return (tempBjet + event->wPairQuarks.first + event->wPairQuarks.second)
            .M();
    }
    else
    {
        return -10;
    }
}

float Plots::fillTopPt(AnalysisEvent* event)
{
    if (event->bTagIndex.size() > 0)
    {
        TLorentzVector tempBjet;
        float smearValue = event->jetSmearValue[event->bTagIndex[0]];
        tempBjet.SetPtEtaPhiE(
            event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
        tempBjet *= smearValue;
        return (tempBjet + event->wPairQuarks.first + event->wPairQuarks.second)
            .Pt();
    }
    else
    {
        return -10;
    }
}

float Plots::fillTopEta(AnalysisEvent* event)
{
    if (event->bTagIndex.size() > 0)
    {
        TLorentzVector tempBjet;
        float smearValue = event->jetSmearValue[event->bTagIndex[0]];
        tempBjet.SetPtEtaPhiE(
            event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
        tempBjet *= smearValue;
        return std::abs(
            (tempBjet + event->wPairQuarks.first + event->wPairQuarks.second)
                .Eta());
    }
    else
    {
        return -10;
    }
}

float Plots::fillTopPhi(AnalysisEvent* event)
{
    if (event->bTagIndex.size() > 0)
    {
        TLorentzVector tempBjet;
        float smearValue = event->jetSmearValue[event->bTagIndex[0]];
        tempBjet.SetPtEtaPhiE(
            event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
            event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
        tempBjet *= smearValue;
        return (tempBjet + event->wPairQuarks.first + event->wPairQuarks.second)
            .Phi();
    }
    else
    {
        return -10;
    }
}
float Plots::fillLepton1D0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return (event->elePF2PATD0PV[event->electronIndexTight[0]]);
    }
    else
    {
        return (event->muonPF2PATDBPV[event->muonIndexTight[0]]);
    }
    return -10;
}
float Plots::fillLepton2D0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return (event->elePF2PATD0PV[event->electronIndexTight[1]]);
    }
    else
    {
        return (event->muonPF2PATDBPV[event->muonIndexTight[1]]);
    }
    return -10;
}
float Plots::fillLepton3D0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 2)
    {
        return (event->elePF2PATD0PV[event->electronIndexTight[2]]);
    }
    if (event->muonIndexTight.size() > 2)
    {
        return event->muonPF2PATDBPV[event->muonIndexTight[2]];
    }
    if (event->electronIndexTight.size() > 1)
    {
        return event->muonPF2PATDBPV[event->muonIndexTight[0]];
    }
    else
    {
        return event->elePF2PATD0PV[event->electronIndexTight[0]];
    }
    return -10;
}

float Plots::fillLepton1DBD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return (event->elePF2PATTrackDBD0[event->electronIndexTight[0]]);
    }
    else
    {
        return (event->muonPF2PATTrackDBD0[event->muonIndexTight[0]]);
    }
    return -10;
}
float Plots::fillLepton2DBD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return (event->elePF2PATTrackDBD0[event->electronIndexTight[1]]);
    }
    else
    {
        return (event->muonPF2PATTrackDBD0[event->muonIndexTight[1]]);
    }
    return -10;
}
float Plots::fillLepton3DBD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 2)
    {
        return (event->elePF2PATTrackDBD0[event->electronIndexTight[2]]);
    }
    if (event->muonIndexTight.size() > 2)
    {
        return event->muonPF2PATTrackDBD0[event->muonIndexTight[2]];
    }
    if (event->electronIndexTight.size() > 1)
    {
        return event->muonPF2PATTrackDBD0[event->muonIndexTight[0]];
    }
    else
    {
        return event->elePF2PATTrackDBD0[event->electronIndexTight[0]];
    }
    return -10;
}

float Plots::fillLepton1BeamSpotCorrectedD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return (event->elePF2PATBeamSpotCorrectedTrackD0
                    [event->electronIndexTight[0]]);
    }
    else
    {
        return (event->muonPF2PATBeamSpotCorrectedD0[event->muonIndexTight[0]]);
    }
    return -10;
}
float Plots::fillLepton2BeamSpotCorrectedD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return (event->elePF2PATBeamSpotCorrectedTrackD0
                    [event->electronIndexTight[1]]);
    }
    else
    {
        return (event->muonPF2PATBeamSpotCorrectedD0[event->muonIndexTight[1]]);
    }
    return -10;
}
float Plots::fillLepton3BeamSpotCorrectedD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 2)
    {
        return (event->elePF2PATBeamSpotCorrectedTrackD0
                    [event->electronIndexTight[2]]);
    }
    if (event->muonIndexTight.size() > 2)
    {
        return event->muonPF2PATBeamSpotCorrectedD0[event->muonIndexTight[2]];
    }
    if (event->electronIndexTight.size() > 1)
    {
        return event->muonPF2PATBeamSpotCorrectedD0[event->muonIndexTight[0]];
    }
    else
    {
        return event
            ->elePF2PATBeamSpotCorrectedTrackD0[event->electronIndexTight[0]];
    }
    return -10;
}

float Plots::fillLepton1InnerTrackD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return -10;
    }
    else
    {
        return (event->muonPF2PATDBInnerTrackD0[event->muonIndexTight[0]]);
    }
    return -10;
}
float Plots::fillLepton2InnerTrackD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 1)
    {
        return -10;
    }
    else
    {
        return (event->muonPF2PATDBInnerTrackD0[event->muonIndexTight[1]]);
    }
    return -10;
}
float Plots::fillLepton3InnerTrackD0(AnalysisEvent* event)
{
    if (event->electronIndexTight.size() > 2)
    {
        return -10;
    }
    if (event->muonIndexTight.size() > 2)
    {
        return event->muonPF2PATDBInnerTrackD0[event->muonIndexTight[2]];
    }
    if (event->electronIndexTight.size() > 1)
    {
        return event->muonPF2PATDBInnerTrackD0[event->muonIndexTight[0]];
    }
    else
    {
        return -10;
    }
    return -10;
}

float Plots::fillwTransverseMass(AnalysisEvent* event)
{
    return std::sqrt(2 * event->wPairQuarks.first.Pt()
                     * event->wPairQuarks.second.Pt()
                     * (1
                        - std::cos(event->wPairQuarks.first.Phi()
                                   - event->wPairQuarks.second.Phi())));
}

float Plots::filljjDelR(AnalysisEvent* event)
{
    TLorentzVector tempJet1;
    TLorentzVector tempJet2;
    if (event->jetIndex.size() < 2)
    {
        return -1.;
    }
    float smearValue1 = event->jetSmearValue[event->jetIndex[0]];
    float smearValue2 = event->jetSmearValue[event->jetIndex[1]];
    tempJet1.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[0]],
                        event->jetPF2PATPy[event->jetIndex[0]],
                        event->jetPF2PATPz[event->jetIndex[0]],
                        event->jetPF2PATE[event->jetIndex[0]]);
    tempJet2.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[1]],
                        event->jetPF2PATPy[event->jetIndex[1]],
                        event->jetPF2PATPz[event->jetIndex[1]],
                        event->jetPF2PATE[event->jetIndex[1]]);
    tempJet1 *= smearValue1;
    tempJet2 *= smearValue2;
    return tempJet1.DeltaR(tempJet2);
}

float Plots::filljjDelPhi(AnalysisEvent* event)
{
    TLorentzVector tempJet1;
    TLorentzVector tempJet2;
    if (event->jetIndex.size() < 2)
    {
        return -1.;
    }
    float smearValue1 = event->jetSmearValue[event->jetIndex[0]];
    float smearValue2 = event->jetSmearValue[event->jetIndex[1]];
    tempJet1.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[0]],
                        event->jetPF2PATPy[event->jetIndex[0]],
                        event->jetPF2PATPz[event->jetIndex[0]],
                        event->jetPF2PATE[event->jetIndex[0]]);
    tempJet2.SetPxPyPzE(event->jetPF2PATPx[event->jetIndex[1]],
                        event->jetPF2PATPy[event->jetIndex[1]],
                        event->jetPF2PATPz[event->jetIndex[1]],
                        event->jetPF2PATE[event->jetIndex[1]]);
    tempJet1 *= smearValue1;
    tempJet2 *= smearValue2;
    return tempJet1.DeltaPhi(tempJet2);
}

float Plots::fillwwDelR(AnalysisEvent* event)
{
    if (event->jetIndex.size() < 3)
    {
        return -1.;
    }
    return event->wPairQuarks.first.DeltaR(event->wPairQuarks.second);
}

float Plots::fillwwDelPhi(AnalysisEvent* event)
{
    if (event->jetIndex.size() < 3)
    {
        return -1.;
    }
    return event->wPairQuarks.first.DeltaPhi(event->wPairQuarks.second);
}

float Plots::filllbDelR(AnalysisEvent* event)
{
    TLorentzVector tempJet1;
    if (event->bTagIndex.size() < 1)
    {
        return -10.;
    }
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempJet1.SetPxPyPzE(
        event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempJet1 *= smearValue;
    return tempJet1.DeltaR(event->wLepton);
}

float Plots::filllbDelPhi(AnalysisEvent* event)
{
    TLorentzVector tempJet1;
    if (event->bTagIndex.size() < 1)
    {
        return -10.;
    }
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempJet1.SetPxPyPzE(
        event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempJet1 *= smearValue;
    return tempJet1.DeltaPhi(event->wLepton);
}

float Plots::fillZLepDelR(AnalysisEvent* event)
{
    return (event->zPairLeptons.first.DeltaR(event->zPairLeptons.second));
}

float Plots::fillZLepDelPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.first.DeltaPhi(event->zPairLeptons.second));
}

float Plots::fillZLep1Quark1DelR(AnalysisEvent* event)
{
    return (event->zPairLeptons.first.DeltaR(event->wPairQuarks.first));
}

float Plots::fillZLep1Quark1DelPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.first.DeltaPhi(event->wPairQuarks.first));
}

float Plots::fillZLep1Quark2DelR(AnalysisEvent* event)
{
    return (event->zPairLeptons.first.DeltaR(event->wPairQuarks.second));
}

float Plots::fillZLep1Quark2DelPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.first.DeltaPhi(event->wPairQuarks.second));
}

float Plots::fillZLep2Quark1DelR(AnalysisEvent* event)
{
    return (event->zPairLeptons.second.DeltaR(event->wPairQuarks.first));
}

float Plots::fillZLep2Quark1DelPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.second.DeltaPhi(event->wPairQuarks.first));
}

float Plots::fillZLep2Quark2DelR(AnalysisEvent* event)
{
    return (event->zPairLeptons.second.DeltaR(event->wPairQuarks.second));
}

float Plots::fillZLep2Quark2DelPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.second.DeltaPhi(event->wPairQuarks.second));
}

float Plots::fillZLep1BjetDelR(AnalysisEvent* event)
{
    if (event->bTagIndex.size() < 1)
    {
        return -10.;
    }
    TLorentzVector tempJet1;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempJet1.SetPxPyPzE(
        event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempJet1 *= smearValue;
    return (event->zPairLeptons.first.DeltaR(tempJet1));
}

float Plots::fillZLep1BjetDelPhi(AnalysisEvent* event)
{
    if (event->bTagIndex.size() < 1)
    {
        return -10.;
    }
    TLorentzVector tempJet1;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempJet1.SetPxPyPzE(
        event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempJet1 *= smearValue;
    return (event->zPairLeptons.first.DeltaPhi(tempJet1));
}

float Plots::fillZLep2BjetDelR(AnalysisEvent* event)
{
    if (event->bTagIndex.size() < 1)
    {
        return -10.;
    }
    TLorentzVector tempJet1;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempJet1.SetPxPyPzE(
        event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempJet1 *= smearValue;
    return (event->zPairLeptons.second.DeltaR(tempJet1));
}

float Plots::fillZLep2BjetDelPhi(AnalysisEvent* event)
{
    if (event->bTagIndex.size() < 1)
    {
        return -10.;
    }
    TLorentzVector tempJet1;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempJet1.SetPxPyPzE(
        event->jetPF2PATPx[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPy[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPz[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempJet1 *= smearValue;
    return (event->zPairLeptons.second.DeltaPhi(tempJet1));
}

float Plots::fillLepHt(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second).Pt();
}

float Plots::fillWquarkHt(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second).Pt();
}

float Plots::fillJetHt(AnalysisEvent* event)
{
    float jetHt(0.0);
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            jetHt += tempJet.Pt();
        }
    }
    return jetHt;
}

float Plots::fillTotHt(AnalysisEvent* event)
{
    float totHt(0.0);
    totHt += (event->zPairLeptons.first + event->zPairLeptons.second).Pt();
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totHt += tempJet.Pt();
        }
    }
    return totHt;
}

float Plots::fillTotHtOverPt(AnalysisEvent* event)
{
    float totHt(0.0);
    totHt += (event->zPairLeptons.first + event->zPairLeptons.second).Pt();
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totHt += tempJet.Pt();
        }
    }
    float totPx(0.0), totPy(0.0);
    totPx += (event->zPairLeptons.first + event->zPairLeptons.second).Px();
    totPy += (event->zPairLeptons.first + event->zPairLeptons.second).Py();
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totPx += tempJet.Px();
            totPy += tempJet.Py();
        }
    }

    return totHt / std::sqrt(totPx * totPx + totPy * totPy);
}

float Plots::fillTotPt(AnalysisEvent* event)
{
    float totPx(0.0), totPy(0.0);
    totPx += (event->zPairLeptons.first + event->zPairLeptons.second).Px();
    totPy += (event->zPairLeptons.first + event->zPairLeptons.second).Py();
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totPx += tempJet.Px();
            totPy += tempJet.Py();
        }
    }
    return std::sqrt(totPx * totPx + totPy * totPy);
}

float Plots::fillTotEta(AnalysisEvent* event)
{
    TLorentzVector totVec;
    totVec = event->zPairLeptons.first + event->zPairLeptons.second;
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
        }
    }
    return std::abs(totVec.Eta());
}

float Plots::fillTotM(AnalysisEvent* event)
{
    TLorentzVector totVec;
    totVec = event->zPairLeptons.first + event->zPairLeptons.second;
    if (event->jetIndex.size() > 0)
    {
        for (auto jetIt = event->jetIndex.begin();
             jetIt != event->jetIndex.end();
             ++jetIt)
        {
            TLorentzVector tempJet;
            float smearValue = event->jetSmearValue[*jetIt];
            tempJet.SetPxPyPzE(event->jetPF2PATPx[*jetIt],
                               event->jetPF2PATPy[*jetIt],
                               event->jetPF2PATPz[*jetIt],
                               event->jetPF2PATE[*jetIt]);
            tempJet *= smearValue;
            totVec += tempJet;
        }
    }
    return totVec.M();
}

float Plots::fillWZdelR(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second)
        .DeltaR(event->wPairQuarks.first + event->wPairQuarks.second);
}

float Plots::fillWZdelPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second)
        .DeltaPhi(event->wPairQuarks.first + event->wPairQuarks.second);
}

float Plots::fillZquark1DelR(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second)
        .DeltaR(event->wPairQuarks.first);
}

float Plots::fillZquark1DelPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second)
        .DeltaPhi(event->wPairQuarks.first);
}

float Plots::fillZquark2DelR(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second)
        .DeltaR(event->wPairQuarks.second);
}

float Plots::fillZquark2DelPhi(AnalysisEvent* event)
{
    return (event->zPairLeptons.first + event->zPairLeptons.second)
        .DeltaPhi(event->wPairQuarks.second);
}

float Plots::fillZtopDelR(AnalysisEvent* event)
{
    TLorentzVector tempBjet;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempBjet.SetPtEtaPhiE(
        event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempBjet *= smearValue;
    return (event->zPairLeptons.first + event->zPairLeptons.second)
        .DeltaR(event->wPairQuarks.first + event->wPairQuarks.second
                + tempBjet);
}

float Plots::fillZtopDelPhi(AnalysisEvent* event)
{
    TLorentzVector tempBjet;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempBjet.SetPtEtaPhiE(
        event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempBjet *= smearValue;
    return (event->zPairLeptons.first + event->zPairLeptons.second)
        .DeltaPhi(event->wPairQuarks.first + event->wPairQuarks.second
                  + tempBjet);
}

float Plots::fillZLep1TopDelR(AnalysisEvent* event)
{
    TLorentzVector tempBjet;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempBjet.SetPtEtaPhiE(
        event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempBjet *= smearValue;
    return (event->zPairLeptons.first)
        .DeltaR(event->wPairQuarks.first + event->wPairQuarks.second
                + tempBjet);
}

float Plots::fillZLep1TopDelPhi(AnalysisEvent* event)
{
    TLorentzVector tempBjet;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempBjet.SetPtEtaPhiE(
        event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempBjet *= smearValue;
    return (event->zPairLeptons.first)
        .DeltaPhi(event->wPairQuarks.first + event->wPairQuarks.second
                  + tempBjet);
}

float Plots::fillZLep2TopDelR(AnalysisEvent* event)
{
    TLorentzVector tempBjet;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempBjet.SetPtEtaPhiE(
        event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempBjet *= smearValue;
    return (event->zPairLeptons.second)
        .DeltaR(event->wPairQuarks.first + event->wPairQuarks.second
                + tempBjet);
}

float Plots::fillZLep2TopDelPhi(AnalysisEvent* event)
{
    TLorentzVector tempBjet;
    float smearValue =
        event->jetSmearValue[event->jetIndex[event->bTagIndex[0]]];
    tempBjet.SetPtEtaPhiE(
        event->jetPF2PATPt[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATEta[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATPhi[event->jetIndex[event->bTagIndex[0]]],
        event->jetPF2PATE[event->jetIndex[event->bTagIndex[0]]]);
    tempBjet *= smearValue;
    return (event->zPairLeptons.second)
        .DeltaPhi(event->wPairQuarks.first + event->wPairQuarks.second
                  + tempBjet);
}

void Plots::fillAllPlots(AnalysisEvent* event, float eventWeight)
{
    for (unsigned i{0}; i < plotPoint.size(); i++)
    {
        if (plotPoint[i].fillPlot)
        {
            plotPoint[i].plotHist->Fill((this->*plotPoint[i].fillExp)(event),
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
