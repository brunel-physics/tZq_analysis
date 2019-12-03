//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 23 14:04:00 2014 by ROOT version 5.32/00
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef _MvaEvent_hpp_
#define _MvaEvent_hpp_

#include "AnalysisEvent.hpp"

#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <array>
#include <iostream>
#include <string>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MvaEvent : public AnalysisEvent
{
    public:
    static constexpr size_t NJETS{15};
    static constexpr size_t NBJETS{10};
    static constexpr size_t NMUONS{2};

    bool isMC;
    Double_t eventWeight;
    Float_t muonMomentumSF[NMUONS];
    Int_t zLep1Index;
    Int_t zLep2Index;
    bool muonLeads;
    Int_t wQuark1Index;
    Int_t wQuark2Index;
    Int_t jetInd[NJETS];
    Int_t bJetInd[NBJETS];
    Float_t jetSmearValue[NJETS];

    // End MVA tree specific

    // List of branches
    TBranch* b_isMC;
    TBranch* b_eventWeight; //!
    TBranch* b_muonMomentumSF; //!
    TBranch* b_zLep1Index; //!
    TBranch* b_zLep2Index; //!
    TBranch* b_muonLeads; //!
    TBranch* b_wQuark1Index; //!
    TBranch* b_wQuark2Index; //!
    TBranch* b_jetInd; //!
    TBranch* b_bJetInd; //!
    TBranch* b_jetSmearValue; //!

    MvaEvent(bool isMC = true,
             TTree* tree = nullptr,
             bool is2016 = false);
    virtual ~MvaEvent();
};

inline MvaEvent::MvaEvent(bool isMC,
                          TTree* tree,
                          bool is2016)
    : AnalysisEvent{isMC, tree, is2016}
{
    fChain->SetBranchAddress("isMC", &isMC, &b_isMC);
    fChain->SetBranchAddress("eventWeight", &eventWeight, &b_eventWeight);
    fChain->SetBranchAddress("muonMomentumSF", &muonMomentumSF,
                             &b_muonMomentumSF);
    fChain->SetBranchAddress("zLep1Index", &zLep1Index, &b_zLep1Index);
    fChain->SetBranchAddress("zLep2Index", &zLep2Index, &b_zLep2Index);
    fChain->SetBranchAddress("muonLeads", &muonLeads, &b_muonLeads);
    fChain->SetBranchAddress("wQuark1Index", &wQuark1Index, &b_wQuark1Index);
    fChain->SetBranchAddress("wQuark2Index", &wQuark2Index, &b_wQuark2Index);
    fChain->SetBranchAddress("jetInd", jetInd, &b_jetInd);
    fChain->SetBranchAddress("bJetInd", bJetInd, &b_bJetInd);
    fChain->SetBranchAddress("jetSmearValue", jetSmearValue, &b_jetSmearValue);
}

inline MvaEvent::~MvaEvent()
{
}

#endif
