//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 23 14:04:00 2014 by ROOT version 5.32/00
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef _MvaEvent_hpp_
#define _MvaEvent_hpp_

#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <array>
#include <iostream>
#include <string>

#include "AnalysisEvent.hpp"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MvaEvent : public AnalysisEvent
{
    public:
    static constexpr size_t NJETS{15};
    static constexpr size_t NBJETS{10};
    static constexpr size_t NMUONS{3};

    bool isMC;
    Float_t eventWeight;
    Int_t zLep1Index;
    Int_t zLep2Index;
    Int_t wQuark1Index;
    Int_t wQuark2Index;
    Int_t jetInd[NJETS];
    Int_t bJetInd[NBJETS];
    Float_t muonMomentumSF[NMUONS];
    Float_t jetSmearValue[NJETS];

    // End MVA tree specific

    // List of branches
    TBranch* b_isMC;
    TBranch* b_eventWeight; //!
    TBranch* b_zLep1Index; //!
    TBranch* b_zLep2Index; //!
    TBranch* b_wQuark1Index; //!
    TBranch* b_wQuark2Index; //!
    TBranch* b_jetInd; //!
    TBranch* b_bJetInd; //!
    TBranch* b_muonMomentumSF; //!
    TBranch* b_jetSmearValue; //!

    MvaEvent(bool isMC = true,
             std::string triggerFlag = "",
             TTree* tree = nullptr,
             bool is2016 = false);
    virtual ~MvaEvent();
};

#endif

#ifdef MvaEvent_cxx
MvaEvent::MvaEvent(bool isMC,
                   std::string triggerFlag,
                   TTree* tree,
                   bool is2016)
    : AnalysisEvent{isMC, triggerFlag, tree, is2016}
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == nullptr)
    {
#ifdef SINGLE_TREE
        // The following code should be used if you want this class to access
        // a single tree instead of a chain
        TFile* f{(TFile*)gROOT->GetListOfFiles()->FindObject(
            "/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_9_0_MJP_skim.root")};
        if (!f || !f->IsOpen())
        {
            f = new TFile{"/data1/tW2012/mc/ttbarInclusive/"
                          "MC_Ntuple_out_9_0_MJP_skim.root"};
        }
        f->GetObject("tree", tree);

#else // SINGLE_TREE

        // The following code should be used if you want this class to access a
        // chain of trees.
        TChain* chain{new TChain{"tree", ""}};
        chain->Add("/data1/tW2012/mc/ttbarInclusive/"
                   "MC_Ntuple_out_100_0_Gu6_skim.root/tree");
        tree = chain;
#endif // SINGLE_TREE
    }

    fChain->SetBranchAddress("isMC", &isMC, &b_isMC);
    fChain->SetBranchAddress("eventWeight", &eventWeight, &b_eventWeight);
    fChain->SetBranchAddress("zLep1Index", &zLep1Index, &b_zLep1Index);
    fChain->SetBranchAddress("zLep2Index", &zLep2Index, &b_zLep2Index);
    fChain->SetBranchAddress("wQuark1Index", &wQuark1Index, &b_wQuark1Index);
    fChain->SetBranchAddress("wQuark2Index", &wQuark2Index, &b_wQuark2Index);
    fChain->SetBranchAddress("jetInd", jetInd, &b_jetInd);
    fChain->SetBranchAddress("bJetInd", bJetInd, &b_bJetInd);
    fChain->SetBranchAddress("muonMomentumSF", muonMomentumSF, &b_muonMomentumSF);
    fChain->SetBranchAddress("jetSmearValue", jetSmearValue, &b_jetSmearValue);
}

inline MvaEvent::~MvaEvent()
{
}

#endif // #ifdef MvaEvent_cxx
