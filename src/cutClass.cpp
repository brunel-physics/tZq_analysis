#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "cutClass.hpp"

#include <boost/functional/hash.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <yaml-cpp/yaml.h>

Cuts::Cuts(const bool doPlots,
           const bool fillCutFlows,
           const bool invertLepCut,
           const bool is2016)
    : doPlots_{doPlots}
    , fillCutFlow_{fillCutFlows}
    , invertLepCut_{invertLepCut}
    , is2016_{is2016}

    , numTightEle_{3}
    , tightElePt_{15.}
    , tightElePtLeading_{35.}
    , tightEleEta_{2.5}
    , tightEleRelIso_{0.107587}

    , numLooseEle_{3}
    , looseElePt_{15.}
    , looseElePtLeading_{35.}
    , looseEleEta_{2.5}
    , looseEleRelIso_{0.15}

    , numTightMu_{0}
    , tightMuonPt_{12.}
    , tightMuonPtLeading_{27.}
    , tightMuonEta_{2.4}
    , tightMuonRelIso_{0.15}

    , numLooseMu_{0}
    , looseMuonPt_{12.}
    , looseMuonPtLeading_{27.}
    , looseMuonEta_{2.4}
    , looseMuonRelIso_{0.25}

    , invZMassCut_{20.}
    , invWMassCut_{20.}

    , numJets_{2}
    , maxJets_{4}
    , jetPt_{30.}
    , jetEta_{5.0}
    , jetIDDo_{true}

    , numbJets_{1}
    , maxbJets_{2}
    , maxbJetEta_{2.5}

    , bDiscCut_{is2016 ? 0.8484f : 0.8838f}

    , numcJets_{1}

    , rc_{is2016 ? "scaleFactors/2016/RoccoR2016.txt"
                 : "scaleFactors/2017/RoccoR2017.txt"}

    , lumiRunsBCDEF_{19713.888}
    , lumiRunsGH_{16146.178}

    , isMC_{true}

    , isNPL_{false}
    , isZplusCR_{false}

    , postLepSelTree_{nullptr}

    // Skips running trigger stuff
    , skipTrigger_{false}

    // Are we making b-tag efficiency plots?
    , makeBTagEffPlots_{false}
    , getBTagWeight_{false}

    // MET and mTW cuts go here.
    , metDileptonCut_{50.0}

{
    std::cout << "\nInitialises fine" << std::endl;
    initialiseJECCors();
    std::cout << "Gets past JEC Cors" << std::endl;

    if (!is2016_)
    {
        std::cout << "lumi set for Runs B-F: " << lumiRunsBCDEF_ << std::endl;
        std::cout << "lumi set for Runs G-H: " << lumiRunsGH_ << std::endl;

        std::cout << "\nLoad 2017 electron SFs from root file ... "
                  << std::endl;

        // Electron tight cut-based tight ID
        electronSFsFile = new TFile("scaleFactors/2017/"
                                    "egammaEffi.txt_EGM2D_runBCDEF_"
                                    "passingTight94X.root");

        // Electron reco SF
        h_eleSFs = dynamic_cast<TH2F*>(electronSFsFile->Get("EGamma_SF2D"));
        electronRecoFile = new TFile{
            "scaleFactors/2017/"
            "egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"}; // Electron Reco

        h_eleReco = dynamic_cast<TH2F*>(electronRecoFile->Get("EGamma_SF2D"));
        std::cout << "Got 2017 electron SFs!\n" << std::endl;

        std::cout << "Load 2017 muon SFs from root file ... " << std::endl;
        muonHltFile1 = new TFile{
            "scaleFactors/2017/HLT_Mu24_EfficienciesAndSF_RunBtoF.root"};
        muonIDsFile1 = new TFile{"scaleFactors/2017/Muon_RunBCDEF_SF_ID.root"};
        muonIsoFile1 = new TFile{"scaleFactors/2017/Muon_RunBCDEF_SF_ISO.root"};

        // Single muon HLT SF
        muonHltFile1->cd("IsoMu27_PtEtaBins");
        h_muonHlt1 = dynamic_cast<TH2F*>(
            muonHltFile1->Get("IsoMu27_PtEtaBins/abseta_pt_ratio"));

        // Hardcoded in 2017
        // // Tight ID
        // h_muonIDs1 = dynamic_cast<TH2D*>(
        //     muonIDsFile1->Get("NUM_TightID_DEN_genTracks_pt_abseta"));
        // h_muonPFiso1 = dynamic_cast<TH2D*>(
        //     muonIsoFile1->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"));
        // std::cout << "Got 2017 muon SFs!\n" << std::endl;
    }
    else
    {
        std::cout << "\nLoad 2016 electron SFs from root file ... "
                  << std::endl;

        // Single electron HLT SF
        electronHltFile =
            new TFile("scaleFactors/2016/"
                      "HLT_Ele32_eta2p1_WPTight_Gsf_FullRunRange.root");

        // Electron cut-based ID
        h_eleHlt = dynamic_cast<TH2F*>(electronHltFile->Get("SF"));
        electronSFsFile =
            new TFile("scaleFactors/2016/egammaEffi_Tight_80X.txt_EGM2D.root");
        h_eleSFs = dynamic_cast<TH2F*>(electronSFsFile->Get("EGamma_SF2D"));

        // Electron reco SF
        electronRecoFile =
            new TFile{"scaleFactors/2016/egammaRecoEffi.txt_EGM2D.root"};
        h_eleReco = dynamic_cast<TH2F*>(electronRecoFile->Get("EGamma_SF2D"));
        std::cout << "Got 2016 electron SFs!\n" << std::endl;

        std::cout << "Load 2016 muon SFs from root file ... " << std::endl;

        // Runs B-F (pre-HIP fix)
        muonHltFile1 = new TFile{"scaleFactors/2016/"
                                 "HLT_Mu24_EfficienciesAndSF_RunBtoF.root"};
        // Runs G-H (post-HIP fix)
        muonHltFile2 = new TFile{"scaleFactors/2016/"
                                 "HLT_Mu24_EfficienciesAndSF_RunGtoH.root"};

        // Runs B-F (pre-HIP fix)
        muonIDsFile1 =
            new TFile{"scaleFactors/2016/MuonID_EfficienciesAndSF_BCDEF.root"};
        // Runs G-H (post-HIP fix)
        muonIDsFile2 =
            new TFile{"scaleFactors/2016/MuonID_EfficienciesAndSF_GH.root"};

        // Runs B-F (pre-HIP fix)
        muonIsoFile1 =
            new TFile{"scaleFactors/2016/MuonISO_EfficienciesAndSF_BCDEF.root"};
        // Runs G-H (post-HIP fix)
        muonIsoFile2 =
            new TFile{"scaleFactors/2016/MuonISO_EfficienciesAndSF_GH.root"};

        // Single Muon HLT SF
        muonHltFile1->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        muonHltFile2->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        h_muonHlt1 = dynamic_cast<TH2F*>(muonHltFile1->Get(
            "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"));
        h_muonHlt2 = dynamic_cast<TH2F*>(muonHltFile2->Get(
            "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"));

        // Tight ID
        muonIDsFile1->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"); // Tight ID
        muonIDsFile2->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"); // Tight ID
        h_muonIDs1 = dynamic_cast<TH2F*>(
            muonIDsFile1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/"
                              "abseta_pt_ratio")); // Tight
                                                   // ID
        h_muonIDs2 = dynamic_cast<TH2F*>(
            muonIDsFile2->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/"
                              "abseta_pt_ratio")); // Tight
                                                   // ID

        // Tight Iso
        muonIsoFile1->cd("TightISO_TightID_pt_eta");
        muonIsoFile2->cd("TightISO_TightID_pt_eta"); // Tight Iso
        h_muonPFiso1 = dynamic_cast<TH2F*>(
            muonIsoFile1->Get("TightISO_TightID_pt_eta/abseta_pt_ratio"));
        h_muonPFiso2 = dynamic_cast<TH2F*>(
            muonIsoFile2->Get("TightISO_TightID_pt_eta/abseta_pt_ratio"));
        std::cout << "Got 2016 muon SFs!\n" << std::endl;
    }
}

Cuts::~Cuts()
{
    electronSFsFile->Close();
    electronRecoFile->Close();
    muonIDsFile1->Close();
    muonIsoFile1->Close();
    muonHltFile1->Close();
    if (is2016_)
    {
        muonHltFile2->Close();
        muonIDsFile2->Close();
        muonIsoFile2->Close();
    }
}

void Cuts::parse_config(const std::string confName)
{
    // Get the configuration file
    const YAML::Node config{YAML::LoadFile(confName)};

    if (config["trigLabel"])
    {
        cutConfTrigLabel_ = config["trigLabel"].as<std::string>();
    }
    if (config["plotPostfix"])
    {
        postfixName_ = config["plotPostfix"].as<std::string>();
    }

    const YAML::Node cuts{config["cuts"]};

    const YAML::Node tight_eles{cuts["tightElectrons"]};
    tightElePt_ = tight_eles["pt"].as<double>();
    tightElePtLeading_ = tight_eles["ptLeading"].as<double>();
    tightEleEta_ = tight_eles["eta"].as<double>();
    tightEleRelIso_ = tight_eles["relIso"].as<double>();
    numTightEle_ = tight_eles["number"].as<unsigned>();

    const YAML::Node loose_eles{cuts["looseElectrons"]};
    looseElePt_ = loose_eles["pt"].as<double>();
    looseElePtLeading_ = loose_eles["ptLeading"].as<double>();
    looseEleEta_ = loose_eles["eta"].as<double>();
    looseEleRelIso_ = loose_eles["relIso"].as<double>();
    numLooseEle_ = loose_eles["number"].as<unsigned>();

    const YAML::Node tight_mus{cuts["tightMuons"]};
    tightMuonPt_ = tight_mus["pt"].as<double>();
    tightMuonPtLeading_ = tight_mus["ptLeading"].as<double>();
    tightMuonEta_ = tight_mus["eta"].as<double>();
    tightMuonRelIso_ = tight_mus["relIso"].as<double>();
    numTightMu_ = tight_mus["number"].as<unsigned>();

    const YAML::Node loose_mus{cuts["looseMuons"]};
    looseMuonPt_ = loose_mus["pt"].as<double>();
    looseMuonPtLeading_ = loose_mus["ptLeading"].as<double>();
    looseMuonEta_ = loose_mus["eta"].as<double>();
    looseMuonRelIso_ = loose_mus["relIso"].as<double>();
    numLooseMu_ = loose_mus["number"].as<unsigned>();

    const YAML::Node jets{cuts["jets"]};
    jetPt_ = jets["pt"].as<double>();
    jetEta_ = jets["eta"].as<double>();
    numJets_ = jets["numJets"].as<unsigned>();
    maxJets_ = jets["maxJets"].as<unsigned>();
    numbJets_ = jets["numbJets"].as<unsigned>();
    maxbJets_ = jets["maxbJets"].as<unsigned>();
    maxbJetEta_ = jets["maxbJetEta"].as<double>();
    // numcJets_ = jets["numcJets"].as<unsigned>();

    std::cerr << "And so it's looking for " << numTightMu_ << " muons and "
              << numTightEle_ << " electrons" << std::endl;
}

bool Cuts::makeCuts(AnalysisEvent& event,
                    double& eventWeight,
                    std::map<std::string, std::shared_ptr<Plots>>& plotMap,
                    TH1D& cutFlow,
                    const int systToRun)
{
    if (!skipTrigger_)
    {
        if (!triggerCuts(event, eventWeight, systToRun))
        {
            return false; // Do trigger cuts
        }
    }

    if (!metFilters(event))
    {
        return false;
    }

    // Make lepton cuts. If the trigLabel contains d, we are in the ttbar CR
    // so the Z mass cut is skipped
    if (!makeLeptonCuts(event, eventWeight, plotMap, cutFlow, systToRun))
    {
        return false;
    }

    std::tie(event.jetIndex, event.jetSmearValue) =
        makeJetCuts(event, systToRun, eventWeight, true);

    if (event.jetIndex.size() < numJets_)
    {
        return false;
    }
    if (event.jetIndex.size() > maxJets_)
    {
        return false;
    }

    event.bTagIndex = makeBCuts(event, event.jetIndex, systToRun);

    if (doPlots_ || fillCutFlow_)
    {
        cutFlow.Fill(2.5, eventWeight);
    }
    if (doPlots_)
    {
        plotMap["jetSel"]->fillAllPlots(event, eventWeight);
    }

    if (event.bTagIndex.size() < numbJets_)
    {
        return false;
    }
    if (event.bTagIndex.size() > maxbJets_)
    {
        return false;
    }
    if (doPlots_)
    {
        plotMap["bTag"]->fillAllPlots(event, eventWeight);
    }
    if (doPlots_ || fillCutFlow_)
    {
        cutFlow.Fill(3.5, eventWeight);
    }

    // Do wMass stuff
    double invWmass{0.};
    invWmass = getWbosonQuarksCand(event, event.jetIndex, systToRun);

    // Debug chi2 cut
    //   double topMass = getTopMass(event);
    //   double topTerm = ( topMass-173.21 )/30.0;
    //   double wTerm = ( (event.wPairQuarks.first +
    //   event.wPairQuarks.second).M() - 80.3585 )/8.0;

    //   double chi2 = topTerm*topTerm + wTerm*wTerm;
    //   if ( chi2 < 2.0 && chi2 > 7.0 ) return false; // control region
    //   if ( chi2 >= 2.0 ) return false; //signal region

    // Signal Region W mass cut
    if (!isZplusCR_)
    {
        if (std::abs(invWmass) > invWMassCut_)
        {
            return false;
        }
    }
    // Z+jets Control Region
    else
    {
        if (std::abs(invWmass) <= invWMassCut_)
        {
            return false;
        }
        if (event.metPF2PATEt >= metDileptonCut_)
        {
            return false;
        }
    }

    if (doPlots_)
    {
        plotMap["wMass"]->fillAllPlots(event, eventWeight);
    }
    if (doPlots_ || fillCutFlow_)
    {
        cutFlow.Fill(4.5, eventWeight);
    }

    return true;
}

std::vector<double> Cuts::getRochesterSFs(const AnalysisEvent& event) const
{
    std::vector<double> SFs{};

    for (auto muonIt = event.muonIndexTight.begin();
         muonIt != event.muonIndexTight.end();
         muonIt++)
    {
        double tempSF{1.0};
        if (isMC_)
        {
            if (event.genMuonPF2PATPT[*muonIt] > 0) // matched gen muon
            {
                tempSF = rc_.kSpreadMC(event.muonPF2PATCharge[*muonIt],
                                       event.muonPF2PATPt[*muonIt],
                                       event.muonPF2PATEta[*muonIt],
                                       event.muonPF2PATPhi[*muonIt],
                                       event.genMuonPF2PATPT[*muonIt]);
            }
            else
            {
                static std::uniform_real_distribution<> u{0, 1};

                // We need a uniformly distributed "random" number, but this
                // should be the same every time, e.g. when we are looking at
                // systematics. So we will seed the random number generator
                // with a hash combining the properties of the muon (and make
                // the hopefully safe assumption no two muons have EXACTLY
                // the same properties within the same event)
                size_t seed{0};
                boost::hash_combine(seed, event.muonPF2PATCharge[*muonIt]);
                boost::hash_combine(seed, event.muonPF2PATPt[*muonIt]);
                boost::hash_combine(seed, event.muonPF2PATEta[*muonIt]);
                boost::hash_combine(seed, event.muonPF2PATPhi[*muonIt]);
                boost::hash_combine(
                    seed, event.muonPF2PATTkLysWithMeasurements[*muonIt]);
                boost::hash_combine(seed, event.eventNum);

                std::mt19937 gen(seed);

                tempSF =
                    rc_.kSmearMC(event.muonPF2PATCharge[*muonIt],
                                 event.muonPF2PATPt[*muonIt],
                                 event.muonPF2PATEta[*muonIt],
                                 event.muonPF2PATPhi[*muonIt],
                                 event.muonPF2PATTkLysWithMeasurements[*muonIt],
                                 u(gen));
            }
        }
        else
        {
            tempSF = rc_.kScaleDT(event.muonPF2PATCharge[*muonIt],
                                  event.muonPF2PATPt[*muonIt],
                                  event.muonPF2PATEta[*muonIt],
                                  event.muonPF2PATPhi[*muonIt]);
        }
        SFs.emplace_back(tempSF);
    }

    return SFs;
}

// Make lepton cuts. Will become customisable in a config later on.
bool Cuts::makeLeptonCuts(
    AnalysisEvent& event,
    double& eventWeight,
    std::map<std::string, std::shared_ptr<Plots>>& plotMap,
    TH1D& cutFlow,
    const int syst,
    const bool skipZCut)
{
    ////Do lepton selection.

    event.electronIndexTight = getTightEles(event);
    if (event.electronIndexTight.size() != numTightEle_)
    {
        return false;
    }
    event.electronIndexLoose = getLooseEles(event);
    if (event.electronIndexLoose.size() != numLooseEle_)
    {
        return false;
    }

    event.muonIndexTight = getTightMuons(event);
    if (event.muonIndexTight.size() != numTightMu_)
    {
        return false;
    }
    event.muonIndexLoose = getLooseMuons(event);
    if (event.muonIndexLoose.size() != numLooseMu_)
    {
        return false;
    }

    // If making NPL shape postLepSkim, MC leptons must BOTH be prompt
    if (isNPL_ && numTightEle_ == 2 && isMC_)
    {
        eventWeight *= -1.0;
        if (!event.genElePF2PATPromptFinalState[event.zPairIndex.first])
        {
            return false;
        }
        if (!event.genElePF2PATPromptFinalState[event.zPairIndex.second])
        {
            return false;
        }
    }

    if (isNPL_ && numTightMu_ == 2 && isMC_)
    { // if mumu channel
        eventWeight *= -1.0;
        if (!event.genMuonPF2PATPromptFinalState[event.zPairIndex.first])
        {
            return false;
        }
        if (!event.genMuonPF2PATPromptFinalState[event.zPairIndex.second])
        {
            return false;
        }
    }

    if (isNPL_ && numTightEle_ == 1 && numTightMu_ == 1 && isMC_)
    { // if emu channel
        eventWeight *= -1.0;
        if (!event.genElePF2PATPromptFinalState[event.zPairIndex.first])
        {
            return false;
        }
        if (!event.genMuonPF2PATPromptFinalState[event.zPairIndex.second])
        {
            return false;
        }
    }

    // This is to make some skims for faster running. Do lepSel and save some
    // files.
    if (postLepSelTree_)
    {
        postLepSelTree_->Fill();
    }

    event.muonMomentumSF = getRochesterSFs(event);

    if (!getDileptonZCand(
            event, event.electronIndexTight, event.muonIndexTight))
    {
        return false;
    }

    eventWeight *= getLeptonWeight(event, syst);

    if (doPlots_ || fillCutFlow_)
    {
        std::tie(event.jetIndex, event.jetSmearValue) =
            makeJetCuts(event, syst, eventWeight, false);
    }
    if (doPlots_)
    {
        plotMap["lepSel"]->fillAllPlots(event, eventWeight);
    }
    if (doPlots_ || fillCutFlow_)
    {
        cutFlow.Fill(0.5, eventWeight);
    }

    if (isNPL_)
    { // if is NPL channel
        double eeWeight{1.0};
        double mumuWeight{1.0};
        double emuWeight{1.0};

        if (invZMassCut_ == 20. && invWMassCut_ == 20.)
        {
            if (is2016_)
            {
                eeWeight = 0.956296241886;
                mumuWeight = 1.01737450337;
                emuWeight = 1.07476479834;
            }
            else
            {
                eeWeight = 1.57352187847;
                mumuWeight = 1.15322585438;
                emuWeight = 1.08021397164;
            }
        }
        if (invZMassCut_ == 20. && invWMassCut_ == 50.)
        {
            eeWeight = 1.12750771638;
            mumuWeight = 0.853155120216;
        }
        if (invZMassCut_ == 50. && invWMassCut_ == 50.)
        {
            eeWeight = 1.2334461839;
            mumuWeight = 0.997331838956;
        }
        if (numTightEle_ == 2)
        {
            eventWeight *= eeWeight;
        }
        if (numTightMu_ == 2)
        {
            eventWeight *= mumuWeight;
        }
        if (numTightEle_ == 1 && numTightMu_ == 1)
        {
            eventWeight *= emuWeight;
        }
    }

    if (std::abs((event.zPairLeptons.first + event.zPairLeptons.second).M()
                 - 91.1)
            > invZMassCut_
        && !skipZCut)
    {
        return false;
    }

    if (doPlots_ || fillCutFlow_)
    {
        std::tie(event.jetIndex, event.jetSmearValue) =
            makeJetCuts(event, syst, eventWeight, false);
    }
    if (doPlots_)
    {
        plotMap["zMass"]->fillAllPlots(event, eventWeight);
    }
    if (doPlots_ || fillCutFlow_)
    {
        cutFlow.Fill(1.5, eventWeight);
    }

    return true;
}

std::vector<int> Cuts::getTightEles(const AnalysisEvent& event) const
{
    std::vector<int> electrons;

    for (int i{0}; i < event.numElePF2PAT; i++)
    {
        if (!event.elePF2PATIsGsf[i])
            continue;

        if (electrons.size() < 1 && event.elePF2PATPT[i] <= tightElePtLeading_)
            continue;
        else if (electrons.size() >= 1 && event.elePF2PATPT[i] <= tightElePt_)
            continue;

        if (std::abs(event.elePF2PATSCEta[i]) > tightEleEta_)
            continue;

        // Ensure we aren't in the barrel/endcap gap and below the max safe eta
        // range
        if ((std::abs(event.elePF2PATSCEta[i]) > 1.4442
             && std::abs(event.elePF2PATSCEta[i]) < 1.566)
            || std::abs(event.elePF2PATSCEta[i]) > 2.50)
            continue;

        // VID cut
        if (event.elePF2PATCutIdTight[i] < 1)
            continue;

        // Cuts not part of the tuned ID
        if (std::abs(event.elePF2PATSCEta[i]) <= 1.479)
        {
            if (std::abs(event.elePF2PATD0PV[i]) >= 0.05)
                continue;
            if (std::abs(event.elePF2PATDZPV[i]) >= 0.10)
                continue;
        }
        else if (std::abs(event.elePF2PATSCEta[i]) > 1.479
                 && std::abs(event.elePF2PATSCEta[i]) < 2.50)
        {
            if (std::abs(event.elePF2PATD0PV[i]) >= 0.10)
                continue;
            if (std::abs(event.elePF2PATDZPV[i]) >= 0.20)
                continue;
        }
        electrons.emplace_back(i);
    }
    return electrons;
}

std::vector<int> Cuts::getLooseEles(const AnalysisEvent& event) const
{
    std::vector<int> electrons;
    for (int i{0}; i < event.numElePF2PAT; i++)
    {
        if (electrons.size() < 1 && event.elePF2PATPT[i] <= looseElePtLeading_)
            continue;
        else if (electrons.size() >= 1 && event.elePF2PATPT[i] <= looseElePt_)
            continue;
        if (std::abs(event.elePF2PATSCEta[i]) > tightEleEta_)
            continue;

        // Ensure we aren't in the barrel/endcap gap and below the max safe
        // eta range
        if ((std::abs(event.elePF2PATSCEta[i]) > 1.4442
             && std::abs(event.elePF2PATSCEta[i]) < 1.566)
            || std::abs(event.elePF2PATSCEta[i]) > 2.50)
            continue;

        // VID cut
        if (!event.elePF2PATCutIdVeto[i])
            continue;

        // Cuts not part of the tuned ID
        if (std::abs(event.elePF2PATSCEta[i]) <= 1.479)
        {
            if (std::abs(event.elePF2PATD0PV[i]) >= 0.05)
                continue;
            if (std::abs(event.elePF2PATDZPV[i]) >= 0.10)
                continue;
        }
        else if (std::abs(event.elePF2PATSCEta[i]) > 1.479
                 && std::abs(event.elePF2PATSCEta[i]) < 2.50)
        {
            if (std::abs(event.elePF2PATD0PV[i]) >= 0.10)
                continue;
            if (std::abs(event.elePF2PATDZPV[i]) >= 0.20)
                continue;
        }
        electrons.emplace_back(i);
    }
    return electrons;
}

std::vector<int> Cuts::getTightMuons(const AnalysisEvent& event) const
{
    std::vector<int> muons;
    if (is2016_)
    {
        for (int i{0}; i < event.numMuonPF2PAT; i++)
        {
            if (!event.muonPF2PATIsPFMuon[i])
                continue;

            if (muons.size() < 1
                && event.muonPF2PATPt[i] <= tightMuonPtLeading_)
                continue;
            else if (muons.size() >= 1 && event.muonPF2PATPt[i] <= tightMuonPt_)
                continue;

            if (std::abs(event.muonPF2PATEta[i]) >= tightMuonEta_)
                continue;
            if (event.muonPF2PATComRelIsodBeta[i] >= tightMuonRelIso_)
                continue;

            // Tight ID Cut
            if (!event.muonPF2PATTrackID[i])
                continue;
            if (!event.muonPF2PATGlobalID[i])
                continue;
            if (event.muonPF2PATGlbTkNormChi2[i] >= 10.)
                continue;
            if (event.muonPF2PATMatchedStations[i] < 2)
                continue; //
            if (std::abs(event.muonPF2PATDBPV[i]) >= 0.2)
                continue;
            if (std::abs(event.muonPF2PATDZPV[i]) >= 0.5)
                continue;
            if (event.muonPF2PATMuonNHits[i] < 1)
                continue;
            if (event.muonPF2PATVldPixHits[i] < 1)
                continue;
            if (event.muonPF2PATTkLysWithMeasurements[i] <= 5)
                continue;
            muons.emplace_back(i);
        }
    }
    else
    {
        for (int i{0}; i < event.numMuonPF2PAT; i++)
        {
            if (event.muonPF2PATIsPFMuon[i] && event.muonPF2PATTightCutId[i]
                && event.muonPF2PATPfIsoTight[i]
                && std::abs(event.muonPF2PATEta[i]) <= tightMuonEta_)
            {
                if (event.muonPF2PATPt[i]
                    >= (muons.empty() ? tightMuonPtLeading_ : tightMuonPt_))
                {
                    muons.emplace_back(i);
                }
            }
        }
    }
    return muons;
}

std::vector<int> Cuts::getLooseMuons(const AnalysisEvent& event) const
{
    std::vector<int> muons;
    if (is2016_)
    {
        for (int i{0}; i < event.numMuonPF2PAT; i++)
        {
            if (!event.muonPF2PATIsPFMuon[i])
                continue;

            if (muons.size() < 1
                && event.muonPF2PATPt[i] <= looseMuonPtLeading_)
                continue;
            else if (muons.size() >= 1 && event.muonPF2PATPt[i] <= looseMuonPt_)
                continue;

            if (std::abs(event.muonPF2PATEta[i]) >= looseMuonEta_)
                continue;
            if (event.muonPF2PATComRelIsodBeta[i] >= looseMuonRelIso_)
                continue;
            if (event.muonPF2PATGlobalID[i] || event.muonPF2PATTrackID[i])
                muons.emplace_back(i);
        }
    }
    else
    {
        for (int i{0}; i < event.numMuonPF2PAT; i++)
        {
            if (event.muonPF2PATIsPFMuon[i] && event.muonPF2PATLooseCutId[i]
                && event.muonPF2PATPfIsoLoose[i]
                && std::abs(event.muonPF2PATEta[i]) < looseMuonEta_)
            {
                if (event.muonPF2PATPt[i]
                    >= (muons.empty() ? looseMuonPtLeading_ : looseMuonPt_))
                {
                    muons.emplace_back(i);
                }
            }
        }
    }
    return muons;
}

bool Cuts::getDileptonZCand(AnalysisEvent& event,
                            const std::vector<int> electrons,
                            const std::vector<int> muons) const
{
    // Check if there are at least two electrons first. Otherwise use muons.

    if (electrons.size() == 2)
    {
        if (!invertLepCut_)
        {
            if (event.elePF2PATCharge[electrons[0]]
                    * event.elePF2PATCharge[electrons[1]]
                >= 0)
            {
                return false; // check electron pair have correct charge.
            }
        }
        else
        {
            if (!(event.elePF2PATCharge[electrons[0]]
                      * event.elePF2PATCharge[electrons[1]]
                  >= 0))
            {
                return false; // check electron pair have correct charge for
                              // same sign control region.
            }
        }

        const TLorentzVector lepton1{event.elePF2PATPX[electrons[0]],
                                     event.elePF2PATPY[electrons[0]],
                                     event.elePF2PATPZ[electrons[0]],
                                     event.elePF2PATE[electrons[0]]};
        const TLorentzVector lepton2{event.elePF2PATPX[electrons[1]],
                                     event.elePF2PATPY[electrons[1]],
                                     event.elePF2PATPZ[electrons[1]],
                                     event.elePF2PATE[electrons[1]]};

        if (lepton1.Pt() > lepton2.Pt())
        {
            event.zPairLeptons.first = lepton1;
            event.zPairIndex.first = electrons[0];

            event.zPairLeptons.second = lepton2;
            event.zPairIndex.second = electrons[1];
        }
        else
        {
            event.zPairLeptons.first = lepton2;
            event.zPairIndex.first = electrons[1];

            event.zPairLeptons.second = lepton1;
            event.zPairIndex.second = electrons[0];
        }

        event.zPairRelIso.first =
            event.elePF2PATComRelIsoRho[event.zPairIndex.first];
        event.zPairRelIso.second =
            event.elePF2PATComRelIsoRho[event.zPairIndex.second];

        return true;
    } // end electron if

    else if (muons.size() == 2)
    {
        if (!invertLepCut_)
        {
            if (event.muonPF2PATCharge[muons[0]]
                    * event.muonPF2PATCharge[muons[1]]
                >= 0)
            {
                return false;
            }
        }
        else
        {
            if (!(event.muonPF2PATCharge[muons[0]]
                      * event.muonPF2PATCharge[muons[1]]
                  >= 0))
            {
                return false;
            }
        }

        TLorentzVector lepton1{event.muonPF2PATPX[muons[0]],
                               event.muonPF2PATPY[muons[0]],
                               event.muonPF2PATPZ[muons[0]],
                               event.muonPF2PATE[muons[0]]};
        TLorentzVector lepton2{event.muonPF2PATPX[muons[1]],
                               event.muonPF2PATPY[muons[1]],
                               event.muonPF2PATPZ[muons[1]],
                               event.muonPF2PATE[muons[1]]};

        lepton1 *= event.muonMomentumSF.at(0);
        lepton2 *= event.muonMomentumSF.at(1);

        event.zPairLeptons.first =
            lepton1.Pt() > lepton2.Pt() ? lepton1 : lepton2;
        event.zPairIndex.first =
            lepton1.Pt() > lepton2.Pt() ? muons[0] : muons[1];
        event.zPairRelIso.first = event.muonPF2PATComRelIsodBeta[muons[0]];
        event.zPairLeptons.second =
            lepton1.Pt() > lepton2.Pt() ? lepton2 : lepton1;
        event.zPairRelIso.second = event.muonPF2PATComRelIsodBeta[muons[1]];
        event.zPairIndex.second =
            lepton1.Pt() > lepton2.Pt() ? muons[1] : muons[0];
        return true;
    }

    else if (electrons.size() == 1 && muons.size() == 1)
    {
        if (!invertLepCut_)
        {
            if (event.elePF2PATCharge[electrons[0]]
                    * event.muonPF2PATCharge[muons[1]]
                >= 0)
            {
                return false;
            }
        }
        else
        {
            if (!(event.elePF2PATCharge[electrons[0]]
                      * event.muonPF2PATCharge[muons[1]]
                  >= 0))
            {
                return false;
            }
        }

        TLorentzVector lepton1{event.elePF2PATPX[electrons[0]],
                               event.elePF2PATPY[electrons[0]],
                               event.elePF2PATPZ[electrons[0]],
                               event.elePF2PATE[electrons[0]]};
        TLorentzVector lepton2{event.muonPF2PATPX[muons[0]],
                               event.muonPF2PATPY[muons[0]],
                               event.muonPF2PATPZ[muons[0]],
                               event.muonPF2PATE[muons[0]]};

        lepton2 *= event.muonMomentumSF.at(0);

        event.zPairLeptons.first = lepton1;
        event.zPairLeptons.second = lepton2;
        return true;
    }
    else
    {
        return false; // Not dilepton candidate if this is the case ...
    }
}

double Cuts::getWbosonQuarksCand(AnalysisEvent& event,
                                 const std::vector<int> jets,
                                 const int syst) const
{
    auto closestWmass{std::numeric_limits<double>::infinity()};
    if (jets.size() > 2)
    {
        for (unsigned k{0}; k < jets.size(); k++)
        {
            for (unsigned l{k + 1}; l < jets.size(); l++)
            {
                // Now ensure that the leading b jet isn't one of these!
                if (event.bTagIndex.size() > 0)
                {
                    if (event
                            .jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                                [jets[k]]
                        > bDiscCut_)
                    {
                        if (event.jetIndex[event.bTagIndex[0]] == jets[k])
                            continue;
                    }
                    else if (
                        event
                            .jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                                [jets[l]]
                        > bDiscCut_)
                    {
                        if (event.jetIndex[event.bTagIndex[0]] == jets[l])
                            continue;
                    }
                }
                const TLorentzVector jetVec1{
                    getJetLVec(event, jets[k], syst, false).first};
                const TLorentzVector jetVec2{
                    getJetLVec(event, jets[l], syst, false).first};

                double invWbosonMass{(jetVec1 + jetVec2).M() - 80.385};

                if (std::abs(invWbosonMass) < std::abs(closestWmass))
                {
                    event.wPairQuarks.first =
                        jetVec1.Pt() > jetVec2.Pt() ? jetVec1 : jetVec2;
                    event.wPairIndex.first =
                        jetVec1.Pt() > jetVec2.Pt() ? jets[k] : jets[l];
                    event.wPairQuarks.second =
                        jetVec1.Pt() > jetVec2.Pt() ? jetVec2 : jetVec1;
                    event.wPairIndex.second =
                        jetVec1.Pt() > jetVec2.Pt() ? jets[l] : jets[k];
                    closestWmass = invWbosonMass;
                }
            }
        }
    }
    return closestWmass;
}

double Cuts::getTopMass(const AnalysisEvent& event) const
{
    TLorentzVector bVec(event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
    return (bVec + event.wPairQuarks.first + event.wPairQuarks.second).M();
}

std::pair<std::vector<int>, std::vector<double>>
    Cuts::makeJetCuts(const AnalysisEvent& event,
                      const int syst,
                      double& eventWeight,
                      const bool isProper) const
{
    std::vector<int> jets;
    std::vector<double> smears;

    double mcTag{1.};
    double mcNoTag{1.};
    double dataTag{1.};
    double dataNoTag{1.};
    // b-tagging errors
    double err1{0.};
    double err2{0.};
    double err3{0.};
    double err4{0.};

    for (int i{0}; i < event.numJetPF2PAT; i++)
    {
        auto [jetVec, smear] = getJetLVec(event, i, syst, true);
        smears.emplace_back(smear);

        if (jetVec.Pt() <= jetPt_ || jetVec.Eta() >= jetEta_)
        {
            continue;
        }

        bool jetId{true};

        if (jetIDDo_ && isProper)
        {
            if (is2016_)
            {
                // Jet ID == loose
                if (std::abs(jetVec.Eta()) <= 2.7)
                { // for cases where jet eta <= 2.7

                    // for all jets with eta <= 2.7
                    if (event.jetPF2PATNeutralHadronEnergyFraction[i] >= 0.99)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNeutralEmEnergyFraction[i] >= 0.99)
                    {
                        jetId = false;
                    }
                    if ((event.jetPF2PATChargedMultiplicity[i]
                         + event.jetPF2PATNeutralMultiplicity[i])
                        <= 1)
                    {
                        jetId = false;
                    }

                    // for jets with eta <= 2.40
                    if (std::abs(jetVec.Eta()) <= 2.40)
                    {
                        if (event.jetPF2PATChargedHadronEnergyFraction[i]
                            <= 0.0)
                        {
                            jetId = false;
                        }
                        if (event.jetPF2PATChargedMultiplicity[i] <= 0.0)
                        {
                            jetId = false;
                        }
                        if (event.jetPF2PATChargedEmEnergyFraction[i] >= 0.99)
                        {
                            jetId = false;
                        }
                    }
                }
                else if (std::abs(jetVec.Eta()) <= 3.0
                         && std::abs(jetVec.Eta()) > 2.70)
                {
                    if (event.jetPF2PATNeutralHadronEnergyFraction[i] >= 0.98)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNeutralEmEnergyFraction[i] <= 0.01)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNeutralMultiplicity[i] <= 2)
                    {
                        jetId = false;
                    }
                }
                else if (std::abs(jetVec.Eta()) > 3.0)
                { // for cases where jet eta > 3.0 and less than 5.0 (or max).
                    if (event.jetPF2PATNeutralEmEnergyFraction[i] >= 0.90)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNeutralMultiplicity[i] <= 10)
                    {
                        jetId = false;
                    }
                }
            }
            else
            {
                // Jet ID == tightLepVeto (loose is deprecated)
                // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
                if (std::abs(jetVec.Eta()) <= 2.7)
                { // for cases where jet eta <= 2.7

                    // for all jets with eta <= 2.7
                    if (event.jetPF2PATNeutralHadronEnergyFraction[i] >= 0.90)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNeutralEmEnergyFraction[i] >= 0.90)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNConstituents[i] <= 1)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATMuonFraction[i] >= 0.8)
                    {
                        jetId = false;
                    }

                    // for jets with eta <= 2.40
                    if (std::abs(jetVec.Eta()) <= 2.40)
                    {
                        if (event.jetPF2PATChargedHadronEnergyFraction[i]
                            <= 0.0)
                        {
                            jetId = false;
                        }
                        if (event.jetPF2PATChargedMultiplicity[i] <= 0.0)
                        {
                            jetId = false;
                        }
                        if (event.jetPF2PATChargedEmEnergyFraction[i] >= 0.8)
                        {
                            jetId = false;
                        }
                    }
                }
                else if (std::abs(jetVec.Eta()) <= 3.0
                         && std::abs(jetVec.Eta()) > 2.70)
                {
                    if (event.jetPF2PATNeutralEmEnergyFraction[i] <= 0.02
                        || event.jetPF2PATNeutralEmEnergyFraction[i] >= 0.99)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNeutralMultiplicity[i] <= 2)
                    {
                        jetId = false;
                    }
                }
                else if (std::abs(jetVec.Eta()) > 3.0)
                { // for cases where jet eta > 3.0 and less than 5.0 (or max).
                    if (event.jetPF2PATNeutralEmEnergyFraction[i] >= 0.90)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNeutralHadronEnergyFraction[i] <= 0.02)
                    {
                        jetId = false;
                    }
                    if (event.jetPF2PATNeutralMultiplicity[i] <= 10)
                    {
                        jetId = false;
                    }
                }
            }
        }

        if (!jetId)
        {
            continue;
        }

        const double deltaLep{std::min(deltaR(event.zPairLeptons.first.Eta(),
                                              event.zPairLeptons.first.Phi(),
                                              jetVec.Eta(),
                                              jetVec.Phi()),
                                       deltaR(event.zPairLeptons.second.Eta(),
                                              event.zPairLeptons.second.Phi(),
                                              jetVec.Eta(),
                                              jetVec.Phi()))};

        if (deltaLep < 0.4 && isProper)
        {
            continue;
        }

        //    if (event.jetPF2PATdRClosestLepton[i] < 0.5) continue;
        if (isMC_ && makeBTagEffPlots_ && isProper)
        {
            // Fill eff info here if needed.
            if (std::abs(event.jetPF2PATPID[i]) == 5)
            { // b-jets
                bTagEffPlots_[0]->Fill(jetVec.Pt(), std::abs(jetVec.Eta()));
                if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                        [i]
                    > bDiscCut_)
                {
                    bTagEffPlots_[4]->Fill(jetVec.Pt(), std::abs(jetVec.Eta()));
                }
            }
            if (std::abs(event.jetPF2PATPID[i]) == 4)
            { // charm
                bTagEffPlots_[1]->Fill(jetVec.Pt(), std::abs(jetVec.Eta()));
                if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                        [i]
                    > bDiscCut_)
                {
                    bTagEffPlots_[5]->Fill(jetVec.Pt(), std::abs(jetVec.Eta()));
                }
            }
            if (std::abs(event.jetPF2PATPID[i]) > 0
                && std::abs(event.jetPF2PATPID[i]) < 4)
            { // light jets
                bTagEffPlots_[2]->Fill(jetVec.Pt(), std::abs(jetVec.Eta()));
                if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                        [i]
                    > bDiscCut_)
                {
                    bTagEffPlots_[6]->Fill(jetVec.Pt(), std::abs(jetVec.Eta()));
                }
            }
            if (std::abs(event.jetPF2PATPID[i]) == 21)
            { // gluons
                bTagEffPlots_[3]->Fill(jetVec.Pt(), std::abs(jetVec.Eta()));
                if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                        [i]
                    > bDiscCut_)
                {
                    bTagEffPlots_[7]->Fill(jetVec.Pt(), std::abs(jetVec.Eta()));
                }
            }
        }

        jets.emplace_back(i);

        if (getBTagWeight_)
        {
            getBWeight(event,
                       jetVec,
                       i,
                       mcTag,
                       mcNoTag,
                       dataTag,
                       dataNoTag,
                       err1,
                       err2,
                       err3,
                       err4);
        }
    }
    // Evaluate b-tag weight for event here.
    if (getBTagWeight_ && isProper)
    {
        double bWeight{(dataNoTag * dataTag) / (mcNoTag * mcTag)};
        if (mcNoTag == 0 || mcTag == 0 || dataNoTag == 0 || dataTag == 0
            || mcNoTag != mcNoTag || mcTag != mcTag || dataTag != dataTag
            || dataNoTag != dataNoTag)
        {
            bWeight = 1.;
        }
        const double bWeightErr{
            std::sqrt(pow(err1 + err2, 2) + pow(err3 + err4, 2)) * bWeight};
        if (syst == 256)
        {
            bWeight += bWeightErr;
        }
        if (syst == 512)
        {
            bWeight -= bWeightErr;
        }

        eventWeight *= bWeight;
    }

    return {jets, smears};
}

std::vector<int> Cuts::makeBCuts(const AnalysisEvent& event,
                                 const std::vector<int> jets,
                                 const int syst) const
{
    std::vector<int> bJets;
    for (unsigned int i = 0; i < jets.size(); i++)
    {
        const TLorentzVector jetVec{
            getJetLVec(event, jets[i], syst, false).first};
        const float bDisc{
            is2016_
                ? event.jetPF2PATBDiscriminator[jets[i]]
                : event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                      [jets[i]]};

        if (bDisc <= bDiscCut_)
        {
            continue;
        }
        if (jetVec.Eta() >= maxbJetEta_)
        {
            continue;
        }
        bJets.emplace_back(i);
    }
    return bJets;
}

bool Cuts::triggerCuts(const AnalysisEvent& event,
                       double& eventWeight,
                       const int syst) const
{
    if (skipTrigger_)
    {
        return true;
    }

    // TRIGGER LOGIC

    // MuEG triggers
    // clang-format off
    const bool muEGTrig{event.muEGTrig()};

    // clang-format on

    // double electron triggers
    const bool eeTrig{event.eeTrig()};

    // double muon triggers
    const bool mumuTrig{event.mumuTrig()};

    // single electron triggers
    const bool eTrig{event.eTrig()};

    // single muon triggers
    const bool muTrig{event.muTrig()};

    std::string channel = "";

    // Dilepton channels
    if (cutConfTrigLabel_.find("e") != std::string::npos)
    {
        channel = "ee";
    }
    if (cutConfTrigLabel_.find("d") != std::string::npos)
    {
        channel = "emu";
    }
    if (cutConfTrigLabel_.find("m") != std::string::npos)
    {
        channel = "mumu";
    }

    double twgt{1.0};

    if (!is2016_ && isMC_)
    {
        // Lepton SFs obtained from triggerScaleFactorsAlgo
        if (channel == "ee")
        {
            if (eTrig || eeTrig)
            {
                twgt = 0.93106;
                if (syst == 1)
                {
                    twgt += 0.01;
                }
                else if (syst == 2)
                {
                    twgt -= 0.01;
                }
            }
        }
        else if (channel == "mumu")
        {
            if (muTrig || mumuTrig)
            {
                twgt = 0.97170;
                if (syst == 1)
                {
                    twgt += 0.01;
                }
                else if (syst == 2)
                {
                    twgt -= 0.01;
                }
            }
        }
        else if (channel == "emu")
        {
            if (muEGTrig)
            {
                twgt = 0.95350;
                if (syst == 1)
                {
                    twgt += 0.02;
                }
                else if (syst == 2)
                {
                    twgt -= 0.02;
                }
            }
        }
        else
        {
            throw std::runtime_error("Unknown channel");
        }
    }
    else if (is2016_ && isMC_)
    { // Apply SFs to MC if 2016
        // Dilepton channels
        if (channel == "ee")
        {
            if (eTrig || eeTrig)
            { // If singleElectron or doubleEG trigger fires ...
                twgt = 0.96917; // 0.97554 for data eff; 0.98715 for SF
                if (syst == 1)
                {
                    twgt += 0.01; // +-/ 0.00138 for eff; 0.00063 for SF
                }
                if (syst == 2)
                {
                    twgt -= 0.01;
                }
            }
        }
        else if (channel == "mumu")
        {
            if (muTrig || mumuTrig)
            { // If doubleMuon or singleMuon trigger fires ...

                // eff pre-HIP fix: 0.98069  +/- -0.00070/0.00073; eff post-HIP
                // fix: 0.99061 +/- -0.00057/0.00061;
                // SF pre-HIP fix 0.98703 +/- 0.00016 and 0.99843 +/- 0.00016
                // for post-HIP fix SF pre-HIP fix 0.98868 +/- 0.00013 and
                // 0.99868 +/- 0.00017  for post-HIP fix

                twgt = (0.97679 * lumiRunsBCDEF_ + 0.98941 * lumiRunsGH_)
                       / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06);

                if (syst == 1)
                {
                    twgt += 0.01;
                }
                if (syst == 2)
                {
                    twgt -= 0.01;
                }
            }
        }
        else if (channel == "emu")
        { // If MuonEG trigger fires, regardless of singleElectron/singleMuon
          // triggers
            if (muEGTrig)
            {
                twgt = 0.98710;
                if (syst == 1)
                {
                    twgt += 0.02; // -0.01220/0.01339 for eff; 0.01018 for SF
                }
                if (syst == 2)
                {
                    twgt -= 0.02;
                }
            }
        }
        else
        {
            std::cout << "Trigger not found!" << std::endl;
            twgt = 0.0; // Return 0.0 if trigger isn't found.
        }
    }

    // Check which trigger fired and if it correctly corresponds to the
    // channel being scanned over.

    if (channel == "emu")
    {
        if (muEGTrig)
        {
            if (isMC_)
            {
                eventWeight *= twgt; // trigger weight should be unchanged
                                     // for data anyway, but good practice to
                                     // explicitly not apply it.
            }
            return true;
        }
    }

    if (channel == "ee")
    {
        //    if ( eeTrig && !(muEGTrig || mumuTrig) ) { // Original trigger
        //    logic, for double triggers only
        if ((eeTrig || eTrig) && !(muEGTrig || mumuTrig || muTrig))
        {
            if (isMC_)
            {
                eventWeight *= twgt; // trigger weight should be unchanged
                                     // for data anyway, but good practice to
                                     // explicitly not apply it.
            }
            return true;
        }
    }

    if (channel == "mumu")
    {
        // Trigger logic for double + single triggers
        if ((mumuTrig || muTrig) && !(eeTrig || muEGTrig || eTrig))
        {
            if (isMC_)
            {
                eventWeight *= twgt; // trigger weight should be unchanged
                                     // for data anyway, but good practice to
                                     // explicitly not apply it.
            }
            return true;
        }
    }
    return false;
}

// Does event pass MET Filter
bool Cuts::metFilters(const AnalysisEvent& event) const
{
    if (event.Flag_HBHENoiseFilter <= 0 || event.Flag_HBHENoiseIsoFilter <= 0
        || event.Flag_globalTightHalo2016Filter <= 0
        || event.Flag_EcalDeadCellTriggerPrimitiveFilter <= 0
        || event.Flag_goodVertices <= 0
        || (!isMC_ && event.Flag_eeBadScFilter <= 0))
    {
        return false;
    }

    if (is2016_
        && (event.Flag_ecalLaserCorrFilter <= 0
            || event.Flag_chargedHadronTrackResolutionFilter <= 0
            || event.Flag_muonBadTrackFilter <= 0
            || (!isMC_ && event.Flag_noBadMuons <= 0)))
    {
        return false;
    }

    if (!is2016_
        && (event.Flag_BadPFMuonFilter <= 0
            || event.Flag_BadChargedCandidateFilter <= 0
            || event.Flag_ecalBadCalibFilter <= 0))
    {
        return false;
    }

    return true;
}

double Cuts::deltaPhi(const double phi1, const double phi2)
{
    return std::atan2(std::sin(phi1 - phi2), std::cos(phi1 - phi2));
}

double Cuts::deltaR(const double eta1,
                    const double phi1,
                    const double eta2,
                    const double phi2)
{
    return std::sqrt(std::pow(eta1 - eta2, 2)
                     + std::pow(deltaPhi(phi1, phi2), 2));
}

double Cuts::getLeptonWeight(const AnalysisEvent& event, const int syst) const
{
    // If number of electrons is > 1  then both z pair are electrons, so get
    // their weight
    if (!isMC_)
    {
        return 1.;
    }

    double leptonWeight{1.};

    if (numTightEle_ == 2)
    {
        leptonWeight *= eleSF(event.zPairLeptons.first.Pt(),
                              event.elePF2PATSCEta[event.zPairIndex.first],
                              syst);
        leptonWeight *= eleSF(event.zPairLeptons.second.Pt(),
                              event.elePF2PATSCEta[event.zPairIndex.second],
                              syst);
    }

    else if (numTightMu_ == 2)
    {
        leptonWeight *= muonSF(event.zPairLeptons.first.Pt(),
                               event.zPairLeptons.first.Eta(),
                               syst);
        leptonWeight *= muonSF(event.zPairLeptons.second.Pt(),
                               event.zPairLeptons.second.Eta(),
                               syst);
    }
    else if (numTightEle_ == 1 && numTightMu_ == 1)
    {
        leptonWeight *= eleSF(event.elePF2PATPT[event.electronIndexTight[0]],
                              event.elePF2PATSCEta[event.electronIndexTight[0]],
                              syst);
        leptonWeight *= muonSF(event.muonPF2PATPt[event.muonIndexTight[0]],
                               event.muonPF2PATEta[event.muonIndexTight[0]],
                               syst);
    }
    return leptonWeight;
}

double Cuts::eleSF(const double pt, const double eta, const int syst) const
{
    const double maxPt{h_eleSFs->GetYaxis()->GetXmax() - 0.1};
    const double minRecoPt{h_eleReco->GetYaxis()->GetXmin() + 0.1};
    int bin1{0};
    int bin2{0};

    // If cut-based, std::abs eta, else just eta
    if (pt <= maxPt)
    {
        bin1 = h_eleSFs->FindBin(eta, pt);
        // Reco SF has a min to consider
        if (pt > minRecoPt)
        {
            bin2 = h_eleReco->FindBin(eta, pt);
        }
        else
        {
            bin2 = h_eleReco->FindBin(eta, minRecoPt);
        }
    }
    else
    {
        bin1 = h_eleSFs->FindBin(eta, maxPt);
        bin2 = h_eleReco->FindBin(eta, maxPt);
    }

    double eleIdSF{h_eleSFs->GetBinContent(bin1)};
    double eleRecoSF{h_eleReco->GetBinContent(bin2)};

    if (syst == 1)
    {
        eleIdSF += h_eleSFs->GetBinError(bin1);
        eleRecoSF += h_eleReco->GetBinError(bin2);
        if (pt > 80.0 || pt <= 20.0)
        {
            eleRecoSF += 0.01;
        }
    }

    if (syst == 2)
    {
        eleIdSF -= h_eleSFs->GetBinError(bin1);
        eleRecoSF -= h_eleReco->GetBinError(bin2);
        if (pt > 80.0 || pt <= 20.0)
        {
            eleRecoSF -= 0.01;
        }
    }

    return eleIdSF * eleRecoSF;
}

double Cuts::muonSF(const double pt, const double eta, const int syst) const
{
    if (!is2016_)
    {
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017

        static constexpr std::array<double, 3> etaBinEdges{0.9, 1.2, 2.1};
        static constexpr std::array<double, 5> ptBinEdges{25, 30, 40, 50, 60};
        // clang-format off
        static constexpr std::array<std::array<std::pair<double, double>, 6>, 4> muonIDSFs{{
                {{{0.9910777627756951, 0.0034967203087024274}, {0.9874104682620840, 0.0018975082960536634}, {0.9907753279135898, 0.0003977503197109705}, {0.9892483588952047, 0.00032329941312374114}, {0.9855545160334763, 0.0008602730194340781}, {0.9898057377093389, 0.0016204327419630327}}},  // η 0-0.9
                {{{0.9927389275515244, 0.0047437370661283660}, {0.9850639397625120, 0.0218787608424192500}, {0.9865359464182247, 0.0006903526352667042}, {0.9849130931014930, 0.02013503091494561000}, {0.9839056384760008, 0.0015917233692683600}, {0.9840604031434680, 0.0121278780491192190}}},  // η 0.9-1.2
                {{{0.9924252719877384, 0.0077859527402420020}, {0.9890884461284933, 0.0149053560477533670}, {0.9946469069883841, 0.0124242236899199730}, {0.9926528825155183, 0.00993167811549696700}, {0.9906364222943529, 0.0009713213798502730}, {0.9920464322143979, 0.0021353964567237746}}},  // η 1.2-2.1111
                {{{0.9758095839531763, 0.0043993151217841040}, {0.9745153594179884, 0.0027111009825340473}, {0.9787410500158746, 0.0010035577872014160}, {0.9781891229195010, 0.00112306059413853970}, {0.9673568416097894, 0.0037006525169638958}, {0.9766311856731202, 0.0086266646688285500}}},  // η 2.1+
            }};
        static constexpr std::array<std::array<std::pair<double, double>, 6>, 4> muonIsoSFs{{
                {{{0.9931118153898082f, 0.0023961047801277897f}, {0.9963959301142516f, 0.0011130262628559282f}, {0.9983181988469478f, 0.00028317863248081760f}, {0.9994347636372417f, 0.00011232714693801277f}, {0.9997680258722230f, 0.00021161394350289537f}, {1.0002119715073154f, 0.00034548765287024270f}}}, // η 0-0.9
                {{{0.9956674008369976f, 0.0040757770757467080f}, {0.9932947173775393f, 0.0020402913397330540f}, {0.9976914943066302f, 0.00052908143248948090f}, {0.9990342041383322f, 0.00020372497501125323f}, {0.9994155529201356f, 0.00041662994217661376f}, {1.0002046210970306f, 0.00066202613921456840f}}}, // η 0.9-1.2
                {{{0.9967414270112102f, 0.0016604845648881112f}, {0.9988489100332861f, 0.0009332375024769512f}, {0.9988826842051872f, 0.00028173209377228390f}, {0.9993872934568914f, 0.00011049799288742951f}, {0.9997410091519127f, 0.00024619498190323770f}, {1.0004725402685148f, 0.00041307624612303565f}}}, // η 1.2-2.1
                {{{0.9973954140213298f, 0.0021903221583127706f}, {0.9987560726170641f, 0.0012087380293472640f}, {0.9990618844158636f, 0.00037710084052130214f}, {0.9997635755214144f, 0.00017351798487330648f}, {1.0002224795137067f, 0.00051831556143466900f}, {0.9993431865975091f, 0.00085142066151194850f}}}, // η 2.1+
            }};
        // clang-format on

        const auto etaBin{std::distance(
            etaBinEdges.begin(),
            std::upper_bound(etaBinEdges.begin(), etaBinEdges.end(), eta))};
        const auto ptBin{std::distance(
            ptBinEdges.begin(),
            std::upper_bound(ptBinEdges.begin(), ptBinEdges.end(), pt))};

        const auto muonIDSF{muonIDSFs.at(etaBin).at(ptBin)};
        const auto muonIsoSF{muonIsoSFs.at(etaBin).at(ptBin)};

        if (syst == 1)
        {
            return (muonIDSF.first + muonIDSF.second)
                   * (muonIsoSF.first + muonIsoSF.second);
        }
        if (syst == 2)
        {
            return (muonIDSF.first - muonIDSF.second)
                   * (muonIsoSF.first - muonIsoSF.second);
        }
        else
        {
            return muonIDSF.first * muonIsoSF.first;
        }
    }
    else
    { // Run2016 needs separate treatments in pre and post HIP eras

        double maxIdPt{h_muonIDs1->GetYaxis()->GetXmax() - 0.1};
        double maxIsoPt{h_muonPFiso1->GetYaxis()->GetXmax() - 0.1};
        double minIdPt{h_muonIDs1->GetYaxis()->GetXmin() + 0.1};
        double minIsoPt{h_muonPFiso1->GetYaxis()->GetXmin() + 0.1};

        int binId1{0}, binIso1{0};
        int binId2{0}, binIso2{0};

        if (pt > maxIdPt)
        {
            binId1 = h_muonIDs1->FindBin(std::abs(eta), maxIdPt);
        }
        else if (pt < minIdPt)
        {
            binId1 = h_muonIDs1->FindBin(std::abs(eta), minIdPt);
        }
        else
        {
            binId1 = h_muonIDs1->FindBin(std::abs(eta), pt);
        }
        if (pt > maxIsoPt)
        {
            binIso1 = h_muonPFiso1->FindBin(std::abs(eta), maxIsoPt);
        }
        else if (pt < minIsoPt)
        {
            binIso1 = h_muonPFiso1->FindBin(std::abs(eta), minIsoPt);
        }
        else
        {
            binIso1 = h_muonPFiso1->FindBin(std::abs(eta), pt);
        }

        if (pt > maxIdPt)
        {
            binId2 = h_muonIDs2->FindBin(std::abs(eta), maxIdPt);
        }
        else if (pt < minIdPt)
        {
            binId2 = h_muonIDs2->FindBin(std::abs(eta), minIdPt);
        }
        else
        {
            binId2 = h_muonIDs2->FindBin(std::abs(eta), pt);
        }
        if (pt > maxIsoPt)
        {
            binIso2 = h_muonPFiso2->FindBin(std::abs(eta), maxIsoPt);
        }
        else if (pt < minIsoPt)
        {
            binIso2 = h_muonPFiso2->FindBin(std::abs(eta), minIsoPt);
        }
        else
        {
            binIso2 = h_muonPFiso2->FindBin(std::abs(eta), pt);
        }

        double muonIdSF{1.0};
        double muonPFisoSF{1.0};
        muonIdSF = (h_muonIDs1->GetBinContent(binId1) * lumiRunsBCDEF_
                    + h_muonIDs2->GetBinContent(binId2) * lumiRunsGH_)
                   / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06);
        muonPFisoSF = (h_muonPFiso1->GetBinContent(binIso1) * lumiRunsBCDEF_
                       + h_muonPFiso2->GetBinContent(binIso2) * lumiRunsGH_)
                      / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06);

        if (syst == 1)
        {
            muonIdSF += (h_muonIDs1->GetBinError(binId1) * lumiRunsBCDEF_
                         + h_muonIDs2->GetBinError(binId2) * lumiRunsGH_)
                            / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06)
                        + 0.01; // Additional 1% uncert for ID and 0.5% for
                                // iso as recommended
            muonPFisoSF += (h_muonPFiso1->GetBinError(binIso1) * lumiRunsBCDEF_
                            + h_muonIDs2->GetBinError(binId2) * lumiRunsGH_)
                               / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06)
                           + 0.005;

            return muonIdSF * muonPFisoSF;
        }
        else if (syst == 2)
        {
            muonIdSF -= (h_muonIDs1->GetBinError(binId1) * lumiRunsBCDEF_
                         + h_muonIDs2->GetBinError(binId2) * lumiRunsGH_)
                            / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06)
                        - 0.01; // Additional 1% uncert for ID and 0.5% for
                                // iso as recommended
            muonPFisoSF -= (h_muonPFiso1->GetBinError(binIso1) * lumiRunsBCDEF_
                            + h_muonIDs2->GetBinError(binId2) * lumiRunsGH_)
                               / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06)
                           - 0.005;
            return muonIdSF * muonPFisoSF;
        }
        else
        {
            return muonIdSF * muonPFisoSF;
        }
    }
}

void Cuts::initialiseJECCors()
{
    std::ifstream jecFile;
    if (!is2016_)
    {
        jecFile.open("scaleFactors/2017/"
                     "Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt",
                     std::ifstream::in);
    }
    else
    {
        jecFile.open("scaleFactors/2016/"
                     "Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt",
                     std::ifstream::in);
    }
    std::string line;
    bool first{true};

    if (!jecFile.is_open())
    {
        std::cout << "Unable to open jecFile." << std::endl;
        exit(0);
    }

    while (getline(jecFile, line))
    {
        std::vector<std::string> tempVec;
        std::stringstream lineStream{line};
        std::string item;
        while (std::getline(lineStream, item, ' '))
        {
            tempVec.emplace_back(item);
        }
        std::vector<double> tempUp;
        std::vector<double> tempDown;

        etaMinJEC_.emplace_back(std::stof(tempVec[0]));
        etaMaxJEC_.emplace_back(std::stof(tempVec[1]));
        for (unsigned i{1}; i < tempVec.size() / 3; i++)
        {
            unsigned ind{i * 3};
            if (first)
            {
                ptMinJEC_.emplace_back(std::stof(tempVec[ind]));
                ptMaxJEC_.emplace_back((ind + 3 >= tempVec.size()
                                            ? 10000.
                                            : std::stof(tempVec[ind + 3])));
            }
            tempUp.emplace_back(std::stof(tempVec[ind + 1]));
            tempDown.emplace_back(std::stof(tempVec[ind + 2]));
        }
        jecSFUp_.emplace_back(tempUp);
        jecSFDown_.emplace_back(tempDown);
        first = false;
    }
}

double Cuts::getJECUncertainty(const double pt,
                               const double eta,
                               const int syst) const
{
    if (!(syst == 4 || syst == 8))
    {
        return 0.;
    }
    unsigned ptBin{0};
    unsigned etaBin{0};
    for (unsigned i{0}; i < ptMinJEC_.size(); i++)
    {
        if (pt > ptMinJEC_[i] && pt < ptMaxJEC_[i])
        {
            ptBin = i;
            break;
        }
    }
    for (unsigned i{0}; i < etaMinJEC_.size(); i++)
    {
        if (eta > etaMinJEC_[i] && eta < etaMaxJEC_[i])
        {
            etaBin = i;
            break;
        }
    }

    const double lowFact{syst == 4 ? jecSFUp_[etaBin][ptBin]
                                   : jecSFDown_[etaBin][ptBin]};
    const double hiFact{syst == 4 ? jecSFUp_[etaBin][ptBin + 1]
                                  : jecSFDown_[etaBin][ptBin + 1]};

    // Now do some interpolation
    const double a{(hiFact - lowFact) / (ptMaxJEC_[ptBin] - ptMinJEC_[ptBin])};
    const double b{(lowFact * (ptMaxJEC_[ptBin]) - hiFact * ptMinJEC_[ptBin])
                   / (ptMaxJEC_[ptBin] - ptMinJEC_[ptBin])};
    return (syst == 4 ? a * pt + b : -(a * pt + b));
}

std::pair<TLorentzVector, double> Cuts::getJetLVec(const AnalysisEvent& event,
                                                   const int index,
                                                   const int syst,
                                                   const bool initialRun) const
{
    static constexpr double MIN_JET_ENERGY{1e-2};
    TLorentzVector returnJet;
    double newSmearValue{1.0};

    if (!initialRun)
    {
        newSmearValue = event.jetSmearValue.at(index);
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);
        returnJet *= newSmearValue;

        if (isMC_)
        {
            double jerUncer{
                getJECUncertainty(returnJet.Pt(), returnJet.Eta(), syst)};
            returnJet *= 1 + jerUncer;
        }

        return {returnJet, newSmearValue};
    }

    if (!isMC_)
    {
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);
        return {returnJet, newSmearValue};
    }

    // Not initial run, smear must be calculated

    // TODO: Check this is correct
    // For now, just leave jets of too large/small pT, large rho, or large η
    // untouched
    const double rho{is2016_ ? event.elePF2PATRhoIso[0]
                             : event.fixedGridRhoFastjetAll};
    if (event.jetPF2PATPtRaw[index] < 15 || event.jetPF2PATPtRaw[index] > 3000
        || rho > (is2016_ ? 40.9 : 42.52)
        || std::abs(event.jetPF2PATEta[index]) > 4.7)
    {
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);
        return {returnJet, newSmearValue};
    }

    // TODO: Should this be gen or reco level?
    // I think reco because gen might not exist? (does not exist when
    // smearing)
    const double ptRes{is2016_ ? jet2016PtSimRes(event.jetPF2PATPtRaw[index],
                                                 event.jetPF2PATEta[index],
                                                 rho)
                               : jet2017PtSimRes(event.jetPF2PATPtRaw[index],
                                                 event.jetPF2PATEta[index],
                                                 rho)};

    auto [jerSF, jerSigma] =
        is2016_ ? jet2016SFs(std::abs(event.jetPF2PATEta[index]))
                : jet2017SFs(std::abs(event.jetPF2PATEta[index]));

    if (syst == 16)
    {
        jerSF += jerSigma;
    }
    else if (syst == 32)
    {
        jerSF -= jerSigma;
    }

    std::optional<size_t> matchingGenIndex{std::nullopt};
    for (size_t genIndex{0}; genIndex < event.NJETSMAX; ++genIndex)
    {
        const double dR{deltaR(event.genJetPF2PATEta[genIndex], event.genJetPF2PATPhi[genIndex], event.jetPF2PATEta[index], event.jetPF2PATPhi[index])};
        const double dPt{event.jetPF2PATPtRaw[index] - event.genJetPF2PATPT[genIndex]};

        if (event.genJetPF2PATPT[genIndex] > 0 && dR < (0.4 / 2.0)
                && std::abs(dPt) < 3.0 * ptRes * event.jetPF2PATPtRaw[index])
        {
            matchingGenIndex = genIndex;
            break;
        }
    }

    if (matchingGenIndex.has_value())
    // If matching from GEN to RECO using dR<Rcone/2 and dPt < 3*sigma,
    // just scale
    {
        const double dPt{event.jetPF2PATPtRaw[index] - event.genJetPF2PATPT[matchingGenIndex.value()]};
        newSmearValue =
            std::max(1. + (jerSF - 1.) * dPt / event.jetPF2PATPtRaw[index], 0.);
    }
    else // If not matched to a gen jet, randomly smear
    {
        std::normal_distribution<> d(
            0, ptRes * std::sqrt(std::max(jerSF * jerSF - 1, 0.)));

        // Like with the Rochester corrections, seed the random number
        // generator with event (jet) properties so that each jet is smeared
        // the same way every time it is processed
        size_t seed{0};
        boost::hash_combine(seed, event.jetPF2PATPtRaw[index]);
        boost::hash_combine(seed, event.jetPF2PATEta[index]);
        boost::hash_combine(seed, event.jetPF2PATPhi[index]);
        boost::hash_combine(seed, event.eventNum);
        std::mt19937 gen(seed);

        newSmearValue = 1.0 + d(gen);
    }

    if (event.jetPF2PATE[index] * newSmearValue < MIN_JET_ENERGY)
    // Negative or too small scale factor
    {
        newSmearValue = MIN_JET_ENERGY / event.jetPF2PATE[index];
    }

    returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                         event.jetPF2PATPy[index],
                         event.jetPF2PATPz[index],
                         event.jetPF2PATE[index]);
    returnJet *= newSmearValue;

    if (isMC_)
    {
        double jerUncer{
            getJECUncertainty(returnJet.Pt(), returnJet.Eta(), syst)};
        returnJet *= 1 + jerUncer;
    }

    return {returnJet, newSmearValue};
}

double
    Cuts::jet2016PtSimRes(const double pt, const double eta, const double rho)
{
    if (pt < 15 || pt > 3000)
    {
        throw std::runtime_error("pT " + std::to_string(pt)
                                 + " out of range to assign resolution");
    }

    static constexpr std::array<double, 14> etaBinEdges{
        0, 0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3, 3.2, 4.7};
    static constexpr std::array<double, 8> rhoBinEdges{
        0, 6.69, 12.39, 18.09, 23.79, 29.49, 35.19, 40.9};
    const auto res = [pt](const double p0,
                          const double p1,
                          const double p2,
                          const double p3) {
        return (
            sqrt(p0 * abs(p0) / (pt * pt) + p1 * p1 * pow(pt, p3) + p2 * p2));
    };
    const auto etaBin{std::distance(etaBinEdges.begin(),
                                    std::upper_bound(etaBinEdges.begin(),
                                                     etaBinEdges.end(),
                                                     std::abs(eta)))};
    const auto rhoBin{std::distance(
        rhoBinEdges.begin(),
        std::upper_bound(rhoBinEdges.begin(), rhoBinEdges.end(), rho))};

    // https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Summer16_25nsV1_MC/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt
    switch (etaBin)
    {
        case 1:
            switch (rhoBin)
            {
                case 1: return res(0.6172, 0.3908, 0.02003, -0.6407);
                case 2: return res(1.775, 0.4231, 0.02199, -0.6701);
                case 3: return res(2.457, 0.4626, 0.02416, -0.7045);
                case 4: return res(2.996, 0.5242, 0.02689, -0.7508);
                case 5: return res(3.623, 0.5591, 0.0288, -0.7747);
                case 6: return res(4.167, 0.6365, 0.03045, -0.8179);
                case 7: return res(4.795, 0.6819, 0.03145, -0.8408);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 2:
            switch (rhoBin)
            {
                case 1: return res(1.003, 0.4142, 0.02486, -0.6698);
                case 2: return res(2.134, 0.3971, 0.02264, -0.6469);
                case 3: return res(2.66, 0.4566, 0.02755, -0.7058);
                case 4: return res(3.264, 0.4799, 0.02702, -0.7156);
                case 5: return res(3.877, 0.5249, 0.02923, -0.7479);
                case 6: return res(4.441, 0.581, 0.03045, -0.7804);
                case 7: return res(4.742, 0.8003, 0.03613, -0.9062);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 3:
            switch (rhoBin)
            {
                case 1: return res(1.423, 0.4736, 0.03233, -0.7093);
                case 2: return res(2.249, 0.5041, 0.03355, -0.7316);
                case 3: return res(2.961, 0.4889, 0.03129, -0.7091);
                case 4: return res(3.4, 0.5757, 0.03541, -0.7742);
                case 5: return res(3.884, 0.6457, 0.03731, -0.8146);
                case 6: return res(4.433, 0.7524, 0.03962, -0.8672);
                case 7: return res(4.681, 0.9075, 0.04182, -0.9304);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 4:
            switch (rhoBin)
            {
                case 1: return res(-0.7275, 0.8099, 0.04885, -0.9097);
                case 2: return res(1.829, 0.8156, 0.04991, -0.9145);
                case 3: return res(2.72, 0.8454, 0.05036, -0.9215);
                case 4: return res(3.07, 0.9201, 0.05067, -0.9439);
                case 5: return res(3.991, 0.8715, 0.05041, -0.9151);
                case 6: return res(4.001, 1.14, 0.05214, -0.9987);
                case 7: return res(4.522, 1.22, 0.05122, -1);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 5:
            switch (rhoBin)
            {
                case 1: return res(-1.692, 1.192, 0.05049, -1.06);
                case 2: return res(-1.804, 1.48, 0.05315, -1.145);
                case 3: return res(1.673, 1.402, 0.0536, -1.116);
                case 4: return res(2.906, 1.305, 0.05377, -1.076);
                case 5: return res(2.766, 1.613, 0.05511, -1.137);
                case 6: return res(3.409, 1.746, 0.05585, -1.143);
                case 7: return res(3.086, 2.034, 0.05795, -1.181);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 6:
            switch (rhoBin)
            {
                case 1: return res(-0.8823, 1.092, 0.03599, -1.062);
                case 2: return res(2.193, 0.9891, 0.03382, -1.012);
                case 3: return res(2.9, 1.043, 0.03477, -1.019);
                case 4: return res(2.371, 1.488, -0.04053, -1.145);
                case 5: return res(3.75, 1.458, 0.04346, -1.122);
                case 6: return res(3.722, 1.808, 0.04668, -1.177);
                case 7: return res(4.836, 1.47, 0.03875, -1.047);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 7:
            switch (rhoBin)
            {
                case 1: return res(1.184, 0.8944, 0.03233, -1.005);
                case 2: return res(1.691, 1.124, 0.03736, -1.094);
                case 3: return res(2.837, 1.077, 0.03437, -1.046);
                case 4: return res(2.459, 1.589, -0.04007, -1.18);
                case 5: return res(4.058, 1.369, -0.03922, -1.087);
                case 6: return res(4.231, 1.679, 0.0432, -1.13);
                case 7: return res(2.635, 2.648, 0.04929, -1.28);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 8:
            switch (rhoBin)
            {
                case 1: return res(0.3022, 1.127, 0.03826, -1.134);
                case 2: return res(2.161, 1.217, 0.03826, -1.142);
                case 3: return res(3.218, 1.21, 0.03662, -1.112);
                case 4: return res(3.328, 1.638, 0.04398, -1.216);
                case 5: return res(5.506, 1.173, 0.04403, -1.054);
                case 6: return res(-2.444, 3.613, 0.05639, -1.437);
                case 7: return res(2.217, 3.133, 0.05032, -1.338);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 9:
            switch (rhoBin)
            {
                case 1: return res(3.125, 0.6026, 0.02576, -0.8702);
                case 2: return res(3.935, 0.6533, 0.02587, -0.889);
                case 3: return res(4.198, 1.024, 0.03618, -1.069);
                case 4: return res(2.948, 2.386, 0.04771, -1.382);
                case 5: return res(4.415, 2.086, 0.04704, -1.294);
                case 6: return res(-3.084, 4.156, 0.05366, -1.503);
                case 7: return res(-6.144, 5.969, 0.05633, -1.602);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 10:
            switch (rhoBin)
            {
                case 1: return res(4.244, 0.2766, -1.86e-08, -0.5068);
                case 2: return res(4.919, 0.3193, 5.463e-06, -0.58);
                case 3: return res(5.909, 0.2752, 4.144e-06, -0.5272);
                case 4: return res(-47.31, 47.18, 0.05853, -1.991);
                case 5: return res(-46.49, 46.33, 0.05698, -1.989);
                case 6: return res(8.651, 0.2522, 6.592e-06, -0.4835);
                case 7: return res(7.716, 2.481, 0.0531, -1.455);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 11:
            switch (rhoBin)
            {
                case 1: return res(4.467, 0.1997, -3.491e-06, -0.2623);
                case 2: return res(4.17, 0.928, 0.07702, -1.063);
                case 3: return res(-0.04491, 3.67, 0.08704, -1.641);
                case 4: return res(5.528, 1.286, 0.07962, -1.187);
                case 5: return res(-78.36, 78.23, 0.08448, -1.996);
                case 6: return res(7.559, 1.147, 0.07023, -1.134);
                case 7: return res(-59.03, 59.03, -0.08184, -1.992);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 12:
            switch (rhoBin)
            {
                case 1: return res(0.0002851, 3.01, 0.1382, -1.702);
                case 2: return res(-33.01, 33.04, 0.1343, -1.991);
                case 3: return res(-67.94, 67.8, 0.1342, -1.996);
                case 4: return res(-47.81, 48, 0.1391, -1.996);
                case 5: return res(7.162, 0.9211, 0.1395, -1.209);
                case 6: return res(8.193, 0.1995, 2.822e-05, -0.132);
                case 7: return res(8.133, 0.9983, 0.1349, -1.181);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 13:
            switch (rhoBin)
            {
                case 1: return res(2.511, 0.3167, 0.09085, -0.7407);
                case 2: return res(3.297, 0.2091, 6.258e-05, -0.2755);
                case 3: return res(1.85, 2.281, 0.1042, -1.635);
                case 4: return res(3.869, 1.001, 0.09955, -1.266);
                case 5: return res(-23.98, 24.11, 0.1057, -1.988);
                case 6: return res(5.403, 0.2371, 1.5e-05, -0.3177);
                case 7: return res(5.753, 0.2337, 0.0002982, -0.3108);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        default:
            throw std::runtime_error("Eta " + std::to_string(eta)
                                     + " out of range");
    }
}

double
    Cuts::jet2017PtSimRes(const double pt, const double eta, const double rho)
{
    if (pt < 15 || pt > 3000)
    {
        throw std::runtime_error("pT " + std::to_string(pt)
                                 + " out of range to assign resolution");
    }

    static constexpr std::array<double, 14> etaBinEdges{
        0, 0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3, 3.2, 4.7};
    static constexpr std::array<double, 8> rhoBinEdges{
        0, 6.37, 12.4, 18.42, 24.45, 30.47, 36.49, 42.52};
    const auto res = [pt](const double p0,
                          const double p1,
                          const double p2,
                          const double p3) {
        return (std::sqrt(p0 * std::abs(p0) / (pt * pt)
                          + p1 * p1 * std::pow(pt, p3) + p2 * p2));
    };
    const auto etaBin{std::distance(etaBinEdges.begin(),
                                    std::upper_bound(etaBinEdges.begin(),
                                                     etaBinEdges.end(),
                                                     std::abs(eta)))};
    const auto rhoBin{std::distance(
        rhoBinEdges.begin(),
        std::upper_bound(rhoBinEdges.begin(), rhoBinEdges.end(), rho))};

    // https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK4PFchs.txt
    switch (etaBin)
    {
        case 1:
            switch (rhoBin)
            {
                case 1: return res(-1.515, 0.5971, 0.03046, -0.7901);
                case 2: return res(-0.7966, 0.6589, 0.03119, -0.8237);
                case 3: return res(1.387, 0.6885, 0.03145, -0.8378);
                case 4: return res(2.151, 0.7185, 0.03168, -0.8502);
                case 5: return res(2.73, 0.7361, 0.03184, -0.8548);
                case 6: return res(3.603, 0.7318, 0.03227, -0.855);
                case 7: return res(3.897, 0.7882, 0.03282, -0.8746);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 2:
            switch (rhoBin)
            {
                case 1: return res(-0.9395, 0.4556, 0.02738, -0.6909);
                case 2: return res(1.339, 0.4621, 0.02785, -0.6965);
                case 3: return res(1.597, 0.5254, 0.02952, -0.7407);
                case 4: return res(2.527, 0.5042, 0.02842, -0.723);
                case 5: return res(2.896, 0.5428, 0.03001, -0.7476);
                case 6: return res(3.514, 0.5437, 0.03055, -0.7486);
                case 7: return res(3.678, 0.6372, 0.03325, -0.8053);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 3:
            switch (rhoBin)
            {
                case 1: return res(-0.8118, 0.491, 0.03583, -0.7149);
                case 2: return res(1.289, 0.49, 0.03539, -0.7073);
                case 3: return res(1.953, 0.5161, 0.03658, -0.7295);
                case 4: return res(2.347, 0.5396, 0.03576, -0.7339);
                case 5: return res(2.794, 0.5687, 0.03825, -0.7602);
                case 6: return res(2.796, 0.7203, 0.04074, -0.8431);
                case 7: return res(3.788, 0.6287, 0.04156, -0.7959);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 4:
            switch (rhoBin)
            {
                case 1: return res(0.6707, 0.5839, 0.04697, -0.752);
                case 2: return res(1.395, 0.6702, 0.0496, -0.8152);
                case 3: return res(2.43, 0.5712, 0.04572, -0.7345);
                case 4: return res(2.439, 0.6623, 0.04496, -0.7771);
                case 5: return res(3.353, 0.5924, 0.04617, -0.7384);
                case 6: return res(3.465, 0.7579, 0.05328, -0.8435);
                case 7: return res(1.982, 1.148, 0.05664, -0.9626);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 5:
            switch (rhoBin)
            {
                case 1: return res(-1.469, 0.9562, 0.05101, -0.955);
                case 2: return res(-1.377, 1.078, 0.05427, -1.003);
                case 3: return res(1.501, 1.072, 0.05498, -1.001);
                case 4: return res(1.53, 1.158, 0.05396, -1.021);
                case 5: return res(1.621, 1.358, 0.0578, -1.078);
                case 6: return res(3.163, 1.131, 0.05725, -0.9809);
                case 7: return res(2.818, 1.326, 0.05893, -0.9977);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 6:
            switch (rhoBin)
            {
                case 1: return res(1.227, 0.8407, -0.0232, -0.9284);
                case 2: return res(-1.339, 1.218, -0.03479, -1.076);
                case 3: return res(-2.011, 1.435, -0.03565, -1.124);
                case 4: return res(3.324, 0.8102, -0.02662, -0.8923);
                case 5: return res(2.188, 1.365, -0.0375, -1.088);
                case 6: return res(2.884, 1.306, 0.03685, -1.038);
                case 7: return res(4.03, 1.141, 0.03059, -0.9262);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 7:
            switch (rhoBin)
            {
                case 1: return res(-1.979, 1.193, -0.03497, -1.109);
                case 2: return res(-2.528, 1.44, -0.03273, -1.143);
                case 3: return res(1.95, 1.118, -0.03202, -1.054);
                case 4: return res(2.377, 1.166, -0.03593, -1.061);
                case 5: return res(3.122, 1.107, -0.0292, -1.005);
                case 6: return res(-1.899, 1.944, 0.03736, -1.185);
                case 7: return res(4.168, 1.452, 0.03836, -1.019);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 8:
            switch (rhoBin)
            {
                case 1: return res(1.947, 0.9639, -0.02799, -1.024);
                case 2: return res(2.643, 0.9054, -0.02701, -0.9753);
                case 3: return res(-3.209, 2.521, -0.04442, -1.385);
                case 4: return res(-5.368, 3.81, -0.04587, -1.525);
                case 5: return res(-2.344, 2.207, 0.03446, -1.265);
                case 6: return res(-11.01, 8.354, 0.05639, -1.706);
                case 7: return res(6.282, 1.064, 8.482e-06, -0.8687);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 9:
            switch (rhoBin)
            {
                case 1: return res(3.639, 0.6502, -0.01427, -0.8624);
                case 2: return res(2.391, 1.635, -0.0378, -1.251);
                case 3: return res(3.431, 1.985, 0.04609, -1.359);
                case 4: return res(5.095, 0.8757, -0.02736, -0.9761);
                case 5: return res(5.034, 1.479, -0.03479, -1.175);
                case 6: return res(6.694, 1.325, 0.03374, -1.101);
                case 7: return res(7.444, 1.137, 4.258e-05, -0.9531);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 10:
            switch (rhoBin)
            {
                case 1: return res(6.114, 0.2385, 1.741e-05, -0.5054);
                case 2: return res(6.931, 0.1964, 7.465e-06, -0.4335);
                case 3: return res(7.858, 0.2435, 6.026e-07, -0.5235);
                case 4: return res(8.713, 0.1314, 8.441e-06, -0.3028);
                case 5: return res(9.413, 0.2792, 1.217e-06, -0.5729);
                case 6: return res(10.51, 0.1659, 1.277e-06, -0.4276);
                case 7: return res(11.77, 8.547e-07, 0.05169, -1.197);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 11:
            switch (rhoBin)
            {
                case 1: return res(6.048, 0.1992, -3.559e-06, -0.2953);
                case 2: return res(6.867, 0.2036, 1.946e-05, -0.3068);
                case 3: return res(8.198, 0.0001314, 0.08772, -1.252);
                case 4: return res(8.756, 0.134, -0.07197, -0.2968);
                case 5: return res(9.615, 0.0001533, -0.08793, -1.445);
                case 6: return res(10.01, 0.1524, 3.815e-05, -0.2422);
                case 7: return res(10.05, 0.1932, 0.0001734, -0.2739);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 12:
            switch (rhoBin)
            {
                case 1: return res(-35.12, 35.21, 0.1466, -1.993);
                case 2: return res(6.573, 0.2026, 6.573e-05, -0.1564);
                case 3: return res(0.004144, 6.019, 0.1549, -1.854);
                case 4: return res(8.341, 0.0001012, 0.1526, -1.689);
                case 5: return res(9.115, 0.0002242, 0.1518, -1.362);
                case 6: return res(9.86, -2.112e-05, 0.1438, -1.114);
                case 7: return res(10.45, 0.0001536, 0.1398, -1.271);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        case 13:
            switch (rhoBin)
            {
                case 1: return res(-29.87, 29.84, 0.1045, -1.995);
                case 2: return res(-23.2, 23.09, 0.1051, -1.987);
                case 3: return res(4.337, 0.2253, 0.06986, -0.4215);
                case 4: return res(4.088, 2.746, 0.1136, -1.959);
                case 5: return res(5.624, 0.1291, 0.002663, -0.04825);
                case 6: return res(6.152, 6.125e-05, 0.1128, -1.319);
                case 7: return res(6.235, 0.1408, 0.0001266, -0.08163);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + " out of range");
            }
        default:
            throw std::runtime_error("Eta " + std::to_string(eta)
                                     + " out of range");
    }
}

std::pair<double, double> Cuts::jet2016SFs(const double eta)
{
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Uncertainty
    constexpr std::array<double, 14> etaBinEdges{0,
                                                 0.522,
                                                 0.783,
                                                 1.131,
                                                 1.305,
                                                 1.740,
                                                 1.930,
                                                 2.043,
                                                 2.322,
                                                 2.5,
                                                 2.853,
                                                 2.964,
                                                 3.319,
                                                 5.191};

    switch (std::distance(
        etaBinEdges.begin(),
        std::upper_bound(etaBinEdges.begin(), etaBinEdges.end(), eta)))
    {
        case 1: return {1.1685, 0.0645};
        case 2: return {1.1948, 0.0652};
        case 3: return {1.1464, 0.0632};
        case 4: return {1.1609, 0.1025};
        case 5: return {1.1278, 0.0986};
        case 6: return {1.1000, 0.1079};
        case 7: return {1.1426, 0.1214};
        case 8: return {1.1512, 0.1440};
        case 9: return {1.2963, 0.2371};
        case 10: return {1.3418, 0.2091};
        case 11: return {1.7788, 0.2008};
        case 12: return {1.1869, 0.1243};
        case 13: return {1.1922, 0.1448};
        default:
            throw std::runtime_error("Eta " + std::to_string(eta)
                                     + " out of range");
    }
}

std::pair<double, double> Cuts::jet2017SFs(const double eta)
{
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Uncertainty
    constexpr std::array<double, 14> etaBinEdges{0,
                                                 0.522,
                                                 0.783,
                                                 1.131,
                                                 1.305,
                                                 1.740,
                                                 1.930,
                                                 2.043,
                                                 2.322,
                                                 2.5,
                                                 2.853,
                                                 2.964,
                                                 3.319,
                                                 5.191};

    switch (std::distance(
        etaBinEdges.begin(),
        std::upper_bound(etaBinEdges.begin(), etaBinEdges.end(), eta)))
    {
        case 1: return {1.1432, 0.0222};
        case 2: return {1.1815, 0.0484};
        case 3: return {1.0989, 0.0456};
        case 4: return {1.1137, 0.1397};
        case 5: return {1.1307, 0.147};
        case 6: return {1.16, 0.0976};
        case 7: return {1.2393, 0.1909};
        case 8: return {1.2604, 0.1501};
        case 9: return {1.4085, 0.2020};
        case 10: return {1.9909, 0.5684};
        case 11: return {2.2923, 0.3743};
        case 12: return {1.2696, 0.1089};
        case 13: return {1.1542, 0.1524};
        default:
            throw std::runtime_error("Eta " + std::to_string(eta)
                                     + " out of range");
    }
}

void Cuts::getBWeight(const AnalysisEvent& event,
                      const TLorentzVector jet,
                      const int index,
                      double& mcTag,
                      double& mcNoTag,
                      double& dataTag,
                      double& dataNoTag,
                      double& err1,
                      double& err2,
                      double& err3,
                      double& err4) const
{
    // Use b-tagging efficiencies and scale factors.
    // Firstly get efficiency for pt/eta bin here.
    double eff{1.};
    int partonFlavour{std::abs(event.jetPF2PATPID[index])};

    if (partonFlavour == 0)
    {
        return;
    }
    if (partonFlavour == 5)
    {
        eff = bTagEffPlots_[4]->GetBinContent(
                  bTagEffPlots_[4]->GetXaxis()->FindBin(jet.Pt()),
                  bTagEffPlots_[4]->GetYaxis()->FindBin(std::abs(jet.Eta())))
              / bTagEffPlots_[0]->GetBinContent(
                    bTagEffPlots_[0]->GetXaxis()->FindBin(jet.Pt()),
                    bTagEffPlots_[0]->GetYaxis()->FindBin(std::abs(jet.Eta())));
    }
    if (partonFlavour == 4)
    {
        eff = bTagEffPlots_[5]->GetBinContent(
                  bTagEffPlots_[5]->GetXaxis()->FindBin(jet.Pt()),
                  bTagEffPlots_[5]->GetYaxis()->FindBin(std::abs(jet.Eta())))
              / bTagEffPlots_[1]->GetBinContent(
                    bTagEffPlots_[1]->GetXaxis()->FindBin(jet.Pt()),
                    bTagEffPlots_[1]->GetYaxis()->FindBin(std::abs(jet.Eta())));
    }
    if (partonFlavour < 4)
    {
        eff = bTagEffPlots_[6]->GetBinContent(
                  bTagEffPlots_[6]->GetXaxis()->FindBin(jet.Pt()),
                  bTagEffPlots_[6]->GetYaxis()->FindBin(std::abs(jet.Eta())))
              / bTagEffPlots_[2]->GetBinContent(
                    bTagEffPlots_[2]->GetXaxis()->FindBin(jet.Pt()),
                    bTagEffPlots_[2]->GetYaxis()->FindBin(std::abs(jet.Eta())));
    }
    if (partonFlavour == 21)
    {
        eff = bTagEffPlots_[7]->GetBinContent(
                  bTagEffPlots_[7]->GetXaxis()->FindBin(jet.Pt()),
                  bTagEffPlots_[7]->GetYaxis()->FindBin(std::abs(jet.Eta())))
              / bTagEffPlots_[3]->GetBinContent(
                    bTagEffPlots_[3]->GetXaxis()->FindBin(jet.Pt()),
                    bTagEffPlots_[3]->GetYaxis()->FindBin(std::abs(jet.Eta())));
    }

    if (std::isnan(eff))
    {
        std::cerr << "WARN: NaN encountered calculating bTag efficiency, "
                     "setting to 1.0. Check efficiency plots."
                  << std::endl;
        eff = 1;
    }

    // Get SF
    // Initalise variables.
    double jet_scalefactor{1.};
    double jet_scalefactor_up{1.};
    double jet_scalefactor_do{1.};

    double SFerr{0.};
    double jetPt{jet.Pt()};
    constexpr double maxBjetPt{670};
    constexpr double maxLjetPt{1000.0};
    bool doubleUncertainty{false};
    // Do some things if it's a b or c

    if (partonFlavour == 5)
    {
        if (jetPt > maxBjetPt)
        {
            jetPt = maxBjetPt;
            doubleUncertainty = true;
        }
        jet_scalefactor = getBSF(0, 0, jetPt);
        jet_scalefactor_up = getBSF(0, 1, jetPt);
        jet_scalefactor_do = getBSF(0, -1, jetPt);
    }

    else if (partonFlavour == 4)
    {
        if (jetPt > maxBjetPt)
        {
            jetPt = maxBjetPt;
            doubleUncertainty = true;
        }
        jet_scalefactor = getBSF(1, 0, jetPt);
        jet_scalefactor_up = getBSF(1, 1, jetPt);
        jet_scalefactor_do = getBSF(1, -1, jetPt);
    }

    // Light jets
    else
    {
        if (jetPt > maxLjetPt)
        {
            jetPt = maxLjetPt;
            doubleUncertainty = true;
        }
        jet_scalefactor = getBSF(2, 0, jetPt);
        jet_scalefactor_up = getBSF(2, 1, jetPt);
        jet_scalefactor_do = getBSF(2, -1, jetPt);
    }

    if (doubleUncertainty)
    {
        jet_scalefactor_up =
            2 * (jet_scalefactor_up - jet_scalefactor) + jet_scalefactor;
        jet_scalefactor_do =
            2 * (jet_scalefactor_do - jet_scalefactor) + jet_scalefactor;
    }

    SFerr = std::abs(jet_scalefactor_up - jet_scalefactor)
                    > std::abs(jet_scalefactor_do - jet_scalefactor)
                ? std::abs(jet_scalefactor_up - jet_scalefactor)
                : std::abs(jet_scalefactor_do - jet_scalefactor);

    // Apply the weight of the jet and set the error
    if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[index]
        > bDiscCut_)
    {
        mcTag *= eff;
        dataTag *= eff * jet_scalefactor;

        if (partonFlavour == 5 || partonFlavour == 4)
        {
            err1 += SFerr / jet_scalefactor;
        }
        else
        {
            err3 += SFerr / jet_scalefactor;
        }
    }
    else
    {
        mcNoTag *= (1 - eff);
        dataNoTag *= (1 - eff * jet_scalefactor);

        if (partonFlavour == 5 || partonFlavour == 4)
        {
            err2 += (-eff * SFerr) / (1 - eff * jet_scalefactor);
        }
        else
        {
            err4 += (-eff * SFerr) / (1 - eff * jet_scalefactor);
        }
    }
}

// Backup temporary method to do Btag Scale Factors whilst debugging is
// ongoing.
// TODO: F1X TH1S

double Cuts::getBSF(const int flavour, const int type, const double pt) const
{
    if (!is2016_)
    { // is 2017
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        // MEDIUM

        static constexpr std::array<double, 10> ptBinEdges{
            20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};
        const auto ptBin{std::distance(
            ptBinEdges.begin(),
            std::upper_bound(ptBinEdges.begin(), ptBinEdges.end(), pt))};
        const auto mujets_sf = [pt, type](const double p0) {
            return (0.941966 * ((1 + 0.24108 * pt) / (1 + 0.248776 * pt)))
                   + type * p0;
        };
        const auto incl_sf = [pt, type]() {
            return (0.949449 + 0.000516201 * pt + 7.13398e-08 * pt * pt
                    + -3.55644e-10 * pt * pt * pt)
                   * (1 + type * 0.082197);
        };

        switch (flavour)
        {
            case 0: // B flavour
                if (type == 0)
                {
                    return mujets_sf(0);
                }
                else if (std::abs(type) == 1)
                {
                    switch (ptBin)
                    {
                        case 1: return mujets_sf(0.051529459655284882);
                        case 2: return mujets_sf(0.017671864479780197);
                        case 3: return mujets_sf(0.022306634113192558);
                        case 4: return mujets_sf(0.023042259737849236);
                        case 5: return mujets_sf(0.039661582559347153);
                        case 6: return mujets_sf(0.061514820903539658);
                        case 7: return mujets_sf(0.071018315851688385);
                        case 8: return mujets_sf(0.054169680923223495);
                        case 9: return mujets_sf(0.063008971512317657);
                        default:
                            throw std::runtime_error(
                                "pT out of range of b tag SFs");
                    }
                }
                else
                {
                    throw std::runtime_error("Unknown b tag systematic type");
                }
            case 1: // C flavour
                if (type == 0)
                {
                    return mujets_sf(0);
                }
                else if (std::abs(type) == 1)
                {
                    switch (ptBin)
                    {
                        case 1: return mujets_sf(0.15458837151527405);
                        case 2: return mujets_sf(0.053015593439340591);
                        case 3: return mujets_sf(0.066919900476932526);
                        case 4: return mujets_sf(0.069126777350902557);
                        case 5: return mujets_sf(0.11898474395275116);
                        case 6: return mujets_sf(0.18454445898532867);
                        case 7: return mujets_sf(0.21305495500564575);
                        case 8: return mujets_sf(0.16250903904438019);
                        case 9: return mujets_sf(0.18902692198753357);
                        default:
                            throw std::runtime_error(
                                "pT out of range of b tag SFs");
                    }
                }
                else
                {
                    throw std::runtime_error("Unknown b tag systematic type");
                }
            case 2: // UDSG flavour
                if (std::abs(type) <= 1)
                {
                    return incl_sf();
                }
                else
                {
                    throw std::runtime_error("Unknown b tag systematic type");
                }
            default:
                throw std::runtime_error("Unknown b tag systematic flavour");
        }
    }
    else
    { // is 2016

        double sf{1.0};
        const double& x{pt};

        // MEDIUM
        if (flavour == 0)
        { // B flavour
            if (type == 0)
                sf = 0.718014
                     * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x)));
            if (type == 1)
            {
                if (pt < 30.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.040554910898208618;
                if (pt < 50.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.01836167648434639;
                if (pt < 70.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.016199169680476189;
                if (pt < 100.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.014634267427027225;
                if (pt < 140.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.014198922552168369;
                if (pt < 200.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.016547618433833122;
                if (pt < 300.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.02140621654689312;
                if (pt < 600.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.023563217371702194;
                else
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.034716218709945679;
            }
            if (type == -1)
            {
                if (pt < 30.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.040554910898208618;
                if (pt < 50.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.01836167648434639;
                if (pt < 70.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.016199169680476189;
                if (pt < 100.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.014634267427027225;
                if (pt < 140.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.014198922552168369;
                if (pt < 200.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.016547618433833122;
                if (pt < 300.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.02140621654689312;
                if (pt < 600.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.023563217371702194;
                else
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.034716218709945679;
            }
        }
        if (flavour == 1)
        { // C flavour
            if (type == 0)
                sf = 0.718014
                     * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x)));
            if (type == 1)
            {
                if (pt < 30.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.12166473269462585;
                if (pt < 50.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.055085029453039169;
                if (pt < 70.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.048597507178783417;
                if (pt < 100.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.043902803212404251;
                if (pt < 140.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.042596768587827682;
                if (pt < 200.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.049642853438854218;
                if (pt < 300.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.06421864777803421;
                if (pt < 600.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.070689648389816284;
                else
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         + 0.10414865612983704;
            }
            if (type == -1)
            {
                if (pt < 30.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.12166473269462585;
                if (pt < 50.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.055085029453039169;
                if (pt < 70.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.048597507178783417;
                if (pt < 100.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.043902803212404251;
                if (pt < 140.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.042596768587827682;
                if (pt < 200.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.049642853438854218;
                if (pt < 300.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.06421864777803421;
                if (pt < 600.0)
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.070689648389816284;
                else
                    sf = (0.718014
                          * ((1. + (0.0685826 * x)) / (1. + (0.0475779 * x))))
                         - 0.10414865612983704;
            }
        }
        if (flavour == 2)
        { // UDSG flavour
            if (type == 0)
                sf = 1.0589 + 0.000382569 * x + -2.4252e-07 * x * x
                     + 2.20966e-10 * x * x * x;
            if (type == 1)
                sf =
                    (1.0589 + 0.000382569 * x + -2.4252e-07 * x * x
                     + 2.20966e-10 * x * x * x)
                    * (1 + (0.100485 + 3.95509e-05 * x + -4.90326e-08 * x * x));
            if (type == -1)
                sf =
                    (1.0589 + 0.000382569 * x + -2.4252e-07 * x * x
                     + 2.20966e-10 * x * x * x)
                    * (1 - (0.100485 + 3.95509e-05 * x + -4.90326e-08 * x * x));
        }

        return sf;
    }
}
