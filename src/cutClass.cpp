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

        // Tight ID
        h_muonIDs1 = dynamic_cast<TH2D*>(
            muonIDsFile1->Get("NUM_TightID_DEN_genTracks_pt_abseta"));
        h_muonPFiso1 = dynamic_cast<TH2D*>(
            muonIsoFile1->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"));
        std::cout << "Got 2017 muon SFs!\n" << std::endl;
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
        muonIDsFile1->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta");
        muonIDsFile2->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta");
        dynamic_cast<TH2F*>(
            muonIDsFile1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/"
                              "abseta_pt_ratio"))
            ->Copy(*h_muonIDs1);
        h_muonIDs2 = dynamic_cast<TH2F*>(
            muonIDsFile2->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/"
                              "abseta_pt_ratio"));

        // Tight Iso
        muonIsoFile1->cd("TightISO_TightID_pt_eta");
        muonIsoFile2->cd("TightISO_TightID_pt_eta");
        dynamic_cast<TH2F*>(
            muonIsoFile1->Get("TightISO_TightID_pt_eta/abseta_pt_ratio"))
            ->Copy(*h_muonPFiso1);
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
        makeJetCuts(event, systToRun, eventWeight, false);

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
                // the same properties
                size_t seed{0};
                boost::hash_combine(seed, event.muonPF2PATCharge[*muonIt]);
                boost::hash_combine(seed, event.muonPF2PATPt[*muonIt]);
                boost::hash_combine(seed, event.muonPF2PATEta[*muonIt]);
                boost::hash_combine(seed, event.muonPF2PATPhi[*muonIt]);
                boost::hash_combine(
                    seed, event.muonPF2PATTkLysWithMeasurements[*muonIt]);

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
                eeWeight = 0.926;
                mumuWeight = 1.114;
                emuWeight = 1.489;
            }
            else
            {
                eeWeight = 1.188;
                mumuWeight = 1.131;
                emuWeight = 1.197;
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
        const TLorentzVector tempVec{event.elePF2PATPX[i],
                                     event.elePF2PATPY[i],
                                     event.elePF2PATPZ[i],
                                     event.elePF2PATE[i]};

        if (electrons.size() < 1 && tempVec.Pt() <= tightElePtLeading_)
            continue;
        else if (electrons.size() >= 1 && tempVec.Pt() <= tightElePt_)
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
        const TLorentzVector tempVec{event.elePF2PATPX[i],
                                     event.elePF2PATPY[i],
                                     event.elePF2PATPZ[i],
                                     event.elePF2PATE[i]};

        if (electrons.size() < 1 && tempVec.Pt() <= looseElePtLeading_)
            continue;
        else if (electrons.size() >= 1 && tempVec.Pt() <= looseElePt_)
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
    return muons;
}

std::vector<int> Cuts::getLooseMuons(const AnalysisEvent& event) const
{
    std::vector<int> muons;
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

        event.zPairLeptons.first =
            lepton1.Pt() > lepton2.Pt() ? lepton1 : lepton2;
        event.zPairIndex.first =
            lepton1.Pt() > lepton2.Pt() ? electrons[0] : electrons[1];
        event.zPairRelIso.first =
            lepton1.Pt() > lepton2.Pt()
                ? event.elePF2PATComRelIsoRho[electrons[0]]
                : event.elePF2PATComRelIsoRho[electrons[1]];
        event.zPairLeptons.second =
            lepton1.Pt() > lepton2.Pt() ? lepton2 : lepton1;
        event.zPairRelIso.second =
            lepton1.Pt() > lepton2.Pt()
                ? event.elePF2PATComRelIsoRho[electrons[1]]
                : event.elePF2PATComRelIsoRho[electrons[0]];
        event.zPairIndex.second =
            lepton1.Pt() > lepton2.Pt() ? electrons[1] : electrons[0];

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
                // Jet ID == tight (loose is deprecated)
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
        if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[jets[i]]
            <= bDiscCut_)
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
                twgt = 0.93124;
                if (syst == 1)
                {
                    twgt += 0.00239;
                }
                else if (syst == 2)
                {
                    twgt -= 0.00239;
                }
            }
        }
        else if (channel == "mumu")
        {
            if (muTrig || mumuTrig)
            {
                twgt = 0.96698;
                if (syst == 1)
                {
                    twgt += 0.00022;
                }
                else if (syst == 2)
                {
                    twgt -= 0.00022;
                }
            }
        }
        else if (channel == "emu")
        {
            if (muEGTrig)
            {
                twgt = 0.94875;
                if (syst == 1)
                {
                    twgt += 0.00425;
                }
                else if (syst == 2)
                {
                    twgt -= 0.00425;
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
                twgt = 0.98715; // 0.97554 for data eff; 0.98715 for SF
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

                twgt = (0.98868 * lumiRunsBCDEF_ + 0.99868 * lumiRunsGH_)
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
                twgt = 0.87661; // 0.87661 for eff; 0.99399 for SF
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
            || event.Flag_muonBadTrackFilter <= 0 || event.Flag_badMuons <= 0
            || event.Flag_duplicateMuons <= 0 || event.Flag_noBadMuons))
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
            if (!is2016_)
            {
                muonIdSF += h_muonIDs1->GetBinError(binId1);
                muonPFisoSF += h_muonPFiso1->GetBinError(binIso1);
            }
            else
            {
                muonIdSF += (h_muonIDs1->GetBinError(binId1) * lumiRunsBCDEF_
                             + h_muonIDs2->GetBinError(binId2) * lumiRunsGH_)
                                / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06)
                            + 0.01; // Additional 1% uncert for ID and 0.5% for
                                    // iso as recommended
                muonPFisoSF +=
                    (h_muonPFiso1->GetBinError(binIso1) * lumiRunsBCDEF_
                     + h_muonIDs2->GetBinError(binId2) * lumiRunsGH_)
                        / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06)
                    + 0.005;
            }
        }
        else if (syst == 2)
        {
            if (!is2016_)
            {
                muonIdSF -= h_muonIDs1->GetBinError(binId1);
                muonPFisoSF -= h_muonPFiso1->GetBinError(binIso1);
            }
            else
            {
                muonIdSF -= (h_muonIDs1->GetBinError(binId1) * lumiRunsBCDEF_
                             + h_muonIDs2->GetBinError(binId2) * lumiRunsGH_)
                                / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06)
                            - 0.01; // Additional 1% uncert for ID and 0.5% for
                                    // iso as recommended
                muonPFisoSF -=
                    (h_muonPFiso1->GetBinError(binIso1) * lumiRunsBCDEF_
                     + h_muonIDs2->GetBinError(binId2) * lumiRunsGH_)
                        / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06)
                    - 0.005;
            }
        }

        return muonIdSF * muonPFisoSF;
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

        if (!isMC_)
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
    if (event.jetPF2PATPtRaw[index] < 15 || event.jetPF2PATPtRaw[index] > 3000
        || event.elePF2PATRhoIso[0] > 42.52
        || std::abs(event.jetPF2PATEta[index]) > 4.7)
    {
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);
        return {returnJet, newSmearValue};
    }

    // TODO: Should this be gen or reco level?
    // I think reco because gen might not exist? (does not exist when smearing)
    const double ptRes{jet2017PtSimRes(event.jetPF2PATPtRaw[index],
                                       event.jetPF2PATEta[index],
                                       event.elePF2PATRhoIso[0])};

    const auto dR{deltaR(event.genJetPF2PATEta[index],
                         event.genJetPF2PATPhi[index],
                         event.jetPF2PATEta[index],
                         event.jetPF2PATPhi[index])};
    const double dPt{event.jetPF2PATPtRaw[index] - event.genJetPF2PATPT[index]};
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

    // If there is no matching gen jet, gen jet pT is -999;
    if (event.genJetPF2PATPT[index] != -999.0 && dR < (0.4 / 2.0)
        && std::abs(dPt) < 3.0 * ptRes * event.jetPF2PATPtRaw[index])
    // If matching from GEN to RECO using dR<Rcone/2 and dPt < 3*sigma,
    // just scale
    {
        newSmearValue =
            std::max(1. + (jerSF - 1.) * dPt / event.jetPF2PATPtRaw[index], 0.);
    }
    else // If not matched to a gen jet, randomly smear
    {
        std::normal_distribution<> d(
            0, ptRes * std::sqrt(std::max(jerSF * jerSF - 1, 0.)));

        // Like with the Rochester corrections, seed the random number generator
        // with event (jet) properties so that each jet is smeared the same
        // way every time it is processed
        size_t seed{0};
        boost::hash_combine(seed, event.jetPF2PATPtRaw[index]);
        boost::hash_combine(seed, event.jetPF2PATEta[index]);
        boost::hash_combine(seed, event.jetPF2PATPhi[index]);
        std::mt19937 gen(rand());

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

    // https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK4PF.txt
    switch (etaBin)
    {
        case 1:
            switch (rhoBin)
            {
                case 1: return res(0.2219, 0.6528, 0.03069, -0.8177);
                case 2: return res(2.561, 0.6994, 0.03112, -0.8421);
                case 3: return res(3.519, 0.7661, 0.03154, -0.8717);
                case 4: return res(4.393, 0.7918, 0.03165, -0.882);
                case 5: return res(5.158, 0.8365, 0.03186, -0.8984);
                case 6: return res(5.771, 0.8938, 0.03207, -0.9182);
                case 7: return res(6.315, 0.9692, 0.03227, -0.9414);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 2:
            switch (rhoBin)
            {
                case 1: return res(1.52, 0.4603, 0.02756, -0.6952);
                case 2: return res(2.975, 0.484, 0.02811, -0.7113);
                case 3: return res(3.843, 0.5315, 0.02946, -0.7462);
                case 4: return res(4.654, 0.5523, 0.02959, -0.7586);
                case 5: return res(5.32, 0.604, 0.03038, -0.7884);
                case 6: return res(5.932, 0.6271, 0.03077, -0.7997);
                case 7: return res(6.416, 0.6823, 0.03098, -0.8249);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 3:
            switch (rhoBin)
            {
                case 1: return res(1.642, 0.4859, 0.03535, -0.7101);
                case 2: return res(2.902, 0.5124, 0.03556, -0.7264);
                case 3: return res(3.767, 0.5547, 0.03671, -0.7556);
                case 4: return res(4.492, 0.6036, 0.03727, -0.7821);
                case 5: return res(5.201, 0.6034, 0.03674, -0.7777);
                case 6: return res(5.755, 0.6417, 0.03702, -0.7945);
                case 7: return res(6.187, 0.7263, 0.03816, -0.8379);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 4:
            switch (rhoBin)
            {
                case 1: return res(0.9151, 0.7297, 0.04994, -0.847);
                case 2: return res(2.699, 0.7325, 0.04927, -0.8412);
                case 3: return res(3.454, 0.8175, 0.05047, -0.8812);
                case 4: return res(4.07, 0.9098, 0.05066, -0.9135);
                case 5: return res(5.068, 0.8488, 0.05008, -0.8833);
                case 6: return res(5.336, 1.033, 0.05161, -0.9531);
                case 7: return res(6.068, 0.9915, 0.05083, -0.9302);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 5:
            switch (rhoBin)
            {
                case 1: return res(-3.251, 1.738, 0.05665, -1.203);
                case 2: return res(-3.048, 2, 0.05759, -1.253);
                case 3: return res(-2.064, 2.059, 0.05732, -1.254);
                case 4: return res(0.8834, 2.116, 0.05735, -1.256);
                case 5: return res(1.75, 2.379, 0.05828, -1.296);
                case 6: return res(4.106, 1.934, 0.05572, -1.203);
                case 7: return res(3.677, 2.399, 0.0572, -1.275);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 6:
            switch (rhoBin)
            {
                case 1: return res(-2.203, 1.393, -0.03217, -1.131);
                case 2: return res(1.102, 1.409, -0.03321, -1.135);
                case 3: return res(1.774, 1.624, -0.03351, -1.182);
                case 4: return res(2.296, 1.838, -0.03475, -1.223);
                case 5: return res(3.541, 1.826, -0.0341, -1.216);
                case 6: return res(2.67, 2.452, -0.03659, -1.319);
                case 7: return res(-0.3526, 2.993, 0.03736, -1.377);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 7:
            return res(-1.066, 1.835, -0.03267, -1.24);
            switch (rhoBin)
            {
                case 1: return res(-2.528, 1.591, -0.03475, -1.215);
                case 2: return res(-1.362, 1.571, -0.03228, -1.192);
                case 3: return res(-1.066, 1.835, -0.03267, -1.24);
                case 4: return res(-0.9215, 2.166, -0.03427, -1.297);
                case 5: return res(2.119, 2.071, -0.03182, -1.262);
                case 6: return res(1.717, 2.571, 0.03417, -1.341);
                case 7: return res(-3.694, 3.636, 0.03391, -1.447);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 8:
            switch (rhoBin)
            {
                case 1: return res(-4.688, 2.887, -0.0398, -1.453);
                case 2: return res(-5.267, 3.392, -0.03859, -1.484);
                case 3: return res(-7.486, 5.453, -0.04113, -1.636);
                case 4: return res(-7.68, 5.674, -0.03997, -1.633);
                case 5: return res(-8.179, 6.134, 0.03822, -1.642);
                case 6: return res(-12.07, 9.833, 0.03918, -1.763);
                case 7: return res(-21.68, 19.59, 0.04008, -1.893);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 9:
            return res(-114.5, 114.1, 0.04377, -1.997);
            switch (rhoBin)
            {
                case 1: return res(-85.96, 85.58, -0.04539, -1.996);
                case 2: return res(-92, 91.58, -0.04383, -1.996);
                case 3: return res(-114.5, 114.1, 0.04377, -1.997);
                case 4: return res(-127.4, 127, 0.04242, -1.997);
                case 5: return res(-100.1, 99.54, 0.04101, -1.995);
                case 6: return res(-171.2, 170.9, 0.04004, -1.998);
                case 7: return res(-125, 124.5, 0.03718, -1.996);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 10:
            return res(-115.9, 115.5, 0.03881, -1.997);
            switch (rhoBin)
            {
                case 1: return res(-80.86, 80.33, 0.04414, -1.994);
                case 2: return res(-89.39, 88.88, 0.04186, -1.995);
                case 3: return res(-115.9, 115.5, 0.03881, -1.997);
                case 4: return res(-108.6, 108.2, 0.0389, -1.996);
                case 5: return res(-135.7, 135.4, 0.03607, -1.997);
                case 6: return res(-150.6, 150.3, -0.03642, -1.998);
                case 7: return res(-137.8, 137.7, 0.03824, -1.998);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 11:
            return res(-11.14, 8.098, -0.06207, -1.683);
            switch (rhoBin)
            {
                case 1: return res(-3.681, 1.906, -0.0582, -1.186);
                case 2: return res(-6.155, 3.845, 0.06349, -1.456);
                case 3: return res(-11.14, 8.098, -0.06207, -1.683);
                case 4: return res(-42.17, 40.57, -0.06272, -1.964);
                case 5: return res(-73.19, 72.17, -0.06209, -1.987);
                case 6: return res(-92.46, 91.64, -0.05969, -1.991);
                case 7: return res(-38.08, 36.32, 0.06421, -1.949);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 12:
            switch (rhoBin)
            {
                case 1: return res(-34.28, 33.8, 0.1526, -1.984);
                case 2: return res(4.892, 0.3432, 0.0002247, -0.3053);
                case 3: return res(-76.26, 75.99, 0.1548, -1.995);
                case 4: return res(-90.46, 90.02, 0.1496, -1.995);
                case 5: return res(-12.56, 10.87, 0.1495, -1.801);
                case 6: return res(-23.44, 21.8, 0.1444, -1.913);
                case 7: return res(7.064, 1.836, 0.1405, -1.217);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        case 13:
            switch (rhoBin)
            {
                case 1: return res(-4.261, 4.085, 0.1075, -1.819);
                case 2: return res(-15.95, 15.87, 0.1088, -1.978);
                case 3: return res(3.137, 1.361, 0.1042, -1.42);
                case 4: return res(-21.12, 21.25, 0.1075, -1.986);
                case 5: return res(5.014, 0.2092, 5.362e-06, -0.2474);
                case 6: return res(-29.77, 29.95, 0.1124, -1.992);
                case 7: return res(5.756, 0.2217, 3.682e-05, -0.2657);
                default:
                    throw std::runtime_error("Rho " + std::to_string(rho)
                                             + "out of range");
            }
        default:
            throw std::runtime_error("Eta " + std::to_string(eta)
                                     + " out of range");
    }
}

std::pair<double, double> Cuts::jet2016SFs(const double eta)
{
    // JER Scaling Factors and uncertainities for 2016
    double jerSF{0.};
    double jerSigma{0.};

    if (eta <= 0.5)
    {
        jerSF = 1.109;
        jerSigma = 0.008;
    }
    else if (eta <= 0.8)
    {
        jerSF = 1.138;
        jerSigma = 0.013;
    }
    else if (eta <= 1.1)
    {
        jerSF = 1.114;
        jerSigma = 0.013;
    }
    else if (eta <= 1.3)
    {
        jerSF = 1.123;
        jerSigma = 0.024;
    }
    else if (eta <= 1.7)
    {
        jerSF = 1.084;
        jerSigma = 0.011;
    }
    else if (eta <= 1.9)
    {
        jerSF = 1.082;
        jerSigma = 0.035;
    }
    else if (eta <= 2.1)
    {
        jerSF = 1.140;
        jerSigma = 0.047;
    }
    else if (eta <= 2.3)
    {
        jerSF = 1.067;
        jerSigma = 0.053;
    }
    else if (eta <= 2.5)
    {
        jerSF = 1.177;
        jerSigma = 0.041;
    }
    else if (eta <= 2.8)
    {
        jerSF = 1.364;
        jerSigma = 0.039;
    }
    else if (eta <= 3.0)
    {
        jerSF = 1.857;
        jerSigma = 0.071;
    }
    else if (eta <= 3.2)
    {
        jerSF = 1.328;
        jerSigma = 0.022;
    }
    else
    {
        jerSF = 1.160;
        jerSigma = 0.029;
    }

    return std::make_pair(jerSF, jerSigma);
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

// Backup temporary method to do Btag Scale Factors whilst debugging is ongoing.
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
