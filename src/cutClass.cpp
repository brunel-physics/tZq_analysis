#include "cutClass.hpp"

#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <yaml-cpp/yaml.h>

Cuts::Cuts(const bool doPlots,
           const bool fillCutFlows,
           const bool invertLepCut,
           const bool lepCutFlow,
           const bool dumpEventNumber,
           const bool is2016,
           const bool isFCNC,
           const bool isCtag)
    :

    // Do plots?
    doPlots_{doPlots}
    , fillCutFlow_{fillCutFlows}
    ,
    // background estimation. May not be possible
    invertLepCut_{invertLepCut}
    ,
    // Synchronisation cut flow.
    synchCutFlow_{lepCutFlow}
    ,
    // Synchronisation cut flow.
    singleEventInfoDump_{false}
    , makeEventDump_{dumpEventNumber}
    , isFCNC_{isFCNC}
    , isCtag_{isCtag}
    , is2016_{is2016}
    ,

    // Set all default parameters. These will be editable later on, probably.
    numTightEle_{3}
    , tightElePt_{15.}
    , tightElePtLeading_{35.}
    , tightEleEta_{2.5}
    , tightEled0_{0.011811}
    , tightEleMissLayers_{3}
    , // Cut based has barrel =2; endcap: veto=3, others=1
    tightEleCheckPhotonVeto_{true}
    , tightEleMVA0_{0.972153}
    , // Medium cut
    tightEleMVA1_{0.922126}
    , // Medium cut
    tightEleMVA2_{0.610764}
    , // Medium cut
    // tightEleMVA0_{0.988153}, // Tight cut
    // tightEleMVA1_{0.967910}, // Tight cut
    // tightEleMVA2_{0.841729}, // Tight cut
    tightEleRelIso_{0.107587}
    ,
    // Loose electron initialisation
    numLooseEle_{3}
    , looseElePt_{15.}
    , looseElePtLeading_{35.}
    , looseEleEta_{2.5}
    , looseEleMVA0_{0.972153}
    , looseEleMVA1_{0.922126}
    , looseEleMVA2_{0.610764}
    , looseEleRelIso_{0.15}
    ,
    // Tight muon initialisation
    numTightMu_{0}
    , tightMuonPt_{12.}
    , tightMuonPtLeading_{27.}
    , tightMuonEta_{2.4}
    , tightMuonRelIso_{0.15}
    ,
    // Loose muons
    numLooseMu_{0}
    , looseMuonPt_{12.}
    , looseMuonPtLeading_{27.}
    , looseMuonEta_{2.4}
    , looseMuonRelIso_{0.25}
    ,
    // zMass cuts
    invZMassCut_{20.}
    , invWMassCut_{20.}
    ,
    // Jet initialisation
    numJets_{2}
    , maxJets_{4}
    , jetPt_{30.}
    , jetEta_{5.0}
    , jetNConsts_{2}
    , jetIDDo_{true}
    ,
    // B-discriminator cut
    numbJets_{1}
    , maxbJets_{2}
    ,
    // bDiscCut_{0.9535}, // Tight cut
    bDiscCut_{is2016 ? 0.8484 : 0.8838}
    , // Medium level
    // bDiscCut_{0.5426}, // Loose cut
    bLooseDiscCut_{is2016 ? 0.5426 : 0.5803}
    , // Loose cut
    bDiscSynchCut_{0.5426}
    ,
    // C-discriminator cut
    numcJets_{1}
    , maxcJets_{1}
    ,
    // cVsLDiscCut_{0.69}, // Tight cut
    cVsLDiscCut_{-0.1}
    , // Medium level
    // cVsLDiscCut_{-0.48}, // Loose cut
    cVsBDiscCut_{0.45}
    , // Tight cut
    // cVsBDiscCut_{0.08}, // Medium level
    // cVsBDiscCut_{-0.17}, // Loose cut

    rc_{"scaleFactors/2017/RoccoR2017v0.txt"}
    , tempSmearValue_{1.0}
    , // Temporary solution to smearing propagation bug fix. A more elegant
      // solution is needed!
    lumiRunsBCDEF_{19713.888}
    , // Lumi for hip era runs
    lumiRunsGH_{16146.178}
    , // Lumi for post-hip era runs

    // Set isMC. Default is true, but it's called everytime a new dataset is
    // processed anyway.
    isMC_{true}
    ,
    // Same for trigger flag.
    triggerFlag_{}
    , isNPL_{false}
    , isZplusCR_{false}
    , doGenMassCuts_{false}
    , doGenPtCuts_{false}
    ,
    // Make cloned lepton sel tree false for now
    postLepSelTree_{nullptr}
    ,
    // Skips running trigger stuff
    skipTrigger_{false}
    ,
    // Are we making b-tag efficiency plots?
    makeBTagEffPlots_{false}
    , getBTagWeight_{false}
    ,

    // MET and mTW cuts go here.
    metCut_{0.0}
    , metDileptonCut_{50.0}
    , mTWCutSynch_{20.0}
    , TopMassCutLower_{95.}
    , TopMassCutUpper_{200.}

{
    // Space here in case other stuff needs to be done.

    std::cout << "lumi set for Runs B-F: " << lumiRunsBCDEF_ << std::endl;
    std::cout << "lumi set for Runs G-H: " << lumiRunsGH_ << std::endl;

    // If doing synchronisation., initialise that here.
    if (synchCutFlow_)
    {
        synchCutFlowHist_ = new TH1D{"synchCutFlow", "synchCutFlow", 11, 0, 11};
        synchNumEles_ = new TH1I{"synchNumEles", "synchNumEles", 11, 0, 11};
        synchNumMus_ = new TH1I{"synchNumMuos", "synchNumMuos", 11, 0, 11};
        synchMuonCutFlow_ =
            new TH1I{"synchMuonCutFlow", "synchMuonCutFlow", 11, 0, 11};
        synchCutTopMassHist_ = new TH1F{
            "synchCutTopMassHist", "synchCutTopMassHist", 200, 0., 200.};
    }

    std::cout << "\nInitialises fine" << std::endl;
    initialiseJECCors();
    std::cout << "Gets past JEC Cors" << std::endl;

    if (!is2016_)
    {
        std::cout << "\nLoad 2017 electron SFs from root file ... "
                  << std::endl;
        electronSFsFile =
            new TFile("scaleFactors/2017/"
                      "egammaEffi.txt_EGM2D_runBCDEF_"
                      "passingTight94X.root"); // Electron cut-based
                                               // Tight ID
        h_eleSFs = dynamic_cast<TH2F*>(electronSFsFile->Get("EGamma_SF2D"));
        electronRecoFile = new TFile{
            "scaleFactors/2017/"
            "egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"}; // Electron Reco
                                                               // SF
        h_eleReco = dynamic_cast<TH2F*>(electronRecoFile->Get("EGamma_SF2D"));
        std::cout << "Got 2017 electron SFs!\n" << std::endl;

        std::cout << "Load 2017 muon SFs from root file ... " << std::endl;
        muonHltFile1 = new TFile{
            "scaleFactors/2017/HLT_Mu24_EfficienciesAndSF_RunBtoF.root"};
        muonIDsFile1 = new TFile{"scaleFactors/2017/Muon_RunBCDEF_SF_ID.root"};
        muonIsoFile1 = new TFile{"scaleFactors/2017/Muon_RunBCDEF_SF_ISO.root"};
        muonHltFile1->cd("IsoMu27_PtEtaBins"); // Single Muon HLT SF
        h_muonHlt1 = dynamic_cast<TH2F*>(muonHltFile1->Get(
            "IsoMu27_PtEtaBins/abseta_pt_ratio")); // Single Muon
                                                   // HLT SF
        h_muonIDs1 = dynamic_cast<TH2D*>(muonIDsFile1->Get(
            "NUM_TightID_DEN_genTracks_pt_abseta")); // Tight ID
        h_muonPFiso1 = dynamic_cast<TH2D*>(muonIsoFile1->Get(
            "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")); // Tight ID
        std::cout << "Got 2017 muon SFs!\n" << std::endl;
    }
    else
    {
        std::cout << "\nLoad 2016 electron SFs from root file ... "
                  << std::endl;
        electronHltFile = new TFile(
            "scaleFactors/2016/"
            "HLT_Ele32_eta2p1_WPTight_Gsf_FullRunRange.root"); // Single
                                                               // Electron
                                                               // HLT SF
        h_eleHlt = dynamic_cast<TH2F*>(electronHltFile->Get("SF"));
        electronSFsFile = new TFile(
            "scaleFactors/2016/egammaEffi_Tight_80X.txt_EGM2D.root"); // Electron
                                                                      // cut-based
                                                                      // Tight
                                                                      // ID
        h_eleSFs = dynamic_cast<TH2F*>(electronSFsFile->Get("EGamma_SF2D"));
        electronRecoFile = new TFile{
            "scaleFactors/2016/egammaRecoEffi.txt_EGM2D.root"}; // Electron Reco
                                                                // SF
        h_eleReco = dynamic_cast<TH2F*>(electronRecoFile->Get("EGamma_SF2D"));
        std::cout << "Got 2016 electron SFs!\n" << std::endl;

        std::cout << "Load 2016 muon SFs from root file ... " << std::endl;
        muonHltFile1 =
            new TFile{"scaleFactors/2016/"
                      "HLT_Mu24_EfficienciesAndSF_RunBtoF.root"}; // RunsB-F
                                                                  // -
                                                                  // pre-HIP
                                                                  // fix
        muonHltFile2 =
            new TFile{"scaleFactors/2016/"
                      "HLT_Mu24_EfficienciesAndSF_RunGtoH.root"}; // RunsB-F
                                                                  // -
                                                                  // pre-HIP
                                                                  // fix
        muonIDsFile1 = new TFile{
            "scaleFactors/2016/MuonID_EfficienciesAndSF_BCDEF.root"}; // RunsB-F
                                                                      // -
                                                                      // pre-HIP
                                                                      // fix
        muonIDsFile2 = new TFile{
            "scaleFactors/2016/MuonID_EfficienciesAndSF_GH.root"}; // RunsG-H -
                                                                   // post-HIP
                                                                   // fix
        muonIsoFile1 = new TFile{
            "scaleFactors/2016/MuonISO_EfficienciesAndSF_BCDEF.root"}; // RunsB-F
                                                                       // -
                                                                       // pre-HIP
                                                                       // fix
        muonIsoFile2 = new TFile{
            "scaleFactors/2016/MuonISO_EfficienciesAndSF_GH.root"}; // RunsG-H -
                                                                    // post-HIP
                                                                    // fix

        muonHltFile1->cd(
            "IsoMu24_OR_IsoTkMu24_PtEtaBins"); // Single Muon HLT SF
        muonHltFile2->cd(
            "IsoMu24_OR_IsoTkMu24_PtEtaBins"); // Single Muon HLT SF
        h_muonHlt1 = dynamic_cast<TH2F*>(muonHltFile1->Get(
            "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio")); // Single Muon
                                                                // HLT SF
        h_muonHlt2 = dynamic_cast<TH2F*>(muonHltFile2->Get(
            "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio")); // Single Muon
                                                                // HLT SF
        muonIDsFile1->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"); // Tight ID
        muonIDsFile2->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"); // Tight ID
        dynamic_cast<TH2F*>(
            muonIDsFile1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/"
                              "abseta_pt_ratio"))
            ->Copy(*h_muonIDs1); // Tight ID
        h_muonIDs2 = dynamic_cast<TH2F*>(
            muonIDsFile2->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/"
                              "abseta_pt_ratio")); // Tight
                                                   // ID
        muonIsoFile1->cd("TightISO_TightID_pt_eta"); // Tight Iso
        muonIsoFile2->cd("TightISO_TightID_pt_eta"); // Tight Iso
        dynamic_cast<TH2F*>(
            muonIsoFile1->Get("TightISO_TightID_pt_eta/abseta_pt_ratio"))
            ->Copy(*h_muonPFiso1); // Tight Iso
        h_muonPFiso2 = dynamic_cast<TH2F*>(muonIsoFile2->Get(
            "TightISO_TightID_pt_eta/abseta_pt_ratio")); // Tight Iso
        std::cout << "Got 2016 muon SFs!\n" << std::endl;
    }

    // Setup bTag calibration code (2016/2017)
    // bTag calib code
    calib2016 = BTagCalibration("CSVv2", "scaleFactors/2016/CSVv2.csv");
    calib2017 =
        BTagCalibration("CSVv2", "scaleFactors/2017/CSVv2_94XSF_V2_B_F.csv");

    // udsg jets
    lightReader = BTagCalibrationReader(
        BTagEntry::OP_TIGHT, "central", {"up", "down"}); // operating point
    // c/b jets
    charmReader = BTagCalibrationReader(
        BTagEntry::OP_TIGHT, "central", {"up", "down"}); // central
    beautyReader = BTagCalibrationReader(
        BTagEntry::OP_TIGHT, "central", {"up", "down"}); // central

    // if doing bTag SFs, load stuff here ...
    // N.B. 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG
    if (getBTagWeight_)
    {
        if (!is2016_)
        {
            lightReader.load(calib2017, BTagEntry::FLAV_UDSG, "incl");
            charmReader.load(calib2017, BTagEntry::FLAV_C, "mujets");
            beautyReader.load(calib2017, BTagEntry::FLAV_B, "mujets");
        }
        else
        {
            lightReader.load(calib2016, BTagEntry::FLAV_UDSG, "incl");
            charmReader.load(calib2016, BTagEntry::FLAV_C, "mujets");
            beautyReader.load(calib2016, BTagEntry::FLAV_B, "mujets");
        }
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

    if (synchCutFlow_)
    {
        delete synchCutFlowHist_;
        delete synchNumEles_;
        delete synchNumMus_;
        delete synchMuonCutFlow_;
        delete synchCutTopMassHist_;
        if (makeEventDump_)
        {
            topMassEventDump_.close();
            step0EventDump_.close();
            step2EventDump_.close();
            step4EventDump_.close();
            step6EventDump_.close();
            step9EventDump_.close();
        }
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
    // numcJets_ = jets["numcJets"].as<unsigned>();
    // maxcJets_ = jets["maxcJets"].as<unsigned>();

    std::cerr << "And so it's looking for " << numTightMu_ << " muons and "
              << numTightEle_ << " electrons" << std::endl;

    if (makeEventDump_ && synchCutFlow_)
    {
        topMassEventDump_.open("topMassEventDump" + postfixName_ + ".txt");
        step0EventDump_.open("step0EventDump" + postfixName_ + ".txt");
        step2EventDump_.open("step2EventDump" + postfixName_ + ".txt");
        step4EventDump_.open("step4EventDump" + postfixName_ + ".txt");
        step6EventDump_.open("step6EventDump" + postfixName_ + ".txt");
        step9EventDump_.open("step9EventDump" + postfixName_ + ".txt");
    }
}

bool Cuts::makeCuts(AnalysisEvent& event,
                    float& eventWeight,
                    std::map<std::string, std::shared_ptr<Plots>>& plotMap,
                    TH1D& cutFlow,
                    const int systToRun)
{
    // If we're doing synchronisation, do this function.
    if (synchCutFlow_)
    {
        return synchCuts(event, eventWeight);
    }

    // For smearing temp solution
    event.jetSmearValue = {1.0};

    // If emu and dilepton - doing ttbar background estimation

    if (cutConfTrigLabel_.find("d") != std::string::npos)
    {
        return ttbarCuts(event, eventWeight, plotMap, cutFlow, systToRun);
    }

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

    // Make lepton cuts. Does the inverted iso cuts if necessary.
    if (!makeLeptonCuts(event, eventWeight, plotMap, cutFlow, systToRun))
    {
        return false;
    }

    std::pair<std::vector<int>, std::vector<float>> jetInfo;
    jetInfo = makeJetCuts(event, systToRun, eventWeight);
    event.jetIndex = jetInfo.first;
    event.jetSmearValue = jetInfo.second;

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

    if (!isFCNC_)
    { // Do wMass stuff
        float invWmass{0.};
        invWmass = getWbosonQuarksCand(event, event.jetIndex, systToRun);

        // Debug chi2 cut
        //   float topMass = getTopMass(event);
        //   float topTerm = ( topMass-173.21 )/30.0;
        //   float wTerm = ( (event.wPairQuarks.first +
        //   event.wPairQuarks.second).M() - 80.3585 )/8.0;

        //   float chi2 = topTerm*topTerm + wTerm*wTerm;
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
    }

    if (isFCNC_) // Do FCNC stuff
    {
        // Leading jet cannot be b-tagged
        if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[0]
            > bDiscCut_)
        {
            return false;
        }

        if (isCtag_) // Do cTagging
        {
            //          event.cTagIndex = makeCCuts(event,event.jetIndex);

            if (event.cTagIndex.size() < numcJets_
                || event.cTagIndex.size() > maxJets_)
            {
                return false;
            }
            if (doPlots_ || fillCutFlow_)
            {
                if (doPlots_)
                {
                    plotMap["cTag"]->fillAllPlots(event, eventWeight);
                }

                cutFlow.Fill(4.5, eventWeight);
            }
        }
    }

    TLorentzVector tempMet;
    tempMet.SetPtEtaPhiE(
        event.metPF2PATPt, 0, event.metPF2PATPhi, event.metPF2PATEt);
    double mtw{
        std::sqrt(2 * event.metPF2PATPt * event.wLepton.Pt()
                  * (1 - std::cos(event.metPF2PATPhi - event.wLepton.Phi())))};
    return true;
}

// Make lepton cuts. Will become customisable in a config later on.
bool Cuts::makeLeptonCuts(
    AnalysisEvent& event,
    float& eventWeight,
    std::map<std::string, std::shared_ptr<Plots>>& plotMap,
    TH1D& cutFlow,
    const int syst,
    const bool isControl)
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
    { // if ee channel
        //        std::cout << "Is ele 1/2 prompt? : " <<
        //        event.genElePF2PATPromptFinalState[event.zPairIndex.first]
        //        << "/" <<
        //        event.genElePF2PATPromptFinalState[event.zPairIndex.second]
        //        << std::endl;
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
    if (postLepSelTree_ && !synchCutFlow_)
    {
        postLepSelTree_->Fill();
    }

    // DO ROCHESTER CORRECTIONS HERE
    std::vector<float> SFs = {1.0};

    for (auto muonIt = event.muonIndexTight.begin();
         muonIt != event.muonIndexTight.end();
         muonIt++)
    {
        float tempSF{1.0};
        if (event.muonPF2PATPt[*muonIt] < 200.)
        { // Only doing Rochester for 2016, method only applicable for pT < 200
            if (isMC_)
            {
                tempSF = rc_.kScaleAndSmearMC(
                    event.muonPF2PATCharge[*muonIt],
                    event.muonPF2PATPt[*muonIt],
                    event.muonPF2PATEta[*muonIt],
                    event.muonPF2PATPhi[*muonIt],
                    event.muonPF2PATTkLysWithMeasurements[*muonIt],
                    gRandom->Rndm(),
                    gRandom->Rndm(),
                    0,
                    0);
            }
            else
            {
                tempSF = rc_.kScaleDT(event.muonPF2PATCharge[*muonIt],
                                      event.muonPF2PATPt[*muonIt],
                                      event.muonPF2PATEta[*muonIt],
                                      event.muonPF2PATPhi[*muonIt],
                                      0,
                                      0);
            }
        }
        SFs.emplace_back(tempSF);
    }

    event.muonMomentumSF = SFs;

    // FINISH ROCHESTER CORRECTIONS BIT

    // Should I make it return which leptons are the zMass candidate? Probably.
    float invZmass{9999.};
    if (!getDileptonZCand(
            event, event.electronIndexTight, event.muonIndexTight))
    {
        return false;
    }

    eventWeight *= getLeptonWeight(event, syst);

    if (doPlots_ || fillCutFlow_)
    {
        event.jetIndex = makeJetCuts(event, syst, eventWeight, false).first;
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
        double eeWeight{1.0}, mumuWeight{1.0}, emuWeight{1.0};
        if (invZMassCut_ == 20. && invWMassCut_ == 20.)
        {
            eeWeight = 0.926;
            mumuWeight = 1.114;
            emuWeight = 1.489;
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
        && !isControl)
    {
        return false;
    }
    //  if (std::abs( (event.zPairLeptons.first +
    //  event.zPairLeptons.second).M() -91.1 ) <= invZMassCut_ && !isControl)
    //  return false;
    if (std::abs(invZmass) < 106 && isControl)
    {
        return false;
    }

    if (doPlots_ || fillCutFlow_)
    {
        event.jetIndex = makeJetCuts(event, syst, eventWeight, false).first;
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
        if (event.elePF2PATCutIdTight[i])
        {
            electrons.emplace_back(i);
        }
    }
    return electrons;
}

std::vector<int> Cuts::getLooseEles(const AnalysisEvent& event) const
{
    std::vector<int> electrons;
    for (int i{0}; i < event.numElePF2PAT; i++)
    {
        if (event.elePF2PATCutIdLoose[i])
        {
            electrons.emplace_back(i);
        }
    }
    return electrons;
}

std::vector<int> Cuts::getTightMuons(const AnalysisEvent& event) const
{
    std::vector<int> muons;
    for (int i{0}; i < event.numMuonPF2PAT; i++)
    {
        if (event.muonPF2PATTightCutId[i])
        {
            muons.emplace_back(i);
        }
    }
    return muons;
}

std::vector<int> Cuts::getLooseMuons(const AnalysisEvent& event) const
{
    std::vector<int> muons;
    for (int i{0}; i < event.numMuonPF2PAT; i++)
    {
        if (event.muonPF2PATLooseCutId[i])
        {
            muons.emplace_back(i);
        }
    }
    return muons;
}

float Cuts::getTrileptonZCand(AnalysisEvent& event,
                              const std::vector<int> electrons,
                              const std::vector<int> muons) const
{
    float closestMass{9999.};
    // Use electrons if there are at least 2. Otherwise use muons.
    if (electrons.size() > 1)
    { // electrons.size() == number of electrons for selected channel
        for (unsigned i{0}; i < electrons.size(); i++)
        {
            for (unsigned j{i + 1}; j < electrons.size(); j++)
            {
                if (event.elePF2PATCharge[electrons[i]]
                        * event.elePF2PATCharge[electrons[j]]
                    > 0)
                {
                    continue;
                }
                TLorentzVector lepton1{event.elePF2PATPX[electrons[0]],
                                       event.elePF2PATPY[electrons[0]],
                                       event.elePF2PATPZ[electrons[0]],
                                       event.elePF2PATE[electrons[0]]};
                TLorentzVector lepton2{event.elePF2PATPX[electrons[1]],
                                       event.elePF2PATPY[electrons[1]],
                                       event.elePF2PATPZ[electrons[1]],
                                       event.elePF2PATE[electrons[1]]};
                // TLorentzVector
                // lepton1{event.elePF2PATGsfPx[electrons[i]],event.elePF2PATGsfPy[electrons[i]],event.elePF2PATGsfPz[electrons[i]],event.elePF2PATGsfE[electrons[i]]};
                // TLorentzVector
                // lepton2{event.elePF2PATGsfPx[electrons[j]],event.elePF2PATGsfPy[electrons[j]],event.elePF2PATGsfPz[electrons[j]],event.elePF2PATGsfE[electrons[j]]};
                double invMass{(lepton1 + lepton2).M() - 91.1};
                if (std::abs(invMass) < std::abs(closestMass))
                {
                    // set up the tlorentz vectors in the event. For plotting
                    // and jazz.
                    event.zPairLeptons.first =
                        lepton1.Pt() > lepton2.Pt() ? lepton1 : lepton2;
                    event.zPairIndex.first = lepton1.Pt() > lepton2.Pt()
                                                 ? electrons[i]
                                                 : electrons[j];
                    event.zPairRelIso.first =
                        lepton1.Pt() > lepton2.Pt()
                            ? event.elePF2PATComRelIsoRho[electrons[i]]
                            : event.elePF2PATComRelIsoRho[electrons[j]];
                    event.zPairRelIso.second =
                        lepton1.Pt() > lepton2.Pt()
                            ? event.elePF2PATComRelIsoRho[electrons[j]]
                            : event.elePF2PATComRelIsoRho[electrons[i]];
                    event.zPairLeptons.second =
                        lepton1.Pt() > lepton2.Pt() ? lepton2 : lepton1;
                    event.zPairIndex.second = lepton1.Pt() > lepton2.Pt()
                                                  ? electrons[j]
                                                  : electrons[i];
                    closestMass = invMass;
                    // Now set up W lepton ...
                    if (electrons.size() == 2)
                    {
                        event.wLepton =
                            TLorentzVector(event.muonPF2PATPX[muons[0]],
                                           event.muonPF2PATPY[muons[0]],
                                           event.muonPF2PATPZ[muons[0]],
                                           event.muonPF2PATE[muons[0]]);
                        event.wLeptonRelIso =
                            event.muonPF2PATComRelIsodBeta[muons[0]];
                        event.wLepIndex = muons[0];
                    }
                    else
                    {
                        for (unsigned k{0}; k < electrons.size(); k++)
                        {
                            if (k == i || k == j)
                            {
                                continue;
                            }
                            event.wLepton =
                                TLorentzVector(event.elePF2PATPX[electrons[k]],
                                               event.elePF2PATPY[electrons[k]],
                                               event.elePF2PATPZ[electrons[k]],
                                               event.elePF2PATE[electrons[k]]);
                            event.wLeptonRelIso =
                                event.elePF2PATComRelIsoRho[electrons[k]];
                            event.wLepIndex = electrons[k];
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (unsigned i{0}; i < muons.size(); i++)
        {
            for (unsigned j{i + 1}; j < muons.size(); j++)
            {
                if (event.muonPF2PATCharge[muons[i]]
                        * event.muonPF2PATCharge[muons[j]]
                    > 0)
                {
                    continue;
                }
                TLorentzVector lepton1{event.muonPF2PATPX[muons[i]],
                                       event.muonPF2PATPY[muons[i]],
                                       event.muonPF2PATPZ[muons[i]],
                                       event.muonPF2PATE[muons[i]]};
                TLorentzVector lepton2{event.muonPF2PATPX[muons[j]],
                                       event.muonPF2PATPY[muons[j]],
                                       event.muonPF2PATPZ[muons[j]],
                                       event.muonPF2PATE[muons[j]]};
                double invMass{(lepton1 + lepton2).M() - 91};
                if (std::abs(invMass) < std::abs(closestMass))
                {
                    // set up the tlorentz vectors in the event. For plotting
                    // and jazz.
                    event.zPairLeptons.first = lepton1;
                    event.zPairIndex.first = muons[i];
                    event.zPairLeptons.second = lepton2;
                    event.zPairIndex.second = muons[j];
                    event.zPairRelIso.first =
                        event.muonPF2PATComRelIsodBeta[muons[i]];
                    event.zPairRelIso.second =
                        event.muonPF2PATComRelIsodBeta[muons[j]];
                    closestMass = invMass;
                    // Now set up W lepton
                    if (muons.size() == 2)
                    {
                        event.wLepton =
                            TLorentzVector(event.elePF2PATPX[electrons[0]],
                                           event.elePF2PATPY[electrons[0]],
                                           event.elePF2PATPZ[electrons[0]],
                                           event.elePF2PATE[electrons[0]]);
                        event.wLeptonRelIso =
                            event.elePF2PATComRelIsoRho[electrons[0]];
                        event.wLepIndex = electrons[0];
                    }
                    else
                    {
                        for (unsigned k{0}; k < muons.size(); k++)
                        {
                            if (k == i || k == j)
                            {
                                continue;
                            }
                            event.wLepton =
                                TLorentzVector(event.muonPF2PATPX[muons[k]],
                                               event.muonPF2PATPY[muons[k]],
                                               event.muonPF2PATPZ[muons[k]],
                                               event.muonPF2PATE[muons[k]]);
                            event.wLeptonRelIso =
                                event.muonPF2PATComRelIsodBeta[muons[k]];
                            event.wLepIndex = muons[k];
                        }
                    }
                }
            }
        }
    }
    return closestMass;
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

        TLorentzVector lepton1{event.elePF2PATPX[electrons[0]],
                               event.elePF2PATPY[electrons[0]],
                               event.elePF2PATPZ[electrons[0]],
                               event.elePF2PATE[electrons[0]]};
        TLorentzVector lepton2{event.elePF2PATPX[electrons[1]],
                               event.elePF2PATPY[electrons[1]],
                               event.elePF2PATPZ[electrons[1]],
                               event.elePF2PATE[electrons[1]]};

        // TLorentzVector
        // lepton1{event.elePF2PATGsfPx[electrons[0]],event.elePF2PATGsfPy[electrons[0]],event.elePF2PATGsfPz[electrons[0]],event.elePF2PATGsfE[electrons[0]]};
        // TLorentzVector
        // lepton2{event.elePF2PATGsfPx[electrons[1]],event.elePF2PATGsfPy[electrons[1]],event.elePF2PATGsfPz[electrons[1]],event.elePF2PATGsfE[electrons[1]]};

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

        lepton1 *= event.muonMomentumSF[0];
        lepton2 *= event.muonMomentumSF[1];

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

        lepton2 *= event.muonMomentumSF[0];

        event.zPairLeptons.first = lepton1;
        event.zPairLeptons.second = lepton2;
        return true;
    }
    else
    {
        return false; // Not dilepton candidate if this is the case ...
    }
}

float Cuts::getWbosonQuarksCand(AnalysisEvent& event,
                                const std::vector<int> jets,
                                const int syst)
{
    float closestWmass{9999.};
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
                TLorentzVector jetVec1{getJetLVec(event, jets[k], syst, false)};
                TLorentzVector jetVec2{getJetLVec(event, jets[l], syst, false)};

                TLorentzVector wQuark1{
                    jetVec1.Px(), jetVec1.Py(), jetVec1.Pz(), jetVec1.E()};
                TLorentzVector wQuark2{
                    jetVec2.Px(), jetVec2.Py(), jetVec2.Pz(), jetVec2.E()};
                double invWbosonMass{(wQuark1 + wQuark2).M() - 80.385};

                if (std::abs(invWbosonMass) < std::abs(closestWmass))
                {
                    event.wPairQuarks.first =
                        wQuark1.Pt() > wQuark2.Pt() ? wQuark1 : wQuark2;
                    event.wPairIndex.first =
                        wQuark1.Pt() > wQuark2.Pt() ? jets[k] : jets[l];
                    event.wPairQuarks.second =
                        wQuark1.Pt() > wQuark2.Pt() ? wQuark2 : wQuark1;
                    event.wPairIndex.second =
                        wQuark1.Pt() > wQuark2.Pt() ? jets[l] : jets[k];
                    // 	    std::cout << "wQuarks pT: " <<
                    // event.wPairQuarks.first.Pt() << "/" <<
                    // event.wPairQuarks.second.Pt() << std::endl;
                    //	    std::cout << "closestMass/invWmass: " <<
                    // closestWmass << "/" << invWbosonMass << std::endl;
                    closestWmass = invWbosonMass;
                }
            }
        }
    }
    return closestWmass;
}

// tW synch - code to get dilepton pairs correct for the tW topology
std::vector<std::pair<int, int>>
    Cuts::getSynchDileptonCandidates(const AnalysisEvent& event,
                                     const std::vector<int> eles,
                                     const std::vector<int> mus) const
{
    std::vector<std::pair<int, int>> dileptonPairs;
    std::vector<int> avaliableElectrons = eles;
    std::vector<int> avaliableMuons = mus;
    bool exhaustedOppSignOptions{false};

    if (numTightEle_ == 2)
    {
        while (avaliableElectrons.size() > 1 && !exhaustedOppSignOptions)
        {
            double leadingPt{-1.};
            bool noneFound = true;
            int lepIndex_1, lepIndex_2;
            for (unsigned i{0}; i < avaliableElectrons.size(); i++)
            {
                for (unsigned j{i + 1}; j < avaliableElectrons.size(); j++)
                {
                    if (event.elePF2PATCharge[avaliableElectrons[i]]
                            * event.elePF2PATCharge[avaliableElectrons[j]]
                        < 0)
                    {
                        noneFound = false;
                    }
                    if (event.elePF2PATCharge[avaliableElectrons[i]]
                            * event.elePF2PATCharge[avaliableElectrons[j]]
                        >= 0)
                    { // Check electron pair have opposite charge
                        if (j == avaliableElectrons.size() - 1
                            && i == avaliableElectrons.size() - 2 && noneFound)
                            exhaustedOppSignOptions =
                                true; // Set flag if reached end of avaliable
                                      // leptons and no opposite sign options
                                      // are left.
                        {
                            continue;
                        }
                    }

                    TLorentzVector lepton1{
                        event.elePF2PATPX[avaliableElectrons[i]],
                        event.elePF2PATPY[avaliableElectrons[i]],
                        event.elePF2PATPZ[avaliableElectrons[i]],
                        event.elePF2PATE[avaliableElectrons[i]]};
                    TLorentzVector lepton2{
                        event.elePF2PATPX[avaliableElectrons[j]],
                        event.elePF2PATPY[avaliableElectrons[j]],
                        event.elePF2PATPZ[avaliableElectrons[j]],
                        event.elePF2PATE[avaliableElectrons[j]]};

                    double invPt{(lepton1 + lepton2).Pt()};
                    if (invPt > leadingPt)
                    {
                        lepIndex_1 = lepton1.Pt() > lepton2.Pt()
                                         ? avaliableElectrons[i]
                                         : avaliableElectrons[j];
                        lepIndex_2 = lepton1.Pt() > lepton2.Pt()
                                         ? avaliableElectrons[j]
                                         : avaliableElectrons[i];
                        leadingPt = invPt;
                    }
                }
            }
            if (!exhaustedOppSignOptions)
            {
                dileptonPairs.emplace_back(
                    std::make_pair(lepIndex_1, lepIndex_2));
                avaliableElectrons.erase(std::remove(avaliableElectrons.begin(),
                                                     avaliableElectrons.end(),
                                                     lepIndex_1),
                                         avaliableElectrons.end());
                avaliableElectrons.erase(std::remove(avaliableElectrons.begin(),
                                                     avaliableElectrons.end(),
                                                     lepIndex_2),
                                         avaliableElectrons.end());
            }
        }
    }

    else if (numTightMu_ == 2)
    {
        while (avaliableMuons.size() > 1 && !exhaustedOppSignOptions)
        {
            double leadingPt{-1.};
            bool noneFound = true;
            int lepIndex_1, lepIndex_2;
            for (unsigned i{0}; i < avaliableMuons.size(); i++)
            {
                for (unsigned j{i + 1}; j < avaliableMuons.size(); j++)
                {
                    if (event.muonPF2PATCharge[avaliableMuons[i]]
                            * event.muonPF2PATCharge[avaliableMuons[j]]
                        < 0)
                    {
                        noneFound = false;
                    }
                    if (event.muonPF2PATCharge[avaliableMuons[i]]
                            * event.muonPF2PATCharge[avaliableMuons[j]]
                        >= 0)
                    { // Chec muon pair have opposite charge
                        if (j == avaliableMuons.size() - 1
                            && i == avaliableMuons.size() - 2 && noneFound)
                            exhaustedOppSignOptions =
                                true; // Set flag if reached end of avaliable
                                      // leptons and no opposite sign options
                                      // are left.
                        {
                            continue;
                        }
                    }

                    TLorentzVector lepton1{
                        event.muonPF2PATPX[avaliableMuons[i]],
                        event.muonPF2PATPY[avaliableMuons[i]],
                        event.muonPF2PATPZ[avaliableMuons[i]],
                        event.muonPF2PATE[avaliableMuons[i]]};
                    TLorentzVector lepton2{
                        event.muonPF2PATPX[avaliableMuons[j]],
                        event.muonPF2PATPY[avaliableMuons[j]],
                        event.muonPF2PATPZ[avaliableMuons[j]],
                        event.muonPF2PATE[avaliableMuons[j]]};

                    double invPt{(lepton1 + lepton2).Pt()};

                    if (invPt > leadingPt)
                    {
                        lepIndex_1 = lepton1.Pt() > lepton2.Pt()
                                         ? avaliableMuons[i]
                                         : avaliableMuons[j];
                        lepIndex_2 = lepton1.Pt() > lepton2.Pt()
                                         ? avaliableMuons[j]
                                         : avaliableMuons[i];
                        leadingPt = invPt;
                    }
                }
            }
            if (!exhaustedOppSignOptions)
            {
                dileptonPairs.emplace_back(
                    std::make_pair(lepIndex_1, lepIndex_2));
                avaliableMuons.erase(std::remove(avaliableMuons.begin(),
                                                 avaliableMuons.end(),
                                                 lepIndex_1),
                                     avaliableMuons.end());
                avaliableMuons.erase(std::remove(avaliableMuons.begin(),
                                                 avaliableMuons.end(),
                                                 lepIndex_2),
                                     avaliableMuons.end());
            }
        }
    }

    else if (numTightEle_ == 1 && numTightMu_ == 1)
    {
        while (avaliableElectrons.size() > 0 && avaliableMuons.size() > 0
               && !exhaustedOppSignOptions)
        {
            double leadingPt{-1.};
            bool noneFound = true;
            int lepIndex_1, lepIndex_2;
            for (unsigned i{0}; i < avaliableElectrons.size(); i++)
            {
                for (unsigned j{0}; j < avaliableMuons.size(); j++)
                {
                    //          std::cout << __LINE__ << " : " << __FILE__ <<
                    //          std::endl; std::cout << "charge i*j = " <<
                    //          event.elePF2PATCharge[avaliableElectrons[i]] <<
                    //          " * " <<
                    //          event.muonPF2PATCharge[avaliableMuons[j]] << "
                    //          = " <<
                    //          event.elePF2PATCharge[avaliableElectrons[i]] *
                    //          event.muonPF2PATCharge[avaliableMuons[j]] <<
                    //          std::endl;
                    if (event.elePF2PATCharge[avaliableElectrons[i]]
                            * event.muonPF2PATCharge[avaliableMuons[j]]
                        < 0)
                    {
                        noneFound = false;
                    }
                    if (event.elePF2PATCharge[avaliableElectrons[i]]
                            * event.muonPF2PATCharge[avaliableMuons[j]]
                        >= 0)
                    { // Chec muon pair have opposite charge
                        if (i == avaliableElectrons.size() - 1
                            && j == avaliableMuons.size() - 1 && noneFound)
                            exhaustedOppSignOptions =
                                true; // Set flag if reached end of avaliable
                                      // leptons and no opposite sign options
                                      // are left.
                        {
                            continue;
                        }
                    }

                    //          std::cout << __LINE__ << " : " << __FILE__ <<
                    //          std::endl;

                    TLorentzVector lepton1{
                        event.elePF2PATPX[avaliableElectrons[i]],
                        event.elePF2PATPY[avaliableElectrons[i]],
                        event.elePF2PATPZ[avaliableElectrons[i]],
                        event.elePF2PATE[avaliableElectrons[i]]};
                    TLorentzVector lepton2{
                        event.muonPF2PATPX[avaliableMuons[j]],
                        event.muonPF2PATPY[avaliableMuons[j]],
                        event.muonPF2PATPZ[avaliableMuons[j]],
                        event.muonPF2PATE[avaliableMuons[j]]};

                    double invPt{(lepton1 + lepton2).Pt()};

                    if (invPt > leadingPt)
                    {
                        lepIndex_1 = avaliableElectrons[i];
                        lepIndex_2 = avaliableMuons[j];
                        leadingPt = invPt;
                    }
                }
            }
            if (!exhaustedOppSignOptions)
            {
                dileptonPairs.emplace_back(
                    std::make_pair(lepIndex_1, lepIndex_2));
                //        std::cout << __LINE__ << " : " << __FILE__ <<
                //        std::endl; std::cout << "lepIndex_1/lepIndex_2: " <<
                //        lepIndex_1 << "/" << lepIndex_2 << std::endl;
                avaliableElectrons.erase(std::remove(avaliableElectrons.begin(),
                                                     avaliableElectrons.end(),
                                                     lepIndex_1),
                                         avaliableElectrons.end());
                avaliableMuons.erase(std::remove(avaliableMuons.begin(),
                                                 avaliableMuons.end(),
                                                 lepIndex_2),
                                     avaliableMuons.end());
            }
        }
    }

    else
    {
        std::cout << "HOW ON EARTH DID YOU MANAGE THIS?" << std::endl;
    }

    return dileptonPairs;
}

float Cuts::getTopMass(const AnalysisEvent& event) const
{
    TLorentzVector metVec{
        event.metPF2PATPx, event.metPF2PATPy, 0, event.metPF2PATEt};
    TLorentzVector bVec(event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);
    float topMass{-1.0};
    topMass = (bVec + event.wPairQuarks.first + event.wPairQuarks.second).M();
    return topMass;
}

std::pair<std::vector<int>, std::vector<float>>
    Cuts::makeJetCuts(AnalysisEvent& event,
                      const int syst,
                      float& eventWeight,
                      const bool isProper)
{
    std::vector<int> jets;
    std::vector<float> SFs;

    float mcTag{1.};
    float mcNoTag{1.};
    float dataTag{1.};
    float dataNoTag{1.};
    // b-tagging errors
    float err1{0.};
    float err2{0.};
    float err3{0.};
    float err4{0.};

    //  std::cout << event.eventNum << std::endl << "Jets: " << std::endl;
    for (int i{0}; i < event.numJetPF2PAT; i++)
    {
        // Check η before retrieving TLorentzVector to prevent runtime error
        // when trying to apply a SF for an out of range η
        if (std::abs(event.jetPF2PATEta[i]) >= jetEta_)
        {
            continue;
        }

        // if (std::sqrt(event.jetPF2PATPx[i] * event.jetPF2PATPx[i] +
        // event.jetPF2PATPy[i] * event.jetPF2PATPy[i]) < jetPt_) continue;
        TLorentzVector jetVec{getJetLVec(event, i, syst, true)};
        // std::cout << getJECUncertainty(sqrt(jetPx*jetPx + jetPy*jetPy),
        // event.jetPF2PATEta[i],syst) << " " << syst << std::endl;
        if (jetVec.Pt() <= jetPt_)
        {
            continue;
        }

        bool jetId{true};

        // Jet ID == loose
        if (jetIDDo_ && isProper)
        {
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
                    if (event.jetPF2PATChargedHadronEnergyFraction[i] <= 0.0)
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

        if (!jetId)
        {
            continue;
        }

        double deltaLep{10000};

        /*
            // Electron cleaning for synch
            for ( unsigned int lep = 0; lep < event.electronIndexTight.size();
           lep++ ) { TLorentzVector tempVec{
           event.elePF2PATPX[event.electronIndexTight[lep]],
           event.elePF2PATPY[event.electronIndexTight[lep]],
           event.elePF2PATPZ[event.electronIndexTight[lep]],
           event.elePF2PATE[event.electronIndexTight[lep]] }; if (deltaLep >
           deltaR( tempVec.Eta(), tempVec.Phi(), jetVec.Eta(),jetVec.Phi() ) )
                deltaLep = deltaR( tempVec.Eta(), tempVec.Phi(),
           jetVec.Eta(),jetVec.Phi());
            }

            // Muon cleaning for synch
            for	( unsigned int lep {0}; lep < event.muonIndexTight.size();
           lep++ ) { TLorentzVector tempVec{
           event.muonPF2PATPX[event.muonIndexTight[lep]],event.muonPF2PATPY[event.muonIndexTight[lep]],event.muonPF2PATPZ[event.muonIndexTight[lep]],event.muonPF2PATE[event.muonIndexTight[lep]]
           }; if (deltaLep > deltaR( tempVec.Eta(), tempVec.Phi(),
           jetVec.Eta(),jetVec.Phi())) deltaLep = deltaR( tempVec.Eta(),
           tempVec.Phi(), jetVec.Eta(),jetVec.Phi());
            }
        */

        if (deltaLep > deltaR(event.zPairLeptons.first.Eta(),
                              event.zPairLeptons.first.Phi(),
                              jetVec.Eta(),
                              jetVec.Phi()))
        {
            deltaLep = deltaR(event.zPairLeptons.first.Eta(),
                              event.zPairLeptons.first.Phi(),
                              jetVec.Eta(),
                              jetVec.Phi());
        }
        if (deltaLep > deltaR(event.zPairLeptons.second.Eta(),
                              event.zPairLeptons.second.Phi(),
                              jetVec.Eta(),
                              jetVec.Phi()))
        {
            deltaLep = deltaR(event.zPairLeptons.second.Eta(),
                              event.zPairLeptons.second.Phi(),
                              jetVec.Eta(),
                              jetVec.Phi());
        }

        // std::cout << event.jetPF2PATPtRaw[i] << " " << deltaLep <<
        // std::endl;
        if (deltaLep < 0.4 && isProper)
        {
            continue; // Only start rejecting things when actually making the
                      // jet cuts!
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
        SFs.emplace_back(tempSmearValue_);
        tempSmearValue_ = 1.0; // Reset temporary smear value. Need a cleaner
                               // solution than this!

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
        float bWeight{(dataNoTag * dataTag) / (mcNoTag * mcTag)};
        if (mcNoTag == 0 || mcTag == 0 || dataNoTag == 0 || dataTag == 0
            || mcNoTag != mcNoTag || mcTag != mcTag || dataTag != dataTag
            || dataNoTag != dataNoTag)
        {
            bWeight = 1.;
        }
        float bWeightErr{float(
            std::sqrt(pow(err1 + err2, 2) + pow(err3 + err4, 2)) * bWeight)};
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

    return std::make_pair(jets, SFs);
}

std::vector<int> Cuts::makeBCuts(AnalysisEvent& event,
                                 const std::vector<int> jets,
                                 const int syst)
{
    std::vector<int> bJets;
    for (unsigned int i = 0; i != jets.size(); i++)
    {
        TLorentzVector jetVec{getJetLVec(event, jets[i], syst, false)};
        if (singleEventInfoDump_)
        {
            std::cout
                << __LINE__ << "/" << __FILE__ << ": "
                << event.jetPF2PATPtRaw[i] << " "
                << event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                       [jets[i]]
                << std::endl;
        }
        if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[jets[i]]
            <= bDiscCut_)
        {
            continue;
        }
        if (jetVec.Eta() >= 2.40)
        {
            continue;
        }
        //    if (event.jetPF2PATEta[ jets[i] ] >= 2.40) continue;
        bJets.emplace_back(i);
    }
    return bJets;
}

std::vector<int> Cuts::makeLooseBCuts(AnalysisEvent& event,
                                      const std::vector<int> jets,
                                      const int syst)
{
    std::vector<int> bJets;

    for (unsigned int i = 0; i != jets.size(); i++)
    {
        TLorentzVector jetVec{getJetLVec(event, jets[i], syst, false)};
        if (singleEventInfoDump_)
        {
            std::cout
                << __LINE__ << "/" << __FILE__ << ": "
                << event.jetPF2PATPtRaw[i] << " "
                << event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                       [jets[i]]
                << std::endl;
        }
        if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[jets[i]]
            <= bLooseDiscCut_)
        {
            continue;
        }
        if (jetVec.Eta() >= 2.40)
        {
            continue;
        }
        //    if (event.jetPF2PATEta[ jets[i] ] >= 2.40) continue;
        bJets.emplace_back(i);
    }
    return bJets;
}

/*
std::vector<int> Cuts::makeCCuts(AnalysisEvent& event, std::vector<int> jets){

std::vector<int> cJets;
for (unsigned i{0}; i < jets.size(); i++){
if (singleEventInfoDump_) std::cout << event.jetPF2PATPtRaw[jets[i]] << " "
<< event.jetPF2PATpfCombinedCvsLJetTags[jets[i]] << std::endl;
//      if (event.jetPF2PATJetCharge[jets[i]] <= 0) continue; // If a
negatively charged jet ... I.e. if not a  u or c ... if
(event.jetPF2PATpfCombinedCvsLJetTags[jets[i]] < cVsLDiscCut_) continue; // If
doesn't pass c vs light discriminator if
(event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[jets[i]] >
bDiscCut_) continue; // If a b jet, continue cJets.emplace_back(i);
}
return cJets;

}
*/

bool Cuts::triggerCuts(const AnalysisEvent& event,
                       float& eventWeight,
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

    // TRIGGER SFs
    // NB, Synch logic doesn't allow for them to be applied currently

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

    float twgt{1.0};

    if (!is2016_ && isMC_)
    { // Placeholder SFs for 2017
      // TODO: Update SFs
        // Dilepton channels
        if (channel == "ee")
        {
            if (eTrig || eeTrig)
            {
                twgt = 0.97397;
                if (syst == 1)
                {
                    twgt += 0.00050;
                }
                else if (syst == 2)
                {
                    twgt += 0.00050;
                }
            }
        }
        else if (channel == "mumu")
        {
            if (muTrig || mumuTrig)
            {
                twgt = 1.06937;
                if (syst == 1)
                {
                    twgt += 0.00097;
                }
                else if (syst == 2)
                {
                    twgt += 0.00097;
                }
            }
        }
        else if (channel == "emu")
        {
            if (muEGTrig)
            {
                twgt = 0.98562;
                if (syst == 1)
                {
                    twgt += 0.00057;
                }
                else if (syst == 2)
                {
                    twgt -= 0.00057;
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

                //        if (syst == 1) twgt += ( 0.00016 * lumiRunsBCDEF_ +
                //        0.00016 * lumiRunsGH_ ) / ( lumiRunsBCDEF_ +
                //        lumiRunsGH_ + 1.0e-06 ); if (syst == 2) twgt -= (
                //        0.00016 * lumiRunsBCDEF_ + 0.00016 * lumiRunsGH_ ) /
                //        ( lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06 );
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
        //    if ( mumuTrig && !(eeTrig || muEGTrig) ) // Original trigger
        //    logic, for double triggers only
        if ((mumuTrig || muTrig)
            && !(eeTrig || muEGTrig
                 || eTrig)) // Trigger logic for double + single triggers
        //    if ( (muTrig) && !(eeTrig || muEGTrig || eTrig) ) // Single
        //    trigger only whilst in debug mode
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

// Does the cuts required for the synchronisation.
bool Cuts::synchCuts(AnalysisEvent& event, float& eventWeight)
{
    // Trigger stuff would go first, but in MC at least (what I'm starting with)
    // I don't have that saved. Idiot.

    // This is to make some skims for faster running. Do lepSel and save some
    // files.
    if (postLepSelTree_ && synchCutFlow_)
    {
        postLepSelTree_->Fill();
    }

    if (singleEventInfoDump_)
    {
        //    std::cout << std::setprecision(6) << std::fixed;
    }

    if (makeEventDump_)
    {
        dumpToFile(event, 9);
    }

    // 2016 tZq dilepton synch

    synchCutFlowHist_->Fill(0.5, eventWeight); // Total events
    //    std::cout << std::setprecision(6) << std::fixed;

    if (!triggerCuts(event, eventWeight, 0))
    {
        return false;
    }
    if (!metFilters(event))
    {
        return false;
    }

    synchCutFlowHist_->Fill(1.5, eventWeight); // Trigger cuts - Step 0
    event.electronIndexTight = getTightEles(event);
    event.muonIndexTight = getTightMuons(event);
    synchNumEles_->Fill(event.electronIndexTight.size());
    synchNumMus_->Fill(event.muonIndexTight.size());

    if (event.electronIndexTight.size() != numTightEle_)
    {
        return false;
    }
    if (event.muonIndexTight.size() != numTightMu_)
    {
        return false;
    }
    if ((event.electronIndexTight.size() + event.muonIndexTight.size()) != 2)
    {
        return false;
    }

    synchCutFlowHist_->Fill(2.5, eventWeight); // 2 Tight Leptons - step 1
    if ((event.electronIndexTight.size() != getLooseEles(event).size())
        || (event.muonIndexTight.size() != getLooseMuons(event).size()))
    {
        return false;
    }
    synchCutFlowHist_->Fill(3.5, eventWeight); // Lepton Veto - step 2

    if (!getDileptonZCand(
            event, event.electronIndexTight, event.muonIndexTight))
    {
        return false;
    }

    if (((event.zPairLeptons.first + event.zPairLeptons.second).M() - 91.1)
        > invZMassCut_)
    {
        return false;
    }
    synchCutFlowHist_->Fill(4.5, eventWeight); // Z veto - step 3

    eventWeight *= getLeptonWeight(event, 0);

    std::pair<std::vector<int>, std::vector<float>> jetInfo;
    jetInfo = makeJetCuts(event, 0, eventWeight);
    event.jetIndex = jetInfo.first;
    event.jetSmearValue = jetInfo.second;

    if (event.jetIndex.size() < numJets_)
    {
        return false;
    }
    if (event.jetIndex.size() > maxJets_)
    {
        return false;
    }

    event.bTagIndex = makeBCuts(event, event.jetIndex, 0);
    synchCutFlowHist_->Fill(5.5, eventWeight); // jet selection - step 4

    if (event.bTagIndex.size() < numbJets_)
    {
        return false;
    }
    if (event.bTagIndex.size() > maxbJets_)
    {
        return false;
    }
    synchCutFlowHist_->Fill(6.5, eventWeight); // b-jet selection - step 5

    synchCutFlowHist_->Fill(7.5, eventWeight); // no Met cut applied

    double wMass = getWbosonQuarksCand(event, event.jetIndex, 0);

    if (std::abs(wMass) > invWMassCut_)
    {
        return false;
    }

    TLorentzVector bVec(event.jetPF2PATPx[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATPy[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATPz[event.jetIndex[event.bTagIndex[0]]],
                        event.jetPF2PATE[event.jetIndex[event.bTagIndex[0]]]);

    double topMass =
        (bVec + event.wPairQuarks.first + event.wPairQuarks.second).M();
    //    double topMass = getTopMass(event);
    //    std::cout << topMass << " / " << getTopMass(event) << std::endl;

    if (topMass > 130)
    {
        if (std::abs(event.jetPF2PATPID[event.jetIndex[event.bTagIndex[0]]])
            == 5)
        {
            std::cout
                << "jet PID/b mass/w Mass/leading b pT: "
                << event.jetPF2PATPID[event.jetIndex[event.bTagIndex[0]]]
                << " / " << bVec.M() << " / "
                << (event.wPairQuarks.first + event.wPairQuarks.second).M()
                << " / "
                << event.jetPF2PATPt[event.jetIndex[event.bTagIndex[0]]]
                << "\t " << std::endl;
        }
    }

    return true;
}

TH1D* Cuts::getSynchCutFlow() const
{
    std::cout << "Eles: " << numTightEle_ << " Muons: " << numTightMu_
              << std::endl;
    char const* names[]{"Total Events",
                        "Trigger",
                        "3 Leptons",
                        "Lepton Veto",
                        "zMass",
                        "1 jet",
                        "1 b-tag",
                        "MET",
                        "mTW",
                        "topMass",
                        "metFilters"};
    for (unsigned i{1}; i < 12; ++i)
    {
        std::cout << names[i - 1] << ": " << synchCutFlowHist_->GetBinContent(i)
                  << std::endl;
    }
    std::cout << "The number of leptons in the passing trigger events:"
              << std::endl;
    std::cout << "Leps\tEle\tMuon" << std::endl;
    for (unsigned i{0}; i < 12; i++)
    {
        std::cout << i << "\t" << synchNumEles_->GetBinContent(i) << "\t"
                  << synchNumMus_->GetBinContent(i) << std::endl;
    }
    char const* labels[] = {"In",
                            "ID",
                            "PtEtaIso",
                            "chi2",
                            "tklay",
                            "DBPV",
                            "TrackHit",
                            "MuHits",
                            "PixHits",
                            "MtStats",
                            "DZPV"};
    for (unsigned i{1}; i < 12; ++i)
    {
        std::cout << labels[i - 1] << ": \t"
                  << synchMuonCutFlow_->GetBinContent(i) << std::endl;
    }

    return synchCutFlowHist_;
    /*
      std::cout << "Eles: " << numTightEle_ << " Muons: " << numTightMu_ <<
      std::endl; char const *names[] {"Total Events","min 1 opp sign lep
      pair","non-Z pairs", "MET", "nJets >= 2","nBjets >= 1","nJets ==
      1","nBjets == 1"}; for (unsigned i{1}; i < 9; ++i){ std::cout <<
      names[i-1] << ": " << synchCutFlowHist_->GetBinContent(i) << std::endl;
      }
      std::cout << "Leps\tEle\tMuon" << std::endl;
      for (unsigned i{0}; i < 9; i++){
        std::cout << i << "\t" << synchNumEles_->GetBinContent(i) << "\t" <<
      synchNumMus_->GetBinContent(i) << std::endl;
      }
      return synchCutFlowHist_;
    */
}

// Method used for the synchronisation. Mimics the IPHC preselection skims.
int Cuts::getLooseLepsNum(const AnalysisEvent& event) const
{
    return getLooseElecs(event) + getLooseMus(event);
}

int Cuts::getLooseElecs(const AnalysisEvent& event) const
{
    int looseLeps{0};
    for (int i{0}; i < event.numElePF2PAT; i++)
    {
        if (!event.elePF2PATIsGsf[i])
        {
            continue;
        }
        TLorentzVector tempVec{event.elePF2PATPX[i],
                               event.elePF2PATPY[i],
                               event.elePF2PATPZ[i],
                               event.elePF2PATE[i]};
        if (tempVec.Pt() < 20)
        {
            continue;
        }
        if (std::abs(tempVec.Eta()) > tightEleEta_)
        {
            continue;
        }
        looseLeps++;
    }
    return looseLeps;
}

int Cuts::getLooseMus(const AnalysisEvent& event) const
{
    int looseLeps{0};
    for (int i{0}; i < event.numMuonPF2PAT; i++)
    {
        if (!event.muonPF2PATGlobalID[i] || !event.muonPF2PATTrackID[i])
        {
            continue;
        }
        if (event.muonPF2PATPt[i] < 20)
        {
            continue;
        }
        if (std::abs(event.muonPF2PATEta[i]) > 2.4)
        {
            continue;
        }
        looseLeps++;
    }
    return looseLeps;
}

std::vector<int> Cuts::getSynchEles(const AnalysisEvent& event) const
{
    std::vector<int> electrons;
    for (int i{0}; i < event.numElePF2PAT; i++)
    {
        if (!event.elePF2PATIsGsf[i])
        {
            continue;
        }
        TLorentzVector tempVec{event.elePF2PATPX[i],
                               event.elePF2PATPY[i],
                               event.elePF2PATPZ[i],
                               event.elePF2PATE[i]};
        if (tempVec.Pt() <= 20.)
        {
            continue;
        }
        if (std::abs(tempVec.Eta()) >= 2.4)
        {
            continue;
        }
        if ((std::abs(event.elePF2PATSCEta[i]) > 1.4442
             && std::abs(event.elePF2PATSCEta[i]) < 1.566))
        {
            continue;
        }
        if (event.elePF2PATCutIdTight[i] < 1)
        {
            continue;
        }
        electrons.emplace_back(i);
    }
    return electrons;
}

std::vector<int> Cuts::getSynchMus(const AnalysisEvent& event) const
{
    std::vector<int> muons;
    for (int i{0}; i < event.numMuonPF2PAT; i++)
    {
        if (!event.muonPF2PATIsPFMuon[i])
        {
            continue;
        }
        if (event.muonPF2PATPt[i] <= 20.)
        {
            continue;
        }
        if (std::abs(event.muonPF2PATEta[i]) >= 2.4)
        {
            continue;
        }

        if (event.muonPF2PATComRelIsodBeta[i] >= tightMuonRelIso_)
        {
            continue;
        }

        // Tight ID Cut
        if (!event.muonPF2PATTrackID[i])
        {
            continue;
        }
        if (!event.muonPF2PATGlobalID[i])
        {
            continue;
        }
        if (event.muonPF2PATGlbTkNormChi2[i] >= 10.)
        {
            continue;
        }
        if (event.muonPF2PATMatchedStations[i] < 2)
        {
            continue;
        }
        if (std::abs(event.muonPF2PATDBPV[i]) >= 0.2)
        {
            continue;
        }
        if (std::abs(event.muonPF2PATDZPV[i]) >= 0.5)
        {
            continue;
        }
        if (event.muonPF2PATMuonNHits[i] < 1)
        {
            continue;
        }
        if (event.muonPF2PATVldPixHits[i] < 1)
        {
            continue;
        }
        if (event.muonPF2PATTkLysWithMeasurements[i] <= 5)
        {
            continue;
        }

        muons.emplace_back(i);
    }
    return muons;
}

// First tentative attempt at doing the background isolation.
bool Cuts::invertIsoCut(AnalysisEvent& event,
                        float& eventWeight,
                        std::map<std::string, std::shared_ptr<Plots>>& plotMap,
                        TH1D& cutFlow)
{
    std::cout << "Invert Iso Cut is not avaliable for the dilepton channel."
              << std::endl;
    return false;
    // Check there are exactly 2 tight leptons with the correct isolation cut.
    event.electronIndexTight = getTightEles(event);
    event.muonIndexTight = getTightMuons(event);

    if ((event.electronIndexTight.size() + event.muonIndexTight.size()) != 2)
    {
        return false;
    }
    if (event.electronIndexTight.size() == 1)
    {
        return false;
    }

    // Check they are a valid zCandidate. (Just call the zCand method here? -
    // No, that won't actually work.)
    float invMass{100.};
    if (numTightEle_ > 1)
    {
        if (event.electronIndexTight.size() < 2)
        {
            return false;
        }
        if (event.elePF2PATCharge[event.electronIndexTight[0]]
                * event.elePF2PATCharge[event.electronIndexTight[1]]
            > 0)
        {
            return false;
        }
        TLorentzVector lep1{event.elePF2PATPX[event.electronIndexTight[0]],
                            event.elePF2PATPY[event.electronIndexTight[0]],
                            event.elePF2PATPZ[event.electronIndexTight[0]],
                            event.elePF2PATE[event.electronIndexTight[0]]};
        TLorentzVector lep2{event.elePF2PATPX[event.electronIndexTight[1]],
                            event.elePF2PATPY[event.electronIndexTight[1]],
                            event.elePF2PATPZ[event.electronIndexTight[1]],
                            event.elePF2PATE[event.electronIndexTight[1]]};
        event.zPairLeptons.first = lep1.Pt() > lep2.Pt() ? lep1 : lep2;
        event.zPairIndex.first = lep1.Pt() > lep2.Pt()
                                     ? event.electronIndexTight[0]
                                     : event.electronIndexTight[1];
        event.zPairLeptons.second = lep1.Pt() > lep2.Pt() ? lep2 : lep1;
        event.zPairIndex.second = lep1.Pt() > lep2.Pt()
                                      ? event.electronIndexTight[1]
                                      : event.electronIndexTight[0];
        invMass =
            (event.zPairLeptons.first + event.zPairLeptons.second).M() - 91.1;
    }
    else
    {
        if (event.muonIndexTight.size() < 2)
        {
            return false;
        }
        if (event.muonPF2PATCharge[event.muonIndexTight[0]]
                * event.muonPF2PATCharge[event.muonIndexTight[1]]
            > 0)
        {
            return false;
        }
        event.zPairLeptons.first =
            TLorentzVector(event.muonPF2PATPX[event.muonIndexTight[0]],
                           event.muonPF2PATPY[event.muonIndexTight[0]],
                           event.muonPF2PATPZ[event.muonIndexTight[0]],
                           event.muonPF2PATE[event.muonIndexTight[0]]);
        event.zPairIndex.first = event.muonIndexTight[0];
        event.zPairLeptons.second =
            TLorentzVector(event.muonPF2PATPX[event.muonIndexTight[1]],
                           event.muonPF2PATPY[event.muonIndexTight[1]],
                           event.muonPF2PATPZ[event.muonIndexTight[1]],
                           event.muonPF2PATE[event.muonIndexTight[1]]);
        event.zPairIndex.second = event.muonIndexTight[1];
        invMass =
            (event.zPairLeptons.first + event.zPairLeptons.second).M() - 91.1;
    }

    // Get rev iso candidates.
    std::vector<int> invIsoEle{getInvIsoEles(event)};
    std::vector<int> invIsoMus{getInvIsoMuons(event)};

    // Check we have the right number of leptons.
    if (event.electronIndexTight.size() + invIsoEle.size() != numTightEle_)
    {
        return false;
    }
    if (event.muonIndexTight.size() + invIsoMus.size() != numTightMu_)
    {
        return false;
    }

    // Debugging
    /*  std::cout << event.numElePF2PAT << " " <<
    event.electronIndexTight.size() << " " << invIsoEle.size(); std::cout << "
    tight index: "; for (unsigned i = 0; i < event.electronIndexTight.size();
    i++) std:: cout << " " << event.electronIndexTight[i]; std::cout << " inv
    index: "; for (unsigned i = 0; i < invIsoEle.size(); i++) std::cout << " "
    << invIsoEle[i] << " " << "relIso: " <<
    event.elePF2PATComRelIsoRho[invIsoEle[i]]/event.elePF2PATPT[invIsoEle[i]]
    ; std::cout << std::endl;*/

    // Put extra lepton into W boson thing.
    if (invIsoEle.size() == 1)
    {
        event.wLepton = TLorentzVector(event.elePF2PATPX[invIsoEle[0]],
                                       event.elePF2PATPY[invIsoEle[0]],
                                       event.elePF2PATPZ[invIsoEle[0]],
                                       event.elePF2PATE[invIsoEle[0]]);
        event.wLepIndex = invIsoEle[0];
        event.wLeptonRelIso = event.elePF2PATComRelIsoRho[invIsoEle[0]];
    }
    else
    {
        event.wLepton = TLorentzVector(event.muonPF2PATPX[invIsoMus[0]],
                                       event.muonPF2PATPY[invIsoMus[0]],
                                       event.muonPF2PATPZ[invIsoMus[0]],
                                       event.muonPF2PATE[invIsoMus[0]]);
        event.wLepIndex = invIsoMus[0];
        event.wLeptonRelIso = event.muonPF2PATComRelIsodBeta[invIsoMus[0]];
    }

    if (doPlots_)
    {
        plotMap["lepSel"]->fillAllPlots(event, eventWeight);
    }
    if (doPlots_ || fillCutFlow_)
    {
        cutFlow.Fill(0.5, eventWeight);
    }

    if (std::abs(invMass) > invZMassCut_)
    {
        return false;
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

std::vector<int> Cuts::getInvIsoEles(const AnalysisEvent& event) const
{
    std::vector<int> electrons;
    for (int i{0}; i < event.numElePF2PAT; i++)
    {
        bool same{false};
        /*    for (int j = 0; j < event.electronIndexTight.size(); j++){
          if (i == event.electronIndexTight[j]) same = true;
          }*/
        if (same)
        {
            continue;
        }
        if (!event.elePF2PATIsGsf[i])
        {
            continue;
        }
        TLorentzVector tempVec{event.elePF2PATPX[i],
                               event.elePF2PATPY[i],
                               event.elePF2PATPZ[i],
                               event.elePF2PATE[i]};
        if (tempVec.Pt() < tightElePt_)
        {
            continue;
        }
        if (std::abs(tempVec.Eta()) > tightEleEta_)
        {
            continue;
        }
        if (std::abs(event.elePF2PATBeamSpotCorrectedTrackD0[i]) > tightEled0_)
        {
            continue;
        }
        if (event.elePF2PATMissingInnerLayers[i] > tightEleMissLayers_)
        {
            continue;
        }
        if (!event.elePF2PATPhotonConversionVeto[i] && tightEleCheckPhotonVeto_)
        {
            continue;
        }
        if (event.elePF2PATComRelIsoRho[i] < tightEleRelIso_)
        {
            continue;
        }

        electrons.emplace_back(i);
    }
    return electrons;
}

std::vector<int> Cuts::getInvIsoMuons(const AnalysisEvent& event) const
{
    std::vector<int> muons;
    for (int i{0}; i < event.numMuonPF2PAT; i++)
    {
        if (!event.muonPF2PATGlobalID[i] && !event.muonPF2PATTrackID[i])
        {
            continue;
        }
        if (event.muonPF2PATPt[i] < tightMuonPt_)
        {
            continue;
        }
        if (std::abs(event.muonPF2PATEta[i]) > tightMuonEta_)
        {
            continue;
        }
        if (event.muonPF2PATComRelIsodBeta[i] < tightMuonRelIso_)
        {
            continue;
        }

        // Do a little test of muon id stuff here.
        if (event.muonPF2PATChi2[i] / event.muonPF2PATNDOF[i] > 10.)
        {
            continue;
        }
        if (std::abs(event.muonPF2PATDBInnerTrackD0[i]) > 0.2)
        {
            continue;
        }
        if (event.muonPF2PATNChambers[i] < 2)
        {
            continue;
        }

        muons.emplace_back(i);
    }
    return muons;
}

bool Cuts::ttbarCuts(AnalysisEvent& event,
                     float& eventWeight,
                     std::map<std::string, std::shared_ptr<Plots>>& plotMap,
                     TH1D& cutFlow,
                     const int systToRun)
{
    if (!skipTrigger_ && !triggerCuts(event, eventWeight, systToRun))
    {
        return false;
    }

    if (!metFilters(event))
    {
        return false;
    }

    if (!(makeLeptonCuts(
            event, eventWeight, plotMap, cutFlow, systToRun, true)))
    {
        return false;
    }

    std::pair<std::vector<int>, std::vector<float>> jetInfo;
    jetInfo = makeJetCuts(event, systToRun, eventWeight);
    event.jetIndex = jetInfo.first;
    event.jetSmearValue = jetInfo.second;

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

    float invWmass{0.};
    invWmass = getWbosonQuarksCand(event, event.jetIndex, systToRun);

    // Debug chi2 cut
    //   float topMass = getTopMass(event);
    //   float topTerm = ( topMass-173.21 )/30.0;
    //   float wTerm = ( (event.wPairQuarks.first +
    //   event.wPairQuarks.second).M() - 80.3585 )/8.0;

    //   float chi2 = topTerm*topTerm + wTerm*wTerm;
    //   if ( chi2 < 2.0 && chi2 > 7.0 ) return false; // control region
    //   if ( chi2 >= 2.0 ) return false; //signal region

    if (std::abs(invWmass) > invWMassCut_)
    {
        return false;
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

// For synchronisation I am dumping the lepton information here.
void Cuts::dumpLeptonInfo(AnalysisEvent& event)
{
    //  std::cout << std::setprecision(6) << std::fixed;
    //  std::cerr << std::setprecision(6) << std::fixed;
    // Dump electron info first
    event.electronIndexTight = getTightEles(event);
    event.muonIndexTight = getTightMuons(event);
    std::cout << "Electrons: found " << event.electronIndexTight.size()
              << std::endl;
    for (unsigned i{0}; i < event.electronIndexTight.size(); i++)
    {
        TLorentzVector tempVec{event.elePF2PATPX[event.electronIndexTight[i]],
                               event.elePF2PATPY[event.electronIndexTight[i]],
                               event.elePF2PATPZ[event.electronIndexTight[i]],
                               event.elePF2PATE[event.electronIndexTight[i]]};
        //    std::cout << i << " | " <<
        //    event.elePF2PATPT[event.electronIndexTight[i]];
        std::cout << " | " << tempVec.Pt();
        // std::cout << " | " <<
        // event.elePF2PATEta[event.electronIndexTight[i]];
        std::cout << " | " << tempVec.Eta();
        std::cout << " | " << tempVec.Phi();
        std::cout << " | " << event.elePF2PATPhi[event.electronIndexTight[i]];
        std::cout << " | " << event.elePF2PATD0PV[event.electronIndexTight[i]];
        std::cout
            << " | "
            << event.elePF2PATMissingInnerLayers[event.electronIndexTight[i]];
        std::cout << " | "
                  << event.elePF2PATComRelIsoRho[event.electronIndexTight[i]];
        std::cout
            << " | "
            << event.elePF2PATPhotonConversionVeto[event.electronIndexTight[i]];
        std::cout << " | "
                  << event.elePF2PATCharge[event.electronIndexTight[i]];
        /*std::cout << " | " <<
        event.elePF2PATRhoIso[event.electronIndexTight[i]]; std::cout << " | "
        << event.elePF2PATAEff03[event.electronIndexTight[i]]; std::cout << "
        | " << event.elePF2PATChHadIso[event.electronIndexTight[i]]; std::cout
        << " | " << event.elePF2PATNtHadIso[event.electronIndexTight[i]];
        std::cout << " | " <<
        event.elePF2PATGammaIso[event.electronIndexTight[i]]; std::cout << " |
        " << event.elePF2PATComRelIsodBeta[event.electronIndexTight[i]];
        std::cout << " | " <<
        event.elePF2PATComRelIso[event.electronIndexTight[i]]; */
        std::cout << std::endl;
    }

    std::cout << "Muons: found " << event.muonIndexTight.size() << std::endl;
    for (unsigned i{0}; i < event.muonIndexTight.size(); i++)
    {
        std::cout << i;
        std::cout << " | " << event.muonPF2PATPt[event.muonIndexTight[i]];
        std::cout << " | " << event.muonPF2PATEta[event.muonIndexTight[i]];
        std::cout << " | " << event.muonPF2PATPhi[event.muonIndexTight[i]];
        std::cout << " | "
                  << event.muonPF2PATChi2[event.muonIndexTight[i]]
                         / event.muonPF2PATNDOF[event.muonIndexTight[i]];
        std::cout
            << " | "
            << event.muonPF2PATTkLysWithMeasurements[event.muonIndexTight[i]];
        std::cout << " | "
                  << event.muonPF2PATMuonNHits[event.muonIndexTight[i]];
        std::cout << " | " << event.muonPF2PATDBPV[event.muonIndexTight[i]];
        std::cout << " | "
                  << event.pvZ - event.muonPF2PATVertZ[event.muonIndexTight[i]];
        std::cout << " | "
                  << event.muonPF2PATVldPixHits[event.muonIndexTight[i]];
        std::cout << " | "
                  << event.muonPF2PATMatchedStations[event.muonIndexTight[i]];
        std::cout << " | "
                  << event.muonPF2PATComRelIsodBeta[event.muonIndexTight[i]];
        std::cout << " | " << event.muonPF2PATCharge[event.muonIndexTight[i]];
        std::cout << std::endl;
    }
    int numbJets{event.numJetPF2PAT > 4 ? 4 : event.numJetPF2PAT};
    std::cout << "Jets: " << event.numJetPF2PAT << std::endl;
    for (int i{0}; i < numbJets; i++)
    {
        std::cout << i;
        //    std::cout << " | " << event.jetPF2PATPt[i];
        // std::cout << " | " << std::sqrt(event.jetPF2PATPx[i] *
        // event.jetPF2PATPx[i] + event.jetPF2PATPy[i] *
        // event.jetPF2PATPy[i]);
        // std::cout << " | " << event.jetPF2PATEt[i];
        std::cout << " | " << event.jetPF2PATPtRaw[i];
        // std::cout << " | " << event.jetPF2PATUnCorPt[i];
        std::cout << " | " << event.jetPF2PATEta[i];
        std::cout << " | " << event.jetPF2PATPhi[i];
        std::cout << " | " << event.jetPF2PATNConstituents[i];
        std::cout << " | "
                  << (event.jetPF2PATNeutralHadronEnergyFraction[i] < 0.99
                      && event.jetPF2PATNeutralEmEnergyFraction[i] < 0.99)
            && ((std::abs(event.jetPF2PATEta[i]) > 2.4)
                || (event.jetPF2PATChargedEmEnergyFraction[i] < 0.99
                    && event.jetPF2PATChargedHadronEnergyFraction[i] > 0.
                    && event.jetPF2PATChargedMultiplicity[i] > 0.));
        std::cout << " | " << event.jetPF2PATdRClosestLepton[i];
        std::cout << std::endl;
    }
    std::cout << "MET: " << event.metPF2PATEt << " | " << event.metPF2PATPt
              << std::endl;

    //  std::cout << std::setprecision(1) << std::fixed;
    //  std::cerr << std::setprecision(1) << std::fixed;
}

void Cuts::dumpLooseLepInfo(const AnalysisEvent& event) const
{
    //  std::cout << std::setprecision(6) << std::fixed;
    //  std::cerr << std::setprecision(6) << std::fixed;

    std::cout << "Electrons: " << event.numElePF2PAT << std::endl;
    for (int i{0}; i < event.numElePF2PAT; i++)
    {
        TLorentzVector tempVec{event.elePF2PATPX[i],
                               event.elePF2PATPY[i],
                               event.elePF2PATPZ[i],
                               event.elePF2PATE[i]};
        std::cout << i;
        std::cout << " | " << tempVec.Pt();
        std::cout << " | " << tempVec.Eta();
        std::cout << " | " << event.elePF2PATSCEta[i];
        std::cout << " | " << event.elePF2PATComRelIsoRho[i];
        std::cout << std::endl;
    }
    std::cout << "Muons: " << event.numMuonPF2PAT << std::endl;
    for (int i{0}; i < event.numMuonPF2PAT; i++)
    {
        std::cout << i;
        std::cout << " | " << event.muonPF2PATPt[i];
        std::cout << " | " << event.muonPF2PATEta[i];
        std::cout << " | " << event.muonPF2PATComRelIsodBeta[i];
        std::cout << " | " << event.muonPF2PATGlobalID[i];
        std::cout << " | " << event.muonPF2PATTrackID[i];
        std::cout << std::endl;
    }

    //  std::cout << std::setprecision(1) << std::fixed;
    //  std::cerr << std::setprecision(1) << std::fixed;
}

double Cuts::deltaR(const float eta1,
                    const float phi1,
                    const float eta2,
                    const float phi2) const
{
    double dEta{eta1 - eta2};
    double dPhi{phi1 - phi2};
    while (std::abs(dPhi) > M_PI)
    {
        dPhi += (dPhi > 0. ? -2 * M_PI : 2 * M_PI);
    }
    //  if(singleEventInfoDump_)  std::cout << eta1 << " " << eta2 << " phi " <<
    //  phi1 << " " << phi2 << " ds: " << eta1-eta2 << " " << phi1-phi2 << " dR:
    //  " << std::sqrt((dEta*dEta)+(dPhi*dPhi)) << std::endl;
    return std::sqrt((dEta * dEta) + (dPhi * dPhi));
}

// For dumping contents of step 4. Bit more complicated than old, so doing it
// elsewhere.
void Cuts::dumpToFile(AnalysisEvent& event, const int step)
{
    std::vector<TLorentzVector> tempLepVec;
    unsigned triggerFlag[3]{};
    std::string channel{"nan"};
    std::pair<int, int> leadingLeptons[3]{}; // Initalise as empty

    // lepton ID for step0?
    event.electronIndexTight = getTightEles(event);
    event.muonIndexTight = getTightMuons(event);

    if (step == 9)
    { // Used for 2016 debug synch
        //    step9EventDump_.precision(3);
        bool synchTrigger{false};
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9 > 0)
        {
            synchTrigger = true;
        }

        //    if ( event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3 > 0
        //    ) synchTrigger = true; if (
        //    event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4 > 0 )
        //    synchTrigger = true; if (
        //    event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5 > 0 )
        //    synchTrigger = true; if (
        //    event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6 > 0 )
        //    synchTrigger = true; if (
        //    event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7 > 0 )
        //    synchTrigger = true; if (
        //    event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8 > 0 )
        //    synchTrigger = true; if (
        //    event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9 > 0 )
        //    synchTrigger = true;

        //    if ( event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1 > 0
        //    ) synchTrigger = true; if (
        //    event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3 > 0 )
        //    synchTrigger = true;

        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0)
        {
            synchTrigger = true;
        }
        if (event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0)
        {
            synchTrigger = true;
        }

        event.electronIndexTight = getTightEles(event);
        event.muonIndexTight = getTightMuons(event);

        float tempWeight{1.};
        std::pair<std::vector<int>, std::vector<float>> jetInfo;
        jetInfo = makeJetCuts(event, 0, tempWeight);
        event.jetIndex = jetInfo.first;

        step9EventDump_ << event.eventNum << " " << event.numElePF2PAT << " "
                        << event.numMuonPF2PAT << " " << event.numJetPF2PAT
                        << " " << event.electronIndexTight.size() << " "
                        << event.muonIndexTight.size() << " "
                        << event.jetIndex.size() << " " << synchTrigger
                        << std::endl;
    }

    if (step == 0)
    { // Used for 2016/2017 synch

        // Get leading 3 leptons pT
        // Search over electrons
        for (auto electronIt = event.electronIndexTight.begin();
             electronIt != event.electronIndexTight.end();
             electronIt++)
        {
            //    for ( int electronIt = 0; electronIt != event.numElePF2PAT;
            //    electronIt++ ) {
            float elePt{event.elePF2PATPT[*electronIt]};
            float itPt[3]{};
            for (unsigned j{0}; j != 3; j++)
            {
                if (leadingLeptons[j].second == 1)
                {
                    itPt[j] = event.elePF2PATPT[leadingLeptons[j].first];
                }
                else if (leadingLeptons[j].second == 2)
                {
                    itPt[j] = event.muonPF2PATPt[leadingLeptons[j].first];
                }
            }

            if (elePt > itPt[2] && elePt <= itPt[1])
            {
                leadingLeptons[2] = std::make_pair(*electronIt, 1);
            }
            else if (elePt > itPt[1] && elePt <= itPt[0])
            {
                leadingLeptons[2] = leadingLeptons[1];
                leadingLeptons[1] = std::make_pair(*electronIt, 1);
            }
            else if (elePt > itPt[0])
            {
                leadingLeptons[2] = leadingLeptons[1];
                leadingLeptons[1] = leadingLeptons[0];
                leadingLeptons[0] = std::make_pair(*electronIt, 1);
            }
        }

        // Search over muons
        for (auto muonIt = event.muonIndexTight.begin();
             muonIt != event.muonIndexTight.end();
             muonIt++)
        {
            //    for ( int  muonIt = 0; muonIt != event.numMuonPF2PAT;
            //    muonIt++ ) {
            float muonPt{event.muonPF2PATPt[*muonIt]};
            float itPt[3]{};
            for (unsigned j{0}; j != 3; j++)
            {
                if (leadingLeptons[j].second == 1)
                {
                    itPt[j] = event.elePF2PATPT[leadingLeptons[j].first];
                }
                else if (leadingLeptons[j].second == 2)
                {
                    itPt[j] = event.muonPF2PATPt[leadingLeptons[j].first];
                }
            }

            if (muonPt > itPt[2] && muonPt <= itPt[1])
            {
                leadingLeptons[2] = std::make_pair(*muonIt, 2);
            }
            else if (muonPt > itPt[1] && muonPt <= itPt[0])
            {
                leadingLeptons[2] = leadingLeptons[1];
                leadingLeptons[1] = std::make_pair(*muonIt, 2);
            }
            else if (muonPt > itPt[0])
            {
                leadingLeptons[2] = leadingLeptons[1];
                leadingLeptons[1] = leadingLeptons[0];
                leadingLeptons[0] = std::make_pair(*muonIt, 2);
            }
        }

        // Setup channel label
        int numEles{0};
        int numMuons{0};
        for (unsigned i{0}; i != 3; ++i)
        {
            if (leadingLeptons[i].second == 1)
            {
                numEles++;
            }
            if (leadingLeptons[i].second == 2)
            {
                numMuons++;
            }
        }

        if (numEles == 3 && numMuons == 0)
        {
            channel = "eee";
        }
        else if (numEles == 2 && numMuons == 1)
        {
            channel = "eem";
        }
        else if (numEles == 1 && numMuons == 2)
        {
            channel = "emm";
        }
        else if (numEles == 0 && numMuons == 3)
        {
            channel = "mmm";
        }
    }

    switch (step)
    {
        case 0:
            step0EventDump_.precision(3);
            step0EventDump_ << "|" << event.eventNum << "|" << triggerFlag[0]
                            << triggerFlag[1] << triggerFlag[2] << "|"
                            << channel << "|";
        case 2:
            step2EventDump_ << event.eventRun << " " << event.eventNum << " ";
            break;
        case 4:
            step4EventDump_ << event.eventRun << " " << event.eventNum << " ";
            break;
        case 6:
            step6EventDump_ << event.eventRun << " " << event.eventNum << " ";
            break;
    }
    for (unsigned i{0}; i < 3; i++)
    {
        switch (step)
        {
            case 0:
                step0EventDump_.precision(3);
                if (leadingLeptons[i].second == 1)
                {
                    step0EventDump_
                        << event.elePF2PATPT[leadingLeptons[i].first] << "|";
                }
                else if (leadingLeptons[i].second == 2)
                {
                    step0EventDump_
                        << event.muonPF2PATPt[leadingLeptons[i].first] << "|";
                }
                else
                {
                    step0EventDump_ << 0 << "|";
                }
                break;
            case 2:
                step2EventDump_ << tempLepVec[i].Pt() << " "
                                << tempLepVec[i].Eta() << " ";
                break;
            case 4:
                step4EventDump_ << tempLepVec[i].Pt() << " "
                                << tempLepVec[i].Eta() << " ";
                break;
            case 6:
                step6EventDump_ << tempLepVec[i].Pt() << " "
                                << tempLepVec[i].Eta() << " ";
                break;
        }
    }

    // Rel iso for step0
    for (unsigned i{0}; i < 3; i++)
    {
        switch (step)
        {
            case 0:
                step0EventDump_.precision(3);
                if (leadingLeptons[i].second == 1)
                {
                    step0EventDump_
                        << event.elePF2PATComRelIsoRho[leadingLeptons[i].first]
                        << "|";
                }
                else if (leadingLeptons[i].second == 2)
                {
                    step0EventDump_
                        << event.muonPF2PATComRelIsodBeta[leadingLeptons[i]
                                                              .first]
                        << "|";
                }
                else
                {
                    step0EventDump_ << 0 << "|";
                }
                break;
        }
    }

    // Id flag superfluous as now we are doing selection on leading leptons
    /*  switch (step) {
      case 0:
        for (unsigned i{0}; i < 3; i++){

          bool IdFlag{false};
          if ( leadingLeptons[i].second == 1 ) {
        for ( unsigned j{0}; j != event.electronIndexTight.size(); j++ ) {
          if ( event.electronIndexTight[j] == leadingLeptons[i].first ) IdFlag
      = 1;
        }
          }
          else if ( leadingLeptons[i].second == 2 ) {
        for ( unsigned j{0}; j != event.muonIndexTight.size(); j++ ) {
          if ( event.muonIndexTight[j] == leadingLeptons[i].first ) IdFlag = 1;
        }
          }
          step0EventDump_ << IdFlag << "|";
        }
        break;
      }*/

    float tempWeight{1.};
    std::pair<std::vector<int>, std::vector<float>> jetInfo;
    jetInfo = makeJetCuts(event, 0, tempWeight);
    event.jetIndex = jetInfo.first;
    event.jetSmearValue = jetInfo.second;

    switch (step)
    {
        case 0:
            step0EventDump_.precision(3);
            step0EventDump_
                << event.jetPF2PATPt[0] << "|"
                << event
                       .jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[0]
                << "|";
            break;
    }

    for (unsigned i{0}; i < 4; i++)
    {
        switch (step)
        {
            case 2:
                step2EventDump_
                    << ((i < event.jetIndex.size())
                            ? event.jetPF2PATPtRaw[event.jetIndex[i]]
                            : -666)
                    << " "
                    << ((i < event.jetIndex.size())
                            ? event.jetPF2PATEta[event.jetIndex[i]]
                            : -666)
                    << " ";
                break;
            case 4:
                step4EventDump_
                    << ((i < event.jetIndex.size())
                            ? event.jetPF2PATPtRaw[event.jetIndex[i]]
                            : -666)
                    << " "
                    << ((i < event.jetIndex.size())
                            ? event.jetPF2PATEta[event.jetIndex[i]]
                            : -666)
                    << " ";
                break;
            case 6:
                step6EventDump_
                    << ((i < event.jetIndex.size())
                            ? event.jetPF2PATPtRaw[event.jetIndex[i]]
                            : -666)
                    << " "
                    << ((i < event.jetIndex.size())
                            ? event.jetPF2PATEta[event.jetIndex[i]]
                            : -666)
                    << " ";
                break;
        }
    }

    switch (step)
    {
        case 0:
            step0EventDump_.precision(3);
            step0EventDump_ << event.metPF2PATPt << "|";
            // Synch Cut Flow stuff
            if (triggerCuts(event, tempWeight, 0))
            {
                step0EventDump_ << "1"; // Trigger Selection - step 0
            }
            else
            {
                step0EventDump_ << "0";
            }
            if ((event.electronIndexTight.size() + event.muonIndexTight.size())
                == 3)
            {
                step0EventDump_ << "1"; // 3 tight lepton selection - step 1
            }
            else
            {
                step0EventDump_ << "0";
            }
            if ((event.electronIndexTight.size() == getLooseEles(event).size())
                && (event.muonIndexTight.size() == getLooseMuons(event).size()))
            {
                step0EventDump_ << "1"; // no additional loose leptons - step 2
            }
            else
            {
                step0EventDump_ << "0";
            }
            if ((event.electronIndexTight.size() + event.muonIndexTight.size())
                < 3)
            {
                step0EventDump_
                    << "0"; // Check to ensure there are at least three leptons
                            // - otherwise memory leak occurs.
            }
            else if (std::abs(getTrileptonZCand(
                         event, event.electronIndexTight, event.muonIndexTight))
                     > invZMassCut_)
            {
                step0EventDump_ << "1"; // Z selection - step 3
            }
            else
            {
                step0EventDump_ << "0";
            }
            jetInfo = makeJetCuts(event, 0, tempWeight);
            event.jetIndex = jetInfo.first;
            event.jetSmearValue = jetInfo.second;
            if (event.jetIndex.size() >= 1)
            {
                step0EventDump_
                    << "1"; // Jet selection, at least one jet - step 4
            }
            else
            {
                step0EventDump_ << "0";
            }
            event.bTagIndex = makeBCuts(event, event.jetIndex);
            if (event.bTagIndex.size() == 1)
            {
                step0EventDump_
                    << "1"; // b-Tag selection, exactly one b-jet - step 5
            }
            else
            {
                step0EventDump_ << "0";
            }
            if (std::sqrt(
                    2 * event.metPF2PATPt * event.wLepton.Pt()
                    * (1 - std::cos(event.metPF2PATPhi - event.wLepton.Phi())))
                > mTWCutSynch_)
            {
                step0EventDump_ << "1"; // MET selection, step 6
            }
            else
            {
                step0EventDump_ << "0";
            }
            //    if ( (getTopMass(event) < TopMassCutUpper_) &&
            //    (getTopMass(event) > TopMassCutLower_) ) step0EventDump_ <<
            //    "1"; // Top Mass cut, step 7 else step0EventDump_ << "0";
            step0EventDump_ << std::endl;
            break;
        case 2:
            step2EventDump_ << event.metPF2PATPt;
            step2EventDump_ << std::endl;
            break;
        case 4:
            step4EventDump_ << event.metPF2PATPt;
            step4EventDump_ << std::endl;
            break;
        case 6:
            step6EventDump_ << event.metPF2PATPt;
            step6EventDump_ << std::endl;
            break;
    }
}

float Cuts::getLeptonWeight(const AnalysisEvent& event, const int syst) const
{
    // If number of electrons is > 1  then both z pair are electrons, so get
    // their weight
    if (!isMC_)
    {
        return 1.;
    }

    float leptonWeight{1.};

    if (numTightEle_ == 2)
    {
        leptonWeight *= eleSF(event.zPairLeptons.first.Pt(),
                              event.elePF2PATSCEta[event.zPairIndex.first],
                              syst);
        leptonWeight *= eleSF(event.zPairLeptons.second.Pt(),
                              event.elePF2PATSCEta[event.zPairIndex.second],
                              syst);
        // Double + Single trigger SFs
        /*
            // Leading lepton is always above 25 GeV, no need to do logic
           for that if ( event.zPairLeptons.second.Pt() < 20 ) leptonWeight
           *= 0.97869; else leptonWeight *= 0.98845;
        */
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

float Cuts::eleSF(const double pt, const double eta, const int syst) const
{
    double maxPt{h_eleSFs->GetYaxis()->GetXmax() - 0.1};
    double minRecoPt{h_eleReco->GetYaxis()->GetXmin() + 0.1};
    unsigned bin1{0}, bin2{0};

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

    float eleIdSF = h_eleSFs->GetBinContent(bin1);
    float eleRecoSF = h_eleReco->GetBinContent(bin2);

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

float Cuts::muonSF(const double pt, const double eta, const int syst) const
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
                {{{0.8320508711624447, 0.0012035114525204718}, {0.8841990841243943, 0.0006456623925952075}, {0.9398541581807597, 0.00019745286244582550}, {0.9715339114232815, 0.00011058622392912809}, {0.9821126684551369, 0.00019326460739818240}, {0.9883449252504323, 0.00025997334405799380}}},  // η 0-1
                {{{0.8268555228705254, 0.0021360555398695150}, {0.8773674976155514, 0.0011967279748543760}, {0.9352046568162169, 0.00037976619824492966}, {0.9693481328970910, 0.00017969930458041310}, {0.9815178116187640, 0.00037593261122859590}, {0.9891116294457432, 0.00048771655437789555}}},  // η 1-1.2
                {{{0.8793298905774019, 0.0009985477167805514}, {0.9186050959597267, 0.0005882991946232953}, {0.9566228717509230, 0.00020888314313107577}, {0.9801908185717985, 0.00006485888238673204}, {0.9881187859538099, 0.00021108637100444955}, {0.9923379543424244, 0.00029549392187439560}}},  // η 1.2-2.1111
                {{{0.9197912050897632, 0.0015156212323584166}, {0.9515318857727514, 0.0008718611746341400}, {0.9761660610840468, 0.00033218732203553470}, {0.9902672402526331, 0.00018127687682536870}, {0.9939443515919232, 0.00037720227885012716}, {0.9953757744989205, 0.00060617533879051520}}},  // η 2.1+
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

        unsigned binId1{0}, binIso1{0};
        unsigned binId2{0}, binIso2{0};

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

        float muonIdSF{1.0};
        float muonPFisoSF{1.0};
        muonIdSF = (h_muonIDs1->GetBinContent(binId1) * lumiRunsBCDEF_
                    + h_muonIDs2->GetBinContent(binId2) * lumiRunsGH_)
                   / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06);
        muonPFisoSF = (h_muonPFiso1->GetBinContent(binIso1) * lumiRunsBCDEF_
                       + h_muonPFiso2->GetBinContent(binIso2) * lumiRunsGH_)
                      / (lumiRunsBCDEF_ + lumiRunsGH_ + 1.0e-06);
        //    muonIdSF = ( h_muonIDs1->GetBinContent(binId1) );
        //    muonPFisoSF = ( h_muonPFiso1->GetBinContent(binIso1) );
        //    muonIdSF = ( h_muonIDs2->GetBinContent(binId2) );
        //    muonPFisoSF = ( h_muonPFiso2->GetBinContent(binIso2) );

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
                //      muonIdSF += ( h_muonIDs1->GetBinError(binId1) );
                //      muonPFisoSF += ( h_muonPFiso1->GetBinError(binIso1) );
                //      muonIdSF += ( h_muonIDs2->GetBinError(binId2) );
                //      muonPFisoSF += ( h_muonPFiso2->GetBinError(binIso2) );
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
                //      muonIdSF -= ( h_muonIDs1->GetBinError(binId1) );
                //      muonPFisoSF -= ( h_muonPFiso1->GetBinError(binIso1) );
                //      muonIdSF -= ( h_muonIDs2->GetBinError(binId2) );
                //      muonPFisoSF -= ( h_muonPFiso2->GetBinError(binIso2) );
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
                     "Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt",
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
        std::vector<float> tempUp;
        std::vector<float> tempDown;

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

float Cuts::getJECUncertainty(const float pt,
                              const float eta,
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

    float lowFact{syst == 4 ? jecSFUp_[etaBin][ptBin]
                            : jecSFDown_[etaBin][ptBin]};
    float hiFact{syst == 4 ? jecSFUp_[etaBin][ptBin + 1]
                           : jecSFDown_[etaBin][ptBin + 1]};
    // Now do some interpolation
    float a{(hiFact - lowFact) / (ptMaxJEC_[ptBin] - ptMinJEC_[ptBin])};
    float b{(lowFact * (ptMaxJEC_[ptBin]) - hiFact * ptMinJEC_[ptBin])
            / (ptMaxJEC_[ptBin] - ptMinJEC_[ptBin])};
    return (syst == 4 ? a * pt + b : -(a * pt + b));
}

TLorentzVector Cuts::getJetLVec(AnalysisEvent& event,
                                const int index,
                                const int syst,
                                const bool initialRun)
{
    TLorentzVector returnJet;
    float newSmearValue{1.0};

    if (!initialRun)
    {
        newSmearValue = event.jetSmearValue[index];
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);
        returnJet *= newSmearValue;
        if (isMC_)
        {
            float jerUncer{
                getJECUncertainty(returnJet.Pt(), returnJet.Eta(), syst)};
            returnJet *= 1 + jerUncer;
        }
        return returnJet;
    }

    if (!isMC_)
    {
        event.jetSmearValue.emplace_back(newSmearValue);
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);
        return returnJet;
    }

    // data is dealt with above - only MC should be dealt with below

    float jerSF{1.};
    float jerSigma{0.};
    std::pair<float, float> jetSFs{};

    if (!is2016_)
    {
        jetSFs = jet2017SFs(std::abs(event.jetPF2PATEta[index]));
    }
    else
    {
        jetSFs = jet2016SFs(std::abs(event.jetPF2PATEta[index]));
    }

    jerSF = jetSFs.first;
    jerSigma = jetSFs.second;

    if (syst == 16)
    {
        jerSF += jerSigma;
    }
    else if (syst == 32)
    {
        jerSF -= jerSigma;
    }

    double dR = deltaR(event.genJetPF2PATEta[index],
                       event.genJetPF2PATPhi[index],
                       event.jetPF2PATEta[index],
                       event.jetPF2PATPhi[index]);
    double min_dR = std::numeric_limits<double>::infinity();
    double dPt = event.jetPF2PATPtRaw[index] - event.genJetPF2PATPT[index];

    if (dR > min_dR)
    { // If dR is greater than infinity ... just return the unsmeared jet
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);
        return returnJet;
    }

    if (isMC_)
    {
        if (event.genJetPF2PATPT[index] > 1e-2 && dR < (0.4 / 2.0)
            && std::abs(dPt) < 3.0 * jerSigma * event.jetPF2PATPtRaw[index])
        { // If matching from GEN to RECO using dR<Rcone/2 and dPt < 3*sigma,
          // just scale, just scale
            newSmearValue =
                1. + (jerSF - 1.) * dPt / (event.jetPF2PATPtRaw[index]);
            returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                                 event.jetPF2PATPy[index],
                                 event.jetPF2PATPz[index],
                                 event.jetPF2PATE[index]);
            returnJet *= newSmearValue;
        }

        else
        { // If not matched to a gen jet, randomly smear
            double sigma = jerSigma * std::sqrt(jerSF * jerSF - 1.0);
            std::normal_distribution<> d(0, sigma);
            std::mt19937 gen(rand());
            newSmearValue = 1.0 + d(gen);
            returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                                 event.jetPF2PATPy[index],
                                 event.jetPF2PATPz[index],
                                 event.jetPF2PATE[index]);
            returnJet *= newSmearValue;
        }

        if (returnJet.E() < 1e-2)
        { // Negative or too small smearFactor. We would change direction of the
          // jet
            double newSmearFactor = 1e-2 / event.jetPF2PATE[index];
            newSmearValue = newSmearFactor;
            returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                                 event.jetPF2PATPy[index],
                                 event.jetPF2PATPz[index],
                                 event.jetPF2PATE[index]);
            returnJet *= newSmearValue;
        }
    }

    else
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);

    if (isMC_)
    {
        float jerUncer{
            getJECUncertainty(returnJet.Pt(), returnJet.Eta(), syst)};
        returnJet *= 1 + jerUncer;
    }

    tempSmearValue_ = newSmearValue;
    return returnJet;
}

std::pair<float, float> Cuts::jet2016SFs(const float eta) const
{
    // JER Scaling Factors and uncertainities for 2016
    float jerSF{0.};
    float jerSigma{0.};

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

std::pair<double, double> Cuts::jet2017SFs(const double eta) const
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
                      float& mcTag,
                      float& mcNoTag,
                      float& dataTag,
                      float& dataNoTag,
                      float& err1,
                      float& err2,
                      float& err3,
                      float& err4) const
{
    // Use b-tagging efficiencies and scale factors.
    // Firstly get efficiency for pt/eta bin here.
    float eff{1.};
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

    float SFerr{0.};
    float jetPt = jet.Pt();
    float maxBjetPt{670};
    float maxLjetPt{1000.0};
    bool doubleUncertainty{false};
    // Do some things if it's a b or c

    if (partonFlavour == 5)
    {
        if (jetPt > maxBjetPt)
        {
            jetPt = maxBjetPt;
            doubleUncertainty = true;
        }
        // jet_scalefactor = beautyReader.eval_auto_bounds("central",
        // BTagEntry::FLAV_B, jet.Eta(), jetPt); jet_scalefactor_up =
        // beautyReader.eval_auto_bounds("up", BTagEntry::FLAV_B, jet.Eta(),
        // jetPt); jet_scalefactor_do = beautyReader.eval_auto_bounds("down",
        // BTagEntry::FLAV_B, jet.Eta(), jetPt);
        jet_scalefactor = getBweight_backup(0, 0, jetPt);
        jet_scalefactor_up = getBweight_backup(0, 1, jetPt);
        jet_scalefactor_do = getBweight_backup(0, -1, jetPt);
    }

    else if (partonFlavour == 4)
    {
        if (jetPt > maxBjetPt)
        {
            jetPt = maxBjetPt;
            doubleUncertainty = true;
        }
        // jet_scalefactor = charmReader.eval_auto_bounds("central",
        // BTagEntry::FLAV_C, jet.Eta(), jetPt); jet_scalefactor_up =
        // charmReader.eval_auto_bounds("up", BTagEntry::FLAV_C, jet.Eta(),
        // jetPt); jet_scalefactor_do = charmReader.eval_auto_bounds("down",
        // BTagEntry::FLAV_C, jet.Eta(), jetPt);
        jet_scalefactor = getBweight_backup(1, 0, jetPt);
        jet_scalefactor_up = getBweight_backup(1, 1, jetPt);
        jet_scalefactor_do = getBweight_backup(1, -1, jetPt);
    }

    // Light jets
    else
    {
        if (jetPt > maxLjetPt)
        {
            jetPt = maxLjetPt;
            doubleUncertainty = true;
        }
        // jet_scalefactor = lightReader.eval_auto_bounds("central",
        // BTagEntry::FLAV_UDSG, jet.Eta(), jetPt); jet_scalefactor_up =
        // lightReader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG, jet.Eta(),
        // jetPt); jet_scalefactor_do = lightReader.eval_auto_bounds("down",
        // BTagEntry::FLAV_UDSG, jet.Eta(), jetPt);
        jet_scalefactor = getBweight_backup(2, 0, jetPt);
        jet_scalefactor_up = getBweight_backup(2, 1, jetPt);
        jet_scalefactor_do = getBweight_backup(2, -1, jetPt);
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

double Cuts::getBweight_backup(const int flavour,
                               const int type,
                               const double pt) const
{
    double sf{1.0};
    const double& x{pt};

    if (!is2016_)
    { // is 2017
        // MEDIUM
        switch (flavour)
        {
            case 0: // B flavour
                switch (type)
                {
                    case 0: // central
                        if (pt > 20 && pt < 1000)
                        {
                            return 0.941966
                                   * ((1. + (0.0241018 * x))
                                      / (1. + (0.0248776 * x)));
                        }
                        else
                        {
                            throw std::runtime_error(
                                "pT out of range of b tag SFs");
                        }
                    case 1: // up
                        if (pt > 20 && pt < 30)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.051529459655284882;
                        }
                        else if (pt < 50)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.017671864479780197;
                        }
                        else if (pt < 70)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.022306634113192558;
                        }
                        else if (pt < 100)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.023042259737849236;
                        }
                        else if (pt < 140)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.039661582559347153;
                        }
                        else if (pt < 200)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.061514820903539658;
                        }
                        else if (pt < 300)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.071018315851688385;
                        }
                        else if (pt < 600)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.054169680923223495;
                        }
                        else if (pt < 1000)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.063008971512317657;
                        }
                        else
                        {
                            throw std::runtime_error(
                                "pT out of range of b tag SFs");
                        }
                    case -1: // down
                        if (pt > 20 && pt < 30)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.051529459655284882;
                        }
                        else if (pt < 50)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.017671864479780197;
                        }
                        else if (pt < 70)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.022306634113192558;
                        }
                        else if (pt < 100)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.023042259737849236;
                        }
                        else if (pt < 140)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.039661582559347153;
                        }
                        else if (pt < 200)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.061514820903539658;
                        }
                        else if (pt < 300)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.071018315851688385;
                        }
                        else if (pt < 600)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.054169680923223495;
                        }
                        else if (pt < 1000)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.063008971512317657;
                        }
                        else
                        {
                            throw std::runtime_error(
                                "pT out of range of b tag SFs");
                        }
                    default:
                        throw std::runtime_error(
                            "Unknown b tag systematic type");
                }
            case 1: // C flavour
                switch (type)
                {
                    case 0: // central
                        if (pt > 20 && pt < 1000)
                        {
                            return 0.941966
                                   * ((1. + (0.0241018 * x))
                                      / (1. + (0.0248776 * x)));
                        }
                        else
                        {
                            throw std::runtime_error(
                                "pT out of range of b tag SFs");
                        }
                    case 1: // up
                        if (pt > 20 && pt < 30)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.15458837151527405;
                        }
                        else if (pt < 50)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.053015593439340591;
                        }
                        else if (pt < 70)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.066919900476932526;
                        }
                        else if (pt < 100)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.069126777350902557;
                        }
                        else if (pt < 140)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.11898474395275116;
                        }
                        else if (pt < 200)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.18454445898532867;
                        }
                        else if (pt < 300)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.21305495500564575;
                        }
                        else if (pt < 600)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.16250903904438019;
                        }
                        else if (pt < 1000)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   + 0.18902692198753357;
                        }
                        else
                        {
                            throw std::runtime_error(
                                "pT out of range of b tag SFs");
                        }
                    case -1: // down
                        if (pt > 20 && pt < 30)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.15458837151527405;
                        }
                        else if (pt < 50)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.053015593439340591;
                        }
                        else if (pt < 70)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.066919900476932526;
                        }
                        else if (pt < 100)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.069126777350902557;
                        }
                        else if (pt < 140)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.11898474395275116;
                        }
                        else if (pt < 200)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.18454445898532867;
                        }
                        else if (pt < 300)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.21305495500564575;
                        }
                        else if (pt < 600)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.16250903904438019;
                        }
                        else if (pt < 1000)
                        {
                            return (0.941966
                                    * ((1. + (0.0241018 * x))
                                       / (1. + (0.0248776 * x))))
                                   - 0.18902692198753357;
                        }
                        else
                        {
                            throw std::runtime_error(
                                "pT out of range of b tag SFs");
                        }
                    default:
                        throw std::runtime_error(
                            "Unknown b tag systematic type");
                }
            case 2: // UDSG flavour
                switch (type)
                {
                    case 0: // central
                        return 0.949449 + 0.000516201 * x + 7.13398e-08 * x * x
                               + -3.55644e-10 * x * x * x;
                    case 1: // up
                        return (0.949449 + 0.000516201 * x + 7.13398e-08 * x * x
                                + -3.55644e-10 * x * x * x)
                               * (1
                                  + (0.115123 + 0.000153114 * x
                                     + -1.72111e-07 * x * x));
                    case -1: // down
                        return (0.949449 + 0.000516201 * x + 7.13398e-08 * x * x
                                + -3.55644e-10 * x * x * x)
                               * (1
                                  - (0.115123 + 0.000153114 * x
                                     + -1.72111e-07 * x * x));
                    default:
                        throw std::runtime_error(
                            "Unknown b tag systematic type");
                }
            default:
                throw std::runtime_error("Unknown b tag systematic flavour");
        }
    }
    else
    { // is 2016
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
    }

    return sf;
}
