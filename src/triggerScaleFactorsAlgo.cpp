#include "AnalysisEvent.hpp"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TMVA/Config.h"
#include "TMVA/Timer.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom.h"
#include "TTree.h"
#include "config_parser.hpp"
#include "triggerScaleFactorsAlgo.hpp"

#include <boost/program_options.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <sys/stat.h>

// Double_t ptBins[] = { 0, 10, 15, 18, 22, 24, 26, 30, 40, 50, 60, 80, 120, 500
// }; Int_t numPt_bins = {13}; Double_t etaBins[] = { -2.4, -2.1, -1.6, -1.2,
// -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4 }; Int_t numEta_bins =
// {13};

Double_t ptBins[]{15, 20, 25, 30, 40, 120, 200};
Int_t numPt_bins{6};
Double_t etaBins[]{-2.4, -1.2, 0.0, 1.2, 2.4};
Int_t numEta_bins{4};

TriggerScaleFactors::TriggerScaleFactors()
    :

    is2016_{false}
    , zCuts_{false}
    , jetCuts_{false}
    , bCuts_{false}
    , applyHltSf_{false}
    , isPart1_{false}
    , isPart2_{false}
    , customElectronCuts_{false}
    , customMuonCuts_{false}
    ,

    // For efficiencies
    numberPassedElectrons()
    , numberTriggeredDoubleElectrons()
    ,

    numberPassedMuons()
    , numberTriggeredDoubleMuons()
    ,

    numberPassedMuonElectrons()
    , numberTriggeredMuonElectrons()
    ,

    // alpha systematics

    numberSelectedElectrons()
    , numberSelectedMuons()
    , numberSelectedMuonElectrons()
    ,

    numberSelectedDoubleElectronsTriggered()
    , numberSelectedDoubleMuonsTriggered()
    , numberSelectedMuonElectronsTriggered()
{
    //// Plots for turn on curve studies

    // MC histos

    p_electron1_pT_MC = new TProfile("electron1_pT_MC",
                                     "p_{T} turn-on curve for leading electron",
                                     numPt_bins,
                                     ptBins);
    p_electron1_eta_MC = new TProfile("electron1_eta_MC",
                                      "#eta turn-on curve for leading",
                                      numEta_bins,
                                      etaBins);
    p_electron2_pT_MC =
        new TProfile("electron2_pT_MC",
                     "p_{T} turn-on curve for subleading electron",
                     numPt_bins,
                     ptBins);
    p_electron2_eta_MC = new TProfile("electron2_eta_MC",
                                      "#eta turn-on curve for subleading",
                                      numEta_bins,
                                      etaBins);

    p_muon1_pT_MC = new TProfile("muon1_pT_MC",
                                 "p_{T} turn-on curve for leading #mu",
                                 numPt_bins,
                                 ptBins);
    p_muon1_eta_MC = new TProfile("muon1_eta_MC",
                                  "#eta turn-on curve for leading #mu",
                                  numEta_bins,
                                  etaBins);
    p_muon2_pT_MC = new TProfile("muon2_pT_MC",
                                 "p_{T} turn-on curve for subleading #mu",
                                 numPt_bins,
                                 ptBins);
    p_muon2_eta_MC = new TProfile("muon2_eta_MC",
                                  "#eta turn-on curve for subleading #mu",
                                  numEta_bins,
                                  etaBins);

    p_muonElectron1_pT_MC =
        new TProfile("muonElectron1_pT_MC",
                     "p_{T} turn-on curve for leading lep e#mu",
                     numPt_bins,
                     ptBins);
    p_muonElectron1_eta_MC =
        new TProfile("muonElectron1_eta_MC",
                     "#eta turn-on curve for leading lep e#mu",
                     numEta_bins,
                     etaBins);
    p_muonElectron2_pT_MC =
        new TProfile("muonElectron2_pT_MC",
                     "p_{T} turn-on curve for subleading lep e#mu",
                     numPt_bins,
                     ptBins);
    p_muonElectron2_eta_MC =
        new TProfile("muonElectron2_eta_MC",
                     "#eta turn-on curve for subleading lep e#mu",
                     numEta_bins,
                     etaBins);

    p_electrons_pT_MC = new TProfile2D(
        "p_electrons_pT_MC", "", numPt_bins, ptBins, numPt_bins, ptBins);
    p_electrons_eta_MC = new TProfile2D(
        "p_electrons_eta_MC", "", numEta_bins, etaBins, numEta_bins, etaBins);
    p_muons_pT_MC = new TProfile2D(
        "p_muons_pT_MC", "", numPt_bins, ptBins, numPt_bins, ptBins);
    p_muons_eta_MC = new TProfile2D(
        "p_muons_eta_MC", "", numEta_bins, etaBins, numEta_bins, etaBins);
    p_muonElectrons_pT_MC = new TProfile2D(
        "p_muonElectrons_pT_MC", "", numPt_bins, ptBins, numPt_bins, ptBins);
    p_muonElectrons_eta_MC = new TProfile2D("p_muonElectrons_eta_MC",
                                            "",
                                            numEta_bins,
                                            etaBins,
                                            numEta_bins,
                                            etaBins);

    // Data histos

    p_electron1_pT_data =
        new TProfile("electron1_pT_data",
                     "p_{T} turn-on curve for leading electron",
                     numPt_bins,
                     ptBins);
    p_electron1_eta_data =
        new TProfile("electron1_eta_data",
                     "#eta turn-on curve for leading electron",
                     numEta_bins,
                     etaBins);
    p_electron2_pT_data =
        new TProfile("electron2_pT_data",
                     "p_{T} turn-on curve for subleading electron",
                     numPt_bins,
                     ptBins);
    p_electron2_eta_data =
        new TProfile("electron2_eta_data",
                     "#eta turn-on curve for subleading electron",
                     numEta_bins,
                     etaBins);

    p_muon1_pT_data = new TProfile("muon1_pT_data",
                                   "p_{T} turn-on curve for leading #mu",
                                   numPt_bins,
                                   ptBins);
    p_muon1_eta_data = new TProfile("muon1_eta_data",
                                    "#eta turn-on curve for leading #mu",
                                    numEta_bins,
                                    etaBins);
    p_muon2_pT_data = new TProfile("muon2_pT_data",
                                   "p_{T} turn-on curve for subleading #mu",
                                   numPt_bins,
                                   ptBins);
    p_muon2_eta_data = new TProfile("muon2_eta_data",
                                    "#eta turn-on curve for subleading #mu",
                                    numEta_bins,
                                    etaBins);

    p_muonElectron1_pT_data =
        new TProfile("muonElectron1_pT_data",
                     "p_{T} turn-on curve for leading lep e#mu",
                     numPt_bins,
                     ptBins);
    p_muonElectron1_eta_data =
        new TProfile("muonElectron1_eta_data",
                     "#eta turn-on curve for leading lep e#mu",
                     numEta_bins,
                     etaBins);
    p_muonElectron2_pT_data =
        new TProfile("muonElectron2_pT_data",
                     "p_{T} turn-on curve for subleading lep e#mu",
                     numPt_bins,
                     ptBins);
    p_muonElectron2_eta_data =
        new TProfile("muonElectron2_eta_data",
                     "#eta turn-on curve for subleading lep e#mu",
                     numEta_bins,
                     etaBins);

    p_electrons_pT_data = new TProfile2D(
        "p_electrons_pT_data", "", numPt_bins, ptBins, numPt_bins, ptBins);
    p_electrons_eta_data = new TProfile2D(
        "p_electrons_eta_data", "", numEta_bins, etaBins, numEta_bins, etaBins);
    p_muons_pT_data = new TProfile2D(
        "p_muons_pT_data", "", numPt_bins, ptBins, numPt_bins, ptBins);
    p_muons_eta_data = new TProfile2D(
        "p_muons_eta_data", "", numEta_bins, etaBins, numEta_bins, etaBins);
    p_muonElectrons_pT_data = new TProfile2D(
        "p_muonElectrons_pT_data", "", numPt_bins, ptBins, numPt_bins, ptBins);
    p_muonElectrons_eta_data = new TProfile2D("p_muonElectrons_eta_data",
                                              "",
                                              numEta_bins,
                                              etaBins,
                                              numEta_bins,
                                              etaBins);

    if (is2016_)
    {
        muonHltFile1 = new TFile{
            "scaleFactors/2016/HLT_Mu24_EfficienciesAndSF_RunBtoF.root"};
        muonHltFile2 = new TFile{
            "scaleFactors/2016/HLT_Mu24_EfficienciesAndSF_RunGtoH.root"};
        muonHltFile1->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        muonHltFile2->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        h_muonHlt1 = dynamic_cast<TH2F*>(muonHltFile1->Get(
            "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"));
        h_muonHlt2 = dynamic_cast<TH2F*>(muonHltFile2->Get(
            "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"));
    }
    else
    {
        muonHltFile1 = new TFile{
            "scaleFactors/2017/HLT_Mu24_EfficienciesAndSF_RunBtoF.root"};
        muonHltFile2 = nullptr;
        muonHltFile1->cd("IsoMu27_PtEtaBins");
        h_muonHlt1 = dynamic_cast<TH2F*>(
            muonHltFile1->Get("IsoMu27_PtEtaBins/abseta_pt_ratio"));
    }
}

TriggerScaleFactors::~TriggerScaleFactors()
{
}

void TriggerScaleFactors::parseCommandLineArguements(int argc, char* argv[])
{
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()("help,h", "Print this message.")(
        "config,c",
        po::value<std::string>(&config)->required(),
        "The configuration file to be used.")(
        "part1",
        po::bool_switch(&isPart1_),
        "Only look at data and triggers from Runs B-F")(
        "part2",
        po::bool_switch(&isPart2_),
        "Only look at data and triggers from Runs B-F")(
        "hlt",
        po::bool_switch(&applyHltSf_),
        "Apply centraly provided HLT SFs onto certian plots")(
        "electronCuts",
        po::value<std::vector<double>>(&electronCutsVars),
        "Set custom electron pT and eta cuts in the format MINELE1PT MAXELE1PT "
        "MINELE2PT MAXLELE2PT MINELE1ETA MAXLELE1ETA MINELE2ETA MAXLELE2ETA. "
        "If a negative value is used for max pT, the cut is skipped. ")(
        "muonCuts",
        po::value<std::vector<double>>(&muonCutsVars),
        "Set custom muon pT and eta cuts in the format MINMU1PT MAXMU1PT "
        "MINMU2PT MAXLMU2PT MINMU1ETA MAXLMU1ETA MINMU2ETA MAXLMU2ETA. If a "
        "negative value is used for max pT, the cut is skipped. ")(
        ",n",
        po::value<long>(&nEvents)->default_value(0),
        "The number of events to be run over. All if set to 0.")(
        "outFolder,o",
        po::value<std::string>(&outFolder)
            ->default_value("plots/scaleFactors/"),
        "The output directory for the plots. Overrides the config file.")(
        "postfix,s",
        po::value<std::string>(&postfix)->default_value("default"),
        "Set postfix for plots. Overrides the config file.")(
        "2016", po::bool_switch(&is2016_), "Use 2016 conditions (SFs, et al.)")(
        "zCuts", po::bool_switch(&zCuts_), "Use z mass veto cut")(
        "jetCuts", po::bool_switch(&jetCuts_), "Use jet cuts")(
        "bCuts", po::bool_switch(&bCuts_), "Use btag cuts")(
        "nFiles,f",
        po::value<int>(&numFiles)->default_value(-1),
        "Number of files to run over. All if set to -1.");
    po::variables_map vm;

    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc;
            std::exit(0);
        }

        po::notify(vm);
    }
    catch (const po::error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << desc;
        std::exit(1);
    }

    gErrorIgnoreLevel = kInfo;

    // Set up environment a little.
    std::cout << std::setprecision(5) << std::fixed;

    if (isPart1_ || isPart2_)
    {
        DO_HIPS = true;
    }
    if (isPart1_)
    {
        HIP_ERA = true;
    }
    if (isPart2_)
    {
        HIP_ERA = false;
    }

    if (vm.count("electronCuts"))
    {
        if (electronCutsVars.size() != 8)
        {
            throw std::logic_error(
                "--electronCuts takes exactly eight arguements.");
        }
        customElectronCuts_ = true;
        std::cout << "CAUTION! Using custom electron pT and eta cuts of "
                  << std::endl;
        std::cout << "min leading electron pT: " << electronCutsVars[0]
                  << std::endl;
        std::cout << "max leading electron pT: " << electronCutsVars[1]
                  << std::endl;
        std::cout << "min subleading electron pT: " << electronCutsVars[2]
                  << std::endl;
        std::cout << "max subleading electron pT: " << electronCutsVars[3]
                  << std::endl;
        std::cout << "min leading electron eta: " << electronCutsVars[4]
                  << std::endl;
        std::cout << "max leading electron eta: " << electronCutsVars[5]
                  << std::endl;
        std::cout << "min subleading electron eta: " << electronCutsVars[6]
                  << std::endl;
        std::cout << "max subleading electron eta: " << electronCutsVars[7]
                  << std::endl;
    }
    if (vm.count("muonCuts"))
    {
        if (muonCutsVars.size() != 8)
        {
            throw std::logic_error(
                "--muonCuts takes exactly eight arguements.");
        }
        customMuonCuts_ = true;
        std::cout << "CAUTION! Using custom electron pT and eta cuts of: "
                  << std::endl;
        std::cout << "min leading muon pT: " << muonCutsVars[0] << std::endl;
        std::cout << "max leading muon pT: " << muonCutsVars[1] << std::endl;
        std::cout << "min subleading muon pT: " << muonCutsVars[2] << std::endl;
        std::cout << "max subleading muon pT: " << muonCutsVars[3] << std::endl;
        std::cout << "min leading muon eta: " << muonCutsVars[4] << std::endl;
        std::cout << "max leading muon eta: " << muonCutsVars[5] << std::endl;
        std::cout << "min subleading muon eta: " << muonCutsVars[6]
                  << std::endl;
        std::cout << "max subleading muon eta: " << muonCutsVars[7]
                  << std::endl;
    }

    // Some vectors that will be filled in the parsing.
    totalLumi = 0;
    try
    {
        Parser::parse_config(config, datasets, totalLumi);
    }
    catch (const std::exception&)
    {
        std::cerr << "ERROR Problem with a confugration file, see previous "
                     "errors for more details. If this is the only error, the "
                     "problem is with the main configuration file."
                  << std::endl;
        throw;
    }
}

void TriggerScaleFactors::runMainAnalysis()
{
    TMVA::gConfig().SetDrawProgressBar(true);

    // PU reweighting
    if (!is2016_)
    {
        // Make pileupReweighting stuff here
        dataPileupFile = new TFile("pileup/2017/truePileupTest.root", "READ");
        dataPU = dynamic_cast<TH1D*>(dataPileupFile->Get("pileup")->Clone());
        mcPileupFile = new TFile("pileup/2017/pileupMC.root", "READ");
        mcPU = dynamic_cast<TH1D*>(mcPileupFile->Get("pileup")->Clone());

        // Get systematic files too.
        systUpFile = new TFile("pileup/2017/truePileupUp.root", "READ");
        pileupUpHist = dynamic_cast<TH1D*>(systUpFile->Get("pileup")->Clone());
        systDownFile = new TFile("pileup/2017/truePileupDown.root", "READ");
        pileupDownHist =
            dynamic_cast<TH1D*>(systDownFile->Get("pileup")->Clone());
    }
    else
    {
        // Make pileupReweighting stuff here
        dataPileupFile = new TFile("pileup/2016/truePileupTest.root", "READ");
        dataPU = dynamic_cast<TH1D*>(dataPileupFile->Get("pileup")->Clone());
        mcPileupFile = new TFile("pileup/2016/pileupMC.root", "READ");
        mcPU = dynamic_cast<TH1D*>(mcPileupFile->Get("pileup")->Clone());

        // Get systematic files too.
        systUpFile = new TFile("pileup/2016/truePileupUp.root", "READ");
        pileupUpHist = dynamic_cast<TH1D*>(systUpFile->Get("pileup")->Clone());
        systDownFile = new TFile("pileup/2016/truePileupDown.root", "READ");
        pileupDownHist =
            dynamic_cast<TH1D*>(systDownFile->Get("pileup")->Clone());
    }

    puReweight = dynamic_cast<TH1D*>(dataPU->Clone());
    puReweight->Scale(1.0 / puReweight->Integral());
    mcPU->Scale(1.0 / mcPU->Integral());
    puReweight->Divide(mcPU);
    puReweight->SetDirectory(nullptr);

    /// And do the same for systematic sampl
    puSystUp = dynamic_cast<TH1D*>(pileupUpHist->Clone());
    puSystUp->Scale(1.0 / puSystUp->Integral());
    puSystUp->Divide(mcPU);
    puSystUp->SetDirectory(nullptr);
    puSystDown = dynamic_cast<TH1D*>(pileupDownHist->Clone());
    puSystDown->Scale(1.0 / puSystDown->Integral());
    puSystDown->Divide(mcPU);
    puSystDown->SetDirectory(nullptr);

    dataPileupFile->Close();
    mcPileupFile->Close();
    systUpFile->Close();
    systDownFile->Close();

    bool datasetFilled = false;

    if (totalLumi == 0.)
        totalLumi = usePreLumi;
    std::cout << "Using lumi: " << totalLumi << std::endl;
    for (auto dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    {
        const std::hash<std::string> hasher;
        srand(hasher(dataset->name()));

        datasetFilled = false;
        TChain* datasetChain = new TChain(dataset->treeName().c_str());

        std::cerr << "Processing dataset " << dataset->name() << std::endl;
        if (!datasetFilled)
        {
            if (!dataset->fillChain(datasetChain, numFiles))
            {
                std::cerr << "There was a problem constructing the chain for "
                          << dataset->name()
                          << ". Continuing with next dataset.\n";
                continue;
            }
            datasetFilled = true;
        }

        else
        {
            std::string inputPostfix{};
            inputPostfix += postfix;
            std::cout << "/scratch/data/TopPhysics/miniSkims2017/"
                             + dataset->name() + inputPostfix + "SmallSkim.root"
                      << std::endl;
            datasetChain->Add(("/scratch/data/TopPhysics/miniSkims2017/"
                               + dataset->name() + inputPostfix
                               + "SmallSkim.root")
                                  .c_str());
            std::ifstream secondTree{"/scratch/data/TopPhysics/miniSkims2017/"
                                     + dataset->name() + inputPostfix
                                     + "SmallSkim1.root"};
            if (secondTree.good())
            {
                datasetChain->Add(("/scratch/data/TopPhysics/miniSkims2017/"
                                   + dataset->name() + inputPostfix
                                   + "SmallSkim1.root")
                                      .c_str());
            }
            std::ifstream thirdTree{"/scratch/data/TopPhysics/miniSkims2017/"
                                    + dataset->name() + inputPostfix
                                    + "SmallSkim2.root"};
            if (thirdTree.good())
            {
                datasetChain->Add(("/scratch/data/TopPhysics/miniSkims2017/"
                                   + dataset->name() + inputPostfix
                                   + "SmallSkim2.root")
                                      .c_str());
            }
        }

        std::cout << "Trigger flag: " << dataset->getTriggerFlag() << std::endl;

        AnalysisEvent event{dataset->isMC(), datasetChain, is2016_};

        int numberOfEvents = datasetChain->GetEntries();
        if (nEvents && nEvents < numberOfEvents)
        {
            numberOfEvents = nEvents;
        }
        auto lEventTimer =
            new TMVA::Timer(numberOfEvents, "Running over dataset ...", false);
        lEventTimer->DrawProgressBar(0, "");
        for (int i = 0; i < numberOfEvents; i++)
        {
            lEventTimer->DrawProgressBar(i);
            event.GetEntry(i);

            double pileupWeight = puReweight->GetBinContent(
                puReweight->GetXaxis()->FindBin(event.numVert));
            double eventWeight = 1.0;
            if (dataset->isMC())
            {
                eventWeight *= pileupWeight;
            }

            //      std::cout << __LINE__ << " : " << __FILE__ << std::endl;

            if (is2016_)
            {
                if (HIP_ERA && event.eventRun >= 278820 && !(dataset->isMC())
                    && DO_HIPS)
                {
                    continue;
                }
                if (!HIP_ERA && event.eventRun < 278820 && !(dataset->isMC())
                    && DO_HIPS)
                {
                    continue;
                }
            }
            if (!metFilters(event, dataset->isMC()))
            {
                continue;
            }

            // If checking impact of jet and bjet cuts add this bool ...
            bool passJetSelection = true;
            if (jetCuts_)
            {
                passJetSelection = makeJetCuts(event, (dataset->isMC()));
            }

            // Does this event pass tight electron cut?
            // Create electron index
            event.electronIndexTight = getTightElectrons(event);
            bool passDoubleElectronSelection(passDileptonSelection(event, 2)
                                             && passJetSelection);
            // Does this event pass tight muon cut?
            // Create muon index
            event.muonIndexTight = getTightMuons(event);
            bool passDoubleMuonSelection(passDileptonSelection(event, 0)
                                         && passJetSelection);

            bool passMuonElectronSelection(passDileptonSelection(event, 1)
                                           && passJetSelection);

            // Triggering stuff
            int triggerDoubleEG(0), triggerDoubleMuon(0),
                triggerMuonElectron(0); // Passes Double Lepton Trigger
            int triggerMetDoubleEG(0), triggerMetDoubleMuon(0),
                triggerMetMuonElectron(
                    0); // Passes Double Lepton and MET triggers

            // Passes event selection and MET triggers
            int triggerMetElectronSelection(0), triggerMetMuonSelection(0),
                triggerMetMuonElectronSelection(
                    0); // Passes lepton selection and MET triggers

            // Does event pass Single/Double EG trigger and the electron
            // selection?
            if (passDoubleElectronSelection)
            {
                triggerDoubleEG = event.eTrig() && event.eeTrig();
                triggerMetDoubleEG = triggerDoubleEG && metTriggerCut(event);
            }
            // Does event pass Single/Double Muon trigger and the muon
            // selection?
            if (passDoubleMuonSelection)
            {
                triggerDoubleMuon = event.muTrig() || event.mumuTrig();
                triggerMetDoubleMuon =
                    triggerDoubleMuon && metTriggerCut(event);
            }
            // Does event pass Single Electron/Single Muon/MuonEG trigger and
            // the muon selection?
            if (passMuonElectronSelection)
            {
                triggerMuonElectron = event.muEGTrig();
                triggerMetMuonElectron =
                    triggerMuonElectron && metTriggerCut(event);
            }
            //
            // Does event pass either double lepton seletion and the MET
            // triggers?
            if (passDoubleElectronSelection)
            {
                triggerMetElectronSelection = (metTriggerCut(event));
            }
            if (passDoubleMuonSelection)
            {
                triggerMetMuonSelection = (metTriggerCut(event));
            }
            if (passMuonElectronSelection)
            {
                triggerMetMuonElectronSelection = (metTriggerCut(event));
            }

            if (dataset->isMC())
            { // If is MC
                // SFs bit

                double SF = 1.0;
                if (applyHltSf_)
                {
                    double maxSfPt = h_muonHlt1->GetYaxis()->GetXmax() - 0.1;
                    double minSfPt = h_muonHlt1->GetYaxis()->GetXmin() + 0.1;
                    unsigned binSf1{0};

                    double pt = event.zPairLeptons.first.Pt();
                    double eta = event.zPairLeptons.first.Eta();

                    if (pt > maxSfPt)
                    {
                        binSf1 = h_muonHlt1->FindBin(std::abs(eta), maxSfPt);
                    }
                    else if (pt < minSfPt)
                    {
                        binSf1 = h_muonHlt1->FindBin(std::abs(eta), minSfPt);
                    }
                    else
                    {
                        binSf1 = h_muonHlt1->FindBin(std::abs(eta), pt);
                    }

                    if (is2016_)
                    {
                        unsigned binSf2{0};
                        if (pt > maxSfPt)
                        {
                            binSf2 =
                                h_muonHlt2->FindBin(std::abs(eta), maxSfPt);
                        }
                        else if (pt < minSfPt)
                        {
                            binSf2 =
                                h_muonHlt2->FindBin(std::abs(eta), minSfPt);
                        }
                        else
                        {
                            binSf2 = h_muonHlt2->FindBin(std::abs(eta), pt);
                        }

                        if (!DO_HIPS)
                        {
                            SF = (h_muonHlt1->GetBinContent(binSf1) * 19648.534
                                  + h_muonHlt2->GetBinContent(binSf2)
                                        * 16144.444)
                                 / (19648.534 + 16144.444 + 1.0e-06);
                        }
                        if (DO_HIPS && HIP_ERA)
                        {
                            SF = h_muonHlt1->GetBinContent(binSf1);
                        }
                        if (DO_HIPS && !HIP_ERA)
                        {
                            SF = h_muonHlt2->GetBinContent(binSf2);
                        }
                    }
                }

                numberPassedElectrons[0] +=
                    triggerMetElectronSelection
                    * eventWeight; // Number of electrons passing the cross
                                   // trigger and electron selection
                numberTriggeredDoubleElectrons[0] +=
                    triggerMetDoubleEG
                    * eventWeight; // Number of electrons passing both cross
                                   // trigger+electron selection AND double EG
                                   // trigger
                numberPassedMuons[0] +=
                    triggerMetMuonSelection
                    * eventWeight; // Number of muons passing the cross trigger
                                   // and muon selection
                numberTriggeredDoubleMuons[0] +=
                    triggerMetDoubleMuon * eventWeight
                    * SF; // Number of muons passing both cross trigger+muon
                          // selection AND double muon trigger
                numberPassedMuonElectrons[0] +=
                    triggerMetMuonElectronSelection
                    * eventWeight; // Number of muonEGs passing cross trigger
                                   // and muonEG selection
                numberTriggeredMuonElectrons[0] +=
                    triggerMetMuonElectron
                    * eventWeight; // Number muonEGs passing both cross
                                   // trigger+muonEG selection AND muonEG
                                   // trigger

                // Systematic stuff
                numberSelectedElectrons[0] +=
                    passDoubleElectronSelection * eventWeight;
                numberSelectedMuons[0] += passDoubleMuonSelection * eventWeight;
                numberSelectedMuonElectrons[0] +=
                    passMuonElectronSelection * eventWeight;

                numberSelectedDoubleElectronsTriggered[0] +=
                    triggerDoubleEG * eventWeight;
                ;
                numberSelectedDoubleMuonsTriggered[0] +=
                    triggerDoubleMuon * eventWeight * SF;
                numberSelectedMuonElectronsTriggered[0] +=
                    triggerMuonElectron * eventWeight;

                // Histos bit
                if (triggerMetElectronSelection > 0)
                { // If passed event selection, then will want to add to
                  // denominator
                    p_electron1_pT_MC->Fill(event.zPairLeptons.first.Pt(),
                                            triggerMetDoubleEG
                                                / triggerMetElectronSelection);
                    p_electron1_eta_MC->Fill(event.zPairLeptons.first.Eta(),
                                             triggerMetDoubleEG
                                                 / triggerMetElectronSelection);
                    p_electron2_pT_MC->Fill(event.zPairLeptons.second.Pt(),
                                            triggerMetDoubleEG
                                                / triggerMetElectronSelection);
                    p_electron2_eta_MC->Fill(event.zPairLeptons.second.Eta(),
                                             triggerMetDoubleEG
                                                 / triggerMetElectronSelection);

                    p_electrons_pT_MC->Fill(event.zPairLeptons.first.Pt(),
                                            event.zPairLeptons.second.Pt(),
                                            triggerMetDoubleEG
                                                / triggerMetElectronSelection);
                    p_electrons_eta_MC->Fill(event.zPairLeptons.first.Eta(),
                                             event.zPairLeptons.second.Eta(),
                                             triggerMetDoubleEG
                                                 / triggerMetElectronSelection);
                }
                if (triggerMetMuonSelection > 0)
                { // If passed event selection, then will want to add to
                  // denominator
                    p_muon1_pT_MC->Fill(event.zPairLeptons.first.Pt(),
                                        triggerMetDoubleMuon * SF
                                            / triggerMetMuonSelection);
                    if (event.zPairLeptons.first.Pt() > 30.)
                    {
                        p_muon1_eta_MC->Fill(event.zPairLeptons.first.Eta(),
                                             triggerMetDoubleMuon * SF
                                                 / triggerMetMuonSelection);
                    }
                    p_muon2_pT_MC->Fill(event.zPairLeptons.second.Pt(),
                                        triggerMetDoubleMuon * SF
                                            / triggerMetMuonSelection);
                    if (event.zPairLeptons.second.Pt() > 30.)
                    {
                        p_muon2_eta_MC->Fill(event.zPairLeptons.second.Eta(),
                                             triggerMetDoubleMuon * SF
                                                 / triggerMetMuonSelection);
                    }

                    p_muons_pT_MC->Fill(event.zPairLeptons.first.Pt(),
                                        event.zPairLeptons.second.Pt(),
                                        triggerMetDoubleMuon * SF
                                            / triggerMetMuonSelection);
                    p_muons_eta_MC->Fill(event.zPairLeptons.first.Eta(),
                                         event.zPairLeptons.second.Eta(),
                                         triggerMetDoubleMuon * SF
                                             / triggerMetMuonSelection);
                }
                if (triggerMetMuonElectronSelection > 0)
                { // If passed event selection, then will want to add to
                  // denominator
                    p_muonElectron1_pT_MC->Fill(
                        event.zPairLeptons.first.Pt(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                    p_muonElectron1_eta_MC->Fill(
                        event.zPairLeptons.first.Eta(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                    p_muonElectron2_pT_MC->Fill(
                        event.zPairLeptons.second.Pt(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                    p_muonElectron2_eta_MC->Fill(
                        event.zPairLeptons.second.Eta(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);

                    p_muonElectrons_pT_MC->Fill(
                        event.zPairLeptons.first.Pt(),
                        event.zPairLeptons.second.Pt(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                    p_muonElectrons_eta_MC->Fill(
                        event.zPairLeptons.first.Eta(),
                        event.zPairLeptons.second.Eta(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                }
            }
            else
            { // Else is data
                // SFs bit
                numberPassedElectrons[1] +=
                    triggerMetElectronSelection
                    * eventWeight; // Number of electrons passing the cross
                                   // trigger and electron selection
                numberTriggeredDoubleElectrons[1] +=
                    triggerMetDoubleEG
                    * eventWeight; // Number of electrons passing both cross
                                   // trigger+electron selection AND double EG
                                   // trigger
                numberPassedMuons[1] +=
                    triggerMetMuonSelection
                    * eventWeight; // Number of muons passing the cross trigger
                                   // and muon selection
                numberTriggeredDoubleMuons[1] +=
                    triggerMetDoubleMuon
                    * eventWeight; // Number of muons passing both cross
                                   // trigger+muon selection AND double muon
                                   // trigger
                numberPassedMuonElectrons[1] +=
                    triggerMetMuonElectronSelection
                    * eventWeight; // Number of muonEGs passing cross trigger
                                   // and muonEG selection
                numberTriggeredMuonElectrons[1] +=
                    triggerMetMuonElectron
                    * eventWeight; // Number muonEGs passing both cross
                                   // trigger+muonEG selection AND muonEG
                                   // trigger

                // NB No systematic stuff required for data

                // Histos bit
                if (triggerMetElectronSelection > 0)
                { // If passed event selection, then will want to add to
                  // denominator
                    p_electron1_pT_data->Fill(
                        event.zPairLeptons.first.Pt(),
                        triggerMetDoubleEG / triggerMetElectronSelection);
                    p_electron1_eta_data->Fill(
                        event.zPairLeptons.first.Eta(),
                        triggerMetDoubleEG / triggerMetElectronSelection);
                    p_electron2_pT_data->Fill(
                        event.zPairLeptons.second.Pt(),
                        triggerMetDoubleEG / triggerMetElectronSelection);
                    p_electron2_eta_data->Fill(
                        event.zPairLeptons.second.Eta(),
                        triggerMetDoubleEG / triggerMetElectronSelection);

                    p_electrons_pT_data->Fill(
                        event.zPairLeptons.first.Pt(),
                        event.zPairLeptons.second.Pt(),
                        triggerMetDoubleEG / triggerMetElectronSelection);
                    p_electrons_eta_data->Fill(
                        event.zPairLeptons.first.Eta(),
                        event.zPairLeptons.second.Eta(),
                        triggerMetDoubleEG / triggerMetElectronSelection);
                }
                if (triggerMetMuonSelection > 0)
                { // If passed event selection, then will want to add to
                  // denominator
                    p_muon1_pT_data->Fill(event.zPairLeptons.first.Pt(),
                                          triggerMetDoubleMuon
                                              / triggerMetMuonSelection);
                    if (event.zPairLeptons.first.Pt() > 30.)
                    {
                        p_muon1_eta_data->Fill(event.zPairLeptons.first.Eta(),
                                               triggerMetDoubleMuon
                                                   / triggerMetMuonSelection);
                    }
                    p_muon2_pT_data->Fill(event.zPairLeptons.second.Pt(),
                                          triggerMetDoubleMuon
                                              / triggerMetMuonSelection);
                    if (event.zPairLeptons.second.Pt() > 30.)
                    {
                        p_muon2_eta_data->Fill(event.zPairLeptons.second.Eta(),
                                               triggerMetDoubleMuon
                                                   / triggerMetMuonSelection);
                    }

                    p_muons_pT_data->Fill(event.zPairLeptons.first.Pt(),
                                          event.zPairLeptons.second.Pt(),
                                          triggerMetDoubleMuon
                                              / triggerMetMuonSelection);
                    p_muons_eta_data->Fill(event.zPairLeptons.first.Eta(),
                                           event.zPairLeptons.second.Eta(),
                                           triggerMetDoubleMuon
                                               / triggerMetMuonSelection);
                }
                if (triggerMetMuonElectronSelection > 0)
                { // If passed event selection, then will want to add to
                  // denominator
                    p_muonElectron1_pT_data->Fill(
                        event.zPairLeptons.first.Pt(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                    p_muonElectron1_eta_data->Fill(
                        event.zPairLeptons.first.Eta(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                    p_muonElectron2_pT_data->Fill(
                        event.zPairLeptons.second.Pt(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                    p_muonElectron2_eta_data->Fill(
                        event.zPairLeptons.second.Eta(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);

                    p_muonElectrons_pT_data->Fill(
                        event.zPairLeptons.first.Pt(),
                        event.zPairLeptons.second.Pt(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                    p_muonElectrons_eta_data->Fill(
                        event.zPairLeptons.first.Eta(),
                        event.zPairLeptons.second.Eta(),
                        triggerMetMuonElectron
                            / triggerMetMuonElectronSelection);
                }
            }
        }

        delete datasetChain;
    } // end dataset loop
}

std::vector<int>
    TriggerScaleFactors::getTightElectrons(const AnalysisEvent& event) const
{
    std::vector<int> electrons;

    for (int i{0}; i < event.numElePF2PAT; i++)
    {
        if (!event.elePF2PATIsGsf[i])
            continue;

        if (electrons.size() < 1 && event.elePF2PATPT[i] <= 38)
            continue;
        else if (electrons.size() >= 1 && event.elePF2PATPT[i] <= 15)
            continue;

        if (std::abs(event.elePF2PATSCEta[i]) > 2.5)
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

std::vector<int>
    TriggerScaleFactors::getTightMuons(const AnalysisEvent& event) const
{
    std::vector<int> muons;
    if (is2016_)
    {
        for (int i{0}; i < event.numMuonPF2PAT; i++)
        {
            if (!event.muonPF2PATIsPFMuon[i])
                continue;

            if (muons.size() < 1 && event.muonPF2PATPt[i] <= 26)
                continue;
            else if (muons.size() >= 1 && event.muonPF2PATPt[i] <= 20)
                continue;

            if (std::abs(event.muonPF2PATEta[i]) >= 2.4)
                continue;
            if (event.muonPF2PATComRelIsodBeta[i] >= 0.15)
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
                && std::abs(event.muonPF2PATEta[i]) <= 2.4)
            {
                if (event.muonPF2PATPt[i] >= (muons.empty() ? 26 : 29))
                {
                    muons.emplace_back(i);
                }
            }
        }
    }
    return muons;
}

bool TriggerScaleFactors::passDileptonSelection(AnalysisEvent& event,
                                                const int nElectrons) const
{
    // Check if there are at least two electrons first. Otherwise use muons.

    float invMass(0.0);
    float pT(0.0);

    // DoubleEG

    if (nElectrons == 2)
    {
        std::vector<int> leptons = event.electronIndexTight;
        for (unsigned i = 0; i < leptons.size(); i++)
        {
            for (unsigned j = i + 1; j < leptons.size(); j++)
            {
                if (event.elePF2PATCharge[leptons[i]]
                        * event.elePF2PATCharge[leptons[j]]
                    >= 0)
                    continue; // check electron pair have correct charge.
                TLorentzVector lepton1 =
                    TLorentzVector(event.elePF2PATPX[leptons[i]],
                                   event.elePF2PATPY[leptons[i]],
                                   event.elePF2PATPZ[leptons[i]],
                                   event.elePF2PATE[leptons[i]]);
                TLorentzVector lepton2 =
                    TLorentzVector(event.elePF2PATPX[leptons[j]],
                                   event.elePF2PATPY[leptons[j]],
                                   event.elePF2PATPZ[leptons[j]],
                                   event.elePF2PATE[leptons[j]]);
                float candidateMass = (lepton1 + lepton2).M();
                if (std::abs((lepton1 + lepton2).Pt()) > std::abs(pT))
                {
                    event.zPairLeptons.first =
                        lepton1.Pt() > lepton2.Pt() ? lepton1 : lepton2;
                    event.zPairIndex.first =
                        lepton1.Pt() > lepton2.Pt() ? leptons[i] : leptons[j];
                    event.zPairLeptons.second =
                        lepton1.Pt() > lepton2.Pt() ? lepton2 : lepton1;
                    event.zPairIndex.second =
                        lepton1.Pt() > lepton2.Pt() ? leptons[j] : leptons[i];
                    invMass = candidateMass;
                    pT = (lepton1 + lepton2).Pt();
                }
            }
        }
    }

    // DoubleMuon
    else if (nElectrons == 0)
    {
        std::vector<int> leptons = event.muonIndexTight;
        for (unsigned i = 0; i < leptons.size(); i++)
        {
            for (unsigned j = i + 1; j < leptons.size(); j++)
            {
                if (event.muonPF2PATCharge[leptons[i]]
                        * event.muonPF2PATCharge[leptons[j]]
                    >= 0)
                {
                    continue;
                }
                TLorentzVector lepton1 =
                    TLorentzVector(event.muonPF2PATPX[leptons[i]],
                                   event.muonPF2PATPY[leptons[i]],
                                   event.muonPF2PATPZ[leptons[i]],
                                   event.muonPF2PATE[leptons[i]]);
                TLorentzVector lepton2 =
                    TLorentzVector(event.muonPF2PATPX[leptons[j]],
                                   event.muonPF2PATPY[leptons[j]],
                                   event.muonPF2PATPZ[leptons[j]],
                                   event.muonPF2PATE[leptons[j]]);
                float candidateMass = (lepton1 + lepton2).M();
                if (std::abs((lepton1 + lepton2).Pt()) > std::abs(pT))
                {
                    event.zPairLeptons.first =
                        lepton1.Pt() > lepton2.Pt() ? lepton1 : lepton2;
                    event.zPairIndex.first =
                        lepton1.Pt() > lepton2.Pt() ? leptons[i] : leptons[j];
                    event.zPairLeptons.second =
                        lepton1.Pt() > lepton2.Pt() ? lepton2 : lepton1;
                    event.zPairIndex.second =
                        lepton1.Pt() > lepton2.Pt() ? leptons[j] : leptons[i];
                    invMass = candidateMass;
                    pT = (lepton1 + lepton2).Pt();
                }
            }
        }
    }

    // MuonEG
    else if (nElectrons == 1)
    {
        std::vector<int> electrons = event.electronIndexTight;
        std::vector<int> muons = event.muonIndexTight;
        if (electrons.empty() || muons.empty())
        {
            return false;
        }
        for (unsigned i = 0; i < electrons.size(); i++)
        {
            for (unsigned j = 0; j < muons.size(); j++)
            {
                if (!(event.elePF2PATCharge[electrons[i]]
                          * event.muonPF2PATCharge[muons[j]]
                      >= 0))
                {
                    continue; // check muon-electron pair have correct (same)
                    // charge.
                }
                TLorentzVector lepton1 =
                    TLorentzVector(event.elePF2PATPX[electrons[i]],
                                   event.elePF2PATPY[electrons[i]],
                                   event.elePF2PATPZ[electrons[i]],
                                   event.elePF2PATE[electrons[i]]);
                TLorentzVector lepton2 =
                    TLorentzVector(event.muonPF2PATPX[muons[j]],
                                   event.muonPF2PATPY[muons[j]],
                                   event.muonPF2PATPZ[muons[j]],
                                   event.muonPF2PATE[muons[j]]);
                float candidateMass = (lepton1 + lepton2).M();
                if (std::abs((lepton1 + lepton2).Pt()) > std::abs(pT))
                {
                    event.zPairLeptons.first =
                        lepton1.Pt() > lepton2.Pt() ? lepton1 : lepton2;
                    event.zPairIndex.first =
                        lepton1.Pt() > lepton2.Pt() ? electrons[i] : muons[j];
                    event.zPairLeptons.second =
                        lepton1.Pt() > lepton2.Pt() ? lepton2 : lepton1;
                    event.zPairIndex.second =
                        lepton1.Pt() > lepton2.Pt() ? muons[j] : electrons[i];
                    invMass = candidateMass;
                    pT = (lepton1 + lepton2).Pt();
                }
            }
        }
    }

    else
    {
        std::cout << "Only dilepton searches currently supported. Exiting ..."
                  << std::endl;
        exit(888);
    }

    if (!zCuts_ && invMass > 0.0)
    {
        return true;
    }
    if (zCuts_)
    {
        double zMass = std::abs(invMass - 91.1);
        if (zMass <= 20.)
        {
            return true;
        }
    }

    return false;
}

bool TriggerScaleFactors::metTriggerCut(const AnalysisEvent& event) const
{
    // clang-format off
    return is2016_
               ? event.HLT_MET200_v1 > 0 || event.HLT_MET200_v2 > 0
                     || event.HLT_MET200_v3 > 0 || event.HLT_MET200_v4 > 0
                     || event.HLT_MET200_v5 > 0 || event.HLT_MET250_v1 > 0
                     || event.HLT_MET250_v2 > 0 || event.HLT_MET250_v3 > 0
                     || event.HLT_MET250_v4 > 0 || event.HLT_MET250_v5 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v2 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v3 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v4 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v5 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v6 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v7 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v8 > 0
                     || event.HLT_PFMET170_HBHECleaned_v2 > 0
                     || event.HLT_PFMET170_HBHECleaned_v3 > 0
                     || event.HLT_PFMET170_HBHECleaned_v4 > 0
                     || event.HLT_PFMET170_HBHECleaned_v5 > 0
                     || event.HLT_PFMET170_HBHECleaned_v6 > 0
                     || event.HLT_PFMET170_HBHECleaned_v7 > 0
                     || event.HLT_PFMET170_HBHECleaned_v8 > 0
                     || event.HLT_PFMET170_HBHECleaned_v9 > 0
                     || event.HLT_PFHT300_PFMET100_v1 > 0
                     || event.HLT_PFHT300_PFMET100_v2 > 0
                     || event.HLT_PFHT300_PFMET100_v3 > 0
                     || event.HLT_PFHT300_PFMET100_v4 > 0
                     || event.HLT_PFHT300_PFMET110_v4 > 0
                     || event.HLT_PFHT300_PFMET110_v5 > 0
                     || event.HLT_PFHT300_PFMET110_v6 > 0
               : event.HLT_MET105_IsoTrk50_v1 > 0
                     || event.HLT_MET105_IsoTrk50_v3 > 0
                     || event.HLT_MET105_IsoTrk50_v4 > 0
                     || event.HLT_MET105_IsoTrk50_v5 > 0
                     || event.HLT_MET105_IsoTrk50_v6 > 0
                     || event.HLT_MET105_IsoTrk50_v7 > 0
                     || event.HLT_MET105_IsoTrk50_v8 > 0
                     || event.HLT_MET120_IsoTrk50_v1 > 0
                     || event.HLT_MET120_IsoTrk50_v3 > 0
                     || event.HLT_MET120_IsoTrk50_v4 > 0
                     || event.HLT_MET120_IsoTrk50_v5 > 0
                     || event.HLT_MET120_IsoTrk50_v6 > 0
                     || event.HLT_MET120_IsoTrk50_v7 > 0
                     || event.HLT_MET120_IsoTrk50_v8 > 0
                     || event.HLT_HT430_DisplacedDijet40_DisplacedTrack_v10 > 0
                     || event.HLT_HT430_DisplacedDijet40_DisplacedTrack_v11 > 0
                     || event.HLT_HT430_DisplacedDijet40_DisplacedTrack_v5 > 0
                     || event.HLT_HT430_DisplacedDijet40_DisplacedTrack_v6 > 0
                     || event.HLT_HT430_DisplacedDijet40_DisplacedTrack_v8 > 0
                     || event.HLT_HT430_DisplacedDijet40_DisplacedTrack_v9 > 0
                     || event.HLT_HT430_DisplacedDijet60_DisplacedTrack_v10 > 0
                     || event.HLT_HT430_DisplacedDijet60_DisplacedTrack_v11 > 0
                     || event.HLT_HT430_DisplacedDijet60_DisplacedTrack_v5 > 0
                     || event.HLT_HT430_DisplacedDijet60_DisplacedTrack_v6 > 0
                     || event.HLT_HT430_DisplacedDijet60_DisplacedTrack_v8 > 0
                     || event.HLT_HT430_DisplacedDijet60_DisplacedTrack_v9 > 0
                     || event.HLT_HT430_DisplacedDijet80_DisplacedTrack_v10 > 0
                     || event.HLT_HT430_DisplacedDijet80_DisplacedTrack_v11 > 0
                     || event.HLT_HT430_DisplacedDijet80_DisplacedTrack_v5 > 0
                     || event.HLT_HT430_DisplacedDijet80_DisplacedTrack_v6 > 0
                     || event.HLT_HT430_DisplacedDijet80_DisplacedTrack_v8 > 0
                     || event.HLT_HT430_DisplacedDijet80_DisplacedTrack_v9 > 0
                     || event.HLT_HT650_DisplacedDijet60_Inclusive_v10 > 0
                     || event.HLT_HT650_DisplacedDijet60_Inclusive_v11 > 0
                     || event.HLT_HT650_DisplacedDijet60_Inclusive_v5 > 0
                     || event.HLT_HT650_DisplacedDijet60_Inclusive_v6 > 0
                     || event.HLT_HT650_DisplacedDijet60_Inclusive_v8 > 0
                     || event.HLT_HT650_DisplacedDijet60_Inclusive_v9 > 0
                     || event.HLT_HT650_DisplacedDijet80_Inclusive_v10 > 0
                     || event.HLT_HT650_DisplacedDijet80_Inclusive_v11 > 0
                     || event.HLT_HT650_DisplacedDijet80_Inclusive_v12 > 0
                     || event.HLT_HT650_DisplacedDijet80_Inclusive_v6 > 0
                     || event.HLT_HT650_DisplacedDijet80_Inclusive_v7 > 0
                     || event.HLT_HT650_DisplacedDijet80_Inclusive_v9 > 0
                     || event.HLT_HT750_DisplacedDijet80_Inclusive_v10 > 0
                     || event.HLT_HT750_DisplacedDijet80_Inclusive_v11 > 0
                     || event.HLT_HT750_DisplacedDijet80_Inclusive_v12 > 0
                     || event.HLT_HT750_DisplacedDijet80_Inclusive_v6 > 0
                     || event.HLT_HT750_DisplacedDijet80_Inclusive_v7 > 0
                     || event.HLT_HT750_DisplacedDijet80_Inclusive_v9 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v1 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v2 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF_v10 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v1 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v2 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_PFHT60_v2 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_PFHT60_v3 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_PFHT60_v4 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_PFHT60_v5 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_PFHT60_v6 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_PFHT60_v7 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v11 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v13 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v14 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v15 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v16 > 0
                     || event.HLT_PFMET120_PFMHT120_IDTight_v9 > 0
                     || event.HLT_PFMET130_PFMHT130_IDTight_v11 > 0
                     || event.HLT_PFMET130_PFMHT130_IDTight_v13 > 0
                     || event.HLT_PFMET130_PFMHT130_IDTight_v14 > 0
                     || event.HLT_PFMET130_PFMHT130_IDTight_v15 > 0
                     || event.HLT_PFMET130_PFMHT130_IDTight_v16 > 0
                     || event.HLT_PFMET130_PFMHT130_IDTight_v9 > 0
                     || event.HLT_PFMET140_PFMHT140_IDTight_v11 > 0
                     || event.HLT_PFMET140_PFMHT140_IDTight_v13 > 0
                     || event.HLT_PFMET140_PFMHT140_IDTight_v14 > 0
                     || event.HLT_PFMET140_PFMHT140_IDTight_v15 > 0
                     || event.HLT_PFMET140_PFMHT140_IDTight_v16 > 0
                     || event.HLT_PFMET140_PFMHT140_IDTight_v17 > 0
                     || event.HLT_PFMET140_PFMHT140_IDTight_v18 > 0
                     || event.HLT_PFMET140_PFMHT140_IDTight_v9 > 0
                     || event.HLT_PFMET200_HBHE_BeamHaloCleaned_v5 > 0
                     || event.HLT_PFMET200_HBHE_BeamHaloCleaned_v6 > 0
                     || event.HLT_PFMET200_HBHE_BeamHaloCleaned_v7 > 0
                     || event.HLT_PFMET250_HBHECleaned_v2 > 0
                     || event.HLT_PFMET250_HBHECleaned_v3 > 0
                     || event.HLT_PFMET250_HBHECleaned_v4 > 0
                     || event.HLT_PFMET250_HBHECleaned_v5 > 0
                     || event.HLT_PFMET250_HBHECleaned_v6 > 0
                     || event.HLT_PFMET250_HBHECleaned_v7 > 0
                     || event.HLT_PFMET300_HBHECleaned_v2 > 0
                     || event.HLT_PFMET300_HBHECleaned_v3 > 0
                     || event.HLT_PFMET300_HBHECleaned_v4 > 0
                     || event.HLT_PFMET300_HBHECleaned_v5 > 0
                     || event.HLT_PFMET300_HBHECleaned_v6 > 0
                     || event.HLT_PFMET300_HBHECleaned_v7 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v1 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v2 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF_v10 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v2 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v3 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v4 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v5 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v6 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v7 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v13 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v14 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v15 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v16 > 0
                     || event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v9 > 0
                     || event.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v10 > 0
                     || event.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v12 > 0
                     || event.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v13 > 0
                     || event.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v14 > 0
                     || event.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v15 > 0
                     || event.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v9 > 0
                     || event.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v10 > 0
                     || event.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v12 > 0
                     || event.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v13 > 0
                     || event.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v14 > 0
                     || event.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v15 > 0
                     || event.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v16 > 0
                     || event.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v17 > 0
                     || event.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v9 > 0
                     || event.HLT_PFHT1050_v11 > 0
                     || event.HLT_PFHT1050_v12 > 0
                     || event.HLT_PFHT1050_v13 > 0
                     || event.HLT_PFHT1050_v14 > 0
                     || event.HLT_PFHT1050_v15 > 0
                     || event.HLT_PFHT1050_v16 > 0
                     || event.HLT_PFHT1050_v7 > 0 || event.HLT_PFHT1050_v9 > 0
                     || event.HLT_PFHT180_v7 > 0
                     || event.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v1 > 0
                     || event.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v10 > 0
                     || event.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v3 > 0
                     || event.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v5 > 0
                     || event.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v6 > 0
                     || event.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v7 > 0
                     || event.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v8 > 0
                     || event.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v9 > 0
                     || event.HLT_PFHT500_PFMET110_PFMHT110_IDTight_v1 > 0
                     || event.HLT_PFHT500_PFMET110_PFMHT110_IDTight_v10 > 0
                     || event.HLT_PFHT500_PFMET110_PFMHT110_IDTight_v3 > 0
                     || event.HLT_PFHT500_PFMET110_PFMHT110_IDTight_v5 > 0
                     || event.HLT_PFHT500_PFMET110_PFMHT110_IDTight_v6 > 0
                     || event.HLT_PFHT500_PFMET110_PFMHT110_IDTight_v7 > 0
                     || event.HLT_PFHT500_PFMET110_PFMHT110_IDTight_v8 > 0
                     || event.HLT_PFHT500_PFMET110_PFMHT110_IDTight_v9 > 0
                     || event.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v1 > 0
                     || event.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v10 > 0
                     || event.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v3 > 0
                     || event.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v5 > 0
                     || event.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v6 > 0
                     || event.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v7 > 0
                     || event.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v8 > 0
                     || event.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v9 > 0
                     || event.HLT_PFHT700_PFMET95_PFMHT95_IDTight_v1 > 0
                     || event.HLT_PFHT700_PFMET95_PFMHT95_IDTight_v10 > 0
                     || event.HLT_PFHT700_PFMET95_PFMHT95_IDTight_v3 > 0
                     || event.HLT_PFHT700_PFMET95_PFMHT95_IDTight_v5 > 0
                     || event.HLT_PFHT700_PFMET95_PFMHT95_IDTight_v6 > 0
                     || event.HLT_PFHT700_PFMET95_PFMHT95_IDTight_v7 > 0
                     || event.HLT_PFHT700_PFMET95_PFMHT95_IDTight_v8 > 0
                     || event.HLT_PFHT700_PFMET95_PFMHT95_IDTight_v9 > 0
                     || event.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v1 > 0
                     || event.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v10 > 0
                     || event.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v3 > 0
                     || event.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v5 > 0
                     || event.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v6 > 0
                     || event.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v7 > 0
                     || event.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v8 > 0
                     || event.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v9 > 0
                     || event.HLT_PFHT800_PFMET85_PFMHT85_IDTight_v1 > 0
                     || event.HLT_PFHT800_PFMET85_PFMHT85_IDTight_v10 > 0
                     || event.HLT_PFHT800_PFMET85_PFMHT85_IDTight_v3 > 0
                     || event.HLT_PFHT800_PFMET85_PFMHT85_IDTight_v5 > 0
                     || event.HLT_PFHT800_PFMET85_PFMHT85_IDTight_v6 > 0
                     || event.HLT_PFHT800_PFMET85_PFMHT85_IDTight_v7 > 0
                     || event.HLT_PFHT800_PFMET85_PFMHT85_IDTight_v8 > 0
                     || event.HLT_PFHT800_PFMET85_PFMHT85_IDTight_v9 > 0;
    // clang-format on
}

bool TriggerScaleFactors::metFilters(const AnalysisEvent& event,
                                     const bool isMC) const
{
    if (event.Flag_HBHENoiseFilter <= 0 || event.Flag_HBHENoiseIsoFilter <= 0
        || event.Flag_globalTightHalo2016Filter <= 0
        || event.Flag_EcalDeadCellTriggerPrimitiveFilter <= 0
        || event.Flag_goodVertices <= 0
        || (!isMC && event.Flag_eeBadScFilter <= 0))
    {
        return false;
    }

    if (is2016_
        && (event.Flag_ecalLaserCorrFilter <= 0
            || event.Flag_chargedHadronTrackResolutionFilter <= 0
            || event.Flag_muonBadTrackFilter <= 0
            || (!isMC && event.Flag_noBadMuons <= 0)))
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

bool TriggerScaleFactors::makeJetCuts(AnalysisEvent& event,
                                      const bool isMC) const
{
    std::vector<int> jets;
    for (int i{0}; i < event.numJetPF2PAT; i++)
    {
        // if (std::sqrt(event.jetPF2PATPx[i] * event.jetPF2PATPx[i] +
        // event.jetPF2PATPy[i] * event.jetPF2PATPy[i]) < jetPt_) continue;
        TLorentzVector jetVec{getJetLVec(event, i, isMC)};
        if (jetVec.Pt() <= 30.)
        {
            continue;
        }
        if (std::abs(jetVec.Eta()) >= 4.70)
        {
            continue;
        }

        bool jetId{true};

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
        else if (std::abs(jetVec.Eta()) <= 3.0 && std::abs(jetVec.Eta()) > 2.70)
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

        if (!jetId)
        {
            continue;
        }
        double deltaLep{10000};

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
        if (deltaLep < 0.4)
        {
            continue; // Only start rejecting things when actually making the
                      // jet cuts!
        }

        jets.emplace_back(i);
    }

    if (jets.size() > 6)
    {
        return false;
    }
    else if (jets.size() < 4)
    {
        return false;
    }

    // if b cuts
    if (bCuts_)
    {
        std::vector<int> bJets;
        for (unsigned int i = 0; i != jets.size(); i++)
        {
            TLorentzVector jetVec{getJetLVec(event, i, isMC)};
            if (event.jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                    [jets[i]]
                <= 0.8484)
            {
                continue;
            }
            if (jetVec.Eta() >= 2.40)
            {
                continue;
            }
            bJets.emplace_back(i);
        }
        if (bJets.size() > 2)
        {
            return false;
        }
        else if (bJets.size() < 1)
        {
            return false;
        }
    }

    return true;
}

TLorentzVector TriggerScaleFactors::getJetLVec(AnalysisEvent& event,
                                               const int index,
                                               const bool isMC_) const
{
    TLorentzVector returnJet;
    float newSmearValue{1.0};
    if (!isMC_)
    {
        event.jetSmearValue.emplace_back(newSmearValue);
        returnJet.SetPxPyPzE(event.jetPF2PATPx[index],
                             event.jetPF2PATPy[index],
                             event.jetPF2PATPz[index],
                             event.jetPF2PATE[index]);
        return returnJet;
    }
    float jerSF{0.};
    float jerSigma{0.};

    double eta = std::abs(event.jetPF2PATEta[index]);

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

    return returnJet;
}

double TriggerScaleFactors::deltaR(const float eta1,
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

void TriggerScaleFactors::savePlots()
{
    double level = 0.60; // ClopperPearson interval level

    // Histos first
    TFile* outFile{
        new TFile{(outFolder + "triggerPlots.root").c_str(), "RECREATE"}};

    // Do pT errors
    for (Int_t bin = 1; bin != numPt_bins + 1; bin++)
    {
        double errUp, errDown, error;

        // electrons MC
        // ele 1 pT MC
        errUp = (p_electron1_pT_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_electron1_pT_MC->GetBinEntries(bin),
                       p_electron1_pT_MC->GetBinEntries(bin)
                           * p_electron1_pT_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_electron1_pT_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_electron1_pT_MC->GetBinEntries(bin),
                         p_electron1_pT_MC->GetBinEntries(bin)
                             * p_electron1_pT_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_electron1_pT_MC->SetBinError(bin, error);
        // ele 2 pT MC
        errUp = (p_electron2_pT_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_electron2_pT_MC->GetBinEntries(bin),
                       p_electron2_pT_MC->GetBinEntries(bin)
                           * p_electron2_pT_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_electron2_pT_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_electron2_pT_MC->GetBinEntries(bin),
                         p_electron2_pT_MC->GetBinEntries(bin)
                             * p_electron2_pT_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_electron2_pT_MC->SetBinError(bin, error);

        // muons MC
        // muon 1 pT MC
        errUp = (p_muon1_pT_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muon1_pT_MC->GetBinEntries(bin),
                       p_muon1_pT_MC->GetBinEntries(bin)
                           * p_muon1_pT_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muon1_pT_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muon1_pT_MC->GetBinEntries(bin),
                         p_muon1_pT_MC->GetBinEntries(bin)
                             * p_muon1_pT_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muon1_pT_MC->SetBinError(bin, error);
        // muon 2 pT MC
        errUp = (p_muon2_pT_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muon2_pT_MC->GetBinEntries(bin),
                       p_muon2_pT_MC->GetBinEntries(bin)
                           * p_muon2_pT_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muon2_pT_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muon2_pT_MC->GetBinEntries(bin),
                         p_muon2_pT_MC->GetBinEntries(bin)
                             * p_muon2_pT_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muon2_pT_MC->SetBinError(bin, error);

        // muonEG MC
        // muonEG 1 pT MC
        errUp = (p_muonElectron1_pT_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muonElectron1_pT_MC->GetBinEntries(bin),
                       p_muonElectron1_pT_MC->GetBinEntries(bin)
                           * p_muonElectron1_pT_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muonElectron1_pT_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muonElectron1_pT_MC->GetBinEntries(bin),
                         p_muonElectron1_pT_MC->GetBinEntries(bin)
                             * p_muonElectron1_pT_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muonElectron1_pT_MC->SetBinError(bin, error);
        // muonEG 2 pT MC
        errUp = (p_muonElectron2_pT_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muonElectron2_pT_MC->GetBinEntries(bin),
                       p_muonElectron2_pT_MC->GetBinEntries(bin)
                           * p_muonElectron2_pT_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muonElectron2_pT_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muonElectron2_pT_MC->GetBinEntries(bin),
                         p_muonElectron2_pT_MC->GetBinEntries(bin)
                             * p_muonElectron2_pT_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muonElectron2_pT_MC->SetBinError(bin, error);

        // electrons data
        // ele 1 pT data
        errUp = (p_electron1_pT_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_electron1_pT_data->GetBinEntries(bin),
                       p_electron1_pT_data->GetBinEntries(bin)
                           * p_electron1_pT_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_electron1_pT_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_electron1_pT_data->GetBinEntries(bin),
                         p_electron1_pT_data->GetBinEntries(bin)
                             * p_electron1_pT_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_electron1_pT_data->SetBinError(bin, error);
        // ele 2 pT data
        errUp = (p_electron2_pT_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_electron2_pT_data->GetBinEntries(bin),
                       p_electron2_pT_data->GetBinEntries(bin)
                           * p_electron2_pT_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_electron2_pT_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_electron2_pT_data->GetBinEntries(bin),
                         p_electron2_pT_data->GetBinEntries(bin)
                             * p_electron2_pT_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_electron2_pT_data->SetBinError(bin, error);

        // muons data
        // muon 1 pT data
        errUp = (p_muon1_pT_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muon1_pT_data->GetBinEntries(bin),
                       p_muon1_pT_data->GetBinEntries(bin)
                           * p_muon1_pT_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muon1_pT_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muon1_pT_data->GetBinEntries(bin),
                         p_muon1_pT_data->GetBinEntries(bin)
                             * p_muon1_pT_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muon1_pT_data->SetBinError(bin, error);
        // muon 2 pT data
        errUp = (p_muon2_pT_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muon2_pT_data->GetBinEntries(bin),
                       p_muon2_pT_data->GetBinEntries(bin)
                           * p_muon2_pT_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muon2_pT_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muon2_pT_data->GetBinEntries(bin),
                         p_muon2_pT_data->GetBinEntries(bin)
                             * p_muon2_pT_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muon2_pT_data->SetBinError(bin, error);

        // muonEG data
        // muonEG 1 pT data
        errUp = (p_muonElectron1_pT_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muonElectron1_pT_data->GetBinEntries(bin),
                       p_muonElectron1_pT_data->GetBinEntries(bin)
                           * p_muonElectron1_pT_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muonElectron1_pT_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muonElectron1_pT_data->GetBinEntries(bin),
                         p_muonElectron1_pT_data->GetBinEntries(bin)
                             * p_muonElectron1_pT_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muonElectron1_pT_data->SetBinError(bin, error);
        // muonEG 2 pT data
        errUp = (p_muonElectron2_pT_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muonElectron2_pT_data->GetBinEntries(bin),
                       p_muonElectron2_pT_data->GetBinEntries(bin)
                           * p_muonElectron2_pT_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muonElectron2_pT_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muonElectron2_pT_data->GetBinEntries(bin),
                         p_muonElectron2_pT_data->GetBinEntries(bin)
                             * p_muonElectron2_pT_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muonElectron2_pT_data->SetBinError(bin, error);
    }

    // Do eta errors
    for (Int_t bin = 1; bin != numEta_bins + 1; bin++)
    {
        double errUp, errDown, error;

        // electrons MC
        // ele 1 eta MC
        errUp = (p_electron1_eta_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_electron1_eta_MC->GetBinEntries(bin),
                       p_electron1_eta_MC->GetBinEntries(bin)
                           * p_electron1_eta_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_electron1_eta_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_electron1_eta_MC->GetBinEntries(bin),
                         p_electron1_eta_MC->GetBinEntries(bin)
                             * p_electron1_eta_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_electron1_eta_MC->SetBinError(bin, error);
        // ele 2 eta MC
        errUp = (p_electron2_eta_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_electron2_eta_MC->GetBinEntries(bin),
                       p_electron2_eta_MC->GetBinEntries(bin)
                           * p_electron2_eta_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_electron2_eta_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_electron2_eta_MC->GetBinEntries(bin),
                         p_electron2_eta_MC->GetBinEntries(bin)
                             * p_electron2_eta_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_electron2_eta_MC->SetBinError(bin, error);

        // muons MC
        // muon 1 eta MC
        errUp = (p_muon1_eta_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muon1_eta_MC->GetBinEntries(bin),
                       p_muon1_eta_MC->GetBinEntries(bin)
                           * p_muon1_eta_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muon1_eta_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muon1_eta_MC->GetBinEntries(bin),
                         p_muon1_eta_MC->GetBinEntries(bin)
                             * p_muon1_eta_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muon1_eta_MC->SetBinError(bin, error);
        // muon 2 eta MC
        errUp = (p_muon2_eta_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muon2_eta_MC->GetBinEntries(bin),
                       p_muon2_eta_MC->GetBinEntries(bin)
                           * p_muon2_eta_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muon2_eta_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muon2_eta_MC->GetBinEntries(bin),
                         p_muon2_eta_MC->GetBinEntries(bin)
                             * p_muon2_eta_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muon2_eta_MC->SetBinError(bin, error);

        // muonEG MC
        // muonEG 1 eta MC
        errUp = (p_muonElectron1_eta_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muonElectron1_eta_MC->GetBinEntries(bin),
                       p_muonElectron1_eta_MC->GetBinEntries(bin)
                           * p_muonElectron1_eta_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muonElectron1_eta_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muonElectron1_eta_MC->GetBinEntries(bin),
                         p_muonElectron1_eta_MC->GetBinEntries(bin)
                             * p_muonElectron1_eta_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muonElectron1_eta_MC->SetBinError(bin, error);
        // muonEG 2 eta MC
        errUp = (p_muonElectron2_eta_MC->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muonElectron2_eta_MC->GetBinEntries(bin),
                       p_muonElectron2_eta_MC->GetBinEntries(bin)
                           * p_muonElectron2_eta_MC->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muonElectron2_eta_MC->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muonElectron2_eta_MC->GetBinEntries(bin),
                         p_muonElectron2_eta_MC->GetBinEntries(bin)
                             * p_muonElectron2_eta_MC->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muonElectron2_eta_MC->SetBinError(bin, error);

        // electrons data
        // ele 1 eta data
        errUp = (p_electron1_eta_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_electron1_eta_data->GetBinEntries(bin),
                       p_electron1_eta_data->GetBinEntries(bin)
                           * p_electron1_eta_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_electron1_eta_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_electron1_eta_data->GetBinEntries(bin),
                         p_electron1_eta_data->GetBinEntries(bin)
                             * p_electron1_eta_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_electron1_eta_data->SetBinError(bin, error);
        // ele 2 eta data
        errUp = (p_electron2_eta_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_electron2_eta_data->GetBinEntries(bin),
                       p_electron2_eta_data->GetBinEntries(bin)
                           * p_electron2_eta_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_electron2_eta_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_electron2_eta_data->GetBinEntries(bin),
                         p_electron2_eta_data->GetBinEntries(bin)
                             * p_electron2_eta_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_electron2_eta_data->SetBinError(bin, error);

        // muons data
        // muon 1 eta data
        errUp = (p_muon1_eta_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muon1_eta_data->GetBinEntries(bin),
                       p_muon1_eta_data->GetBinEntries(bin)
                           * p_muon1_eta_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muon1_eta_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muon1_eta_data->GetBinEntries(bin),
                         p_muon1_eta_data->GetBinEntries(bin)
                             * p_muon1_eta_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muon1_eta_data->SetBinError(bin, error);
        // muon 2 eta data
        errUp = (p_muon2_eta_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muon2_eta_data->GetBinEntries(bin),
                       p_muon2_eta_data->GetBinEntries(bin)
                           * p_muon2_eta_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muon2_eta_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muon2_eta_data->GetBinEntries(bin),
                         p_muon2_eta_data->GetBinEntries(bin)
                             * p_muon2_eta_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muon2_eta_data->SetBinError(bin, error);

        // muonEG data
        // muonEG 1 eta data
        errUp = (p_muonElectron1_eta_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muonElectron1_eta_data->GetBinEntries(bin),
                       p_muonElectron1_eta_data->GetBinEntries(bin)
                           * p_muonElectron1_eta_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muonElectron1_eta_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muonElectron1_eta_data->GetBinEntries(bin),
                         p_muonElectron1_eta_data->GetBinEntries(bin)
                             * p_muonElectron1_eta_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muonElectron1_eta_data->SetBinError(bin, error);
        // muonEG 2 eta data
        errUp = (p_muonElectron2_eta_data->GetBinContent(bin)
                 - TEfficiency::ClopperPearson(
                       p_muonElectron2_eta_data->GetBinEntries(bin),
                       p_muonElectron2_eta_data->GetBinEntries(bin)
                           * p_muonElectron2_eta_data->GetBinContent(bin),
                       level,
                       true));
        errDown = (p_muonElectron2_eta_data->GetBinContent(bin)
                   - TEfficiency::ClopperPearson(
                         p_muonElectron2_eta_data->GetBinEntries(bin),
                         p_muonElectron2_eta_data->GetBinEntries(bin)
                             * p_muonElectron2_eta_data->GetBinContent(bin),
                         level,
                         false));
        error = errUp > errDown ? errUp : errDown;
        p_muonElectron2_eta_data->SetBinError(bin, error);
    }

    // Turn on curves ratio histos
    TProfile* p_electron1_pT_data_MC = dynamic_cast<TProfile*>(
        p_electron1_pT_data->Clone("p_electron1_pT_data_MC"));
    p_electron1_pT_data_MC->Divide(p_electron1_pT_MC);

    TProfile* p_electron1_eta_data_MC = dynamic_cast<TProfile*>(
        p_electron1_eta_data->Clone("p_electron1_eta_data_MC"));
    p_electron1_eta_data_MC->Divide(p_electron1_eta_MC);

    TProfile* p_electron2_pT_data_MC = dynamic_cast<TProfile*>(
        p_electron2_pT_data->Clone("p_electron2_pT_data_MC"));
    p_electron2_pT_data_MC->Divide(p_electron2_pT_MC);

    TProfile* p_electron2_eta_data_MC = dynamic_cast<TProfile*>(
        p_electron2_eta_data->Clone("p_electron2_eta_data_MC"));
    p_electron2_eta_data_MC->Divide(p_electron2_eta_MC);

    TProfile* p_muon1_pT_data_MC =
        dynamic_cast<TProfile*>(p_muon1_pT_data->Clone("p_muon1_pT_data_MC"));
    p_muon1_pT_data_MC->Divide(p_muon1_pT_MC);

    TProfile* p_muon1_eta_data_MC =
        dynamic_cast<TProfile*>(p_muon1_eta_data->Clone("p_muon1_eta_data_MC"));
    p_muon1_eta_data_MC->Divide(p_muon1_eta_MC);

    TProfile* p_muon2_pT_data_MC =
        dynamic_cast<TProfile*>(p_muon2_pT_data->Clone("p_muon2_pT_data_MC"));
    p_muon2_pT_data_MC->Divide(p_muon2_pT_MC);

    TProfile* p_muon2_eta_data_MC =
        dynamic_cast<TProfile*>(p_muon2_eta_data->Clone("p_muon2_eta_data_MC"));
    p_muon2_eta_data_MC->Divide(p_muon2_eta_MC);

    TProfile* p_muonElectron1_pT_data_MC = dynamic_cast<TProfile*>(
        p_muonElectron1_pT_data->Clone("p_muonElectron1_pT_data_MC"));
    p_muonElectron1_pT_data_MC->Divide(p_muonElectron1_pT_MC);

    TProfile* p_muonElectron1_eta_data_MC = dynamic_cast<TProfile*>(
        p_muonElectron1_eta_data->Clone("p_muonElectron1_eta_data_MC"));
    p_muonElectron1_eta_data_MC->Divide(p_muonElectron1_eta_MC);

    TProfile* p_muonElectron2_pT_data_MC = dynamic_cast<TProfile*>(
        p_muonElectron2_pT_data->Clone("p_muonElectron2_pt_data_MC"));
    p_muonElectron2_pT_data_MC->Divide(p_muonElectron2_pT_MC);

    TProfile* p_muonElectron2_eta_data_MC = dynamic_cast<TProfile*>(
        p_muonElectron2_eta_data->Clone("p_muonElectron2_eta_data_MC"));
    p_muonElectron2_eta_data_MC->Divide(p_muonElectron2_eta_MC);

    // SF histos

    TProfile2D* electronPtSF =
        dynamic_cast<TProfile2D*>(p_electrons_pT_data->Clone("electronPtSF"));
    electronPtSF->Divide(p_electrons_pT_MC);

    TProfile2D* electronEtaSF =
        dynamic_cast<TProfile2D*>(p_electrons_eta_data->Clone("electronEtaSF"));
    electronEtaSF->Divide(p_electrons_eta_MC);

    TProfile2D* muonPtSF =
        dynamic_cast<TProfile2D*>(p_muons_pT_data->Clone("muonPtSF"));
    muonPtSF->Divide(p_muons_pT_MC);

    TProfile2D* muonEtaSF =
        dynamic_cast<TProfile2D*>(p_muons_eta_data->Clone("muonEtaSF"));
    muonEtaSF->Divide(p_muons_eta_MC);

    TProfile2D* muonElectronPtSF = dynamic_cast<TProfile2D*>(
        p_muonElectrons_pT_data->Clone("muonElectronPtSF"));
    muonElectronPtSF->Divide(p_muonElectrons_pT_MC);

    TProfile2D* muonElectronEtaSF = dynamic_cast<TProfile2D*>(
        p_muonElectrons_eta_data->Clone("muonElectronEtaSF"));
    muonElectronEtaSF->Divide(p_muonElectrons_eta_MC);

    // Write Histos

    p_electron1_pT_MC->Write();
    p_electron1_eta_MC->Write();
    p_electron2_pT_MC->Write();
    p_electron2_eta_MC->Write();
    p_muon1_pT_MC->Write();
    p_muon1_eta_MC->Write();
    p_muon2_pT_MC->Write();
    p_muon2_eta_MC->Write();
    p_muonElectron1_pT_MC->Write();
    p_muonElectron1_eta_MC->Write();
    p_muonElectron2_pT_MC->Write();
    p_muonElectron2_eta_MC->Write();

    p_electron1_pT_data->Write();
    p_electron1_eta_data->Write();
    p_electron2_pT_data->Write();
    p_electron2_eta_data->Write();
    p_muon1_pT_data->Write();
    p_muon1_eta_data->Write();
    p_muon2_pT_data->Write();
    p_muon2_eta_data->Write();
    p_muonElectron1_pT_data->Write();
    p_muonElectron1_eta_data->Write();
    p_muonElectron2_pT_data->Write();
    p_muonElectron2_eta_data->Write();

    p_electrons_pT_MC->Write();
    p_electrons_eta_MC->Write();
    p_muons_pT_MC->Write();
    p_muons_eta_MC->Write();
    p_muonElectrons_pT_MC->Write();
    p_muonElectrons_eta_MC->Write();

    p_electrons_pT_data->Write();
    p_electrons_eta_data->Write();
    p_muons_pT_data->Write();
    p_muons_eta_data->Write();
    p_muonElectrons_pT_data->Write();
    p_muonElectrons_eta_data->Write();

    p_electron1_pT_data_MC->Write();
    p_electron1_eta_data_MC->Write();
    p_electron2_pT_data_MC->Write();
    p_electron2_eta_data_MC->Write();

    p_muon1_pT_data_MC->Write();
    p_muon1_eta_data_MC->Write();
    p_muon2_pT_data_MC->Write();
    p_muon2_eta_data_MC->Write();

    p_muonElectron1_pT_data_MC->Write();
    p_muonElectron1_eta_data_MC->Write();
    p_muonElectron2_pT_data_MC->Write();
    p_muonElectron2_eta_data_MC->Write();

    electronPtSF->Write();
    electronEtaSF->Write();
    muonPtSF->Write();
    muonEtaSF->Write();
    muonElectronPtSF->Write();
    muonElectronEtaSF->Write();

    TCanvas* lCanvasEle1PtEff =
        new TCanvas("lCanvasEle1PtEff", "lCanvasEle1PtEff");
    lCanvasEle1PtEff->cd(1);
    p_electron1_pT_data->SetStats(false);
    p_electron1_pT_data->Draw();
    p_electron1_pT_data->SetLineColor(kBlue);
    p_electron1_pT_MC->SetStats(false);
    p_electron1_pT_MC->Draw("same");
    p_electron1_pT_MC->SetLineColor(kRed);
    lCanvasEle1PtEff->Write();

    TCanvas* lCanvasEle2PtEff =
        new TCanvas("lCanvasEle2PtEff", "lCanvasEle2PtEff");
    lCanvasEle2PtEff->cd(1);
    p_electron2_pT_data->SetStats(false);
    p_electron2_pT_data->Draw();
    p_electron2_pT_data->SetLineColor(kBlue);
    p_electron2_pT_MC->SetStats(false);
    p_electron2_pT_MC->Draw("same");
    p_electron2_pT_MC->SetLineColor(kRed);
    lCanvasEle2PtEff->Write();

    TCanvas* lCanvasEle1EtaEff =
        new TCanvas("lCanvasEle1EtaEff", "lCanvasEle1EtaEff");
    lCanvasEle1EtaEff->cd(1);
    p_electron1_eta_data->SetStats(false);
    p_electron1_eta_data->Draw();
    p_electron1_eta_data->SetLineColor(kBlue);
    p_electron1_eta_MC->SetStats(false);
    p_electron1_eta_MC->Draw("same");
    p_electron1_eta_MC->SetLineColor(kRed);
    lCanvasEle1EtaEff->Write();

    TCanvas* lCanvasEle2EtaEff =
        new TCanvas("lCanvasEle2EtaEff", "lCanvasEle2EtaEff");
    lCanvasEle2EtaEff->cd(1);
    p_electron2_eta_data->SetStats(false);
    p_electron2_eta_data->Draw();
    p_electron2_eta_data->SetLineColor(kBlue);
    p_electron2_eta_MC->SetStats(false);
    p_electron2_eta_MC->Draw("same");
    p_electron2_eta_MC->SetLineColor(kRed);
    lCanvasEle2EtaEff->Write();

    TCanvas* lCanvasMuon1PtEff =
        new TCanvas("lCanvasMuon1PtEff", "lCanvasMuon1PtEff");
    lCanvasMuon1PtEff->cd(1);
    p_muon1_pT_data->SetStats(false);
    p_muon1_pT_data->Draw();
    p_muon1_pT_data->SetLineColor(kBlue);
    p_muon1_pT_MC->SetStats(false);
    p_muon1_pT_MC->Draw("same");
    p_muon1_pT_MC->SetLineColor(kRed);
    lCanvasMuon1PtEff->Write();

    TCanvas* lCanvasMuon2PtEff =
        new TCanvas("lCanvasMuon2PtEff", "lCanvasMuon2PtEff");
    lCanvasMuon2PtEff->cd(1);
    p_muon2_pT_data->SetStats(false);
    p_muon2_pT_data->Draw();
    p_muon2_pT_data->SetLineColor(kBlue);
    p_muon2_pT_MC->SetStats(false);
    p_muon2_pT_MC->Draw("same");
    p_muon2_pT_MC->SetLineColor(kRed);
    lCanvasMuon2PtEff->Write();

    TCanvas* lCanvasMuon1EtaEff =
        new TCanvas("lCanvasMuon1EtaEff", "lCanvasMuon1EtaEff");
    lCanvasMuon1EtaEff->cd(1);
    p_muon1_eta_data->SetStats(false);
    p_muon1_eta_data->Draw();
    p_muon1_eta_data->SetLineColor(kBlue);
    p_muon1_eta_MC->SetStats(false);
    p_muon1_eta_MC->Draw("same");
    p_muon1_eta_MC->SetLineColor(kRed);
    lCanvasMuon1EtaEff->Write();

    TCanvas* lCanvasMuon2EtaEff =
        new TCanvas("lCanvasMuon2EtaEff", "lCanvasMuon2EtaEff");
    lCanvasMuon2EtaEff->cd(1);
    p_muon2_eta_data->SetStats(false);
    p_muon2_eta_data->Draw();
    p_muon2_eta_data->SetLineColor(kBlue);
    p_muon1_eta_MC->SetStats(false);
    p_muon2_eta_MC->Draw("same");
    p_muon2_eta_MC->SetLineColor(kRed);
    lCanvasMuon2EtaEff->Write();

    outFile->Close();

    // Calculate MC efficiency

    //// LeptonTriggers
    double doubleElectronEfficiencyMC =
        numberTriggeredDoubleElectrons[0] / (numberPassedElectrons[0] + 1.0e-6);
    double doubleMuonEfficiencyMC =
        numberTriggeredDoubleMuons[0] / (numberPassedMuons[0] + 1.0e-6);
    double muonElectronEfficiencyMC = numberTriggeredMuonElectrons[0]
                                      / (numberPassedMuonElectrons[0] + 1.0e-6);

    // Calculate Data efficiency

    //// DoubleLeptonTriggers
    double doubleElectronEfficiencyData =
        numberTriggeredDoubleElectrons[1] / (numberPassedElectrons[1] + 1.0e-6);
    double doubleMuonEfficiencyData =
        numberTriggeredDoubleMuons[1] / (numberPassedMuons[1] + 1.0e-6);
    double muonElectronEfficiencyData =
        numberTriggeredMuonElectrons[1]
        / (numberPassedMuonElectrons[1] + 1.0e-6);

    // Calculate SF

    //// LeptonTriggers
    double doubleElectronSF =
        doubleElectronEfficiencyData / (doubleElectronEfficiencyMC + 1.0e-6);
    double doubleMuonSF =
        doubleMuonEfficiencyData / (doubleMuonEfficiencyMC + 1.0e-6);
    double muonElectronSF =
        muonElectronEfficiencyData / (muonElectronEfficiencyMC + 1.0e-6);

    // Calculate alphas
    double alphaDoubleElectron =
        ((numberSelectedDoubleElectronsTriggered[0]
          / numberSelectedElectrons[0])
         * (numberPassedElectrons[0] / numberSelectedElectrons[0]))
        / (numberTriggeredDoubleElectrons[0] / numberSelectedElectrons[0]
           + 1.0e-6);
    double alphaDoubleMuon =
        ((numberSelectedDoubleMuonsTriggered[0] / numberSelectedMuons[0])
         * (numberPassedMuons[0] / numberSelectedMuons[0]))
        / (numberTriggeredDoubleMuons[0] / numberSelectedMuons[0] + 1.0e-6);
    double alphaMuonElectron =
        ((numberSelectedMuonElectronsTriggered[0]
          / numberSelectedMuonElectrons[0])
         * (numberPassedMuonElectrons[0] / numberSelectedMuonElectrons[0]))
        / (numberTriggeredMuonElectrons[0] / numberSelectedMuonElectrons[0]
           + 1.0e-6);

    // Calculate uncertainities

    //// LeptonTriggers

    double doubleElectronDataUpperUncert =
        doubleElectronEfficiencyData
        - TEfficiency::ClopperPearson(numberPassedElectrons[1],
                                      numberTriggeredDoubleElectrons[1],
                                      level,
                                      true);
    double doubleElectronMcUpperUncert =
        doubleElectronEfficiencyMC
        - TEfficiency::ClopperPearson(numberPassedElectrons[0],
                                      numberTriggeredDoubleElectrons[0],
                                      level,
                                      true);
    double doubleElectronDataLowerUncert =
        doubleElectronEfficiencyData
        - TEfficiency::ClopperPearson(numberPassedElectrons[1],
                                      numberTriggeredDoubleElectrons[1],
                                      level,
                                      false);
    double doubleElectronMcLowerUncert =
        doubleElectronEfficiencyMC
        - TEfficiency::ClopperPearson(numberPassedElectrons[0],
                                      numberTriggeredDoubleElectrons[0],
                                      level,
                                      false);

    double doubleMuonDataUpperUncert =
        doubleMuonEfficiencyData
        - TEfficiency::ClopperPearson(
              numberPassedMuons[1], numberTriggeredDoubleMuons[1], level, true);
    double doubleMuonMcUpperUncert =
        doubleMuonEfficiencyMC
        - TEfficiency::ClopperPearson(
              numberPassedMuons[0], numberTriggeredDoubleMuons[0], level, true);
    double doubleMuonDataLowerUncert =
        doubleMuonEfficiencyData
        - TEfficiency::ClopperPearson(numberPassedMuons[1],
                                      numberTriggeredDoubleMuons[1],
                                      level,
                                      false);
    double doubleMuonMcLowerUncert =
        doubleMuonEfficiencyMC
        - TEfficiency::ClopperPearson(numberPassedMuons[0],
                                      numberTriggeredDoubleMuons[0],
                                      level,
                                      false);

    double muonElectronDataUpperUncert =
        muonElectronEfficiencyData
        - TEfficiency::ClopperPearson(numberPassedMuonElectrons[1],
                                      numberTriggeredMuonElectrons[1],
                                      level,
                                      true);
    double muonElectronMcUpperUncert =
        muonElectronEfficiencyMC
        - TEfficiency::ClopperPearson(numberPassedMuonElectrons[0],
                                      numberTriggeredMuonElectrons[0],
                                      level,
                                      true);
    double muonElectronDataLowerUncert =
        muonElectronEfficiencyData
        - TEfficiency::ClopperPearson(numberPassedMuonElectrons[1],
                                      numberTriggeredMuonElectrons[1],
                                      level,
                                      false);
    double muonElectronMcLowerUncert =
        muonElectronEfficiencyMC
        - TEfficiency::ClopperPearson(numberPassedMuonElectrons[0],
                                      numberTriggeredMuonElectrons[0],
                                      level,
                                      false);

    double doubleEleSfUp =
        (doubleElectronEfficiencyData + doubleElectronDataUpperUncert)
            / (doubleElectronEfficiencyMC - doubleElectronMcLowerUncert
               + 1.0e-6)
        - doubleElectronSF;
    double doubleEleSfDown =
        (doubleElectronEfficiencyData + doubleElectronDataLowerUncert)
            / (doubleElectronEfficiencyMC - doubleElectronMcUpperUncert
               + 1.0e-6)
        - doubleElectronSF;
    double doubleEleSfUncert = 0.0;
    if (doubleEleSfUp > doubleEleSfDown)
    {
        doubleEleSfUncert = doubleEleSfUp;
    }
    else
    {
        doubleEleSfUncert = doubleEleSfDown;
    }

    double doubleMuonSfUp =
        (doubleMuonEfficiencyData + doubleMuonDataUpperUncert)
            / (doubleMuonEfficiencyMC - doubleMuonMcLowerUncert + 1.0e-6)
        - doubleMuonSF;
    double doubleMuonSfDown =
        (doubleMuonEfficiencyData + doubleMuonDataLowerUncert)
            / (doubleMuonEfficiencyMC - doubleMuonMcUpperUncert + 1.0e-6)
        - doubleMuonSF;
    double doubleMuonSfUncert = 0.0;
    if (doubleMuonSfUp > doubleMuonSfDown)
    {
        doubleMuonSfUncert = doubleMuonSfUp;
    }
    else
    {
        doubleMuonSfUncert = doubleMuonSfDown;
    }

    double muonElectronSfUp =
        (muonElectronEfficiencyData + muonElectronDataUpperUncert)
            / (muonElectronEfficiencyMC - muonElectronMcLowerUncert + 1.0e-6)
        - muonElectronSF;
    double muonElectronSfDown =
        (muonElectronEfficiencyData + muonElectronDataLowerUncert)
            / (muonElectronEfficiencyMC - muonElectronMcUpperUncert + 1.0e-6)
        - muonElectronSF;
    double muonElectronSfUncert = 0.0;
    if (muonElectronSfUp > muonElectronSfDown)
    {
        muonElectronSfUncert = muonElectronSfUp;
    }
    else
    {
        muonElectronSfUncert = muonElectronSfDown;
    }

    // Print output

    std::cout << "-----------------------------------------------------------"
              << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
    std::cout << "Double Electron data efficiency: "
              << doubleElectronEfficiencyData << " +/- "
              << doubleElectronDataUpperUncert << "/"
              << doubleElectronDataLowerUncert << std::endl;
    std::cout << "Double Electron MC efficiency: " << doubleElectronEfficiencyMC
              << " +/- " << doubleElectronMcUpperUncert << "/"
              << doubleElectronMcLowerUncert << std::endl;
    std::cout << "Double Electron trigger SF: " << doubleElectronSF << " +/- "
              << doubleEleSfUncert << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
    std::cout << "Double Muon data efficiency: " << doubleMuonEfficiencyData
              << " +/- " << doubleMuonDataUpperUncert << "/"
              << doubleMuonDataLowerUncert << std::endl;
    std::cout << "Double Muon MC efficiency: " << doubleMuonEfficiencyMC
              << " +/- " << doubleMuonMcUpperUncert << "/"
              << doubleMuonMcLowerUncert << std::endl;
    std::cout << "Double Muon trigger SF: " << doubleMuonSF << " +/- "
              << doubleMuonSfUncert << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
    std::cout << "MuonEG data efficiency: " << muonElectronEfficiencyData
              << " +/- " << muonElectronDataUpperUncert << "/"
              << muonElectronDataLowerUncert << std::endl;
    std::cout << "MuonEG MC efficiency: " << muonElectronEfficiencyMC << " +/- "
              << muonElectronMcUpperUncert << " / " << muonElectronMcLowerUncert
              << std::endl;
    std::cout << "MuonEG trigger SF: " << muonElectronSF << " +/- "
              << muonElectronSfUncert << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
    std::cout << "alpha for DoubleEG/DoubleMuon/MuonEG Triggers: "
              << alphaDoubleElectron << "/" << alphaDoubleMuon << "/"
              << alphaMuonElectron << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
    std::cout << "-----------------------------------------------------------"
              << std::endl;
}
