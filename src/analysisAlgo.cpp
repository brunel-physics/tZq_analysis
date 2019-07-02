#include "AnalysisEvent.hpp"
#include "Compression.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TMVA/Config.h"
#include "TMVA/Timer.h"
#include "TPad.h"
#include "TTree.h"
#include "analysisAlgo.hpp"
#include "config_parser.hpp"

#include <LHAPDF/LHAPDF.h>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

AnalysisAlgo::AnalysisAlgo()
    : plots{false}
    , makeHistos{false}
    , useHistos{false}
    , channel{}
    , cutConfName{}
    , plotConfName{}
    , customJetRegion{false}
    , is2016_{false}
    , doNPLs_{false}
    , doZplusCR_{false}
{
}

AnalysisAlgo::~AnalysisAlgo()
{
}

void AnalysisAlgo::parseCommandLineArguements(int argc, char* argv[])
{
    std::stringstream events;
    std::stringstream jetRegion;

    gErrorIgnoreLevel = kInfo;
    // Set up environment a little.
    std::cerr << std::setprecision(6) << std::fixed;
    std::cout << std::setprecision(6) << std::fixed;

    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()("help,h", "Print this message.")(
        "config,c",
        po::value<std::string>(&config)->required(),
        "The configuration file to be used.")(
        "2016",
        po::bool_switch(&is2016_),
        "Use 2016 conditions (SFs, et al.).")(
        ",n",
        po::value<long>(&nEvents)->default_value(0),
        "The number of events to be run over. All if set to 0.")(
        "allPlots,p", po::bool_switch(&plots), "Make all plots")(
        "makeHistos",
        po::bool_switch(&makeHistos),
        "Make histos to be used in future plots")(
        "useHistos",
        po::bool_switch(&useHistos),
        "Use saved histos to make plots")(
        "histoDir",
        po::value<std::string>(&histoDir)->default_value("histos/mz20mw50/"),
        "The output directory for the histos used to make the plots.")(
        "outFolder,o",
        po::value<std::string>(&outFolder)->default_value("plots/"),
        "The output directory for the plots. Overrides the config file.")(
        "postfix,s",
        po::value<std::string>(&postfix)->default_value("default"),
        "Set postfix for plots. Overrides the config file.")(
        "lumi,l",
        po::value<double>(&usePreLumi)->default_value(35860.066),
        "Lumi to scale MC plots to.")(
        "cutConf,x",
        po::value<std::string>(&cutConfName),
        "Override the cut configuration given in the config file.")(
        "plotConf",
        po::value<std::string>(&plotConfName),
        "Override the plot configuration given in the config file. Sets "
        "--allPlots.")("invert,i",
                       po::bool_switch(&invertLepCut),
                       "Inverts the different charge cut for leptons.")(
        "MC,m",
        po::bool_switch(&skipData),
        "Monte Carlo only mode. Ignores all data in the config file.")(
        "data,b",
        po::bool_switch(&skipMC),
        "Data only mode. Ignores all data in the config file.")(
        "bTag,t",
        po::bool_switch(&usebTagWeight),
        "Use b-tagging efficiencies to reweight the Monte Carlo. Currently "
        "requires -u.")(
        "NPLs", po::bool_switch(&doNPLs_), "Make or use NPL shapes")(
        "zPlus",
        po::bool_switch(&doZplusCR_),
        "Use Z+jets CR for dilepton channel. Region mwCut and metCut set by "
        "--mwCut and --metCut.")(
        "nFiles,f",
        po::value<int>(&numFiles)->default_value(-1),
        "Number of files to run over. All if set to -1.")(
        "events,e",
        po::value<std::vector<int>>(&eventNumbers)->multitoken(),
        "Specify a space-separated list of events to run over.")(
        ",g",
        po::bool_switch(&makePostLepTree),
        "Make post lepton selection trees and bTag efficiencies.")(
        ",u",
        po::bool_switch(&usePostLepTree),
        "Use post lepton selection trees.")(
        "makeMVATree,z",
        po::bool_switch(&makeMVATree),
        "Produce trees after event selection for multivariate analysis.")(
        "syst,v",
        po::value<int>(&systToRun)->default_value(0),
        "Mask for systematics to be run. 65535 enables all systematics.")(
        "channels,k",
        po::value<int>(&channelsToRun)->default_value(0),
        "Mask describing the channels to be run over. The mask "
        "is the sum of each channel's mask, which are:\n"
        "    ee        - 1\n"
        "    mumu      - 2\n"
        "    ee ss     - 4\n"
        "    mumu ss   - 8\n"
        "    emu       - 16\n"
        "    emu ss    - 32\n"
        "0 runs the channels specified in the config file")(
        "skipTrig", po::bool_switch(&skipTrig), "Skip running triggers.")(
        "mvaDir",
        po::value<std::string>(&mvaDir),
        "Output directory for the MVA files.")(
        "jetRegion",
        po::value<std::vector<unsigned>>(&jetRegVars),
        "Set a sustom jet region in the format NJETS NBJETS MAXJETS MAXBJETS.")(
        "metCut",
        po::value<float>(&metCut)->default_value(0),
        "Apply an MET cut. Dilepton Z+jets CR.")(
        "mzCut",
        po::value<float>(&mzCut)->default_value(20.),
        "Apply an mZ cut. Dilepton only.")(
        "mwCut",
        po::value<float>(&mwCut)->default_value(20.),
        "Apply an mW cut. Dilepton only.");
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

        if (vm.count("channels") && !vm.count("config"))
        {
            throw std::logic_error(
                "--channels requires --config to be specified");
        }
        if (vm.count("jetRegion"))
        {
            if (jetRegVars.size() != 4)
            {
                throw std::logic_error(
                    "--jetRegion takes exactly four arguments.");
            }
            customJetRegion = true;
            std::cout << "CAUTION! Using a custom jet region of "
                      << jetRegVars[0] << "-" << jetRegVars[2] << " jets, and "
                      << jetRegVars[1] << "-" << jetRegVars[3] << " b-jets"
                      << std::endl;
        }
        if (usebTagWeight && !usePostLepTree)
        {
            throw std::logic_error(
                "Currently bTag weights can only be retrieved "
                "from post lepton selection trees. Please set -u.");
        }
    }
    catch (const std::logic_error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "Use -h or --help for help." << std::endl;
        std::exit(1);
    }

    if (vm.count("plotConf"))
    {
        plots = true;
    }

    // Some vectors that will be filled in the parsing.
    totalLumi = 0;

    try
    {
        Parser::parse_config(config,
                             datasets,
                             totalLumi,
                             plotTitles,
                             plotNames,
                             xMin,
                             xMax,
                             nBins,
                             fillExp,
                             xAxisLabels,
                             cutStage,
                             cutConfName,
                             plotConfName,
                             outFolder,
                             postfix,
                             channel,
                             doNPLs_);
    }
    catch (const std::exception)
    {
        std::cerr << "ERROR Problem with a confugration file, see previous "
                     "errors for more details. If this is the only error, the "
                     "problem is with the main configuration file."
                  << std::endl;
        throw;
    }

    if (channelsToRun)
    {
        std::cout << "Running over the channels: " << std::endl;
        for (unsigned channelInd = 1; channelInd != 32;
             channelInd = channelInd << 1)
        {
            if (!(channelInd & channelsToRun) && channelsToRun)
                continue;
            if (channelInd & 5)
            { // ee channel
                std::cout << "ee ";
            }
            if (channelInd & 10)
            { // mumu channel
                std::cout << "mumu ";
            }
            if (channelInd & 3)
            { // nominal samples
                std::cout << "nominal" << std::endl;
            }
            if (channelInd & 12)
            { // same sign samples
                std::cout << "same lepton charge" << std::endl;
            }
            if (channelInd & 16)
            { // nominal samples and emu
                std::cout << "emu - used only for ttbar control region "
                          << std::endl;
            }
            if (channelInd & 32)
            { // nominal samples and emu
                std::cout
                    << "emu - used only for ttbar same sign control region "
                    << std::endl;
            }
        }
    }
}

void AnalysisAlgo::setupSystematics()
{
    systNames.emplace_back("");
    systNames.emplace_back("__trig__plus");
    systNames.emplace_back("__trig__minus");
    systNames.emplace_back("__jer__plus");
    systNames.emplace_back("__jer__minus");
    systNames.emplace_back("__jes__plus");
    systNames.emplace_back("__jes__minus");
    systNames.emplace_back("__pileup__plus");
    systNames.emplace_back("__pileup__minus");
    systNames.emplace_back("__bTag__plus");
    systNames.emplace_back("__bTag__minus");
    systNames.emplace_back("__pdf__plus");
    systNames.emplace_back("__pdf__minus");
    systNames.emplace_back("__ME__plus");
    systNames.emplace_back("__ME__minus");
    systNames.emplace_back("__alphaS__plus");
    systNames.emplace_back("__alphaS__minus");
    systNames.emplace_back("__isr__plus");
    systNames.emplace_back("__isr__minus");
    systNames.emplace_back("__fsr__plus");
    systNames.emplace_back("__fsr__minus");

    if (!is2016_)
    { // If 2017 mode, get 2017 PU
        // Make pileupReweighting stuff here
        dataPileupFile = new TFile{"pileup/2017/truePileupTest.root", "READ"};
        dataPU = dynamic_cast<TH1D*>(dataPileupFile->Get("pileup")->Clone());
        mcPileupFile = new TFile{"pileup/2017/pileupMC.root", "READ"};
        mcPU = dynamic_cast<TH1D*>(mcPileupFile->Get("pileup")->Clone());

        // Get systematic files too.
        systUpFile = new TFile{"pileup/2017/truePileupUp.root", "READ"};
        pileupUpHist = dynamic_cast<TH1D*>(systUpFile->Get("pileup")->Clone());
        systDownFile = new TFile{"pileup/2017/truePileupDown.root", "READ"};
        pileupDownHist =
            dynamic_cast<TH1D*>(systDownFile->Get("pileup")->Clone());
    }
    else
    {
        // Make pileupReweighting stuff here
        dataPileupFile = new TFile{"pileup/2016/truePileupTest.root", "READ"};
        dataPU = (TH1D*)(dataPileupFile->Get("pileup")->Clone());
        mcPileupFile = new TFile{"pileup/2016/pileupMC.root", "READ"};
        mcPU = (TH1D*)(mcPileupFile->Get("pileup")->Clone());

        // Get systematic files too.
        systUpFile = new TFile{"pileup/2016/truePileupUp.root", "READ"};
        pileupUpHist = (TH1D*)(systUpFile->Get("pileup")->Clone());
        systDownFile = new TFile{"pileup/2016/truePileupDown.root", "READ"};
        pileupDownHist = (TH1D*)(systDownFile->Get("pileup")->Clone());
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

    // Initialise PDFs
    if ((systToRun & 1024 || systToRun & 2048) && is2016_)
    {
        LHAPDF::initPDFSet(1, "NNPDF30_nlo_nf_5_pdfas.LHgrid");
        //    LHAPDF::initPDFSet(1, "cteq6ll.LHpdf");
        //    LHAPDF::initPDFSet(1, "cteq6lg.LHgrid");
    }
}

void AnalysisAlgo::setupCuts()
{
    // Make cuts object. The methods in it should perhaps just be i nthe
    // AnalysisEvent class....
    cutObj = new Cuts{plots, plots, invertLepCut, is2016_};

    try
    {
        cutObj->parse_config(cutConfName);
    }
    catch (const std::exception)
    {
        std::cerr << "ERROR parsing cut configuration file" << std::endl;
        throw;
    }

    // For studying some trigger things. Default is false.
    cutObj->setSkipTrig(skipTrig);
    if (customJetRegion)
    {
        cutObj->setJetRegion(
            jetRegVars[0], jetRegVars[1], jetRegVars[2], jetRegVars[3]);
    }
    cutObj->setMetCut(metCut);
    cutObj->setMWCut(mwCut);
    cutObj->setMZCut(mzCut);
    if (doZplusCR_)
    {
        cutObj->setZplusControlRegionFlag(true);
    }
}

void AnalysisAlgo::setupPlots()
{
    // Do a little initialisation for the plots here. Will later on be done in a
    // config file. Initialise plot stage names.
    stageNames.emplace_back(std::make_pair("lepSel", "Lepton Cuts"));
    stageNames.emplace_back(std::make_pair("zMass", "Z Mass Cuts"));
    stageNames.emplace_back(std::make_pair("jetSel", "Jet Cuts"));
    stageNames.emplace_back(std::make_pair("bTag", "b-tag Cuts"));
    stageNames.emplace_back(std::make_pair("wMass", "W Mass Cuts"));
}

void AnalysisAlgo::runMainAnalysis()
{
    TMVA::gConfig().SetDrawProgressBar(true);

    if (totalLumi == 0.)
    {
        totalLumi = usePreLumi;
    }
    std::cout << "Using lumi: " << totalLumi << std::endl;

    bool datasetFilled{false};

    const std::string postLepSelSkimDir{
        std::string{"/scratch2/data/TopPhysics/postLepSelSkims"}
        + (is2016_ ? "2016" : "2017") + "/"};

    // Begin to loop over all datasets
    for (auto dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    {
        datasetFilled = false;
        TChain* datasetChain{new TChain{dataset->treeName().c_str()}};
        datasetChain->SetAutoSave(0);
        unsigned channelIndMax{256};

        const std::hash<std::string> hasher;
        srand(hasher(dataset->name()));

        channelIndMax = 64;
        for (unsigned channelInd{1}; channelInd != channelIndMax;
             channelInd = channelInd << 1)
        {
            if (!(channelInd & channelsToRun) && channelsToRun)
            {
                continue;
            }

            std::string chanName = channelSetup(channelInd);

            if (dataset->isMC() && skipMC)
            {
                continue;
            }
            if (!dataset->isMC() && skipData)
            {
                continue;
            }

            if (plots)
            { // Initialise a load of stuff that's required by the plotting
              // macro.

                // Gather all variables for plotting to make it easier to follow
                std::string histoName{dataset->getFillHisto()},
                    plotLabel{dataset->getPlotLabel()},
                    plotType{dataset->getPlotType()};
                int plotColour{dataset->getColour()};

                int systMask{1};
                for (unsigned systInd{0}; systInd < systNames.size(); systInd++)
                {
                    if (systInd > 0 && !(systToRun & systMask))
                    {
                        systMask = systMask << 1;
                        continue;
                    }
                    if (cutFlowMap.find(histoName + systNames[systInd])
                        == cutFlowMap.end())
                    {
                        const size_t numCutFlowBins{stageNames.size()};
                        cutFlowMap[histoName] = new TH1D{
                            (histoName + systNames[systInd] + "cutFlow")
                                .c_str(),
                            (histoName + systNames[systInd] + "cutFlow")
                                .c_str(),
                            boost::numeric_cast<int>(numCutFlowBins),
                            0,
                            boost::numeric_cast<double>(numCutFlowBins)};
                        if (systInd == 0
                            && datasetInfos.find(histoName)
                                   == datasetInfos.end())
                        {
                            legOrder.emplace_back(histoName);
                            plotOrder.emplace_back(histoName);
                            datasetInfos[histoName] = datasetInfo();
                            datasetInfos[histoName].colour = plotColour;
                            datasetInfos[histoName].legLabel = plotLabel;
                            datasetInfos[histoName].legType = plotType;
                        }
                        if (plots)
                        { // Only make all the plots if it's entirely necessary.
                            std::cout << "Made plots under " << histoName
                                      << " : " << systNames[systInd] + channel
                                      << std::endl;
                            if (plotsMap.find(channel) == plotsMap.end())
                            {
                                plotsVec.emplace_back(systNames[systInd]
                                                      + channel);
                            }
                            plotsMap[systNames[systInd] + channel]
                                    [(histoName)] = {};
                            for (unsigned j{0}; j < stageNames.size(); j++)
                            {
                                plotsMap[systNames[systInd] + channel]
                                        [histoName][stageNames[j].first] =
                                            std::make_shared<Plots>(
                                                plotTitles,
                                                plotNames,
                                                xMin,
                                                xMax,
                                                nBins,
                                                fillExp,
                                                xAxisLabels,
                                                cutStage,
                                                j,
                                                histoName + "_"
                                                    + stageNames[j].first
                                                    + systNames[systInd] + "_"
                                                    + channel);
                            }
                        }
                    } // end cutFlow find loop
                    if (systInd > 0)
                    {
                        systMask = systMask << 1;
                    }
                } // end systematic loop

            } // end plots if statement

            // If making plots and using saved histos, skip running over the
            // datasets ...
            if (plots && useHistos)
            {
                continue;
            }

            // If making either plots, make cut flow object.
            std::cerr << "Processing dataset " << dataset->name() << std::endl;
            if (!usePostLepTree)
            {
                if (!datasetFilled)
                {
                    if (!dataset->fillChain(datasetChain, numFiles))
                    {
                        std::cerr
                            << "There was a problem constructing the chain for "
                            << dataset->name() << " made of " << numFiles
                            << " files. Continuing with next dataset.\n";
                        continue;
                    }
                    datasetFilled = true;
                }
            }
            else
            {
                std::string inputPostfix{};
                inputPostfix += postfix;
                if (invertLepCut)
                {
                    inputPostfix += "invLep";
                }
                if (doNPLs_ && dataset->getPlotLabel() == "NPL")
                {
                    inputPostfix +=
                        "invLep"; // If plotting non-prompt leptons for this
                                  // dataset, be sure to read in the same sign
                                  // lepton post lepton skim!
                    // If making plots which include non-prompt leptons, then
                    // when running over "non-prompt" samples (as determined by
                    // their label), set the cutClass object NPL flag to true
                    // and invert charge seletion criteria, i.e. choose same
                    // sign leptons
                    cutObj->setNplFlag(true);
                    cutObj->setInvLepCut(true);
                }
                else if (doNPLs_ && dataset->getPlotLabel() != "NPL")
                {
                    cutObj->setNplFlag(false);
                    cutObj->setInvLepCut(false);
                }
                std::cout << postLepSelSkimDir + dataset->name() + inputPostfix
                                 + "SmallSkim.root"
                          << std::endl;
                datasetChain->Add((postLepSelSkimDir + dataset->name()
                                   + inputPostfix + "SmallSkim.root")
                                      .c_str());
            }

            cutObj->setMC(dataset->isMC());
            cutObj->setTriggerFlag(dataset->getTriggerFlag());
            std::cout << "Trigger flag: " << dataset->getTriggerFlag()
                      << std::endl;

            // Here we will initialise the b-tag eff plots if we are doing b-tag
            // efficiencies
            std::vector<TH2D*> bTagEffPlots;
            std::vector<std::string> denomNum{"Denom", "Num"};
            std::vector<std::string> typesOfEff{"b", "c", "uds", "g"};
            if (makePostLepTree && dataset->isMC())
            {
                int ptBins{4};
                int etaBins{4};
                float ptMin{0};
                float ptMax{200};
                float etaMin{0};
                float etaMax{2.4};
                for (unsigned denNum{0}; denNum < denomNum.size(); denNum++)
                {
                    for (unsigned type{0}; type < typesOfEff.size(); type++)
                    {
                        bTagEffPlots.emplace_back(
                            new TH2D{("bTagEff_" + denomNum[denNum] + "_"
                                      + typesOfEff[type])
                                         .c_str(),
                                     ("bTagEff_" + denomNum[denNum] + "_"
                                      + typesOfEff[type])
                                         .c_str(),
                                     ptBins,
                                     ptMin,
                                     ptMax,
                                     etaBins,
                                     etaMin,
                                     etaMax});
                    }
                }
                cutObj->setBTagPlots(bTagEffPlots, true);
            } // end btag eff plots.
            if (usePostLepTree && usebTagWeight && dataset->isMC())
            {
                // Get efficiency plots from the file. Will have to be from
                // post-lep sel trees I guess.
                std::string inputPostfix{};
                inputPostfix += postfix;
                if (doNPLs_ && dataset->getPlotLabel() == "NPL")
                {
                    inputPostfix +=
                        "invLep"; // If plotting non-prompt leptons for this
                                  // dataset, be sure to read in the same sign
                                  // lepton post lepton $
                }
                if (invertLepCut)
                {
                    inputPostfix += "invLep";
                }
                TFile* datasetFileForHists;
                datasetFileForHists =
                    new TFile((postLepSelSkimDir + dataset->name()
                               + inputPostfix + "SmallSkim.root")
                                  .c_str(),
                              "READ");
                for (unsigned denNum{0}; denNum < denomNum.size(); denNum++)
                {
                    for (unsigned eff{0}; eff < typesOfEff.size(); eff++)
                    {
                        bTagEffPlots.emplace_back(dynamic_cast<TH2D*>(
                            datasetFileForHists
                                ->Get(("bTagEff_" + denomNum[denNum] + "_"
                                       + typesOfEff[eff])
                                          .c_str())
                                ->Clone()));
                    }
                }
                for (unsigned plotIt{0}; plotIt < bTagEffPlots.size(); plotIt++)
                {
                    bTagEffPlots[plotIt]->SetDirectory(nullptr);
                }
                cutObj->setBTagPlots(bTagEffPlots, false);
                datasetFileForHists->Close();
            }

            // Here we will initialise the generator level weight histograms
            TH1I* generatorWeightPlot{nullptr};
            if (dataset->isMC())
            {
                if (usePostLepTree)
                {
                    std::string inputPostfix{};
                    inputPostfix += postfix;
                    if (invertLepCut)
                    {
                        inputPostfix += "invLep";
                    }
                    if (doNPLs_ && dataset->getPlotLabel() == "NPL")
                    {
                        inputPostfix +=
                            "invLep"; // If plotting non-prompt leptons for this
                                      // dataset, be sure to read in the same
                                      // sign lepton post lepton $
                    }
                    TFile* datasetFileForHists;
                    datasetFileForHists =
                        new TFile((postLepSelSkimDir + dataset->name()
                                   + inputPostfix + "SmallSkim.root")
                                      .c_str(),
                                  "READ");
                    generatorWeightPlot = dynamic_cast<TH1I*>(
                        datasetFileForHists->Get("sumNumPosMinusNegWeights")
                            ->Clone());
                    generatorWeightPlot->SetDirectory(nullptr);
                    datasetFileForHists->Close();
                }
                else
                {
                    generatorWeightPlot = dynamic_cast<TH1I*>(
                        dataset->getGeneratorWeightHistogram(numFiles)
                            ->Clone());
                    generatorWeightPlot->SetDirectory(nullptr);
                }
            }

            // extract the dataset weight. MC =
            // (lumi*crossSection)/(totalEvents), data = 1.0
            float datasetWeight{dataset->getDatasetWeight(totalLumi)};

            std::cout << datasetChain->GetEntries()
                      << " number of items in tree. Dataset weight: "
                      << datasetWeight << std::endl;
            if (datasetChain->GetEntries() == 0)
            {
                std::cout << "No entries in tree, skipping..." << std::endl;
                continue;
            }
            AnalysisEvent event{dataset->isMC(), datasetChain, is2016_};

            // Adding in some stuff here to make a skim file out of post lep sel
            // stuff
            TFile* outFile1{nullptr};
            TTree* cloneTree{nullptr};

            // If we're making the post lepton selection trees, set them up
            // here.
            if (makePostLepTree)
            {
                std::string invPostFix;
                if (invertLepCut)
                {
                    invPostFix = "invLep";
                }

                outFile1 = new TFile{(postLepSelSkimDir + dataset->name()
                                      + postfix + invPostFix + "SmallSkim.root")
                                         .c_str(),
                                     "RECREATE"};
                cloneTree = datasetChain->CloneTree(0);
                cloneTree->SetDirectory(outFile1);
                cutObj->setCloneTree(cloneTree);
            }

            // If we're making the MVA tree, set it up here.
            TFile* mvaOutFile{nullptr};
            std::vector<TTree*> mvaTree;
            // Add a few variables into the MVA tree for easy access of stuff
            // like lepton index etc
            double eventWeight{0.};
            int zLep1Index{
                -1}; // Addresses in elePF2PATWhatever of the z lepton
            int zLep2Index{-1};
            int wQuark1Index{-1};
            int wQuark2Index{-1};
            std::vector<size_t> jetInd{}; // The index of the selected jets;
            std::vector<size_t> bJetInd{}; // Index of selected b-jets;
            std::vector<double> muonMomentumSF{};
            float jetSmearValue[15]{};
            int isMC{dataset->isMC()}; // isMC flag for debug purposes
            event.isMC_ = (dataset->isMC());
            // Now add in the branches:

            if (makeMVATree)
            {
                boost::filesystem::create_directories(mvaDir);
                std::string invPostFix{};
                if (invertLepCut)
                {
                    invPostFix = "invLep";
                }
                mvaOutFile = new TFile{(mvaDir + dataset->name() + postfix
                                        + (invertLepCut ? invPostFix : "")
                                        + "mvaOut.root")
                                           .c_str(),
                                       "RECREATE"};
                mvaOutFile->SetCompressionSettings(ROOT::CompressionSettings(ROOT::kLZ4, 4));
                if (!mvaOutFile->IsOpen())
                {
                    throw std::runtime_error(
                        "MVA Tree TFile could not be opened!");
                }
                int systMask{1};
                // std::cout << "Making systematic trees for " <<
                // dataset->name() << ": ";
                for (unsigned systIn{0}; systIn < systNames.size(); systIn++)
                {
                    // std::cout << systNames[systIn] << " ";
                    //  	std::cout << "Making systs: " << systMask << " " <<
                    //  systToRun << " " << systIn << " " << (systMask &
                    //  systToRun) << std::endl;
                    /*  	if (systIn > 0 && !(systMask & systToRun)){
                      if (systIn > 0) systMask = systMask << 1;
                      continue;
                      }*/
                    mvaTree.emplace_back(datasetChain->CloneTree(0));
                    mvaTree[systIn]->SetDirectory(mvaOutFile);
                    mvaTree[systIn]->SetName(
                        (mvaTree[systIn]->GetName() + systNames[systIn])
                            .c_str());
                    mvaTree[systIn]->Branch(
                        "eventWeight", &eventWeight, "eventWeight/F");
                    mvaTree[systIn]->Branch(
                        "zLep1Index", &zLep1Index, "zLep1Index/I");
                    mvaTree[systIn]->Branch(
                        "zLep2Index", &zLep2Index, "zLep2Index/I");
                    mvaTree[systIn]->Branch(
                        "wQuark1Index", &wQuark1Index, "wQuark1Index/I");
                    mvaTree[systIn]->Branch(
                        "wQuark2Index", &wQuark2Index, "wQuark2Index/I");
                    mvaTree[systIn]->Branch("jetInd", &jetInd);
                    mvaTree[systIn]->Branch("muonMomentumSF", &muonMomentumSF);
                    mvaTree[systIn]->Branch(
                        "jetSmearValue", &jetSmearValue, "jetSmearValue[15]/F");
                    mvaTree[systIn]->Branch("bJetInd", &bJetInd);
                    mvaTree[systIn]->Branch("isMC", &isMC, "isMC/I");
                    if (systIn > 0)
                    {
                        systMask = systMask << 1;
                    }
                }
                std::cout << std::endl;
            }

            long long numberOfEvents{datasetChain->GetEntries()};
            if (nEvents && nEvents < numberOfEvents)
            {
                numberOfEvents = nEvents;
            }
            //    datasetChain->Draw("numElePF2PAT","numMuonPF2PAT > 2");
            //    TH1F * htemp = (TH1F*)gPad->GetPrimitive("htemp");
            //    htemp->SaveAs("tempCanvas.png");
            int foundEvents{0};
            double foundEventsNorm{0.0};

            // If event is amc@nlo, need to sum number of positive and negative
            // weights first.
            if (dataset->isMC())
            {
                // Load in plots
                sumPositiveWeights_ = dataset->getTotalEvents();
                sumNegativeWeights_ = generatorWeightPlot->GetBinContent(4);
                sumNegativeWeightsScaleUp_ = generatorWeightPlot->GetBinContent(
                    7); // Systematic Scale up
                sumNegativeWeightsScaleDown_ =
                    generatorWeightPlot->GetBinContent(
                        1); // Systematic Scale down
                if (sumNegativeWeights_ > sumPositiveWeights_)
                {
                    std::cout << "Something SERIOUSLY went wrong here - the "
                                 "number of postitive weights minus negative "
                                 "ones is greater than their sum?!"
                              << std::endl;
                    std::cout << "number of postitive weights: "
                              << sumPositiveWeights_ << std::endl;
                    std::cout
                        << "number of negative weights: " << sumNegativeWeights_
                        << std::endl;
                    exit(999);
                }
            }

            TMVA::Timer* lEventTimer{
                new TMVA::Timer{boost::numeric_cast<int>(numberOfEvents),
                                "Running over dataset ...",
                                false}};
            lEventTimer->DrawProgressBar(0, "");
            for (int i{0}; i < numberOfEvents; i++)
            {
                std::stringstream lSStrFoundEvents;
                lSStrFoundEvents << foundEvents;
                lEventTimer->DrawProgressBar(
                    i, ("Found " + lSStrFoundEvents.str() + " events."));
                event.GetEntry(i);
                // Do the systematics indicated by the systematic flag, oooor
                // just do data if that's your thing. Whatevs.
                int systMask{1};
                for (unsigned systInd{0}; systInd < systNames.size(); systInd++)
                {
                    if (!dataset->isMC() && systInd > 0)
                    {
                        break;
                    }
                    //	std::cout << systInd << " " << systMask << std::endl;
                    if (systInd > 0 && !(systMask & systToRun))
                    {
                        if (systInd > 0)
                        {
                            systMask = systMask << 1;
                        }
                        continue;
                    }
                    eventWeight = 1;

                    // apply generator weights here.
                    double generatorWeight{1.0};
                    if (dataset->isMC() && sumNegativeWeights_ >= 0)
                    {
                        if (systMask == 4096)
                        {
                            generatorWeight =
                                (sumPositiveWeights_)
                                / (sumNegativeWeightsScaleUp_)
                                * (event.weight_muF2muR2
                                   / std::abs(event.origWeightForNorm));
                        }
                        else if (systMask == 8192)
                        {
                            generatorWeight =
                                (sumPositiveWeights_)
                                / (sumNegativeWeightsScaleDown_)
                                * (event.weight_muF0p5muR0p5
                                   / std::abs(event.origWeightForNorm));
                        }
                        else
                        {
                            generatorWeight =
                                (sumPositiveWeights_) / (sumNegativeWeights_)
                                * (event.origWeightForNorm
                                   / std::abs(event.origWeightForNorm));
                        }
                        //	    	      std::cout << std::setprecision(5) <<
                        // std::fixed; 	                std::cout <<
                        // sumPositiveWeights_ << "/" << sumNegativeWeights_ <<
                        // "*" << event.origWeightForNorm
                        //<< "/" << std::abs(event.origWeightForNorm) <<
                        // std::endl; 	                std::cout << "generator
                        // level SF = " << generatorWeight << std::endl;
                        // std::cout << "NB. This should only not be 1.0 for
                        // aMC@NLO." << std::endl;
                    }
                    eventWeight *= generatorWeight;
                    // apply pileup weights here.
                    if (dataset->isMC())
                    { // no weights applied for synchronisation
                        double pileupWeight{puReweight->GetBinContent(
                            puReweight->GetXaxis()->FindBin(event.numVert))};
                        if (systMask == 64)
                        {
                            pileupWeight = puSystUp->GetBinContent(
                                puSystUp->GetXaxis()->FindBin(event.numVert));
                        }
                        if (systMask == 128)
                        {
                            pileupWeight = puSystDown->GetBinContent(
                                puSystDown->GetXaxis()->FindBin(event.numVert));
                        }
                        eventWeight *= pileupWeight;
                        // std::cout << "pileupWeight: " <<  pileupWeight <<
                        // std::endl;
                    }

                    // Scale according to lumi
                    eventWeight *= datasetWeight;

                    // apply negative weighting for SameSign MC lepton samples
                    // so that further downstream
                    if (dataset->isMC() && invertLepCut && !plots)
                    {
                        eventWeight *=
                            -1.0; // Should NOT be done when plotting
                                  // non-prompts - separate code for that
                    }

                    // Apply in cutClass, as the RATIO weight of OS/SS
                    // non-prompts cannot be applied before charge cuts (Z cand
                    // cuts) are applied If NPLs shape (for plotting purposes)
                    // apply OS/SS ratio SF
                    // if ( plots && doNpls_ && dataset->getPlotLabel() == "NPL"
                    // && !trileptonChannel_ ) { if ( channel == "ee" )
                    // eventWeight *= 1.24806; if ( channel == "mumu" )
                    // eventWeight *= 1.03226; if ( dataset->isMC() )
                    // eventWeight
                    // *= -1.0;
                    //}

                    // If ttbar, do reweight
                    //          std::cout << "eventWeight: " << eventWeight <<
                    //          std::endl;
                    if (dataset->name() == "ttbarInclusivePowerheg"
                        || dataset->name()
                               == "ttbarInclusivePowerheg_colourFlip"
                        || dataset->name() == "ttbarInclusivePowerheg_hdampUP"
                        || dataset->name() == "ttbarInclusivePowerheg_hdampDown"
                        || dataset->name() == "ttbarInclusivePowerheg_fsrup"
                        || dataset->name() == "ttbarInclusivePowerheg_fsrdown"
                        || dataset->name() == "ttbarInclusivePowerheg_isrup"
                        || dataset->name() == "ttbarInclusivePowerheg_isrdown"
                        || dataset->name() == "ttbar_2l2v"
                        || dataset->name() == "ttbar_hadronic"
                        || dataset->name() == "ttbar_semileptonic")
                    {
                        eventWeight *= event.topPtReweight;
                    }
                    //	  std::cout << "event.topPtReweight: " <<
                    // event.topPtReweight << std::endl;
                    //          std::cout << "eventWeight: " << eventWeight <<
                    //          std::endl;

                    //	  std::cout << "channel: " << channel << std::endl;
                    std::string histoName{dataset->getFillHisto()};

                    if (!cutObj->makeCuts(
                            event,
                            eventWeight,
                            plotsMap[systNames[systInd] + channel][histoName],
                            *cutFlowMap[histoName + systNames[systInd]],
                            systInd ? systMask : systInd))
                    {
                        if (systInd)
                        {
                            systMask = systMask << 1;
                        }
                        continue;
                    }

                    // Do Run 1 style PDF reweighting things for tW samples as
                    // they use Powerheg V1 Everything else uses LHE event
                    // weights
                    if (systMask == 1024 || systMask == 2048)
                    {
                        if (is2016_
                            && (dataset->name() == "tWInclusive"
                                || dataset->name() == "tbarWInclusive"
                                || dataset->name() == "tWInclusive_scaleup"
                                || dataset->name() == "tWInclusive_scaledown"
                                || dataset->name() == "tbarWInclusive_scaleup"
                                || dataset->name()
                                       == "tbarWInclusive_scaledown"))
                        {
                            // std::cout << std::setprecision(15) << eventWeight
                            // << " ";
                            LHAPDF::usePDFMember(1, 0);
                            float q{event.genPDFScale};
                            float x1{event.genPDFx1};
                            float x2{event.genPDFx2};
                            int id1{event.genPDFf1};
                            int id2{event.genPDFf2};
                            if (id2 == 21)
                            {
                                id2 = 0;
                            }
                            if (id1 == 21)
                            {
                                id1 = 0;
                            }
                            double xpdf1{LHAPDF::xfx(1, x1, q, id1)};
                            double xpdf2{LHAPDF::xfx(1, x2, q, id2)};
                            std::vector<float> pdf_weights;
                            // std::cout << q << " " << x1 << " " << x2 << " "
                            // << id1 << " " << id2 << " "; std::cout << xpdf1
                            // << " " << xpdf2 << " " << xpdf1 * xpdf2 << " ";
                            float min{1};
                            float max{1};
                            float pdfWeightUp{0};
                            float pdfWeightDown{0};
                            for (int j{1}; j <= 100; j++)
                            {
                                LHAPDF::usePDFMember(1, j);
                                double xpdf1_new{LHAPDF::xfx(1, x1, q, id1)};
                                double xpdf2_new{LHAPDF::xfx(1, x2, q, id2)};
                                // std::cout << " " << x1 << " " << id1 << " "
                                // << x2 << " " << id2 << " " << q << " "
                                // <<xpdf1
                                // << " " << xpdf2 << " " << xpdf1_new << " " <<
                                // xpdf2_new << " ";
                                double weight{1};
                                if ((xpdf1 * xpdf2) > 0.00001)
                                {
                                    weight =
                                        xpdf1_new * xpdf2_new / (xpdf1 * xpdf2);
                                }
                                pdf_weights.emplace_back(weight);
                                if (weight > 1.0)
                                {
                                    pdfWeightUp += (1 - weight) * (1 - weight);
                                }
                                if (weight < 1.0)
                                {
                                    pdfWeightDown +=
                                        (1 - weight) * (1 - weight);
                                }
                                if (weight > max)
                                {
                                    max = weight;
                                }
                                if (weight < min)
                                {
                                    min = weight;
                                }
                                //	      std::cout << " " << xpdf1_new << " "
                                //<< xpdf2_new << " " << weight << " ";
                            }
                            if (systMask == 1024)
                            {
                                eventWeight *= max;
                            }
                            if (systMask == 2048)
                            {
                                eventWeight *= min;
                            }
                            // std::cout << eventWeight << std::setprecision(4)
                            // << max << " " << min << " " <<
                            // 1+std::sqrt(pdfWeightUp) << " " <<
                            // 1-std::sqrt(pdfWeightDown) << std::endl;
                            // std::cout
                            // << std::setprecision(9) << " " << min << " " <<
                            // max << " " << eventWeight << std::endl;
                        }
                        // LHE event weights for everything else
                        else
                        {
                            if (systMask == 1024)
                            {
                                eventWeight *= event.weight_pdfMax; // Max
                            }
                            if (systMask == 2048)
                            {
                                eventWeight *= event.weight_pdfMin; // Min
                            }
                        }
                    }
                    if (systMask == 16384 || systMask == 32768)
                    {
                        if (systMask == 16384)
                        {
                            eventWeight *=
                                event.weight_alphaMin; // Max, but incorrectly
                                                       // named branch
                        }
                        if (systMask == 32768)
                        {
                            eventWeight *=
                                event.weight_alphaMax; // Min, but incorrectly
                                                       // named branch
                        }
                    }

                    // PSWeights
                    if (systMask == 65536)
                    {
                        eventWeight *= event.isrDefLo;
                    }
                    if (systMask == 131072)
                    {
                        eventWeight *= event.isrDefHi;
                    }
                    if (systMask == 262144)
                    {
                        eventWeight *= event.fsrDefLo;
                    }
                    if (systMask == 524288)
                    {
                        eventWeight *= event.fsrDefHi;
                    }

                    // Do the Zpt reweighting here
                    if (makeMVATree)
                    {
                        jetInd.clear();
                        bJetInd.clear();
                        zLep1Index = event.zPairIndex.first;
                        zLep2Index = event.zPairIndex.second;
                        muonMomentumSF = event.muonMomentumSF;
                        wQuark1Index = event.wPairIndex.first;
                        wQuark2Index = event.wPairIndex.second;
                        for (unsigned jetIndexIt{0}; jetIndexIt < 15;
                             jetIndexIt++)
                        {
                            if (jetIndexIt < event.jetIndex.size())
                            {
                                jetInd.emplace_back(event.jetIndex[jetIndexIt]);
                                jetSmearValue[jetIndexIt] =
                                    event.jetSmearValue[jetIndexIt];
                            }
                            else
                            {
                                jetSmearValue[jetIndexIt] = 0.0;
                            }
                        }
                        for (unsigned bJetIt{0}; bJetIt < 10; bJetIt++)
                        {
                            if (bJetIt < event.bTagIndex.size())
                            {
                                bJetInd.emplace_back(event.bTagIndex[bJetIt]);
                            }
                        }
                        mvaTree[systInd]->Fill();
                    }

                    foundEvents++;
                    foundEventsNorm += eventWeight;
                    if (systInd > 0)
                    {
                        systMask = systMask << 1;
                    }
                } // End systematics loop.
            } // end event loop

            // If we're making post lepSel skims save the tree here
            if (makePostLepTree)
            {
                outFile1->cd();
                std::cout << "\nPrinting some info on the tree "
                          << dataset->name() << " " << cloneTree->GetEntries()
                          << std::endl;
                std::cout << "But there were :" << datasetChain->GetEntries()
                          << " entries in the original tree" << std::endl;
                cloneTree->Write();
                // Write out mc generator level info
                if (dataset->isMC())
                {
                    generatorWeightPlot->Write();
                }
                for (unsigned i{0}; i < bTagEffPlots.size(); i++)
                {
                    bTagEffPlots[i]->Write();
                }

                delete cloneTree;
                cloneTree = nullptr;
                outFile1->Write();
                outFile1->Close();
                outFile1 = nullptr;
            }

            // Save mva outputs
            if (makeMVATree)
            {
                std::string invPostFix{};
                if (invertLepCut)
                {
                    invPostFix = "invLep";
                }

                std::cout << (mvaDir + dataset->name() + postfix
                              + (invertLepCut ? invPostFix : "")
                              + "mvaOut.root")
                          << std::endl;
                mvaOutFile->cd();
                std::cout << std::endl;
                int systMask{1};
                std::cout << "Saving Systematics: ";
                for (unsigned systInd{0}; systInd < systNames.size(); systInd++)
                {
                    if (systInd > 0 && !(systToRun & systMask))
                    {
                        systMask = systMask << 1;
                        continue;
                    }
                    std::cout << systNames[systInd] << ": "
                              << mvaTree[systInd]->GetEntriesFast() << " "
                              << std::flush;
                    mvaTree[systInd]->FlushBaskets();
                    if (systInd > 0)
                    {
                        systMask = systMask << 1;
                    }
                    if (!dataset->isMC())
                    {
                        break;
                    }
                }
                std::cout << std::endl;
                // Save the efficiency plots for b-tagging here if we're doing
                // that.
                if (makePostLepTree)
                {
                    for (unsigned i{0}; i < bTagEffPlots.size(); i++)
                    {
                        bTagEffPlots[i]->Write();
                    }
                }
                mvaOutFile->Write();
                for (unsigned i{0}; i < mvaTree.size(); i++)
                {
                    delete mvaTree[i];
                }
                mvaOutFile->Close();
            }
            std::cerr << "\nFound " << foundEvents << " in " << dataset->name()
                      << std::endl;
            std::cerr << "Found " << foundEventsNorm
                      << " after normalisation in " << dataset->name()
                      << std::endl;
            std::cerr << "\n\n";
            // Delete generator level plot. Avoid memory leaks, kids.
            delete generatorWeightPlot;
            generatorWeightPlot = nullptr;
            // Delete plots from out btag vector. Avoid memory leaks, kids.
            if (makePostLepTree)
            {
                for (unsigned i{0}; i < bTagEffPlots.size(); i++)
                {
                    delete bTagEffPlots[i];
                }
            }

            // datasetChain->MakeClass("AnalysisEvent");
        } // end channel loop.
        delete datasetChain;
    } // end dataset loop
}

void AnalysisAlgo::savePlots()
{
    // Save all plot objects. For testing purposes.

    // Now test out the histogram plotter class I just wrote.
    // Make the plotting object.
    if (plots)
    {
        HistogramPlotter plotObj =
            HistogramPlotter(legOrder, plotOrder, datasetInfos, is2016_);

        // If either making or reading in histos, then set the correct read in
        // directory
        if ((makeHistos || useHistos) && plots)
        {
            plotObj.setHistogramFolder(histoDir);
        }

        // If making histos, save the output!
        if (makeHistos && plots)
        {
            std::cout << "Saving histograms for later use ..." << std::endl;
            for (unsigned i{0}; i < plotsVec.size(); i++)
            {
                plotObj.saveHistos(plotsMap[plotsVec[i]]);
            }
            plotObj.saveHistos(
                cutFlowMap,
                "cutFlow",
                channel); // Don't forget to save the cutflow too!
        }

        if (!makeHistos)
        {
            if (useHistos)
            {
                plotObj.loadHistos(); // If using saved histos, read them in ...
            }
            plotObj.setLabelOne("CMS Preliminary");
            plotObj.setLabelTwo("Some amount of lumi");
            plotObj.setPostfix("");
            plotObj.setOutputFolder(outFolder);

            for (unsigned i{0}; i < plotsVec.size(); i++)
            {
                std::cout << plotsVec[i] << std::endl;
                if (plots)
                {
                    plotObj.plotHistos(plotsMap[plotsVec[i]]);
                }
            }

            // cut flow x axis labels
            std::vector<std::string> cutFlowLabels;
            for (std::vector<
                     std::pair<std::string, std::string>>::const_iterator lIt =
                     stageNames.begin();
                 lIt != stageNames.end();
                 ++lIt)
            {
                cutFlowLabels.emplace_back((*lIt).second);
            }
            if (useHistos)
            {
                cutFlowMap = plotObj.loadCutFlowMap("cutFlow", channel);
            }
            plotObj.makePlot(
                cutFlowMap, "data/MC Yield", "cutFlow", cutFlowLabels);
        }
    }

    // Delete all the plot objects.

    std::cerr << "Gets to the delete bit" << std::endl;
    std::cerr << "But not past it" << std::endl;
}

std::string AnalysisAlgo::channelSetup(unsigned channelInd)
{
    std::string chanName{};

    if (channelsToRun)
    {
        if (channelInd & 5)
        { // ee channels
            cutObj->setNumLeps(0, 0, 2, 2);
            cutObj->setCutConfTrigLabel("e");
            channel = "ee";
            postfix = "ee";
            chanName += "ee";
        }
        if (channelInd & 10)
        { // mumu channels
            cutObj->setNumLeps(2, 2, 0, 0);
            cutObj->setCutConfTrigLabel("m");
            channel = "mumu";
            postfix = "mumu";
            chanName += "mumu";
        }
        if (channelInd & 3)
        { // nominal samples
            cutObj->setInvLepCut(false);
            invertLepCut = false;
            chanName += "nom";
        }
        if (channelInd & 12)
        { // same sign samples
            cutObj->setInvLepCut(true);
            invertLepCut = true;
            chanName += "inv";
        }
        if (channelInd & 16)
        { // emu channel for ttbar background estimation
            cutObj->setNumLeps(1, 1, 1, 1);
            cutObj->setCutConfTrigLabel("d");
            channel = "emu";
            postfix = "emu";
            chanName += "emu";
        }
        if (channelInd & 32)
        { // same signemu channel for NPL ttbar background estimation
            cutObj->setNumLeps(1, 1, 1, 1);
            cutObj->setCutConfTrigLabel("d");
            channel = "emu";
            postfix = "emu";
            cutObj->setInvLepCut(true);
            invertLepCut = true;
            chanName += "invemu";
        }
    }
    return chanName;
}
