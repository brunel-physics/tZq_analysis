#include "debugInfo.hpp"

#include "AnalysisEvent.hpp"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLatex.h"
#include "TMVA/Timer.h"
#include "TTree.h"
#include "config_parser.hpp"

#include <boost/program_options.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

int main(int argc, char* argv[])
{
    DebugInfo debugInfo;

    debugInfo.parseCommandLineArguements(argc, argv);
    debugInfo.runMainAnalysis();
}

DebugInfo::DebugInfo()
{
}

DebugInfo::~DebugInfo()
{
}

void DebugInfo::parseCommandLineArguements(int argc, char* argv[])
{
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()("help,h", "Print this message.")(
        "config,c",
        po::value<std::string>(&config)->required(),
        "The configuration file to be used.")(
        ",n",
        po::value<long>(&nEvents)->default_value(0),
        "The number of events to be run over. All if set to 0.")(
        "outFolder,o",
        po::value<std::string>(&outFolder)->default_value("plots/debug/"),
        "The output directory for the plots. Overrides the config file.")(
        "postfix,s",
        po::value<std::string>(&postfix)->default_value("default"),
        "Set postfix for plots. Overrides the config file.")(
        "2016", po::bool_switch(&is2016_), "Use 2016 conditions (SFs, et al.)")(
        "MC,m", po::bool_switch(&isMC_), "Running over Monte Carlo.")(
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
    std::cout << std::setprecision(3) << std::fixed;

    // Some vectors that will be filled in the parsing.
    totalLumi = 0;
    if (!Parser::parse_config(config, datasets, totalLumi))
    {
        std::cerr << "There was an error parsing the config file.\n";
        exit(0);
    }
}

void DebugInfo::runMainAnalysis()
{
    TH1F* histMuChannel =
        new TH1F("histMuChannel", "mu trigger yields vs Run;", 7, 0.5, 7.5);
    TH1F* histMuMuChannel =
        new TH1F("histMuMuChannel", "mu trigger yields vs Run;", 7, 0.5, 7.5);

    TH1F* histChannel =
        new TH1F("histChannel",
                 "channel yields vs RunH/RunsB-G;channel;RunH/RunsB-G",
                 2,
                 0.5,
                 2.5);
    TH1F* histTrileptonChannel = new TH1F(
        "histTrileptonChannel",
        "trilepton channel yields vs RunH/RunsB-G;channel;RunH/RunsB-G",
        4,
        0.5,
        4.5);
    TH1F* histNumJets =
        new TH1F("histNumJets",
                 "number of jets vs RunH/RunsB-G;# jets;RunH/RunsB-G",
                 10,
                 -0.5,
                 9.5);
    TH1F* histNumBJets =
        new TH1F("histNumBJets",
                 "number of b-jets vs RunH/RunsB-G;# b-jets;RunH/RunsB-G",
                 5,
                 -0.5,
                 4.5);

    std::pair<int, int> numElectrons{0, 0};
    std::pair<int, int> numMuons{0, 0};
    std::pair<int, int> numEEE{0, 0};
    std::pair<int, int> numEEMU{0, 0};
    std::pair<int, int> numEMUMU{0, 0};
    std::pair<int, int> numMUMUMU{0, 0};
    std::pair<int, int> numJets[10]{{0, 0}};
    std::pair<int, int> numBJets[5]{{0, 0}};

    bool datasetFilled = false;

    if (totalLumi == 0.)
    {
        totalLumi = usePreLumi;
    }
    std::cout << "Using lumi: " << totalLumi << std::endl;
    for (auto dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    {
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

        AnalysisEvent* event = new AnalysisEvent(
            dataset->isMC(), dataset->getTriggerFlag(), datasetChain, is2016_);

        int numberOfEvents = datasetChain->GetEntries();
        if (nEvents && nEvents < numberOfEvents)
        {
            numberOfEvents = nEvents;
        }
        auto lEventTimer =
            new TMVA::Timer(numberOfEvents, "Running over dataset ...", false);
        lEventTimer->DrawProgressBar(0, "");
        if (numberOfEvents < 1)
        {
            continue;
        }
        for (int i = 0; i < numberOfEvents; i++)
        {
            lEventTimer->DrawProgressBar(i);
            event->GetEntry(i);

            // MET Triggers
            /*
                  bool metTrig = false;
                  if( event->HLT_MET250_v2 > 0 ) metTrig = true;
                  if( event->HLT_MET250_v3 > 0 ) metTrig = true;
                  if( event->HLT_MET250_v4 > 0 ) metTrig = true;
                  if( event->HLT_MET250_v5 > 0 ) metTrig = true;

                  if( event->HLT_PFHT300_PFMET100_v1 > 0 ) > 0 ) metTrig = true;
                  if( event->HLT_PFHT300_PFMET100_v2 > 0 ) metTrig = true;
                  if( event->HLT_PFHT300_PFMET100_v3 > 0 ) metTrig = true;
                  if( event->HLT_PFHT300_PFMET100_v4 > 0 ) metTrig = true;
                  if( event->HLT_PFHT300_PFMET110_v4 > 0 ) metTrig = true;
                  if( event->HLT_PFHT300_PFMET110_v5 > 0 ) metTrig = true;
                  if( event->HLT_PFHT300_PFMET110_v6 > 0 ) metTrig = true;

                  if( event->HLT_PFMET120_PFMHT120_IDTight_v2 > 0 ) metTrig =
               true; if( event->HLT_PFMET120_PFMHT120_IDTight_v3 > 0 ) metTrig =
               true; if( event->HLT_PFMET120_PFMHT120_IDTight_v4 > 0 ) metTrig =
               true; if( event->HLT_PFMET120_PFMHT120_IDTight_v5 > 0 ) metTrig =
               true; if( event->HLT_PFMET120_PFMHT120_IDTight_v6 > 0 ) metTrig =
               true; if( event->HLT_PFMET120_PFMHT120_IDTight_v7 > 0 ) metTrig =
               true; if( event->HLT_PFMET120_PFMHT120_IDTight_v8 > 0 ) metTrig =
               true; if( event->HLT_PFMET170_HBHECleaned_v2 > 0 ) metTrig =
               true; if( event->HLT_PFMET170_HBHECleaned_v3 > 0 ) metTrig =
               true; if( event->HLT_PFMET170_HBHECleaned_v4 > 0 ) metTrig =
               true; if( event->HLT_PFMET170_HBHECleaned_v5 > 0 ) metTrig =
               true; if( event->HLT_PFMET170_HBHECleaned_v6 > 0 ) metTrig =
               true; if( event->HLT_PFMET170_HBHECleaned_v7 > 0 ) metTrig =
               true; if( event->HLT_PFMET170_HBHECleaned_v8 > 0 ) metTrig =
               true; if( event->HLT_PFMET170_HBHECleaned_v9 > 0 ) metTrig =
               true;

                  if ( !metTrig ) continue;
            */

            // Muon Triggers
            bool muTrig{false};

            if (event->HLT_IsoMu24_v1 > 0)
            {
                muTrig = true;
            }
            if (event->HLT_IsoMu24_v2 > 0)
            {
                muTrig = true;
            }
            if (event->HLT_IsoMu24_v3 > 0)
            {
                muTrig = true;
            }
            if (event->HLT_IsoMu24_v4 > 0)
            {
                muTrig = true;
            }
            if (event->HLT_IsoTkMu24_v1 > 0)
            {
                muTrig = true;
            }
            if (event->HLT_IsoTkMu24_v2 > 0)
            {
                muTrig = true;
            }
            if (event->HLT_IsoTkMu24_v3 > 0)
            {
                muTrig = true;
            }
            if (event->HLT_IsoTkMu24_v4 > 0)
            {
                muTrig = true;
            }

            // Double Muon Trigger
            bool mumuTrig{false};

            // non-DZ legs are prescaled for Run2016H
            if (event->eventRun < 280919 && !isMC_)
            {
                if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2 > 0)
                {
                    mumuTrig = true;
                }
                if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3 > 0)
                {
                    mumuTrig = true;
                }
                if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4 > 0)
                {
                    mumuTrig = true;
                }
                if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6 > 0)
                {
                    mumuTrig = true;
                }

                if (event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2 > 0)
                {
                    mumuTrig = true;
                }
                if (event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3 > 0)
                {
                    mumuTrig = true;
                }
                if (event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5 > 0)
                {
                    mumuTrig = true;
                }
            }
            // non-DZ legs in MC
            else if (isMC_)
            {
                if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6 > 0)
                {
                    mumuTrig = true;
                }
                if (event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5 > 0)
                {
                    mumuTrig = true;
                }
            }

            // DZ legs avaliable all the time but inefficient in data for Runs
            // B-F -> hence uses of non-DZ legs
            if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0)
            {
                mumuTrig = true;
            }
            if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 > 0)
            {
                mumuTrig = true;
            }
            if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0)
            {
                mumuTrig = true;
            }
            if (event->HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 > 0)
            {
                mumuTrig = true;
            }
            if (event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0)
            {
                mumuTrig = true;
            }
            if (event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0)
            {
                mumuTrig = true;
            }
            if (event->HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6 > 0)
            {
                mumuTrig = true;
            }

            int runBin{-1};
            if (event->eventRun >= 272007 && event->eventRun < 275657)
            {
                runBin = 1; // Run B
            }
            if (event->eventRun >= 275657 && event->eventRun < 276315)
            {
                runBin = 2; // Run C
            }
            if (event->eventRun >= 276315 && event->eventRun < 276831)
            {
                runBin = 3; // Run D
            }
            if (event->eventRun >= 276831 && event->eventRun < 277772)
            {
                runBin = 4; // Run E
            }
            if (event->eventRun >= 277772 && event->eventRun < 278820)
            {
                runBin = 5; // Run F
            }
            if (event->eventRun >= 278820 && event->eventRun < 280919)
            {
                runBin = 6; // Run G
            }
            if (event->eventRun >= 280919)
            {
                runBin = 7; // Run H
            }

            double lumi = 0.0;
            if (runBin == 1)
            {
                lumi = 5784.596;
            }
            if (runBin == 2)
            {
                lumi = 2573.399;
            }
            if (runBin == 3)
            {
                lumi = 4248.384;
            }
            if (runBin == 4)
            {
                lumi = 4009.132;
            }
            if (runBin == 5)
            {
                lumi = 3101.618;
            }
            if (runBin == 6)
            {
                std::cout << "runBin: " << runBin << std::endl;
                lumi = 7540.488;
            }
            if (runBin == 7)
            {
                std::cout << "runBin: " << runBin << std::endl;
                lumi = 8605.690;
            }

            histMuChannel->Fill(runBin, muTrig * lumi);
            histMuMuChannel->Fill(runBin, mumuTrig * lumi);

            if (event->numElePF2PAT == 2 && event->numMuonPF2PAT == 0)
            {
                if (event->elePF2PATPT[0] < 15.0)
                {
                    continue;
                }
                if (event->elePF2PATPT[1] < 15.0)
                {
                    continue;
                }
                if (event->elePF2PATCharge[0] * event->elePF2PATCharge[1] >= 0)
                {
                    continue; // check electron pair have correct charge.
                }
                TLorentzVector lepton1{event->elePF2PATGsfPx[0],
                                       event->elePF2PATGsfPy[0],
                                       event->elePF2PATGsfPz[0],
                                       event->elePF2PATGsfE[0]};
                TLorentzVector lepton2{event->elePF2PATGsfPx[1],
                                       event->elePF2PATGsfPy[1],
                                       event->elePF2PATGsfPz[1],
                                       event->elePF2PATGsfE[1]};
                double invMass{(lepton1 + lepton2).M() - 91.1};
                if (std::abs(invMass) > 30.0)
                {
                    continue;
                }
                if (!isMC_)
                {
                    if (event->eventRun <= 280385)
                    {
                        numElectrons.first += 1; // If Runs B-G
                    }
                    else
                    {
                        numElectrons.second += 1; // else if Run H
                    }
                }
                else
                {
                    numElectrons.first += 1; // just MC
                }
            }

            if (event->numMuonPF2PAT == 2 && event->numElePF2PAT == 0)
            {
                if (event->muonPF2PATPt[0] < 15.0)
                {
                    continue;
                }
                if (event->muonPF2PATPt[1] < 15.0)
                {
                    continue;
                }
                if (event->muonPF2PATCharge[0] * event->muonPF2PATCharge[1]
                    >= 0)
                {
                    continue;
                }
                TLorentzVector lepton1{event->muonPF2PATPX[0],
                                       event->muonPF2PATPY[0],
                                       event->muonPF2PATPZ[0],
                                       event->muonPF2PATE[0]};
                TLorentzVector lepton2{event->muonPF2PATPX[1],
                                       event->muonPF2PATPY[1],
                                       event->muonPF2PATPZ[1],
                                       event->muonPF2PATE[1]};
                double invMass{(lepton1 + lepton2).M() - 91.1};
                if (std::abs(invMass) > 30.0)
                {
                    continue;
                }
                if (!isMC_)
                {
                    if (event->eventRun <= 280385)
                    {
                        numMuons.first += 1; // If Runs B-G
                    }
                    else
                    {
                        numMuons.second += 1; // else if Run H
                    }
                }
                else
                {
                    numMuons.first += 1; // just MC
                }
            }
            /*
                  if ( event->numElePF2PAT == 3 && event->numMuonPF2PAT == 0 ) {
            //        if ( event->elePF2PATPT[0] < 15.0 ) continue;
            //        if ( event->elePF2PATPT[1] < 15.0 ) continue;
            //        if ( event->elePF2PATPT[2] < 15.0 ) continue;
                    if (!isMC_){
                      if ( event->eventRun <= 280385 ) numEEE.first += 1; // If
            Runs B-G else numEEE.second += 1; // else if Run H
                    }
                    else numEEE.first += 1; // just MC
                  }
                  if ( event->numElePF2PAT == 2 && event->numMuonPF2PAT == 1 ) {
            //        if ( event->elePF2PATPT[0] < 15.0 ) continue;
            //        if ( event->elePF2PATPT[1] < 15.0 ) continue;
            //        if ( event->muonPF2PATPt[0] < 15.0 ) continue;
            //        if ( event->elePF2PATCharge[0] * event->elePF2PATCharge[1]
            >= 0 )  continue; // check electron pair have correct charge. if
            (!isMC_){ if ( event->eventRun <= 280385 ) numEEMU.first += 1; // If
            Runs B-G else numEEMU.second += 1; // else if Run H
                    }
                    else numEEMU.first += 1; // just MC
                  }
                  if ( event->numElePF2PAT == 1 && event->numMuonPF2PAT == 2 ) {
            //        if ( event->elePF2PATPT[0] < 15.0 ) continue;
            //        if ( event->muonPF2PATPt[0] < 15.0 ) continue;
            //        if ( event->muonPF2PATPt[1] < 15.0 ) continue;
            //        if ( event->muonPF2PATCharge[0] *
            event->muonPF2PATCharge[1] >= 0 ) continue; if (!isMC_){ if (
            event->eventRun <= 280385 ) numEMUMU.first += 1; // If Runs B-G else
            numEMUMU.second += 1; // else if Run H
                    }
                    else numEMUMU.first += 1; // just MC
                  }
                  if ( event->numElePF2PAT == 0 && event->numMuonPF2PAT == 3 ) {
            //        if ( event->muonPF2PATPt[0] < 15.0 ) continue;
            //        if ( event->muonPF2PATPt[1] < 15.0 ) continue;
            //        if ( event->muonPF2PATPt[2] < 15.0 ) continue;
                    if (!isMC_){
                      if ( event->eventRun <= 280385 ) numMUMUMU.first += 1; //
            If Runs B-G else numMUMUMU.second += 1; // else if Run H
                    }
                    else numMUMUMU.first += 1; // just MC
                  }
            */
            if (event->numJetPF2PAT < 0)
            {
                continue;
            }
            if (event->numJetPF2PAT < 10)
            {
                if (!isMC_)
                {
                    if (event->eventRun <= 280385)
                    {
                        numJets[event->numJetPF2PAT].first += 1;
                    }
                    else
                    {
                        numJets[event->numJetPF2PAT].second += 1;
                    }
                }
                else
                {
                    numJets[event->numJetPF2PAT].first += 1; // just MC
                }
            }
            unsigned int bJets{0};
            if (event->numJetPF2PAT < 1)
            {
                continue;
            }
            for (int j = 0; j < event->numJetPF2PAT; j++)
            {
                if (event->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                        [i]
                    > 0.5426)
                {
                    bJets += 1;
                }
            }
            if (bJets < 5)
            {
                if (!isMC_)
                {
                    if (event->eventRun <= 280385)
                    {
                        numBJets[bJets].first += 1;
                    }
                    else
                    {
                        numBJets[bJets].second += 1;
                    }
                }
                else
                {
                    numBJets[bJets].first += 1; // just MC
                }
            }
        }
        delete datasetChain;
    } // end dataset loop

    //  std::cout << "numElectrons.first/numElectrons.second : " <<
    //  numElectrons.first << "/" << numElectrons.second << std::endl; std::cout
    //  << "numMuons.first/numMuons.second : " << numMuons.first << "/" <<
    //  numMuons.second << std::endl;

    float electronFraction =
        float(numElectrons.second) / (float(numElectrons.first) + 1.0e-06);
    float muonFraction =
        float(numMuons.second) / (float(numMuons.first) + 1.0e-06);

    if (!isMC_)
    {
        histChannel->Fill(1, electronFraction);
    }
    else
    {
        histChannel->Fill(1, numElectrons.first);
    }
    if (!isMC_)
    {
        histChannel->Fill(2, muonFraction);
    }
    else
    {
        histChannel->Fill(2, numMuons.first);
    }

    float eeeFraction = float(numEEE.second) / (float(numEEE.first) + 1.0e-06);
    float eemuFraction =
        float(numEEMU.second) / (float(numEEMU.first) + 1.0e-06);
    float emumuFraction =
        float(numEMUMU.second) / (float(numEMUMU.first) + 1.0e-06);
    float mumumuFraction =
        float(numMUMUMU.second) / (float(numMUMUMU.first) + 1.0e-06);

    if (!isMC_)
    {
        histTrileptonChannel->Fill(1, eeeFraction);
    }
    else
    {
        histTrileptonChannel->Fill(1, numEEE.first);
    }
    if (!isMC_)
    {
        histTrileptonChannel->Fill(2, eemuFraction);
    }
    else
    {
        histTrileptonChannel->Fill(2, numEEMU.first);
    }
    if (!isMC_)
    {
        histTrileptonChannel->Fill(3, emumuFraction);
    }
    else
    {
        histTrileptonChannel->Fill(3, numEMUMU.first);
    }
    if (!isMC_)
    {
        histTrileptonChannel->Fill(4, mumumuFraction);
    }
    else
    {
        histTrileptonChannel->Fill(4, numMUMUMU.first);
    }

    for (int i = 0; i < 10; i++)
    {
        if (!isMC_)
        {
            histNumJets->Fill(i,
                              numJets[i].second / (numJets[i].first + 1.0e-06));
        }
        else
        {
            histNumJets->Fill(i, numJets[i].first);
        }
    }
    for (int i = 0; i < 5; i++)
    {
        if (!isMC_)
        {
            histNumBJets->Fill(
                i, numBJets[i].second / (numBJets[i].first + 1.0e-06));
        }
        else
        {
            histNumBJets->Fill(i, numBJets[i].first);
        }
    }

    histMuChannel->GetXaxis()->SetBinLabel(1, "B");
    histMuChannel->GetXaxis()->SetBinLabel(2, "C");
    histMuChannel->GetXaxis()->SetBinLabel(3, "D");
    histMuChannel->GetXaxis()->SetBinLabel(4, "E");
    histMuChannel->GetXaxis()->SetBinLabel(5, "F");
    histMuChannel->GetXaxis()->SetBinLabel(6, "G");
    histMuChannel->GetXaxis()->SetBinLabel(7, "H");

    histMuMuChannel->GetXaxis()->SetBinLabel(1, "B");
    histMuMuChannel->GetXaxis()->SetBinLabel(2, "C");
    histMuMuChannel->GetXaxis()->SetBinLabel(3, "D");
    histMuMuChannel->GetXaxis()->SetBinLabel(4, "E");
    histMuMuChannel->GetXaxis()->SetBinLabel(5, "F");
    histMuMuChannel->GetXaxis()->SetBinLabel(6, "G");
    histMuMuChannel->GetXaxis()->SetBinLabel(7, "H");

    histChannel->GetXaxis()->SetBinLabel(1, "ee");
    histChannel->GetXaxis()->SetBinLabel(2, "#mu#mu");

    histTrileptonChannel->GetXaxis()->SetBinLabel(1, "eee");
    histTrileptonChannel->GetXaxis()->SetBinLabel(2, "ee#mu");
    histTrileptonChannel->GetXaxis()->SetBinLabel(3, "e#mu#mu");
    histTrileptonChannel->GetXaxis()->SetBinLabel(4, "#mu#mu#mu");

    histChannel->SetMinimum(0.0);
    histNumJets->SetMinimum(0.0);
    histNumBJets->SetMinimum(0.0);

    mkdir((outFolder).c_str(), 0700);
    TFile* outFile =
        new TFile((outFolder + postfix + ".root").c_str(), "RECREATE");

    histMuChannel->Write();
    histMuMuChannel->Write();
    histChannel->Write();
    histTrileptonChannel->Write();
    histNumJets->Write();
    histNumBJets->Write();
    outFile->Close();
}
