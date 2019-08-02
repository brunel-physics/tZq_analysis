#include "MvaEvent.hpp"
#include "TLorentzVector.h"
#include "TMVA/Config.h"
#include "TMVA/Timer.h"
#include "TTree.h"
#include "config_parser.hpp"
#include "makeMVAinputAlgo.hpp"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <limits>
#include <memory>

MakeMvaInputs::MakeMvaInputs()
    : inputVars{}
    , oldMetFlag{false}
    , is2016{false}
    , era{}
    , ttbarControlRegion{false}
    , useSidebandRegion{false}
    , doMC{false}
    , doSysts{false}
    , doData{false}
    , doFakes{false}
    , inputDir{"mvaTest/"}
    , outputDir{"mvaInputs/"}
{
}

MakeMvaInputs::~MakeMvaInputs()
{
}

void MakeMvaInputs::parseCommandLineArguements(const int argc, char* argv[])
{
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()("help,h", "Print this message.")(
        "met,m", po::bool_switch(&oldMetFlag), "Use old MET uncerts recipe")(
        "ttbar", po::bool_switch(&ttbarControlRegion), "Make ttbar CR stuff")(
        "2016",
        po::bool_switch(&is2016),
        "Process 2016 data (default to 2017)")(
        "inputDir,i",
        po::value<std::string>(&inputDir),
        "Mva skims input directory")("outputDir,o",
                                     po::value<std::string>(&outputDir),
                                     "Mva inputs output directory")(
        "sideband,s",
        po::bool_switch(&useSidebandRegion),
        "Make side band CR plots")(
        "data,D", po::bool_switch(&doData), "Run the data analysis")(
        "systs,S",
        po::bool_switch(&doSysts),
        "Run dedicated systematic analysis")(
        "MC,M", po::bool_switch(&doMC), "Run MC analysis")(
        "fakes,F", po::bool_switch(&doFakes), "Run fakes analysis");

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

    catch (const std::logic_error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "Use -h or --help for help." << std::endl;
        std::exit(1);
    }
}

void MakeMvaInputs::runMainAnalysis()
{
    TMVA::gConfig().SetDrawProgressBar(true);

    era = is2016 ? "2016" : "2017";

    const auto getListOfMCs{[this]() -> std::map<std::string, std::string> {
        if (is2016)
        { // 2016
            return {{"ttHTobb", "TTHbb"},
                    {"ttHToNonbb", "TTHnonbb"},
                    {"WWW", "WWW"},
                    {"WWZ", "WWZ"},
                    {"WZZ", "WZZ"},
                    {"ZZZ", "ZZZ"},
                    {"WW1l1nu2q", "WW1l1nu2q"},
                    {"WW2l2nu", "WW2l2nu"},
                    {"ZZ4l", "ZZ4l"},
                    {"ZZ2l2nu", "ZZ2l2nu"},
                    {"ZZ2l2q", "ZZ2l2q"},
                    {"WZjets", "WZ3l1nu"},
                    {"WZ2l2q", "WZ2l2q"},
                    {"WZ1l1nu2q", "WZ1l1nu2q"},
                    {"sChannel", "TsChan"},
                    {"tChannel", "TtChan"},
                    {"tbarChannel", "TbartChan"},
                    {"tWInclusive", "TW"},
                    {"tbarWInclusive", "TbarW"},
                    {"tZq", "TZQ"},
                    {"tHq", "THQ"},
                    {"ttWlnu", "TTWlnu"},
                    {"ttW2q", "TTW2q"},
                    {"ttZ2l2nu", "TTZ2l2nu"},
                    {"ttZ2q", "TTZ2q"},
                    {"ttbarInclusivePowerheg", "TT"},
                    {"tWZ", "TWZ"},
                    {"wPlusJets", "Wjets"},
                    {"DYJetsToLL_Pt-0To50", "DYJetsLLPt0To50"},
                    {"DYJetsToLL_Pt-50To100", "DYJetsLLPt50To100"},
                    {"DYJetsToLL_Pt-100To250", "DYJetsLLPt100To250"},
                    {"DYJetsToLL_Pt-250To400", "DYJetsLLPt250To400"},
                    {"DYJetsToLL_Pt-400To650", "DYJetsLLPt400To650"},
                    {"DYJetsToLL_Pt-650ToInf", "DYJetsLLPt650ToInf"}};
        }
        else
        { // 2017
            return {{"ttH_bb", "TTHbb"},
                    {"ttH_nonbb", "TTHnonbb"},
                    {"WWW", "WWW"},
                    {"WWZ", "WWZ"},
                    {"WZZ", "WZZ"},
                    {"ZZZ", "ZZZ"},
                    {"WW_1l1nu2q", "WW1l1nu2q"},
                    {"WW_2l2nu", "WW2l2nu"},
                    {"ZZ_4l", "ZZ4l"},
                    {"ZZ_2l2nu", "ZZ2l2nu"},
                    {"ZZ_2l2q", "ZZ2l2q"},
                    {"WZ_3l1nu", "WZ3l3nu"},
                    {"WZ_2l2q", "WZ2l2q"},
                    {"WZ_1l1nu2q", "WZ1l1nu2q"},
                    {"WG_lnug", "WG11nu1g"},
                    {"ZG_llg", "ZG2l1g"},
                    {"t_s_channel", "TsChan"},
                    {"t_t_channel", "TtChan"},
                    {"tbar_t_channel", "TbartChan"},
                    {"tW", "TW"},
                    {"tbarW", "TbarW"},
                    {"tZq", "TZQ"},
                    {"tHq", "THQ"},
                    {"ttW_lnu", "TTWlnu"},
                    {"ttW_2q", "TTW2q"},
                    {"ttZ_2l2nu", "TTZ2l2nu"},
                    {"ttZ_2q", "TTZ2q"},
                    {"ttgamma", "TTG"},
                    {"ttbar_2l2v", "TT2l2v"},
                    {"ttbar_hadronic", "TTjets"},
                    {"ttbar_semileptonic", "TT2l2q"},
                    {"tWZ", "TWZ"},
                    {"wPlusJets", "Wjets"},
                    {"DYJetsToLL_M-10to50", "DYJetsToLLM10to50"},
                    {"DYJetsToLL_M-50", "DYJetsToLLM50"}};
        }
    }};

    static const auto listOfMCs = getListOfMCs();
    const auto channels{[=]() -> std::vector<std::string> {
        if (ttbarControlRegion)
        {
            return {"emu"};
        }
        else
        {
            return {"ee", "mumu"};
        }
    }()};

    const std::vector<std::string> systs = {"",
                                            "__trig__plus",
                                            "__trig__minus",
                                            "__jer__plus",
                                            "__jer__minus",
                                            "__jes__plus",
                                            "__jes__minus",
                                            "__pileup__plus",
                                            "__pileup__minus",
                                            "__bTag__plus",
                                            "__bTag__minus",
                                            "__met__plus",
                                            "__met__minus",
                                            "__pdf__plus",
                                            "__pdf__minus",
                                            "__ME__plus",
                                            "__ME__minus"};

    if (doMC)
    {
        standardAnalysis(listOfMCs, systs, channels, useSidebandRegion);
    }
    if (doSysts)
    {
        standardAnalysis({{"tChannel_scaleUp", "TtChan__scaleUp"},
                          {"tChannel_scaleDown", "TtChan__scaleDown"},
                          {"tChannel_hdampUp", "TtChan__hdampUp"},
                          {"tChannel_hdampDown", "TtChan__hdampDown"},
                          {"tbarChannel_scaleUp", "TbartChan__scaleUp"},
                          {"tbarChannel_scaleDown", "TbartChan__scaleDown"},
                          {"tbarChannel_hdampUp", "TbartChan__hdampUp"},
                          {"tbarChannel_hdampDown", "TbartChan__hdampDown"},
                          {"ttbarInclusivePowerheg_hdampUp", "TT__hdampUp"},
                          {"ttbarInclusivePowerheg_hdampDown", "TT__hdampDown"},
                          {"ttbarInclusivePowerheg_fsrUp", "TT__fsrUp"},
                          {"ttbarInclusivePowerheg_fsrDown", "TT__fsrDown"},
                          {"ttbarInclusivePowerheg_isrUp", "TT__isrUp"},
                          {"ttbarInclusivePowerheg_isrDown", "TT__isrDown"},
                          {"tWInclusive_scaleUp", "TtW__scaleUp"},
                          {"tWInclusive_scaleDown", "TtW__scaleDown"},
                          {"tbarWInclusive_scaleUp", "TbartW__scaleUp"},
                          {"tbarWInclusive_scaleDown", "TbartW__scaleDown"},
                          {"tZq_scaleUp", "tZq__scaleUp"},
                          {"tZq_scaleDown", "tZq__scaleDown"}},
                         {""},
                         channels,
                         useSidebandRegion);
    }
    if (doData)
    {
        dataAnalysis(channels, useSidebandRegion);
    }
    if (doFakes)
    {
        sameSignAnalysis(listOfMCs, channels, useSidebandRegion);
    }
}

void MakeMvaInputs::standardAnalysis(
    const std::map<std::string, std::string>& listOfMCs,
    const std::vector<std::string>& systs,
    const std::vector<std::string>& channels,
    const bool useSidebandRegion)
{
    std::string treeNamePostfixSig{""};
    std::string treeNamePostfixSB{""};
    if (useSidebandRegion)
    {
        std::cout << "Using control region stuff" << std::endl;
        treeNamePostfixSig = "sig_";
        treeNamePostfixSB = "ctrl_";
    }

    // loop over nominal samples
    for (const auto& mc : listOfMCs)
    {
        const std::string sample{mc.first};
        const std::string outSample{mc.second};

        const auto longest_string{[](std::vector<std::string> v) {
            return std::max_element(v.begin(),
                                    v.end(),
                                    [](std::string a, std::string b) {
                                        return a.size() < b.size();
                                    })
                ->size();
        }};
        boost::format systFormat{
            "%-" + (std::to_string(longest_string(channels))) + "s    %-"
            + std::to_string(longest_string(systs))
            + "s    %12.2f %+8.2f %+10.2f%%"};

        std::cout << "Doing " << sample << " : " << std::endl;

        auto outFile{new TFile{
            (outputDir + "histofile_" + listOfMCs.at(sample) + ".root").c_str(),
            "RECREATE"}};

        // loop over systematics
        std::unordered_map<std::string, long double> nominalEvents{};
        for (const auto& syst : systs)
        {
            auto outTreeSig{new TTree{
                ("Ttree_" + treeNamePostfixSig + outSample + syst).c_str(),
                ("Ttree_" + treeNamePostfixSig + outSample + syst).c_str()}};
            TTree* outTreeSdBnd{};
            setupBranches(outTreeSig);

            if (useSidebandRegion)
            {
                outTreeSdBnd = new TTree{
                    ("Ttree_" + treeNamePostfixSB + outSample + syst).c_str(),
                    ("Ttree_" + treeNamePostfixSB + outSample + syst).c_str()};
                setupBranches(outTreeSdBnd);
            }

            // loop over channels
            for (const auto& channel : channels)
            {
                auto inFile{new TFile{
                    (inputDir + sample + channel + "mvaOut.root").c_str(),
                    "READ"}};
                TTree* tree;
                if (syst == "__met__plus" || syst == "__met__minus")
                {
                    tree = dynamic_cast<TTree*>(inFile->Get("tree"));
                }
                else
                {
                    tree = dynamic_cast<TTree*>(
                        inFile->Get(("tree" + syst).c_str()));
                }

                //        TChain* tree;
                //        if ( syst == "__met__plus" || syst ==
                //        "__met__minus" ) tree = new TChain("tree"); else
                //        tree = new TChain(("tree"+syst).c_str());
                //        tree->Add((inputDir+sample+channel+"mvaOut.root").c_str());
                const long long numberOfEvents{tree->GetEntries()};
                auto event{new MvaEvent{true, tree, true}};

                // loop over events
                long double nEvents{0};
                for (long long i{0}; i < numberOfEvents; i++)
                {
                    event->GetEntry(i);

                    fillTree(outTreeSig,
                             outTreeSdBnd,
                             event,
                             outSample + syst,
                             channel,
                             false);

                    nEvents += event->eventWeight;
                } // end event loop

                if (syst.empty())
                {
                    nominalEvents.emplace(channel, nEvents);
                }
                std::cout << systFormat % channel % syst % nEvents
                                 % (nEvents - nominalEvents[channel])
                                 % (((nEvents - nominalEvents[channel])
                                    / nominalEvents[channel]) * 100)
                          << std::endl;

                inFile->Close();
            } // end channel loop
            outFile->cd();
            outTreeSig->Write();
            delete outTreeSig;
            delete outTreeSdBnd;
            if (useSidebandRegion)
            {
                outTreeSdBnd->Write();
            }
        } // end systematic loop
        outFile->Write();
        outFile->Close();
    } // end sample loop
}

void MakeMvaInputs::dataAnalysis(const std::vector<std::string>& channels,
                                 const bool useSidebandRegion)
{
    const std::unordered_map<std::string, std::string> outChanToData = {
        {"ee", "DataEG"}, {"mumu", "DataMu"}, {"emu", "MuonEG"}};

    std::string treeNamePostfixSig{""};
    std::string treeNamePostfixSB{""};
    if (useSidebandRegion)
    {
        std::cout << "Using control region stuff" << std::endl;
        treeNamePostfixSig = "sig_";
        treeNamePostfixSB = "ctrl_";
    }

    for (const auto& channel : channels)
    {
        std::cout << "Data " << channel << std::endl;
        const std::string outChan{outChanToData.at(channel)};

        auto outTreeSig{
            new TTree{("Ttree_" + treeNamePostfixSig + outChan).c_str(),
                      ("Ttree_" + treeNamePostfixSig + outChan).c_str()}};
        setupBranches(outTreeSig);
        TTree* outTreeSdBnd{};
        if (useSidebandRegion)
        {
            outTreeSdBnd =
                new TTree{("Ttree_" + treeNamePostfixSB + outChan).c_str(),
                          ("Ttree_" + treeNamePostfixSB + outChan).c_str()};
            setupBranches(outTreeSdBnd);
        }
        TFile outFile{(outputDir + "histofile_" + outChan + ".root").c_str(),
                      "RECREATE"};
        TChain dataChain{"tree"};
        dataChain.Add(
            (inputDir + channel + "Run" + era + channel + "mvaOut.root")
                .c_str());

        MvaEvent event{false, &dataChain, true};
        const long long numberOfEvents{dataChain.GetEntries()};
        TMVA::Timer lEventTimer{boost::numeric_cast<int>(numberOfEvents),
                                "Running over dataset ...",
                                false};
        for (long long i{0}; i < numberOfEvents; i++)
        {
            lEventTimer.DrawProgressBar(i);
            event.GetEntry(i);
            fillTree(outTreeSig, outTreeSdBnd, &event, outChan, channel, false);
        }
        outFile.cd();
        outFile.Write();
        outTreeSig->Write();
        if (useSidebandRegion)
        {
            outTreeSdBnd->Write();
        }
        outFile.Close();
    }
}

void MakeMvaInputs::sameSignAnalysis(
    const std::map<std::string, std::string>& listOfMCs,
    const std::vector<std::string>& channels,
    const bool useSidebandRegion)
{
    std::vector<std::string> outFakeChannels{"FakeEG", "FakeMu"};
    const std::unordered_map<std::string, std::string> outFakeChanToData{
        {"FakeEG", "ee"}, {"FakeMu", "mumu"}};
    const std::unordered_map<std::string, std::string> chanMap{
        {"ee", "eeRun" + era}, {"mumu", "mumuRun" + era}};

    std::string treeNamePostfixSig;
    std::string treeNamePostfixSB;
    if (useSidebandRegion)
    {
        std::cout << "Using control region stuff" << std::endl;
        treeNamePostfixSig = "sig_";
        treeNamePostfixSB = "ctrl_";
    }

    for (const auto& outChan : outFakeChannels)
    {
        // Get same sign data
        const std::string chan{outFakeChanToData.at(outChan)};

        TChain dataChain{"tree"};
        // Get expected real SS events from MC
        for (const auto& mc : listOfMCs)
        {
            const std::string sample{mc.first};

            // std::cout << "Doing SS fakes " << sample << std::endl;
            if (!dataChain.Add(
                    (inputDir + sample + chan + "invLepmvaOut.root").c_str()))
                abort();
        }

        if (!dataChain.AddFile(
                (inputDir + chanMap.at(chan) + chan + "invLepmvaOut.root")
                    .c_str()))
            abort();

        auto outFile{
            new TFile{(outputDir + "histofile_" + outChan + ".root").c_str(),
                      "RECREATE"}};
        auto outTreeSig{
            new TTree{("Ttree_" + treeNamePostfixSig + outChan).c_str(),
                      ("Ttree_" + treeNamePostfixSig + outChan).c_str()}};
        setupBranches(outTreeSig);
        TTree* outTreeSdBnd;
        if (useSidebandRegion)
        {
            outTreeSdBnd =
                new TTree{("Ttree_" + treeNamePostfixSB + outChan).c_str(),
                          ("Ttree_" + treeNamePostfixSB + outChan).c_str()};
            setupBranches(outTreeSdBnd);
        }
        MvaEvent event{false, &dataChain, true};

        const long long numberOfEvents{dataChain.GetEntries()};
        TMVA::Timer lEventTimer{boost::numeric_cast<int>(numberOfEvents),
                                "Running over dataset ...",
                                false};
        for (long long i{0}; i < numberOfEvents; i++)
        {
            lEventTimer.DrawProgressBar(i);
            event.GetEntry(i);
            fillTree(outTreeSig, outTreeSdBnd, &event, outChan, chan, true);
        } // end event loop

        outFile->cd();
        outFile->Write();
        outTreeSig->Write();
        if (useSidebandRegion)
        {
            outTreeSdBnd->Write();
        }
        outFile->Close();
    }
}

std::pair<TLorentzVector, TLorentzVector>
    MakeMvaInputs::sortOutLeptons(const MvaEvent* tree,
                                  const std::string& channel) const
{
    TLorentzVector zLep1;
    TLorentzVector zLep2;

    const int zlep1Index{tree->zLep1Index};
    const int zlep2Index{tree->zLep2Index};

    if (channel == "ee")
    {
        zLep1.SetPxPyPzE(tree->elePF2PATPX[zlep1Index],
                         tree->elePF2PATPY[zlep1Index],
                         tree->elePF2PATPZ[zlep1Index],
                         tree->elePF2PATE[zlep1Index]);
        zLep2.SetPxPyPzE(tree->elePF2PATPX[zlep2Index],
                         tree->elePF2PATPY[zlep2Index],
                         tree->elePF2PATPZ[zlep2Index],
                         tree->elePF2PATE[zlep2Index]);
    }
    if (channel == "mumu")
    {
        zLep1.SetPxPyPzE(tree->muonPF2PATPX[zlep1Index],
                         tree->muonPF2PATPY[zlep1Index],
                         tree->muonPF2PATPZ[zlep1Index],
                         tree->muonPF2PATE[zlep1Index]);
        zLep2.SetPxPyPzE(tree->muonPF2PATPX[zlep2Index],
                         tree->muonPF2PATPY[zlep2Index],
                         tree->muonPF2PATPZ[zlep2Index],
                         tree->muonPF2PATE[zlep2Index]);
    }
    if (channel == "emu")
    {
        zLep1.SetPxPyPzE(tree->elePF2PATPX[zlep1Index],
                         tree->elePF2PATPY[zlep1Index],
                         tree->elePF2PATPZ[zlep1Index],
                         tree->elePF2PATE[zlep1Index]);
        zLep2.SetPxPyPzE(tree->muonPF2PATPX[zlep2Index],
                         tree->muonPF2PATPY[zlep2Index],
                         tree->muonPF2PATPZ[zlep2Index],
                         tree->muonPF2PATE[zlep2Index]);
    }

    return {zLep1, zLep2};
}

std::pair<TLorentzVector, TLorentzVector>
    MakeMvaInputs::sortOutHadronicW(const MvaEvent* tree) const
{
    TLorentzVector wQuark1;
    TLorentzVector wQuark2;
    wQuark1.SetPxPyPzE(tree->jetPF2PATPx[tree->wQuark1Index],
                       tree->jetPF2PATPy[tree->wQuark1Index],
                       tree->jetPF2PATPz[tree->wQuark1Index],
                       tree->jetPF2PATE[tree->wQuark1Index]);
    wQuark2.SetPxPyPzE(tree->jetPF2PATPx[tree->wQuark2Index],
                       tree->jetPF2PATPy[tree->wQuark2Index],
                       tree->jetPF2PATPz[tree->wQuark2Index],
                       tree->jetPF2PATE[tree->wQuark2Index]);

    return {wQuark1, wQuark2};
}

std::pair<std::vector<int>, std::vector<TLorentzVector>> MakeMvaInputs::getJets(
    const MvaEvent* tree, const int syst, TLorentzVector met) const
{
    std::vector<int> jetList{};
    std::vector<TLorentzVector> jetVecList{};

    for (int i{0}; i != tree->NJETS; i++)
    {
        if (tree->jetInd[i] > -1)
        {
            jetList.emplace_back(tree->jetInd[i]);
            jetVecList.emplace_back(getJetVec(tree,
                                              tree->jetInd[i],
                                              tree->jetSmearValue[i],
                                              met,
                                              syst,
                                              true));
        }
        else
        {
            continue;
        }
    }

    return {jetList, jetVecList};
}

std::pair<std::vector<int>, std::vector<TLorentzVector>>
    MakeMvaInputs::getBjets(const MvaEvent* tree,
                            const int syst,
                            TLorentzVector met,
                            const std::vector<int>& jets) const
{
    std::vector<int> bJetList{};
    std::vector<TLorentzVector> bJetVecList{};

    for (int i{0}; i != tree->NBJETS; i++)
    {
        if (tree->bJetInd[i] > -1)
        {
            bJetList.emplace_back(tree->bJetInd[i]);
            bJetVecList.emplace_back(
                getJetVec(tree,
                          jets.at(tree->bJetInd[i]),
                          tree->jetSmearValue[tree->bJetInd[i]],
                          met,
                          syst,
                          false));
        }
        else
        {
            continue;
        }
    }

    return {bJetList, bJetVecList};
}

TLorentzVector MakeMvaInputs::getJetVec(const MvaEvent* tree,
                                        const int index,
                                        const float smearValue,
                                        TLorentzVector& metVec,
                                        const int syst,
                                        const bool doMetSmear) const
{
    TLorentzVector returnJet;
    returnJet.SetPxPyPzE(tree->jetPF2PATPx[index],
                         tree->jetPF2PATPy[index],
                         tree->jetPF2PATPz[index],
                         tree->jetPF2PATE[index]);
    returnJet *= smearValue;

    const static JetCorrectionUncertainty jetUnc{
        is2016 ? "scaleFactors/2016/"
                 "Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt"
               : "scaleFactors/2017/"
                 "Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt"};

    if (syst == 16)
    {
        returnJet *=
            1 + jetUnc.getUncertainty(returnJet.Pt(), returnJet.Eta(), 1);
    }
    else if (syst == 32)
    {
        returnJet *=
            1 + jetUnc.getUncertainty(returnJet.Pt(), returnJet.Eta(), 2);
    }

    if (doMetSmear && smearValue > 0.01)
    {
        metVec.SetPx(metVec.Px() + tree->jetPF2PATPx[index]);
        metVec.SetPy(metVec.Py() + tree->jetPF2PATPy[index]);

        metVec.SetPx(metVec.Px() - returnJet.Px());
        metVec.SetPy(metVec.Py() - returnJet.Py());
    }

    return returnJet;
}

TLorentzVector
    MakeMvaInputs::doUncMet(TLorentzVector met,
                            const TLorentzVector& zLep1,
                            const TLorentzVector& zLep2,
                            const std::vector<TLorentzVector>& jetVecs,
                            const unsigned syst) const
{
    double uncMetX{met.Px() + zLep1.Px() + zLep2.Px()};
    double uncMetY{met.Py() + zLep1.Py() + zLep2.Py()};

    for (const auto& jetVec : jetVecs)
    {
        uncMetX += jetVec.Px();
        uncMetY += jetVec.Py();
    }

    if (syst == 1024)
    {
        met.SetPx(met.Px() + 0.1 * uncMetX);
        met.SetPy(met.Py() + 0.1 * uncMetY);
    }

    else if (syst == 2048)
    {
        met.SetPx(met.Px() - 0.1 * uncMetX);
        met.SetPy(met.Py() - 0.1 * uncMetY);
    }

    return met;
}

void MakeMvaInputs::setupBranches(TTree* tree)
{
    tree->Branch("EvtWeight", &inputVars["eventWeight"], "EvtWeight/F");
    tree->Branch("EvtNumber", &inputVars["eventNumber"], "EvtNumber/F");
    tree->Branch("mTW", &inputVars["mTW"], "mTW/F");
    tree->Branch("wQuark1Pt", &inputVars["wQuark1Pt"], "wQuark1Pt/F");
    tree->Branch("wQuark1Eta", &inputVars["wQuark1Eta"], "wQuark1Eta/F");
    tree->Branch("wQuark1Phi", &inputVars["wQuark1Phi"], "wQuark1Phi/F");
    tree->Branch("wQuark2Pt", &inputVars["wQuark2Pt"], "wQuark2Pt/F");
    tree->Branch("wQuark2Eta", &inputVars["wQuark2Eta"], "wQuark2Eta/F");
    tree->Branch("wQuark2Phi", &inputVars["wQuark2Phi"], "wQuark2Phi/F");
    tree->Branch("wPairMass", &inputVars["wPairMass"], "wPairMass/F");
    tree->Branch("wPairPt", &inputVars["wPairPt"], "wPairPt/F");
    tree->Branch("wPairEta", &inputVars["wPairEta"], "wPairEta/F");
    tree->Branch("wPairPhi", &inputVars["wPairPhi"], "wPairPhi/F");
    tree->Branch("met", &inputVars["met"], "met/F");
    tree->Branch("nJets", &inputVars["nJets"], "nJets/F");
    tree->Branch("leadJetPt", &inputVars["leadJetPt"], "leadJetPt/F");
    tree->Branch("leadJetEta", &inputVars["leadJetEta"], "leadJetEta/F");
    tree->Branch("leadJetPhi", &inputVars["leadJetPhi"], "leadJetPhi/F");
    tree->Branch("leadJetbTag", &inputVars["leadJetbTag"], "leadJetbTag/F");
    tree->Branch("secJetPt", &inputVars["secJetPt"], "secJetPt/F");
    tree->Branch("secJetEta", &inputVars["secJetEta"], "secJetEta/F");
    tree->Branch("secJetPhi", &inputVars["secJetPhi"], "secJetPhi/F");
    tree->Branch("secJetbTag", &inputVars["secJetbTag"], "secJetbTag/F");
    tree->Branch("thirdJetPt", &inputVars["thirdJetPt"], "thirdJetPt/F");
    tree->Branch("thirdJetEta", &inputVars["thirdJetEta"], "thirdJetEta/F");
    tree->Branch("thirdJetPhi", &inputVars["thirdJetPhi"], "thirdJetPhi/F");
    tree->Branch("thirdJetbTag", &inputVars["thirdJetbTag"], "thirdJetbTag/F");
    tree->Branch("fourthJetPt", &inputVars["fourthJetPt"], "fourthJetPt/F");
    tree->Branch("fourthJetEta", &inputVars["fourthJetEta"], "fourthJetEta/F");
    tree->Branch("fourthJetPhi", &inputVars["fourthJetPhi"], "fourthJetPhi/F");
    tree->Branch(
        "fourthJetbTag", &inputVars["fourthJetbTag"], "fourthJetbTag/F");
    tree->Branch("nBjets", &inputVars["nBjets"], "nBjets/F");
    tree->Branch("bTagDisc", &inputVars["bTagDisc"], "bTagDisc/F");
    tree->Branch("lep1Pt", &inputVars["lep1Pt"], "lep1Pt/F");
    tree->Branch("lep1Eta", &inputVars["lep1Eta"], "lep1Eta/F");
    tree->Branch("lep1Phi", &inputVars["lep1Phi"], "lep1Phi/F");
    tree->Branch("lep1RelIso", &inputVars["lep1RelIso"], "lep1RelIso/F");
    tree->Branch("lep1D0", &inputVars["lep1D0"], "lep1D0/F");
    tree->Branch("lep2Pt", &inputVars["lep2Pt"], "lep2Pt/F");
    tree->Branch("lep2Eta", &inputVars["lep2Eta"], "lep2Eta/F");
    tree->Branch("lep2Phi", &inputVars["lep2Phi"], "lep2Phi/F");
    tree->Branch("lep2RelIso", &inputVars["lep2RelIso"], "lep2RelIso/F");
    tree->Branch("lep2D0", &inputVars["lep2D0"], "lep2D0/F");
    tree->Branch("lepMass", &inputVars["lepMass"], "lepMass/F");
    tree->Branch("lepPt", &inputVars["lepPt"], "lepPt/F");
    tree->Branch("lepEta", &inputVars["lepEta"], "lepEta/F");
    tree->Branch("lepPhi", &inputVars["lepPhi"], "lepPhi/F");
    tree->Branch("zMass", &inputVars["zMass"], "zMass/F");
    tree->Branch("zPt", &inputVars["zPt"], "zPt/F");
    tree->Branch("zEta", &inputVars["zEta"], "zEta/F");
    tree->Branch("zPhi", &inputVars["zPhi"], "zPhi/F");
    tree->Branch("topMass", &inputVars["topMass"], "topMass/F");
    tree->Branch("topPt", &inputVars["topPt"], "topPt/F");
    tree->Branch("topEta", &inputVars["topEta"], "topEta/F");
    tree->Branch("topPhi", &inputVars["topPhi"], "topPhi/F");
    tree->Branch("jjdelR", &inputVars["j1j2delR"], "jjdelR/F");
    tree->Branch("jjdelPhi", &inputVars["j1j2delPhi"], "jjdelPhi/F");
    tree->Branch("wwdelR", &inputVars["w1w2delR"], "wwdelR/F");
    tree->Branch("wwdelPhi", &inputVars["w1w2delPhi"], "wwdelPhi/F");
    tree->Branch("zLepdelR", &inputVars["zLepdelR"], "zLepdelR/F");
    tree->Branch("zLepdelPhi", &inputVars["zLepdelPhi"], "zLepdelPhi/F");
    tree->Branch(
        "zl1Quark1DelR", &inputVars["zl1Quark1DelR"], "zl1Quark1DelR/F");
    tree->Branch(
        "zl1Quark1DelPhi", &inputVars["zl1Quark1DelPhi"], "zl1Quark1DelPhi/F");
    tree->Branch(
        "zl1Quark2DelR", &inputVars["zl1Quark2DelR"], "zl1Quark2DelR/F");
    tree->Branch(
        "zl1Quark2DelPhi", &inputVars["zl1Quark2DelPhi"], "zl1Quark2DelPhi/F");
    tree->Branch(
        "zl2Quark1DelR", &inputVars["zl2Quark1DelR"], "zl2Quark1DelR/F");
    tree->Branch(
        "zl2Quark1DelPhi", &inputVars["zl2Quark1DelPhi"], "zl2Quark1DelPhi/F");
    tree->Branch(
        "zl2Quark2DelR", &inputVars["zl2Quark2DelR"], "zl2Quark2DelR/F");
    tree->Branch(
        "zl2Quark2DelPhi", &inputVars["zl2Quark2DelPhi"], "zl2Quark2DelPhi/F");
    tree->Branch("zlb1DelR", &inputVars["zlb1DelR"], "zlb1DelR/F");
    tree->Branch("zlb1DelPhi", &inputVars["zlb1DelPhi"], "zlb1DelPhi/F");
    tree->Branch("zlb2DelR", &inputVars["zlb2DelR"], "zlb2DelR/F");
    tree->Branch("zlb2DelPhi", &inputVars["zlb2DelPhi"], "zlb2DelPhi/F");
    tree->Branch("lepHt", &inputVars["lepHt"], "lepHt/F");
    tree->Branch("wQuarkHt", &inputVars["wQuarkHt"], "wQuarkHt/F");
    tree->Branch("totPt", &inputVars["totPt"], "totPt/F");
    tree->Branch("totEta", &inputVars["totEta"], "totEta/F");
    tree->Branch("totPhi", &inputVars["totPhi"], "totPhi/F");
    tree->Branch("totPtVec", &inputVars["totPtVec"], "totPtVec/F");
    tree->Branch("totVecM", &inputVars["totVecM"], "totVecM/F");
    tree->Branch("Channel", &inputVars["chan"], "Channel/F");
    tree->Branch("totPt2Jet", &inputVars["totPt2Jet"], "totPt2Jet/F");
    tree->Branch("wzdelR", &inputVars["wZdelR"], "wzdelR/F");
    tree->Branch("wzdelPhi", &inputVars["wZdelPhi"], "wzdelPhi/F");
    tree->Branch("zQuark1DelR", &inputVars["zQuark1DelR"], "zQuark1DelR/F");
    tree->Branch(
        "zQuark1DelPhi", &inputVars["zQuark1DelPhi"], "zQuark1DelPhi/F");
    tree->Branch("zQuark2DelR", &inputVars["zQuark2DelR"], "zQuark2DelR/F");
    tree->Branch(
        "zQuark2DelPhi", &inputVars["zQuark2DelPhi"], "zQuark2DelPhi/F");
    tree->Branch("zTopDelR", &inputVars["zTopDelR"], "zTopDelR/F");
    tree->Branch("zTopDelPhi", &inputVars["zTopDelPhi"], "zTopDelPhi/F");
    tree->Branch("zl1TopDelR", &inputVars["zl1TopDelR"], "zl1TopDelR/F");
    tree->Branch("zl1TopDelPhi", &inputVars["zl1TopDelPhi"], "zl1TopDelPhi/F");
    tree->Branch("zl2TopDelR", &inputVars["zl2TopDelR"], "zl2TopDelR/F");
    tree->Branch("zl2TopDelPhi", &inputVars["zl2TopDelPhi"], "zl2TopDelPhi/F");
    tree->Branch("wTopDelR", &inputVars["wTopDelR"], "wTopDelR/F");
    tree->Branch("wTopDelPhi", &inputVars["wTopDelPhi"], "wTopDelPhi/F");
    tree->Branch("w1TopDelR", &inputVars["w1TopDelR"], "w1TopDelR/F");
    tree->Branch("w1TopDelPhi", &inputVars["w1TopDelPhi"], "w1TopDelPhi/F");
    tree->Branch("w2TopDelR", &inputVars["w2TopDelR"], "w2TopDelR/F");
    tree->Branch("w2TopDelPhi", &inputVars["w2TopDelPhi"], "w2TopDelPhi/F");
    tree->Branch("zjminR", &inputVars["minZJetR"], "zjminR/F");
    tree->Branch("zjminPhi", &inputVars["minZJetPhi"], "zjminPhi/F");
    tree->Branch("totHt", &inputVars["totHt"], "totHt/F");
    tree->Branch("jetHt", &inputVars["jetHt"], "jetHt/F");
    tree->Branch("jetMass", &inputVars["jetMass"], "jetMass/F");
    tree->Branch("jetPt", &inputVars["jetPt"], "jetPt/F");
    tree->Branch("jetEta", &inputVars["jetEta"], "jetEta/F");
    tree->Branch("jetPhi", &inputVars["jetPhi"], "jetPhi/F");
    tree->Branch("jetMass3", &inputVars["jetMass3"], "jetMass3/F");
    tree->Branch("totHtOverPt", &inputVars["totHtOverPt"], "totHtOverPt/F");
    tree->Branch("chi2", &inputVars["chi2"], "chi2/F");
}

void MakeMvaInputs::fillTree(TTree* outTreeSig,
                             TTree* outTreeSdBnd,
                             MvaEvent* tree,
                             const std::string& label,
                             const std::string& channel,
                             const bool SameSignMC)
{
    unsigned syst{0};
    const double NaN{std::numeric_limits<double>::quiet_NaN()};

    if (label.find("__met__plus") != std::string::npos)
    {
        syst = 1024;
    }
    if (label.find("__met__minus") != std::string::npos)
    {
        syst = 2048;
    }

    if (channel == "emu")
    {
        inputVars.at("chan") = 2.;
    }
    if (channel == "ee")
    {
        inputVars.at("chan") = 1.;
    }
    if (channel == "mumu")
    {
        inputVars.at("chan") = 0.;
    }

    inputVars.at("eventNumber") = tree->eventNum;

    const std::pair<TLorentzVector, TLorentzVector> zPairLeptons{
        sortOutLeptons(tree, channel)};
    const TLorentzVector zLep1{zPairLeptons.first};
    const TLorentzVector zLep2{zPairLeptons.second};

    TLorentzVector metVec;

    if (oldMetFlag)
    {
        metVec.SetPtEtaPhiE(
            tree->metPF2PATEt, 0, tree->metPF2PATPhi, tree->metPF2PATEt);
    }
    else
    {
        if (syst == 1024)
        {
            metVec.SetPtEtaPhiE(tree->metPF2PATUnclusteredEnUp,
                                0,
                                tree->metPF2PATPhi,
                                tree->metPF2PATUnclusteredEnUp);
        }
        else if (syst == 2048)
        {
            metVec.SetPtEtaPhiE(tree->metPF2PATUnclusteredEnDown,
                                0,
                                tree->metPF2PATPhi,
                                tree->metPF2PATUnclusteredEnDown);
        }
        else
        {
            metVec.SetPtEtaPhiE(
                tree->metPF2PATEt, 0, tree->metPF2PATPhi, tree->metPF2PATEt);
        }
    }

    const std::pair<std::vector<int>, std::vector<TLorentzVector>> jetPair{
        getJets(tree, syst, metVec)};
    const std::vector<int> jets{jetPair.first};
    const std::vector<TLorentzVector> jetVecs{jetPair.second};

    const std::pair<std::vector<int>, std::vector<TLorentzVector>> bJetPair{
        getBjets(tree, syst, metVec, jets)};
    const std::vector<int> bJets{bJetPair.first};
    const std::vector<TLorentzVector> bJetVecs{bJetPair.second};

    const std::pair<TLorentzVector, TLorentzVector> wQuarkPair{
        sortOutHadronicW(tree)};
    const TLorentzVector wQuark1{wQuarkPair.first};
    const TLorentzVector wQuark2{wQuarkPair.second};

    // Do unclustered met stuff here now that we have all of the objects, all
    // corrected for their various SFs etc ...
    if (oldMetFlag && (syst == 1024 || syst == 2048))
    {
        metVec = doUncMet(metVec, zLep1, zLep2, jetVecs, syst);
    }

    // SFs for NPL lepton estimation normilisation
    // mz20 mw 20, ee = 1.06773225071; mumu = 1.02492608673;
    constexpr double SF_EE{1.068};
    constexpr double SF_MUMU{1.114};

    if (SameSignMC == true && channel == "ee")
    {
        inputVars.at("eventWeight") = tree->eventWeight * SF_EE;
    }
    else if (SameSignMC == true && channel == "mumu")
    {
        inputVars.at("eventWeight") = tree->eventWeight * SF_MUMU;
    }
    else
    {
        inputVars.at("eventWeight") = tree->eventWeight;
    }

    inputVars.at("leadJetPt") = jetVecs[0].Pt();
    inputVars.at("leadJetEta") = jetVecs[0].Eta();
    inputVars.at("leadJetPhi") = jetVecs[0].Phi();

    float totPx{0.0};
    float totPy{0.0};

    totPx += zLep1.Px() + zLep2.Px();
    totPy += zLep1.Py() + zLep2.Py();
    inputVars.at("lep1Pt") = zLep1.Pt();
    inputVars.at("lep1Eta") = zLep1.Eta();
    inputVars.at("lep1Phi") = zLep1.Phi();
    inputVars.at("lep2Pt") = zLep2.Pt();
    inputVars.at("lep2Eta") = zLep2.Eta();
    inputVars.at("lep2Phi") = zLep2.Phi();

    if (channel == "ee")
    {
        inputVars.at("lep1RelIso") =
            tree->elePF2PATComRelIsoRho[tree->zLep1Index];
        inputVars.at("lep1D0") = tree->elePF2PATD0PV[tree->zLep1Index];
        inputVars.at("lep2RelIso") =
            tree->elePF2PATComRelIsoRho[tree->zLep2Index];
        inputVars.at("lep2D0") = tree->elePF2PATD0PV[tree->zLep2Index];
    }
    if (channel == "mumu")
    {
        inputVars.at("lep1RelIso") =
            tree->muonPF2PATComRelIsodBeta[tree->zLep1Index];
        inputVars.at("lep1D0") = tree->muonPF2PATDBPV[tree->zLep1Index];
        inputVars.at("lep2RelIso") =
            tree->muonPF2PATComRelIsodBeta[tree->zLep2Index];
        inputVars.at("lep2D0") = tree->muonPF2PATDBPV[tree->zLep2Index];
    }
    if (channel == "emu")
    {
        inputVars.at("lep1RelIso") =
            tree->elePF2PATComRelIsoRho[tree->zLep1Index];
        inputVars.at("lep1D0") = tree->elePF2PATD0PV[tree->zLep1Index];
        inputVars.at("lep2RelIso") =
            tree->muonPF2PATComRelIsodBeta[tree->zLep2Index];
        inputVars.at("lep2D0") = tree->muonPF2PATDBPV[tree->zLep2Index];
    }

    inputVars.at("lepMass") = (zLep1 + zLep2).M();
    inputVars.at("lepPt") = std::sqrt(totPx * totPx + totPy * totPy);
    inputVars.at("lepEta") = (zLep1 + zLep2).Eta();
    inputVars.at("lepPhi") = (zLep1 + zLep2).Phi();
    inputVars.at("wQuark1Pt") = wQuark1.Pt();
    inputVars.at("wQuark1Eta") = wQuark1.Eta();
    inputVars.at("wQuark1Phi") = wQuark1.Phi();
    inputVars.at("wQuark2Pt") = wQuark2.Pt();
    inputVars.at("wQuark2Eta") = wQuark2.Eta();
    inputVars.at("wQuark2Phi") = wQuark2.Phi();

    const double wPairMass{(wQuark1 + wQuark2).M()};
    inputVars.at("wPairMass") = wPairMass;
    inputVars.at("wPairPt") = (wQuark1 + wQuark2).Pt();
    inputVars.at("wPairEta") = (wQuark1 + wQuark2).Eta();
    inputVars.at("wPairPhi") = (wQuark1 + wQuark2).Phi();
    totPx += jetVecs[0].Px();
    totPy += jetVecs[0].Py();

    if (jetVecs.size() > 1)
    {
        totPx += jetVecs[1].Px();
        totPy += jetVecs[1].Py();
    }
    inputVars.at("totPt2Jet") = std::sqrt(totPx * totPx + totPy * totPy);

    for (unsigned i{2}; i != jetVecs.size(); i++)
    {
        totPx += jetVecs[i].Px();
        totPy += jetVecs[i].Py();
    }

    inputVars.at("totPt") = std::sqrt(totPx * totPx + totPy * totPy);
    TLorentzVector totVec{zLep1 + zLep2};

    for (const auto& jetVec : jetVecs)
    {
        totVec += jetVec;
    }

    inputVars.at("totEta") = totVec.Eta();
    inputVars.at("totEta") = totVec.Phi();
    inputVars.at("totPtVec") = totVec.Pt();
    inputVars.at("totVecM") = totVec.M();
    inputVars.at("mTW") =
        std::sqrt(2 * tree->jetPF2PATPt[tree->wQuark1Index]
                  * tree->jetPF2PATPt[tree->wQuark2Index]
                  * (1
                     - cos(tree->jetPF2PATPhi[tree->wQuark1Index]
                           - tree->jetPF2PATPhi[tree->wQuark2Index])));
    inputVars.at("nJets") = boost::numeric_cast<float>(jets.size());
    inputVars.at("nBjets") = boost::numeric_cast<float>(bJets.size());
    inputVars.at("met") = metVec.Pt();
    inputVars.at("bTagDisc") = tree->jetPF2PATBDiscriminator[jets[bJets[0]]];
    inputVars.at("leadJetbTag") = tree->jetPF2PATBDiscriminator[jets[0]];
    inputVars.at("secJetbTag") = NaN;
    inputVars.at("secJetPt") = NaN;
    inputVars.at("secJetEta") = NaN;
    inputVars.at("secJetPhi") = NaN;
    inputVars.at("thirdJetbTag") = NaN;
    inputVars.at("thirdJetPt") = NaN;
    inputVars.at("thirdJetEta") = NaN;
    inputVars.at("thirdJetPhi") = NaN;
    inputVars.at("fourthJetbTag") = NaN;
    inputVars.at("fourthJetPt") = NaN;
    inputVars.at("fourthJetEta") = NaN;
    inputVars.at("fourthJetPhi") = NaN;

    if (jetVecs.size() > 1)
    {
        inputVars.at("secJetPt") = jetVecs[1].Pt();
        inputVars.at("secJetEta") = jetVecs[1].Eta();
        inputVars.at("secJetPhi") = jetVecs[1].Phi();
        inputVars.at("secJetbTag") = tree->jetPF2PATBDiscriminator[jets[1]];
    }

    if (jetVecs.size() > 2)
    {
        inputVars.at("thirdJetPt") = jetVecs[2].Pt();
        inputVars.at("thirdJetEta") = jetVecs[2].Eta();
        inputVars.at("thirdJetPhi") = jetVecs[2].Phi();
        inputVars.at("thirdJetbTag") = tree->jetPF2PATBDiscriminator[jets[2]];
    }

    if (jetVecs.size() > 3)
    {
        inputVars.at("fourthJetPt") = jetVecs[3].Pt();
        inputVars.at("fourthJetEta") = jetVecs[3].Eta();
        inputVars.at("fourthJetPhi") = jetVecs[3].Phi();
        inputVars.at("fourthJetbTag") = tree->jetPF2PATBDiscriminator[jets[3]];
    }

    const double topMass{(bJetVecs[0] + wQuark1 + wQuark2).M()};
    inputVars.at("topMass") = topMass;
    inputVars.at("topPt") = (bJetVecs[0] + wQuark1 + wQuark2).Pt();
    inputVars.at("topEta") = (bJetVecs[0] + wQuark1 + wQuark2).Eta();
    inputVars.at("topPhi") = (bJetVecs[0] + wQuark1 + wQuark2).Phi();
    inputVars.at("wZdelR") = (zLep2 + zLep1).DeltaR(wQuark1 + wQuark2);
    inputVars.at("wZdelPhi") = (zLep2 + zLep1).DeltaPhi(wQuark1 + wQuark2);

    inputVars.at("zQuark1DelR") = (zLep2 + zLep1).DeltaR(wQuark1);
    inputVars.at("zQuark1DelPhi") = (zLep2 + zLep1).DeltaPhi(wQuark1);
    inputVars.at("zQuark2DelR") = (zLep2 + zLep1).DeltaR(wQuark2);
    inputVars.at("zQuark2DelPhi") = (zLep2 + zLep1).DeltaPhi(wQuark2);

    inputVars.at("zTopDelR") =
        (zLep2 + zLep1).DeltaR(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("zTopDelPhi") =
        (zLep2 + zLep1).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("zl1TopDelR") =
        (zLep1).DeltaR(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("zl1TopDelPhi") =
        (zLep1).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("zl2TopDelR") =
        (zLep2).DeltaR(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("zl2TopDelPhi") =
        (zLep2).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2);

    inputVars.at("wTopDelR") =
        (wQuark1 + wQuark2).DeltaR(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("wTopDelPhi") =
        (wQuark1 + wQuark2).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("w1TopDelR") =
        (wQuark1).DeltaR(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("w1TopDelR") =
        (wQuark1).DeltaR(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("w1TopDelPhi") =
        (wQuark1).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("w2TopDelR") =
        (wQuark2).DeltaR(bJetVecs[0] + wQuark1 + wQuark2);
    inputVars.at("w2TopDelPhi") =
        (wQuark2).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2);

    inputVars.at("j1j2delR") = NaN;
    inputVars.at("j1j2delPhi") = NaN;

    if (jetVecs.size() > 1)
    {
        inputVars.at("j1j2delR") = jetVecs[0].DeltaR(jetVecs[1]);
        inputVars.at("j1j2delPhi") = jetVecs[0].DeltaPhi(jetVecs[1]);
    }

    inputVars.at("w1w2delR") = (wQuark1).DeltaR(wQuark2);
    inputVars.at("w1w2delPhi") = (wQuark1).DeltaPhi(wQuark2);
    inputVars.at("zLepdelR") = (zLep1).DeltaR(zLep2);
    inputVars.at("zLepdelPhi") = (zLep1).DeltaPhi(zLep2);
    inputVars.at("zl1Quark1DelR") = (zLep1).DeltaR(wQuark1);
    inputVars.at("zl1Quark1DelPhi") = (zLep1).DeltaPhi(wQuark1);
    inputVars.at("zl1Quark2DelR") = (zLep1).DeltaR(wQuark2);
    inputVars.at("zl1Quark2DelPhi") = (zLep1).DeltaPhi(wQuark2);
    inputVars.at("zl2Quark1DelR") = (zLep2).DeltaR(wQuark1);
    inputVars.at("zl2Quark1DelPhi") = (zLep2).DeltaPhi(wQuark1);
    inputVars.at("zl2Quark2DelR") = (zLep2).DeltaR(wQuark2);
    inputVars.at("zl2Quark2DelPhi") = (zLep2).DeltaPhi(wQuark2);

    float jetHt{0.};
    TLorentzVector jetVector;
    inputVars.at("minZJetR") = std::numeric_limits<double>::infinity();
    inputVars.at("minZJetPhi") = std::numeric_limits<double>::infinity();

    for (const auto& jetVec : jetVecs)
    {
        jetHt += jetVec.Pt();
        jetVector += jetVec;
        if (jetVec.DeltaR(zLep2 + zLep1) < inputVars.at("minZJetR"))
        {
            inputVars.at("minZJetR") = jetVec.DeltaR(zLep2 + zLep1);
        }
        if (jetVec.DeltaPhi(zLep2 + zLep1) < inputVars.at("minZJetPhi"))
        {
            inputVars.at("minZJetPhi") = jetVec.DeltaPhi(zLep2 + zLep1);
        }
    }

    inputVars.at("zlb1DelR") = zLep1.DeltaR(bJetVecs[0]);
    inputVars.at("zlb1DelPhi") = zLep1.DeltaPhi(bJetVecs[0]);
    inputVars.at("zlb2DelR") = zLep2.DeltaR(bJetVecs[0]);
    inputVars.at("zlb2DelPhi") = zLep2.DeltaPhi(bJetVecs[0]);

    const double lepHt{zLep1.Pt() + zLep2.Pt()};

    inputVars.at("lepHt") = lepHt;
    inputVars.at("jetHt") = jetHt;
    inputVars.at("jetMass") = jetVector.M();
    inputVars.at("jetPt") = jetVector.Pt();
    inputVars.at("jetEta") = jetVector.Eta();
    inputVars.at("jetPhi") = jetVector.Phi();

    if (channel != "emu")
    {
        inputVars.at("jetMass3") = (jetVecs[0] + jetVecs[1] + jetVecs[2]).M();
    }
    else
    {
        inputVars.at("jetMass3") = (jetVecs[0] + jetVecs[1]).M();
    }

    inputVars.at("wQuarkHt") = wQuark1.Pt() + wQuark2.Pt();

    inputVars.at("totHt") = lepHt + jetHt;
    inputVars.at("totHtOverPt") =
        inputVars.at("totHt") / std::sqrt(totPx * totPx + totPy * totPy);
    inputVars.at("zMass") = (zLep1 + zLep2).M();
    inputVars.at("zPt") = (zLep2 + zLep1).Pt();
    inputVars.at("zEta") = (zLep2 + zLep1).Eta();
    inputVars.at("zPhi") = (zLep2 + zLep1).Phi();

    constexpr double W_MASS{80.385};
    constexpr double TOP_MASS{173.1};

    // from scripts/plotMassPeaks.py
    constexpr double W_SIGMA{8};
    constexpr double TOP_SIGMA{30};

    const double wChi2Term{(wPairMass - W_MASS) / W_SIGMA};
    const double topChi2Term{(topMass - TOP_MASS) / TOP_SIGMA};
    inputVars.at("chi2") = std::pow(wChi2Term, 2) + std::pow(topChi2Term, 2);

    constexpr double MIN_SIDEBAND_CHI2{40};
    constexpr double MAX_SIDEBAND_CHI2{150};

    if (useSidebandRegion)
    {
        if (inputVars.at("chi2") >= MIN_SIDEBAND_CHI2
            and inputVars.at("chi2") < MAX_SIDEBAND_CHI2)
        {
            outTreeSdBnd->Fill();
        }
        if (inputVars.at("chi2") < MIN_SIDEBAND_CHI2)
        {
            outTreeSig->Fill();
        }
    }
    else
    {
        outTreeSig->Fill();
    }
}
