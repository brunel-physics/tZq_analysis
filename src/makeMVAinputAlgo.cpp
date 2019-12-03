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
                    {"WW_lnu2q", "WW1l1nu2q"},
                    {"WW_2l2nu", "WW2l2nu"},
                    {"ZZ_4l", "ZZ4l"},
                    {"ZZ_2l2nu", "ZZ2l2nu"},
                    {"ZZ_2l2q", "ZZ2l2q"},
                    {"WZ_3lnu", "WZ3l3nu"},
                    {"WZ_2l2q", "WZ2l2q"},
                    {"WZ_lnu2q", "WZ1l1nu2q"},
                    {"WG_lnug", "WG11nu1g"},
                    {"ZG_lnug", "ZG2l1g"},
                    {"t_s_channel", "TsChan"},
                    {"t_t_channel", "TtChan"},
                    {"tbar_t_channel", "TbartChan"},
                    {"tW", "TW"},
                    {"tbarW", "TbarW"},
                    {"tZq", "TZQ"},
                    {"tHq", "THQ"},
                    {"ttWTolnu", "TTWlnu"},
                    {"ttW2q", "TTW2q"},
                    {"ttZToll", "TTZ2l2nu"},
                    {"ttZ2q", "TTZ2q"},
                    {"ttgamma", "TTG"},
                    {"ttbar_2l2v", "TT2l2v"},
                    {"ttbar_hadronic", "TTjets"},
                    {"ttbar_semileptonic", "TT1l1v2q"},
                    {"tWZ", "TWZ"},
                    {"wPlusJets", "Wjets"},
                    {"DYJetsToLL_M-10to50", "DYJetsToLLM10to50"},
                    {"DYJetsToLL_M-50", "DYJetsToLLM50"}};
        }
    }};

    const auto getListOfSysts{[this]() -> std::map<std::string, std::string> {
        if (is2016)
        { // 2016
            return {{"tChannel_scaleUp", "TtChan__scaleUp"},
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
                    {"tWInclusive_scaleUp", "TW__scaleUp"},
                    {"tWInclusive_scaleDown", "TW__scaleDown"},
                    {"tbarWInclusive_scaleUp", "TbarW__scaleUp"},
                    {"tbarWInclusive_scaleDown", "TbarW__scaleDown"},
                    {"tZq_scaleUp", "TZQ__scaleUp"},
                    {"tZq_scaleDown", "TZQ__scaleDown"}};
        }
        else
        { // 2017
            return {{"ttbar_2l2v_hdampUp", "TT2l2v__hdampUp"},
                    {"ttbar_2l2v_hdampDown", "TT2l2v__hdampDown"},
                    {"ttbar_semileptonic_hdampUp", "TT1l1v2q__hdampUp"},
                    {"ttbar_semileptonic_hdampDown", "TT1l1v2q__hdampDown"},
                    {"ttbar_semileptonic_hdampUp", "TTjets__hdampUp"},
                    {"ttbar_semileptonic_hdampDown", "TTjets__hdampDown"}};
        }
    }};

    static const auto listOfMCs = getListOfMCs();
    static const auto listOfSysts = getListOfSysts();
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

    std::vector<std::string> systs = {"",
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
    if (!is2016)
    {
        systs.emplace_back("__isr__plus");
        systs.emplace_back("__isr__minus");
        systs.emplace_back("__fsr__plus");
        systs.emplace_back("__fsr__minus");
    }

    if (doMC)
    {
        standardAnalysis(listOfMCs, systs, channels, useSidebandRegion);
    }
    if (doSysts)
    {
        standardAnalysis(listOfSysts, {""}, channels, useSidebandRegion);
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
                auto event{new MvaEvent{true, tree, is2016}};

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
                                     / nominalEvents[channel])
                                    * 100)
                          << std::endl;

                inFile->Close();
            } // end channel loop
            outFile->cd();
            outTreeSig->SetDirectory(outFile);
            outTreeSig->FlushBaskets();
            if (useSidebandRegion)
            {
                outTreeSdBnd->SetDirectory(outFile);
                outTreeSdBnd->FlushBaskets();
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

        MvaEvent event{false, &dataChain, is2016};
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
        outTreeSig->SetDirectory(&outFile);
        outTreeSig->FlushBaskets();
        if (useSidebandRegion)
        {
            outTreeSdBnd->SetDirectory(&outFile);
            outTreeSdBnd->FlushBaskets();
        }
        outFile.Write();
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
        MvaEvent event{false, &dataChain, is2016};

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
        outTreeSig->SetDirectory(outFile);
        outTreeSig->FlushBaskets();
        if (useSidebandRegion)
        {
            outTreeSdBnd->SetDirectory(outFile);
            outTreeSdBnd->FlushBaskets();
        }
        outFile->Write();
        outFile->Close();
    }
}

std::pair<TLorentzVector, TLorentzVector>
    MakeMvaInputs::sortOutLeptons(const MvaEvent* tree,
                                  const std::string& channel) const
{
    TLorentzVector zLep1;
    TLorentzVector zLep2;

    const int zl1Index{tree->zLep1Index};
    const int zl2Index{tree->zLep2Index};

    if (channel == "ee")
    {
        zLep1.SetPxPyPzE(tree->elePF2PATPX[zl1Index],
                         tree->elePF2PATPY[zl1Index],
                         tree->elePF2PATPZ[zl1Index],
                         tree->elePF2PATE[zl1Index]);
        zLep2.SetPxPyPzE(tree->elePF2PATPX[zl2Index],
                         tree->elePF2PATPY[zl2Index],
                         tree->elePF2PATPZ[zl2Index],
                         tree->elePF2PATE[zl2Index]);
    }
    else if (channel == "mumu")
    {
        zLep1.SetPxPyPzE(tree->muonPF2PATPX[zl1Index],
                         tree->muonPF2PATPY[zl1Index],
                         tree->muonPF2PATPZ[zl1Index],
                         tree->muonPF2PATE[zl1Index]);
        zLep2.SetPxPyPzE(tree->muonPF2PATPX[zl2Index],
                         tree->muonPF2PATPY[zl2Index],
                         tree->muonPF2PATPZ[zl2Index],
                         tree->muonPF2PATE[zl2Index]);

        zLep1 *= tree->muonMomentumSF[0];
        zLep2 *= tree->muonMomentumSF[1];
    }
    else if (channel == "emu")
    {
        if (tree->muonLeads)
        {
            zLep1.SetPxPyPzE(tree->muonPF2PATPX[zl1Index],
                    tree->muonPF2PATPY[zl1Index],
                    tree->muonPF2PATPZ[zl1Index],
                    tree->muonPF2PATE[zl1Index]);
            zLep2.SetPxPyPzE(tree->elePF2PATPX[zl2Index],
                    tree->elePF2PATPY[zl2Index],
                    tree->elePF2PATPZ[zl2Index],
                    tree->elePF2PATE[zl2Index]);

            zLep1 *= tree->muonMomentumSF[0];
        }
        else
        {
            zLep1.SetPxPyPzE(tree->elePF2PATPX[zl1Index],
                    tree->elePF2PATPY[zl1Index],
                    tree->elePF2PATPZ[zl1Index],
                    tree->elePF2PATE[zl1Index]);
            zLep2.SetPxPyPzE(tree->muonPF2PATPX[zl2Index],
                    tree->muonPF2PATPY[zl2Index],
                    tree->muonPF2PATPZ[zl2Index],
                    tree->muonPF2PATE[zl2Index]);

            zLep2 *= tree->muonMomentumSF[0];
        }
    }

    return {zLep1, zLep2};
}

std::pair<TLorentzVector, TLorentzVector>
    MakeMvaInputs::sortOutHadronicW(const MvaEvent* tree,
                                    const int syst,
                                    TLorentzVector met,
                                    const std::vector<int>& jets) const
{
    const auto wQuark1{getJetVec(tree,
                                 jets.at(tree->wQuark1Index),
                                 tree->jetSmearValue[tree->wQuark1Index],
                                 met,
                                 syst,
                                 false)};
    const auto wQuark2{getJetVec(tree,
                                 jets.at(tree->wQuark2Index),
                                 tree->jetSmearValue[tree->wQuark2Index],
                                 met,
                                 syst,
                                 false)};

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
    tree->Branch("Channel", &inputVars["chan"], "Channel/F");
    tree->Branch("EvtNumber", &inputVars["eventNumber"], "EvtNumber/F");
    tree->Branch("EvtWeight", &inputVars["eventWeight"], "EvtWeight/F");
    tree->Branch("bEta", &inputVars["bEta"], "bEta/F");
    tree->Branch("bPhi", &inputVars["bPhi"], "bPhi/F");
    tree->Branch("bPt", &inputVars["bPt"], "bPt/F");
    tree->Branch("bbTag", &inputVars["bbTag"], "bbTag/F");
    tree->Branch("chi2", &inputVars["chi2"], "chi2/F");
    tree->Branch("j1Eta", &inputVars["j1Eta"], "j1Eta/F");
    tree->Branch("j1Phi", &inputVars["j1Phi"], "j1Phi/F");
    tree->Branch("j1Pt", &inputVars["j1Pt"], "j1Pt/F");
    tree->Branch("j1bDelR", &inputVars["j1bDelR"], "j1bDelR/F");
    tree->Branch("j1bTag", &inputVars["j1bTag"], "j1bTag/F");
    tree->Branch("j1j2DelR", &inputVars["j1j2DelR"], "j1j2DelR/F");
    tree->Branch("j1j3DelR", &inputVars["j1j3DelR"], "j1j3DelR/F");
    tree->Branch("j1j4DelR", &inputVars["j1j4DelR"], "j1j4DelR/F");
    tree->Branch("j1l1DelR", &inputVars["j1l1DelR"], "j1l1DelR/F");
    tree->Branch("j1l2DelR", &inputVars["j1l2DelR"], "j1l2DelR/F");
    tree->Branch("j1tDelR", &inputVars["j1tDelR"], "j1tDelR/F");
    tree->Branch("j1wDelR", &inputVars["j1wDelR"], "j1wDelR/F");
    tree->Branch("j1wj1DelR", &inputVars["j1wj1DelR"], "j1wj1DelR/F");
    tree->Branch("j1wj2DelR", &inputVars["j1wj2DelR"], "j1wj2DelR/F");
    tree->Branch("j1zDelR", &inputVars["j1zDelR"], "j1zDelR/F");
    tree->Branch("j2Eta", &inputVars["j2Eta"], "j2Eta/F");
    tree->Branch("j2Phi", &inputVars["j2Phi"], "j2Phi/F");
    tree->Branch("j2Pt", &inputVars["j2Pt"], "j2Pt/F");
    tree->Branch("j2bDelR", &inputVars["j2bDelR"], "j2bDelR/F");
    tree->Branch("j2bTag", &inputVars["j2bTag"], "j2bTag/F");
    tree->Branch("j2j3DelR", &inputVars["j2j3DelR"], "j2j3DelR/F");
    tree->Branch("j2j4DelR", &inputVars["j2j4DelR"], "j2j4DelR/F");
    tree->Branch("j2l1DelR", &inputVars["j2l1DelR"], "j2l1DelR/F");
    tree->Branch("j2l2DelR", &inputVars["j2l2DelR"], "j2l2DelR/F");
    tree->Branch("j2tDelR", &inputVars["j2tDelR"], "j2tDelR/F");
    tree->Branch("j2wDelR", &inputVars["j2wDelR"], "j2wDelR/F");
    tree->Branch("j2wj1DelR", &inputVars["j2wj1DelR"], "j2wj1DelR/F");
    tree->Branch("j2wj2DelR", &inputVars["j2wj2DelR"], "j2wj2DelR/F");
    tree->Branch("j2zDelR", &inputVars["j2zDelR"], "j2zDelR/F");
    tree->Branch("j3Eta", &inputVars["j3Eta"], "j3Eta/F");
    tree->Branch("j3Phi", &inputVars["j3Phi"], "j3Phi/F");
    tree->Branch("j3Pt", &inputVars["j3Pt"], "j3Pt/F");
    tree->Branch("j3bDelR", &inputVars["j3bDelR"], "j3bDelR/F");
    tree->Branch("j3bTag", &inputVars["j3bTag"], "j3bTag/F");
    tree->Branch("j3j4DelR", &inputVars["j3j4DelR"], "j3j4DelR/F");
    tree->Branch("j3l1DelR", &inputVars["j3l1DelR"], "j3l1DelR/F");
    tree->Branch("j3l2DelR", &inputVars["j3l2DelR"], "j3l2DelR/F");
    tree->Branch("j3tDelR", &inputVars["j3tDelR"], "j3tDelR/F");
    tree->Branch("j3wDelR", &inputVars["j3wDelR"], "j3wDelR/F");
    tree->Branch("j3wj1DelR", &inputVars["j3wj1DelR"], "j3wj1DelR/F");
    tree->Branch("j3wj2DelR", &inputVars["j3wj2DelR"], "j3wj2DelR/F");
    tree->Branch("j3zDelR", &inputVars["j3zDelR"], "j3zDelR/F");
    tree->Branch("j4Eta", &inputVars["j4Eta"], "j4Eta/F");
    tree->Branch("j4Phi", &inputVars["j4Phi"], "j4Phi/F");
    tree->Branch("j4Pt", &inputVars["j4Pt"], "j4Pt/F");
    tree->Branch("j4bDelR", &inputVars["j4bDelR"], "j4bDelR/F");
    tree->Branch("j4bTag", &inputVars["j4bTag"], "j4bTag/F");
    tree->Branch("j4l1DelR", &inputVars["j4l1DelR"], "j4l1DelR/F");
    tree->Branch("j4l2DelR", &inputVars["j4l2DelR"], "j4l2DelR/F");
    tree->Branch("j4tDelR", &inputVars["j4tDelR"], "j4tDelR/F");
    tree->Branch("j4wDelR", &inputVars["j4wDelR"], "j4wDelR/F");
    tree->Branch("j4wj1DelR", &inputVars["j4wj1DelR"], "j4wj1DelR/F");
    tree->Branch("j4wj2DelR", &inputVars["j4wj2DelR"], "j4wj2DelR/F");
    tree->Branch("j4zDelR", &inputVars["j4zDelR"], "j4zDelR/F");
    tree->Branch("j5Eta", &inputVars["j5Eta"], "j5Eta/F");
    tree->Branch("j5Phi", &inputVars["j5Phi"], "j5Phi/F");
    tree->Branch("j5Pt", &inputVars["j5Pt"], "j5Pt/F");
    tree->Branch("j5bTag", &inputVars["j5bTag"], "j5bTag/F");
    tree->Branch("j6Eta", &inputVars["j6Eta"], "j6Eta/F");
    tree->Branch("j6Phi", &inputVars["j6Phi"], "j6Phi/F");
    tree->Branch("j6Pt", &inputVars["j6Pt"], "j6Pt/F");
    tree->Branch("j6bTag", &inputVars["j6bTag"], "j6bTag/F");
    tree->Branch("jetMass", &inputVars["jetMass"], "jetMass/F");
    tree->Branch("jetMass3", &inputVars["jetMass3"], "jetMass3/F");
    tree->Branch("jetMt", &inputVars["jetMt"], "jetMt/F");
    tree->Branch("jetPt", &inputVars["jetPt"], "jetPt/F");
    tree->Branch("l1D0", &inputVars["l1D0"], "l1D0/F");
    tree->Branch("l1Eta", &inputVars["l1Eta"], "l1Eta/F");
    tree->Branch("l1Phi", &inputVars["l1Phi"], "l1Phi/F");
    tree->Branch("l1Pt", &inputVars["l1Pt"], "l1Pt/F");
    tree->Branch("l1RelIso", &inputVars["l1RelIso"], "l1RelIso/F");
    tree->Branch("l1bDelR", &inputVars["l1bDelR"], "l1bDelR/F");
    tree->Branch("l1tDelR", &inputVars["l1tDelR"], "l1tDelR/F");
    tree->Branch("l1wj1DelR", &inputVars["l1wj1DelR"], "l1wj1DelR/F");
    tree->Branch("l1wj2DelR", &inputVars["l1wj2DelR"], "l1wj2DelR/F");
    tree->Branch("l2D0", &inputVars["l2D0"], "l2D0/F");
    tree->Branch("l2DelR", &inputVars["l2DelR"], "l2DelR/F");
    tree->Branch("l2Eta", &inputVars["l2Eta"], "l2Eta/F");
    tree->Branch("l2Phi", &inputVars["l2Phi"], "l2Phi/F");
    tree->Branch("l2Pt", &inputVars["l2Pt"], "l2Pt/F");
    tree->Branch("l2RelIso", &inputVars["l2RelIso"], "l2RelIso/F");
    tree->Branch("l2bDelR", &inputVars["l2bDelR"], "l2bDelR/F");
    tree->Branch("l2tDelR", &inputVars["l2tDelR"], "l2tDelR/F");
    tree->Branch("l2wj1DelR", &inputVars["l2wj1DelR"], "l2wj1DelR/F");
    tree->Branch("l2wj2DelR", &inputVars["l2wj2DelR"], "l2wj2DelR/F");
    tree->Branch("met", &inputVars["met"], "met/F");
    tree->Branch("nBjets", &inputVars["nBjets"], "nBjets/F");
    tree->Branch("nJets", &inputVars["nJets"], "nJets/F");
    tree->Branch("tEta", &inputVars["tEta"], "tEta/F");
    tree->Branch("tMass", &inputVars["tMass"], "tMass/F");
    tree->Branch("tMt", &inputVars["tMt"], "tMt/F");
    tree->Branch("tPhi", &inputVars["tPhi"], "tPhi/F");
    tree->Branch("tPt", &inputVars["tPt"], "tPt/F");
    tree->Branch("tbDelR", &inputVars["tbDelR"], "tbDelR/F");
    tree->Branch("totMass", &inputVars["totMass"], "toMass/F");
    tree->Branch("totMt", &inputVars["totMt"], "totMt/F");
    tree->Branch("totPt", &inputVars["totPt"], "totPt/F");
    tree->Branch("wEta", &inputVars["wEta"], "wEta/F");
    tree->Branch("wMass", &inputVars["wMass"], "wMass/F");
    tree->Branch("wMt", &inputVars["wMt"], "wMt/F");
    tree->Branch("wPhi", &inputVars["wPhi"], "wPhi/F");
    tree->Branch("wPt", &inputVars["wPt"], "wPt/F");
    tree->Branch("wbDelR", &inputVars["wbDelR"], "wbDelR/F");
    tree->Branch("wj1DelR", &inputVars["wj1DelR"], "wj1DelR/F");
    tree->Branch("wj1Eta", &inputVars["wj1Eta"], "wj1Eta/F");
    tree->Branch("wj1Phi", &inputVars["wj1Phi"], "wj1Phi/F");
    tree->Branch("wj1Pt", &inputVars["wj1Pt"], "wj1Pt/F");
    tree->Branch("wj1bDelR", &inputVars["wj1bDelR"], "wj1bDelR/F");
    tree->Branch("wj1tDelR", &inputVars["wj1tDelR"], "wj1tDelR/F");
    tree->Branch("wj2DelR", &inputVars["wj2DelR"], "wj2DelR/F");
    tree->Branch("wj2Eta", &inputVars["wj2Eta"], "wj2Eta/F");
    tree->Branch("wj2Phi", &inputVars["wj2Phi"], "wj2Phi/F");
    tree->Branch("wj2Pt", &inputVars["wj2Pt"], "wj2Pt/F");
    tree->Branch("wj2bDelR", &inputVars["wj2bDelR"], "wj2bDelR/F");
    tree->Branch("wj2tDelR", &inputVars["wj2tDelR"], "wj2tDelR/F");
    tree->Branch("wtDelR", &inputVars["wtDelR"], "wtDelR/F");
    tree->Branch("wwDelR", &inputVars["w1w2DelR"], "wwDelR/F");
    tree->Branch("wzDelR", &inputVars["wZDelR"], "wzDelR/F");
    tree->Branch("zEta", &inputVars["zEta"], "zEta/F");
    tree->Branch("zMass", &inputVars["zMass"], "zMass/F");
    tree->Branch("zMt", &inputVars["zMt"], "zMt/F");
    tree->Branch("zPhi", &inputVars["zPhi"], "zPhi/F");
    tree->Branch("zPt", &inputVars["zPt"], "zPt/F");
    tree->Branch("zbDelR", &inputVars["zbDelR"], "zbDelR/F");
    tree->Branch("zjMaxR", &inputVars["zjMaxR"], "zjMaxR/F");
    tree->Branch("zjMinR", &inputVars["zjMinR"], "zjMinR/F");
    tree->Branch("ztDelR", &inputVars["ztDelR"], "ztDelR/F");
    tree->Branch("zwj1DelR", &inputVars["zwj1DelR"], "zwj1DelR/F");
    tree->Branch("zwj2DelR", &inputVars["zwj2DelR"], "zwj2DelR/F");
    tree->Branch("zzDelR", &inputVars["zzDelR"], "zzDelR/F");
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
        sortOutHadronicW(tree, syst, metVec, jets)};
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
    const double SF_EE{is2016 ? 0.956296241886 : 1.57352187847};
    const double SF_MUMU{is2016 ? 1.01737450337 : 1.15322585438};

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

    inputVars.at("j1Pt") = jetVecs[0].Pt();
    inputVars.at("j1Eta") = jetVecs[0].Eta();
    inputVars.at("j1Phi") = jetVecs[0].Phi();

    inputVars.at("l1Pt") = zLep1.Pt();
    inputVars.at("l1Eta") = zLep1.Eta();
    inputVars.at("l1Phi") = zLep1.Phi();
    inputVars.at("l2Pt") = zLep2.Pt();
    inputVars.at("l2Eta") = zLep2.Eta();
    inputVars.at("l2Phi") = zLep2.Phi();

    const TLorentzVector zVec{zLep1 + zLep2};
    inputVars.at("zMass") = zVec.M();
    // if (abs(zVec.M() - 91.1876) > 100)
    // {
    //     std::cout << tree->muonLeads << '\t' << zVec.M() << std::endl;;
    // }
    inputVars.at("zPt") = zVec.Pt();
    inputVars.at("zEta") = zVec.Eta();
    inputVars.at("zPhi") = zVec.Phi();
    inputVars.at("zMt") = zVec.Mt();

    const TLorentzVector wVec{wQuark1 + wQuark2};
    const double wMass{(wQuark1 + wQuark2).M()};
    inputVars.at("wMass") = wMass;
    inputVars.at("wPt") = wVec.Pt();
    inputVars.at("wEta") = wVec.Eta();
    inputVars.at("wPhi") = wVec.Phi();

    const TLorentzVector tVec{bJetVecs[0] + wVec};
    const double topMass{tVec.M()};
    inputVars.at("tMass") = topMass;
    inputVars.at("tMt") = tVec.Mt();
    inputVars.at("tPt") = tVec.Pt();
    inputVars.at("tEta") = tVec.Eta();
    inputVars.at("tPhi") = tVec.Phi();

    if (channel == "ee")
    {
        inputVars.at("l1RelIso") =
            tree->elePF2PATComRelIsoRho[tree->zLep1Index];
        inputVars.at("l1D0") = tree->elePF2PATD0PV[tree->zLep1Index];
        inputVars.at("l2RelIso") =
            tree->elePF2PATComRelIsoRho[tree->zLep2Index];
        inputVars.at("l2D0") = tree->elePF2PATD0PV[tree->zLep2Index];
    }
    if (channel == "mumu")
    {
        inputVars.at("l1RelIso") =
            tree->muonPF2PATComRelIsodBeta[tree->zLep1Index];
        inputVars.at("l1D0") = tree->muonPF2PATDBPV[tree->zLep1Index];
        inputVars.at("l2RelIso") =
            tree->muonPF2PATComRelIsodBeta[tree->zLep2Index];
        inputVars.at("l2D0") = tree->muonPF2PATDBPV[tree->zLep2Index];
    }
    if (channel == "emu")
    {
        inputVars.at("l1RelIso") =
            tree->elePF2PATComRelIsoRho[tree->zLep1Index];
        inputVars.at("l1D0") = tree->elePF2PATD0PV[tree->zLep1Index];
        inputVars.at("l2RelIso") =
            tree->muonPF2PATComRelIsodBeta[tree->zLep2Index];
        inputVars.at("l2D0") = tree->muonPF2PATDBPV[tree->zLep2Index];
    }

    inputVars.at("wj1Pt") = wQuark1.Pt();
    inputVars.at("wj1Eta") = wQuark1.Eta();
    inputVars.at("wj1Phi") = wQuark1.Phi();
    inputVars.at("wj2Pt") = wQuark2.Pt();
    inputVars.at("wj2Eta") = wQuark2.Eta();
    inputVars.at("wj2Phi") = wQuark2.Phi();

    TLorentzVector totVec{zVec};
    for (const auto& jetVec : jetVecs)
    {
        totVec += jetVec;
    }

    inputVars.at("totPt") = totVec.Pt();
    inputVars.at("totMass") = totVec.M();
    inputVars.at("totMt") = totVec.Mt();
    inputVars.at("wMt") = wVec.Mt();
    inputVars.at("nJets") = boost::numeric_cast<float>(jets.size());
    inputVars.at("nBjets") = boost::numeric_cast<float>(bJets.size());
    inputVars.at("met") = metVec.Et();

    inputVars.at("bbTag") =
        tree->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
            [jets[bJets[0]]];
    inputVars.at("bPt") = bJetVecs[0].Pt();
    inputVars.at("bEta") = bJetVecs[0].Eta();
    inputVars.at("bPhi") = bJetVecs[0].Phi();

    inputVars.at("j1bTag") =
        tree->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[jets[0]];
    inputVars.at("j2bTag") = 0;
    inputVars.at("j2Pt") = 0;
    inputVars.at("j2Eta") = 0;
    inputVars.at("j2Phi") = 0;
    inputVars.at("j3bTag") = 0;
    inputVars.at("j3Pt") = 0;
    inputVars.at("j3Eta") = 0;
    inputVars.at("j3Phi") = 0;
    inputVars.at("j4bTag") = 0;
    inputVars.at("j4Pt") = 0;
    inputVars.at("j4Eta") = 0;
    inputVars.at("j4Phi") = 0;

    if (jetVecs.size() > 1)
    {
        inputVars.at("j2Pt") = jetVecs[1].Pt();
        inputVars.at("j2Eta") = jetVecs[1].Eta();
        inputVars.at("j2Phi") = jetVecs[1].Phi();
        inputVars.at("j2bTag") =
            tree->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                [jets[1]];
    }
    if (jetVecs.size() > 2)
    {
        inputVars.at("j3Pt") = jetVecs[2].Pt();
        inputVars.at("j3Eta") = jetVecs[2].Eta();
        inputVars.at("j3Phi") = jetVecs[2].Phi();
        inputVars.at("j3bTag") =
            tree->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                [jets[2]];
    }
    if (jetVecs.size() > 3)
    {
        inputVars.at("j4Pt") = jetVecs[3].Pt();
        inputVars.at("j4Eta") = jetVecs[3].Eta();
        inputVars.at("j4Phi") = jetVecs[3].Phi();
        inputVars.at("j4bTag") =
            tree->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                [jets[3]];
    }
    if (jetVecs.size() > 4)
    {
        inputVars.at("j5Pt") = jetVecs[4].Pt();
        inputVars.at("j5Eta") = jetVecs[4].Eta();
        inputVars.at("j5Phi") = jetVecs[4].Phi();
        inputVars.at("j5bTag") =
            tree->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                [jets[4]];
    }
    if (jetVecs.size() > 5)
    {
        inputVars.at("j6Pt") = jetVecs[5].Pt();
        inputVars.at("j6Eta") = jetVecs[5].Eta();
        inputVars.at("j6Phi") = jetVecs[5].Phi();
        inputVars.at("j6bTag") =
            tree->jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags
                [jets[5]];
    }

    inputVars.at("wZDelR") = zVec.DeltaR(wVec);

    inputVars.at("zwj1DelR") = zVec.DeltaR(wQuark1);
    inputVars.at("zwj2DelR") = zVec.DeltaR(wQuark2);

    inputVars.at("ztDelR") = zVec.DeltaR(tVec);
    inputVars.at("l1tDelR") = zLep1.DeltaR(tVec);
    inputVars.at("l2tDelR") = zLep2.DeltaR(tVec);

    inputVars.at("wtDelR") = wVec.DeltaR(tVec);
    inputVars.at("wj1tDelR") = wQuark1.DeltaR(tVec);
    inputVars.at("wj1tDelR") = wQuark1.DeltaPhi(tVec);
    inputVars.at("wj2tDelR") = wQuark2.DeltaR(tVec);

    inputVars.at("wtDelR") = wVec.DeltaR(tVec);

    inputVars.at("zbDelR") = zVec.DeltaR(bJetVecs[0]);
    inputVars.at("tbDelR") = tVec.DeltaR(bJetVecs[0]);
    inputVars.at("wbDelR") = wVec.DeltaR(bJetVecs[0]);
    inputVars.at("wj1bDelR") = wQuark1.DeltaR(bJetVecs[0]);
    inputVars.at("wj2bDelR") = wQuark2.DeltaR(bJetVecs[0]);
    inputVars.at("l1bDelR") = zLep1.DeltaR(bJetVecs[0]);
    inputVars.at("l2bDelR") = zLep2.DeltaR(bJetVecs[0]);

    for (unsigned i{0}; i < 4; i++)
    {
        for (unsigned j{i + 1}; j < 4; j++)
        {
            inputVars.at("j" + std::to_string(i + 1) + "j"
                         + std::to_string(j + 1) + "DelR") =
                jetVecs.at(i).DeltaR(jetVecs.at(j));
        }
        inputVars.at("j" + std::to_string(i + 1) + "bDelR") =
            jetVecs.at(i).DeltaR(bJetVecs.at(0));
        inputVars.at("j" + std::to_string(i + 1) + "tDelR") =
            jetVecs.at(i).DeltaR(tVec);

        inputVars.at("j" + std::to_string(i + 1) + "l1DelR") =
            jetVecs.at(i).DeltaR(zLep1);
        inputVars.at("j" + std::to_string(i + 1) + "l2DelR") =
            jetVecs.at(i).DeltaR(zLep2);
        inputVars.at("j" + std::to_string(i + 1) + "zDelR") =
            jetVecs.at(i).DeltaR(zVec);

        inputVars.at("j" + std::to_string(i + 1) + "wj1DelR") =
            jetVecs.at(i).DeltaR(wQuark1);
        inputVars.at("j" + std::to_string(i + 1) + "wj2DelR") =
            jetVecs.at(i).DeltaR(wQuark2);
        inputVars.at("j" + std::to_string(i + 1) + "wDelR") =
            jetVecs.at(i).DeltaR(wVec);
    }

    inputVars.at("w1w2DelR") = wQuark1.DeltaR(wQuark2);
    inputVars.at("zzDelR") = zLep1.DeltaR(zLep2);
    inputVars.at("l1wj1DelR") = zLep1.DeltaR(wQuark1);
    inputVars.at("l1wj2DelR") = zLep1.DeltaR(wQuark2);
    inputVars.at("l2wj1DelR") = zLep2.DeltaR(wQuark1);
    inputVars.at("l2wj2DelR") = zLep2.DeltaR(wQuark2);

    TLorentzVector jetVector;
    inputVars.at("zjMinR") = std::numeric_limits<float>::infinity();
    inputVars.at("zjMaxR") = -std::numeric_limits<float>::infinity();

    for (const auto& jetVec : jetVecs)
    {
        jetVector += jetVec;
        if (jetVec.DeltaR(zVec) < inputVars.at("zjMinR"))
        {
            inputVars.at("zjMinR") = jetVec.DeltaR(zVec);
        }
        if (jetVec.DeltaR(zVec) > inputVars.at("zjMaxR"))
        {
            inputVars.at("zjMaxR") = jetVec.DeltaR(zVec);
        }
    }

    inputVars.at("jetMass") = jetVector.M();
    inputVars.at("jetMt") = jetVector.Mt();
    inputVars.at("jetPt") = jetVector.Pt();

    inputVars.at("jetMass3") = (jetVecs[0] + jetVecs[1] + jetVecs[2]).M();

    constexpr double W_MASS{80.385};
    constexpr double TOP_MASS{173.1};

    // from scripts/plotMassPeaks.py
    constexpr double W_SIGMA{8};
    constexpr double TOP_SIGMA{30};

    const double wChi2Term{(wMass - W_MASS) / W_SIGMA};
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
