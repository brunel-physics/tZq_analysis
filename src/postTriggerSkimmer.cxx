#include "AnalysisEvent.hpp"

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <array>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/range/iterator_range.hpp>
#include <iostream>
#include <regex>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std::string_literals;
namespace fs = boost::filesystem;

int main(int argc, char* argv[])
{
    std::vector<std::string> dileptonDirs;
    std::vector<std::string> singleLeptonDirs;
    std::string datasetName;
    std::string channel;
    bool is2016;

    int singleElectron{0};
    int dupElectron{0};
    int singleMuon{0};
    int dupMuon{0};

    const std::string postTriggerSkimDir{
        "/data0/data/TopPhysics/postTriggerSkims201"s + (is2016 ? "6/" : "7/")};

    // Define command-line flags
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()("help,h", "Print this message.")(
        "channel,c",
        po::value<std::string>(&channel)->required(),
        "Channel to operate over. Either ee, emu or mumu.")(
        "2016", po::bool_switch(&is2016), "Use 2016 conditions (SFs, et al.).")(
        "dileptonDirs,d",
        po::value<std::vector<std::string>>(&dileptonDirs)
            ->multitoken()
            ->required(),
        "Directories in which to look for double lepton datasets.")(
        "singleLeptonDirs,s",
        po::value<std::vector<std::string>>(&singleLeptonDirs)
            ->multitoken()
            ->required(),
        "Directories in which to look for single lepton datasets.")(
        "datasetName,o",
        po::value<std::string>(&datasetName)->required(),
        "Output dataset name.");
    po::variables_map vm;

    // Parse arguments
    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc;
            return 0;
        }

        po::notify(vm);
    }
    catch (const po::error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    const std::regex mask{".*\\.root"};
    int fileNum{0};

    std::unordered_set<std::pair<int, int>,
                       boost::hash<std::pair<Int_t, Int_t>>>
        triggerDoubleCountCheck;

    for (const auto& dileptonDir : dileptonDirs)
    { // for each dilepton input directory
        for (const auto& file : boost::make_iterator_range(
                 fs::directory_iterator{dileptonDir}, {}))
        { // for each file in directory
            const std::string path{file.path().string()};
            std::cout << "path: " << path << std::endl;

            if (!fs::is_regular_file(file.status())
                || !std::regex_match(path, mask))
            {
                continue; // skip if not a root file
            }

            const std::string numName{std::to_string(fileNum)};
            const std::string numNamePlus{std::to_string(fileNum + 2)};

            if (fs::is_regular_file(postTriggerSkimDir + datasetName
                                    + "/triggerSkim" + numNamePlus + ".root"))
            {
                // don't overwrite existing skim files, except for the last two
                fileNum++;
                continue;
            }

            TChain datasetChain{"tree"};
            datasetChain.Add(path.c_str());

            TTree* const outTree = datasetChain.CloneTree(0);

            std::string outFilePath{postTriggerSkimDir + datasetName
                                    + "/triggerSkim" + numName + ".root"};
            TFile outFile{outFilePath.c_str(), "RECREATE"};

            const long long int numberOfEvents{datasetChain.GetEntries()};

            boost::progress_display progress{
                boost::numeric_cast<unsigned long>(numberOfEvents),
                std::cout,
                outFilePath + "\n"};
            AnalysisEvent event{false, &datasetChain, is2016};

            for (long long int i{0}; i < numberOfEvents; i++)
            {
                ++progress; // update progress bar (++ must be prefix)
                event.GetEntry(i);

                if (channel == "ee")
                {
                    // clang-format off
                    const bool eeTrig{event.eeTrig()};
                    // clang-format on

                    if (eeTrig)
                    {
                        triggerDoubleCountCheck.emplace(event.eventRun,
                                                        event.eventNum);
                        outTree->Fill();
                    }
                }

                if (channel == "mumu")
                {
                    const bool mumuTrig{event.mumuTrig()};

                    if (mumuTrig)
                    {
                        triggerDoubleCountCheck.emplace(event.eventRun,
                                                        event.eventNum);
                        outTree->Fill();
                    }
                }

                if (channel == "emu")
                {
                    const bool muEGTrig{event.muEGTrig()};

                    if (muEGTrig)
                    {
                        triggerDoubleCountCheck.emplace(event.eventRun,
                                                        event.eventNum);
                        outTree->Fill();
                    }
                } // end emu
            }
            outFile.cd();
            outTree->Write();

            outFile.Write();
            outFile.Close();

            fileNum++;
            std::cout << std::endl;

            delete outTree;
        }
    }

    for (const auto& singleLeptonDir : singleLeptonDirs)
    { // for each single lepton input directory
        for (const auto& file : boost::make_iterator_range(
                 fs::directory_iterator{singleLeptonDir}, {}))
        { // for each file in directory
            const std::string path{file.path().string()};
            std::cout << "path: " << path << std::endl;

            if (!fs::is_regular_file(file.status())
                || !std::regex_match(path, mask))
            {
                continue; // skip if not a root file
            }

            const std::string numName{std::to_string(fileNum)};
            const std::string numNamePlus{std::to_string(fileNum + 2)};

            if (fs::is_regular_file(postTriggerSkimDir + datasetName
                                    + "/triggerSkim" + numNamePlus + ".root"))
            {
                // don't overwrite existing skim files, except for the last two
                fileNum++;
                continue;
            }

            TChain datasetChain{"tree"};
            datasetChain.Add(path.c_str());
            TTree* const outTree = datasetChain.CloneTree(0);

            std::string outFilePath{postTriggerSkimDir + datasetName
                                    + "/triggerSkim" + numName + ".root"};
            TFile outFile{outFilePath.c_str(), "RECREATE"};

            const long long int numberOfEvents{datasetChain.GetEntries()};

            boost::progress_display progress{
                boost::numeric_cast<unsigned long>(numberOfEvents),
                std::cout,
                outFilePath + "\n"};
            AnalysisEvent event{false, &datasetChain, is2016};

            for (long long int i{0}; i < numberOfEvents; i++)
            {
                ++progress; // update progress bar (++ must be prefix)
                event.GetEntry(i);
                if (channel == "ee")
                {
                    const bool eTrig{event.eTrig()};
                    if (eTrig)
                    {
                        auto it{triggerDoubleCountCheck.find(
                            {event.eventRun, event.eventNum})};
                        singleElectron++;
                        // If event has already been found ... skip event
                        if (it != triggerDoubleCountCheck.end())
                        {
                            dupElectron++;
                        }
                        // If event has not already been found, add to new skim
                        else
                        {
                            // triggerDoubleCountCheck.emplace(
                            // event.eventRun, event.eventNum);
                            outTree->Fill();
                        }
                    }
                } // end single electron check for ee

                if (channel == "mumu")
                {
                    const bool muTrig{event.muTrig()};

                    // If single Muon triggered fired, check to see if event
                    // also fired a DoubleMuon trigger
                    if (muTrig)
                    {
                        singleMuon++;
                        auto it{triggerDoubleCountCheck.find(
                            {event.eventRun, event.eventNum})};
                        // If event has already been found ... skip event
                        if (it != triggerDoubleCountCheck.end())
                        {
                            dupMuon++;
                        }
                        // If event has not already been found, add to new skim
                        else
                        {
                            // triggerDoubleCountCheck.emplace(
                            // event.eventRun, event.eventNum);
                            outTree->Fill();
                        }
                    }

                } // end single muon check for mumu

                if (channel == "emu")
                {
                    // check eTrigger for emu first
                    // clang-format off
                    const bool eTrig{event.eTrig()};
                    // then check muTrigger for emu
                    const bool muTrig{event.muTrig()};
                    // clang-format on

                    // If either single lepton  triggered fired, check to see if
                    // event also fired a MuonEG trigger
                    if (eTrig || muTrig)
                    {
                        if (eTrig)
                        {
                            singleElectron++;
                        }
                        if (muTrig)
                        {
                            singleMuon++;
                        }
                        auto it{triggerDoubleCountCheck.find(
                            {event.eventRun, event.eventNum})};
                        // If event has already been found ... skip event
                        if (it != triggerDoubleCountCheck.end())
                        {
                            if (eTrig)
                            {
                                dupElectron++;
                            }
                            if (muTrig)
                            {
                                dupMuon++;
                            }
                        }
                        // If event has not already been found, add to new skim
                        else
                        {
                            // triggerDoubleCountCheck.emplace(
                            // event.eventRun, event.eventNum);
                            outTree->Fill();
                        }
                    }
                } // end single lepton check for emu
            }
            outFile.cd();
            outTree->Write();

            outFile.Write();
            outFile.Close();

            fileNum++;
            std::cout << std::endl;

            delete outTree;
        }
    }

    if (channel == "ee" || channel == "emu")
    {
        std::cout << "Single electron trigger fired with double lepton "
                     "trigger/Total single electron triggers fired: "
                  << dupElectron << " / " << singleElectron << std::endl;
    }
    if (channel == "mumu" || channel == "emu")
    {
        std::cout << "Single muon trigger fired with double lepton "
                     "trigger/Total single muon triggers fired: "
                  << dupMuon << " / " << singleMuon << std::endl;
    }
}
