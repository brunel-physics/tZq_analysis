#include "dataset.hpp"

#include "TChain.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <fstream>
#include <iostream>
#include <regex>

namespace fs = boost::filesystem;

Dataset::Dataset(std::string name,
                 float lumi,
                 bool isMC,
                 float crossSection,
                 std::vector<std::string> locations,
                 std::string histoName,
                 std::string treeName,
                 long long totalEvents,
                 std::string colourHex,
                 std::string plotLabel,
                 std::string plotType,
                 std::string triggerFlag)
    : colour_{TColor::GetColor(colourHex.c_str())}
{
    name_ = name;
    lumi_ = lumi;
    isMC_ = isMC;
    crossSection_ = crossSection;
    // TODO: Do something with the fileList here. This will build the TChain.
    fillName_ = histoName;
    treeName_ = treeName;
    locations_ = locations;
    totalEvents_ = totalEvents;
    plotType_ = plotType;
    plotLabel_ = plotLabel;
    triggerFlag_ = triggerFlag;
    std::cout << "For dataset " << name_ << " trigger flag is " << triggerFlag_
              << std::endl;

    for (auto& location : locations_)
    {
        if (location.back() != '/')
        {
            location += '/';
        }
    }
}

// Method that fills a TChain with the files that will be used for the analysis.
// Returns 1 if succesful, otherwise returns 0. This can probably be largely
// ignored.
int Dataset::fillChain(TChain* chain, int nFiles)
{
    for (const auto& location : locations_)
    {
        const fs::path dir{location};
        if (fs::is_directory(dir))
        {
            chain->Add(TString{location + "*.root"});
        }
        else
        {
            std::cout << "ERROR: " << location << "is not a valid directory"
                      << std::endl;
            return 0;
        }
    }
    return 1;
}

// Function that returns the weight of a dataset. This is 1 is the dataset is
// data, but varies by lumi, number of events and cross section otherwise.
float Dataset::getDatasetWeight(double lumi)
{
    if (!isMC_)
    {
        return 1.;
    }
    return (lumi * crossSection_) / totalEvents_;
}

// Function to weight each event in a dataset
float Dataset::getEventWeight()
{
    return 1.;
}

// Function that constructs a histogram of all the generator level weights from
// across the entire dataset
TH1I* Dataset::getGeneratorWeightHistogram(int nFiles)
{
    TH1I* generatorWeightPlot{nullptr};
    const std::regex mask{R"(\.root$)"};
    bool firstFile{true};
    for (const auto& location : locations_)
    {
        for (const auto& file :
             boost::make_iterator_range(fs::directory_iterator{location}, {}))
        {
            const std::string path{file.path().string()};
            if (!fs::is_regular_file(file.status())
                || !std::regex_search(path, mask))
            {
                continue;
            }

            TFile* tempFile{new TFile{path.c_str(), "READ"}};
            if (firstFile)
            {
                generatorWeightPlot = dynamic_cast<TH1I*>(
                    (TH1I*)tempFile->Get("sumNumPosMinusNegWeights")->Clone());
                firstFile = false;
            }
            else
            {
                generatorWeightPlot->Add(
                    (TH1I*)tempFile->Get("sumNumPosMinusNegWeights"));
                tempFile->Close();
                delete tempFile;
            }
        }
    }

    return generatorWeightPlot;
}
