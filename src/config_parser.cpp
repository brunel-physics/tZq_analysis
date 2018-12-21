// config_parser.cpp
#include "config_parser.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <yaml-cpp/yaml.h>

#include "Rtypes.h"

int Parser::parse_config(std::string conf,
                         std::vector<Dataset>* datasets,
                         double* lumi)
{
    const YAML::Node root{YAML::LoadFile(conf)};

    auto datasetConf{root["datasets"].as<std::string>()};

    if (!parse_files(datasetConf, datasets, lumi))
    {
        std::cerr << "Dataset parsing failed!" << std::endl;
        return 0;
    }

    return 1;
}

int Parser::parse_config(std::string conf,
                         std::vector<Dataset>* datasets,
                         double* lumi,
                         std::vector<std::string>* plotTitles,
                         std::vector<std::string>* plotNames,
                         std::vector<float>* xMin,
                         std::vector<float>* xMax,
                         std::vector<int>* nBins,
                         std::vector<std::string>* fillExp,
                         std::vector<std::string>* xAxisLabels,
                         std::vector<int>* cutStage,
                         std::string* cutsConfName,
                         std::string* plotConfName,
                         std::string* outFolder,
                         std::string* postfix,
                         std::string* channel)
{
    const YAML::Node root{YAML::LoadFile(conf)};

    auto datasetConf{root["datasets"].as<std::string>()};

    if (!parse_files(datasetConf, datasets, lumi))
    {
        std::cerr << "Dataset parsing failed!" << std::endl;
        return 0;
    }

    if (*cutsConfName == "")
    {
        *cutsConfName = root["cuts"].as<std::string>();
    }
    if (*plotConfName == "") // If you haven't already chosen the plots to use,
                             // use the ones here.
    {
        *plotConfName = root["plots"].as<std::string>();
    }

    if (!parse_plots(*plotConfName,
                     plotTitles,
                     plotNames,
                     xMin,
                     xMax,
                     nBins,
                     fillExp,
                     xAxisLabels,
                     cutStage))
    {
        std::cerr << "There was a problem parsing the plots" << std::endl;
        return 0;
    }

    if (*outFolder == "plots/" && root["outputFolder"])
    {
        *outFolder = root["outputFolder"].as<std::string>();
    }
    if (*postfix == "default" && root["outputPostfix"])
    {
        *postfix = root["outputPostfix"].as<std::string>();
    }
    if (root["channelName"])
    {
        *channel = root["channelName"].as<std::string>();
    }

    return 1;
}

// For reading the file config.
int Parser::parse_files(std::string fileConf,
                        std::vector<Dataset>* datasets,
                        double* totalLumi)
{
    const YAML::Node root{YAML::LoadFile(fileConf)};
    const std::unordered_map<std::string, int> colourMap{
        {"kAzure", kAzure},
        {"kBlack", kRed},
        {"kBlue", kBlue},
        {"kCyan", kCyan},
        {"kGray", kGray},
        {"kGreen", kGreen},
        {"kMagenta", kMagenta},
        {"kOrange", kOrange},
        {"kPink", kPink},
        {"kPurple", kMagenta + 3},
        {"kRed", kRed},
        {"kRed1", kRed + 5},
        {"kSpring", kSpring},
        {"kTeal", kTeal},
        {"kViolet", kViolet},
        {"kYellow", kYellow}};

    std::cerr << "Adding datasets:" << std::endl;

    for (YAML::const_iterator it = root.begin(); it != root.end(); ++it)
    {
        const bool isMC{it->second["mc"].as<bool>()};

        datasets->emplace_back(
            it->first.as<std::string>(),
            isMC ? 0 : it->second["luminosity"].as<double>(),
            isMC,
            isMC ? it->second["cross_section"].as<double>() : 0,
            it->second["file_list"].as<std::string>(),
            it->second["histogram"].as<std::string>(),
            "tree",
            isMC ? it->second["total_events"].as<long>() : 0,
            colourMap.at(it->second["colour"].as<std::string>()),
            it->second["label"].as<std::string>(),
            it->second["plot_type"].as<std::string>(),
            isMC ? "" : it->second["trigger_flag"].as<std::string>());

        std::cerr << datasets->back().name() << "\t(" << (isMC ? "MC" : "Data")
                  << ')' << std::endl;
    }

    return 1;
}

int Parser::parse_plots(std::string plotConf,
                        std::vector<std::string>* plotTitles,
                        std::vector<std::string>* plotNames,
                        std::vector<float>* xMin,
                        std::vector<float>* xMax,
                        std::vector<int>* nBins,
                        std::vector<std::string>* fillExp,
                        std::vector<std::string>* xAxisLabels,
                        std::vector<int>* cutStage)
{
    const YAML::Node root{YAML::LoadFile(plotConf)};
    const YAML::Node plots{root["plots"]};

    for (YAML::const_iterator it = plots.begin(); it != plots.end(); ++it)
    {
        plotTitles->emplace_back((*it)["title"].as<std::string>());
        plotNames->emplace_back((*it)["name"].as<std::string>());
        xMin->emplace_back((*it)["xMin"].as<float>());
        xMax->emplace_back((*it)["xMax"].as<float>());
        nBins->emplace_back((*it)["nBins"].as<int>());
        fillExp->emplace_back((*it)["fillExp"].as<std::string>());
        xAxisLabels->emplace_back((*it)["xAxisLabel"].as<std::string>());
        cutStage->emplace_back((*it)["cutStage"].as<int>());
    }

    return 1;
}
