// config_parser.cpp
#include "config_parser.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <yaml-cpp/yaml.h>

void Parser::parse_config(const std::string conf,
                          std::vector<Dataset>& datasets,
                          double& lumi)
{
    const YAML::Node root{YAML::LoadFile(conf)};
    auto datasetConfs{root["datasets"].as<std::vector<std::string>>()};
    try
    {
        parse_files(datasetConfs, datasets, lumi);
    }
    catch (const std::exception)
    {
        std::cerr << "ERROR while parsing dataset file" << std::endl;
        throw;
    }
}

void Parser::parse_config(const std::string conf,
                          std::vector<Dataset>& datasets,
                          double& lumi,
                          std::vector<std::string>& plotTitles,
                          std::vector<std::string>& plotNames,
                          std::vector<float>& xMin,
                          std::vector<float>& xMax,
                          std::vector<int>& nBins,
                          std::vector<std::string>& fillExp,
                          std::vector<std::string>& xAxisLabels,
                          std::vector<int>& cutStage,
                          std::string& cutsConfName,
                          std::string& plotConfName,
                          std::string& outFolder,
                          std::string& postfix,
                          std::string& channel)
{
    const YAML::Node root{YAML::LoadFile(conf)};

    auto datasetConfs{root["datasets"].as<std::vector<std::string>>()};

    try
    {
        parse_files(datasetConfs, datasets, lumi);
    }
    catch (const std::exception)
    {
        std::cerr << "ERROR while parsing dataset file" << std::endl;
        throw;
    }

    if (cutsConfName == "")
    {
        cutsConfName = root["cuts"].as<std::string>();
    }
    if (plotConfName == "") // If you haven't already chosen the plots to use,
                             // use the ones here.
    {
        plotConfName = root["plots"].as<std::string>();
    }

    try
    {
        parse_plots(plotConfName,
                    plotTitles,
                    plotNames,
                    xMin,
                    xMax,
                    nBins,
                    fillExp,
                    xAxisLabels,
                    cutStage);
    }
    catch (const std::exception)
    {
        std::cerr << "ERROR while parsing plot configuration" << std::endl;
        throw;
    }

    if (outFolder == "plots/" && root["outputFolder"])
    {
        outFolder = root["outputFolder"].as<std::string>();
    }
    if (postfix == "default" && root["outputPostfix"])
    {
        postfix = root["outputPostfix"].as<std::string>();
    }
    if (root["channelName"])
    {
        channel = root["channelName"].as<std::string>();
    }
}

// For reading the file config.
void Parser::parse_files(const std::vector<std::string> files,
                         std::vector<Dataset>& datasets,
                         double& totalLumi)
{
    std::cerr << "Adding datasets:" << std::endl;

    for (const auto& file: files)
    {
        const YAML::Node root{YAML::LoadFile(file)};
        const bool isMC{root["mc"].as<bool>()};
        datasets.emplace_back(root["name"].as<std::string>(),
                              isMC ? 0 : root["luminosity"].as<double>(),
                              isMC,
                              isMC ? root["cross_section"].as<double>() : 0,
                              root["locations"].as<std::vector<std::string>>(),
                              root["histogram"].as<std::string>(),
                              "tree",
                              isMC ? root["total_events"].as<long>() : 0,
                              root["colour"].as<std::string>(),
                              root["label"].as<std::string>(),
                              root["plot_type"].as<std::string>(),
                              isMC ? ""
                                   : root["trigger_flag"].as<std::string>());

        if (root["luminosity"])
        {
            totalLumi += root["luminosity"].as<double>();
        }

        std::cerr << datasets.back().name() << "\t(" << (isMC ? "MC" : "Data")
                  << ')' << std::endl;
    }
}

void Parser::parse_plots(const std::string plotConf,
                         std::vector<std::string>& plotTitles,
                         std::vector<std::string>& plotNames,
                         std::vector<float>& xMin,
                         std::vector<float>& xMax,
                         std::vector<int>& nBins,
                         std::vector<std::string>& fillExp,
                         std::vector<std::string>& xAxisLabels,
                         std::vector<int>& cutStage)
{
    const YAML::Node root{YAML::LoadFile(plotConf)};
    const YAML::Node plots{root["plots"]};

    for (YAML::const_iterator it = plots.begin(); it != plots.end(); ++it)
    {
        plotTitles.emplace_back((*it)["title"].as<std::string>());
        plotNames.emplace_back((*it)["name"].as<std::string>());
        xMin.emplace_back((*it)["xMin"].as<float>());
        xMax.emplace_back((*it)["xMax"].as<float>());
        nBins.emplace_back((*it)["nBins"].as<int>());
        fillExp.emplace_back((*it)["fillExp"].as<std::string>());
        xAxisLabels.emplace_back((*it)["xAxisLabel"].as<std::string>());
        cutStage.emplace_back((*it)["cutStage"].as<int>());
    }
}
