// config_parser.hpp

#ifndef _config_parser_hpp_
#define _config_parser_hpp_

#include "TColor.h"
#include "dataset.hpp"

#include <map>
#include <string>
#include <vector>

namespace Parser
{
    void parse_config(const std::string conf,
                      std::vector<Dataset>& datasets,
                      double& lumi);
    void parse_config(const std::string conf,
                      std::vector<Dataset>& datasets,
                      double& lumi,
                      std::vector<std::string>&,
                      std::vector<std::string>&,
                      std::vector<float>&,
                      std::vector<float>&,
                      std::vector<int>&,
                      std::vector<std::string>&,
                      std::vector<std::string>&,
                      std::vector<int>&,
                      std::string&,
                      std::string&,
                      std::string&,
                      std::string&,
                      std::string&);
    void parse_files(const std::vector<std::string> files,
                     std::vector<Dataset>& datasets,
                     double& lumi);
    void parse_plots(const std::string plotConf,
                     std::vector<std::string>&,
                     std::vector<std::string>&,
                     std::vector<float>&,
                     std::vector<float>&,
                     std::vector<int>&,
                     std::vector<std::string>&,
                     std::vector<std::string>&,
                     std::vector<int>&);
} // namespace Parser

#endif
