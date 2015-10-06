//config_parser.hpp

#ifndef _config_parser_hpp_
#define _config_parser_hpp_

#include <string>
#include "dataset.hpp"
#include <map>
#include "TColor.h"
#include <vector>

namespace Parser {
  int parse_config(std::string conf,std::vector<Dataset> * datasets,double* lumi,std::vector<std::string>*,std::vector<float>*,std::vector<float>*,std::vector<int>*,std::vector<std::string>*,std::vector<std::string>*,std::vector<int>*,std::string*,std::string*,std::string*,std::string*,std::string*);
  int parse_files(std::string FileConf, std::vector<Dataset> * datasets,double* lumi);
  int parse_plots(std::string plotConf,std::vector<std::string>*,std::vector<float>*,std::vector<float>*,std::vector<int>*,std::vector<std::string>*,std::vector<std::string>*,std::vector<int>*);
  std::map<std::string,int> getColourMap();
}

#endif
