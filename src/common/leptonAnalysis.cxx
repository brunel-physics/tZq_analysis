#include "AnalysisEvent.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <string>
#include <stdlib.h> 
#include <iostream>

int main(int argc, char* argv[]) {

  std::string config = "";

  // Loop for parsing command line arguments.
  for (int i = 1; i < argc; ++i){
    std::string arg = argv[i];
    if ((arg=="-c")||(arg=="--config")){ // Sets configuration file - Required!
      if (i + 1 < argc) {
      config = argv[++i];
      } else{
	std::cerr << "--config requires an argument!";
	return 0;
      }
    }

  }

  if (config == ""){
    std::cerr << "We need a configuration file! Type -h for usage. Error";
    return 0;
  }

}
