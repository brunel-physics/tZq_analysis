#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h> 
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"

#include "analysisAlgo.hpp"
#include <libconfig.h++>

#include <iomanip>
#include <map>
#include <math.h>
#include <LHAPDF/LHAPDF.h>

int main(int argc, char* argv[]){

  analysisAlgo analysisMain;

  gErrorIgnoreLevel = kInfo;
  //Set up environment a little.
  std::cerr << std::setprecision(1) << std::fixed;
  std::cout << std::setprecision(1) << std::fixed;
  // "This is the main function. It basically just loads a load of other stuff.";
  //Parse command line arguments - looking for config file.
  if (argc < 3){
    show_usage(argv[0]);
    return 1;
  }

  analysisMain.parseCommandLineArguements(argc, argv);

  analysisMain.setupSystematics();
  analysisMain.setupCuts();
  analysisMain.setupPlots();
  analysisMain.runMainAnalysis();
  analysisMain.savePlots();
}

