#include <iostream>

#include "triggerScaleFactorsAlgo.hpp"

#include <iomanip>

int main(int argc, char* argv[]){

  TriggerScaleFactors triggerScaleFactors;

  triggerScaleFactors.parseCommandLineArguements(argc, argv);
  triggerScaleFactors.runMainAnalysis();
  triggerScaleFactors.savePlots();
}
