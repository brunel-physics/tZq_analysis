#include "controlRegions.hpp"
#include <iomanip>

int main(int argc, char* argv[]){

  ControlRegions controlRegions;

  controlRegions.parseCommandLineArguements(argc, argv);
  controlRegions.runMainAnalysis();
  controlRegions.savePlots();
}
