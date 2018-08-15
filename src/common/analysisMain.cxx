#include "analysisAlgo.hpp"

#include <iomanip>
#include <iostream>

int main(int argc, char* argv[])
{
    AnalysisAlgo analysisMain;

    analysisMain.parseCommandLineArguements(argc, argv);
    analysisMain.setupSystematics();
    analysisMain.setupCuts();
    analysisMain.setupPlots();
    analysisMain.runMainAnalysis();
    analysisMain.savePlots();
}
