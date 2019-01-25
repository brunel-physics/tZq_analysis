#include "analysisAlgo.hpp"

#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <limits>

int main(int argc, char* argv[])
{
    TTree::SetMaxTreeSize(std::numeric_limits<long long>::max());

    AnalysisAlgo analysisMain;

    analysisMain.parseCommandLineArguements(argc, argv);
    analysisMain.setupSystematics();
    analysisMain.setupCuts();
    analysisMain.setupPlots();
    analysisMain.runMainAnalysis();
    analysisMain.savePlots();
}
