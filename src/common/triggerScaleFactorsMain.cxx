#include "triggerScaleFactorsAlgo.hpp"

#include <iomanip>
#include <iostream>

int main(int argc, char* argv[])
{
    TriggerScaleFactors triggerScaleFactors;

    triggerScaleFactors.parseCommandLineArguements(argc, argv);
    triggerScaleFactors.runMainAnalysis();
    triggerScaleFactors.savePlots();
}
