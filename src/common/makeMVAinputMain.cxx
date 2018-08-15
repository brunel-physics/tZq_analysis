#include "makeMVAinputAlgo.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
    MakeMvaInputs makeMvaInputs;

    makeMvaInputs.parseCommandLineArguements(argc, argv);
    makeMvaInputs.runMainAnalysis();
}
