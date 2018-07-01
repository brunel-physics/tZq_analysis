#include <iostream>

#include "makeMVAinputAlgo.hpp"

int main(int argc, char* argv[]){

  MakeMvaInputs makeMvaInputs;

  makeMvaInputs.parseCommandLineArguements(argc, argv);
  makeMvaInputs.runMainAnalysis();

}


