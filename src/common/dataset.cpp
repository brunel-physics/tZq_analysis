#include "dataset.hpp"
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"

Dataset::Dataset(std::string name, float lumi, bool isMC, float crossSection, std::string fileList, std::string histoName, std::string treeName, long long totalEvents, int colourInt, std::string plotLabel, std::string plotType, std::string triggerFlag){
  name_ = name;
  lumi_ = lumi;
  isMC_ = isMC;
  crossSection_ = crossSection;
  //TODO: Do something with the fileList here. This will build the TChain.
  fillName_ = histoName;
  treeName_= treeName;
  fileList_ = fileList;
  totalEvents_ = totalEvents;
  colour_ = colourInt;
  plotType_ = plotType;
  plotLabel_ = plotLabel;
  triggerFlag_ = triggerFlag;
  std::cout << "For dataset " << name_ << " trigger flag is " << triggerFlag_ << std::endl;
}

//Method that fills a TChain with the files that will be used for the analysis. Returns 1 if succesful, otherwise returns 0. This can probably be largely ignored.
int Dataset::fillChain(TChain * chain , int nFiles){
  std::cerr <<  fileList_ << "\n";
  std::ifstream fileList(fileList_);
  if (!fileList.is_open()){
    std::cerr << "Couldn't read file list for " << name_ <<  std::endl;
    return 0;
  }
  std::string line;
  int files{0};
  while(getline(fileList,line)){
    chain->Add(line.c_str());
    if (nFiles > 0) {
      files ++;
      if (files > nFiles) break;
    }
  }
  return 1;
}

// Function that returns the weight of a dataset. This is 1 is the dataset is data, but varies by lumi, number of events and cross section otherwise.
float Dataset::getDatasetWeight(double lumi){
  if (!isMC_) return 1.;
  return (lumi * crossSection_)/totalEvents_;
}

//Function to weight each event in a dataset
float Dataset::getEventWeight(){
  return 1.;

}

//Function that constructs a histogram of all the generator level weights from across the entire dataset
TH1I* Dataset::getGeneratorWeightHistogram( int nFiles ){
//  std::cerr <<  fileList_ << "\n";
  std::ifstream fileList(fileList_);

  TH1I* generatorWeightPlot{nullptr};

  if (!fileList.is_open()){
    std::cerr << "Couldn't read file list for " << name_ <<  std::endl;
    return generatorWeightPlot;
  }

  std::string line;
  int files{0};
  bool firstFile{true};

  while(getline(fileList,line)) {
    TFile* tempFile {new TFile {line.c_str(), "READ"} };
    if (firstFile){
      generatorWeightPlot = dynamic_cast<TH1I*>( (TH1I*)tempFile->Get("sumNumPosMinusNegWeights")->Clone() );
      firstFile = false;
    }
    else {
      generatorWeightPlot->Add( (TH1I*)tempFile->Get("sumNumPosMinusNegWeights") );
      tempFile->Close();
      delete tempFile;
    }
    if (nFiles > 0) {
      files ++;
      if (files > nFiles) break;
    }
  }
  return generatorWeightPlot;

}
