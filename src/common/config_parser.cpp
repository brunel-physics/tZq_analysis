//config_parser.cpp
#include "config_parser.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <libconfig.h++>

int Parser::parse_config(std::string conf, std::vector<Dataset> * datasets, double * lumi,std::vector<std::string>*plotNames,std::vector<float>*xMin,std::vector<float>*xMax,std::vector<int>*nBins,std::vector<std::string>*fillExp,std::vector<std::string>*xAxisLabels,std::vector<int>*cutStage,std::string* cutsConfName,std::string* plotConfName, std::string* outFolder, std::string* postfix, std::string* channel){
  //Re-write config parser here.
  libconfig::Config config;

  try{ //Attempt to open the configuration file.
    config.readFile(conf.c_str());
  }
  catch (const libconfig::FileIOException &exep){ //No file error
    std::cerr << "Error opening file" << std::endl;
    return 0;
  }
  catch(const libconfig::ParseException &e){ //Parsing the file error
    std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine() << " - " << e.getError() << std::endl;
    return 0;
  }
  libconfig::Setting &root = config.getRoot();
  std::string datasetConf;
  root.lookupValue("datasets",datasetConf);
  std::cerr << root.exists("datasets");
  if (!parse_files(datasetConf,datasets,lumi)){
    std::cerr << "Dataset parsing failed!" << std::endl;
    return 0;
  }
  if (*cutsConfName == "")
    root.lookupValue("cuts",*cutsConfName);
  if (*plotConfName == "") // If you haven't already chosen the plots to use, use the ones here.
    root.lookupValue("plots",*plotConfName);
  if (!parse_plots(*plotConfName,plotNames,xMin,xMax,nBins,fillExp,xAxisLabels,cutStage)){
    std::cerr << "There was a problem parsing the plots" << std::endl;
    return 0;
  }
  if (*outFolder == "plots/" && root.exists("outputFolder")) root.lookupValue("outputFolder",*outFolder);
  if (*postfix == "default" && root.exists("outputPostfix")) root.lookupValue("outputPostfix",*postfix);
  if (root.exists("channelName")) root.lookupValue("channelName",*channel);
  //Succesfully parsed everything.
  return 1;
}

//For reading the file config. 
int Parser::parse_files(std::string fileConf, std::vector<Dataset> * datasets, double * totalLumi){
  std::cout << "file config: " << fileConf << std::endl;
  std::ifstream conf(fileConf);
  if (!conf.is_open()){
    std::cerr << "Couldn't open file config! Exiting!\n";
    return 0;
  }
  std::map<std::string,int> colourMap = getColourMap();
  std::string line;
  //Main config loop. Every time it finds a new section name it should add a dataset to the datasets vector passed.
  std::string datasetName= "";
  std::string file = "";
  bool isMC = false;
  long totalEvents = 0;
  std::string histoFill = "";
  double crossSection = 0.;
  double lumi = 0.;
  std::string treeLabel = "tree";
  std::string plotLabel = "";
  std::string plotType = "";
  int colourInt = 0;
  std::string triggerFlag = "";
  std::cerr << "Adding datasets: \n";
  while (getline(conf,line)){
    if (line.find("[")==0) {
      if (datasetName != ""){
	std::cerr << "  " << datasetName << std::endl;
	std::cerr << "   " << isMC << "\t" << file << "\n   " << crossSection << "\t" << totalEvents << "\t" << histoFill << std::endl;
	if (isMC) datasets->push_back(Dataset(datasetName, 0.,isMC,crossSection,file,histoFill,treeLabel,totalEvents,colourInt,plotLabel,plotType,triggerFlag));
	else datasets->push_back(Dataset(datasetName, lumi,isMC,0.,file,histoFill,treeLabel,totalEvents,colourInt,plotLabel,plotType,triggerFlag));
      }
      treeLabel = "tree";
      triggerFlag = "";
      datasetName = line.substr(1,line.find("]")-1);
    }
    //Add in file name that will be used to construct TChain
    if (line.find("fileName") == 0){
      file = line.substr(9);
      file = file.substr(file.find_first_not_of(" ="));
    }
    //Add lumi (if data)
    if (line.find("luminosity") == 0){
      std::string temp = line.substr(11);
      temp = temp.substr(temp.find_first_not_of(" ="));
      //      lumi = temp.c_str().std::stof();
      lumi = std::stof(temp);
      *totalLumi+=lumi;
    }
    //Get cross section info (if mc)
    if (line.find("crossSection") == 0){
      std::string temp = line.substr(13);
      temp = temp.substr(temp.find_first_not_of(" ="));
      crossSection = std::stof(temp);
    }
    //MC or data
    if (line.find("runType") == 0){
      if (line.find("mc") < 20) isMC = true;
      else { isMC = false;
	crossSection = 0.;
      }
    }
    //get total events in MC
    if (line.find("totalEvents") == 0){
      std::string temp = line.substr(12);
      temp = temp.substr(temp.find_first_not_of(" ="));
      //      std::cout << temp;
      totalEvents = std::stol(temp);
    }
    if (line.find("histoName") == 0){
      histoFill = line.substr(10);
      histoFill = histoFill.substr(histoFill.find_first_not_of(" ="));
    }
    if (line.find("treeName") == 0){
      treeLabel = line.substr(9);
      treeLabel = treeLabel.substr(treeLabel.find_first_not_of(" ="));
    }
    if (line.find("colour") == 0){
      std::string tempString = line.substr(7);
      tempString = tempString.substr(tempString.find_first_not_of(" ="));
      colourInt = colourMap[tempString];
    }
    if (line.find("label") == 0){
      plotLabel = line.substr(6);
      plotLabel = plotLabel.substr(plotLabel.find_first_not_of(" ="));
    }
    if (line.find("plotType") == 0){
      plotType = line.substr(9);
      plotType = plotType.substr(plotType.find_first_not_of(" ="));
    }
    if (line.find("triggerFlag") == 0){
      triggerFlag = line.substr(12);
      triggerFlag = triggerFlag.substr(triggerFlag.find_first_not_of(" ="));
    }
    if (line.find("pileupDistro") == 0){
      
    }
      
    
    
  }
  conf.close();
  return 1;
}

int Parser::parse_plots(std::string plotConf,std::vector<std::string> *plotNames,std::vector<float> *xMin,std::vector<float> *xMax,std::vector<int> *nBins,std::vector<std::string> *fillExp,std::vector<std::string> *xAxisLabels, std::vector<int>*cutStage){
  //make the object
  libconfig::Config plotCfg;

  //apparently this is done inside a try catch loop. I don't really like this, but the example does it so...
  try {
    plotCfg.readFile(plotConf.c_str());
  }
  catch (const libconfig::FileIOException &exep){
    std::cerr << "Error opening file" << std::endl;
    return 0;
  }
  catch(const libconfig::ParseException &e){
    std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine() << " - " << e.getError() << std::endl;
    return 0;
  }
  const libconfig::Setting& root = plotCfg.getRoot();
  
  if (! root.exists("plots")){
    std::cerr << "No plots setting in the config file! What the hell are you doing?!?" << std::endl;
    return 0;
  }
  libconfig::Setting& plots = root["plots"];
  std::string nameT,fillExpT,xAxisLabelT;
  float xMinT,xMaxT;
  int nBinsT,cutStageT;
  for (int i = 0; i < plots.getLength(); i++){
    const libconfig::Setting &plot = plots[i];
    plot.lookupValue("name",nameT);
    plotNames->push_back(nameT);
    plot.lookupValue("xMin",xMinT);
    xMin->push_back(xMinT);
    plot.lookupValue("xMax",xMaxT);
    xMax->push_back(xMaxT);
    plot.lookupValue("nBins",nBinsT);
    nBins->push_back(nBinsT);
    plot.lookupValue("fillExp",fillExpT);
    fillExp->push_back(fillExpT);
    plot.lookupValue("xAxisLabel",xAxisLabelT);
    xAxisLabels->push_back(xAxisLabelT);
    plot.lookupValue("cutStage",cutStageT);
    cutStage->push_back(cutStageT);
  }


  return 1;
}


std::map<std::string,int> Parser::getColourMap(){
  std::map<std::string,int> colourMap;
  colourMap["kRed"] = kRed;
  colourMap["kBlue"] = kBlue;
  colourMap["kGreen"] = kGreen;
  colourMap["kYellow"] = kYellow;
  colourMap["kTeal"] = kTeal;
  colourMap["kMagenta"] = kMagenta;
  colourMap["kCyan"] = kCyan;
  colourMap["kAzure"] = kAzure;
  colourMap["kSpring"] = kSpring;
  colourMap["kPink"] = kPink;
  colourMap["kOrange"] = kOrange;
  colourMap["kRed1"] = kRed + 5;
  colourMap["kBlack"] = kBlack;
  return colourMap;
}
