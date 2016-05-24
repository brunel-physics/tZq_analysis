#include "histogramPlotter.hpp"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include <sys/stat.h>

//For debugging. *sigh*
#include <iostream>


HistogramPlotter::HistogramPlotter(std::vector<std::string> legOrder, std::vector<std::string> plotOrder, std::map<std::string,datasetInfo> dsetMap):
  //Initialise a load of variables. Labels are empty by default, but this can be changed by calling set label routines.
  lumiStr_(""),
  outputFolder_(""),
  postfix_("defaultPostfix"),
  

  //Some things that actually need to be set. plot order, legend order and dataset information map.
  plotOrder_(plotOrder),
  legOrder_(legOrder),
  dsetMap_(dsetMap)
{

  gErrorIgnoreLevel = kInfo;

  extensions_.push_back(".root");
  //  extensions_.push_back(".png");

  //Make three labels but don't put anything in them just yet. This will be called on the plotting object if we want a label.
  labelOne_ = new TPaveText(0.16,0.88,0.5,0.94,"NDCBR");
  labelOne_->SetTextAlign(12);
  labelOne_->SetTextSize(0.045);
  labelOne_->SetFillColor(kWhite);
  labelOne_->SetBorderSize(0);

  labelTwo_ = new TPaveText(0.16,0.85,0.5,0.88,"NDCBR");
  labelTwo_->SetTextAlign(12);
  labelTwo_->SetTextSize(0.045);
  labelTwo_->SetFillColor(kWhite);
  labelTwo_->SetBorderSize(0);

  labelThree_ = new TPaveText(0.16,0.88,0.5,0.94,"NDCBR");
  labelThree_->SetTextAlign(12);
  labelThree_->SetTextSize(0.045);
  labelThree_->SetFillColor(kWhite);
  labelThree_->SetBorderSize(0);

  gStyle->SetLabelFont(18,"");
  
  //for (auto leg_iter = legOrder.begin(); leg_iter != legOrder.end(); leg_iter++){
  //  
  //}
  //I guess I should put in the plot style stuff here? Worry about that later on. Get it plotting something at all first I guess.
}

HistogramPlotter::~HistogramPlotter(){
  delete labelOne_;
  delete labelTwo_;
  delete labelThree_;
}

void HistogramPlotter::plotHistos(std::map<std::string, std::map<std::string, Plots*> > plotMap){
  //Get a list of keys from the map.
  auto firstIt = plotMap.begin();
  std::vector<std::string> stageNameVec;
  for (auto stageNameIt = firstIt->second.begin(); stageNameIt != firstIt->second.end(); stageNameIt++){
    stageNameVec.push_back(stageNameIt->first);
  }
  //Loop over all the plots, for each stage name. Then create a map for each with all datasets in it.
  int plotNumb = firstIt->second.begin()->second->getPlotPoint().size();
  for (int i = 0; i < plotNumb; i++){
    for (auto stageIt = stageNameVec.begin(); stageIt != stageNameVec.end(); stageIt++){
      std::map<std::string, TH1F*> tempPlotMap;
      for (auto mapIt = plotMap.begin(); mapIt != plotMap.end(); mapIt++){
	tempPlotMap[mapIt->first] = mapIt->second[*stageIt]->getPlotPoint()[i].plotHist;
      }
      makePlot(tempPlotMap,firstIt->second[*stageIt]->getPlotPoint()[i].name,"");
    }
  }
}

void HistogramPlotter::makePlot(std::map<std::string, TH1F*> plotMap, std::string plotName){
  std::vector<std::string> blankLabels;
  makePlot(plotMap,plotName,"",blankLabels);
}

void HistogramPlotter::makePlot(std::map<std::string, TH1F*> plotMap, std::string plotName, std::vector<std::string> xAxisLabels){
  makePlot(plotMap,plotName,"",xAxisLabels);
}

void HistogramPlotter::makePlot(std::map<std::string, TH1F*> plotMap, std::string plotName, std::string subLabel){
  std::vector<std::string> blankLabels;
  makePlot(plotMap,plotName,subLabel,blankLabels);
}  

void HistogramPlotter::makePlot(std::map<std::string, TH1F*> plotMap, std::string plotName, std::string subLabel, std::vector<std::string> xAxisLabels){
  std::cerr << "Making a plot called: " << plotName << std::endl;
  //Make the legend. This is clearly the first thing I should do.
  TLegend* legend_ = new TLegend(0.7,0.7,0.94,0.94);
  legend_->SetFillStyle(1001);
  legend_->SetBorderSize(1);
  legend_->SetFillColor(kWhite);
  for (auto leg_iter = legOrder_.begin(); leg_iter != legOrder_.end(); leg_iter++){
    legend_->AddEntry(plotMap[*leg_iter], dsetMap_[*leg_iter].legLabel.c_str(), dsetMap_[*leg_iter].legType.c_str());
  }
    
  //Initialise the stack
  THStack* mcStack = new THStack(plotName.c_str(),plotName.c_str());
  //Do a few colour changing things and add MC to the stack.
  for (auto plot_iter = plotOrder_.rbegin(); plot_iter != plotOrder_.rend(); plot_iter++){
    plotMap[*plot_iter]->SetFillColor(dsetMap_[*plot_iter].colour);
    plotMap[*plot_iter]->SetLineColor(kBlack);
    plotMap[*plot_iter]->SetLineWidth(1);
    if( *plot_iter == "data"){
      plotMap["data"]->SetMarkerStyle(20);
      plotMap["data"]->SetMarkerSize(1.2);
      plotMap["data"]->SetMarkerColor(kBlack);
      continue;
    }
    mcStack->Add(plotMap[*plot_iter]);
  }

  //Initialise ratio plots
  TH1F* ratioHisto = (TH1F*) plotMap["data"]->Clone();
  ratioHisto->GetTitle();
  ratioHisto->Divide( (TH1F*)(mcStack->GetStack()->Last()) );
  ratioHisto->SetMarkerStyle(20);
  ratioHisto->SetMarkerSize(1.2);
  ratioHisto->SetMarkerColor(kBlack);

  // Set up canvas
  TCanvas * canvy = new TCanvas((plotName + subLabel + postfix_).c_str(), (plotName + subLabel + postfix_).c_str());
  canvy->cd();

  // Top Histogram
  TPad* canvy_1 = new TPad("canvy_1", "newpad",0.01,0.255,0.99,0.99);
  canvy_1->Draw(); 
  canvy_1->cd();
  canvy_1->SetTopMargin(0.1);
  canvy_1->SetBottomMargin(0.01);
  canvy_1->SetRightMargin(0.1);
  canvy_1->SetFillStyle(0);

  mcStack->Draw("");

  if (xAxisLabels.size() > 0){
    for (unsigned i = 1; i <= xAxisLabels.size(); i++){
      mcStack->GetXaxis()->SetBinLabel(i,xAxisLabels[i-1].c_str());
    }
  }
  
  setLabelThree(subLabel);
  //labelThree_->Draw();
  //  labelTwo_->Draw();
  //  labelOne_->Draw();

  float max = mcStack->GetMaximum();
  if (plotMap.find("data") != plotMap.end()){
    max = TMath::Max(mcStack->GetMaximum(),plotMap["data"]->GetMaximum());
    plotMap["data"]->Draw("e x0, same");
  }

  mcStack->SetMaximum(max*1.3);
      
  legend_->Draw();
  
  // Bottom ratio plots
  canvy->cd();
  TPad* canvy_2 = new TPad("canvy_2", "newpad2",0.01,0.01,0.99,0.245);
  canvy_2->Draw(); 
  canvy_2->cd();
  canvy_2->SetTopMargin(0.01);
  canvy_2->SetBottomMargin(0.3);
  canvy_2->SetRightMargin(0.1);
  canvy_2->SetFillStyle(0);
  
  ratioHisto->SetMinimum(0.5);
  ratioHisto->SetMaximum(1.5);
  ratioHisto->Draw("e x0");

  // Save the plots.
  for (unsigned ext_it = 0; ext_it < extensions_.size(); ext_it++){
    canvy->SaveAs((outputFolder_ + plotName + subLabel + postfix_ + extensions_[ext_it]).c_str());
  }

/*  //Make log plots
  canvy->SetLogy();
  mcStack.SetMaximum(max*10);

  //Save the log plots
  for (unsigned ext_it = 0; ext_it < extensions_.size(); ext_it++){
    canvy->SaveAs((outputFolder_ + plotName + subLabel + postfix_ + "_log" + extensions_[ext_it]).c_str());
  }*/
  
  delete canvy;
}

//Resets the label size on the plot.
void HistogramPlotter::setLabelTextSize(float size){
  labelOne_->SetTextSize(size);
  labelTwo_->SetTextSize(size);
  labelThree_->SetTextSize(size);
}

void HistogramPlotter::setOutputFolder(std::string output){
  outputFolder_ = output;

  mkdir(outputFolder_.c_str(),0700);

}
