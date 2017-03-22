#include "histogramPlotter.hpp"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"
#include "TH1F.h"

// For CMS Guideline styling
#include "TLatex.h"
#include "TASImage.h"

#include <sys/stat.h>
#include <boost/filesystem.hpp>

//For debugging. *sigh*
#include <iostream>

const bool BLIND_PLOTS( true );
const bool writeExtraText( true );

HistogramPlotter::HistogramPlotter(std::vector<std::string> legOrder, std::vector<std::string> plotOrder, std::map<std::string,datasetInfo> dsetMap, const bool is2016):
  //Initialise a load of variables. Labels are empty by default, but this can be changed by calling set label routines.
  lumiStr_{},
  outputFolder_{},
  postfix_{"defaultPostfix"},
  is2016_{is2016},

  //Some things that actually need to be set. plot order, legend order and dataset information map.
  plotOrder_{plotOrder},
  legOrder_{legOrder},
  dsetMap_{dsetMap}
{

  gErrorIgnoreLevel = kInfo;

  extensions_.emplace_back(".root");

  //Make three labels but don't put anything in them just yet. This will be called on the plotting object if we want a label.
  labelOne_ = new TPaveText{0.16,0.88,0.5,0.94,"NDCBR"};
  labelOne_->SetTextAlign(12);
  labelOne_->SetTextSize(0.045);
  labelOne_->SetFillColor(kWhite);
  labelOne_->SetBorderSize(0);

  labelTwo_ = new TPaveText{0.16,0.85,0.5,0.88,"NDCBR"};
  labelTwo_->SetTextAlign(12);
  labelTwo_->SetTextSize(0.045);
  labelTwo_->SetFillColor(kWhite);
  labelTwo_->SetBorderSize(0);

  labelThree_ = new TPaveText{0.16,0.88,0.5,0.94,"NDCBR"};
  labelThree_->SetTextAlign(12);
  labelThree_->SetTextSize(0.045);
  labelThree_->SetFillColor(kWhite);
  labelThree_->SetBorderSize(0);

  setTDRStyle();

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
    stageNameVec.emplace_back(stageNameIt->first);
  }
  //Loop over all the plots, for each stage name. Then create a map for each with all datasets in it.
  unsigned long plotNumb{firstIt->second.begin()->second->getPlotPoint().size()};
  for (unsigned i{0}; i < plotNumb; i++){
    for (auto stageIt = stageNameVec.begin(); stageIt != stageNameVec.end(); stageIt++){
      std::map<std::string, TH1F*> tempPlotMap;
      for (auto mapIt = plotMap.begin(); mapIt != plotMap.end(); mapIt++){
	tempPlotMap[mapIt->first] = mapIt->second[*stageIt]->getPlotPoint()[i].plotHist;
      }
      std::vector<std::string> xAxisLabel = {firstIt->second[*stageIt]->getPlotPoint()[i].xAxisLabel};
      makePlot(tempPlotMap,firstIt->second[*stageIt]->getPlotPoint()[i].title,firstIt->second[*stageIt]->getPlotPoint()[i].name,xAxisLabel);
    }
  }
}

void HistogramPlotter::makePlot(std::map<std::string, TH1F*> plotMap, std::string plotTitle, std::string plotName){
  std::vector<std::string> blankLabels;
  makePlot(plotMap,plotTitle,plotName,"",blankLabels);
}

void HistogramPlotter::makePlot(std::map<std::string, TH1F*> plotMap, std::string plotTitle, std::string plotName, std::vector<std::string> xAxisLabels){
  makePlot(plotMap,plotTitle,plotName,"",xAxisLabels);
}

void HistogramPlotter::makePlot(std::map<std::string, TH1F*> plotMap, std::string plotTitle, std::string plotName, std::string subLabel){
  std::vector<std::string> blankLabels;
  makePlot(plotMap,plotTitle,plotName,subLabel,blankLabels);
}

void HistogramPlotter::makePlot(std::map<std::string, TH1F*> plotMap, std::string plotTitle, std::string plotName, std::string subLabel, std::vector<std::string> xAxisLabels){
  std::cerr << "Making a plot called: " << plotName << std::endl;

  //Make the legend. This is clearly the first thing I should do.
  TLegend* legend_{new TLegend{0.82,0.72,0.98,0.90}};
//  legend_->SetFillStyle(1001);
  legend_->SetFillStyle(0);
  legend_->SetBorderSize(0);
//  legend_->SetFillColor(kWhite);
  for (auto leg_iter = legOrder_.begin(); leg_iter != legOrder_.end(); leg_iter++){
    legend_->AddEntry(plotMap[*leg_iter], dsetMap_[*leg_iter].legLabel.c_str(), dsetMap_[*leg_iter].legType.c_str());
  }

  std::string xAxisLabel {""};
  if ( xAxisLabels.size() < 2 ) xAxisLabel = xAxisLabels[0]; 

  //Initialise the stack
  THStack* mcStack{new THStack{plotName.c_str(),(plotTitle+";"+xAxisLabel+";Events").c_str()}};
  //Do a few colour changing things and add MC to the stack.
  for (auto plot_iter = plotOrder_.rbegin(); plot_iter != plotOrder_.rend(); plot_iter++){
    plotMap[*plot_iter]->SetFillColor(dsetMap_[*plot_iter].colour);
    plotMap[*plot_iter]->SetLineColor(kBlack);
    plotMap[*plot_iter]->SetLineWidth(1);
    if( *plot_iter == "data"){
      plotMap["data"]->SetMarkerStyle(20);
      plotMap["data"]->SetMarkerSize(0.9);
      plotMap["data"]->SetMarkerColor(kBlack);
      continue;
    }
    mcStack->Add(plotMap[*plot_iter]);
  }

  if ( !BLIND_PLOTS ) {
    //Initialise ratio plots
    ratioHisto = dynamic_cast<TH1F*>(plotMap["data"]->Clone());
    ratioHisto->Sumw2();
    ratioHisto->Divide(ratioHisto, dynamic_cast<TH1F*>(mcStack->GetStack()->Last()),1,1 );

    ratioHisto->SetMarkerStyle(20);
    ratioHisto->SetMarkerSize(0.85);
    ratioHisto->SetMarkerColor(kBlack);
  }

  // Set up canvas
  int W = 800;
  int H = 600;
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  std::string canvasName = plotName + subLabel + postfix_;
  int period;
  int pos = 0;
  if( writeExtraText ) canvasName += "-prelim";
  if( !is2016_ ) {
    canvasName += "-2015";
    period = 1;
  }
  else {
    canvasName += "-2016";
    period = 2;
  }

  TCanvas * canvy = new TCanvas((canvasName).c_str(), (canvasName).c_str(), 50, 50, W, H);
  canvy->cd();

  canvy->SetFillColor(0);
  canvy->SetBorderMode(0);
  canvy->SetFrameFillStyle(0);
  canvy->SetFrameBorderMode(0);
  canvy->SetLeftMargin( L/W );
  canvy->SetRightMargin( R/W );
  canvy->SetTopMargin( T/H );
  canvy->SetBottomMargin( B/H );
  canvy->SetTickx(0);
  canvy->SetTicky(0);

  if( !BLIND_PLOTS ) {
    // Top Histogram
    canvy_1 = new TPad("canvy_1", "newpad",0.01,0.315,0.99,0.99);
    canvy_1->Draw();
    canvy_1->cd();
    canvy_1->SetFillColor(0);
    canvy_1->SetBorderMode(0);
    canvy_1->SetFrameFillStyle(0);
    canvy_1->SetFrameBorderMode(0);
    canvy_1->SetLeftMargin( L/W );
    canvy_1->SetRightMargin( R/W );
    canvy_1->SetTopMargin( T/H );
    canvy_1->SetBottomMargin( B/H * 0.3 );
    canvy_1->SetTickx(0);
    canvy_1->SetTicky(0);
  }

  mcStack->Draw("");

  if ( !BLIND_PLOTS ) mcStack->GetHistogram()->GetXaxis()->SetLabelOffset(999);

  setLabelThree("");
  //labelThree_->Draw();
  //  labelTwo_->Draw();
  //  labelOne_->Draw();

  float max = mcStack->GetMaximum();
  if (plotMap.find("data") != plotMap.end()){
    max = TMath::Max(mcStack->GetMaximum(),plotMap["data"]->GetMaximum());
    plotMap["data"]->Draw("e x0, same");
  }

  mcStack->SetMaximum(max*1.1);
  mcStack->GetXaxis()->SetNdivisions(6,5,0);
  mcStack->GetYaxis()->SetNdivisions(6,5,0);
  mcStack->GetYaxis()->SetTitleOffset( L/W * 8.6 );

  if (xAxisLabels.size() > 1){
    for (unsigned i{1}; i <= xAxisLabels.size(); i++){
      if ( !BLIND_PLOTS ) mcStack->GetXaxis()->SetBinLabel(i,"");
      else mcStack->GetXaxis()->SetBinLabel(i,xAxisLabels[i-1].c_str());
    }
  }

  legend_->Draw();

  if ( !BLIND_PLOTS ) {
    // Bottom ratio plots
    canvy->cd();
    canvy_2 = new TPad("canvy_2", "newpad2",0.01,0.01,0.99,0.315);
  //  canvy_2->SetOptStat(0);
    canvy_2->Draw();
    canvy_2->cd();

    canvy_2->SetFillColor(0);
    canvy_2->SetBorderMode(0);
    canvy_2->SetFrameFillStyle(0);
    canvy_2->SetFrameBorderMode(0);
    canvy_2->SetLeftMargin( L/W );
    canvy_2->SetRightMargin( R/W );
    canvy_2->SetTopMargin( T/H );
    canvy_2->SetBottomMargin( B/H );
    canvy_2->SetTickx(0);
    canvy_2->SetTicky(0);
    canvy_2->SetGridy(1);

    ratioHisto->SetStats(false);
    ratioHisto->SetMinimum(0.5);
    ratioHisto->SetMaximum(1.5);
    ratioHisto->SetTitle(";;data/MC");
    ratioHisto->GetXaxis()->SetNdivisions(6,5,0);
    ratioHisto->GetYaxis()->SetNdivisions(6,5,0);

    ratioHisto->GetXaxis()->SetLabelSize(0.12);
    ratioHisto->GetXaxis()->SetLabelOffset(0.018);
    ratioHisto->GetYaxis()->SetLabelSize(0.06);
    ratioHisto->GetYaxis()->SetTitleSize(0.12);
    ratioHisto->GetYaxis()->SetTitleOffset( L/W * 3. );
    ratioHisto->GetYaxis()->CenterTitle();

    if (xAxisLabels.size() > 1){
      for (unsigned i{1}; i <= xAxisLabels.size(); i++){
        ratioHisto->GetXaxis()->SetBinLabel(i,xAxisLabels[i-1].c_str());
      }
    }

    ratioHisto->Draw("e x0, SCAT");
  }

  // writing the lumi information and the CMS "logo"
  if ( BLIND_PLOTS ) {
    CMS_lumi ( canvy, period, pos );
    canvy->Update();
    canvy->RedrawAxis();
    canvy->Draw();
  }
  else {
    CMS_lumi ( canvy_1, period, pos );
    canvy_1->Update();
    canvy_1->RedrawAxis();
    canvy_1->Draw();
  }
  // Save the plots.
  for (unsigned ext_it = 0; ext_it < extensions_.size(); ext_it++){
    canvy->SaveAs((outputFolder_ + plotName + extensions_[ext_it]).c_str());
  }

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

  boost::filesystem::create_directory(outputFolder_.c_str());

}

void HistogramPlotter::CMS_lumi(TPad* pad, int period, int posX) {

  TString cmsText     = "CMS";
  float cmsTextFont   = 61;  // default is helvetic-bold
  TString extraText   = "Preliminary";
  float extraTextFont = 52;  // default is helvetica-italics

  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextSize     = 0.6;
  float lumiTextOffset   = 0.2;
  float cmsTextSize      = 0.75;
  float cmsTextOffset    = 0.1;  // only used in outOfFrame version

  float relPosX    = 0.045;
  float relPosY    = 0.035;
  float relExtraDY = 1.2;

  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;

  TString lumi_2016 = "35.86 fb^{-1}";
  TString lumi_2015 = "2.26 fb^{-1}";
  TString lumi_sqrtS = "";

  bool drawLogo      = false;

  bool outOfFrame    = false;
  if( posX/10==0 ) 
    {
      outOfFrame = true;
    }
  int alignY_=3;
  int alignX_=2;
  if( posX/10==0 ) alignX_=1;
  if( posX==0    ) alignX_=1;
  if( posX==0    ) alignY_=1;
  if( posX/10==1 ) alignX_=1;
  if( posX/10==2 ) alignX_=2;
  if( posX/10==3 ) alignX_=3;
  //if( posX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();

  TString lumiText;
  if( period==1 )
    {
      lumiText += lumi_2015;
      lumiText += " (13 TeV)";
    }
  else if ( period==2 )
    {
      lumiText += lumi_2016;
      lumiText += " (13 TeV)";
    }
  else if ( period==0 )
    {
      lumiText += lumi_sqrtS;
    }
   
//  std::cout << lumiText << std::endl;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  if( outOfFrame )
    {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
    }
  
  pad->cd();

  float posX_=0;
  if( posX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( posX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( posX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
	{
	  posX_ =   l + 0.045*(1-l-r)*W/H;
	  posY_ = 1-t - 0.045*(1-t-b);
	  float xl_0 = posX_;
	  float yl_0 = posY_ - 0.15;
	  float xl_1 = posX_ + 0.15*H/W;
	  float yl_1 = posY_;
	  TASImage* CMS_logo = new TASImage("scripts/CMS-BW-label.png");
	  TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
	  pad_logo->Draw();
	  pad_logo->cd();
	  CMS_logo->Draw("X");
	  pad_logo->Modified();
	  pad->cd();
	}
      else
	{
	  latex.SetTextFont(cmsTextFont);
	  latex.SetTextSize(cmsTextSize*t);
	  latex.SetTextAlign(align_);
	  latex.DrawLatex(posX_, posY_, cmsText);
	  if( writeExtraText ) 
	    {
	      latex.SetTextFont(extraTextFont);
	      latex.SetTextAlign(align_);
	      latex.SetTextSize(extraTextSize*t);
	      latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
	    }
	}
    }
  else if( writeExtraText )
    {
      if( posX==0) 
	{
	  posX_ =   l +  relPosX*(1.75-l-r);
	  posY_ =   1-t+lumiTextOffset*t;
	}
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, extraText);      
    }
}

void HistogramPlotter::setTDRStyle() { 

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

}
