#ifndef _dataset_hpp_
#define _dataset_hpp_

#include <string>
#include <vector>

class TChain;
class TH1I;

// class to hold information about the dataset. This will be extracted during
// config parsing.
class Dataset
{
    std::string name_;
    float lumi_;
    bool isMC_;
    float crossSection_;
    std::string fillName_;
    std::string treeName_;
    std::vector<std::string> locations_;
    long long totalEvents_;
    const int colour_;
    std::string plotLabel_;
    std::string plotType_;
    std::string triggerFlag_;

    public:
    Dataset(std::string name,
            float lumi,
            bool isMC,
            float crossSection,
            std::vector<std::string> locations,
            std::string histoName,
            std::string treeName,
            long long,
            std::string,
            std::string,
            std::string,
            std::string);
    std::string name()
    {
        return name_;
    }
    float lumi()
    {
        return lumi_;
    }
    bool isMC()
    {
        return isMC_;
    }
    std::string getFillHisto()
    {
        return fillName_;
    }
    std::string treeName()
    {
        return treeName_;
    }
    int getColour()
    {
        return colour_;
    }
    std::string getPlotLabel()
    {
        return plotLabel_;
    }
    std::string getPlotType()
    {
        return plotType_;
    }
    int fillChain(TChain*, int);
    float getDatasetWeight(double);
    float getEventWeight();
    std::string getTriggerFlag()
    {
        return triggerFlag_;
    }
    long long getTotalEvents()
    {
        return totalEvents_;
    } // Function to return total number of events
    TH1I* getGeneratorWeightHistogram(int nFiles);
};

#endif
