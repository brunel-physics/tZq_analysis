#ifndef _plots_hpp_
#define _plots_hpp_

#include "AnalysisEvent.hpp"

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

typedef struct plot plot;

class TH1D;

class Plots
{
    private:
    std::vector<plot> plotPoint;

    public:
    Plots(const std::vector<std::string> titles,
          const std::vector<std::string> names,
          const std::vector<float> xMins,
          const std::vector<float> xMaxs,
          const std::vector<int> nBins,
          const std::vector<std::string> fillExps,
          const std::vector<std::string> xAxisLabels,
          const std::vector<int> cutStage,
          const unsigned thisCutStage,
          const std::string postfixName);
    ~Plots();
    void fillAllPlots(const AnalysisEvent& event, const double eventWeight);
    void saveAllPlots();
    void fillOnePlot(std::string, AnalysisEvent&, float);
    void saveOnePlots(int);
    std::vector<plot> getPlotPoint()
    {
        return plotPoint;
    }
    std::unordered_map<std::string,
                       std::function<std::vector<float>(const AnalysisEvent&)>>
        getFncMap() const;
};

struct plot
{
    std::string name;
    std::string title;
    TH1D* plotHist;
    std::function<std::vector<float>(const AnalysisEvent&)> fillExp;
    std::string xAxisLabel;
    bool fillPlot;
};

#endif // _plots_hpp_ endif
