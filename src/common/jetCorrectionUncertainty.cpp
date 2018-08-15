#include "jetCorrectionUncertainty.hpp"

#include "MvaEvent.hpp"
#include "TLorentzVector.h"
#include "TTree.h"
#include "config_parser.hpp"
#include "jetCorrectionUncertainty.hpp"

JetCorrectionUncertainty::JetCorrectionUncertainty(std::string dataFile)
    : ptMinJEC_{}
    , ptMaxJEC_{}
    , etaMinJEC_{}
    , etaMaxJEC_{}
    , jecSFUp_{}
    , jecSFDown_{}

{
    std::ifstream jecFile;
    jecFile.open(dataFile, std::ifstream::in);
    //    jecFile.open(
    //    "scaleFactors/2016/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt",
    //    std::ifstream::in );
    std::string line;
    bool first{true};

    if (!jecFile.is_open())
    {
        std::cout << "Unable to open jecFile." << std::endl;
        exit(0);
    }

    while (getline(jecFile, line))
    {
        std::vector<std::string> tempVec;
        std::stringstream lineStream{line};
        std::string item;
        while (std::getline(lineStream, item, ' '))
        {
            tempVec.emplace_back(item);
        }
        std::vector<float> tempUp;
        std::vector<float> tempDown;

        etaMinJEC_.emplace_back(std::stof(tempVec[0]));
        etaMaxJEC_.emplace_back(std::stof(tempVec[1]));
        for (unsigned i{1}; i < tempVec.size() / 3; i++)
        {
            unsigned ind{i * 3};
            if (first)
            {
                ptMinJEC_.emplace_back(std::stof(tempVec[ind]));
                ptMaxJEC_.emplace_back((ind + 3 >= tempVec.size()
                                            ? 10000.
                                            : std::stof(tempVec[ind + 3])));
            }
            tempUp.emplace_back(std::stof(tempVec[ind + 1]));
            tempDown.emplace_back(std::stof(tempVec[ind + 2]));
        }
        jecSFUp_.emplace_back(tempUp);
        jecSFDown_.emplace_back(tempDown);
        first = false;
    }
}

JetCorrectionUncertainty::~JetCorrectionUncertainty()
{
}

double
    JetCorrectionUncertainty::getUncertainty(double pt, double eta, int jesUD)
{
    if (jesUD == 0)
    {
        return 1.0;
    }

    uint ptBin{0}, etaBin{0};

    for (uint i = 0; i != ptMinJEC_.size(); i++)
    {
        if (pt > ptMinJEC_[i] && pt <= ptMaxJEC_[i])
        {
            ptBin = i;
            break;
        }
    }

    for (uint i = 0; i != etaMinJEC_.size(); i++)
    {
        if (eta > etaMinJEC_[i] && eta <= etaMaxJEC_[i])
        {
            etaBin = i;
            break;
        }
    }

    double lowFact{0.}, highFact{0.};

    if (jesUD == 1)
    {
        lowFact = jecSFUp_[etaBin][ptBin];
        highFact = jecSFUp_[etaBin][ptBin + 1];
    }
    else
    {
        lowFact = jecSFDown_[etaBin][ptBin];
        highFact = jecSFDown_[etaBin][ptBin + 1];
    }

    double a = (highFact - lowFact) / (ptMaxJEC_[ptBin] - ptMinJEC_[ptBin]);
    double b = (lowFact * ptMaxJEC_[ptBin] - highFact * ptMinJEC_[ptBin])
               / (ptMaxJEC_[ptBin] - ptMinJEC_[ptBin]);

    return a * pt + b;
}

std::pair<double, double> JetCorrectionUncertainty::getMetAfterJESUnc(
    double metPx, double metPy, MvaEvent tree, int jesUD)
{
    for (int i = 0; i != tree.numJetPF2PAT; i++)
    {
        metPx += tree.jetPF2PATPx[i];
        metPy += tree.jetPF2PATPy[i];

        double uncertainty =
            getUncertainty(tree.jetPF2PATPt[i], tree.jetPF2PATEta[i], jesUD);

        if (jesUD == 1)
        {
            metPx -= (1 + uncertainty) * tree.jetPF2PATPx[i];
            metPy -= (1 + uncertainty) * tree.jetPF2PATPy[i];
        }
        else
        {
            metPx -= (1 - uncertainty) * tree.jetPF2PATPx[i];
            metPy -= (1 - uncertainty) * tree.jetPF2PATPy[i];
        }
    }
    return std::make_pair(metPx, metPy);
}
