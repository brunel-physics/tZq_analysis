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

    if (!jecFile.is_open())
    {
        std::cout << "Unable to open jecFile." << std::endl;
        exit(0);
    }

    bool first{true};
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
            const unsigned ind{i * 3};
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

double JetCorrectionUncertainty::getUncertainty(const double pt,
                                                const double eta,
                                                const int jesUD) const
{
    if (jesUD == 0)
    {
        return 1.0;
    }

    unsigned ptBin{0};
    unsigned etaBin{0};

    for (size_t i{0}; i != ptMinJEC_.size(); i++)
    {
        if (pt > ptMinJEC_[i] && pt <= ptMaxJEC_[i])
        {
            ptBin = i;
            break;
        }
    }

    for (size_t i{0}; i != etaMinJEC_.size(); i++)
    {
        if (eta > etaMinJEC_[i] && eta <= etaMaxJEC_[i])
        {
            etaBin = i;
            break;
        }
    }

    double lowFact{0.};
    double highFact{0.};

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

    const double a{(highFact - lowFact)
                   / (ptMaxJEC_[ptBin] - ptMinJEC_[ptBin])};
    const double b{(lowFact * ptMaxJEC_[ptBin] - highFact * ptMinJEC_[ptBin])
                   / (ptMaxJEC_[ptBin] - ptMinJEC_[ptBin])};

    return a * pt + b;
}

std::pair<double, double>
    JetCorrectionUncertainty::getMetAfterJESUnc(double metPx,
                                                double metPy,
                                                const MvaEvent& tree,
                                                const int jesUD) const
{
    for (int i{0}; i != tree.numJetPF2PAT; i++)
    {
        metPx += tree.jetPF2PATPx[i];
        metPy += tree.jetPF2PATPy[i];

        const double uncertainty{
            getUncertainty(tree.jetPF2PATPt[i], tree.jetPF2PATEta[i], jesUD)};

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
    return {metPx, metPy};
}
