#ifndef _jetCorrectionUncertainty_hpp_
#define _jetCorrectionUncertainty_hpp_

#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>

class MvaEvent;

class JetCorrectionUncertainty{

        public:
        // Constructor
        JetCorrectionUncertainty( std::string dataFile );
        ~JetCorrectionUncertainty();

        double getUncertainty( double pt, double eta, int jesUD );
        std::pair< double, double> getMetAfterJESUnc( double metPx, double metPy, MvaEvent tree, int jesUD );

        private:
        
        std::vector<double> ptMinJEC_;
        std::vector<double> ptMaxJEC_;

        std::vector<double> etaMinJEC_;
        std::vector<double> etaMaxJEC_;

        std::vector<std::vector <float> > jecSFUp_;
        std::vector<std::vector <float> > jecSFDown_;
};

#endif

