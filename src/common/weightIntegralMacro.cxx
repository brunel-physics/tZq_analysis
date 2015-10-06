
#include <stdlib.h> 
#include <iostream>
#include <sstream>

#include "SkimFileEvent.hpp"

int main(int argc, char *argv[])
{

  std::cout << "NOTE: This works only with skim files made with the eventWeight in tree branch!!!" << std::endl;

  Int_t lNumberOfFiles = argc;
  //  std::cout << "Number of Files: " << argc-1 << std::endl;

  for ( Int_t i = 1; i < lNumberOfFiles; i++ )
  {
    std::stringstream lSStr;
    //    std::cout << lSStr.str().c_str() << std::endl;
    lSStr << argv[i];

    TFile *lFile = new TFile ( (lSStr.str()).c_str() );
    TTree *lTree = (TTree*)lFile->Get("tree");  

    SkimFileEvent* lEvent = new SkimFileEvent(lTree);
    Float_t lEventWeightIntegral = 0.0;	
    Int_t lNumEvents = lTree->GetEntries();

    std::cout << "Number of Events for "+lSStr.str()+": " << lNumEvents << std::endl;

    for ( Int_t j = 0; j < lNumEvents; j++ )
      {
	lTree->GetEvent(j);
	Float_t lTempWeight = lEvent->eventWeight;
	//	std::cout << "lTempWeight: " << lTempWeight << std::endl;
	lEventWeightIntegral += lTempWeight;
	//	std::cout << "lEventWeightIntegral: " << lEventWeightIntegral << std::endl;
      }

    std::cout << "EventWeightIntegral for "+lSStr.str()+": " << lEventWeightIntegral << std::endl;
  }
}
