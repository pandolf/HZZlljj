#include "Ntp1Analyzer_HWWlvjj.h"
#include <stdlib.h>



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_HWWlvjj [dataset] [inputfile=""] [flags=""]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_HWWlvjj* na;

  if( argc<4 ) {
    na = new Ntp1Analyzer_HWWlvjj(dataset);
  } else {
    std::string flags(argv[3]);
    na = new Ntp1Analyzer_HWWlvjj(dataset, flags);
  }

//na->AddRequiredTrigger( "HLT_Ele15_LW_L1R" );
//na->AddRequiredTrigger( "HLT_Mu9" );
//na->AddRequiredTrigger( "HLT_Mu15" );
na->ReadJSONFile("Cert_160404-167913_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");

  if( argc==2 ) {
    na->LoadInput();
  } else {
    std::string inputfile(argv[2]);
    na->LoadInputFromFile(inputfile.c_str());
  }

  na->Loop();

  delete na;

  return 0;

}
