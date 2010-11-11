#include "Ntp1Analyzer_TMVA.h"
#include <stdlib.h>



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_TMVA [dataset] [inputfile=""] [flags=""]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_TMVA* na;

  std::string flags = "";
  if( argc<4 ) {
    na = new Ntp1Analyzer_TMVA(dataset, "BESTZ");
  } else {
    std::string flags_str(argv[3]);
    flags = flags_str;
    na = new Ntp1Analyzer_TMVA(dataset, "BESTZ", flags);
  }


//na->AddRequiredTrigger( "HLT_Ele15_LW_L1R" );
//na->AddRequiredTrigger( "HLT_Mu9" );
//na->AddRequiredTrigger( "HLT_Mu15" );


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
