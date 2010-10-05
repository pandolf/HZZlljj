#include "Ntp1Analyzer_HZZlljj.h"
#include <stdlib.h>



int main( int argc, char* argv[]) {

  if( argc < 2 ) {
    std::cout << "USAGE: ./do2ndLevel_Hzzlljj [dataset]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_HZZlljj* na = new Ntp1Analyzer_HZZlljj(dataset);

  na->AddRequiredTrigger( "HLT_Ele15_LW_L1R" );
  na->AddRequiredTrigger( "HLT_Mu9" );
  na->AddRequiredTrigger( "HLT_Mu15" );
  na->LoadInput();
  na->Loop();

  delete na;

  return 0;

}
