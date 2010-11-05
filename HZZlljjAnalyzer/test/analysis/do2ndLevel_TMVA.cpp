#include "Ntp1Analyzer_TMVA.h"
#include <stdlib.h>



int main( int argc, char* argv[]) {

  if( argc < 2 ) {
    std::cout << "USAGE: ./do2ndLevel_TMVA [dataset]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_TMVA* na = new Ntp1Analyzer_TMVA(dataset);

  na->LoadInput();
  na->Loop();

  delete na;

  return 0;

}
