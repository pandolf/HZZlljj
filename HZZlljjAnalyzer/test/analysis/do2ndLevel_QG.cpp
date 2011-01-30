#include "Ntp1Analyzer_QG.h"
#include <stdlib.h>



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_QG [dataset] [inputfile=""] [flags=""]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_QG* na;

  if( argc<4 ) {
    na = new Ntp1Analyzer_QG(dataset);
  } else {
    std::string flags(argv[3]);
    na = new Ntp1Analyzer_QG(dataset, flags);
  }


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
