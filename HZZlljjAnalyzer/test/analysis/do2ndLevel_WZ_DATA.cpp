#include "Ntp1Analyzer_WZ.h"
#include <stdlib.h>



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_WZ [dataset] [inputfile=""] [flags=""]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_WZ* na;

  if( argc<4 ) {
    na = new Ntp1Analyzer_WZ(dataset);
  } else {
    std::string flags(argv[3]);
    na = new Ntp1Analyzer_WZ(dataset, flags);
  }

  na->AddRequiredTrigger( "HLT_Mu9" );
  na->AddRequiredTrigger( "HLT_Mu15_v1" );
  na->AddRequiredTrigger( "HLT_Ele10_LW_L1R" );
  na->AddRequiredTrigger( "HLT_Ele15_LW_L1R" );
  na->AddRequiredTrigger( "HLT_Ele15_SW_L1R" );
  na->AddRequiredTrigger( "HLT_Ele15_SW_CaloEleId_L1R" );
  na->AddRequiredTrigger( "HLT_Ele17_SW_CaloEleId_L1R" );
  na->AddRequiredTrigger( "HLT_Ele17_SW_TightEleId_L1R" );
  na->AddRequiredTrigger( "HLT_Ele17_SW_TighterEleIdIsol_L1R_v2" );
  na->AddRequiredTrigger( "HLT_Ele17_SW_TighterEleIdIsol_L1R_v3" );

  //na->ReadJSONFile("Cert_132440-149442_7TeV_StreamExpress_Collisions10_CMSSWConfig_v3.txt");
  na->ReadJSONFile("Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_CMSSWConfig.txt");

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
