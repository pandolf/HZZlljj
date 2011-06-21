#include "Ntp1Analyzer_HZZlljj.h"
#include <stdlib.h>



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_HZZlljj [dataset] [inputfile=""] [flags=""]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_HZZlljj* na;

  if( argc<4 ) {
    na = new Ntp1Analyzer_HZZlljj(dataset);
  } else {
    std::string flags(argv[3]);
    na = new Ntp1Analyzer_HZZlljj(dataset, flags);
  }

  //na->AddRequiredTrigger( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2" );
  //na->AddRequiredTrigger( "HLT_DoubleMu7_v1" );

//na->AddRequiredTrigger( "HLT_Mu9" );
//na->AddRequiredTrigger( "HLT_Mu15_v1" );
//na->AddRequiredTrigger( "HLT_Ele10_LW_L1R" );
//na->AddRequiredTrigger( "HLT_Ele15_LW_L1R" );
//na->AddRequiredTrigger( "HLT_Ele15_SW_L1R" );
//na->AddRequiredTrigger( "HLT_Ele15_SW_CaloEleId_L1R" );
//na->AddRequiredTrigger( "HLT_Ele17_SW_CaloEleId_L1R" );
//na->AddRequiredTrigger( "HLT_Ele17_SW_TightEleId_L1R" );
//na->AddRequiredTrigger( "HLT_Ele17_SW_TighterEleIdIsol_L1R_v2" );
//na->AddRequiredTrigger( "HLT_Ele17_SW_TighterEleIdIsol_L1R_v3" );

  //na->ReadJSONFile("Cert_160404-163369_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
  //na->ReadJSONFile("Cert_160404-163757_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
  //na->ReadJSONFile("Cert_160404-165121_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
  na->ReadJSONFile("Cert_160404-165542_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");

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
