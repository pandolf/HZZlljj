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


  TString dataset_tstr(dataset);
//if( dataset_tstr.BeginsWith("DoubleMu") ) {
//  na->AddRequiredTrigger( "HLT_Mu13_Mu8" );
//  na->AddRequiredTrigger( "HLT_DoubleMu7" );
//  //na->AddRequiredTriggerNOT( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL" );
//} else if( dataset_tstr.BeginsWith("SingleMu") ) {
//  na->AddRequiredTrigger( "HLT_IsoMu24" );
//  na->AddRequiredTriggerNOT( "HLT_Mu13_Mu8" );
//  na->AddRequiredTriggerNOT( "HLT_DoubleMu7" );
//  na->AddRequiredTriggerNOT( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL" );
//} else if( dataset_tstr.BeginsWith("DoubleElectron") ) {
//  //na->AddRequiredTrigger( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL" );
//  ////na->AddRequiredTriggerNOT( "HLT_Mu13_Mu8" );
//  ////na->AddRequiredTriggerNOT( "HLT_DoubleMu7" );
//}


  //na->ReadJSONFile("Cert_160404-167913_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
  //na->ReadJSONFile("Cert_goodCollisions2011_CMSSWConfig.txt");
//if( dataset_tstr.BeginsWith("DoubleElectron") )
//  na->ReadJSONFile("Cert_160404-172802_7TeV_PromptReco_Collisions11_DoubleElectron_plusMay10addon_CMSSWConfig.txt");
//else
//  na->ReadJSONFile("Cert_160404-172802_7TeV_PromptReco_Collisions11_plusMay10addon_CMSSWConfig.txt");
////na->ReadJSONFile("Cert_160404-171116_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");

  //na->ReadJSONFile("Cert_160404-177053_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
  na->ReadJSONFile("Cert_160404-179431_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");

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
