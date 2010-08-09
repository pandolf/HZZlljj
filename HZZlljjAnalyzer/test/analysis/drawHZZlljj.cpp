#include <stdlib.h>
#include "DrawBase.h"
#include "fitTools.h"





int main(int argc, char* argv[]) {

//if( argc != 7 && argc!=8 ) {
//  std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_SIGNAL_dataset] [mc_BG_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')] [flags=\"\"]" << std::endl;
//  exit(23);
//}


  DrawBase* db = new DrawBase("HZZlljj");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "HZZlljjPlots_AllM_vs_Zjets";
  db->set_outputdir(outputdir_str);

  db->set_lumiNormalization( 1000. ); //1 fb-1


  std::string mcSignal130FileName = "HZZlljj_HZZ_qqll_gluonfusion_M130.root";
  TFile* mcSignal130File = TFile::Open(mcSignal130FileName.c_str());
  std::cout << "Opened mc file '" << mcSignal130FileName << "'." << std::endl;
  db->add_mcFile( mcSignal130File, "HZZlljj (130)", kRed+1);

  std::string mcSignal150FileName = "HZZlljj_HZZ_qqll_gluonfusion_M150.root";
  TFile* mcSignal150File = TFile::Open(mcSignal150FileName.c_str());
  std::cout << "Opened mc file '" << mcSignal150FileName << "'." << std::endl;
  db->add_mcFile( mcSignal150File, "HZZlljj (150)", kOrange+1);

  std::string mcSignal500FileName = "HZZlljj_HZZ_qqll_gluonfusion_M500.root";
  TFile* mcSignal500File = TFile::Open("HZZlljj_HZZ_qqll_gluonfusion_M500.root");
  std::cout << "Opened mc file '" << mcSignal500FileName << "'." << std::endl;
  db->add_mcFile( mcSignal500File, "HZZlljj (500)", kRed+3);

  std::string mcZJetsFileName = "HZZlljj_ZJets_madgraph.root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  std::cout << "Opened mc file '" << mcZJetsFileName << "'." << std::endl;
  db->add_mcFile( mcZJetsFile, "Z + Jets", 38, 3001);

  bool log = true;

  db->drawHisto( "ZZInvMass", "", "",1, log);

  db->set_noStack( (bool)true );
  db->set_shapeNormalization();
  
  db->drawHisto("etaJet", "", "", 1);
  db->drawHisto("RchJet", "", "", 1);

  delete db;
  db = 0;

  return 0;

}  


