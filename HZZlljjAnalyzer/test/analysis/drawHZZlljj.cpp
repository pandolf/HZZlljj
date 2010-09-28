#include <stdlib.h>
#include <iostream>
#include <string>
#include "DrawBase.h"
#include "fitTools.h"





int main(int argc, char* argv[]) {

  if( argc != 2 ) {
    std::cout << "USAGE: ./drawHZZlljj [(string) LO/HI/MED]" << std::endl;
    exit(23);
  }

  std::string lohi(argv[1]);
  if( lohi != "LO" && lohi != "HI" && lohi != "MED") {
    std::cout << "LO/HI/MED must be set to 'LO' or 'HI' or 'MED'. Exiting." << std::endl;
    exit(33);
  }

  DrawBase* db = new DrawBase("HZZlljj");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "HZZlljjPlots_"+lohi+"mass_vs_Zjets";
  db->set_outputdir(outputdir_str);


  std::string flags;
  if( lohi=="LO" ) {
    
    flags = "loMass";

    std::string mcSignal130FileName = "HZZlljj_HZZ_qqll_gluonfusion_M130.root";
    TFile* mcSignal130File = TFile::Open(mcSignal130FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal130FileName << "'." << std::endl;
    db->add_mcFile( mcSignal130File, "HZZ_qqll_gluonfusion_M130", "HZZlljj (130)", kRed+1);

    std::string mcSignal150FileName = "HZZlljj_HZZ_qqll_gluonfusion_M150.root";
    TFile* mcSignal150File = TFile::Open(mcSignal150FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal150FileName << "'." << std::endl;
    db->add_mcFile( mcSignal150File, "HZZ_qqll_gluonfusion_M150", "HZZlljj (150)", kOrange+1);

  } else if( lohi=="MED" ) {

    flags = "medMass";

    std::string mcSignal200FileName = "HZZlljj_HZZ_qqll_gluonfusion_M200.root";
    TFile* mcSignal200File = TFile::Open("HZZlljj_HZZ_qqll_gluonfusion_M200.root");
    std::cout << "Opened mc file '" << mcSignal200FileName << "'." << std::endl;
    db->add_mcFile( mcSignal200File, "HZZ_qqll_gluonfusion_M200", "HZZlljj (200)", kRed+1);

  } else { //HI

    flags = "hiMass";

    std::string mcSignal300FileName = "HZZlljj_HZZ_qqll_gluonfusion_M300.root";
    TFile* mcSignal300File = TFile::Open("HZZlljj_HZZ_qqll_gluonfusion_M300.root");
    std::cout << "Opened mc file '" << mcSignal300FileName << "'." << std::endl;
    db->add_mcFile( mcSignal300File, "HZZ_qqll_gluonfusion_M300", "HZZlljj (300)", kRed+3);

    std::string mcSignal400FileName = "HZZlljj_HZZ_qqll_gluonfusion_M400.root";
    TFile* mcSignal400File = TFile::Open("HZZlljj_HZZ_qqll_gluonfusion_M400.root");
    std::cout << "Opened mc file '" << mcSignal400FileName << "'." << std::endl;
    db->add_mcFile( mcSignal400File, "HZZ_qqll_gluonfusion_M400", "HZZlljj (400)", kOrange);

    std::string mcSignal500FileName = "HZZlljj_HZZ_qqll_gluonfusion_M500.root";
    TFile* mcSignal500File = TFile::Open("HZZlljj_HZZ_qqll_gluonfusion_M500.root");
    std::cout << "Opened mc file '" << mcSignal500FileName << "'." << std::endl;
    db->add_mcFile( mcSignal500File, "HZZ_qqll_gluonfusion_M500", "HZZlljj (500)", kRed+1);

  }

  std::string mcZJetsFileName = "HZZlljj_ZJets_madgraph.root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  std::cout << "Opened mc file '" << mcZJetsFileName << "'." << std::endl;
  db->add_mcFile( mcZJetsFile, "ZJets_madgraph", "Z + Jets", 38, 3001);


  bool log = true;

  db->set_noStack( (bool)true );
  db->set_shapeNormalization();
  
  db->drawHisto("ptJet1", "", "", "Leading Jet p_{T} [GeV/c]", 1);
  db->drawHisto("ptJet2", "", "", "Subleading Jet p_{T} [GeV/c]", 1);
  db->drawHisto("ptLept1", "", "", "Leading Lepton p_{T} [GeV/c]", 1);
  db->drawHisto("ptLept2", "", "", "Subleading Lepton p_{T} [GeV/c]", 1);
  db->drawHisto("etaJet1", "", "", "Leading Jet Pseudorapidity", 1);
  db->drawHisto("etaJet2", "", "", "Subleading Jet Pseudorapidity", 1);
  db->drawHisto("RchJet1", "", "", "Leading Jet R_{ch}", 1);
  db->drawHisto("RchJet2", "", "", "Subleading Jet R_{ch}", 1);
  db->drawHisto("JetJetPt", "", "", "Dijet Transverse Momentum [GeV/c]", 1, log);
  db->drawHisto("LeptLeptPt", "", "", "Dilepton Transverse Momentum [GeV/c]", 1, log);
  db->drawHisto("ptHardestZ", "", "", "Hardest Z p_{T} [GeV/c]", 1, log);
  db->drawHisto("JetJetInvMass", "", "", "Dijet Invariant Mass [GeV/c^{2}]", 1, log);
  db->drawHisto("LeptLeptInvMass", "", "", "Dilepton Invariant Mass [GeV/c^2]", 1, log);
  db->drawHisto("deltaRjj", "", "", "Jet-Jet #Delta R", 1, log);
  db->drawHisto("ptZZ", "", "", "ZZ Transverse Momentum [GeV/c]", 1, log);
  db->drawHisto("deltaRZZ", "", "", "Z-Z #Delta R", 1, log);


  db->set_noStack( (bool)false );
  db->set_lumiNormalization( 1000. ); //1 fb-1


  db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  flags =  flags + "_ZjjTag";
  db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  if( lohi=="LO" ) {
    flags =  flags + "_ZllAntiTag";
    db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    flags =  flags + "_Rch40";
    db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  } else if( lohi=="MED" ) {
    flags = flags + "_kinem";
    db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  } else if( lohi=="HI" ) {
    db->drawHisto( "ZZInvMass", "", "hiMass_fullSelection_tight", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    db->drawHisto( "ZZInvMass", "", "hiMass_fullSelection_medium", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    db->drawHisto( "ZZInvMass", "", "hiMass_fullSelection_loose", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  }

  delete db;
  db = 0;

  return 0;

}  


