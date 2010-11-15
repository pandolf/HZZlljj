#include <stdlib.h>
#include <iostream>
#include <string>
#include "DrawBase.h"
#include "fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    //std::cout << "USAGE: ./drawHZZlljj [(string) LO/HI/MED] [(string) ZJets dataset=\"ZJets_alpgen\"]" << std::endl;
    std::cout << "USAGE: ./drawHZZlljj [(string) LO/HI/MED] [(string) LEPT_TYPE=\"ALL\"]" << std::endl;
    exit(23);
  }

  std::string lohi(argv[1]);
  if( lohi != "LO" && lohi != "HI" && lohi != "MED" && lohi!="400" && lohi!="500") {
    std::cout << "LO/HI/MED must be set to 'LO' or 'HI' or 'MED' or '400' or '500'. Exiting." << std::endl;
    exit(33);
  }

  std::string zJets_dataset = "ZJets_alpgen";
  //std::string zJets_dataset = "ZJets_madgraph";
//if( argc==3 ) {
//  std::string zJets_dataset_str(argv[2]);
//  zJets_dataset = zJets_dataset_str;
//}

  std::string lept_type = "ALL";
  if( argc==3 ) {
    std::string lept_type_str(argv[2]);
    lept_type = lept_type_str;
  }

  DrawBase* db = new DrawBase("HZZlljj");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "HZZlljjPlots_"+lohi+"mass_vs_"+zJets_dataset;
  //if( lept_type!="ALL" ) outputdir_str += "_" + lept_type;
  outputdir_str += "_" + lept_type;
  db->set_outputdir(outputdir_str);


  std::string flags;
  if( lohi=="LO" ) {
    
    flags = "loMass";

    std::string mcSignal130FileName = "HZZlljj_HZZ_qqll_gluonfusion_M130";
    //if( lept_type!="ALL" ) mcSignal130FileName += "_" + lept_type;
    mcSignal130FileName += "_" + lept_type;
    mcSignal130FileName += ".root";
    TFile* mcSignal130File = TFile::Open(mcSignal130FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal130FileName << "'." << std::endl;
    db->add_mcFile( mcSignal130File, "HZZ_qqll_gluonfusion_M130", "HZZlljj (130)", kRed+1);

    std::string mcSignal150FileName = "HZZlljj_HZZ_qqll_gluonfusion_M150";
    //if( lept_type!="ALL" ) mcSignal150FileName += "_" + lept_type;
    mcSignal150FileName += "_" + lept_type;
    mcSignal150FileName += ".root";
    TFile* mcSignal150File = TFile::Open(mcSignal150FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal150FileName << "'." << std::endl;
    db->add_mcFile( mcSignal150File, "HZZ_qqll_gluonfusion_M150", "HZZlljj (150)", kOrange+1);

  } else if( lohi=="MED" ) {

    flags = "medMass";

    std::string mcSignal200FileName = "HZZlljj_HZZ_qqll_gluonfusion_M200";
    //if( lept_type!="ALL" ) mcSignal200FileName += "_" + lept_type;
    mcSignal200FileName += "_" + lept_type;
    mcSignal200FileName += ".root";
    TFile* mcSignal200File = TFile::Open("HZZlljj_HZZ_qqll_gluonfusion_M200.root");
    std::cout << "Opened mc file '" << mcSignal200FileName << "'." << std::endl;
    db->add_mcFile( mcSignal200File, "HZZ_qqll_gluonfusion_M200", "HZZlljj (200) x40", kRed+1);

  } else { //HI or 400 or 500

    flags = "hiMass";
    if( lohi=="400" || lohi=="500" ) flags=lohi;

    //std::string mcSignal300FileName = "HZZlljj_HZZ_qqll_gluonfusion_M300";
    std::string mcSignal300FileName = "HZZlljj_JHUgen_HiggsSM300_2l2j_FASTSIM";
    mcSignal300FileName += "_" + lept_type;
    mcSignal300FileName += ".root";
  //TFile* mcSignal300File = TFile::Open(mcSignal300FileName.c_str());
  //std::cout << "Opened mc file '" << mcSignal300FileName << "'." << std::endl;
  //db->add_mcFile( mcSignal300File, "HZZ_qqll_gluonfusion_M300", "HZZlljj (300)", kRed+3);

    //std::string mcSignal400FileName = "HZZlljj_HZZ_qqll_gluonfusion_M400";
    std::string mcSignal400FileName = "HZZlljj_JHUgen_HiggsSM400_2l2j_FASTSIM";
    mcSignal400FileName += "_" + lept_type;
    mcSignal400FileName += ".root";
    TFile* mcSignal400File = TFile::Open(mcSignal400FileName.c_str());
    if( lohi=="400" || lohi=="HI" ) {
      std::cout << "Opened mc file '" << mcSignal400FileName << "'." << std::endl;
      db->add_mcFile( mcSignal400File, "HZZ_qqll_gluonfusion_M400", "HZZlljj (400)", kOrange);
    }

    //std::string mcSignal500FileName = "HZZlljj_HZZ_qqll_gluonfusion_M500";
    std::string mcSignal500FileName = "HZZlljj_JHUgen_HiggsSM500_2l2j_FASTSIM";
    mcSignal500FileName += "_" + lept_type;
    mcSignal500FileName += ".root";
    TFile* mcSignal500File = TFile::Open(mcSignal500FileName.c_str());
    if( lohi=="500" || lohi=="HI" ) {
      std::cout << "Opened mc file '" << mcSignal500FileName << "'." << std::endl;
      db->add_mcFile( mcSignal500File, "HZZ_qqll_gluonfusion_M500", "HZZlljj (500)", kRed+1);
    }

  }

  db->set_noStack( (bool)true );
  db->set_shapeNormalization();
  
/*
  // do Z->JetJet comparisons on signal only:
  std::vector< HistoAndName > massZjj_Rch;
  HistoAndName hn1;
  hn1.histoName = "massZjj_RchHIHI";
  hn1.legendName = "Both R_{ch} > 70%";
  massZjj_Rch.push_back( hn1 );
  HistoAndName hn2;
  hn2.histoName = "massZjj_RchLOLO";
  hn2.legendName = "Both R_{ch} < 70%";
  massZjj_Rch.push_back( hn2 );
  db->compareDifferentHistos( massZjj_Rch, "DiJet Invariant Mass [GeV/c^{2}]", "massZjj_vs_Rch");

  // do Z->JetJet comparisons on signal only:
  std::vector< HistoAndName > massZZ_kinfit;
  HistoAndName hn1_ZZ;
  hn1_ZZ.histoName = "ZZInvMass_hiMass_MCassoc";
  hn1_ZZ.legendName = "Preselection";
  massZZ_kinfit.push_back( hn1_ZZ );
  HistoAndName hn2_ZZ;
  hn2_ZZ.histoName = "ZZInvMass_hiMass_MCassoc_ZjjMassConstr";
  hn2_ZZ.legendName = "M_{Z} Constraint";
  massZZ_kinfit.push_back( hn2_ZZ );
  HistoAndName hn3_ZZ;
  hn3_ZZ.histoName = "ZZInvMass_hiMass_MCassoc_kinfit_jets";
  hn3_ZZ.legendName = "Kin. Fit (PFJets)";
  massZZ_kinfit.push_back( hn3_ZZ );
  HistoAndName hn4_ZZ;
  hn4_ZZ.histoName = "ZZInvMass_hiMass_MCassoc_kinfit_cands";
  hn4_ZZ.legendName = "Kin. Fit (PFCands)";
  massZZ_kinfit.push_back( hn4_ZZ );
  db->compareDifferentHistos( massZZ_kinfit, "ZZ Invariant Mass [GeV/c^{2}]", "ZZInvMass_kinfit");
*/

  // then add bg:

  std::string mcZJetsFileName = "HZZlljj_" + zJets_dataset;
  //if( lept_type!="ALL" ) mcZJetsFileName += "_" + lept_type;
  mcZJetsFileName += "_" + lept_type;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  std::cout << "Opened mc file '" << mcZJetsFileName << "'." << std::endl;
  db->add_mcFile( mcZJetsFile, zJets_dataset, "Z + Jets", 38, 3001);

  std::string mcZZFileName = "HZZlljj_ZZ_Spring10";
  //if( lept_type!="ALL" ) mcZZFileName += "_" + lept_type;
  mcZZFileName += "_" + lept_type;
  mcZZFileName += ".root";
  TFile* mcZZFile = TFile::Open(mcZZFileName.c_str());
  std::cout << "Opened mc file '" << mcZZFileName << "'." << std::endl;
  db->add_mcFile( mcZZFile, "ZZ_Spring10", "ZZ", kCyan+1, 3003);

  std::string mcTTbarFileName = "HZZlljj_TTbar_2l_Spring10";
  //if( lept_type!="ALL" ) mcTTbarFileName += "_" + lept_type;
  mcTTbarFileName += "_" + lept_type;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  std::cout << "Opened mc file '" << mcTTbarFileName << "'." << std::endl;
  db->add_mcFile( mcTTbarFile, "TTbar_2l_Spring10", "t#bar{t}", 30, 3002);



  bool log = true;

/*
  db->drawHisto("ptJetLead", "", "", "p_{T} of Leading Jet in the Event [GeV/c]", 1);
  db->drawHisto("ptJetLead2", "", "", "p_{T} of Second Leading Jet in the Event [GeV/c]", 1);
  db->drawHisto("ptJetLead3", "", "", "p_{T} of Third Leading Jet in the Event [GeV/c]", 1);
  db->drawHisto("ptJetRecoil", "", "", "Recoil Jet p_{T} [GeV/c]", 1);
  db->drawHisto("ptJet1", "", "", "First Jet p_{T} [GeV/c]", 1);
  db->drawHisto("ptJet2", "", "", "Second Jet p_{T} [GeV/c]", 1);
  db->drawHisto("ptJet2Rel", "", "", "Relative Second Jet p_{T}", 1);
  db->drawHisto("ptJet2OverLead", "", "", "Second Jet p_{T} / Lead Jet p_{T}", 1);
  db->drawHisto("iJet1", "", "", "Ranking of First Jet", 1);
  db->drawHisto("iJet2", "", "", "Ranking of Second Jet", 1);
  db->drawHisto("iJet2MinusiJet1", "", "", "Ranking of Second Jet Minus First Jet", 1);
  db->drawHisto("iJet2PlusiJet1", "", "", "Ranking of Second Jet Plus First Jet", 1);
  db->drawHisto("ptJet2OverJet1", "", "", "Second Jet p_{T} / First Jet p_{T}", 1);
  db->drawHisto("ptLept1", "", "", "Leading Lepton p_{T} [GeV/c]", 1);
  db->drawHisto("ptLept2", "", "", "Subleading Lepton p_{T} [GeV/c]", 1);
  db->drawHisto("ptLept2OverLept1", "", "", "Subleading Lepton p_{T} / Leading Lepton p_{T}", 1);
  db->drawHisto("etaJet1", "", "", "Leading Jet Pseudorapidity", 1);
  db->drawHisto("etaJet2", "", "", "Subleading Jet Pseudorapidity", 1);
  db->drawHisto("RchJet1", "", "", "Leading Jet R_{ch}", 1);
  db->drawHisto("RchJet2", "", "", "Subleading Jet R_{ch}", 1);
  db->drawHisto("massJet1", "", "", "Leading Jet Mass [GeV/c^{2}]", 1, log);
  db->drawHisto("massJet2", "", "", "Subleading Jet Mass [GeV/c^{2}]", 1, log);
  db->drawHisto("ptZjj", "", "", "Dijet Transverse Momentum [GeV/c]", 1, log);
  db->drawHisto("pzZjj", "", "", "Dijet Longitudinal Momentum [GeV/c]", 1, log);
  db->drawHisto("ptZll", "", "", "Dilepton Transverse Momentum [GeV/c]", 1, log);
  db->drawHisto("pzZll", "", "", "Dilepton Longitudinal Momentum [GeV/c]", 1, log);
  db->drawHisto("ptHardestZ", "", "", "Hardest Z p_{T} [GeV/c]", 1, log);
  db->drawHisto("massZjj", "", "", "Dijet Invariant Mass [GeV/c^{2}]", 1, log);
  db->drawHisto("massZll", "", "", "Dilepton Invariant Mass [GeV/c^{2}]", 1, log);
  db->drawHisto("deltaRll", "", "", "Lepton-Lepton #Delta R", 1, log);
  db->drawHisto("deltaRjj", "", "", "Jet-Jet #Delta R", 1, log);
  db->drawHisto("ptHiggs", "", "", "ZZ Transverse Momentum [GeV/c]", 1, log);
  db->drawHisto("pzHiggs", "", "", "ZZ Longitudinal Momentum [GeV/c]", 1, log);
  db->drawHisto("etaHiggs", "", "", "ZZ Pseudorapidity", 1);
  db->drawHisto("deltaRZZ", "", "", "Z-Z #DeltaR", 1, log);
  db->drawHisto("deltaEtaZZ", "", "", "Z-Z #Delta#eta", 1);
  db->drawHisto("deltaEtaAbsZZ", "", "", "Z-Z #Delta|#eta|", 1);
  db->drawHisto("deltaPhiZZ", "", "", "Z-Z #Delta#phi [rad]", 1);
  db->drawHisto("deltaPtZZ", "", "", "Z-Z #Delta p_{T} [GeV/c]", 1);
  db->drawHisto("pfMet", "", "", "PF Missing E_{T} [GeV]");
  db->drawHisto("pfMet_minusHiggs", "", "", "PF Missing E_{T} - Higgs p_{T} [GeV]");
  db->drawHisto("deltaPhi_ZllRecoil", "", "", "#Delta#phi Between Leptonic Z and Recoil [rad]", 1);
  db->drawHisto("deltaPhi_ZjjRecoil", "", "", "#Delta#phi Between Hadronic Z and Recoil [rad]", 1);
  db->drawHisto("deltaR_ZjjRecoil", "", "", "#DeltaR Between Hadronic Z and Recoil", 1);
  db->drawHisto("deltaPhi_HiggsRecoil", "", "", "#Delta#phi Between Higgs and Recoil [rad]", 1);
  db->drawHisto("ptFullSystem", "", "", "p_{T} of Higgs + Recoil [GeV/c]", 1);
  db->drawHisto("ptHiggsOverRecoil", "", "", "Higgs p_{T} / Recoil p_{T}", 1);
  db->drawHisto("ptRecoilOverJet2", "", "", "Recoil p_{T} / Second Jet p_{T}", 1);

*/
  db->set_noStack( (bool)false );
  db->set_lumiNormalization( 1000. ); //1 fb-1
  //db->set_lumiNormalization( 7. ); //1 fb-1



  db->set_rebin(2);

  db->drawHisto("ZZInvMass_hiMass_loose", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);
  db->drawHisto("ZZInvMass_hiMass_tight", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);
  db->drawHisto("ZZInvMass_hiMass_opt400_LowEff", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);
  db->drawHisto("ZZInvMass_hiMass_opt400_HighEff", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);
  db->drawHisto("ZZInvMass_hiMass_opt500", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);

/*
  db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  flags =  flags + "_ZjjTag";
  db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  if( lohi=="LO" ) {
    flags =  flags + "_ZllAntiTag";
    db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    flags =  flags + "_Rch40";
    db->drawHisto( "ZZInvMass", "", flags, "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  } else if( lohi=="MED" ) {
    db->drawHisto( "ZZInvMass", "", "medMass_fullSelection_nokin", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    db->drawHisto( "ZZInvMass", "", "medMass_fullSelection", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  } else if( lohi=="HI" ) {
    db->drawHisto( "ZZInvMass", "", "hiMass_fullSelection_tightOLD", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    db->drawHisto( "ZZInvMass", "", "hiMass_fullSelection_tight", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    db->drawHisto( "ZZInvMass", "", "hiMass_fullSelection_medium", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    db->drawHisto( "ZZInvMass", "", "hiMass_fullSelection_medium_ZjjMassConstr", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
    db->drawHisto( "ZZInvMass", "", "hiMass_fullSelection_loose", "ZZ Invariant Mass [GeV/c^{2}]", 1, log);
  }
*/


  delete db;
  db = 0;

  return 0;

}  


