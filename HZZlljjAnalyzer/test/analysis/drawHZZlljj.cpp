#include <stdlib.h>
#include <iostream>
#include <string>
#include "DrawBase.h"
#include "fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    //std::cout << "USAGE: ./drawHZZlljj [(string) LO/HI/MED] [(string) ZJets dataset=\"ZJets_alpgen\"]" << std::endl;
    std::cout << "USAGE: ./drawHZZlljj [(string) selectionType] [(string) LEPT_TYPE=\"ALL\"]" << std::endl;
    exit(23);
  }

  std::string selType(argv[1]);
//if( selType != "LO" && selType != "HI" && selType != "MED" && selType!="400" && selType!="500") {
//  std::cout << "LO/HI/MED must be set to 'LO' or 'HI' or 'MED' or '400' or '500'. Exiting." << std::endl;
//  exit(33);
//}

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

  std::string outputdir_str = "HZZlljjPlots_"+selType;
  //if( lept_type!="ALL" ) outputdir_str += "_" + lept_type;
  outputdir_str += "_" + lept_type;
  db->set_outputdir(outputdir_str);

  std::string lohi;
  if( selType=="opt400" || selType=="opt500" || selType=="opt600" || selType=="tight" ) lohi="HI"; //temp solution
  else lohi="LO";

  std::string flags;
  if( lohi=="LO" ) {
    
    flags = "loMass";

    std::string mcSignal130FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-130_7TeV-jhu-pythia6";
    mcSignal130FileName += "_" + selType;
    mcSignal130FileName += "_" + lept_type;
    mcSignal130FileName += ".root";
    TFile* mcSignal130File = TFile::Open(mcSignal130FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal130FileName << "'." << std::endl;
    db->add_mcFile( mcSignal130File, "HZZ_qqll_gluonfusion_M130", "HZZlljj (130)", kRed+1);

    std::string mcSignal150FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-150_7TeV-jhu-pythia6";
    mcSignal150FileName += "_" + selType;
    mcSignal150FileName += "_" + lept_type;
    mcSignal150FileName += ".root";
    TFile* mcSignal150File = TFile::Open(mcSignal150FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal150FileName << "'." << std::endl;
    db->add_mcFile( mcSignal150File, "HZZ_qqll_gluonfusion_M150", "HZZlljj (150)", kOrange+1);

  } else if( lohi=="MED" ) {

    flags = "medMass";

    std::string mcSignal200FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-200_7TeV-jhu-pythia6";
    //if( lept_type!="ALL" ) mcSignal200FileName += "_" + lept_type;
    mcSignal200FileName += "_" + selType;
    mcSignal200FileName += "_" + lept_type;
    mcSignal200FileName += ".root";
    TFile* mcSignal200File = TFile::Open("HZZlljj_HZZ_qqll_gluonfusion_M200.root");
    std::cout << "Opened mc file '" << mcSignal200FileName << "'." << std::endl;
    db->add_mcFile( mcSignal200File, "HZZ_qqll_gluonfusion_M200", "HZZlljj (200) x40", kRed+1);

  } else { //HI or 400 or 500

    flags = "hiMass";
    if( lohi=="400" || lohi=="500" ) flags=lohi;

    std::string mcSignal300FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-300_7TeV-jhu-pythia6";
    mcSignal300FileName += "_" + selType;
    mcSignal300FileName += "_" + lept_type;
    mcSignal300FileName += ".root";
  //TFile* mcSignal300File = TFile::Open(mcSignal300FileName.c_str());
  //std::cout << "Opened mc file '" << mcSignal300FileName << "'." << std::endl;
  //db->add_mcFile( mcSignal300File, "HZZ_qqll_gluonfusion_M300", "HZZlljj (300)", kRed+3);

    std::string mcSignal400FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6";
    mcSignal400FileName += "_" + selType;
    mcSignal400FileName += "_" + lept_type;
    mcSignal400FileName += ".root";
    TFile* mcSignal400File = TFile::Open(mcSignal400FileName.c_str());
    if( lohi=="400" || lohi=="HI" ) {
      std::cout << "Opened mc file '" << mcSignal400FileName << "'." << std::endl;
      db->add_mcFile( mcSignal400File, "HZZ_qqll_gluonfusion_M400", "HZZlljj (400)", kOrange);
    }

    std::string mcSignal500FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-500_7TeV-jhu-pythia6";
    mcSignal500FileName += "_" + selType;
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

  std::string mcZJetsFileName = "HZZlljj_ZJets_alpgen_TuneZ2_Fall10";
  //if( lept_type!="ALL" ) mcZJetsFileName += "_" + lept_type;
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + lept_type;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  std::cout << "Opened mc file '" << mcZJetsFileName << "'." << std::endl;
  db->add_mcFile( mcZJetsFile, zJets_dataset, "Z + Jets", 38, 3001);

//std::string mcZZFileName = "HZZlljj_ZZ_Spring10";
////if( lept_type!="ALL" ) mcZZFileName += "_" + lept_type;
//mcZZFileName += "_" + lept_type;
//mcZZFileName += ".root";
//TFile* mcZZFile = TFile::Open(mcZZFileName.c_str());
//std::cout << "Opened mc file '" << mcZZFileName << "'." << std::endl;
//db->add_mcFile( mcZZFile, "ZZ_Spring10", "ZZ", kCyan+1, 3003);

//std::string mcTTbarFileName = "HZZlljj_TTbar_2l_Spring10";
////if( lept_type!="ALL" ) mcTTbarFileName += "_" + lept_type;
//mcTTbarFileName += "_" + lept_type;
//mcTTbarFileName += ".root";
//TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
//std::cout << "Opened mc file '" << mcTTbarFileName << "'." << std::endl;
//db->add_mcFile( mcTTbarFile, "TTbar_2l_Spring10", "t#bar{t}", 30, 3002);



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
  //db->set_noStack( (bool)false );
  db->set_shapeNormalization();
  //db->set_lumiNormalization( 1000. ); //1 fb-1
  //db->set_lumiNormalization( 7. ); //1 fb-1



  //db->set_rebin(2);

  db->drawHisto("ptJet_all_presel", "Jet Transverse Momentum", "GeV/c", "Jets", log);
  db->drawHisto("ptDJet_all_presel", "ptD", "", "Jets", log);
  db->drawHisto("rmsCandJet_all_presel", "PFCandidate RMS", "", "Jets", log);
  db->drawHisto("nChargedJet_all_presel", "Charged Multiplicity", "", "Jets", log);
  db->drawHisto("nNeutralJet_all_presel", "Neutral Multiplicity", "", "Jets", log);

//db->drawHisto("ptJetBest1", "Jet Transverse Momentum", "GeV/c", "Jets", log);
//db->drawHisto("ptDJetBest1", "ptD", "", "Jets", log);
//db->drawHisto("rmsCandJetBest1", "PFCandidate RMS", "", "Jets", log);
//db->drawHisto("nChargedJetBest1", "Charged Multiplicity", "", "Jets", log);
//db->drawHisto("nNeutralJetBest1", "Neutral Multiplicity", "", "Jets", log);

//db->drawHisto("ptJetBest2", "Jet Transverse Momentum", "GeV/c", "Jets", log);
//db->drawHisto("ptDJetBest2", "ptD", "", "Jets", log);
//db->drawHisto("rmsCandJetBest2", "PFCandidate RMS", "", "Jets", log);
//db->drawHisto("nChargedJetBest2", "Charged Multiplicity", "", "Jets", log);
//db->drawHisto("nNeutralJetBest2", "Neutral Multiplicity", "", "Jets", log);

  db->drawHisto("ptJetRecoil", "Recoil Jet Jet Transverse Momentum", "GeV/c", "Jets", log);
  db->drawHisto("ptDJetRecoil", "Recoil Jet ptD", "", "Jets", log);
  db->drawHisto("rmsCandJetRecoil", "Recoil Jet PFCandidate RMS", "", "Jets", log);
  db->drawHisto("nChargedJetRecoil", "Recoil Jet Charged Multiplicity", "", "Jets", log);
  db->drawHisto("nNeutralJetRecoil", "Recoil Jet Neutral Multiplicity", "", "Jets", log);

  db->drawHisto("ptDJet1_partMatched", "First Jet ptD", "", "Jets", log);
  db->drawHisto("rmsCandJet1_partMatched", "First Jet PFCandidate RMS", "", "Jets", log);
  db->drawHisto("nChargedJet1_partMatched", "First Jet Charged Multiplicity", "", "Jets", log);
  db->drawHisto("nNeutralJet1_partMatched", "First Jet Neutral Multiplicity", "", "Jets", log);

  db->drawHisto("ptDJet2_partMatched", "Second Jet ptD", "", "Jets", log);
  db->drawHisto("rmsCandJet2_partMatched", "Second Jet PFCandidate RMS", "", "Jets", log);
  db->drawHisto("nChargedJet2_partMatched", "Second Jet Charged Multiplicity", "", "Jets", log);
  db->drawHisto("nNeutralJet2_partMatched", "Second Jet Neutral Multiplicity", "", "Jets", log);

  db->drawHisto("partFlavorJet1_partMatched", "First Jet Parton Flavour", "", "Jets", log);
  db->drawHisto("partFlavorJet2_partMatched", "Second Jet Parton Flavour", "", "Jets", log);


//db->drawHisto("ZZInvMass_hiMass_loose", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);
//db->drawHisto("ZZInvMass_hiMass_tight", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);
//db->drawHisto("ZZInvMass_hiMass_opt400_LowEff", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);
//db->drawHisto("ZZInvMass_hiMass_opt400_HighEff", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);
//db->drawHisto("ZZInvMass_hiMass_opt500", "", "", "ZZ Invariant Mass [GeV/c^{2}]", 1);

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


