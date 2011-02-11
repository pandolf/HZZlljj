#include <stdlib.h>
#include <iostream>
#include <string>
#include "DrawBase.h"
#include "fitTools.h"



void draw_vs_pt_plots( DrawBase* db, int nPtBins, Double_t* ptBins, const std::string& histoName, const std::string& axisName, const std::string& units="", const std::string& instanceName="Entries", bool log=false );


int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    //std::cout << "USAGE: ./drawHZZlljj [(string) LO/HI/MED] [(string) ZJets dataset=\"ZJets_alpgen\"]" << std::endl;
    std::cout << "USAGE: ./drawHZZlljj [(string) selectionType] [(string) leptType=\"ALL\"]" << std::endl;
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

  std::string leptType = "ALL";
  if( argc==3 ) {
    std::string leptType_str(argv[2]);
    leptType = leptType_str;
  }

  DrawBase* db = new DrawBase("HZZlljj");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "HZZlljjPlots_"+selType;
  //if( leptType!="ALL" ) outputdir_str += "_" + leptType;
  outputdir_str += "_" + leptType;
  db->set_outputdir(outputdir_str);

  std::string lohi;
  if( selType=="opt400" || selType=="opt500" || selType=="opt600" || selType=="tight" ) lohi="HI"; //temp solution
  else if( selType=="loose" ) lohi="MED";
  else lohi="LO";

  std::string flags;
  if( lohi=="LO" ) {
    
    flags = "loMass";

    std::string mcSignal130FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-130_7TeV-jhu-pythia6";
    mcSignal130FileName += "_" + selType;
    mcSignal130FileName += "_" + leptType;
    mcSignal130FileName += ".root";
    TFile* mcSignal130File = TFile::Open(mcSignal130FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal130FileName << "'." << std::endl;
    db->add_mcFile( mcSignal130File, "HZZ_qqll_gluonfusion_M130", "HZZlljj (130)", kRed+1);

    std::string mcSignal150FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-150_7TeV-jhu-pythia6";
    mcSignal150FileName += "_" + selType;
    mcSignal150FileName += "_" + leptType;
    mcSignal150FileName += ".root";
    TFile* mcSignal150File = TFile::Open(mcSignal150FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal150FileName << "'." << std::endl;
    db->add_mcFile( mcSignal150File, "HZZ_qqll_gluonfusion_M150", "HZZlljj (150)", kOrange+1);

  } else if( lohi=="MED" ) {

    flags = "medMass";

    std::string mcSignal200FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-200_7TeV-jhu-pythia6";
    //if( leptType!="ALL" ) mcSignal200FileName += "_" + leptType;
    mcSignal200FileName += "_" + selType;
    mcSignal200FileName += "_" + leptType;
    mcSignal200FileName += ".root";
    TFile* mcSignal200File = TFile::Open(mcSignal200FileName.c_str());
    std::cout << "Opened mc file '" << mcSignal200FileName << "'." << std::endl;
    db->add_mcFile( mcSignal200File, "HZZ_qqll_gluonfusion_M200", "HZZlljj (200)", kRed+1);

  } else { //HI or 400 or 500

    flags = "hiMass";
    if( lohi=="400" || lohi=="500" ) flags=lohi;

    std::string mcSignal300FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-300_7TeV-jhu-pythia6";
    mcSignal300FileName += "_" + selType;
    mcSignal300FileName += "_" + leptType;
    mcSignal300FileName += ".root";
  //TFile* mcSignal300File = TFile::Open(mcSignal300FileName.c_str());
  //std::cout << "Opened mc file '" << mcSignal300FileName << "'." << std::endl;
  //db->add_mcFile( mcSignal300File, "HZZ_qqll_gluonfusion_M300", "HZZlljj (300)", kRed+3);

    std::string mcSignal400FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6";
    mcSignal400FileName += "_" + selType;
    mcSignal400FileName += "_" + leptType;
    mcSignal400FileName += ".root";
    TFile* mcSignal400File = TFile::Open(mcSignal400FileName.c_str());
    if( lohi=="400" || lohi=="HI" ) {
      std::cout << "Opened mc file '" << mcSignal400FileName << "'." << std::endl;
      db->add_mcFile( mcSignal400File, "HZZ_qqll_gluonfusion_M400", "HZZlljj (400)", kOrange);
    }

    std::string mcSignal500FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-500_7TeV-jhu-pythia6";
    mcSignal500FileName += "_" + selType;
    mcSignal500FileName += "_" + leptType;
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
*/

  std::vector< HistoAndName > massZZ_QG;
  HistoAndName hn1_ZZ;
  hn1_ZZ.histoName = "ZZInvMass_hiMass";
  hn1_ZZ.legendName = "Kinematic Selection";
  massZZ_QG.push_back( hn1_ZZ );
  HistoAndName hn2_ZZ;
  hn2_ZZ.histoName = "ZZInvMass_hiMass_QGlikeli";
  hn2_ZZ.legendName = "Kin. + Q-G Selection";
  massZZ_QG.push_back( hn2_ZZ );
  db->compareDifferentHistos( massZZ_QG, "ZZInvMass_QG", "ZZ Invariant Mass",  "GeV/c^{2}");

/*
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
  //if( leptType!="ALL" ) mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  std::cout << "Opened mc file '" << mcZJetsFileName << "'." << std::endl;
  db->add_mcFile( mcZJetsFile, zJets_dataset, "Z + Jets", 38, 3001);

//std::string mcZZFileName = "HZZlljj_ZZ_Spring10";
////if( leptType!="ALL" ) mcZZFileName += "_" + leptType;
//mcZZFileName += "_" + leptType;
//mcZZFileName += ".root";
//TFile* mcZZFile = TFile::Open(mcZZFileName.c_str());
//std::cout << "Opened mc file '" << mcZZFileName << "'." << std::endl;
//db->add_mcFile( mcZZFile, "ZZ_Spring10", "ZZ", kCyan+1, 3003);

//std::string mcTTbarFileName = "HZZlljj_TTbar_2l_Spring10";
////if( leptType!="ALL" ) mcTTbarFileName += "_" + leptType;
//mcTTbarFileName += "_" + leptType;
//mcTTbarFileName += ".root";
//TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
//std::cout << "Opened mc file '" << mcTTbarFileName << "'." << std::endl;
//db->add_mcFile( mcTTbarFile, "TTbar_2l_Spring10", "t#bar{t}", 30, 3002);



  bool log = true;

  //db->set_noStack( (bool)false );
  db->set_shapeNormalization();
  //db->set_lumiNormalization( 1000. ); //1 fb-1
  //db->set_lumiNormalization( 7. ); //1 fb-1

  db->drawHisto("nJets_presel", "Jet Multiplicity (p_{T} > 30 GeV/c)", "", "Events", log);
  db->drawHisto("ptJet_all_presel", "Jet Transverse Momentum", "GeV/c", "Jets", log);
  db->drawHisto("mZjj_all_presel", "DiJet Invariant Mass", "GeV/c^{2}", "Jet Pairs", log);

  db->drawHisto("deltaRjj_all_presel", "#DeltaR Between Jets (p_{T} > 30 GeV/c)", "", "Jet Pairs");
  db->drawHisto("deltaRll_presel", "#DeltaR Between Leptons", "", "Lepton Pairs");
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("etaLept1_presel", "Lead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaLept2_presel", "Sublead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaJet_all_presel", "Jet Pseudorapidity", "", "Jets");

  db->set_yAxisMaxScale( 1.4 );
  db->drawHisto("ptLept1_presel", "Lead Lepton Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptLept2_presel", "Sublead Lepton Transverse Momentum", "GeV/c", "Events", log);

  db->drawHisto("mZll_presel", "Dilepton Invariant Mass", "GeV/c^{2}", "Events", log);
  if( leptType=="ALL" || leptType=="MU" )
    db->drawHisto("mZmumu_presel", "DiMuon Invariant Mass", "GeV/c^{2}", "Events", log);
  if( leptType=="ALL" || leptType=="ELE" )
    db->drawHisto("mZee_presel", "DiElectron Invariant Mass", "GeV/c^{2}", "Events", log);

  db->drawHisto("ptZll_presel", "Dilepton Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZjj_all_presel", "Dijet Transverse Momentum", "GeV/c", "Events", log);


  db->drawHisto("ptLept1", "Leading Lepton Transverse Momentum", "GeV/c", "Events");
  db->drawHisto("ptLept2", "Subleading Lepton Transverse Momentum", "GeV/c", "Events");
  db->drawHisto("ptJet1", "Leading Jet Transverse Momentum", "GeV/c", "Events");
  db->drawHisto("ptJet2", "Subleading Jet Transverse Momentum", "GeV/c", "Events");

  db->drawHisto("partFlavorJet1", "Leading Jet Parton Flavor", "", "Events");
  db->drawHisto("partFlavorJet2", "Subleading Jet Parton Flavor", "", "Events");
  db->drawHisto("QGLikelihoodJet1", "Leading Jet Q-G Likelihood", "", "Events");
  db->drawHisto("QGLikelihoodJet2", "Subleading Jet Q-G Likelihood", "", "Events");
  db->drawHisto("QGLikelihoodJet1_norms", "Leading Jet Q-G Likelihood", "", "Events");
  db->drawHisto("QGLikelihoodJet2_norms", "Subleading Jet Q-G Likelihood", "", "Events");
  db->drawHisto("QGLikelihoodSum", "Leading + Subleading Jet Q-G Likelihood", "", "Events");

  db->drawHisto("QGLikelihoodSum", "Leading + Subleading Jet Q-G Likelihood", "", "Events");


  const int nPtBins = 20;
  Double_t ptBins[nPtBins+1];
  fitTools::getBins_int( nPtBins+1, ptBins, 15., 1000. );

//draw_vs_pt_plots( db, nPtBins, ptBins, "QGLikelihoodJet1", "Leading Jet Q-G Likelihood");
//draw_vs_pt_plots( db, nPtBins, ptBins, "QGLikelihoodJet2", "Subleading Jet Q-G Likelihood");


  std::string mcVVFileName = "HZZlljj_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10";
  //if( leptType!="ALL" ) mcZme += "_" + selType;
  mcVVFileName += "_" + selType;
  mcVVFileName += "_" + leptType;
  mcVVFileName += ".root";
  TFile* mcVVFile = TFile::Open(mcVVFileName.c_str());
  std::cout << "Opened mc file '" << mcVVFileName << "'." << std::endl;
  db->add_mcFile( mcVVFile, "VVtoAnything_TuneZ2", "ZZ/WZ/WW", kCyan+1, 3003);

  std::string mcTTbarFileName = "HZZlljj_TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10";
  //if( leptType!="ALL" ) mcZZFileName += "_" + leptType;
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  std::cout << "Opened mc file '" << mcTTbarFileName << "'." << std::endl;
  db->add_mcFile( mcTTbarFile, "TTJets_TuneD6T", "t#bar{t}", 30, 3002);


  db->set_legendTitle("");
  db->set_rebin(2);
  db->set_lumiNormalization(1000.);
  db->set_noStack((bool)false);
  db->drawHisto("ZZInvMass_hiMass", "ZZ Invariant Mass", "GeV/c^{2}", "Events");
  db->drawHisto("ZZInvMass_hiMass_QGlikeli", "ZZ Invariant Mass", "GeV/c^{2}", "Events");


  delete db;
  db = 0;

  return 0;

}  


void draw_vs_pt_plots( DrawBase* db, int nPtBins, Double_t* ptBins, const std::string& histoName, const std::string& axisName, const std::string& units, const std::string& instanceName, bool log ) {

  for( unsigned iPtBin=0; iPtBin<nPtBins-1; ++iPtBin ) {

    Double_t ptMin = ptBins[iPtBin];
    Double_t ptMax = ptBins[iPtBin+1];

    char histoName_thisBin[250];
    sprintf( histoName_thisBin, "%s_pt_%.0lf_%.0lf", histoName.c_str(), ptMin, ptMax);

    std::string histoName_thisBin_str(histoName_thisBin);

    char legendTitle[200];
    sprintf( legendTitle, "%.0lf < pt < %.0lf GeV/c", ptMin, ptMax);
    std::string legendTitle_str(legendTitle);
    db->set_legendTitle(legendTitle_str);

    db->drawHisto( histoName_thisBin_str, axisName, units, instanceName, log );

  } //for bins

} //draw_vs_pt_plots
