#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawHZZlljj [(string)selType] [(string) leptType=\"ALL\"]" << std::endl;
    exit(23);
  }


  std::string selType(argv[1]);

  std::string leptType = "ALL";
  if( argc==3 ) {
    std::string leptType_str(argv[2]);
    leptType = leptType_str;
  }

  if( leptType!="ELE" && leptType!="MU" && leptType!="ALL" ) {
    std::cout << "Unknown leptType '" << leptType << "'. Only 'ELE', 'MU' and 'ALL' supported. Exiting." << std::endl;
    exit(17);
  }

  DrawBase* db = new DrawBase("HZZlljj");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "HZZlljjPlots_DATA_" + selType + "_" + leptType;
  db->set_outputdir(outputdir_str);


  if( leptType=="MU" ) {
    std::string dataFileName = "HZZlljj_Mu_Nov4ReReco_PU_"+selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "Muon_38x_35pb" );
  } else if( leptType=="ELE" ) { 
    std::string dataFileName = "HZZlljj_Electron_Nov4ReReco_PU_"+selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "Electron_38x_35pb" );
  } else {
    std::string dataFileName = "HZZlljj_EleMu_Nov4ReReco_PU_"+selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "EleMu_38x_35pb" );
  }

//std::string mcSignal200FileName = "HZZlljj_SMHiggsToZZTo2L2Q_M-200_7TeV-jhu-pythia6";
//mcSignal200FileName += "_" + selType;
//mcSignal200FileName += "_" + leptType;
//mcSignal200FileName += ".root";
//TFile* mcSignal200File = TFile::Open(mcSignal200FileName.c_str());
//std::cout << "Opened mc file '" << mcSignal200FileName << "'." << std::endl;
//db->add_mcFile( mcSignal200File, "HZZ_qqll_gluonfusion_M200", "HZZlljj (200)", kRed+1);

  std::string mcZJetsFileName = "HZZlljj_ZJets_alpgen_TuneZ2_Fall10";
  //if( leptType!="ALL" ) mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  std::cout << "Opened mc file '" << mcZJetsFileName << "'." << std::endl;
  db->add_mcFile( mcZJetsFile, "ZJets_alpgen_TuneZ2_Fall10", "Z + Jets", 38, 3001);

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



  if( leptType=="ELE" )
    db->set_lumiNormalization(35.2);
  else if( leptType=="MU" )
    db->set_lumiNormalization(33.3);
  else if( leptType=="ALL" )
    db->set_lumiNormalization(34.3); //average

  bool log = true;

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
  db->drawHisto("ptLept1", "Lead Lepton Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptLept2", "Sublead Lepton Transverse Momentum", "GeV/c", "Events", log);

  db->drawHisto("mZll_presel", "Dilepton Invariant Mass", "GeV/c^{2}", "Events", log);
  if( leptType=="ALL" || leptType=="MU" )
    db->drawHisto("mZmumu_presel", "DiMuon Invariant Mass", "GeV/c^{2}", "Events", log);
  if( leptType=="ALL" || leptType=="ELE" )
    db->drawHisto("mZee_presel", "DiElectron Invariant Mass", "GeV/c^{2}", "Events", log);

  db->set_rebin(2);
  db->drawHisto("ptZll_presel", "Dilepton Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZjj_all_presel", "Dijet Transverse Momentum", "GeV/c", "Events", log);

  db->drawHisto("mZll", "Dilepton Invariant Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(4);
  db->drawHisto("mZjj", "Dijet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("ptZll", "Dilepton Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZjj", "Dijet Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("deltaRjj", "#Delta R Between Jets", "", "Events", log);
  db->set_rebin(4);
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("QGLikelihoodJet1", "Leading Jet Q-G Likelihood", "", "Events");
  db->drawHisto("QGLikelihoodJet2", "Subleading Jet Q-G Likelihood", "", "Events");
  db->drawHisto("QGLikelihood_normsJet1", "Leading Jet Q-G Likelihood", "", "Events");
  db->drawHisto("QGLikelihood_normsJet2", "Subleading Jet Q-G Likelihood", "", "Events");
  db->set_yAxisMaxScale( 1.4 );
  db->drawHisto("QGLikelihoodRevProd_norms", "Q-G Likelihood Product", "", "Events");

  db->set_yAxisMaxScale( 1.6 );
  db->set_rebin(3);
  db->drawHisto("cosThetaStar", "cos(#theta^{*})", "", "Events");
  db->drawHisto("cosTheta1", "cos(#theta_{1})", "", "Events");
  db->drawHisto("cosTheta2", "cos(#theta_{2})", "", "Events");
  db->drawHisto("phi", "#phi", "rad", "Events");
  db->drawHisto("phi1", "#phi_{1}", "rad", "Events");


  db->set_rebin(2);
  db->drawHisto("mZZ_medMass", "ZZ Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZZ_hiMass", "ZZ Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZZ_hiMass_QGlikeliRevProd_norms", "ZZ Invariant Mass", "GeV/c^{2}", "Events", log);

  delete db;
  db = 0;

  return 0;

}  


