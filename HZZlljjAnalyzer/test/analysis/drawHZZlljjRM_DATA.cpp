#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"


bool withSignal_=false;



int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 && argc !=4 ) {
    std::cout << "USAGE: ./drawHZZlljjRM [(string)selType] [(string)data_dataset=\"HR11\"] [(string) normType=\"LUMI\"]" << std::endl;
    exit(23);
  }

  std::string leptType = "ALL";

  std::string selType(argv[1]);

  std::string ZJetsMC = "madgraph";

//if( argc>=3 ) {
//  std::string PUType_str(argv[2]);
//  PUType = PUType_str;
//}

  std::string data_dataset = "HR11";
  if( argc>=3 ) {
    std::string data_dataset_str(argv[2]);
    data_dataset = data_dataset_str;
  }

  std::string normType = "LUMI";
  if( argc>=4 ) {
    std::string normType_str(argv[3]);
    normType = normType_str;
  }

  if( normType!="LUMI" && normType!="SHAPE" ) {
    std::cout << "Unknown normalization type: '" << normType << "'. Exiting." << std::endl;
    exit(191919);
  }


  DrawBase* db = new DrawBase("HZZlljjRM");
  db->set_pdf_aussi((bool)false);


  //std::string data_dataset = "DATA_Run2011A_v2_Sub2";
  //std::string data_dataset = "DATA_1fb";
  //std::string data_dataset = "DATA_EPS_FINAL";
  //std::string data_dataset = "DATA_EPS_FINAL_FULL";
  //std::string data_dataset = "DATA_LP11";
  //std::string data_dataset = "DATA_Run2011A_FULL";
  //std::string data_dataset = "DATA_Run2011B_v1";
  //std::string data_dataset = "DATA_HR11";

  

  std::string PUType = "HR11";
  if( data_dataset=="LP11" || data_dataset=="Run2011A_FULL" )
    PUType = "Run2011A";
  else if( data_dataset=="HR11" )
    PUType = "HR11";
  else if( data_dataset=="Run2011B_v1" )
    //PUType = "HR11";
    //PUType = "Run2011B_73pb";
    PUType = "Run2011B";
  else {
    std::cout << "Don't knwo what PUType to choose for data dataset '" << data_dataset << "'. Choosing HR11 default." << std::endl;
    PUType = "HR11";
  }

    


  std::string outputdir_str = "HZZlljjRMPlots_DATA_" + data_dataset + "_" + ZJetsMC;
  if( withSignal_ ) outputdir_str += "_plusSignal";
  outputdir_str += "_" + selType + "_PU" + PUType + "_" + leptType;
  if( normType=="SHAPE" ) outputdir_str += "_SHAPE";
  //outputdir_str += "_OLDPU";
  db->set_outputdir(outputdir_str);


  if( leptType=="MU" ) {
    std::string dataFileName = "HZZlljjRM_DoubleMu_Run2011A_v2_Sub2_"+selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "DoubleMu_Run2011A" );
  } else if( leptType=="ELE" ) { 
    std::string dataFileName = "HZZlljjRM_DoubleElectron_Run2011A_v2_Sub2_"+selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "DoubleElectron_Run2011A" );
  } else {
    std::string dataFileName = "HZZlljjRM_DATA_" + data_dataset + "_"+selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "DATA_Run2011A" );
  }

  std::string signalFileName = "HZZlljjRM_GluGluToHToZZTo2L2Q_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1";
  signalFileName += "_" + selType;
  signalFileName += "_PU" + PUType;
  signalFileName += "_" + leptType;
  signalFileName += ".root";
  TFile* signalFile = TFile::Open(signalFileName.c_str());
  if( selType=="presel" )
    db->add_mcFile_superimp( signalFile, "H400", "H(400) #times 100", 100., kRed-1);

  if( withSignal_ )
    db->add_mcFile( signalFile, "H400", "H(400)", kRed+2, 3004);

  std::string mcZJetsFileName;
  if( ZJetsMC=="alpgen" )
    mcZJetsFileName = "HZZlljjRM_ZJets_alpgen_TuneZ2_Spring11_v2";
  else if( ZJetsMC=="madgraph" )
    mcZJetsFileName = "HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
    //mcZJetsFileName = "HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_OLDPU";
  else {
    std::cout << "Unknown ZJetsMC '" << ZJetsMC << "'. Exiting." << std::endl;
    exit(13);
  }
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_PU" + PUType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  //db->add_mcFile( mcZJetsFile, "ZJets", "Z + Jets", 38, 3001);
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + Jets", 30, 3001);


  if( ZJetsMC=="alpgen" ) {
    std::string mcZBBFileName = "HZZlljjRM_ZBB_alpgen_TuneZ2_Spring11_v2";
    mcZBBFileName += "_" + selType;
    mcZBBFileName += "_PU" + PUType;
    mcZBBFileName += "_" + leptType;
    mcZBBFileName += ".root";
    TFile* mcZBBFile = TFile::Open(mcZBBFileName.c_str());
    db->add_mcFile( mcZBBFile, "ZBB_alpgen_TuneZ2_Spring11", "Z + bb", 40, 3003);

  //std::string mcZCCFileName = "HZZlljjRM_ZCC_alpgen_TuneZ2_Spring11_v2";
  //mcZCCFileName += "_" + selType;
  //mcZCCFileName += "_" + leptType;
  //mcZCCFileName += ".root";
  //TFile* mcZCCFile = TFile::Open(mcZCCFileName.c_str());
  //db->add_mcFile( mcZCCFile, "ZCC_alpgen_TuneZ2_Spring11", "Z + cc", 40, 3003);
  }

  //std::string mcVVFileName = "HZZlljjRM_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_OLDPU";
  //std::string mcVVFileName = "HZZlljjRM_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2";
  std::string mcVVFileName = "HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1";
  mcVVFileName += "_" + selType;
  mcVVFileName += "_PU" + PUType;
  mcVVFileName += "_" + leptType;
  mcVVFileName += ".root";
  TFile* mcVVFile = TFile::Open(mcVVFileName.c_str());
  //db->add_mcFile( mcVVFile, "VVtoAnything_TuneZ2", "ZZ/WZ/WW", kCyan+1, 3003);
  db->add_mcFile( mcVVFile, "VVtoAnything_TuneZ2", "ZZ/WZ/WW", 38, 3003);

  //std::string mcTTbarFileName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_OLDPU";
  //std::string mcTTbarFileName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2";
  std::string mcTTbarFileName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1";
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_PU" + PUType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  //db->add_mcFile( mcTTbarFile, "TTtW", "tt/tW", 30, 3002);
  db->add_mcFile( mcTTbarFile, "TTtW", "tt/tW", 39, 3002);






  if( normType=="LUMI" ) {

    if( data_dataset=="Run2011A_v2_Sub2" )
      db->set_lumiNormalization(175.);
    else if( data_dataset=="1fb" )
      db->set_lumiNormalization(859.);
    else if( data_dataset=="EPS" )
      db->set_lumiNormalization(960.); 
    else if( data_dataset=="EPS_FINAL" )
      db->set_lumiNormalization(1000.); 
    else if( data_dataset=="EPS_FINAL_FULL" )
      db->set_lumiNormalization(1143.); 
    else if( data_dataset=="EPS_FINAL_FULL_plusSingleMu" )
      db->set_lumiNormalization(1143.); 
    else if( data_dataset=="DoubleElectron_Aug05ReReco" )
      db->set_lumiNormalization(227.); 
    else if( data_dataset=="DoubleMu_Aug05ReReco" )
      db->set_lumiNormalization(285.); 
    else if( data_dataset=="EPS_FINAL_plusSingleMu" )
      db->set_lumiNormalization(1143.);
    else if( data_dataset=="LP11" )
      //db->set_lumiNormalization(1500.);
      db->set_lumiNormalization(1580.);
    else if( data_dataset=="Run2011A_FULL" )
      db->set_lumiNormalization(2100.);
    else if( data_dataset=="HR11" )
      db->set_lumiNormalization(4200.);
    else if( data_dataset=="Run2011B_v1" )
      db->set_lumiNormalization(2100.);

  } else { //shape

    if( data_dataset=="Run2011A_v2_Sub2" )
      db->set_lumi(175.);
    else if( data_dataset=="1fb" )
      db->set_lumi(859.);
    else if( data_dataset=="EPS" )
      db->set_lumi(960.); 
    else if( data_dataset=="EPS_FINAL" )
      db->set_lumi(1000.); 
  
    db->set_shapeNormalization();

  }



  bool log = true;

  db->drawHisto("nvertex", "Number of Reconstructed Vertexes", "", "Events", log);
  db->drawHisto("nvertex_PUW", "Number of Reconstructed Vertexes", "", "Events", log);
  db->drawHisto("nvertex_PUW_ave", "Number of Reconstructed Vertexes", "", "Events", log);

  db->set_getBinLabels(true);
  db->set_yAxisMaxScaleLog(50.);
  db->drawHisto("nEventsCategories_presel", "", "", "Events", log);
  db->set_yAxisMaxScaleLog(5.);
  db->set_getBinLabels(false);
  db->set_xAxisMin();

  db->drawHisto("rhoPF", "Particle Flow Energy Density (#rho)", "GeV", "Events", log);

  db->drawHisto("nJets_presel", "Jet Multiplicity (p_{T} > 30 GeV)", "", "Events", log);
  db->set_rebin(4);
  db->drawHisto("ptJet_all_presel", "Jet Transverse Momentum", "GeV", "Jets", log);

  db->set_rebin(10);
  db->drawHisto("mZjj_all_presel", "DiJet Invariant Mass", "GeV", "Jet Pairs", log);

  db->set_rebin(2);
  db->drawHisto("deltaRjj", "#DeltaR Between Jets (p_{T} > 30 GeV)", "", "Jet Pairs");
  db->set_rebin(1);
  db->drawHisto("deltaRjj_all_presel", "#DeltaR Between Jets (p_{T} > 30 GeV)", "", "Jet Pairs");
  db->drawHisto("deltaRll_presel", "#DeltaR Between Leptons", "", "Lepton Pairs");
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("etaLept1_presel", "Lead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaLept2_presel", "Sublead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaJet_all_presel", "Jet Pseudorapidity", "", "Jets");

  db->set_yAxisMaxScale( 1.1 );
  db->set_rebin(5);
  db->set_xAxisMax(250.);
  db->drawHisto("ptLept1_presel", "Lead Lepton p_{T}", "GeV", "Events", log);
  db->drawHisto("ptLept1", "Lead Lepton p_{T}", "GeV", "Events", log);
  db->set_xAxisMax(150.);
  db->drawHisto("ptLept2_presel", "Sublead Lepton p_{T}", "GeV", "Events", log);
  db->drawHisto("ptLept2", "Sublead Lepton p_{T}", "GeV", "Events", log);

  db->set_xAxisMax(250.);
  db->drawHisto("ptJet1", "Lead Jet p_{T}", "GeV", "Events", log);
  db->drawHisto("ptJet1_prekin", "Lead Jet p_{T}", "GeV", "Events", log);
  db->set_xAxisMax(150.);
  db->drawHisto("ptJet2", "Sublead Jet p_{T}", "GeV", "Events", log);
  db->drawHisto("ptJet2_prekin", "Sublead Jet p_{T}", "GeV", "Events", log);
  db->set_xAxisMax();
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("etaJet1", "Lead Jet Pseudorapidity", "", "Events", log);
  db->drawHisto("etaJet2", "Sublead Jet Pseudorapidity", "", "Events", log);
  db->drawHisto("tcheJet", "TCHE", "", "Events", log);
  db->drawHisto("tcheJet1", "Lead Jet TCHE", "", "Events", log);
  db->drawHisto("tcheJet2", "Sublead Jet TCHE", "", "Events", log);

  db->set_rebin(10);
  db->drawHisto("ptZll_presel", "Dilepton Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("ptZjj_all_presel", "Dijet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("ptZll", "Dilepton Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("ptZjj", "Dijet Transverse Momentum", "GeV", "Events", log);

  db->set_rebin(5);
  db->drawHisto("mZll", "m_{ll}", "GeV", "Events", log);
  db->set_yAxisMaxScale( 1.3 );
  db->set_rebin(10);
  db->drawHisto("mZjj", "m_{jj}", "GeV", "Events", log);
  //db->drawHisto("mZjj_nogluetag", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZjj_0btag", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZjj_1btag", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZjj_2btag", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("150 < m_{lljj} < 250 GeV");
  db->drawHisto("mZjj_loMass", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("250 < m_{lljj} < 400 GeV");
  db->drawHisto("mZjj_medMass", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("m_{lljj} > 400 GeV");
  db->drawHisto("mZjj_hiMass", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("");

  db->set_rebin(1);
  db->drawHisto("mZll_presel", "m_{ll}", "GeV", "Events", log);

  db->set_xAxisMax(1200.);
  db->set_rebin(20);
  db->drawHisto("mZZ_kinfit_hiMass_all", "m_{lljj}", "GeV", "Events", log);
  db->set_xAxisMax(800.);
  db->set_legendTitle("Gluon-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_gluetag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag", "m_{lljj}", "GeV", "Events", log);

  db->set_legendTitle("0 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_0btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_1btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_2btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("");
  db->set_xAxisMax();

  db->set_rebin(1);
  db->drawHisto("pfMet", "Particle Flow Missing E_{T}", "GeV", "Events", log);
  db->drawHisto("mEtSig", "ME_{T} / Sum E_{T}", "", "Events", log);
  db->drawHisto("metSignificance", "ME_{T} Significance", "", "Events", log);
  db->drawHisto("metSignificance_2btag", "ME_{T} Significance", "", "Events", log);

  db->set_rebin(1);
  db->drawHisto("nChargedJet1", "Leading Jet Charged Multiplicity", "", "Events");
  db->drawHisto("nNeutralJet1", "Leading Jet Neutral Multiplicity", "", "Events");
  db->drawHisto("ptDJet1", "Leading Jet p_{T}D", "", "Events");
  db->drawHisto("nChargedJet2", "Subleading Jet Charged Multiplicity", "", "Events");
  db->drawHisto("nNeutralJet2", "Subleading Jet Neutral Multiplicity", "", "Events");
  db->drawHisto("ptDJet2", "Subleading Jet p_{T}D", "", "Events");

  db->set_rebin(4);
  db->set_yAxisMaxScale(1.6);
  db->drawHisto("QGLikelihoodJet1", "Leading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodJet2", "Subleading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodProd", "Q-G Likelihood Product", "", "Events");
  db->set_legendTitle("282 < m_{lljj} < 330 GeV");
  db->drawHisto("QGLikelihoodJet1_MW300", "Leading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodJet2_MW300", "Subleading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodProd_MW300", "Q-G Likelihood Product", "", "Events");
  db->set_legendTitle("376 < m_{lljj} < 440 GeV");
  db->drawHisto("QGLikelihoodJet1_MW400", "Leading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodJet2_MW400", "Subleading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodProd_MW400", "Q-G Likelihood Product", "", "Events");
  db->set_legendTitle("");
  db->drawHisto("QGLikelihoodNoPUJet1", "Leading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodNoPUJet2", "Subleading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodNoPUProd", "Q-G Likelihood Product", "", "Events");
  db->set_yAxisMaxScale();

  db->set_yAxisMaxScale( 1.6 );
  db->set_rebin(3);
  db->drawHisto("cosThetaStar", "cos(#theta^{*})", "", "Events");
  db->drawHisto("cosTheta2", "cos(#theta_{2})", "", "Events");
  db->set_yAxisMaxScale( 1.8 );
  db->drawHisto("cosTheta1", "cos(#theta_{1})", "", "Events");
  db->drawHisto("phi", "#phi", "rad", "Events");
  db->drawHisto("phi1", "#phi_{1}", "rad", "Events");
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("helicityLD", "Angular Likelihood Discriminant", "", "Events");
  //db->drawHisto("helicityLD_nogluetag", "Angular Likelihood Discriminant", "", "Events");
  db->set_yAxisMaxScale();

  db->set_rebin(1);
  db->set_legendTitle("");


  //------------
  // MUON PLOTS:
  //------------

  if( data_dataset=="LP11" ) db->set_lumiNormalization(1530.);
  //if( data_dataset=="DATA_LP11" ) db->set_lumiNormalization(1615.);

  db->set_yAxisMaxScale( 1.6 );
  if( leptType=="ALL" || leptType=="MU" )
    db->drawHisto("mZmumu_presel", "m_{#mu#mu}", "GeV", "Events", log);
  db->set_yAxisMaxScale();

  db->set_rebin(10);
  db->set_legendTitle("Dimuon channel");
  db->drawHisto("mZjj_MU", "m_{jj}", "GeV", "Events", log);
  //db->drawHisto("mZjj_nogluetag_MU", "m_{jj}", "GeV", "Events", log);

  db->set_rebin(20);
  db->set_xAxisMax(800.);
  db->set_legendTitle("Gluon-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_gluetag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->drawHisto_fromTree("tree_passedEvents", "mZZ", "eventWeight*(nBTags==-1 && leptType==0 && !isSidebands)", 1500, 150., 1650., "mZZ_kinfit_hiMass_gluetag_MU_PROVA", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);

  db->set_legendTitle("0 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_0btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_1btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_2btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_xAxisMax();

  db->set_rebin(1);
  db->set_legendTitle("");



  //----------------
  // ELECTRON PLOTS:
  //----------------

  //if( data_dataset=="DATA_LP11" ) db->set_lumiNormalization(1556.);
  if( data_dataset=="LP11" ) db->set_lumiNormalization(1480.);

  
  db->set_yAxisMaxScale( 1.6 );
  if( leptType=="ALL" || leptType=="ELE" )
    db->drawHisto("mZee_presel", "m_{ee}", "GeV", "Events", log);
  db->set_yAxisMaxScale( );

  db->set_rebin(10);
  db->set_legendTitle("Dielectron channel");
  db->drawHisto("mZjj_ELE", "m_{jj}", "GeV", "Events", log);
  //db->drawHisto("mZjj_nogluetag_ELE", "m_{jj}", "GeV", "Events", log);

  db->set_rebin(20);
  db->set_xAxisMax(800.);
  db->set_legendTitle("Gluon-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_gluetag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag_ELE", "m_{eejj}", "GeV", "Events", log);

  db->set_legendTitle("0 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_0btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_1btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_2btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_xAxisMax();




  delete db;
  db = 0;

  return 0;

}  


