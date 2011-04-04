#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawWZ [(string)selType] [data=\"Nov4\"]" << std::endl;
    exit(23);
  }


  std::string selType(argv[1]);

  std::string leptType = "ALL";
//if( argc==3 ) {
//  std::string leptType_str(argv[2]);
//  leptType = leptType_str;
//}

  //std::string dataType="Nov4ReReco_PU";
  std::string dataType="Dec22ReReco";
  if( argc==3 ) {
    std::string dataType_str(argv[2]);
    dataType = dataType_str;
  }

  if( leptType!="ELE" && leptType!="MU" && leptType!="ALL" ) {
    std::cout << "Unknown leptType '" << leptType << "'. Only 'ELE', 'MU' and 'ALL' supported. Exiting." << std::endl;
    exit(17);
  }

  DrawBase* db = new DrawBase("WZ");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "WZPlots_DATA_" + dataType + "_" + selType + "_" + leptType;
  db->set_outputdir(outputdir_str);


  if( leptType=="MU" ) {
    std::string dataFileName = "WZ_Mu_" + dataType + "_" +selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "Muon_38x_35pb" );
  } else if( leptType=="ELE" ) { 
    std::string dataFileName = "WZ_Electron_" + dataType + "_" +selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "Electron_38x_35pb" );
  } else {
    std::string dataFileName = "WZ_EleMu_" + dataType + "_" +selType+"_"+leptType+".root";
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    db->add_dataFile( dataFile, "EleMu_38x_35pb" );
  }


  std::string mcZJetsFileName = "WZ_ZJets_alpgen_TuneZ2_Fall10";
//  std::string mcZJetsFileName = "WZ_DYJetsToLL_TuneZ2_M-50_madgraph_Fall10";
  //if( leptType!="ALL" ) mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  std::cout << "Opened mc file '" << mcZJetsFileName << "'." << std::endl;
  db->add_mcFile( mcZJetsFile, "ZJets_alpgen_TuneZ2_Fall10", "Z + Jets", 38, 3001);

  std::string mcVVFileName = "WZ_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10";
  //if( leptType!="ALL" ) mcZme += "_" + selType;
  mcVVFileName += "_" + selType;
  mcVVFileName += "_" + leptType;
  mcVVFileName += ".root";
  TFile* mcVVFile = TFile::Open(mcVVFileName.c_str());
  std::cout << "Opened mc file '" << mcVVFileName << "'." << std::endl;
  db->add_mcFile( mcVVFile, "VVtoAnything_TuneZ2", "ZZ/WZ/WW", kCyan+1, 3003);

  std::string mcTTbarFileName = "WZ_TTJets_TuneZ2_7TeV-madgraph-tauola_Fall10";
  //if( leptType!="ALL" ) mcZZFileName += "_" + leptType;
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  std::cout << "Opened mc file '" << mcTTbarFileName << "'." << std::endl;
  db->add_mcFile( mcTTbarFile, "TTJets_TuneD6T", "t#bar{t}", 30, 3002);



//if( leptType=="ELE" )
//  db->set_lumiNormalization(35.2);
//else if( leptType=="MU" )
//  db->set_lumiNormalization(33.3);
//else if( leptType=="ALL" )
//  db->set_lumiNormalization(34.3); //average

  db->set_lumiNormalization(34.);

  bool log = true;

//db->drawHisto("etaMu", "Muon Pseudorapidity", "", "Events");
//db->drawHisto("phiMu", "Muon Azimuth", "", "Events");
//db->drawHisto("ptEle", "Electron Transverse Momentum", "GeV/c", "Events");
//db->drawHisto("etaEle", "Electron Pseudorapidity", "", "Events");
//db->drawHisto("phiEle", "Electron Azimuth", "", "Events");

  db->set_rebin(5);
  db->drawHisto("ptZll", "Z Transverse Momentum", "GeV/c", "Events", log );
  db->drawHisto("ptZmumu", "Z#mu#mu Transverse Momentum", "GeV/c", "Events", log );
  db->drawHisto("ptZee", "Zee Transverse Momentum", "GeV/c", "Events", log );
  db->set_rebin(1);
  db->drawHisto("etaZll", "Z Pseudorapidity", "", "Events");
  db->drawHisto("phiZll", "Z Azimuth", "", "Events");
  db->drawHisto("etaZmumu", "Z#mu#mu Pseudorapidity", "", "Events");
  db->drawHisto("phiZmumu", "Z#mu#mu Azimuth", "", "Events");
  db->drawHisto("etaZee", "Zee Pseudorapidity", "", "Events");
  db->drawHisto("phiZee", "Zee Azimuth", "", "Events");

  db->set_rebin(20);
  db->drawHisto("mZllW", "Z(ll)W(jj) Mass", "GeV/c^{2}", "Events");
  db->drawHisto("mZmumuW", "Z(#mu#mu)W(jj) Mass", "GeV/c^{2}", "Events");
  db->drawHisto("mZeeW", "Z(ee)W(jj) Mass", "GeV/c^{2}", "Events");
  db->drawHisto("mZllW_kinfit", "Z(ll)W(jj) Mass", "GeV/c^{2}", "Events");
  db->drawHisto("mZmumuW_kinfit", "Z(#mu#mu)W(jj) Mass", "GeV/c^{2}", "Events");
  db->drawHisto("mZeeW_kinfit", "Z(ee)W(jj) Mass", "GeV/c^{2}", "Events");

  db->set_rebin(1);
  db->drawHisto("nJets", "Jet Multiplicity", "", "Events", log);
  db->drawHisto("nJets_inBump", "Jet Multiplicity", "", "Events", log);
  db->drawHisto("nJets_outBump", "Jet Multiplicity", "", "Events", log);

  db->set_rebin(10);
  db->drawHisto("ptMu1_inBump", "Lead Muon Transverse Momentum", "GeV/c", "Events");
  db->drawHisto("ptMu2_inBump", "Sublead Muon Transverse Momentum", "GeV/c", "Events");
  db->drawHisto("etaMu1_inBump", "Lead Muon Pseudorapidity", "", "Events");
  db->drawHisto("etaMu2_inBump", "Sublead Muon Pseudorapidity", "", "Events");
  db->drawHisto("phiMu1_inBump", "Lead Muon Azimuth", "", "Events");
  db->drawHisto("phiMu2_inBump", "Sublead Muon Azimuth", "", "Events");
  db->set_rebin(5);
  db->drawHisto("ptMu1_outBump", "Lead Muon Transverse Momentum", "GeV/c", "Events");
  db->drawHisto("ptMu2_outBump", "Sublead Muon Transverse Momentum", "GeV/c", "Events");
  db->drawHisto("etaMu1_outBump", "Lead Muon Pseaudorapidity", "", "Events");
  db->drawHisto("etaMu2_outBump", "Sublead Muon Pseaudorapidity", "", "Events");
  db->drawHisto("phiMu1_outBump", "Lead Muon Azimuth", "", "Events");
  db->drawHisto("phiMu2_outBump", "Sublead Muon Azimuth", "", "Events");
  db->set_rebin(2);
  db->drawHisto("mZll_inBump", "DiLepton Mass", "GeV/c^{2}", "Events", log );
  db->drawHisto("mZmumu_inBump", "DiMuon Mass", "GeV/c^{2}", "Events", log );
  db->drawHisto("mZee_inBump", "DiElectron Mass", "GeV/c^{2}", "Events", log );
  db->drawHisto("mZll_outBump", "DiLepton Mass", "GeV/c^{2}", "Events", log );
  db->drawHisto("mZmumu_outBump", "DiMuon Mass", "GeV/c^{2}", "Events", log );
  db->drawHisto("mZee_outBump", "DiElectron Mass", "GeV/c^{2}", "Events", log );

  db->set_rebin(10);
  db->drawHisto("deltaPhiMu1Met_inBump", "Lead Muon-ME_{T} #Delta#phi", "rad", "Events", log);
  db->drawHisto("deltaPhiMu1Met_outBump", "Lead Muon-ME_{T} #Delta#phi", "rad", "Events", log);

  db->drawHisto("eleEnergyFractionJet_inBump", "Electron Energy Fraction", "", "Events", log);
  db->drawHisto("eleEnergyFractionJet_outBump", "Electron Energy Fraction", "", "Events", log);
  db->drawHisto("muonEnergyFractionJet_inBump", "Muon Energy Fraction", "", "Events", log);
  db->drawHisto("muonEnergyFractionJet_outBump", "Muon Energy Fraction", "", "Events", log);
  db->drawHisto("mJet_inBump", "Jet Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mJet_outBump", "Jet Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZllJet_inBump", "Z-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZllJet_outBump", "Z-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZllJet_inBump", "Z-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZllJet_outBump", "Z-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZllMet_inBump", "Z-ME_{T} Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZllMet_outBump", "Z-ME_{T} Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZeeJet_inBump", "Zee-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZeeJet_outBump", "Zee-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZeeJet_inBump", "Zee-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZeeJet_outBump", "Zee-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZeeMet_inBump", "Zee-ME_{T} Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZeeMet_outBump", "Zee-ME_{T} Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZmumuJet_inBump", "Z#mu#mu-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZmumuJet_outBump", "Z#mu#mu-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZmumuJet_inBump", "Z#mu#mu-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZmumuJet_outBump", "Z#mu#mu-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZmumuMet_inBump", "Z#mu#mu-ME_{T} Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZmumuMet_outBump", "Z#mu#mu-ME_{T} Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("met_inBump", "Missing Transverse Energy", "GeV", "Events", log);
  db->drawHisto("met_outBump", "Missing Transverse Energy", "GeV", "Events", log);

  db->set_legendTitle("nJets = 0");
  db->drawHisto("ptZll_0jet", "Z Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZmumu_0jet", "Z#mu#mu Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZee_0jet", "Zee Transverse Momentum", "GeV/c", "Events", log);

  db->set_legendTitle("nJets = 1");
  db->drawHisto("ptZll_1jet", "Z Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZmumu_1jet", "Z#mu#mu Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZee_1jet", "Zee Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("mJet_1jet", "Jet Mass", "GeV/c^{2}", "Events", log);
//db->set_rebin(5);
  db->drawHisto("mZllJet_1jet", "Z-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZllJet_1jet", "Z-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZeeJet_1jet", "Zee-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZeeJet_1jet", "Zee-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZmumuJet_1jet", "Z#mu#mu-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZmumuJet_1jet", "Z#mu#mu-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
//db->set_rebin(2);
  db->drawHisto("met_1jet", "Missing Transverse Energy", "GeV", "Events", log);

  db->set_legendTitle("nJets = 2");
  db->drawHisto("ptZll_2jet", "Z Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZmumu_2jet", "Z#mu#mu Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZee_2jet", "Zee Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("mJet_2jet", "Jet Mass", "GeV/c^{2}", "Events", log);
//db->set_rebin(5);
  db->drawHisto("mZllJet_2jet", "Z-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZllJet_2jet", "Z-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZeeJet_2jet", "Zee-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZeeJet_2jet", "Zee-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZmumuJet_2jet", "Z#mu#mu-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZmumuJet_2jet", "Z#mu#mu-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
//db->set_rebin(2);
  db->drawHisto("met_2jet", "Missing Transverse Energy", "GeV", "Events", log);

/*
  db->set_legendTitle("nJets = 3");
  db->drawHisto("ptZll_3jet", "Z Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZmumu_3jet", "Z#mu#mu Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZee_3jet", "Zee Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("mJet_3jet", "Jet Mass", "GeV/c^{2}", "Events", log);
//db->set_rebin(5);
  db->drawHisto("mZllJet_3jet", "Z-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZllJet_3jet", "Z-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZeeJet_3jet", "Zee-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZeeJet_3jet", "Zee-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZmumuJet_3jet", "Z#mu#mu-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZmumuJet_3jet", "Z#mu#mu-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
//db->set_rebin(2);
  db->drawHisto("met_3jet", "Missing Transverse Energy", "GeV", "Events", log);

  db->set_legendTitle("nJets = 4");
  db->drawHisto("ptZll_4jet", "Z Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZmumu_4jet", "Z#mu#mu Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("ptZee_4jet", "Zee Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("mJet_4jet", "Jet Mass", "GeV/c^{2}", "Events", log);
//db->set_rebin(5);
  db->drawHisto("mZllJet_4jet", "Z-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZllJet_4jet", "Z-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZeeJet_4jet", "Zee-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZeeJet_4jet", "Zee-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mZmumuJet_4jet", "Z#mu#mu-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtZmumuJet_4jet", "Z#mu#mu-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
//db->set_rebin(2);
  db->drawHisto("met_4jet", "Missing Transverse Energy", "GeV", "Events", log);
*/

  delete db;
  db = 0;

  return 0;

}  


