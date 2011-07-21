#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"



void draw_vs_pt_plots( DrawBase* db, int nPtBins, Double_t* ptBins, const std::string& histoName, const std::string& axisName, const std::string& units="", const std::string& instanceName="Entries", bool log=false );


int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawHWWlvjj [(string) selectionType] [(string) leptType=\"ALL\"]" << std::endl;
    exit(23);
  }

  std::string selType(argv[1]);

 // std::string zJets_dataset = "ZJets_alpgen";

  std::string leptType = "ALL";
  if( argc==3 ) {
    std::string leptType_str(argv[2]);
    leptType = leptType_str;
  }

  DrawBase* db = new DrawBase("HWWlvjj");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "HWWlvjjPlots_"+selType;
  //if( leptType!="ALL" ) outputdir_str += "_" + leptType;
  outputdir_str += "_" + leptType;
  db->set_outputdir(outputdir_str);

  //TFile* mcSignal_WW200 = TFile::Open("HWWlvjj_WW200_loose_ALL.root");
  //db->add_mcFile(mcSignal_WW200, "WW200", "HWW (200)", 2);
  //TFile* mcSignal_WW300 = TFile::Open("HWWlvjj_WW300_loose_ALL.root");
  //db->add_mcFile(mcSignal_WW300, "WW300", "HWW (300)", 46); 
  //TFile* mcSignal_WW400 = TFile::Open("HWWlvjj_WW400_loose_ALL.root");
  //db->add_mcFile(mcSignal_WW400, "WW400", "HWW (400)", 93);
  //TFile* mcSignal_WW500 = TFile::Open("HWWlvjj_WW500_helicity_ALL.root");
  //db->add_mcFile(mcSignal_WW500, "WW500", "HWW (500)", 98);


  TFile* mcBkg_WJets = TFile::Open("HWWlvjj_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  db->add_mcFile(mcBkg_WJets, "WJets", "W + Jets", 38);
  TFile* mcBkg_VV = TFile::Open("HWWlvjj_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  db->add_mcFile(mcBkg_VV, "VV", " Diboson ", 65);
  TFile* mcBkg_TT = TFile::Open("HWWlvjj_TT_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_3_helicity_ALL.root");
  db->add_mcFile(mcBkg_TT, "tt", " tt ", 30);
  TFile* mcBkg_DY = TFile::Open("HWWlvjj_DY_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  db->add_mcFile(mcBkg_DY, "DY", " DY ", 39);
  TFile* mcBkg_T = TFile::Open("HWWlvjj_TToBLNu_TuneZ2_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  db->add_mcFile(mcBkg_T, "top", " top ", 50);
  TFile* mcBkg_QCD = TFile::Open("HWWlvjj_QCD_EMEnriched_TuneZ2_7TeV-pythia6_3_helicity_ALL.root");
  db->add_mcFile(mcBkg_QCD, "QCD", " QCD ", 7);


  TFile* Data_ = TFile::Open("HWWlvjj_DATA_6july_helicity_ALL.root");
  db->add_dataFile(Data_, "Data", " Data ", 1);

  db->set_lumiNormalization( 1000. ); //1 fb-1 li somma e norm alla lumi
  //db->set_noStack( false);//(bool)true );
  //db->set_shapeNormalization();

  db->drawHisto( "lept1Eta", "Lepton Eta", "", "Events", false);
  db->drawHisto( "lept2Eta", "Neutrino Eta", "", "Events", false);
  db->drawHisto( "lept1Pt", "Lepton Pt", "GeV/c", "Events", false);
  db->drawHisto( "lept2Pt", "Neutrino Pt", "GeV/c", "Events", false);
  db->drawHisto( "Mdilept", "Dilepton invariant Mass", "GeV/c^{2}", "Events", false);
  db->drawHisto( "Jet1Eta", "Leading jet Eta", "", "Events", false);
  db->drawHisto( "Jet2Eta", "Subleading jet Eta", "", "Events", false);
  db->drawHisto( "Jet1Pt", "Leading jet Pt", "GeV/c", "Events", false);
  db->drawHisto( "Jet2Pt", "Subleading jet Pt", "GeV/c", "Events", false);
  db->drawHisto( "Mdijet", "Dijet invariant Mass", "GeV/c^{2}", "Events", false);
  db->drawHisto( "energyMet", "Missing Energy", "GeV", "Events", true);
   
  db->drawHisto( "lept1Eta_mu", "(Mu) Lepton Eta", "", "Events", false);
  db->drawHisto( "lept2Eta_mu", "(Mu) Neutrino Eta", "", "Events", false);
  db->drawHisto( "lept1Pt_mu", "(Mu) Lepton Pt", "GeV/c", "Events", false);
  db->drawHisto( "lept2Pt_mu", "(Mu) Neutrino Pt", "GeV/c", "Events", false);
  db->drawHisto( "Mdilept_mu", "(Mu) Dilepton invariant Mass", "GeV/c^{2}", "Events", false);
  db->drawHisto( "Jet1Eta_mu", "(Mu) Leading jet Eta", "", "Events", false);
  db->drawHisto( "Jet2Eta_mu", "(Mu) Subleading jet Eta", "", "Events", false);
  db->drawHisto( "Jet1Pt_mu", "(Mu) Leading jet Pt", "GeV/c", "Events", false);
  db->drawHisto( "Jet2Pt_mu", "(Mu) Subleading jet Pt", "GeV/c", "Events", false);
  db->drawHisto( "Mdijet_mu", "(Mu) Dijet invariant Mass", "GeV/c^{2}", "Events", false);
  db->drawHisto( "energyMet_mu", "(Mu) Missing Energy", "GeV", "Events", true);

  db->drawHisto( "lept1Eta_e", "(Ele) Lepton Eta", "", "Events", false);
  db->drawHisto( "lept2Eta_e", "(Ele) Neutrino Eta", "", "Events", false);
  db->drawHisto( "lept1Pt_e", "(Ele) Lepton Pt", "GeV/c", "Events", false);
  db->drawHisto( "lept2Pt_e", "(Ele) Neutrino Pt", "GeV/c", "Events", false);
  db->drawHisto( "Mdilept_e", "(Ele) Dilepton invariant Mass", "GeV/c^{2}", "Events", false);
  db->drawHisto( "Jet1Eta_e", "(Ele) Leading jet Eta", "", "Events", false);
  db->drawHisto( "Jet2Eta_e", "(Ele) Subleading jet Eta", "", "Events", false);
  db->drawHisto( "Jet1Pt_e", "(Ele) Leading jet Pt", "GeV/c", "Events", false);
  db->drawHisto( "Jet2Pt_e", "(Ele) Subleading jet Pt", "GeV/c", "Events", false);
  db->drawHisto( "Mdijet_e", "(Ele) Dijet invariant Mass", "GeV/c^{2}", "Events", false);
  db->drawHisto( "energyMet_e", "(Ele) Missing Energy", "GeV", "Events", true);

  db->set_rebin(2);
  db->set_xAxisMin(400.);
  db->set_yAxisMax(500);
  db->drawHisto( "mWW_kinfit", "WW Inv. Mass After all cuts", "GeV/c^{2}", "Events", false, 2);

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

// To Have this Label without this code
//  TPaveText* cmsLabel = db->get_labelCMS();
//  TPaveText* sqrtLabel = db->get_labelSqrt();
//  cmsLabel->Draw("same"); 
//  sqrtLabel->Draw("same");
