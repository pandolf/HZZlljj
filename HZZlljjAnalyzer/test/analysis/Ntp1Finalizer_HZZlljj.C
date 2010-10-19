#include <TH2F.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <TRegexp.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "fitTools.C"
//#include "TFitConstraintM.h"
//#include "TFitParticleEtEtaPhi.h"
//#include "TKinFitter.h"
//#include "TFitConstraintM.cc"
//#include "TFitParticleEtEtaPhi.cc"
//#include "TKinFitter.cc"


std::string LEPT_TYPE;


double delta_phi(double phi1, double phi2) {

  double dphi = fabs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}


float delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}




bool DEBUG_ = false;

Double_t totalLumi=0.;
TChain* tree;


void addFile(const std::string& dataset);
std::vector<TH1F*> getResponseHistos(const std::string& name, int binArraySize, Double_t* ptBins);
void writeResponseHistos( TFile* file, std::vector<TH1F*> h1_response, std::string dirName );
Double_t ErrEt(Float_t Et, Float_t Eta);
Double_t ErrEta(Float_t Et, Float_t Eta);
Double_t ErrPhi(Float_t Et, Float_t Eta);


void finalize(const std::string& dataset, const std::string lept_type="ALL") {

  tree = new TChain("reducedTree");

  if( lept_type!="ALL" && lept_type!="ELE" && lept_type!="MU" ) {
    std::cout << "Lepton type: '" << lept_type << "' currently not supported. Please use 'ALL', 'ELE' or 'MU'." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9181);
  }


  LEPT_TYPE = lept_type;

  std::string infileName, treeName;


  if( dataset=="DATA_EG_37X" ) {

    addFile( "EG_Run2010A_Jul15thReReco_v1" );
    addFile( "EG_Run2010A_Jul26thReReco_v1" );

  } else if( dataset=="Electron_Run2010B" ) {

    addFile( "Electron_Run2010B_PromptReco_v2_runs146240_146733_prod2" );
    addFile( "Electron_Run2010B_PromptReco_v2_runs146733_147111" );

  } else if( dataset=="Run2010B_runs146240_146733" ) {

    addFile( "Electron_Run2010B_PromptReco_v2_runs146240_146733_prod2" );
    addFile( "MU_Run2010B_PromptReco_v2_runs146240_146733" );

  } else if( dataset=="Mu_upto147589" ) {

    addFile("Mu_Run2010B_PromptReco_v2_runs146734_147589");
    addFile("Mu_Run2010A-Sep17ReReco_v2_runs135821_144114");
    addFile("Mu_Run2010B_PromptReco_v2_runs146240_146733");

  } else if( dataset=="ELEMUCombined" ) {

    addFile("Mu_Run2010B_PromptReco_v2_runs146734_147589");
    addFile("Mu_Run2010A-Sep17ReReco_v2_runs135821_144114");
    addFile("Mu_Run2010B_PromptReco_v2_runs146240_146733");
    addFile("EG_upto146724");

  } else if( dataset=="ZJets_alpgen" ) {

    addFile( "Z0Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z1Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z1Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z1Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z1Jets_Pt800to1600-alpgen_Spring10" );
    addFile( "Z2Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z2Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z2Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z2Jets_Pt800to1600-alpgen_Spring10" );
    addFile( "Z3Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z3Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z3Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z3Jets_Pt800to1600-alpgen_Spring10" );
    addFile( "Z4Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z4Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z4Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z4Jets_Pt800to1600-alpgen_Spring10" );
    addFile( "Z5Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z5Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z5Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z5Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z1Jets_alpgen_Spring10" ) {

    addFile( "Z1Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z1Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z1Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z1Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z2Jets_alpgen_Spring10" ) {

    addFile( "Z2Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z2Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z2Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z2Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z3Jets_alpgen_Spring10" ) {

    addFile( "Z3Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z3Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z3Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z3Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z4Jets_alpgen_Spring10" ) {

    addFile( "Z4Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z4Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z4Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z4Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="Z5Jets_alpgen_Spring10" ) {

    addFile( "Z5Jets_Pt0to100-alpgen_Spring10" );
    addFile( "Z5Jets_Pt100to300-alpgen_Spring10" );
    addFile( "Z5Jets_Pt300to800-alpgen_Spring10" );
    addFile( "Z5Jets_Pt800to1600-alpgen_Spring10" );

  } else if( dataset=="all" ) {

    if( lept_type=="ELE" ) {
      finalize( "EG_upto146724", lept_type );
    } else if( lept_type=="MU") {
      //finalize( "Mu_upto147589", lept_type);
      finalize("Mu_upto147589", lept_type);
    } else if( lept_type=="ALL" ) {
      finalize("ELEMUCombined", lept_type);
    }

    finalize( "HZZ_qqll_gluonfusion_M130", lept_type );
    finalize( "HZZ_qqll_gluonfusion_M150", lept_type );
    finalize( "HZZ_qqll_gluonfusion_M200", lept_type );
    finalize( "HZZ_qqll_gluonfusion_M300", lept_type );
    finalize( "HZZ_qqll_gluonfusion_M400", lept_type );
    finalize( "HZZ_qqll_gluonfusion_M500", lept_type );
    finalize( "TTbar_2l_Spring10", lept_type );
    finalize( "ZZ_Spring10", lept_type );
    finalize( "ZJets_alpgen", lept_type );
    finalize( "ZJets_madgraph", lept_type );
    return;

  } else {
  
    addFile( dataset );

  }



  std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  if( dataset=="Run2010B_runs146240_146733" )
    h1_totalLumi->SetBinContent(1, 1220000.);
  else
    h1_totalLumi->SetBinContent(1, totalLumi);

  TH1F* h1_run = new TH1F("run", "", 15149, 132440, 147589);

  TH1D* h1_ptLept1 = new TH1D("ptLept1", "", 50, 10., 220.);
  h1_ptLept1->Sumw2();
  TH1D* h1_ptLept2 = new TH1D("ptLept2", "", 50, 10., 220.);
  h1_ptLept2->Sumw2();
  TH1D* h1_ptLept2OverLept1 = new TH1D("ptLept2OverLept1", "", 50, 0., 1.001);
  h1_ptLept2OverLept1->Sumw2();

  TH1D* h1_ptJetLead = new TH1D("ptJetLead", "", 50, 20., 220.);
  h1_ptJetLead->Sumw2();
  TH1D* h1_ptJetLead2 = new TH1D("ptJetLead2", "", 50, 20., 220.);
  h1_ptJetLead2->Sumw2();
  TH1D* h1_ptJetLead3 = new TH1D("ptJetLead3", "", 50, 20., 220.);
  h1_ptJetLead3->Sumw2();
  TH1D* h1_ptJetRecoil = new TH1D("ptJetRecoil", "", 50, 0., 200.);
  h1_ptJetRecoil->Sumw2();
  TH1D* h1_iJet1 = new TH1D("iJet1", "", 8, -0.5, 7.5);
  h1_iJet1->Sumw2();
  TH1D* h1_iJet2 = new TH1D("iJet2", "", 8, -0.5, 7.5);
  h1_iJet2->Sumw2();
  TH1D* h1_iJet2MinusiJet1 = new TH1D("iJet2MinusiJet1", "", 8, -0.5, 7.5);
  h1_iJet2MinusiJet1->Sumw2();
  TH1D* h1_iJet2PlusiJet1 = new TH1D("iJet2PlusiJet1", "", 8, -0.5, 7.5);
  h1_iJet2PlusiJet1->Sumw2();
  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 50, 20., 220.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 50, 20., 140.);
  h1_ptJet2->Sumw2();
  TH1D* h1_ptJet2Rel = new TH1D("ptJet2Rel", "", 50, 0., 0.1);
  h1_ptJet2Rel->Sumw2();
  TH1D* h1_ptJet2OverLead = new TH1D("ptJet2OverLead", "", 50, 0., 1.001);
  h1_ptJet2OverLead->Sumw2();
  TH1D* h1_ptJet2OverJet1 = new TH1D("ptJet2OverJet1", "", 50, 0., 1.001);
  h1_ptJet2OverJet1->Sumw2();
  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 50, -5., 5.);
  h1_etaJet1->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 50, -5., 5.);
  h1_etaJet2->Sumw2();
  TH1D* h1_RchJet1 = new TH1D("RchJet1", "", 20, 0., 1.001);
  h1_RchJet1->Sumw2();
  TH1D* h1_RchJet2 = new TH1D("RchJet2", "", 20, 0., 1.001);
  h1_RchJet2->Sumw2();
  TH1D* h1_massJet1 = new TH1D("massJet1", "", 50, 0., 100.);
  h1_massJet1->Sumw2();
  TH1D* h1_massJet2 = new TH1D("massJet2", "", 50, 0., 100.);
  h1_massJet2->Sumw2();

  TH2D* h2_etaPhi_map = new TH2D("etaPhi_map", "", 50, -5., 5., 50, -3.149, 3.149 );
  h2_etaPhi_map->Sumw2();
  TH2D* h2_etaPhi_map_cutOnH = new TH2D("etaPhi_map_cutOnH", "", 50, -5., 5., 50, -3.149, 3.149 );
  h2_etaPhi_map_cutOnH->Sumw2();

  TH1D* h1_ptZll = new TH1D("ptZll", "", 50, 0., 300.);
  h1_ptZll->Sumw2();
  TH1D* h1_pzZll = new TH1D("pzZll", "", 50, 0., 300.);
  h1_pzZll->Sumw2();
  TH1D* h1_ptZjj = new TH1D("ptZjj", "", 50, 0., 300.);
  h1_ptZjj->Sumw2();
  TH1D* h1_pzZjj = new TH1D("pzZjj", "", 50, 0., 300.);
  h1_pzZjj->Sumw2();
  TH1D* h1_ptHardestZ = new TH1D("ptHardestZ", "", 50, 0., 240.);
  h1_ptHardestZ->Sumw2();

  TH1D* h1_deltaRll = new TH1D("deltaRll", "", 20, 0.5, 5.);
  h1_deltaRll->Sumw2();
  TH1D* h1_deltaRll_Nm1 = new TH1D("deltaRll_Nm1", "", 20, 0.5, 5.);
  h1_deltaRll_Nm1->Sumw2();
  TH1D* h1_deltaRjj = new TH1D("deltaRjj", "", 20, 0.5, 5.);
  h1_deltaRjj->Sumw2();
  TH1D* h1_deltaRjj_Nm1 = new TH1D("deltaRjj_Nm1", "", 20, 0.5, 5.);
  h1_deltaRjj_Nm1->Sumw2();
  TH1D* h1_deltaRZZ = new TH1D("deltaRZZ", "", 50, 0.5, 5.);
  h1_deltaRZZ->Sumw2();
  TH1D* h1_deltaRZZ_Nm1 = new TH1D("deltaRZZ_Nm1", "", 50, 0.5, 5.);
  h1_deltaRZZ_Nm1->Sumw2();
  TH1D* h1_deltaEtaZZ = new TH1D("deltaEtaZZ", "", 50, -10., 10.);
  h1_deltaEtaZZ->Sumw2();
  TH1D* h1_deltaEtaAbsZZ = new TH1D("deltaEtaAbsZZ", "", 50, -5., 5.);
  h1_deltaEtaAbsZZ->Sumw2();
  TH1D* h1_deltaPhiZZ = new TH1D("deltaPhiZZ", "", 50, 1.5, 3.1416);
  h1_deltaPhiZZ->Sumw2();
  TH1D* h1_deltaPtZZ = new TH1D("deltaPtZZ", "", 50, 0., 100.);
  h1_deltaPtZZ->Sumw2();
  TH1D* h1_ptHiggs = new TH1D("ptHiggs", "", 50, 0., 100.);
  h1_ptHiggs->Sumw2();
  TH1D* h1_pzHiggs = new TH1D("pzHiggs", "", 50, 0., 1000.);
  h1_pzHiggs->Sumw2();
  TH1D* h1_etaHiggs = new TH1D("etaHiggs", "", 50, -5., 5.);
  h1_etaHiggs->Sumw2();

  TH1D* h1_ptFullSystem = new TH1D("ptFullSystem", "", 50, 0., 100.);
  h1_ptFullSystem->Sumw2();
  TH1D* h1_ptRecoilOverJet2 = new TH1D("ptRecoilOverJet2", "", 50, 0., 3.);
  h1_ptRecoilOverJet2->Sumw2();
  TH1D* h1_ptHiggsOverRecoil = new TH1D("ptHiggsOverRecoil", "", 50, 0., 3.);
  h1_ptHiggsOverRecoil->Sumw2();

  TH1D* h1_deltaPhi_HiggsRecoil = new TH1D("deltaPhi_HiggsRecoil", "", 50, 1.5, 3.1416);
  h1_deltaPhi_HiggsRecoil->Sumw2();
  TH1D* h1_deltaPhi_ZllRecoil = new TH1D("deltaPhi_ZllRecoil", "", 50, 0., 3.1416);
  h1_deltaPhi_ZllRecoil->Sumw2();
  TH1D* h1_deltaPhi_ZjjRecoil = new TH1D("deltaPhi_ZjjRecoil", "", 50, 0., 3.1416);
  h1_deltaPhi_ZjjRecoil->Sumw2();
  TH1D* h1_deltaR_ZjjRecoil = new TH1D("deltaR_ZjjRecoil", "", 50, 0., 5.);
  h1_deltaR_ZjjRecoil->Sumw2();

  TH1F* h1_pfMet = new TH1F("pfMet", "", 50, 0., 120.);
  h1_pfMet->Sumw2();
  TH1F* h1_pfMet_minusHiggs = new TH1F("pfMet_minusHiggs", "", 50, 0., 120.);
  h1_pfMet_minusHiggs->Sumw2();

  int nBins_invMass = 20;
  float invMassMin = 0.;
  float invMassMax = 120.;

  TH1D* h1_massZll = new TH1D("massZll", "", nBins_invMass, 55., invMassMax);
  h1_massZll->Sumw2();
  TH1D* h1_massZll_Nm1 = new TH1D("massZll_Nm1", "", nBins_invMass, 55., invMassMax);
  h1_massZll_Nm1->Sumw2();
  TH1D* h1_MuMuInvMass = new TH1D("MuMuInvMass", "", nBins_invMass, 40., invMassMax);
  h1_MuMuInvMass->Sumw2();
  TH1D* h1_EleEleInvMass = new TH1D("EleEleInvMass", "", nBins_invMass, 40., invMassMax);
  h1_EleEleInvMass->Sumw2();
  TH1D* h1_massZjj = new TH1D("massZjj", "", nBins_invMass, 20., 200.);
  h1_massZjj->Sumw2();
  TH1D* h1_massZjj_cutOnH = new TH1D("massZjj_cutOnH", "", nBins_invMass, 20., 200.);
  h1_massZjj_cutOnH->Sumw2();
  TH1D* h1_massZjj_Nm1 = new TH1D("massZjj_Nm1", "", nBins_invMass, 20., 200.);
  h1_massZjj_Nm1->Sumw2();
  TH1D* h1_massZjj_RchHIHI = new TH1D("massZjj_RchHIHI", "", nBins_invMass, 20., 200.);
  h1_massZjj_RchHIHI->Sumw2();
  TH1D* h1_massZjj_RchHILO = new TH1D("massZjj_RchHILO", "", nBins_invMass, 20., 200.);
  h1_massZjj_RchHILO->Sumw2();
  TH1D* h1_massZjj_RchLOLO = new TH1D("massZjj_RchLOLO", "", nBins_invMass, 20., 200.);
  h1_massZjj_RchLOLO->Sumw2();
  TH1D* h1_massZjj_Rch2_050 = new TH1D("massZjj_Rch2_050", "", nBins_invMass, 20., 200.);
  h1_massZjj_Rch2_050->Sumw2();
  TH1D* h1_massZjj_Rch2_5070 = new TH1D("massZjj_Rch2_5070", "", nBins_invMass, 20., 200.);
  h1_massZjj_Rch2_5070->Sumw2();
  TH1D* h1_massZjj_Rch2_70100 = new TH1D("massZjj_Rch2_70100", "", nBins_invMass, 20., 200.);
  h1_massZjj_Rch2_70100->Sumw2();

  TH1D* h1_ZZInvMass_loMass = new TH1D("ZZInvMass_loMass", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass->Sumw2();
  TH1D* h1_ZZInvMass_medMass = new TH1D("ZZInvMass_medMass", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass->Sumw2();
  TH1D* h1_ZZInvMass_medMass_FINEBINNING = new TH1D("ZZInvMass_medMass_FINEBINNING", "", 20000, 100., 350.);
  h1_ZZInvMass_medMass_FINEBINNING->Sumw2();
  TH1D* h1_ZZInvMass_hiMass = new TH1D("ZZInvMass_hiMass", "", nBins_invMass, 400., 600.);
  h1_ZZInvMass_hiMass->Sumw2();

  TH1D* h1_ZZInvMass_loMass_ZjjTag = new TH1D("ZZInvMass_loMass_ZjjTag", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass_ZjjTag->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_ZjjTag = new TH1D("ZZInvMass_hiMass_ZjjTag", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_ZjjTag->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_tightOLD = new TH1D("ZZInvMass_hiMass_fullSelection_tightOLD", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_fullSelection_tightOLD->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_tight = new TH1D("ZZInvMass_hiMass_fullSelection_tight", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_fullSelection_tight->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_medium = new TH1D("ZZInvMass_hiMass_fullSelection_medium", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_fullSelection_medium->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr = new TH1D("ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_loose = new TH1D("ZZInvMass_hiMass_fullSelection_loose", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_fullSelection_loose->Sumw2();

  TH1D* h1_ZZInvMass_medMass_fullSelection = new TH1D("ZZInvMass_medMass_fullSelection", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_FINEBINNING = new TH1D("ZZInvMass_medMass_fullSelection_FINEBINNING", "", 20000, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_FINEBINNING->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_nokin = new TH1D("ZZInvMass_medMass_fullSelection_nokin", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_nokin->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass = new TH1D("ZZInvMass_medMass_fullSelection_nokin_lowInvMass", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_lept20 = new TH1D("ZZInvMass_medMass_fullSelection_nokin_lowInvMass_lept20", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_lept20->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetLead = new TH1D("ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetLead", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetLead->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets = new TH1D("ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets_Rch = new TH1D("ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets_Rch", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets_Rch->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetPt = new TH1D("ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetPt", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetPt->Sumw2();
  TH1D* h1_ZZInvMass_medMass_fullSelection_tight = new TH1D("ZZInvMass_medMass_fullSelection_tight", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass_fullSelection_tight->Sumw2();

  TH1D* h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag = new TH1D("ZZInvMass_loMass_ZjjTag_ZllAntiTag", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag = new TH1D("ZZInvMass_hiMass_ZjjTag_ZllAntiTag", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag->Sumw2();

  TH1D* h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40 = new TH1D("ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40 = new TH1D("ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40->Sumw2();


  Double_t ptBins[16];
  fitTools::getBins_int( 16, ptBins, 20., 500.);

  TProfile* hp_ptJetGenMean = new TProfile("ptJetGenMean", "", 15, ptBins);
  std::vector<TH1F*> h1_response_vs_pt = getResponseHistos("response", 16, ptBins);
  std::vector<TH1F*> h1_response_vs_pt_Rch050 = getResponseHistos("response_Rch050", 16, ptBins);
  std::vector<TH1F*> h1_response_vs_pt_Rch5070 = getResponseHistos("response_Rch5070", 16, ptBins);
  std::vector<TH1F*> h1_response_vs_pt_Rch70100 = getResponseHistos("response_Rch70100", 16, ptBins);


  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Int_t LS;
  tree->SetBranchAddress("LS", &LS);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

  Float_t pfMet;
  tree->SetBranchAddress("epfMet", &pfMet);
  Float_t phiMet;
  tree->SetBranchAddress("phipfMet", &phiMet);


  int leptType;
  tree->SetBranchAddress("leptType", &leptType);

  Float_t eLept1;
  tree->SetBranchAddress("eLept1", &eLept1);
  Float_t ptLept1;
  tree->SetBranchAddress("ptLept1", &ptLept1);
  Float_t etaLept1;
  tree->SetBranchAddress("etaLept1", &etaLept1);
  Float_t phiLept1;
  tree->SetBranchAddress("phiLept1", &phiLept1);

  Float_t eLept2;
  tree->SetBranchAddress("eLept2", &eLept2);
  Float_t ptLept2;
  tree->SetBranchAddress("ptLept2", &ptLept2);
  Float_t etaLept2;
  tree->SetBranchAddress("etaLept2", &etaLept2);
  Float_t phiLept2;
  tree->SetBranchAddress("phiLept2", &phiLept2);

  Float_t eJetLead;
  tree->SetBranchAddress("eJetLead", &eJetLead);
  Float_t ptJetLead;
  tree->SetBranchAddress("ptJetLead", &ptJetLead);
  Float_t etaJetLead;
  tree->SetBranchAddress("etaJetLead", &etaJetLead);
  Float_t phiJetLead;
  tree->SetBranchAddress("phiJetLead", &phiJetLead);

  Float_t eJetLead2;
  tree->SetBranchAddress("eJetLead2", &eJetLead2);
  Float_t ptJetLead2;
  tree->SetBranchAddress("ptJetLead2", &ptJetLead2);
  Float_t etaJetLead2;
  tree->SetBranchAddress("etaJetLead2", &etaJetLead2);
  Float_t phiJetLead2;
  tree->SetBranchAddress("phiJetLead2", &phiJetLead2);

  Float_t eJetLead3;
  tree->SetBranchAddress("eJetLead3", &eJetLead3);
  Float_t ptJetLead3;
  tree->SetBranchAddress("ptJetLead3", &ptJetLead3);
  Float_t etaJetLead3;
  tree->SetBranchAddress("etaJetLead3", &etaJetLead3);
  Float_t phiJetLead3;
  tree->SetBranchAddress("phiJetLead3", &phiJetLead3);

  Float_t eJetRecoil;
  tree->SetBranchAddress("eJetRecoil", &eJetRecoil);
  Float_t ptJetRecoil;
  tree->SetBranchAddress("ptJetRecoil", &ptJetRecoil);
  Float_t etaJetRecoil;
  tree->SetBranchAddress("etaJetRecoil", &etaJetRecoil);
  Float_t phiJetRecoil;
  tree->SetBranchAddress("phiJetRecoil", &phiJetRecoil);

  Int_t iJet1;
  tree->SetBranchAddress("iJet1", &iJet1);
  Float_t eJet1;
  tree->SetBranchAddress("eJet1", &eJet1);
  Float_t ptJet1;
  tree->SetBranchAddress("ptJet1", &ptJet1);
  Float_t etaJet1;
  tree->SetBranchAddress("etaJet1", &etaJet1);
  Float_t phiJet1;
  tree->SetBranchAddress("phiJet1", &phiJet1);
  Float_t eJetGen1;
  tree->SetBranchAddress("eJetGen1", &eJetGen1);
  Float_t ptJetGen1;
  tree->SetBranchAddress("ptJetGen1", &ptJetGen1);
  Float_t etaJetGen1;
  tree->SetBranchAddress("etaJetGen1", &etaJetGen1);
  Float_t phiJetGen1;
  tree->SetBranchAddress("phiJetGen1", &phiJetGen1);
  Float_t ePart1;
  tree->SetBranchAddress("ePart1", &ePart1);
  Float_t ptPart1;
  tree->SetBranchAddress("ptPart1", &ptPart1);
  Float_t etaPart1;
  tree->SetBranchAddress("etaPart1", &etaPart1);
  Float_t phiPart1;
  tree->SetBranchAddress("phiPart1", &phiPart1);
  Float_t eChargedHadronsJet1;
  tree->SetBranchAddress("eChargedHadronsJet1", &eChargedHadronsJet1);

  Int_t iJet2;
  tree->SetBranchAddress("iJet2", &iJet2);
  Float_t eJet2;
  tree->SetBranchAddress("eJet2", &eJet2);
  Float_t ptJet2;
  tree->SetBranchAddress("ptJet2", &ptJet2);
  Float_t etaJet2;
  tree->SetBranchAddress("etaJet2", &etaJet2);
  Float_t phiJet2;
  tree->SetBranchAddress("phiJet2", &phiJet2);
  Float_t eJetGen2;
  tree->SetBranchAddress("eJetGen2", &eJetGen2);
  Float_t ptJetGen2;
  tree->SetBranchAddress("ptJetGen2", &ptJetGen2);
  Float_t etaJetGen2;
  tree->SetBranchAddress("etaJetGen2", &etaJetGen2);
  Float_t phiJetGen2;
  tree->SetBranchAddress("phiJetGen2", &phiJetGen2);
  Float_t ePart2;
  tree->SetBranchAddress("ePart2", &ePart2);
  Float_t ptPart2;
  tree->SetBranchAddress("ptPart2", &ptPart2);
  Float_t etaPart2;
  tree->SetBranchAddress("etaPart2", &etaPart2);
  Float_t phiPart2;
  tree->SetBranchAddress("phiPart2", &phiPart2);
  Float_t eChargedHadronsJet2;
  tree->SetBranchAddress("eChargedHadronsJet2", &eChargedHadronsJet2);



  float nEvents400_pre=0.;
  float nEvents400_pre_leptPt=0.;
  float nEvents400_pre_leptPt_leptMass =0.;
  float nEvents400_pre_leptPt_leptMass_jetPt=0.;
  float nEvents400_pre_leptPt_leptMass_jetPt_jetMass=0.;
  float nEvents400_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj=0.;

  float nEvents500_pre=0.;
  float nEvents500_pre_leptPt=0.;
  float nEvents500_pre_leptPt_leptMass =0.;
  float nEvents500_pre_leptPt_leptMass_jetPt=0.;
  float nEvents500_pre_leptPt_leptMass_jetPt_jetMass=0.;
  float nEvents500_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj=0.;


  int nEntries = tree->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    if( LEPT_TYPE!="ALL" ) {
      if( LEPT_TYPE=="ELE" && leptType==0 ) continue;
      if( LEPT_TYPE=="MU" && leptType==1 ) continue;
    }


    if( run>5 ) {
    
      std::map<int, std::map<int, std::vector<int> > >::iterator it;

      it = run_lumi_ev_map.find(run);


      if( it==run_lumi_ev_map.end() ) {

        std::vector<int> events;
        events.push_back(event);
        std::map<int, std::vector<int> > lumi_ev_map;
        lumi_ev_map.insert( std::pair<int,std::vector<int> >(LS, events));
        run_lumi_ev_map.insert( std::pair<int, std::map<int, std::vector<int> > > (run, lumi_ev_map) );

      } else { //run exists, look for LS


        std::map<int, std::vector<int> >::iterator it_LS;
        it_LS = it->second.find( LS );

        if( it_LS==(it->second.end())  ) {

          std::vector<int> events;
          events.push_back(event);
          it->second.insert( std::pair<int, std::vector<int> > (LS, events) );

        } else { //LS exists, look for event

          std::vector<int>::iterator ev;
          for( ev=it_LS->second.begin(); ev!=it_LS->second.end(); ++ev )
            if( *ev==event ) break;


          if( ev==it_LS->second.end() ) {

            it_LS->second.push_back(event);

          } else {

            std::cout << "DISCARDING DUPLICATE EVENT!! Run: " << run << " LS: " << LS << " event: " << event << std::endl;

            continue;

          }
        }
      }


      h1_run->Fill( run, eventWeight );
    
    } //if is not mc



    TLorentzVector lept1, lept2;
    lept1.SetPtEtaPhiE( ptLept1, etaLept1, phiLept1, eLept1 );
    lept2.SetPtEtaPhiE( ptLept2, etaLept2, phiLept2, eLept2 );

    if( lept1.Pt() < 15. || lept2.Pt()<15. ) continue;

    h1_ptLept1->Fill( lept1.Pt(), eventWeight );
    h1_ptLept2->Fill( lept2.Pt(), eventWeight );

    h1_ptLept2OverLept1->Fill( lept2.Pt()/lept1.Pt(), eventWeight );

    float deltaRll = lept1.DeltaR(lept2);
    //h1_deltaRll->Fill( deltaRll, eventWeight );


    TLorentzVector Zll = ( lept1 + lept2 );
    h1_ptZll->Fill( Zll.Pt(), eventWeight );
    h1_pzZll->Fill( Zll.Pz(), eventWeight );


    TLorentzVector jetLead;
    jetLead.SetPtEtaPhiE( ptJetLead, etaJetLead, phiJetLead, eJetLead );
    TLorentzVector jetLead2;
    jetLead2.SetPtEtaPhiE( ptJetLead2, etaJetLead2, phiJetLead2, eJetLead2 );

    TLorentzVector jetRecoil;
    jetRecoil.SetPtEtaPhiE( ptJetRecoil, etaJetRecoil, phiJetRecoil, eJetRecoil );

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiE( ptJet1, etaJet1, phiJet1, eJet1 );
    jet2.SetPtEtaPhiE( ptJet2, etaJet2, phiJet2, eJet2 );

    bool jet1_is_lead = false;
    if( iJet1==0 ) jet1_is_lead=true;
    float Rch1 = eChargedHadronsJet1/eJet1;
    float Rch2 = eChargedHadronsJet2/eJet2;
    if( fabs(jet1.Eta())<2.5 && Rch1<0.1 ) continue;
    if( fabs(jet2.Eta())<2.5 && Rch2<0.1 ) continue;
  //float Rnh1 = eNeutralHadronsJet1/eJet1;
  //float Rnh2 = eNeutralHadronsJet2/eJet2;
  //float Rgamma1 = ePhotonsJet1/eJet1;
  //float Rgamma2 = ePhotonsJet2/eJet2;

    TLorentzVector jetgen1, jetgen2;
    jetgen1.SetPtEtaPhiE( ptJetGen1, etaJetGen1, phiJetGen1, eJetGen1 );
    jetgen2.SetPtEtaPhiE( ptJetGen2, etaJetGen2, phiJetGen2, eJetGen2 );

    h1_ptJetLead->Fill( ptJetLead, eventWeight );
    h1_ptJetLead2->Fill( ptJetLead2, eventWeight );
    h1_ptJetLead3->Fill( ptJetLead3, eventWeight );

    h1_ptJetRecoil->Fill( jetRecoil.Pt(), eventWeight );


    h1_ptJet1->Fill( ptJet1, eventWeight );
    h1_ptJet2->Fill( ptJet2, eventWeight );

    float ptJet2Rel = ptJet2/ptJetLead;
    if( iJet2>1 ) ptJet2Rel /= ptJetLead2; 
    if( iJet2>2 ) ptJet2Rel /= ptJetLead3; 
    h1_ptJet2Rel->Fill( ptJet2Rel, eventWeight );
    h1_ptJet2OverJet1->Fill( ptJet2/ptJet1, eventWeight );
    h1_ptJet2OverLead->Fill( ptJet2/ptJetLead, eventWeight );

    h1_massJet1->Fill( jet1.M(), eventWeight );
    h1_massJet2->Fill( jet2.M(), eventWeight );

    h1_iJet1->Fill( iJet1, eventWeight );
    h1_iJet2->Fill( iJet2, eventWeight );
    h1_iJet2MinusiJet1->Fill( iJet2-iJet1, eventWeight );
    h1_iJet2PlusiJet1->Fill( iJet2+iJet1, eventWeight );


    int binMax = hp_ptJetGenMean->GetXaxis()->GetNbins();
 
    int theBin1 = hp_ptJetGenMean->FindBin( ptJetGen1 );
    if( theBin1>0 && theBin1 <= binMax ) {
      hp_ptJetGenMean->Fill( ptJetGen1, ptJetGen1, eventWeight );
      h1_response_vs_pt[theBin1-1]->Fill( ptJet1/ptJetGen1, eventWeight );
      if( Rch1 < 0.5 ) h1_response_vs_pt_Rch050[theBin1-1]->Fill( ptJet1/ptJetGen1, eventWeight );
      if( Rch1 < 0.7 ) h1_response_vs_pt_Rch5070[theBin1-1]->Fill( ptJet1/ptJetGen1, eventWeight );
      else h1_response_vs_pt_Rch70100[theBin1-1]->Fill( ptJet1/ptJetGen1, eventWeight );
    }

    int theBin2 = hp_ptJetGenMean->FindBin( ptJetGen2 );
    if( theBin2>0 && theBin2 <= binMax ) {
      hp_ptJetGenMean->Fill( ptJetGen2, ptJetGen2, eventWeight );
      h1_response_vs_pt[theBin2-1]->Fill( ptJet2/ptJetGen2, eventWeight );
      if( Rch2 < 0.5 ) h1_response_vs_pt_Rch050[theBin2-1]->Fill( ptJet2/ptJetGen2, eventWeight );
      if( Rch2 < 0.7 ) h1_response_vs_pt_Rch5070[theBin2-1]->Fill( ptJet2/ptJetGen2, eventWeight );
      else h1_response_vs_pt_Rch70100[theBin2-1]->Fill( ptJet2/ptJetGen2, eventWeight );
    }


    float deltaRjj = jet1.DeltaR(jet2);
    //h1_deltaRjj->Fill( deltaRjj, eventWeight );

    TLorentzVector Zjj = ( jet1 + jet2 );
    h1_ptZjj->Fill( Zjj.Pt(), eventWeight );
    h1_pzZjj->Fill( Zjj.Pz(), eventWeight );
    if( Rch2 < 0.5 ) h1_massZjj_Rch2_050->Fill( Zjj.M(), eventWeight );
    else if( Rch2 < 0.7 ) h1_massZjj_Rch2_5070->Fill( Zjj.M(), eventWeight );
    else h1_massZjj_Rch2_70100->Fill( Zjj.M(), eventWeight );

    bool isSignal=false;
    TString dataset_tstr(dataset);
    TRegexp re("HZZ_qqll");
    if( dataset_tstr.Contains(re) ) isSignal=true;
    if( isSignal ) {
      TLorentzVector part1, part2;
      if( ptPart1!=0. ) part1.SetPtEtaPhiE( ptPart1, etaPart1, phiPart1, ePart1 );
      else part1.SetPtEtaPhiE(0., 20., 0., 0.);
      if( ptPart2!=0. ) part2.SetPtEtaPhiE( ptPart2, etaPart2, phiPart2, ePart2 );
      else part2.SetPtEtaPhiE(0., 20., 0., 0.);

      if( jet1.DeltaR(part1) < 0.25 && jet2.DeltaR(part2) < 0.25 ) {
        if( Rch1>0.7 && Rch2>0.7 ) h1_massZjj_RchHIHI->Fill( Zjj.M(), eventWeight );
        else if( Rch1>0.7 || Rch2>0.7 ) h1_massZjj_RchHILO->Fill( Zjj.M(), eventWeight );
        else  h1_massZjj_RchLOLO->Fill( Zjj.M(), eventWeight );
      }
    } //if signal


    h1_deltaPhi_ZllRecoil->Fill( Zll.DeltaPhi(jetRecoil), eventWeight );
    h1_deltaPhi_ZjjRecoil->Fill( Zjj.DeltaPhi(jetRecoil), eventWeight );
    h1_deltaR_ZjjRecoil->Fill( Zjj.DeltaR(jetRecoil), eventWeight );


    //constrain the Z->jj mass to the pdg value
    TLorentzVector Zjj_constr;
    Zjj_constr.SetXYZM( Zjj.Px(), Zjj.Py(), Zjj.Pz(), 91.1876);
 
    float ptHardestZ;
    if( Zjj.Pt() > Zll.Pt() )
      ptHardestZ = Zjj.Pt();
    else
      ptHardestZ = Zll.Pt();

    h1_ptHardestZ->Fill( ptHardestZ, eventWeight );

    float ptHardestZ_constr;
    if( Zjj_constr.Pt() > Zll.Pt() )
      ptHardestZ_constr = Zjj_constr.Pt();
    else
      ptHardestZ_constr = Zll.Pt();


    float deltaRZZ = Zll.DeltaR(Zjj);
    //h1_deltaRZZ->Fill( deltaRZZ, eventWeight );
    float deltaRZZ_constr = Zll.DeltaR(Zjj_constr);
    h1_deltaEtaZZ->Fill( Zll.Eta() - Zjj.Eta(), eventWeight );
    h1_deltaEtaAbsZZ->Fill( fabs(Zll.Eta()) - fabs(Zjj.Eta()), eventWeight );
    h1_deltaPhiZZ->Fill( Zll.DeltaPhi(Zjj), eventWeight );
    h1_deltaPtZZ->Fill( Zll.Pt()-Zjj.Pt(), eventWeight );

    TLorentzVector ZZ = Zll + Zjj;
    TLorentzVector ZZ_constr = Zll + Zjj_constr;

    float ptHiggs = ZZ.Pt();
    h1_ptHiggs->Fill( ptHiggs, eventWeight );
    h1_pzHiggs->Fill( ZZ.Pz(), eventWeight );
    float ptHiggs_constr = ZZ_constr.Pt();
    h1_etaHiggs->Fill( ZZ.Eta(), eventWeight );

    // compute met using leptons and jets and recoil jet:
    TLorentzVector fullSystem = lept1 + lept2 + jet1 + jet2 + jetRecoil;
  
    h1_ptFullSystem->Fill( fullSystem.Pt(), eventWeight);
    h1_ptRecoilOverJet2->Fill( jetRecoil.Pt()/jet2.Pt(), eventWeight);
    h1_ptHiggsOverRecoil->Fill( ZZ.Pt()/jetRecoil.Pt(), eventWeight);
    h1_deltaPhi_HiggsRecoil->Fill( ZZ.DeltaPhi(jetRecoil), eventWeight );

    // ------------------------
    //   KINEMATIC FIT: BEGIN
    // ------------------------
//  TLorentzVector chjet1( jet1*Rch1 );
//  TLorentzVector nhjet1( jet1*Rnh1 );
//  TLorentzVector gammajet1( jet1*Rgamma1 );
//  TLorentzVector chjet2( jet2*Rch2 );
//  TLorentzVector nhjet2( jet2*Rnh2 );
//  TLorentzVector gammajet2( jet2*Rgamma2 );

//  TMatrixD m_chjet1(3,3);
//  TMatrixD m_nhjet1(3,3);
//  TMatrixD m_gammajet1(3,3);
//  TMatrixD m_chjet2(3,3);
//  TMatrixD m_nhjet2(3,3);
//  TMatrixD m_gammajet2(3,3);

//  m_chjet1.Zero();
//  m_nhjet1.Zero();
//  m_gammajet1.Zero();
//  m_chjet2.Zero();
//  m_nhjet2.Zero();
//  m_gammajet2.Zero();

//  m_chJet1(0,0) = ErrEt (chjet1.Et(), chjet1.Eta()); // et
//  m_chJet1(1,1) = ErrEta(chjet1.Et(), chjet1.Eta()); // eta
//  m_chJet1(2,2) = ErrPhi(chjet1.Et(), chjet1.Eta()); // phi

//  TMatrixD m_jet1(3,3);
//  TMatrixD m_jet2(3,3);

//  m_jet1(0,0) = ErrEt (jet1.Et(), jet1.Eta()); // et
//  m_jet1(1,1) = ErrEta(jet1.Et(), jet1.Eta()); // eta
//  m_jet1(2,2) = ErrPhi(jet1.Et(), jet1.Eta()); // phi
//  m_jet2(0,0) = ErrEt (jet2.Et(), jet2.Eta()); // et
//  m_jet2(1,1) = ErrEta(jet2.Et(), jet2.Eta()); // eta
//  m_jet2(2,2) = ErrPhi(jet2.Et(), jet2.Eta()); // phi

//  TFitParticleEtEtaPhi *fitJet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &jet1, &m_jet2 );
//  TFitParticleEtEtaPhi *fitJet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &jet2, &m_jet2 );
    


    // ------------------------
    //   KINEMATIC FIT: END
    // ------------------------


    if( leptType==0 )
      h1_MuMuInvMass->Fill( Zll.M(), eventWeight );
    else if( leptType==1 )
      h1_EleEleInvMass->Fill( Zll.M(), eventWeight );
    else
      std::cout << "WARNING!! found incredible leptType: '" << leptType << "'." << std::endl;
    h1_massZll->Fill( Zll.M(), eventWeight );

    // require Z->ll mass:
    if( Zll.M()<86. || Zll.M()>96. ) continue;

    h2_etaPhi_map->Fill( etaJet1, phiJet1, eventWeight );
    h2_etaPhi_map->Fill( etaJet2, phiJet2, eventWeight );
    if( ZZ.M() < 200. && ZZ.M()>180. ) {
      h2_etaPhi_map_cutOnH->Fill( etaJet1, phiJet1, eventWeight );
      h2_etaPhi_map_cutOnH->Fill( etaJet2, phiJet2, eventWeight );
    }

    h1_pfMet->Fill( pfMet, eventWeight );
    float pxMet = pfMet*cos(phiMet);
    float pyMet = pfMet*sin(phiMet);
    float pxHiggs = ptHiggs*cos(ZZ.Phi());
    float pyHiggs = ptHiggs*sin(ZZ.Phi());
    float pxMet_minusHiggs = pxMet-pxHiggs; 
    float pyMet_minusHiggs = pyMet-pyHiggs; 
    float pfMet_minusHiggs = sqrt( pxMet_minusHiggs*pxMet_minusHiggs + pyMet_minusHiggs*pyMet_minusHiggs);
    h1_pfMet_minusHiggs->Fill( pfMet_minusHiggs, eventWeight );


    // ----------------------------------
    //   LOW MASS ANALYSIS (mH < 180)
    // ----------------------------------

    h1_ZZInvMass_loMass->Fill( ZZ.M(), eventWeight );

    if( Zjj.M() > 80. && Zjj.M()<110. && ptJet1>35.) {

      h1_ZZInvMass_loMass_ZjjTag->Fill( ZZ.M(), eventWeight );

      if( Zll.M() < 70. ) {

        h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag->Fill( ZZ.M(), eventWeight );

        if( Rch1>0.4 && Rch2>0.4 ) {

          h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40->Fill( ZZ.M(), eventWeight );
 
        }
      }
    }


    // ----------------------------------
    //   HIGH MASS ANALYSIS (mH > 180)
    // ----------------------------------

    h1_ZZInvMass_medMass->Fill( ZZ.M(), eventWeight );
    h1_ZZInvMass_medMass_FINEBINNING->Fill( ZZ.M(), eventWeight );
    h1_ZZInvMass_hiMass->Fill( ZZ.M(), eventWeight );


    h1_massZjj->Fill( Zjj.M(), eventWeight );
    h1_deltaRll->Fill( deltaRll, eventWeight );
    h1_deltaRll->Fill( deltaRjj, eventWeight );
    h1_deltaRZZ->Fill( deltaRjj, eventWeight );

    if( deltaRjj>2. && deltaRll > 2. && deltaRZZ < 3.5 ) h1_massZjj_Nm1->Fill( Zjj.M(), eventWeight );
    if( Zjj.M()>75. && Zjj.M()<100. && deltaRjj>2. && deltaRll > 2. && deltaRZZ < 3.5 ) h1_massZll_Nm1->Fill( Zll.M(), eventWeight );
    if( Zjj.M()>75. && Zjj.M()<100. && deltaRjj>2. && deltaRZZ < 3.5 ) h1_deltaRll_Nm1->Fill( deltaRll, eventWeight );
    if( Zjj.M()>75. && Zjj.M()<100. && deltaRll>2. && deltaRZZ < 3.5 ) h1_deltaRll_Nm1->Fill( deltaRjj, eventWeight );
    if( Zjj.M()>75. && Zjj.M()<100. && deltaRll>2. && deltaRjj > 2. ) h1_deltaRZZ_Nm1->Fill( deltaRjj, eventWeight );


    if( ZZ.M()<205. && ZZ.M()>180. ) {
      h1_massZjj_cutOnH->Fill( Zjj.M(), eventWeight );
    }

    if( Zjj.M()<150. ) {

      if( jet1.Pt()>30. && jet2.Pt()>30. )
        h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetPt->Fill( ZZ.M(), eventWeight );

      h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass->Fill( ZZ.M(), eventWeight );
      if( fabs(etaJet1)<2. && fabs(etaJet2)<2. ) {
        h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets->Fill( ZZ.M(), eventWeight );
        if( Rch1>0.4 && Rch2>0.4 ) {
          h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets_Rch->Fill( ZZ.M(), eventWeight );
        }
      }

      if( lept1.Pt()>20. && lept2.Pt()>20. ) h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_lept20->Fill( ZZ.M(), eventWeight );
      
      TLorentzVector Zjj_lead = jetLead + jetLead2;
      TLorentzVector ZZ_lead = Zjj_lead + Zll;
      h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetLead->Fill( ZZ_lead.M(), eventWeight );


    }

    if( Zjj.M() > 70. && Zjj.M()<110. ) {

//std::cout << "run: " << run << " event: " << event << std::endl;
      h1_ZZInvMass_medMass_fullSelection_nokin->Fill( ZZ.M(), eventWeight );
      h1_etaJet1->Fill( etaJet1, eventWeight );
      h1_etaJet2->Fill( etaJet2, eventWeight );


      if( deltaRjj>2. && deltaRll > 2. && deltaRZZ < 3.5 ) {

        h1_ZZInvMass_medMass_fullSelection->Fill( ZZ.M(), eventWeight );
        h1_ZZInvMass_medMass_fullSelection_FINEBINNING->Fill( ZZ.M(), eventWeight );
        h1_RchJet1->Fill( Rch1, eventWeight );
        h1_RchJet2->Fill( Rch2, eventWeight );
  ////if( Zjj.M() > 80. && Zjj.M() < 100 ) {
    if( ZZ.M() > 190. && ZZ.M() < 200 ) {
      std::cout << "Run: " << run << " LS: " << LS << " event: " << event << std::endl;
      std::cout << " pt1: " << jet1.Pt() << " eta1: " << jet1.Eta() << " phi1: " << jet1.Phi() << " Rch1: " << Rch1 << std::endl;
      std::cout << " pt2: " << jet2.Pt() << " eta2: " << jet2.Eta() << " phi2: " << jet2.Phi() << " Rch2: " << Rch2 << std::endl;
      std::cout << " M(jj): " << Zjj.M() << std::endl;
    }

        if( ptJet1>45. && ptJet2>30. ) {

          h1_ZZInvMass_medMass_fullSelection_tight->Fill( ZZ.M(), eventWeight );

        }

      }

    }


    if( Zjj.M() > 80. && Zjj.M()<110. && ptJet1>115. && ptJet2>55. && deltaRjj<1.5 && ptHiggs>30.) {

      h1_ZZInvMass_hiMass_ZjjTag->Fill( ZZ.M(), eventWeight );

    }

    if( ptJetLead > 100. &&
        ptJetLead2 > 50. &&
        ptLept1 > 80. &&
        ptLept2 > 50. &&
        Zjj.M() > 80. &&
        Zjj.M() < 110. &&
        Zll.M() > 86. &&
        Zll.M() < 96. &&
        ptHardestZ > 120. &&
        ZZ.Pt() > 20. &&
        deltaRZZ > 2.8 ) {

      h1_ZZInvMass_hiMass_fullSelection_tightOLD->Fill( ZZ.M(), eventWeight );

    }

    if( ZZ.M()>390. && ZZ.M()<420. ) nEvents400_pre += eventWeight;
    if( ZZ.M()>460. && ZZ.M()<540. ) nEvents500_pre += eventWeight;

    if( ptLept1 > 100. &&
        ptLept2 > 50. ) {

      if( ZZ.M()>390. && ZZ.M()<420. ) nEvents400_pre_leptPt += eventWeight;
      if( ZZ.M()>460. && ZZ.M()<540. ) nEvents500_pre_leptPt += eventWeight;

      if( Zll.M() > 86. &&
          Zll.M() < 96. ) {

        if( ZZ.M()>390. && ZZ.M()<420. ) nEvents400_pre_leptPt_leptMass += eventWeight;
        if( ZZ.M()>460. && ZZ.M()<540. ) nEvents500_pre_leptPt_leptMass += eventWeight;

        if( ptJetLead > 100. &&
            ptJetLead2 > 60. ) {
      
          if( ZZ.M()>390. && ZZ.M()<420. ) nEvents400_pre_leptPt_leptMass_jetPt += eventWeight;
          if( ZZ.M()>460. && ZZ.M()<540. ) nEvents500_pre_leptPt_leptMass_jetPt += eventWeight;
      
            //ptJet2 > 0.15*ptJetLead &&
            //ptJetLead3 > 30. &&
            //iJet2 > 1 &&
            //ptJetRecoil > 20.&&
      
      
          if( Zjj.M() > 80. &&
              Zjj.M() < 105. ) {
        
            if( ZZ.M()>390. && ZZ.M()<420. ) nEvents400_pre_leptPt_leptMass_jetPt_jetMass += eventWeight;
            if( ZZ.M()>460. && ZZ.M()<540. ) nEvents500_pre_leptPt_leptMass_jetPt_jetMass += eventWeight;

              if( deltaRjj < 1.5 ) {

                if( ZZ.M()>390. && ZZ.M()<420. ) nEvents400_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj += eventWeight;
                if( ZZ.M()>460. && ZZ.M()<540. ) nEvents500_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj += eventWeight;

                  h1_ZZInvMass_hiMass_fullSelection_tight->Fill( ZZ.M(), eventWeight );

            }
          }
        }
      }
    }


    if( ptJetLead > 80. &&
        ptJetRecoil > 15. &&
        ptLept1 > 80. &&
        ptLept2 > 50. &&
        Zjj.M() > 80. &&
        Zjj.M() < 105. &&
        Zll.M() > 86. &&
        Zll.M() < 96. &&
        deltaRll < 1.5 &&
        deltaRjj < 1.5 &&
        deltaRZZ < 3.5 ) {

      h1_ZZInvMass_hiMass_fullSelection_medium->Fill( ZZ.M(), eventWeight );

    }

    if( ptJet1 > 80. &&
        ptJet2 > 40. &&
        ptLept1 > 80. &&
        ptLept2 > 50. &&
        Zjj.M() > 80. &&
        Zjj.M() < 110. &&
        Zll.M() > 86. &&
        Zll.M() < 96. &&
        ptHardestZ_constr > 100. &&
        deltaRjj < 2. &&
        ptHiggs_constr > 20. &&
        deltaRZZ_constr > 2.8 ) {

      h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr->Fill( ZZ_constr.M(), eventWeight );

    }


    if( ptJetLead > 60. &&
        ptJetRecoil > 15. &&
        ptLept1 > 65. &&
        ptLept2 > 30. &&
        Zjj.M() > 80. &&
        Zjj.M() < 105. &&
        Zll.M() > 86. &&
        Zll.M() < 96. &&
        deltaRll < 1.5 &&
        deltaRZZ < 3.5 ) {

      h1_ZZInvMass_hiMass_fullSelection_loose->Fill( ZZ.M(), eventWeight );

    }


  } //for entries


  std::string ofs400_name = "effTable400_tight_"+dataset+".txt";
  ofstream ofs400(ofs400_name.c_str());
  ofs400 << "DATASET\tPreselection\tLepton pt\tLepton mass\tjet pt\tdijet mass\tjet deltaR" << std::endl;
  ofs400 << dataset << "\t"
      << nEvents400_pre*1000. << "\t"
      << nEvents400_pre_leptPt*1000. << "\t"
      << nEvents400_pre_leptPt_leptMass*1000. << "\t"
      << nEvents400_pre_leptPt_leptMass_jetPt*1000. << "\t"
      << nEvents400_pre_leptPt_leptMass_jetPt_jetMass*1000. << "\t"
      << nEvents400_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj*1000. << "\t"
      << std::endl;

  std::string ofs500_name = "effTable500_tight_"+dataset+".txt";
  ofstream ofs500(ofs500_name.c_str());
  ofs500 << "DATASET\tPreselection\tLepton pt\tLepton mass\tjet pt\tdijet mass\tjet deltaR" << std::endl;
  ofs500 << dataset << "\t"
      << nEvents500_pre*1000. << "\t"
      << nEvents500_pre_leptPt*1000. << "\t"
      << nEvents500_pre_leptPt_leptMass*1000. << "\t"
      << nEvents500_pre_leptPt_leptMass_jetPt*1000. << "\t"
      << nEvents500_pre_leptPt_leptMass_jetPt_jetMass*1000. << "\t"
      << nEvents500_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj*1000. << "\t"
      << std::endl;


  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaHZZlljj_"+dataset;
  else {
   if(dataset!="") outfileName = "HZZlljj_"+dataset;
   else outfileName = "HZZlljj";
  }

  if( LEPT_TYPE!="ALL" ) outfileName += "_" + LEPT_TYPE;
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();

  h1_totalLumi->Write();
  h1_run->Write();

  h1_ptLept1->Write();
  h1_ptLept2->Write();
  h1_ptLept2OverLept1->Write();

  h1_ptJetLead->Write();
  h1_ptJetLead2->Write();
  h1_ptJetLead3->Write();
  h1_ptJetRecoil->Write();
  h1_ptJet1->Write();
  h1_ptJet2->Write();
  h1_ptJet2Rel->Write();
  h1_ptJet2OverLead->Write();
  h1_iJet1->Write();
  h1_iJet2->Write();
  h1_iJet2MinusiJet1->Write();
  h1_iJet2PlusiJet1->Write();
  h1_ptJet2OverJet1->Write();
  h1_etaJet1->Write();
  h1_etaJet2->Write();
  h1_RchJet1->Write();
  h1_RchJet2->Write();
  h1_massJet1->Write();
  h1_massJet2->Write();

  h2_etaPhi_map->Write();
  h2_etaPhi_map_cutOnH->Write();

  h1_deltaRll->Write();
  h1_deltaRll_Nm1->Write();
  h1_deltaRjj->Write();
  h1_deltaRjj_Nm1->Write();
  h1_ptHiggs->Write();
  h1_pzHiggs->Write();
  h1_etaHiggs->Write();
  h1_deltaRZZ->Write();
  h1_deltaRZZ_Nm1->Write();
  h1_deltaEtaZZ->Write();
  h1_deltaEtaAbsZZ->Write();
  h1_deltaPhiZZ->Write();
  h1_deltaPtZZ->Write();

  h1_ptFullSystem->Write();
  h1_ptRecoilOverJet2->Write();
  h1_ptHiggsOverRecoil->Write();
  h1_deltaPhi_ZllRecoil->Write();
  h1_deltaPhi_ZjjRecoil->Write();
  h1_deltaR_ZjjRecoil->Write();
  h1_deltaPhi_HiggsRecoil->Write();

  h1_ptZll->Write();
  h1_pzZll->Write();
  h1_ptZjj->Write();
  h1_pzZjj->Write();
  h1_ptHardestZ->Write();

  h1_massZll->Write();
  h1_massZll_Nm1->Write();
  h1_MuMuInvMass->Write();
  h1_EleEleInvMass->Write();

  h1_massZjj->Write();
  h1_massZjj_cutOnH->Write();
  h1_massZjj_Nm1->Write();
  h1_massZjj_Rch2_050->Write();
  h1_massZjj_Rch2_5070->Write();
  h1_massZjj_Rch2_70100->Write();
  h1_massZjj_RchHIHI->Write();
  h1_massZjj_RchHILO->Write();
  h1_massZjj_RchLOLO->Write();

  h1_pfMet->Write();
  h1_pfMet_minusHiggs->Write();

  h1_ZZInvMass_loMass->Write();
  h1_ZZInvMass_medMass->Write();
  h1_ZZInvMass_medMass_FINEBINNING->Write();
  h1_ZZInvMass_hiMass->Write();

  h1_ZZInvMass_loMass_ZjjTag->Write();
  h1_ZZInvMass_medMass->Write();
  h1_ZZInvMass_hiMass_ZjjTag->Write();
  h1_ZZInvMass_hiMass_fullSelection_tightOLD->Write();
  h1_ZZInvMass_hiMass_fullSelection_tight->Write();
  h1_ZZInvMass_hiMass_fullSelection_medium->Write();
  h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr->Write();
  h1_ZZInvMass_hiMass_fullSelection_loose->Write();

  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass->Write();
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_lept20->Write();
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetLead->Write();
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets->Write();
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets_Rch->Write();
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetPt->Write();
  h1_ZZInvMass_medMass_fullSelection_nokin->Write();
  h1_ZZInvMass_medMass_fullSelection->Write();
  h1_ZZInvMass_medMass_fullSelection_FINEBINNING->Write();
  h1_ZZInvMass_medMass_fullSelection_tight->Write();

  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag->Write();
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag->Write();

  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40->Write();
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40->Write();

  hp_ptJetGenMean->Write();

  writeResponseHistos( outFile, h1_response_vs_pt, "response" );
  writeResponseHistos( outFile, h1_response_vs_pt_Rch050, "response_Rch050" );
  writeResponseHistos( outFile, h1_response_vs_pt_Rch5070, "response_Rch5070" );
  writeResponseHistos( outFile, h1_response_vs_pt_Rch70100, "response_Rch70100" );


  outFile->Close();

  delete h1_ptLept1;
  h1_ptLept1 = 0;
  delete h1_ptLept2;
  h1_ptLept2 = 0;
  delete h1_ptLept2OverLept1;
  h1_ptLept2OverLept1 = 0;

  delete h1_ptJetLead;
  h1_ptJetLead = 0;
  delete h1_ptJetLead2;
  h1_ptJetLead2 = 0;
  delete h1_ptJetLead3;
  h1_ptJetLead3 = 0;
  delete h1_ptRecoilOverJet2;
  h1_ptRecoilOverJet2 = 0;
  delete h1_ptHiggsOverRecoil;
  h1_ptHiggsOverRecoil = 0;
  delete h1_deltaPhi_HiggsRecoil;
  h1_deltaPhi_HiggsRecoil = 0;
  delete h1_deltaPhi_ZllRecoil;
  h1_deltaPhi_ZllRecoil = 0;
  delete h1_deltaPhi_ZjjRecoil;
  h1_deltaPhi_ZjjRecoil = 0;
  delete h1_deltaR_ZjjRecoil;
  h1_deltaR_ZjjRecoil = 0;
  delete h1_ptFullSystem;
  h1_ptFullSystem = 0;
  delete h1_ptJetRecoil;
  h1_ptJetRecoil = 0;

  delete h1_iJet1;
  h1_iJet1 = 0;
  delete h1_iJet2;
  h1_iJet2 = 0;
  delete h1_iJet2MinusiJet1;
  h1_iJet2MinusiJet1 = 0;
  delete h1_iJet2PlusiJet1;
  h1_iJet2PlusiJet1 = 0;
  delete h1_ptJet1;
  h1_ptJet1 = 0;
  delete h1_ptJet2;
  h1_ptJet2 = 0;
  delete h1_ptJet2Rel;
  h1_ptJet2Rel = 0;
  delete h1_ptJet2OverLead;
  h1_ptJet2OverLead = 0;
  delete h1_ptJet2OverJet1;
  h1_ptJet2OverJet1 = 0;
  delete h1_etaJet1;
  h1_etaJet1 = 0;
  delete h1_etaJet2;
  h1_etaJet2 = 0;
  delete h1_RchJet1;
  h1_RchJet1 = 0;
  delete h1_RchJet2;
  h1_RchJet2 = 0;
  delete h1_massJet1;
  h1_massJet1 = 0;
  delete h1_massJet2;
  h1_massJet2 = 0;

  delete h2_etaPhi_map;
  h2_etaPhi_map = 0;
  delete h2_etaPhi_map_cutOnH;
  h2_etaPhi_map_cutOnH = 0;

  delete h1_deltaRll;
  h1_deltaRll = 0;
  delete h1_deltaRll_Nm1;
  h1_deltaRll_Nm1 = 0;
  delete h1_deltaRjj;
  h1_deltaRjj = 0;
  delete h1_deltaRjj_Nm1;
  h1_deltaRjj_Nm1 = 0;
  delete h1_ptHiggs;
  h1_ptHiggs = 0;
  delete h1_pzHiggs;
  h1_pzHiggs = 0;
  delete h1_etaHiggs;
  h1_etaHiggs = 0;
  delete h1_deltaRZZ;
  h1_deltaRZZ = 0;
  delete h1_deltaRZZ_Nm1;
  h1_deltaRZZ_Nm1 = 0;
  delete h1_deltaEtaZZ;
  h1_deltaEtaZZ = 0;
  delete h1_deltaEtaAbsZZ;
  h1_deltaEtaAbsZZ = 0;
  delete h1_deltaPhiZZ;
  h1_deltaPhiZZ = 0;
  delete h1_deltaPtZZ;
  h1_deltaPtZZ = 0;

  delete h1_pfMet;
  h1_pfMet = 0;
  delete h1_pfMet_minusHiggs;
  h1_pfMet_minusHiggs = 0;

  delete h1_totalLumi;
  h1_totalLumi = 0;
  delete h1_run;
  h1_run = 0;

  delete h1_ptZjj;
  h1_ptZjj = 0;
  delete h1_pzZjj;
  h1_pzZjj = 0;
  delete h1_ptZll;
  h1_ptZll = 0;
  delete h1_pzZll;
  h1_pzZll = 0;

  delete h1_massZll;
  h1_massZll = 0;
  delete h1_massZll_Nm1;
  h1_massZll_Nm1 = 0;
  delete h1_MuMuInvMass;
  h1_MuMuInvMass = 0;
  delete h1_EleEleInvMass;
  h1_EleEleInvMass = 0;
  delete h1_massZjj;
  h1_massZjj = 0;
  delete h1_massZjj_cutOnH;
  h1_massZjj_cutOnH = 0;
  delete h1_massZjj_Nm1;
  h1_massZjj_Nm1 = 0;
  delete h1_massZjj_Rch2_050;
  h1_massZjj_Rch2_050 = 0;
  delete h1_massZjj_Rch2_5070;
  h1_massZjj_Rch2_5070 = 0;
  delete h1_massZjj_Rch2_70100;
  h1_massZjj_Rch2_70100 = 0;
  delete h1_massZjj_RchHIHI;
  h1_massZjj_RchHIHI = 0;
  delete h1_massZjj_RchHILO;
  h1_massZjj_RchHILO = 0;
  delete h1_massZjj_RchLOLO;
  h1_massZjj_RchLOLO = 0;
  delete h1_ptHardestZ;
  h1_ptHardestZ = 0;

  delete h1_ZZInvMass_loMass;
  h1_ZZInvMass_loMass = 0;
  delete h1_ZZInvMass_medMass;
  h1_ZZInvMass_medMass = 0;
  delete h1_ZZInvMass_medMass_FINEBINNING;
  h1_ZZInvMass_medMass_FINEBINNING = 0;
  delete h1_ZZInvMass_hiMass;
  h1_ZZInvMass_hiMass = 0;

  delete h1_ZZInvMass_loMass_ZjjTag;
  h1_ZZInvMass_loMass_ZjjTag = 0;
  delete h1_ZZInvMass_hiMass_ZjjTag;
  h1_ZZInvMass_hiMass_ZjjTag = 0;
  delete h1_ZZInvMass_hiMass_fullSelection_tightOLD;
  h1_ZZInvMass_hiMass_fullSelection_tightOLD = 0;
  delete h1_ZZInvMass_hiMass_fullSelection_tight;
  h1_ZZInvMass_hiMass_fullSelection_tight = 0;
  delete h1_ZZInvMass_hiMass_fullSelection_medium;
  h1_ZZInvMass_hiMass_fullSelection_medium = 0;
  delete h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr;
  h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr = 0;
  delete h1_ZZInvMass_hiMass_fullSelection_loose;
  h1_ZZInvMass_hiMass_fullSelection_loose = 0;

  delete h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass;
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass = 0;
  delete h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_lept20;
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_lept20 = 0;
  delete h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetLead;
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetLead = 0;
  delete h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets;
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets = 0;
  delete h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets_Rch;
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_etaJets_Rch = 0;
  delete h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetPt;
  h1_ZZInvMass_medMass_fullSelection_nokin_lowInvMass_jetPt = 0;
  delete h1_ZZInvMass_medMass_fullSelection_nokin;
  h1_ZZInvMass_medMass_fullSelection_nokin = 0;
  delete h1_ZZInvMass_medMass_fullSelection;
  h1_ZZInvMass_medMass_fullSelection = 0;
  delete h1_ZZInvMass_medMass_fullSelection_FINEBINNING;
  h1_ZZInvMass_medMass_fullSelection_FINEBINNING = 0;
  delete h1_ZZInvMass_medMass_fullSelection_tight;
  h1_ZZInvMass_medMass_fullSelection_tight = 0;

  delete h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag;
  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag = 0;
  delete h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag;
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag = 0;

  delete h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40;
  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40 = 0;
  delete h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40;
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40 = 0;
  delete tree;
  tree = 0;


  delete hp_ptJetGenMean;
  hp_ptJetGenMean = 0;

  int nHistos = h1_response_vs_pt.size();
  for( unsigned iHisto=0; iHisto<nHistos; ++iHisto ) {
    delete h1_response_vs_pt[nHistos-iHisto-1];
    h1_response_vs_pt[nHistos-iHisto-1]=0;
    delete h1_response_vs_pt_Rch050[nHistos-iHisto-1];
    h1_response_vs_pt_Rch050[nHistos-iHisto-1]=0;
    delete h1_response_vs_pt_Rch5070[nHistos-iHisto-1];
    h1_response_vs_pt_Rch5070[nHistos-iHisto-1]=0;
    delete h1_response_vs_pt_Rch70100[nHistos-iHisto-1];
    h1_response_vs_pt_Rch70100[nHistos-iHisto-1]=0;
  }


  totalLumi = 0.;

}


void addFile( const std::string& dataset ) {

  std::string infileName = "HZZlljj_2ndLevelTreeW_" + dataset + ".root"; //the W is important: means that files have passed treatment (merging and weights)
  std::string treeName = infileName +"/reducedTree";
  tree->Add(treeName.c_str());
  std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
  TFile* infile = TFile::Open(infileName.c_str(), "READ");
  TH1F* h1_lumi = (TH1F*)infile->Get("lumi");
  if( h1_lumi!=0 ) {
    totalLumi += h1_lumi->GetBinContent(1);
    std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
  } else {
    std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
  }
  infile->Close();


} //addFile


std::vector<TH1F*> getResponseHistos(const std::string& name, int binArraySize, Double_t* ptBins) {

  std::vector<TH1F*> returnVector;

  for( unsigned i=0; i<(binArraySize-1); ++i ) {
    char histoName[100];
    sprintf( histoName, "%s_ptBin_%.0f_%.0f", name.c_str(), ptBins[i], ptBins[i+1]);
    int nbins = 50;
    float xmin = 0.5;
    float xmax = 2.;
    TH1F* newHisto = new TH1F(histoName, "", nbins, xmin, xmax);
    newHisto->Sumw2();
    returnVector.push_back(newHisto);
  }

  return returnVector;

}


void writeResponseHistos( TFile* file, std::vector<TH1F*> h1_response, std::string dirName ) {

  file->mkdir( dirName.c_str() );
  file->cd( dirName.c_str() );

  for( unsigned iHisto=0; iHisto<h1_response.size(); ++iHisto ) h1_response[iHisto]->Write();

  file->cd();

}

Double_t ErrEt(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 5.6;
    b = 1.25;
    c = 0.033;
  }
  else{
    a = 4.8;
    b = 0.89;
    c = 0.043;
  }
  InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
  return InvPerr2;
}

Double_t ErrEta(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 1.215;
    b = 0.037;
    c = 7.941 * 0.0001;
  }
  else{
    a = 1.773;
    b = 0.034;
    c = 3.56 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

Double_t ErrPhi(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 6.65;
    b = 0.04;
    c = 8.49 * 0.00001;
  }
  else{
    a = 2.908;
    b = 0.021;
    c = 2.59 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

