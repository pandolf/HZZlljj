#include "Ntp1Finalizer_WZ.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"

#include "TFitConstraintM.h"
#include "TFitParticleEtEtaPhi.h"
#include "TKinFitter.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "CommonTools/fitTools.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"



float BUMP_MIN = 90.;
float BUMP_MAX = 100.;



class AnalysisJet : public TLorentzVector {

 public:

  AnalysisJet( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {
    rmsCand=0.;
    ptD=0.;
    nCharged=0;
    nNeutral=0;
  }

  AnalysisJet( const TLorentzVector &v) : TLorentzVector( v ) {
    rmsCand=0.;
    ptD=0.;
    nCharged=0;
    nNeutral=0;
  }

  float rmsCand;
  float ptD;
  int nCharged;
  int nNeutral;

  float muonEnergyFraction;
  float electronEnergyFraction;

};





HelicityLikelihoodDiscriminant::HelicityAngles computeHelicityAngles(TLorentzVector leptMinus, TLorentzVector leptPlus, TLorentzVector jet1, TLorentzVector jet2 );



void print(TKinFitter *fitter);
Double_t ErrEt(Float_t Et, Float_t Eta);
Double_t ErrEta(Float_t Et, Float_t Eta);
Double_t ErrPhi(Float_t Et, Float_t Eta);
Double_t ErrEt(Float_t Et, Float_t Eta, int particleType);
Double_t ErrEta(Float_t Et, Float_t Eta, int particleType);
Double_t ErrPhi(Float_t Et, Float_t Eta, int particleType);


std::vector<TH1D*> getHistoVector(int nPtBins, Double_t *ptBins, std::string histoName, int nBins, float xMin, float xMax );


// constructor:

Ntp1Finalizer_WZ::Ntp1Finalizer_WZ( const std::string& dataset, const std::string& selectionType, const std::string& leptType ) : Ntp1Finalizer( "WZ", dataset, leptType ) {

  if( leptType!="ALL" && leptType!="MU" && leptType!="ELE" ) {
    std::cout << "Lept type '" << leptType << "' currently not supported. Exiting." << std::endl;
    exit(9176);
  }

  leptType_ = leptType;

  setSelectionType(selectionType);

  std::string fullFlags = selectionType_ + "_" + leptType_;
  this->set_flags(fullFlags); //this is for the outfile name

}





void Ntp1Finalizer_WZ::finalize() {

  if( outFile_==0 ) this->createOutputFile();


  TH1F* h1_run = new TH1F("run", "", 15149, 132440, 147589);


  TH1D* h1_ptMu1 = new TH1D("ptMu1", "", 250, 0., 250);
  h1_ptMu1->Sumw2();
  TH1D* h1_ptMu2 = new TH1D("ptMu2", "", 250, 0., 250);
  h1_ptMu2->Sumw2();
  TH1D* h1_etaMu = new TH1D("etaMu", "", 100, -2.5, 2.5 );
  h1_etaMu->Sumw2();
  TH1D* h1_phiMu = new TH1D("phiMu", "", 100, -3.1416, 3.1416);
  h1_phiMu->Sumw2();

  TH1D* h1_ptMu1_inBump = new TH1D("ptMu1_inBump", "", 250, 0., 250);
  h1_ptMu1_inBump->Sumw2();
  TH1D* h1_ptMu2_inBump = new TH1D("ptMu2_inBump", "", 250, 0., 250);
  h1_ptMu2_inBump->Sumw2();
  TH1D* h1_etaMu1_inBump = new TH1D("etaMu1_inBump", "", 100, -2.5, 2.5);
  h1_etaMu1_inBump->Sumw2();
  TH1D* h1_etaMu2_inBump = new TH1D("etaMu2_inBump", "", 100, -2.5, 2.5);
  h1_etaMu2_inBump->Sumw2();
  TH1D* h1_phiMu1_inBump = new TH1D("phiMu1_inBump", "", 100, -3.1416, 3.1416);
  h1_phiMu1_inBump->Sumw2();
  TH1D* h1_phiMu2_inBump = new TH1D("phiMu2_inBump", "", 100, -3.1416, 3.1416);
  h1_phiMu2_inBump->Sumw2();

  TH1D* h1_ptMu1_outBump = new TH1D("ptMu1_outBump", "", 250, 0., 250);
  h1_ptMu1_outBump->Sumw2();
  TH1D* h1_ptMu2_outBump = new TH1D("ptMu2_outBump", "", 250, 0., 250);
  h1_ptMu2_outBump->Sumw2();
  TH1D* h1_etaMu1_outBump = new TH1D("etaMu1_outBump", "", 100, -2.5, 2.5);
  h1_etaMu1_outBump->Sumw2();
  TH1D* h1_etaMu2_outBump = new TH1D("etaMu2_outBump", "", 100, -2.5, 2.5);
  h1_etaMu2_outBump->Sumw2();
  TH1D* h1_phiMu1_outBump = new TH1D("phiMu1_outBump", "", 100, -3.1416, 3.1416);
  h1_phiMu1_outBump->Sumw2();
  TH1D* h1_phiMu2_outBump = new TH1D("phiMu2_outBump", "", 100, -3.1416, 3.1416);
  h1_phiMu2_outBump->Sumw2();

  TH1D* h1_deltaPhiMu1Met_inBump = new TH1D("deltaPhiMu1Met_inBump", "", 100, 0, 3.14159);
  h1_deltaPhiMu1Met_inBump->Sumw2();
  TH1D* h1_deltaPhiMu1Met_outBump = new TH1D("deltaPhiMu1Met_outBump", "", 100, 0, 3.14159);
  h1_deltaPhiMu1Met_outBump->Sumw2();

  TH2D* h2_etaMu_vs_ptMu = new TH2D("etaMu_vs_ptMu", "", 100, 0., 250., 100, -2.5, 2.5 );
  h2_etaMu_vs_ptMu->Sumw2();
  TH2D* h2_phiMu_vs_ptMu = new TH2D("phiMu_vs_ptMu", "", 100, 0., 250., 100, -3.1416, 3.1416);
  h2_phiMu_vs_ptMu->Sumw2();
  TH2D* h2_phiMu_vs_etaMu = new TH2D("phiMu_vs_etaMu", "", 100, -2.5, 2.5, 100, -3.1416, 3.1416);
  h2_phiMu_vs_etaMu->Sumw2();

  TH2D* h2_etaMu_vs_ptZll = new TH2D("etaMu_vs_ptZll", "", 100, 0., 250., 100, -2.5, 2.5 );
  h2_etaMu_vs_ptZll->Sumw2();
  TH2D* h2_phiMu_vs_ptZll = new TH2D("phiMu_vs_ptZll", "", 100, 0., 250., 100, -3.1416, 3.1416);
  h2_phiMu_vs_ptZll->Sumw2();
  TH2D* h2_phiMu_vs_etaZll = new TH2D("phiMu_vs_etaZll", "", 100, -2.5, 2.5, 100, -3.1416, 3.1416);
  h2_phiMu_vs_etaZll->Sumw2();


  TH1D* h1_ptEle1 = new TH1D("ptEle1", "", 250, 0., 250);
  h1_ptEle1->Sumw2();
  TH1D* h1_ptEle2 = new TH1D("ptEle2", "", 250, 0., 250);
  h1_ptEle2->Sumw2();
  TH1D* h1_etaEle = new TH1D("etaEle", "", 100, -2.5, 2.5 );
  h1_etaEle->Sumw2();
  TH1D* h1_phiEle = new TH1D("phiEle", "", 100, -3.1416, 3.1416);
  h1_phiEle->Sumw2();

  TH2D* h2_etaEle_vs_ptEle = new TH2D("etaEle_vs_ptEle", "", 100, 0., 250., 100, -2.5, 2.5 );
  h2_etaEle_vs_ptEle->Sumw2();
  TH2D* h2_phiEle_vs_ptEle = new TH2D("phiEle_vs_ptEle", "", 100, 0., 250., 100, -3.1416, 3.1416);
  h2_phiEle_vs_ptEle->Sumw2();
  TH2D* h2_phiEle_vs_etaEle = new TH2D("phiEle_vs_etaEle", "", 100, -2.5, 2.5, 100, -3.1416, 3.1416);
  h2_phiEle_vs_etaEle->Sumw2();

  TH2D* h2_etaEle_vs_ptZll = new TH2D("etaEle_vs_ptZll", "", 100, 0., 250., 100, -2.5, 2.5 );
  h2_etaEle_vs_ptZll->Sumw2();
  TH2D* h2_phiEle_vs_ptZll = new TH2D("phiEle_vs_ptZll", "", 100, 0., 250., 100, -3.1416, 3.1416);
  h2_phiEle_vs_ptZll->Sumw2();
  TH2D* h2_phiEle_vs_etaZll = new TH2D("phiEle_vs_etaZll", "", 100, -2.5, 2.5, 100, -3.1416, 3.1416);
  h2_phiEle_vs_etaZll->Sumw2();


  TH1D* h1_ptZll= new TH1D("ptZll", "", 250, 0., 250.);
  h1_ptZll->Sumw2();
  TH1D* h1_etaZll= new TH1D("etaZll", "", 15, -1.3, 1.3);
  h1_etaZll->Sumw2();
  TH1D* h1_phiZll= new TH1D("phiZll", "", 15, -3.1416, 3.1416);
  h1_phiZll->Sumw2();
  TH1D* h1_ptZmumu= new TH1D("ptZmumu", "", 250, 0., 250.);
  h1_ptZmumu->Sumw2();
  TH1D* h1_etaZmumu= new TH1D("etaZmumu", "", 15, -1.3, 1.3);
  h1_etaZmumu->Sumw2();
  TH1D* h1_phiZmumu= new TH1D("phiZmumu", "", 15, -3.1416, 3.1416);
  h1_phiZmumu->Sumw2();
  TH1D* h1_ptZee= new TH1D("ptZee", "", 250, 0., 250.);
  h1_ptZee->Sumw2();
  TH1D* h1_etaZee= new TH1D("etaZee", "", 15, -1.3, 1.3);
  h1_etaZee->Sumw2();
  TH1D* h1_phiZee= new TH1D("phiZee", "", 15, -3.1416, 3.1416);
  h1_phiZee->Sumw2();

  TH1D* h1_nJets= new TH1D("nJets", "", 7, -0.5, 6.5);
  h1_nJets->Sumw2();
  TH1D* h1_nJets_inBump = new TH1D("nJets_inBump", "", 7, -0.5, 6.5);
  h1_nJets_inBump->Sumw2();
  TH1D* h1_nJets_outBump = new TH1D("nJets_outBump ", "", 7, -0.5, 6.5);
  h1_nJets_outBump ->Sumw2();

  TH1D* h1_ptZll_0jet = new TH1D("ptZll_0jet", "", 250, 0., 250.);
  h1_ptZll_0jet->Sumw2();
  TH1D* h1_ptZll_1jet = new TH1D("ptZll_1jet", "", 250, 0., 250.);
  h1_ptZll_1jet->Sumw2();
  TH1D* h1_ptZll_2jet = new TH1D("ptZll_2jet", "", 250, 0., 250.);
  h1_ptZll_2jet->Sumw2();
  TH1D* h1_ptZll_3jet = new TH1D("ptZll_3jet", "", 250, 0., 250.);
  h1_ptZll_3jet->Sumw2();
  TH1D* h1_ptZll_4jet = new TH1D("ptZll_4jet", "", 250, 0., 250.);
  h1_ptZll_4jet->Sumw2();
  TH1D* h1_ptZmumu_0jet = new TH1D("ptZmumu_0jet", "", 250, 0., 250.);
  h1_ptZmumu_0jet->Sumw2();
  TH1D* h1_ptZmumu_1jet = new TH1D("ptZmumu_1jet", "", 250, 0., 250.);
  h1_ptZmumu_1jet->Sumw2();
  TH1D* h1_ptZmumu_2jet = new TH1D("ptZmumu_2jet", "", 250, 0., 250.);
  h1_ptZmumu_2jet->Sumw2();
  TH1D* h1_ptZmumu_3jet = new TH1D("ptZmumu_3jet", "", 250, 0., 250.);
  h1_ptZmumu_3jet->Sumw2();
  TH1D* h1_ptZmumu_4jet = new TH1D("ptZmumu_4jet", "", 250, 0., 250.);
  h1_ptZmumu_4jet->Sumw2();
  TH1D* h1_ptZee_0jet = new TH1D("ptZee_0jet", "", 250, 0., 250.);
  h1_ptZee_0jet->Sumw2();
  TH1D* h1_ptZee_1jet = new TH1D("ptZee_1jet", "", 250, 0., 250.);
  h1_ptZee_1jet->Sumw2();
  TH1D* h1_ptZee_2jet = new TH1D("ptZee_2jet", "", 250, 0., 250.);
  h1_ptZee_2jet->Sumw2();
  TH1D* h1_ptZee_3jet = new TH1D("ptZee_3jet", "", 250, 0., 250.);
  h1_ptZee_3jet->Sumw2();
  TH1D* h1_ptZee_4jet = new TH1D("ptZee_4jet", "", 250, 0., 250.);
  h1_ptZee_4jet->Sumw2();

  TH1D* h1_eleEnergyFractionJet_inBump = new TH1D("eleEnergyFractionJet_inBump", "", 100, 0., 1.0001);
  h1_eleEnergyFractionJet_inBump->Sumw2();
  TH1D* h1_eleEnergyFractionJet_outBump = new TH1D("eleEnergyFractionJet_outBump", "", 100, 0., 1.0001);
  h1_eleEnergyFractionJet_outBump->Sumw2();
  TH1D* h1_muonEnergyFractionJet_inBump = new TH1D("muonEnergyFractionJet_inBump", "", 100, 0., 1.0001);
  h1_muonEnergyFractionJet_inBump->Sumw2();
  TH1D* h1_muonEnergyFractionJet_outBump = new TH1D("muonEnergyFractionJet_outBump", "", 100, 0., 1.0001);
  h1_muonEnergyFractionJet_outBump->Sumw2();

  TH1D* h1_mJet_inBump = new TH1D("mJet_inBump", "", 150, 0., 150.);
  h1_mJet_inBump->Sumw2();
  TH1D* h1_mJet_outBump = new TH1D("mJet_outBump", "", 150, 0., 150.);
  h1_mJet_outBump->Sumw2();
  TH1D* h1_mJet_1jet = new TH1D("mJet_1jet", "", 150, 0., 150.);
  h1_mJet_1jet->Sumw2();
  TH1D* h1_mJet_2jet = new TH1D("mJet_2jet", "", 150, 0., 150.);
  h1_mJet_2jet->Sumw2();
  TH1D* h1_mJet_3jet = new TH1D("mJet_3jet", "", 150, 0., 150.);
  h1_mJet_3jet->Sumw2();
  TH1D* h1_mJet_4jet = new TH1D("mJet_4jet", "", 150, 0., 150.);
  h1_mJet_4jet->Sumw2();

  TH1D* h1_mZll_inBump = new TH1D("mZll_inBump", "", 100, 30., 130.);
  h1_mZll_inBump->Sumw2();
  TH1D* h1_mZmumu_inBump = new TH1D("mZmumu_inBump", "", 100, 30., 130.);
  h1_mZmumu_inBump->Sumw2();
  TH1D* h1_mZee_inBump = new TH1D("mZee_inBump", "", 100, 30., 130.);
  h1_mZee_inBump->Sumw2();
  TH1D* h1_mZll_outBump = new TH1D("mZll_outBump", "", 100, 30., 130.);
  h1_mZll_outBump->Sumw2();
  TH1D* h1_mZmumu_outBump = new TH1D("mZmumu_outBump", "", 100, 30., 130.);
  h1_mZmumu_outBump->Sumw2();
  TH1D* h1_mZee_outBump = new TH1D("mZee_outBump", "", 100, 30., 130.);
  h1_mZee_outBump->Sumw2();

  float mZllJet_min = 50.;
  float mZllJet_max = 750.;

  TH1D* h1_mtZllMet_inBump = new TH1D("mtZllMet_inBump", "", 450, 50., 500);
  h1_mtZllMet_inBump->Sumw2();
  TH1D* h1_mtZllMet_outBump = new TH1D("mtZllMet_outBump", "", 450, 50., 500);
  h1_mtZllMet_outBump->Sumw2();
  TH1D* h1_mtZmumuMet_inBump = new TH1D("mtZmumuMet_inBump", "", 450, 50., 500);
  h1_mtZmumuMet_inBump->Sumw2();
  TH1D* h1_mtZmumuMet_outBump = new TH1D("mtZmumuMet_outBump", "", 450, 50., 500);
  h1_mtZmumuMet_outBump->Sumw2();
  TH1D* h1_mtZeeMet_inBump = new TH1D("mtZeeMet_inBump", "", 450, 50., 500);
  h1_mtZeeMet_inBump->Sumw2();
  TH1D* h1_mtZeeMet_outBump = new TH1D("mtZeeMet_outBump", "", 450, 50., 500);
  h1_mtZeeMet_outBump->Sumw2();

  TH1D* h1_mZllJet_inBump = new TH1D("mZllJet_inBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZllJet_inBump->Sumw2();
  TH1D* h1_mZllJet_outBump = new TH1D("mZllJet_outBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZllJet_outBump->Sumw2();
  TH1D* h1_mZllJet_1jet = new TH1D("mZllJet_1jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZllJet_1jet->Sumw2();
  TH1D* h1_mZllJet_2jet = new TH1D("mZllJet_2jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZllJet_2jet->Sumw2();
  TH1D* h1_mZllJet_3jet = new TH1D("mZllJet_3jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZllJet_3jet->Sumw2();
  TH1D* h1_mZllJet_4jet = new TH1D("mZllJet_4jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZllJet_4jet->Sumw2();

  TH1D* h1_mZmumuJet_inBump = new TH1D("mZmumuJet_inBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZmumuJet_inBump->Sumw2();
  TH1D* h1_mZmumuJet_outBump = new TH1D("mZmumuJet_outBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZmumuJet_outBump->Sumw2();
  TH1D* h1_mZmumuJet_1jet = new TH1D("mZmumuJet_1jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZmumuJet_1jet->Sumw2();
  TH1D* h1_mZmumuJet_2jet = new TH1D("mZmumuJet_2jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZmumuJet_2jet->Sumw2();
  TH1D* h1_mZmumuJet_3jet = new TH1D("mZmumuJet_3jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZmumuJet_3jet->Sumw2();
  TH1D* h1_mZmumuJet_4jet = new TH1D("mZmumuJet_4jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZmumuJet_4jet->Sumw2();

  TH1D* h1_mZeeJet_inBump = new TH1D("mZeeJet_inBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZeeJet_inBump->Sumw2();
  TH1D* h1_mZeeJet_outBump = new TH1D("mZeeJet_outBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZeeJet_outBump->Sumw2();
  TH1D* h1_mZeeJet_1jet = new TH1D("mZeeJet_1jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZeeJet_1jet->Sumw2();
  TH1D* h1_mZeeJet_2jet = new TH1D("mZeeJet_2jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZeeJet_2jet->Sumw2();
  TH1D* h1_mZeeJet_3jet = new TH1D("mZeeJet_3jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZeeJet_3jet->Sumw2();
  TH1D* h1_mZeeJet_4jet = new TH1D("mZeeJet_4jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mZeeJet_4jet->Sumw2();

  TH1D* h1_mtZllJet_inBump = new TH1D("mtZllJet_inBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZllJet_inBump->Sumw2();
  TH1D* h1_mtZllJet_outBump = new TH1D("mtZllJet_outBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZllJet_outBump->Sumw2();
  TH1D* h1_mtZllJet_1jet = new TH1D("mtZllJet_1jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZllJet_1jet->Sumw2();
  TH1D* h1_mtZllJet_2jet = new TH1D("mtZllJet_2jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZllJet_2jet->Sumw2();
  TH1D* h1_mtZllJet_3jet = new TH1D("mtZllJet_3jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZllJet_3jet->Sumw2();
  TH1D* h1_mtZllJet_4jet = new TH1D("mtZllJet_4jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZllJet_4jet->Sumw2();

  TH1D* h1_mtZmumuJet_inBump = new TH1D("mtZmumuJet_inBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZmumuJet_inBump->Sumw2();
  TH1D* h1_mtZmumuJet_outBump = new TH1D("mtZmumuJet_outBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZmumuJet_outBump->Sumw2();
  TH1D* h1_mtZmumuJet_1jet = new TH1D("mtZmumuJet_1jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZmumuJet_1jet->Sumw2();
  TH1D* h1_mtZmumuJet_2jet = new TH1D("mtZmumuJet_2jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZmumuJet_2jet->Sumw2();
  TH1D* h1_mtZmumuJet_3jet = new TH1D("mtZmumuJet_3jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZmumuJet_3jet->Sumw2();
  TH1D* h1_mtZmumuJet_4jet = new TH1D("mtZmumuJet_4jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZmumuJet_4jet->Sumw2();

  TH1D* h1_mtZeeJet_inBump = new TH1D("mtZeeJet_inBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZeeJet_inBump->Sumw2();
  TH1D* h1_mtZeeJet_outBump = new TH1D("mtZeeJet_outBump", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZeeJet_outBump->Sumw2();
  TH1D* h1_mtZeeJet_1jet = new TH1D("mtZeeJet_1jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZeeJet_1jet->Sumw2();
  TH1D* h1_mtZeeJet_2jet = new TH1D("mtZeeJet_2jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZeeJet_2jet->Sumw2();
  TH1D* h1_mtZeeJet_3jet = new TH1D("mtZeeJet_3jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZeeJet_3jet->Sumw2();
  TH1D* h1_mtZeeJet_4jet = new TH1D("mtZeeJet_4jet", "", (int)mZllJet_max-mZllJet_min, mZllJet_min, mZllJet_max);
  h1_mtZeeJet_4jet->Sumw2();

  TH1D* h1_met_inBump = new TH1D("met_inBump", "", 250, 0., 250.);
  h1_met_inBump->Sumw2();
  TH1D* h1_met_outBump = new TH1D("met_outBump", "", 250, 0., 250.);
  h1_met_outBump->Sumw2();
  TH1D* h1_met_1jet = new TH1D("met_1jet", "", 250, 0., 250.);
  h1_met_1jet->Sumw2();
  TH1D* h1_met_2jet = new TH1D("met_2jet", "", 250, 0., 250.);
  h1_met_2jet->Sumw2();
  TH1D* h1_met_3jet = new TH1D("met_3jet", "", 250, 0., 250.);
  h1_met_3jet->Sumw2();
  TH1D* h1_met_4jet = new TH1D("met_4jet", "", 250, 0., 250.);
  h1_met_4jet->Sumw2();

  TH1D* h1_mZllW= new TH1D("mZllW", "", 500, 150., 550.);
  h1_mZllW->Sumw2();
  TH1D* h1_mZllW_kinfit = new TH1D("mZllW_kinfit", "", 500, 150., 550.);
  h1_mZllW_kinfit->Sumw2();
  TH1D* h1_mZmumuW= new TH1D("mZmumuW", "", 500, 150., 550.);
  h1_mZmumuW->Sumw2();
  TH1D* h1_mZmumuW_kinfit = new TH1D("mZmumuW_kinfit", "", 500, 150., 550.);
  h1_mZmumuW_kinfit->Sumw2();
  TH1D* h1_mZeeW= new TH1D("mZeeW", "", 500, 150., 550.);
  h1_mZeeW->Sumw2();
  TH1D* h1_mZeeW_kinfit = new TH1D("mZeeW_kinfit", "", 500, 150., 550.);
  h1_mZeeW_kinfit->Sumw2();



  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Int_t LS;
  tree_->SetBranchAddress("LS", &LS);
  Int_t event;
  tree_->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);

  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Float_t pfMet;
  tree_->SetBranchAddress("epfMet", &pfMet);
  Float_t phiMet;
  tree_->SetBranchAddress("phipfMet", &phiMet);


  int leptType;
  tree_->SetBranchAddress("leptType", &leptType);

  Float_t eLept1;
  tree_->SetBranchAddress("eLept1", &eLept1);
  Float_t ptLept1;
  tree_->SetBranchAddress("ptLept1", &ptLept1);
  Float_t etaLept1;
  tree_->SetBranchAddress("etaLept1", &etaLept1);
  Float_t phiLept1;
  tree_->SetBranchAddress("phiLept1", &phiLept1);
  Int_t chargeLept1;
  tree_->SetBranchAddress("chargeLept1", &chargeLept1);

  Float_t eLept2;
  tree_->SetBranchAddress("eLept2", &eLept2);
  Float_t ptLept2;
  tree_->SetBranchAddress("ptLept2", &ptLept2);
  Float_t etaLept2;
  tree_->SetBranchAddress("etaLept2", &etaLept2);
  Float_t phiLept2;
  tree_->SetBranchAddress("phiLept2", &phiLept2);
  Int_t chargeLept2;
  tree_->SetBranchAddress("chargeLept2", &chargeLept2);


  Int_t nJet;
  tree_->SetBranchAddress("nJet", &nJet);

  Int_t iJet[50];
  tree_->SetBranchAddress("iJet", iJet);
  Float_t eJet[50];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[50];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t etaJet[50];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[50];
  tree_->SetBranchAddress("phiJet", phiJet);
  Float_t eChargedHadronsJet[50];
  tree_->SetBranchAddress("eChargedHadronsJet", eChargedHadronsJet);
  Float_t rmsCandJet[50];
  tree_->SetBranchAddress("rmsCandJet", rmsCandJet);
  Float_t ptDJet[50];
  tree_->SetBranchAddress("ptDJet", ptDJet);
  Int_t nChargedJet[50];
  tree_->SetBranchAddress("nChargedJet", nChargedJet);
  Int_t nNeutralJet[50];
  tree_->SetBranchAddress("nNeutralJet", nNeutralJet);
  Float_t eMuonsJet[50];
  tree_->SetBranchAddress("eMuonsJet", eMuonsJet);
  Float_t eElectronsJet[50];
  tree_->SetBranchAddress("eElectronsJet", eElectronsJet);

  Int_t nPart;
  tree_->SetBranchAddress("nPart", &nPart);
  Float_t ePart[20];
  tree_->SetBranchAddress("ePart", ePart);
  Float_t ptPart[20];
  tree_->SetBranchAddress("ptPart", ptPart);
  Float_t etaPart[20];
  tree_->SetBranchAddress("etaPart", etaPart);
  Float_t phiPart[20];
  tree_->SetBranchAddress("phiPart", phiPart);
  Int_t pdgIdPart[20];
  tree_->SetBranchAddress("pdgIdPart", pdgIdPart);



  int nEntries = tree_->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_3_8_7/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root");


  ofstream ofs_events("events.txt");


  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    if( leptType_!="ALL" ) {
      if( leptType_=="ELE" && leptType==0 ) continue;
      if( leptType_=="MU" && leptType==1 ) continue;
    }



    if( run>5 ) { //is not MC:
    

      // remove duplicate events:

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

    TLorentzVector diLepton = lept1+lept2;

    if( lept1.Pt()<ptLept1_thresh_ ) continue;
    if( lept2.Pt()<ptLept2_thresh_ ) continue;
    if( fabs(lept1.Eta())>etaLept1_thresh_ ) continue;
    if( fabs(lept2.Eta())>etaLept2_thresh_ ) continue;

    if( fabs(diLepton.Eta())>etaZll_thresh_ ) continue;

    // fill this before mZ cut:
    if( diLepton.Pt()<BUMP_MAX && diLepton.Pt()>BUMP_MIN ) {
      h1_mZll_inBump->Fill( diLepton.M(), eventWeight );
      if( leptType==0 ) h1_mZmumu_inBump->Fill( diLepton.M(), eventWeight );
      else h1_mZee_inBump->Fill( diLepton.M(), eventWeight );
    } else {
      h1_mZll_outBump->Fill( diLepton.M(), eventWeight );
      if( leptType==0 ) h1_mZmumu_outBump->Fill( diLepton.M(), eventWeight );
      else h1_mZee_outBump->Fill( diLepton.M(), eventWeight );
    }

    if( diLepton.M() < mZll_threshLo_ || diLepton.M() > mZll_threshHi_ ) continue;



    // passed selection

    std::vector<AnalysisJet> selectedJets;

    for( unsigned iJet=0; iJet<nJet; ++iJet) {

      AnalysisJet jet;
      jet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      jet.rmsCand = rmsCandJet[iJet];
      jet.ptD = ptDJet[iJet];
      jet.nCharged = nChargedJet[iJet];
      jet.nNeutral = nNeutralJet[iJet];
      jet.muonEnergyFraction = eMuonsJet[iJet]/jet.Energy();
      jet.electronEnergyFraction = eElectronsJet[iJet]/jet.Energy();


      if( jet.Pt()>ptJet_thresh_ && fabs(jet.Eta())<etaJet_thresh_ )
        selectedJets.push_back( jet );

    }


    TLorentzVector v_met;
    v_met.SetPtEtaPhiE( pfMet, 0., phiMet, pfMet );
    TLorentzVector Zmet = diLepton+v_met;

    if( diLepton.Pt()<BUMP_MAX && diLepton.Pt()>BUMP_MIN ) {

ofs_events << run << ":" << LS << ":" << event << std::endl;

      if( leptType==0 ) {
        h1_ptMu1_inBump->Fill( lept1.Pt(), eventWeight );
        h1_ptMu2_inBump->Fill( lept2.Pt(), eventWeight );
        h1_etaMu1_inBump->Fill( lept1.Eta(), eventWeight );
        h1_etaMu2_inBump->Fill( lept2.Eta(), eventWeight );
        h1_phiMu1_inBump->Fill( lept1.Phi(), eventWeight );
        h1_phiMu2_inBump->Fill( lept2.Phi(), eventWeight );
        h1_mtZmumuMet_inBump->Fill( Zmet.Mt(), eventWeight );
        h1_deltaPhiMu1Met_inBump->Fill( fabs(v_met.DeltaPhi(lept1)), eventWeight );
      } else {
        h1_mtZeeMet_inBump->Fill( Zmet.Mt(), eventWeight );
      }

      h1_nJets_inBump->Fill( selectedJets.size(), eventWeight );
      h1_met_inBump->Fill( pfMet, eventWeight );
      h1_mtZllMet_inBump->Fill( Zmet.Mt(), eventWeight );

      for( unsigned iJet=0; iJet<selectedJets.size(); ++iJet ) {
        AnalysisJet thisJet = selectedJets[iJet];
        h1_mJet_inBump->Fill( thisJet.M(), eventWeight );
        h1_eleEnergyFractionJet_inBump->Fill( thisJet.electronEnergyFraction, eventWeight );
        h1_muonEnergyFractionJet_inBump->Fill( thisJet.muonEnergyFraction, eventWeight );
        TLorentzVector ZJet = diLepton+thisJet;
        h1_mZllJet_inBump->Fill( ZJet.M(), eventWeight );
        h1_mtZllJet_inBump->Fill( ZJet.Mt(), eventWeight );
        if( leptType==0 ) {
          h1_mZmumuJet_inBump->Fill( diLepton.Pt(), eventWeight );
          h1_mtZmumuJet_inBump->Fill( diLepton.Pt(), eventWeight );
        } else {
          h1_mZeeJet_inBump->Fill( diLepton.Pt(), eventWeight );
          h1_mtZeeJet_inBump->Fill( diLepton.Pt(), eventWeight );
        }
      }

    } else {

      if( leptType==0 ) {
        h1_ptMu1_outBump->Fill( lept1.Pt(), eventWeight );
        h1_ptMu2_outBump->Fill( lept2.Pt(), eventWeight );
        h1_etaMu1_outBump->Fill( lept1.Eta(), eventWeight );
        h1_etaMu2_outBump->Fill( lept2.Eta(), eventWeight );
        h1_phiMu1_outBump->Fill( lept1.Phi(), eventWeight );
        h1_phiMu2_outBump->Fill( lept2.Phi(), eventWeight );
        h1_mtZmumuMet_outBump->Fill( Zmet.Mt(), eventWeight );
        h1_deltaPhiMu1Met_outBump->Fill( fabs(v_met.DeltaPhi(lept1)), eventWeight );
      } else {
        h1_mtZeeMet_outBump->Fill( Zmet.Mt(), eventWeight );
      }

      h1_nJets_outBump->Fill( selectedJets.size(), eventWeight );
      h1_met_outBump->Fill( pfMet, eventWeight );
      h1_mtZllMet_outBump->Fill( Zmet.Mt(), eventWeight );
      if( leptType==0 ) h1_mtZmumuMet_outBump->Fill( Zmet.Mt(), eventWeight );
      else h1_mtZeeMet_outBump->Fill( Zmet.Mt(), eventWeight );

      for( unsigned iJet=0; iJet<selectedJets.size(); ++iJet ) {
        AnalysisJet thisJet = selectedJets[iJet];
        h1_mJet_outBump->Fill( thisJet.M(), eventWeight );
        h1_eleEnergyFractionJet_outBump->Fill( thisJet.electronEnergyFraction, eventWeight );
        h1_muonEnergyFractionJet_outBump->Fill( thisJet.muonEnergyFraction, eventWeight );
        TLorentzVector ZJet = diLepton+thisJet;
        h1_mZllJet_outBump->Fill( ZJet.M(), eventWeight );
        h1_mtZllJet_outBump->Fill( ZJet.Mt(), eventWeight );
        if( leptType==0 ) {
          h1_mZmumuJet_outBump->Fill( diLepton.Pt(), eventWeight );
          h1_mtZmumuJet_outBump->Fill( diLepton.Pt(), eventWeight );
        } else {
          h1_mZeeJet_outBump->Fill( diLepton.Pt(), eventWeight );
          h1_mtZeeJet_outBump->Fill( diLepton.Pt(), eventWeight );
        }
      }
    }




    h1_nJets->Fill( selectedJets.size(), eventWeight );

    h1_ptZll->Fill( diLepton.Pt(), eventWeight );
    h1_etaZll->Fill( diLepton.Eta(), eventWeight );
    h1_phiZll->Fill( diLepton.Phi(), eventWeight );

    if( leptType==0 ) {

      h1_ptZmumu->Fill( diLepton.Pt(), eventWeight );
      h1_etaZmumu->Fill( diLepton.Eta(), eventWeight );
      h1_phiZmumu->Fill( diLepton.Phi(), eventWeight );

      h1_ptMu1->Fill( lept1.Pt(), eventWeight );
      h1_ptMu2->Fill( lept2.Pt(), eventWeight );
      h1_etaMu->Fill( lept1.Eta(), eventWeight );
      h1_etaMu->Fill( lept2.Eta(), eventWeight );
      h1_phiMu->Fill( lept1.Phi(), eventWeight );
      h1_phiMu->Fill( lept2.Phi(), eventWeight );

      h2_etaMu_vs_ptMu->Fill( lept1.Pt(), lept1.Eta(), eventWeight );
      h2_etaMu_vs_ptMu->Fill( lept2.Pt(), lept2.Eta(), eventWeight );
      h2_phiMu_vs_ptMu->Fill( lept1.Pt(), lept1.Phi(), eventWeight );
      h2_phiMu_vs_ptMu->Fill( lept2.Pt(), lept2.Phi(), eventWeight );
      h2_phiMu_vs_etaMu->Fill( lept1.Eta(), lept1.Phi(), eventWeight );
      h2_phiMu_vs_etaMu->Fill( lept2.Eta(), lept2.Phi(), eventWeight );

      h2_etaMu_vs_ptZll->Fill( diLepton.Pt(), lept1.Eta(), eventWeight );
      h2_etaMu_vs_ptZll->Fill( diLepton.Pt(), lept2.Eta(), eventWeight );
      h2_phiMu_vs_ptZll->Fill( diLepton.Pt(), lept1.Phi(), eventWeight );
      h2_phiMu_vs_ptZll->Fill( diLepton.Pt(), lept2.Phi(), eventWeight );
      h2_phiMu_vs_etaZll->Fill( diLepton.Eta(), lept1.Phi(), eventWeight );
      h2_phiMu_vs_etaZll->Fill( diLepton.Eta(), lept2.Phi(), eventWeight );

    } else {

      h1_ptZee->Fill( diLepton.Pt(), eventWeight );
      h1_etaZee->Fill( diLepton.Eta(), eventWeight );
      h1_phiZee->Fill( diLepton.Phi(), eventWeight );

      h1_ptEle1->Fill( lept1.Eta(), eventWeight );
      h1_ptEle2->Fill( lept2.Eta(), eventWeight );
      h1_etaEle->Fill( lept1.Eta(), eventWeight );
      h1_etaEle->Fill( lept2.Eta(), eventWeight );
      h1_phiEle->Fill( lept1.Eta(), eventWeight );
      h1_phiEle->Fill( lept2.Eta(), eventWeight );

      h2_etaEle_vs_ptEle->Fill( lept1.Pt(), lept1.Eta(), eventWeight );
      h2_etaEle_vs_ptEle->Fill( lept2.Pt(), lept2.Eta(), eventWeight );
      h2_phiEle_vs_ptEle->Fill( lept1.Pt(), lept1.Phi(), eventWeight );
      h2_phiEle_vs_ptEle->Fill( lept2.Pt(), lept2.Phi(), eventWeight );
      h2_phiEle_vs_etaEle->Fill( lept1.Eta(), lept1.Phi(), eventWeight );
      h2_phiEle_vs_etaEle->Fill( lept2.Eta(), lept2.Phi(), eventWeight );

      h2_etaEle_vs_ptZll->Fill( diLepton.Pt(), lept1.Eta(), eventWeight );
      h2_etaEle_vs_ptZll->Fill( diLepton.Pt(), lept2.Eta(), eventWeight );
      h2_phiEle_vs_ptZll->Fill( diLepton.Pt(), lept1.Phi(), eventWeight );
      h2_phiEle_vs_ptZll->Fill( diLepton.Pt(), lept2.Phi(), eventWeight );
      h2_phiEle_vs_etaZll->Fill( diLepton.Eta(), lept1.Phi(), eventWeight );
      h2_phiEle_vs_etaZll->Fill( diLepton.Eta(), lept2.Phi(), eventWeight );

    }

    
    if( selectedJets.size()==0 ) {

      h1_ptZll_0jet->Fill( diLepton.Pt(), eventWeight );
      if( leptType==0 ) h1_ptZmumu_0jet->Fill( diLepton.Pt(), eventWeight );
      else if( leptType==1 ) h1_ptZee_0jet->Fill( diLepton.Pt(), eventWeight );


    } else if( selectedJets.size()==1 ) {

      h1_ptZll_1jet->Fill( diLepton.Pt(), eventWeight );
      if( leptType==0 ) h1_ptZmumu_1jet->Fill( diLepton.Pt(), eventWeight );
      else if( leptType==1 ) h1_ptZee_1jet->Fill( diLepton.Pt(), eventWeight );

      for( unsigned iJet=0; iJet<selectedJets.size(); ++iJet ) {
        AnalysisJet thisJet = selectedJets[iJet];
        h1_mJet_1jet->Fill( thisJet.M(), eventWeight );
        TLorentzVector ZJet = diLepton+thisJet;
        h1_mZllJet_1jet->Fill( ZJet.M(), eventWeight );
        h1_mtZllJet_1jet->Fill( ZJet.Mt(), eventWeight );
        if( leptType==0 ) {
          h1_mZmumuJet_1jet->Fill( diLepton.Pt(), eventWeight );
          h1_mtZmumuJet_1jet->Fill( diLepton.Pt(), eventWeight );
        } else {
          h1_mZeeJet_1jet->Fill( diLepton.Pt(), eventWeight );
          h1_mtZeeJet_1jet->Fill( diLepton.Pt(), eventWeight );
        }
        h1_met_1jet->Fill( pfMet, eventWeight );
      }

    } else if( selectedJets.size()==2 ) {

      h1_ptZll_2jet->Fill( diLepton.Pt(), eventWeight );
      if( leptType==0 ) h1_ptZmumu_2jet->Fill( diLepton.Pt(), eventWeight );
      else if( leptType==1 ) h1_ptZee_2jet->Fill( diLepton.Pt(), eventWeight );

      for( unsigned iJet=0; iJet<selectedJets.size(); ++iJet ) {
        AnalysisJet thisJet = selectedJets[iJet];
        h1_mJet_2jet->Fill( thisJet.M(), eventWeight );
        TLorentzVector ZJet = diLepton+thisJet;
        h1_mZllJet_2jet->Fill( ZJet.M(), eventWeight );
        h1_mtZllJet_2jet->Fill( ZJet.Mt(), eventWeight );
        if( leptType==0 ) {
          h1_mZmumuJet_2jet->Fill( diLepton.Pt(), eventWeight );
          h1_mtZmumuJet_2jet->Fill( diLepton.Pt(), eventWeight );
        } else {
          h1_mZeeJet_2jet->Fill( diLepton.Pt(), eventWeight );
          h1_mtZeeJet_2jet->Fill( diLepton.Pt(), eventWeight );
        }
        h1_met_2jet->Fill( pfMet, eventWeight );
      }

    } else if( selectedJets.size()==3 ) {

      h1_ptZll_3jet->Fill( diLepton.Pt(), eventWeight );
      if( leptType==0 ) h1_ptZmumu_3jet->Fill( diLepton.Pt(), eventWeight );
      else if( leptType==1 ) h1_ptZee_3jet->Fill( diLepton.Pt(), eventWeight );

      for( unsigned iJet=0; iJet<selectedJets.size(); ++iJet ) {
        AnalysisJet thisJet = selectedJets[iJet];
        h1_mJet_3jet->Fill( thisJet.M(), eventWeight );
        TLorentzVector ZJet = diLepton+thisJet;
        h1_mZllJet_3jet->Fill( ZJet.M(), eventWeight );
        h1_mtZllJet_3jet->Fill( ZJet.Mt(), eventWeight );
        if( leptType==0 ) {
          h1_mZmumuJet_3jet->Fill( diLepton.Pt(), eventWeight );
          h1_mtZmumuJet_3jet->Fill( diLepton.Pt(), eventWeight );
        } else {
          h1_mZeeJet_3jet->Fill( diLepton.Pt(), eventWeight );
          h1_mtZeeJet_3jet->Fill( diLepton.Pt(), eventWeight );
        }
        h1_met_3jet->Fill( pfMet, eventWeight );
      }

    } else if( selectedJets.size()==4 ) {

      h1_ptZll_4jet->Fill( diLepton.Pt(), eventWeight );
      if( leptType==0 ) h1_ptZmumu_4jet->Fill( diLepton.Pt(), eventWeight );
      else if( leptType==1 ) h1_ptZee_4jet->Fill( diLepton.Pt(), eventWeight );

      for( unsigned iJet=0; iJet<selectedJets.size(); ++iJet ) {
        AnalysisJet thisJet = selectedJets[iJet];
        h1_mJet_4jet->Fill( thisJet.M(), eventWeight );
        TLorentzVector ZJet = diLepton+thisJet;
        h1_mZllJet_4jet->Fill( ZJet.M(), eventWeight );
        h1_mtZllJet_4jet->Fill( ZJet.Mt(), eventWeight );
        if( leptType==0 ) {
          h1_mZmumuJet_4jet->Fill( diLepton.Pt(), eventWeight );
          h1_mtZmumuJet_4jet->Fill( diLepton.Pt(), eventWeight );
        } else {
          h1_mZeeJet_4jet->Fill( diLepton.Pt(), eventWeight );
          h1_mtZeeJet_4jet->Fill( diLepton.Pt(), eventWeight );
        }
        h1_met_4jet->Fill( pfMet, eventWeight );
      }

    }


    // now look for best W->jj
    if( selectedJets.size()>1 ) {

      AnalysisJet bestJet1, bestJet2;

      float Wmass = 80.399;
      bool firstOne=true;
      float bestMass=99999;

      for( unsigned iJet=0; iJet<selectedJets.size()-1; ++iJet ) {
        for( unsigned jJet=iJet+1; jJet<selectedJets.size(); ++jJet ) {
          TLorentzVector diJet = selectedJets[iJet] + selectedJets[jJet];
          if( fabs(diJet.M()-Wmass)<fabs(bestMass-Wmass) || firstOne ) {
            firstOne=false;
            bestJet1 = selectedJets[iJet];
            bestJet2 = selectedJets[jJet];
            bestMass = diJet.M();
          }
        } //for j
      } //for i



      // ------------------------
      //   KINEMATIC FIT: BEGIN
      // ------------------------


      TMatrixD m_jet1(3,3);
      TMatrixD m_jet2(3,3);

      m_jet1(0,0) = 0.5*ErrEt (bestJet1.Et(), bestJet1.Eta()); // et
      m_jet1(1,1) = 0.5*ErrEta(bestJet1.Et(), bestJet1.Eta()); // eta
      m_jet1(2,2) = 0.5*ErrPhi(bestJet1.Et(), bestJet1.Eta()); // phi
      m_jet2(0,0) = 0.5*ErrEt (bestJet2.Et(), bestJet2.Eta()); // et
      m_jet2(1,1) = 0.5*ErrEta(bestJet2.Et(), bestJet2.Eta()); // eta
      m_jet2(2,2) = 0.5*ErrPhi(bestJet2.Et(), bestJet2.Eta()); // phi

      TFitParticleEtEtaPhi *fitJet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &bestJet1, &m_jet1 );
      TFitParticleEtEtaPhi *fitJet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &bestJet2, &m_jet2 );
      
      TFitConstraintM *mCons_jets = new TFitConstraintM( "ZMassConstraint_jets", "ZMass-Constraint", 0, 0 , 91.19);
      mCons_jets->addParticles1( fitJet1, fitJet2 );

      TKinFitter* fitter_jets = new TKinFitter("fitter_jets", "fitter_jets");
      fitter_jets->addMeasParticle( fitJet1 );
      fitter_jets->addMeasParticle( fitJet2 );
      fitter_jets->addConstraint( mCons_jets );

      //Set convergence criteria
      fitter_jets->setMaxNbIter( 30 );
      fitter_jets->setMaxDeltaS( 1e-2 );
      fitter_jets->setMaxF( 1e-1 );
      fitter_jets->setVerbosity(0);

      //Perform the fit
      fitter_jets->fit();


      TLorentzVector jet1_kinfit(*fitJet1->getCurr4Vec());
      TLorentzVector jet2_kinfit(*fitJet2->getCurr4Vec());
      TLorentzVector Wjj_kinfit = jet1_kinfit + jet2_kinfit;

      TLorentzVector ZW= diLepton + bestJet1 + bestJet2;
      TLorentzVector ZW_kinfit = diLepton + Wjj_kinfit;

      h1_mZllW->Fill( ZW.M(), eventWeight );
      h1_mZllW_kinfit->Fill( ZW_kinfit.M(), eventWeight );
      if( leptType==0 ) {
        h1_mZmumuW->Fill( ZW.M(), eventWeight );
        h1_mZmumuW_kinfit->Fill( ZW_kinfit.M(), eventWeight );
      } else {
        h1_mZeeW->Fill( ZW.M(), eventWeight );
        h1_mZeeW_kinfit->Fill( ZW_kinfit.M(), eventWeight );
      }

    } // if 2 or more jets


/*
      vh1_rmsCandJet1[jet1PtBin]->Fill( jet1.rmsCand, eventWeight );
      vh1_ptDJet1[jet1PtBin]->Fill( jet1.ptD, eventWeight );
      vh1_nChargedJet1[jet1PtBin]->Fill( jet1.nCharged, eventWeight );
      vh1_nNeutralJet1[jet1PtBin]->Fill( jet1.nNeutral, eventWeight );
      float QGLikelihoodJet1 = qglikeli->computeQGLikelihood( jet1.Pt(), jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );

      vh1_QGLikelihoodJet1[jet1PtBin]->Fill( QGLikelihoodJet1, eventWeight );
      h1_QGLikelihoodJet1->Fill( QGLikelihoodJet1, eventWeight );
      if( fabs(jet1.Eta())<2. ) h1_QGLikelihoodJet1_eta2->Fill(QGLikelihoodJet1, eventWeight);
      
      vh1_rmsCandJet2[jet2PtBin]->Fill( jet2.rmsCand, eventWeight );
      vh1_ptDJet2[jet2PtBin]->Fill( jet2.ptD, eventWeight );
      vh1_nChargedJet2[jet2PtBin]->Fill( jet2.nCharged, eventWeight );
      vh1_nNeutralJet2[jet2PtBin]->Fill( jet2.nNeutral, eventWeight );
      float QGLikelihoodJet2 = qglikeli->computeQGLikelihood( jet2.Pt(), jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );

      vh1_QGLikelihoodJet2[jet2PtBin]->Fill( QGLikelihoodJet2, eventWeight );
      h1_QGLikelihoodJet2->Fill( QGLikelihoodJet2, eventWeight );
      if( fabs(jet2.Eta())<2. ) h1_QGLikelihoodJet2_eta2->Fill(QGLikelihoodJet2, eventWeight);

      float QGLikelihoodProd = QGLikelihoodJet1*QGLikelihoodJet2;
      float QGLikelihoodRevProd = (1.-QGLikelihoodJet1)*(1.-QGLikelihoodJet2);

      h1_QGLikelihoodProd->Fill( QGLikelihoodProd, eventWeight );
      
      h1_QGLikelihoodRevProd->Fill( QGLikelihoodRevProd, eventWeight );
        
      h2_QGLikelihoodJet1_vs_Jet2->Fill( QGLikelihoodJet2, QGLikelihoodJet1, eventWeight );

      if( QGLikelihoodJet1>0.93 && QGLikelihoodJet2>0.93 ) h1_QGLikelihoodProd_hi->Fill(QGLikelihoodProd, eventWeight); 

      // last step of selection:
      // QG and helicity LD's

      if( QGLikelihoodProd < QGLikelihoodProd_thresh_ ) continue;
      if( helicityLD < helicityLD_thresh_ ) continue;


      // *****************************************
      // *****  PASSED ANALYSIS SELECTION ********
      // *****************************************


        // fill histograms:

        if( jet1.Pt()>jet2.Pt() ) {
          h1_ptJet1->Fill( jet1.Pt(), eventWeight );
          h1_ptJet2->Fill( jet2.Pt(), eventWeight );
        } else {
          h1_ptJet1->Fill( jet2.Pt(), eventWeight );
          h1_ptJet2->Fill( jet1.Pt(), eventWeight );
        }
        h1_eMuonsJet1->Fill( jet1.muonEnergyFraction, eventWeight );
        h1_eMuonsJet2->Fill( jet2.muonEnergyFraction, eventWeight );
        h1_eElectronsJet1->Fill( jet1.electronEnergyFraction, eventWeight );
        h1_eElectronsJet2->Fill( jet2.electronEnergyFraction, eventWeight );

        h1_ptLept1->Fill( lept1.Pt(), eventWeight );
        h1_ptLept2->Fill( lept2.Pt(), eventWeight );
        h1_deltaRjj->Fill( jet1.DeltaR(jet2), eventWeight);
        h1_ptZll->Fill( diLepton.Pt(), eventWeight);
        h1_ptWjj->Fill( bestZDiJet.Pt(), eventWeight);
        if( leptType==0 )
          h1_mZmumu->Fill( diLepton.M(), eventWeight );
        else
          h1_mZee->Fill( diLepton.M(), eventWeight );
        h1_mZll->Fill( diLepton.M(), eventWeight);
        h1_mWjj->Fill( bestZDiJet.M(), eventWeight);
        h1_mZZ_UL->Fill(ZZ.M(), eventWeight);
        h1_mZZ_hiMass->Fill(ZZ.M(), eventWeight);
      //  h1_mZZ_highestMass->Fill(ZZ.M(), eventWeight);
        h1_mZZ_medMass->Fill(ZZ.M(), eventWeight);
        h1_mZZ_300Mass->Fill(ZZ.M(), eventWeight);
        h1_mZZ_WjjMassConstr_hiMass->Fill(ZZ_constr.M(), eventWeight);
        h1_mZZ_WjjMassConstr_medMass->Fill(ZZ_constr.M(), eventWeight);
        h1_mZZ_WjjMassConstr_300Mass->Fill(ZZ_constr.M(), eventWeight);
        h1_mZZ_UL_kinfit->Fill(ZZ_kinfit.M(), eventWeight);
        h1_mZZ_kinfit_hiMass->Fill(ZZ_kinfit.M(), eventWeight);
      //  h1_mZZ_kinfit_highestMass->Fill(ZZ_kinfit.M(), eventWeight);
        h1_mZZ_kinfit_medMass->Fill(ZZ_kinfit.M(), eventWeight);
        h1_mZZ_kinfit_300Mass->Fill(ZZ_kinfit.M(), eventWeight);
        h2_mWjj_vs_mZZ->Fill( ZZ.M(), bestZDiJet.M() );
        h2_mWjj_vs_mZZ_kinfit->Fill( ZZ_kinfit.M(), bestZDiJet.M() );


        h1_deltaRZZ->Fill(bestZDiJet.DeltaR(diLepton), eventWeight);

        if( ZZ_kinfit.M()>390. && ZZ_kinfit.M()<460. ) {
        h1_ptZZ->Fill( ZZ.Pt(), eventWeight );
        h1_ptZZ_kinfit->Fill( ZZ_kinfit.Pt(), eventWeight );
        h1_etaZZ->Fill( ZZ.Eta(), eventWeight );
        h1_etaZZ_kinfit->Fill( ZZ_kinfit.Eta(), eventWeight );
        }

        h1_cosThetaStar->Fill(hangles.helCosThetaStar, eventWeight);
        h1_cosTheta1->Fill(hangles.helCosTheta1, eventWeight);
        h1_cosTheta2->Fill(hangles.helCosTheta2, eventWeight);
        h1_phi->Fill(hangles.helPhi, eventWeight);
        h1_phi1->Fill(hangles.helPhi1, eventWeight);

        h1_cosThetaStar_kinfit->Fill(hangles_kinfit.helCosThetaStar, eventWeight);
        h1_cosTheta1_kinfit->Fill(hangles_kinfit.helCosTheta1, eventWeight);
        h1_cosTheta2_kinfit->Fill(hangles_kinfit.helCosTheta2, eventWeight);
        h1_phi_kinfit->Fill(hangles_kinfit.helPhi, eventWeight);
        h1_phi1_kinfit->Fill(hangles_kinfit.helPhi1, eventWeight);



        //match to partons:
        int partFlavor1=0;
        float deltaRmin1=999.;
        for(unsigned iPart=0; iPart<nPart; ++iPart ) {
          TLorentzVector thisPart;
          thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
          float thisDeltaR = jet1.DeltaR(thisPart);
          if( thisDeltaR<deltaRmin1 ) {
            partFlavor1 = pdgIdPart[iPart];
            deltaRmin1 = thisDeltaR;
          }
        }
        h1_deltaR_part1->Fill(deltaRmin1, eventWeight);
        h1_partFlavorJet1->Fill( partFlavor1, eventWeight );

        float deltaRmin2=999.;
        int partFlavor2=0;
        for(unsigned iPart=0; iPart<nPart; ++iPart ) {
          TLorentzVector thisPart;
          thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
          float thisDeltaR = jet2.DeltaR(thisPart);
          if( thisDeltaR<deltaRmin2 ) {
            partFlavor2 = pdgIdPart[iPart];
            deltaRmin2 = thisDeltaR;
          }
        }
        h1_deltaR_part2->Fill(deltaRmin2, eventWeight);
        h1_partFlavorJet2->Fill( partFlavor2, eventWeight );

      } //if passes selection

    } //if selected jet pairs
*/




  } //for entries
 
  ofs_events.close();

  outFile_->cd();

  h1_run->Write();

  h1_ptMu1->Write();
  h1_ptMu2->Write();
  h1_etaMu->Write();
  h1_phiMu->Write();

  h1_ptMu1_inBump->Write();
  h1_ptMu2_inBump->Write();
  h1_etaMu1_inBump->Write();
  h1_etaMu2_inBump->Write();
  h1_phiMu1_inBump->Write();
  h1_phiMu2_inBump->Write();

  h1_ptMu1_outBump->Write();
  h1_ptMu2_outBump->Write();
  h1_etaMu1_outBump->Write();
  h1_etaMu2_outBump->Write();
  h1_phiMu1_outBump->Write();
  h1_phiMu2_outBump->Write();

  h1_deltaPhiMu1Met_inBump->Write();
  h1_deltaPhiMu1Met_outBump->Write();

  h2_etaMu_vs_ptMu->Write();
  h2_phiMu_vs_ptMu->Write();
  h2_phiMu_vs_etaMu->Write();

  h2_etaMu_vs_ptZll->Write();
  h2_phiMu_vs_ptZll->Write();
  h2_phiMu_vs_etaZll->Write();

  h1_ptEle1->Write();
  h1_ptEle2->Write();
  h1_etaEle->Write();
  h1_phiEle->Write();

  h2_etaEle_vs_ptEle->Write();
  h2_phiEle_vs_ptEle->Write();
  h2_phiEle_vs_etaEle->Write();

  h2_etaEle_vs_ptZll->Write();
  h2_phiEle_vs_ptZll->Write();
  h2_phiEle_vs_etaZll->Write();


  h1_ptZll->Write();
  h1_etaZll->Write();
  h1_phiZll->Write();
  h1_ptZmumu->Write();
  h1_etaZmumu->Write();
  h1_phiZmumu->Write();
  h1_ptZee->Write();
  h1_etaZee->Write();
  h1_phiZee->Write();

  h1_nJets->Write();
  h1_nJets_inBump->Write();
  h1_nJets_outBump->Write();

  h1_ptZll_0jet->Write();
  h1_ptZll_1jet->Write();
  h1_ptZll_2jet->Write();
  h1_ptZll_3jet->Write();
  h1_ptZll_4jet->Write();
  h1_ptZmumu_0jet->Write();
  h1_ptZmumu_1jet->Write();
  h1_ptZmumu_2jet->Write();
  h1_ptZmumu_3jet->Write();
  h1_ptZmumu_4jet->Write();
  h1_ptZee_0jet->Write();
  h1_ptZee_1jet->Write();
  h1_ptZee_2jet->Write();
  h1_ptZee_3jet->Write();
  h1_ptZee_4jet->Write();

  h1_eleEnergyFractionJet_inBump->Write();
  h1_muonEnergyFractionJet_inBump->Write();
  h1_eleEnergyFractionJet_outBump->Write();
  h1_muonEnergyFractionJet_outBump->Write();

  h1_mJet_inBump->Write();
  h1_mJet_outBump->Write();
  h1_mJet_1jet->Write();
  h1_mJet_2jet->Write();
  h1_mJet_3jet->Write();
  h1_mJet_4jet->Write();

  h1_mZll_inBump->Write();
  h1_mZmumu_inBump->Write();
  h1_mZee_inBump->Write();
  h1_mZll_outBump->Write();
  h1_mZmumu_outBump->Write();
  h1_mZee_outBump->Write();

  h1_mZllJet_inBump->Write();
  h1_mZllJet_outBump->Write();
  h1_mZllJet_1jet->Write();
  h1_mZllJet_2jet->Write();
  h1_mZllJet_3jet->Write();
  h1_mZllJet_4jet->Write();

  h1_mZmumuJet_inBump->Write();
  h1_mZmumuJet_outBump->Write();
  h1_mZmumuJet_1jet->Write();
  h1_mZmumuJet_2jet->Write();
  h1_mZmumuJet_3jet->Write();
  h1_mZmumuJet_4jet->Write();

  h1_mZeeJet_inBump->Write();
  h1_mZeeJet_outBump->Write();
  h1_mZeeJet_1jet->Write();
  h1_mZeeJet_2jet->Write();
  h1_mZeeJet_3jet->Write();
  h1_mZeeJet_4jet->Write();
  
  h1_mtZllMet_inBump->Write();
  h1_mtZllMet_outBump->Write();
  h1_mtZmumuMet_inBump->Write();
  h1_mtZmumuMet_outBump->Write();
  h1_mtZeeMet_inBump->Write();
  h1_mtZeeMet_outBump->Write();
  
  h1_mtZllJet_inBump->Write();
  h1_mtZllJet_outBump->Write();
  h1_mtZllJet_1jet->Write();
  h1_mtZllJet_2jet->Write();
  h1_mtZllJet_3jet->Write();
  h1_mtZllJet_4jet->Write();

  h1_mtZmumuJet_inBump->Write();
  h1_mtZmumuJet_outBump->Write();
  h1_mtZmumuJet_1jet->Write();
  h1_mtZmumuJet_2jet->Write();
  h1_mtZmumuJet_3jet->Write();
  h1_mtZmumuJet_4jet->Write();

  h1_mtZeeJet_inBump->Write();
  h1_mtZeeJet_outBump->Write();
  h1_mtZeeJet_1jet->Write();
  h1_mtZeeJet_2jet->Write();
  h1_mtZeeJet_3jet->Write();
  h1_mtZeeJet_4jet->Write();
  
  h1_met_inBump->Write();
  h1_met_outBump->Write();
  h1_met_1jet->Write();
  h1_met_2jet->Write();
  h1_met_3jet->Write();
  h1_met_4jet->Write();

  h1_mZllW->Write();
  h1_mZllW_kinfit->Write();
  h1_mZmumuW->Write();
  h1_mZmumuW_kinfit->Write();
  h1_mZeeW->Write();
  h1_mZeeW_kinfit->Write();

  outFile_->Close();

} // finalize()



void Ntp1Finalizer_WZ::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;

  if( selectionType_=="presel" ) {

    ptLept1_thresh_ = 10.;
    ptLept2_thresh_ = 10.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;
    mZll_threshLo_ = 60.;
    mZll_threshHi_ = 120.;
    etaZll_thresh_ = 10.;

  } else if( selectionType_=="loose" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet_thresh_ = 30.;
    etaJet_thresh_ = 2.4;
    mZll_threshLo_ = 60.;
    mZll_threshHi_ = 120.;
    etaZll_thresh_ = 10.;

  } else if( selectionType_=="MIT" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet_thresh_ = 30.;
    etaJet_thresh_ = 5.;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    etaZll_thresh_ = 1.;

  } else if( selectionType_=="adish" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 2.1;
    etaLept2_thresh_ = 2.1;
    ptJet_thresh_ = 30.;
    etaJet_thresh_ = 5.;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    etaZll_thresh_ = 10.;

  } else if( selectionType_=="lept40" ) {

    ptLept1_thresh_ = 30.;
    ptLept2_thresh_ = 30.;
    etaLept1_thresh_ = 2.1;
    etaLept2_thresh_ = 2.1;
    ptJet_thresh_ = 30.;
    etaJet_thresh_ = 5.;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    etaZll_thresh_ = 10.;

  } else if( selectionType_=="lept40" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 40.;
    etaLept1_thresh_ = 2.1;
    etaLept2_thresh_ = 2.1;
    ptJet_thresh_ = 30.;
    etaJet_thresh_ = 5.;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    etaZll_thresh_ = 10.;

  } else if( selectionType_=="lept50" ) {

    ptLept1_thresh_ = 50.;
    ptLept2_thresh_ = 50.;
    etaLept1_thresh_ = 2.1;
    etaLept2_thresh_ = 2.1;
    ptJet_thresh_ = 30.;
    etaJet_thresh_ = 5.;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    etaZll_thresh_ = 10.;

  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType


// error functions for jets:

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



// error functions for PFCands:

Double_t ErrEt(Float_t Et, Float_t Eta, int particleType) {
  Double_t InvPerr2, S, N, C;

  // following values taken from PFMet significance code: 
  // RecoMET/METProducers/python/METSigParams_cfi.py
  if( particleType==1 ) { // charged hadrons
    N = 0.05;
    S = 0.;
    C = 0.;
  } else if( particleType==2 ) { // electrons
    N = 0.05;
    S = 0.;
    C = 0.;
  } else if( particleType==3 ) { // muons
    N = 0.05;
    S = 0.;
    C = 0.;
  } else if( particleType==4 ) { // photons
    N = 0.25; //slightly more conservative:
    S = 0.03;
    C = 0.01;
  } else if( particleType==5 ) { // neutral hadrons
    N = 0.;
    S = 1.22;
    C = 0.05;
  } else if( particleType==6 ) { // HF hadrons
    N = 0.;
    S = 1.22;
    C = 0.05;
  } else if( particleType==7 ) { // HF EM
    N = 0.;
    S = 1.22;
    C = 0.05;
  } else {
    std::cout << "IMPROBABLE particleType '" << particleType << "' found! PROBLEM!!" << std::endl;
    return 0.;
  }
  InvPerr2 = (N * N) + (S * S) * Et + (C * C) * Et * Et;
  return InvPerr2;
}





Double_t ErrEta(Float_t Et, Float_t Eta, int particleType) {

  Double_t E = Et*cosh(Eta);
  Double_t InvPerr2, a, b, c;

  InvPerr2 = ErrPhi( Et, Eta, particleType );
  // ErrPhi is in radians^2:
  InvPerr2 = sqrt(InvPerr2);
  // have to convert it to pseudorapidity:
  InvPerr2 = -log( tan (0.5*(InvPerr2) ) );
  // go back to sqyared:
  InvPerr2 *= InvPerr2;

//if( particleType==1 ) { // charged hadrons
//} else if( particleType==2 ) { // electrons
//} else if( particleType==3 ) { // muons
//} else if( particleType==4 ) { // photons
//} else if( particleType==5 ) { // neutral hadrons
//} else if( particleType==6 ) { // HF hadrons
//} else if( particleType==7 ) { // HF EM
//}

  return InvPerr2;
}

Double_t ErrPhi(Float_t Et, Float_t Eta, int particleType) {

  Double_t E = Et*cosh(Eta);
  Double_t InvPerr2, a, b, c;


  if( particleType==1 ) { // charged hadrons

    InvPerr2 = 0.002*Et;

  } else if( particleType==2 ) { // electrons

    InvPerr2 = 0.002*Et;

  } else if( particleType==3 ) { // muons

    InvPerr2 = 0.002*Et;

  } else if( particleType==4 ) { // photons

    bool inBarrel = fabs(Eta) < 1.475;

    a = 4.; // see DN 2007-011
    b = 11.;
    c = 0.6;
    if(!inBarrel) {
      float crystalScaleFactor = 28.6/22.;
      a *= crystalScaleFactor;
      b *= crystalScaleFactor;
      c *= crystalScaleFactor;
    }
    // this is error in mm's:
    InvPerr2 = a*a/E + b*b/(E*E) + c*c;
    // convert to radians:
    float crystalSize = (inBarrel) ? 22. : 28.6;
    InvPerr2 /= (crystalSize*crystalSize);
    // (one crystal is one degree):
    float degToRad = 3.14159/180.;
    InvPerr2 *= (degToRad*degToRad);

  //// have to convert it to radians:
  //float ECAL_radius = 1300.; //in mm
  //InvPerr2 = TMath::Pi()/2. - atan( InvPerr2/ECAL_radius );

  } else if( particleType==5 ) { // neutral hadrons

    InvPerr2 = 0.0251*Et;

  } else if( particleType==6 ) { // HF hadrons

    InvPerr2 = 0.0251*Et;

  } else if( particleType==7 ) { // HF EM

    InvPerr2 = 0.0251*Et;

  } else {
    std::cout << "IMPROBABLE particleType '" << particleType << "' found! PROBLEM!!" << std::endl;
  }


  return InvPerr2;
}


void print(TKinFitter *fitter)
{
  std::cout << "=============================================" << std ::endl;
  std::cout << "-> Number of measured Particles  : " << fitter->nbMeasParticles() << std::endl;
  std::cout << "-> Number of unmeasured particles: " << fitter->nbUnmeasParticles() << std::endl;
  std::cout << "-> Number of constraints         : " << fitter->nbConstraints() << std::endl;
  std::cout << "-> Number of degrees of freedom  : " << fitter->getNDF() << std::endl;
  std::cout << "-> Number of parameters A        : " << fitter->getNParA() << std::endl;
  std::cout << "-> Number of parameters B        : " << fitter->getNParB() << std::endl;
  std::cout << "-> Maximum number of iterations  : " << fitter->getMaxNumberIter() << std::endl;
  std::cout << "-> Maximum deltaS                : " << fitter->getMaxDeltaS() << std::endl;
  std::cout << "-> Maximum F                     : " << fitter->getMaxF() << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std ::endl;
  std::cout << "-> Status                        : " << fitter->getStatus() << std::endl;
  std::cout << "-> Number of iterations          : " << fitter->getNbIter() << std::endl;
  std::cout << "-> S                             : " << fitter->getS() << std::endl;
  std::cout << "-> F                             : " << fitter->getF() << std::endl;
  std::cout << "=============================================" << std ::endl;
}


int getNJets( int nPairs ) {

  int nJets = 0;

  int i=1;
  bool found=false;
  while( !found) {
    nJets += i;
    i++;
    if( nJets== nPairs ) found=true;
  }

  return i;

}


std::vector<TH1D*> getHistoVector(int nPtBins, Double_t *ptBins, std::string histoName, int nBins, float xMin, float xMax ) {

  std::vector<TH1D*> returnVector;

  for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

    Double_t ptMin = ptBins[iPtBin];
    Double_t ptMax = ptBins[iPtBin+1];

    char histoName_thisBin[200];
    sprintf( histoName_thisBin, "%s_pt_%.0lf_%.0lf", histoName.c_str(), ptMin, ptMax);

    TH1D* newHisto = new TH1D(histoName_thisBin, "", nBins, xMin, xMax);
    newHisto->Sumw2();

    returnVector.push_back(newHisto);

  } //for pt bins

  return returnVector;

}



HelicityLikelihoodDiscriminant::HelicityAngles computeHelicityAngles(TLorentzVector leptMinus, TLorentzVector leptPlus, TLorentzVector jet1, TLorentzVector jet2 ) {

  HelicityLikelihoodDiscriminant::HelicityAngles returnAngles;

  TLorentzVector Zll = leptPlus + leptMinus;
  TLorentzVector Wjj = jet1 + jet2;

  TLorentzVector Higgs = Wjj + Zll;

  
  // define lept1 as negatively charged lepton:
  TLorentzVector lept1 = leptMinus;
  TLorentzVector lept2 = leptPlus;

  // no charge for jets (well, not really)
  // so choose jet with positive scalar product with Wjj 
  // in its restframe:
  TLorentzVector jet1_Wjjstar_tmp(jet1);
  jet1_Wjjstar_tmp.Boost(-Wjj.BoostVector());
  if( jet1_Wjjstar_tmp.Phi()<0. ) { //swap them
    TLorentzVector jet1_tmp(jet1);
    TLorentzVector jet2_tmp(jet2);
    jet1 = jet2_tmp;
    jet2 = jet1_tmp;
  }


  //     BOOSTS:

  // boosts in Higgs CoM frame:
  TLorentzVector lept1_Hstar(lept1);
  lept1_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector lept2_Hstar(lept2);
  lept2_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector jet1_Hstar(jet1);
  jet1_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector jet2_Hstar(jet2);
  jet2_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector Zll_Hstar(Zll);
  Zll_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector Wjj_Hstar(Wjj);
  Wjj_Hstar.Boost(-Higgs.BoostVector());

  // boosts in Zll CoM frame:
  TLorentzVector lept1_Zllstar(lept1);
  lept1_Zllstar.Boost(-Zll.BoostVector());
  TLorentzVector H_Zllstar(Higgs);
  H_Zllstar.Boost(-Zll.BoostVector());

  // boosts in Wjj CoM frame:
  TLorentzVector jet1_Wjjstar(jet1);
  jet1_Wjjstar.Boost(-Wjj.BoostVector());
  TLorentzVector H_Wjjstar(Higgs);
  H_Wjjstar.Boost(-Wjj.BoostVector());


  returnAngles.helCosThetaStar = Zll_Hstar.CosTheta();


  TVector3 v_pbeamLAB( 0.0, 0.0, 1.0 );

  //cross prod beam x Zll
  TVector3 v_1 = (v_pbeamLAB.Cross(  (Zll_Hstar.Vect()).Unit()) ).Unit();//versor normal to z-z' plane


  //v_2 = cross prod l1 x l2 = versor normal to Zll decay plane
  // careful to the order: L1, the z-axis and Z->ll make a right-handed (non-orthogonal) frame (xyz); at the end we want the angle btw x and y
  TVector3 v_2((Zll_Hstar.Vect().Cross(lept1_Hstar.Vect().Unit())).Unit());


  //v_3 = similar to v_2, BUT
  //now, if we want a right-handed set of Unit-vectors, keeping the same direction of the z-axis
  //we must swap the direction of one of the other two vectors of the Z bosons. 
  //Keeping the same direction of the z-axis
  //means measuring phiZll and phiWjj w.r.t. to the same v_1 vector (i.e. w.r.t. the same z'-Zll plane)
  TVector3 v_3(((-1.0*Wjj_Hstar.Vect()).Cross(jet1_Hstar.Vect().Unit())).Unit()) ;

  //in other terms: we can define v_3 as above and then do the crss prod with v_1
  //or define v_3 in a way consistent with v_2 and then do the cross product with a newly defined
  //Unit vector v_4 =  (v_pbeamLAB.Cross(  (WjjboostedX->momentum()).Unit()) ).Unit();//versor normal to z-Wjj plane
 
  // helphiZll:
  float phiZll = fabs( acos(v_1.Dot(v_2)) );
  if(v_pbeamLAB.Dot(v_2)>0.0)phiZll=-1.0*phiZll;
  else phiZll=+1.0*phiZll;

  // helphiWjj:
  float phiWjj = fabs( acos(v_1.Dot(v_3)) );
  if(v_pbeamLAB.Dot(v_3)>0.0)phiWjj=+1.0*phiWjj; 
  else phiWjj=-1.0*phiWjj;


  float phi1 = phiZll;


  //phi
  float phi = fabs( acos(v_2.Dot(v_3)) );//two-fold ambiguity when doing the acos + pi ambiguity from sign of v_3 
  if(lept1_Hstar.Vect().Dot(v_3)>0.0)phi= +1.0 * phi;
  else phi= -1.0 * phi;

  returnAngles.helPhi1 = phi1;
  returnAngles.helPhi = phi;


  returnAngles.helCosTheta1 =  (-1.0*(lept1_Zllstar.X()* H_Zllstar.X()+
                                   lept1_Zllstar.Y()* H_Zllstar.Y()+
                                   lept1_Zllstar.Z()* H_Zllstar.Z())/
                                  (lept1_Zllstar.Vect().Mag()* H_Zllstar.Vect().Mag())  );


  returnAngles.helCosTheta2 =  fabs( (jet1_Wjjstar.X()* H_Wjjstar.X()+
                                   jet1_Wjjstar.Y()* H_Wjjstar.Y()+
                                   jet1_Wjjstar.Z()* H_Wjjstar.Z())/
                                  (jet1_Wjjstar.Vect().Mag()* H_Wjjstar.Vect().Mag())  );

  returnAngles.mzz = Higgs.M();


  return returnAngles;

}
