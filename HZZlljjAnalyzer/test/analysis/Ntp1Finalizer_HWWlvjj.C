#include "Ntp1Finalizer_HWWlvjj.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "CommonTools/fitTools.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"
#include "KinematicFit/DiJetKinFitter.h"
#include "KinematicFit/LeptonNeutrinoKinFitter.h"


inline double delta_phi(double phi1, double phi2);
float getPzRight( TLorentzVector lepton, float pxPFMet, float pyPFMet,  TLorentzVector neuMC );
float getPzWrong( TLorentzVector lepton, float pxPFMet, float pyPFMet,  TLorentzVector neuMC );
float getPz( TLorentzVector lepton, float pxPFMet, float pyPFMet, TLorentzVector jet1, TLorentzVector jet2);
float getPzMH( TLorentzVector lepton, float pxPFMet, float pyPFMet, TLorentzVector jet1, TLorentzVector jet2);
float getPzATLAS( TLorentzVector lepton, float pxPFMet, float pyPFMet, TLorentzVector jet1, TLorentzVector jet2);

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

  float ptGen;
  float etaGen;
float phiGen;
  float eGen;

  //btags:
  float trackCountingHighEffBJetTag;
  float trackCountingHighPurBJetTag;
  float simpleSecondaryVertexHighEffBJetTag;
  float simpleSecondaryVertexHighPurBJetTag;
  float jetBProbabilityBJetTag;
  float jetProbabilityBJetTag;

};



HelicityLikelihoodDiscriminant::HelicityAngles computeHelicityAngles(TLorentzVector leptMinus, TLorentzVector leptPlus, TLorentzVector jet1, TLorentzVector jet2 );



int getNJets( int nPairs );

std::vector<TH1D*> getHistoVector(int nPtBins, Double_t *ptBins, std::string histoName, int nBins, float xMin, float xMax );


// constructor:

Ntp1Finalizer_HWWlvjj::Ntp1Finalizer_HWWlvjj( const std::string& dataset, const std::string& selectionType, const std::string& leptType ) : Ntp1Finalizer( "HWWlvjj", dataset, leptType ) {

  if( leptType!="ALL" && leptType!="MU" && leptType!="ELE" ) {
    std::cout << "Lept type '" << leptType << "' currently not supported. Exiting." << std::endl;
    exit(9176);
  }

  leptType_ = leptType;

  setSelectionType(selectionType);

  std::string fullFlags = selectionType_ + "_" + leptType_;
  this->set_flags(fullFlags); //this is for the outfile name

}





void Ntp1Finalizer_HWWlvjj::finalize() {

  if( outFile_==0 ) this->createOutputFile();


  //  std::string btagFileName = "btagfile_" + selectionType_ + ".root";
  //  TFile* btagFile = TFile::Open(btagFileName.c_str(), "RECREATE");

  TH1D* h1_Ev_Presel = new TH1D("Ev_Presel","",1,0.,1.);
  TH1D* h1_Ev_Jet1 = new TH1D("Ev_Jet1","",1,0.,1.);
  TH1D* h1_Ev_Jet2 = new TH1D("Ev_Jet2","",1,0.,1.);
  TH1D* h1_Ev_EtaJet = new TH1D("Ev_EtaJet","",1,0.,1.);
  TH1D* h1_Ev_DrJet = new TH1D("Ev_DrJet","",1,0.,1.);
  TH1D* h1_Ev_mWjj = new TH1D("Ev_mWjj","",1,0.,1.);
  TH1D* h1_Ev_PtWjj = new TH1D("Ev_PtWjj","",1,0.,1.);
  TH1D* h1_Ev_DileptonPt = new TH1D("Ev_DileptonPt","",1,0.,1.);
  TH1D* h1_Ev_Lept1 = new TH1D("Ev_Lept1","",1,0.,1.);
  TH1D* h1_Ev_Lept2 = new TH1D("Ev_Lept2","",1,0.,1.);
  TH1D* h1_Ev_EtaLept = new TH1D("Ev_EtaLept","",1,0.,1.);
  TH1D* h1_Ev_mtW = new TH1D("Ev_mtW","",1,0.,1.);
  TH1D* h1_Ev_DrLept = new TH1D("Ev_DrLept","",1,0.,1.);
  TH1D* h1_Ev_mWW = new TH1D("Ev_mWW","",1,0.,1.);
  TH1D* h1_Ev_nCounter = new TH1D("Ev_nCounter","",1,0.,1.);

  TH1D* h1_Ev_nCounterW = new TH1D("Ev_nCounterW","",1,0.,1.);
  TH1D* h1_Ev_nEventsPassed_fb_kinfit = new TH1D("Ev_nEventsPassed_fb_kinfit","",1,0.,1.);
  TH1D* h1_Ev_nEventsPassed_kinfit = new TH1D("Ev_nEventsPassed_kinfit","",1,0.,1.);
  // TH1D* h1_Ev_nEventsPassed_kinfit_antiBtag = new TH1D("Ev_nEventsPassed_kinfit_antiBtag","",1,0.,1.);
  TH1D* h1_Ev_nEventsPassed_fb_nokinfit = new TH1D("Ev_nEventsPassed_fb_nokinfit","",1,0.,1.);
  TH1D* h1_Ev_nEventsPassed_nokinfit = new TH1D("Ev_nEventsPassed_nokinfit","",1,0.,1.);

  TH1F* h1_run = new TH1F("run", "", 15149, 132440, 147589);


  TH1D* h1_ptLept1_JustPresel = new TH1D("ptLept1_JustPresel", "", 50, 10., 220.);
  h1_ptLept1_JustPresel->Sumw2();
  TH1D* h1_ptLept2_JustPresel = new TH1D("ptLept2_JustPresel", "", 50, 10., 220.);
  h1_ptLept2_JustPresel->Sumw2();
  TH1D* h1_etaLept1_presel = new TH1D("etaLept1_presel", "", 25, -2.5, 2.5);
  h1_etaLept1_presel->Sumw2();
  TH1D* h1_etaLept2_presel = new TH1D("etaLept2_presel", "", 25, -2.5, 2.5);
  h1_etaLept2_presel->Sumw2();
  TH1D* h1_mtW_JustPresel = new TH1D("mtW_JustPresel", "", 50, 0.,120. );
  h1_mtW_JustPresel->Sumw2();

  TH1D* h1_ptJet_all_presel = new TH1D("ptJet_all_presel", "", 27, 30., 400.);
  h1_ptJet_all_presel->Sumw2();
  TH1D* h1_etaJet_all_presel = new TH1D("etaJet_all_presel", "", 25, -5., 5.);
  h1_etaJet_all_presel->Sumw2();

  const int nPtBins = 20;
  Double_t ptBins[nPtBins+1];
  fitTools::getBins_int( nPtBins+1, ptBins, 15., 1000. );

  std::vector<TH1D*> vh1_ptDJet_all_presel = getHistoVector(nPtBins, ptBins, "ptDJet_all_presel", 60, 0., 1.);
  std::vector<TH1D*> vh1_rmsCandJet_all_presel = getHistoVector(nPtBins, ptBins, "rmsCandJet_all_presel", 60, 0., 0.07);
  std::vector<TH1D*> vh1_nChargedJet_all_presel = getHistoVector(nPtBins, ptBins, "nChargedJet_all_presel", 41, -0.5, 40.5);
  std::vector<TH1D*> vh1_nNeutralJet_all_presel = getHistoVector(nPtBins, ptBins, "nNeutralJet_all_presel", 41, -0.5, 40.5);

  TH1D* h1_nJets_presel = new TH1D("nJets_presel", "", 7, 1.5, 8.5);
  h1_nJets_presel->Sumw2();
  TH1D* h1_nPairs_presel = new TH1D("nPairs_presel", "", 21, 0.5, 21.5);
  h1_nPairs_presel->Sumw2();

  TH1D* h1_ptLept1= new TH1D("ptLept1", "", 25, 20., 300.);
  h1_ptLept1->Sumw2();
  TH1D* h1_ptLept2= new TH1D("ptLept2", "", 25, 20., 150.);
  h1_ptLept2->Sumw2();
  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 50, 20., 220.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 50, 20., 140.);
  h1_ptJet2->Sumw2();
  TH1D* h1_mtW= new TH1D("mtW", "", 50, 0., 120.);
  h1_mtW->Sumw2();
  TH1D* h1_eElectronsJet1 = new TH1D("eElectronsJet1", "", 60, 0., 1.0001);
  h1_eElectronsJet1->Sumw2();
  TH1D* h1_eElectronsJet2 = new TH1D("eElectronsJet2", "", 60, 0., 1.0001);
  h1_eElectronsJet2->Sumw2();
  TH1D* h1_eMuonsJet1 = new TH1D("eMuonsJet1", "", 60, 0., 1.0001);
  h1_eMuonsJet1->Sumw2();
  TH1D* h1_eMuonsJet2 = new TH1D("eMuonsJet2", "", 60, 0., 1.0001);
  h1_eMuonsJet2->Sumw2();

  TH1D* h1_ptResoJet1_beforeKin = new TH1D("ptResoJet1_beforeKin", "", 100, -1., 1.);
  h1_ptResoJet1_beforeKin->Sumw2();
  TH1D* h1_ptResoJet2_beforeKin = new TH1D("ptResoJet2_beforeKin", "", 100, -1., 1.);
  h1_ptResoJet2_beforeKin->Sumw2();

  TH1D* h1_ptResoJet1_afterKin = new TH1D("ptResoJet1_afterKin", "", 100, -1., 1.);
  h1_ptResoJet1_afterKin->Sumw2();
  TH1D* h1_ptResoJet2_afterKin = new TH1D("ptResoJet2_afterKin", "", 100, -1., 1.);
  h1_ptResoJet2_afterKin->Sumw2();

  TH1D* h1_ptWreso_beforeKin = new TH1D("ptWreso_beforeKin", "", 100, -1., 1.);
  h1_ptWreso_beforeKin->Sumw2();
  TH1D* h1_ptWreso_afterKin = new TH1D("ptWreso_afterKin", "", 100, -1., 1.);
  h1_ptWreso_afterKin->Sumw2();

  TH1D* h1_pzResoNeut_GetPz = new TH1D("pzResoNeut_GetPz", "", 100, -5., 5.);
  h1_pzResoNeut_GetPz->Sumw2();
  TH1D* h1_pzResoNeut_KinFit = new TH1D("pzResoNeut_KinFit", "", 100, -5., 5.);
  h1_pzResoNeut_KinFit->Sumw2();

//TH1D* h1_ptJetRecoil = new TH1D("ptJetRecoil", "", 27, 30., 400.);
//h1_ptJetRecoil->Sumw2();
//TH1D* h1_ptDJetRecoil = new TH1D("ptDJetRecoil", "", 60, 0., 1.);
//h1_ptDJetRecoil->Sumw2();
//TH1D* h1_rmsCandJetRecoil = new TH1D("rmsCandJetRecoil", "", 60, 0., 0.07);
//h1_rmsCandJetRecoil->Sumw2();
//TH1D* h1_nChargedJetRecoil = new TH1D("nChargedJetRecoil", "", 41, -0.5, 40.5);
//h1_nChargedJetRecoil->Sumw2();
//TH1D* h1_nNeutralJetRecoil = new TH1D("nNeutralJetRecoil", "", 41, -0.5, 40.5);
//h1_nNeutralJetRecoil->Sumw2();

  int nBins_invMass = 80;
  float invMassMin = 30.;
  float invMassMax = 120.;
  float invMassMin_ll = 60.;

  TH1D* h1_mWjj_all_presel = new TH1D("mWjj_all_presel", "", 50, 20, 200.);
  h1_mWjj_all_presel->Sumw2();

  TH1D* h1_deltaRjj_all_presel = new TH1D("deltaRjj_all_presel", "", 50, 0.5, 5.);
  h1_deltaRjj_all_presel->Sumw2();
  TH1D* h1_deltaRll_JustPresel = new TH1D("deltaRll_JustPresel", "", 20, 0.5, 5.);
  h1_deltaRll_JustPresel->Sumw2();

  TH1D* h1_deltaRjj_JustPresel= new TH1D("deltaRjj_JustPresel", "", 50, 0.5, 5.);
  h1_deltaRjj_JustPresel->Sumw2();
  TH1D* h1_ptJet1_JustPresel= new TH1D("ptJet1_JustPresel", "", 50, 20., 220.);
  h1_ptJet1_JustPresel->Sumw2();	
  TH1D* h1_ptJet2_JustPresel= new TH1D("ptJet2_JustPresel", "", 50, 20., 140.);
  h1_ptJet2_JustPresel->Sumw2();
	
  TH1D* h1_mWjj_JustPresel= new TH1D("mWjj_JustPresel", "", 50, 20., 200.);
   h1_mWjj_JustPresel->Sumw2();

  TH1D* h1_mWll = new TH1D("mWll", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWll->Sumw2();
  TH1D* h1_mWll_presel = new TH1D("mWll_presel", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWll_presel->Sumw2();
  TH1D* h1_mWll_presel_0jets = new TH1D("mWll_presel_0jets", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWll_presel_0jets->Sumw2();

  TH1D* h1_mWmumu = new TH1D("mWmumu", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWmumu->Sumw2();
  TH1D* h1_mWmumu_presel = new TH1D("mWmumu_presel", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWmumu_presel->Sumw2();
  TH1D* h1_mWmumu_presel_0jets = new TH1D("mWmumu_presel_0jets", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWmumu_presel_0jets->Sumw2();
  TH1D* h1_mWee = new TH1D("mWee", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWee->Sumw2();
  TH1D* h1_mWee_presel = new TH1D("mWee_presel", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWee_presel->Sumw2();
  TH1D* h1_mWee_presel_0jets = new TH1D("mWee_presel_0jets", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mWee_presel_0jets->Sumw2();

  TH1D* h1_ptDJet1 = new TH1D("ptDJet1", "", 60, 0., 1.0001);
  h1_ptDJet1->Sumw2();
  TH1D* h1_ptDJet2 = new TH1D("ptDJet2", "", 60, 0., 1.0001);
  h1_ptDJet2->Sumw2();

  TH1D* h1_deltaR_part1 = new TH1D("deltaR_part1", "", 50, 0., 0.8);
  h1_deltaR_part1->Sumw2();
  //TH1D* h1_partFlavorJet1= new TH1D("partFlavorJet1", "", 30, -7.5, 22.5);
  TH1D* h1_partFlavorJet1= new TH1D("partFlavorJet1", "", 38, -15.5, 22.5);
  h1_partFlavorJet1->Sumw2();
  std::vector<TH1D*> vh1_ptDJet1 = getHistoVector(nPtBins, ptBins, "ptDJet1", 60, 0., 1.);
  std::vector<TH1D*> vh1_rmsCandJet1 = getHistoVector(nPtBins, ptBins, "rmsCandJet1", 60, 0., 0.1);
  std::vector<TH1D*> vh1_nChargedJet1 = getHistoVector(nPtBins, ptBins, "nChargedJet1", 51, -0.5, 50.5);
  std::vector<TH1D*> vh1_nNeutralJet1 = getHistoVector(nPtBins, ptBins, "nNeutralJet1", 51, -0.5, 50.5);
  std::vector<TH1D*> vh1_QGLikelihoodJet1 = getHistoVector(nPtBins, ptBins, "QGLikelihoodJet1", 60, 0., 1.);
  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet1->Sumw2();
  //@ TH1D* h1_QGLikelihoodJet1_antiBtag_SSVhighEff = new TH1D("QGLikelihoodJet1_antiBtag_SSVhighEff", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet1_antiBtag_SSVhighEff->Sumw2();
  // TH1D* h1_QGLikelihoodJet1_antiBtag_SSVhighPur = new TH1D("QGLikelihoodJet1_antiBtag_SSVhighPur", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet1_antiBtag_SSVhighPur->Sumw2();
  // TH1D* h1_QGLikelihoodJet1_antiBtag_SSVhighPurhighEff = new TH1D("QGLikelihoodJet1_antiBtag_SSVhighPurhighEff", "", 60, 0., 1.0001);
  //h1_QGLikelihoodJet1_antiBtag_SSVhighPurhighEff->Sumw2();
  // TH1D* h1_QGLikelihoodJet1_antiBtag_TChighEff = new TH1D("QGLikelihoodJet1_antiBtag_TChighEff", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet1_antiBtag_TChighEff->Sumw2();
  //TH1D* h1_QGLikelihoodJet1_antiBtag_TChighPur = new TH1D("QGLikelihoodJet1_antiBtag_TChighPur", "", 60, 0., 1.0001);
  //h1_QGLikelihoodJet1_antiBtag_TChighPur->Sumw2();
  //TH1D* h1_QGLikelihoodJet1_antiBtag_jetBProb = new TH1D("QGLikelihoodJet1_antiBtag_jetBProb", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet1_antiBtag_jetBProb->Sumw2();
  // TH1D* h1_QGLikelihoodJet1_antiBtag_jetProb = new TH1D("QGLikelihoodJet1_antiBtag_jetProb", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet1_antiBtag_jetProb->Sumw2();

  TH1D* h1_simpleSecondaryVertexHighEffBJetTagJet1 = new TH1D("simpleSecondaryVertexHighEffBJetTagJet1", "", 50, -1.5, 4.);
  h1_simpleSecondaryVertexHighEffBJetTagJet1->Sumw2();
  TH1D* h1_simpleSecondaryVertexHighPurBJetTagJet1 = new TH1D("simpleSecondaryVertexHighPurBJetTagJet1", "", 50, -1.5, 4.);
  h1_simpleSecondaryVertexHighPurBJetTagJet1->Sumw2();
  TH1D* h1_jetBProbabilityBJetTagJet1 = new TH1D("jetBProbabilityBJetTagJet1", "", 50, 0., 8.);
  h1_jetBProbabilityBJetTagJet1->Sumw2();
  TH1D* h1_jetProbabilityBJetTagJet1 = new TH1D("jetProbabilityBJetTagJet1", "", 50, 0., 2.5);
  h1_jetProbabilityBJetTagJet1->Sumw2();


  TH1D* h1_deltaR_part2 = new TH1D("deltaR_part2", "", 60, 0., 0.8);
  h1_deltaR_part2->Sumw2();
  //TH1D* h1_partFlavorJet2= new TH1D("partFlavorJet2", "", 30, -7.5, 22.5);
  TH1D* h1_partFlavorJet2= new TH1D("partFlavorJet2", "", 38, -15.5, 22.5);
  h1_partFlavorJet2->Sumw2();

  std::vector<TH1D*> vh1_ptDJet2 = getHistoVector(nPtBins, ptBins, "ptDJet2", 60, 0., 1.);
  std::vector<TH1D*> vh1_rmsCandJet2 = getHistoVector(nPtBins, ptBins, "rmsCandJet2", 60, 0., 0.1);
  std::vector<TH1D*> vh1_nChargedJet2 = getHistoVector(nPtBins, ptBins, "nChargedJet2", 51, -0.5, 50.5);
  std::vector<TH1D*> vh1_nNeutralJet2 = getHistoVector(nPtBins, ptBins, "nNeutralJet2", 51, -0.5, 50.5);
  std::vector<TH1D*> vh1_QGLikelihoodJet2 = getHistoVector(nPtBins, ptBins, "QGLikelihoodJet2", 60, 0., 1.);
  TH1D* h1_QGLikelihoodJet2 = new TH1D("QGLikelihoodJet2", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet2->Sumw2();
  //  TH1D* h1_QGLikelihoodJet2_antiBtag_SSVhighEff = new TH1D("QGLikelihoodJet2_antiBtag_SSVhighEff", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet2_antiBtag_SSVhighEff->Sumw2();
  // TH1D* h1_QGLikelihoodJet2_antiBtag_SSVhighPur = new TH1D("QGLikelihoodJet2_antiBtag_SSVhighPur", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet2_antiBtag_SSVhighPur->Sumw2();
  // TH1D* h1_QGLikelihoodJet2_antiBtag_SSVhighPurhighEff = new TH1D("QGLikelihoodJet2_antiBtag_SSVhighPurhighEff", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet2_antiBtag_SSVhighPurhighEff->Sumw2();
  // TH1D* h1_QGLikelihoodJet2_antiBtag_TChighEff = new TH1D("QGLikelihoodJet2_antiBtag_TChighEff", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet2_antiBtag_TChighEff->Sumw2();
  // TH1D* h1_QGLikelihoodJet2_antiBtag_TChighPur = new TH1D("QGLikelihoodJet2_antiBtag_TChighPur", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet2_antiBtag_TChighPur->Sumw2();
  // TH1D* h1_QGLikelihoodJet2_antiBtag_jetBProb = new TH1D("QGLikelihoodJet2_antiBtag_jetBProb", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet2_antiBtag_jetBProb->Sumw2();
  // TH1D* h1_QGLikelihoodJet2_antiBtag_jetProb = new TH1D("QGLikelihoodJet2_antiBtag_jetProb", "", 60, 0., 1.0001);
  // h1_QGLikelihoodJet2_antiBtag_jetProb->Sumw2();

  TH1D* h1_simpleSecondaryVertexHighEffBJetTagJet2 = new TH1D("simpleSecondaryVertexHighEffBJetTagJet2", "", 50, -1.5, 4.);
  h1_simpleSecondaryVertexHighEffBJetTagJet2->Sumw2();
  TH1D* h1_simpleSecondaryVertexHighPurBJetTagJet2 = new TH1D("simpleSecondaryVertexHighPurBJetTagJet2", "", 50, -1.5, 4.);
  h1_simpleSecondaryVertexHighPurBJetTagJet2->Sumw2();
  TH1D* h1_jetBProbabilityBJetTagJet2 = new TH1D("jetBProbabilityBJetTagJet2", "", 50, 0., 8.);
  h1_jetBProbabilityBJetTagJet2->Sumw2();
  TH1D* h1_jetProbabilityBJetTagJet2 = new TH1D("jetProbabilityBJetTagJet2", "", 50, 0., 2.5);
  h1_jetProbabilityBJetTagJet2->Sumw2();

  TH1D* h1_QGLikelihoodProd = new TH1D("QGLikelihoodProd", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd->Sumw2();
  // TH1D* h1_QGLikelihoodProd_antiBtag = new TH1D("QGLikelihoodProd_antiBtag", "", 60, 0., 1.0001);
  // h1_QGLikelihoodProd_antiBtag->Sumw2();
  TH1D* h1_QGLikelihoodProd_hi = new TH1D("QGLikelihoodProd_hi", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_hi->Sumw2();
  TH1D* h1_QGLikelihoodRevProd = new TH1D("QGLikelihoodRevProd", "", 60, 0., 1.0001);
  h1_QGLikelihoodRevProd->Sumw2();
  TH2D* h2_QGLikelihoodJet1_vs_Jet2 = new TH2D("QGLikelihoodJet1_vs_Jet2", "", 60, 0., 1.0001, 60, 0., 1.0001);
  h2_QGLikelihoodJet1_vs_Jet2->Sumw2();


  TH1D* h1_mWjj= new TH1D("mWjj", "", nBins_invMass, 20., 200.);
  h1_mWjj->Sumw2();
  TH1D* h1_mWjj_loChiSquareProb= new TH1D("mWjj_loChiSquareProb", "", 100, 30., 200.);
  h1_mWjj_loChiSquareProb->Sumw2();
  TH1D* h1_mWjj_hiChiSquareProb= new TH1D("mWjj_hiChiSquareProb", "", 100, 30., 200.);
  h1_mWjj_hiChiSquareProb->Sumw2();

  TH1D* h1_ptWll_JustPresel = new TH1D("ptWll_JustPresel", "", 50, 0., 300.);
  h1_ptWll_JustPresel->Sumw2();
  TH1D* h1_ptWjj_all_presel = new TH1D("ptWjj_all_presel", "", 50, 0., 3000.);
  h1_ptWjj_all_presel->Sumw2();

  TH1D* h1_ptWll = new TH1D("ptWll", "", 50, 140., 340.);
  h1_ptWll->Sumw2();
  TH1D* h1_ptWjj = new TH1D("ptWjj", "", 40, 140., 340.);
  h1_ptWjj->Sumw2();

  TH1D* h1_deltaRjj= new TH1D("deltaRjj", "", 36, 0.5, 5.);
  h1_deltaRjj->Sumw2();

  TH1D* h1_cosThetaStar = new TH1D("cosThetaStar", "", 45, -1.001, 1.001);
  h1_cosThetaStar->Sumw2();
  TH1D* h1_cosTheta1 = new TH1D("cosTheta1", "", 45, -1.001, 1.001);
  h1_cosTheta1->Sumw2();
  TH1D* h1_cosTheta2 = new TH1D("cosTheta2", "", 45, 0., 1.001);
  h1_cosTheta2->Sumw2();
  TH1D* h1_phi = new TH1D("phi", "", 45, -3.1416, 3.1416);
  h1_phi->Sumw2();
  TH1D* h1_phi1 = new TH1D("phi1", "", 45, -3.1416, 3.1416);
  h1_phi1->Sumw2();

  TH1D* h1_cosThetaStar_kinfit = new TH1D("cosThetaStar_kinfit", "", 45, -1.0001, 1.001);
  h1_cosThetaStar_kinfit->Sumw2();
  TH1D* h1_cosTheta1_kinfit = new TH1D("cosTheta1_kinfit", "", 45, -1.001, 1.001);
  h1_cosTheta1_kinfit->Sumw2();
  TH1D* h1_cosTheta2_kinfit = new TH1D("cosTheta2_kinfit", "", 45, 0., 1.001);
  h1_cosTheta2_kinfit->Sumw2();
  TH1D* h1_phi_kinfit = new TH1D("phi_kinfit", "", 45, -3.1416, 3.1416);
  h1_phi_kinfit->Sumw2();
  TH1D* h1_phi1_kinfit = new TH1D("phi1_kinfit", "", 45, -3.1416, 3.1416);
  h1_phi1_kinfit->Sumw2();

  
  TH1D* h1_kinfit_chiSquare = new TH1D("kinfit_chiSquare", "", 60, 0., 10.);
  h1_kinfit_chiSquare->Sumw2();
  TH1D* h1_kinfit_chiSquareProb = new TH1D("kinfit_chiSquareProb", "", 60, 0., 1.0001);
  h1_kinfit_chiSquareProb->Sumw2();

  TH1D* h1_helicityLD = new TH1D("helicityLD", "", 60, 0., 1.);
  h1_helicityLD->Sumw2();
  TH1D* h1_helicityLD_kinfit = new TH1D("helicityLD_kinfit", "", 60, 0., 1.);
  h1_helicityLD_kinfit->Sumw2();

  TH1D* h1_deltaRWW= new TH1D("deltaRWW", "", 60, 0., 6.);
  h1_deltaRWW->Sumw2();

 
  TH1D* h1_mWW_GetPz = new TH1D("mWW_GetPz", "", 100, 200., 600.);
  h1_mWW_GetPz->Sumw2();//###
  TH1D* h1_mWW_kinLept = new TH1D("mWW_kinLept", "", 100, 200., 600.);
  h1_mWW_kinLept->Sumw2();

  TH1D* h1_mWW_nokinfit = new TH1D("mWW_nokinfit", "", 100, 200., 700.);
  h1_mWW_nokinfit->Sumw2();//###
  TH1D* h1_mWW_kinfit = new TH1D("mWW_kinfit", "", 100, 200., 700.);
  h1_mWW_kinfit->Sumw2();
  TH1D* h1_mWW_kinfit_cut = new TH1D("mWW_kinfit_cut", "", 100, 200., 700.);
  h1_mWW_kinfit_cut->Sumw2();
  TH1D* h1_mWW_hiChiSquareProb = new TH1D("mWW_hiChiSquareProb", "", 200, 100., 700.);
  h1_mWW_hiChiSquareProb->Sumw2();
  TH1D* h1_mWW_loChiSquareProb = new TH1D("mW_loChiSquareProb", "", 200, 100., 700.);
  h1_mWW_loChiSquareProb->Sumw2();
  TH1D* h1_mWW_mWjj_cut = new TH1D("mWW_mWjj_cut", "", 200, 100., 700.);
  h1_mWW_mWjj_cut->Sumw2();
  TH1D* h1_mWW_mWjj_notcut = new TH1D("mWW_mWjj_notcut", "", 200, 100., 700.);
  h1_mWW_mWjj_notcut->Sumw2();
  TH1D* h1_mWW_UL = new TH1D("mWW_UL", "", 900, 100., 1000.);
  h1_mWW_UL->Sumw2();
  TH1D* h1_mWW_UL_kinfit = new TH1D("mWW_UL_kinfit", "", 900, 100., 1000.);
  h1_mWW_UL_kinfit->Sumw2();
  TH1D* h1_mWW_hiMass= new TH1D("mWW_hiMass", "", 90, 250., 700.);
  h1_mWW_hiMass->Sumw2();
  TH1D* h1_mWW_kinfit_hiMass= new TH1D("mWW_kinfit_hiMass", "", 90, 250., 700.);
  h1_mWW_kinfit_hiMass->Sumw2();
//TH1D* h1_mWW_highestMass= new TH1D("mWW_highestMass", "", 70, 350., 700.);
//h1_mWW_highestMass->Sumw2();
//TH1D* h1_mWW_kinfit_highestMass= new TH1D("mWW_kinfit_highestMass", "", 70, 350., 700.);
//h1_mWW_kinfit_highestMass->Sumw2();
  TH1D* h1_mWW_WjjMassConstr_hiMass  = new TH1D("mWW_WjjMassConstr_hiMass", "", 70, 250., 600.);
  h1_mWW_WjjMassConstr_hiMass->Sumw2();
  TH1D* h1_mWW_300Mass= new TH1D("mWW_300Mass", "", 50, 250., 500.);
  h1_mWW_300Mass->Sumw2();
  TH1D* h1_mWW_WjjMassConstr_300Mass  = new TH1D("mWW_WjjMassConstr_300Mass", "", 50, 250., 500.);
  h1_mWW_WjjMassConstr_300Mass->Sumw2();
  TH1D* h1_mWW_kinfit_300Mass= new TH1D("mWW_kinfit_300Mass", "", 50, 250., 500.);
  h1_mWW_kinfit_300Mass->Sumw2();
  TH1D* h1_mWW_medMass= new TH1D("mWW_medMass", "", 70, 150., 350.);
  h1_mWW_medMass->Sumw2();
  TH1D* h1_mWW_WjjMassConstr_medMass  = new TH1D("mWW_WjjMassConstr_medMass", "", 70, 150., 350.);
  h1_mWW_WjjMassConstr_medMass->Sumw2();
  TH1D* h1_mWW_kinfit_medMass= new TH1D("mWW_kinfit_medMass", "", 70, 150., 350.);
  h1_mWW_kinfit_medMass->Sumw2();

  TH1D* h1_ptWW  = new TH1D("ptWW", "", 100, 0., 300.);
  h1_ptWW->Sumw2();
  TH1D* h1_ptWW_kinfit  = new TH1D("ptWW_kinfit", "", 100, 0., 300.);
  h1_ptWW_kinfit->Sumw2();
  TH1D* h1_etaWW  = new TH1D("etaWW", "", 100, -5.5, 5.5);
  h1_etaWW->Sumw2();
  TH1D* h1_etaWW_kinfit  = new TH1D("etaWW_kinfit", "", 100, -5.5, 5.5);
  h1_etaWW_kinfit->Sumw2();
  

  TH1D* h1_mWW_MCassoc  = new TH1D("mWW_MCassoc", "", 100, 200., 600.);
  h1_mWW_MCassoc->Sumw2();
  TH1D* h1_mWW_MCassoc_WjjMassConstr  = new TH1D("mWW_MCassoc_WjjMassConstr", "", 100, 200., 600.);
  h1_mWW_MCassoc_WjjMassConstr->Sumw2();
  TH1D* h1_mWW_MCassoc_kinfit  = new TH1D("mWW_MCassoc_kinfit", "", 100, 200., 600.);
  h1_mWW_MCassoc_kinfit->Sumw2();
  TH1D* h1_mWW_MCassoc_kinfit_cands  = new TH1D("mWW_MCassoc_kinfit_cands", "", 100, 200., 600.);
  h1_mWW_MCassoc_kinfit_cands->Sumw2();

  TH1D* h1_partFlavor_tight = new TH1D("partFlavor_tight", "", 31, -8.5, 22.5);

  TH2D* h2_mWjj_vs_mWW = new TH2D("mWjj_vs_mWW", "", 100, 200., 600., 100, 60., 110.);//traslati di 10 per m(W)
  h2_mWjj_vs_mWW->Sumw2();
  TH2D* h2_mWjj_vs_mWW_kinfit = new TH2D("mWjj_vs_mWW_kinfit", "", 100, 200., 600., 100, 60., 110.);
  h2_mWjj_vs_mWW_kinfit->Sumw2();


  TH1D* h1_deltaE_ch = new TH1D("deltaE_ch", "", 100, -1., 1.);
  TH1D* h1_deltaE_gamma = new TH1D("deltaE_gamma", "", 100, -1., 1.);
  TH1D* h1_deltaE_nh = new TH1D("deltaE_nh", "", 100, -1., 1.);

  TH1D* h1_deltaEta_ch    = new TH1D("deltaEta_ch", "", 100, -1., 1.);
  TH1D* h1_deltaEta_gamma = new TH1D("deltaEta_gamma", "", 100, -1., 1.);
  TH1D* h1_deltaEta_nh    = new TH1D("deltaEta_nh", "", 100, -1., 1.);

  TH1D* h1_deltaPhi_ch    = new TH1D("deltaPhi_ch", "", 100, -1., 1.);
  TH1D* h1_deltaPhi_gamma = new TH1D("deltaPhi_gamma", "", 100, -1., 1.);
  TH1D* h1_deltaPhi_nh    = new TH1D("deltaPhi_nh", "", 100, -1., 1.);

  TH1D* h1_deltaPt_ch = new TH1D("deltaPt_ch", "", 100, -1., 1.);
  TH1D* h1_deltaPt_gamma = new TH1D("deltaPt_gamma", "", 100, -1., 1.);
  TH1D* h1_deltaPt_nh = new TH1D("deltaPt_nh", "", 100, -1., 1.);

  TH1D* h1_energyMet= new TH1D("energyMet","",100,0.,400);
//@@@
  TH1D* h1_resoPz = new TH1D("resoPz", "", 100, -5., 5.);
  TH1D* h1_resoPt = new TH1D("resoPt", "", 100, -1., 1.);

  TH1D* h1_resoPzMH = new TH1D("resoPzMH", "", 100, -5., 5.);
  TH1D* h1_resoPzATLAS = new TH1D("resoPzATLAS", "", 100, -5., 5.);
  TH1D* h1_resoPzRight = new TH1D("resoPzRight", "", 100, -5., 5.);
  TH1D* h1_resoPzWrong = new TH1D("resoPzWrong", "", 100, -5., 5.);
  TH1D* h1_resoPtW = new TH1D("resoPtW", "", 100, -5., 5.);
  TH1D* h1_resomH = new TH1D("resomH", "", 100, -1.5, 1.5);
  TH1D* h1_resomHATLAS = new TH1D("resomHATLAS", "", 100, -1.5, 1.5);
  TH1D* h1_resoPzW = new TH1D("resoPzW", "", 100, -5., 5.);
  TH1D* h1_resoPzWMH = new TH1D("resoPzWMH", "", 100, -5., 5.);
  TH1D* h1_resoPzWATLAS = new TH1D("resoPzWATLAS", "", 100, -5., 5.);

  TH1D* h1_FindPz_EtaR =  new TH1D("FindPz_EtaR", "", 100, -5., 5.);
  TH1D* h1_FindPz_EtaW =  new TH1D("FindPz_EtaW", "", 100, -5., 5.);
  TH1D* h1_FindPz_DRR =  new TH1D("FindPz_DRR", "", 100, -5., 5.);
  TH1D* h1_FindPz_DRW =  new TH1D("FindPz_DRW", "", 100, -5., 5.);
  TH1D* h1_FindPz_EtaWnW = new TH1D("FindPz_EtaWnW", "", 100, -5., 5.);
  TH1D* h1_FindPz_EtaWnR = new TH1D("FindPz_EtaWnR", "", 100, -5., 5.);
  TH1D* h1_FindPz_PzW = new TH1D("FindPz_PzW", "", 100, 0., 500.);
  TH1D* h1_FindPz_PzR = new TH1D("FindPz_PzR", "", 100, 0., 500.);   
  TH1D* h1_FindPz_DeltaEta = new TH1D("FindPz_DeltaEta", "", 100, -5., 5.);
  TH1D* h1_FindPz_DeltaEtaWn = new TH1D("FindPz_DeltaEtaWn", "", 100, -4., 4.);
  TH1D* h1_FindPz_DeltaPz = new TH1D("FindPz_DeltaPz", "", 100, 5., 300.);
  TH1D* h1_FindPz_WNeutWW = new TH1D("FindPz_WNeutWW", "", 100, -1., 1.);
  TH1D* h1_FindPz_WNeutWR = new TH1D("FindPz_WNeutWR", "", 100, -1., 1.);


//Double_t ptBins[16];
//fitTools::getBins_int( 16, ptBins, 20., 500.);

//TProfile* hp_ptJetGenMean = new TProfile("ptJetGenMean", "", 15, ptBins);
//std::vector<TH1F*> h1_response_vs_pt = getResponseHistos("response", 16, ptBins);
//std::vector<TH1F*> h1_response_vs_pt_Rch050 = getResponseHistos("response_Rch050", 16, ptBins);
//std::vector<TH1F*> h1_response_vs_pt_Rch5070 = getResponseHistos("response_Rch5070", 16, ptBins);
//std::vector<TH1F*> h1_response_vs_pt_Rch70100 = getResponseHistos("response_Rch70100", 16, ptBins);



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

  Float_t energyPFMet;
  tree_->SetBranchAddress("energyPFMet", &energyPFMet);
  Float_t phiPFMet;
  tree_->SetBranchAddress("phiPFMet", &phiPFMet);
  Float_t pxPFMet;
  tree_->SetBranchAddress("pxPFMet", &pxPFMet);
  Float_t pyPFMet;
  tree_->SetBranchAddress("pyPFMet", &pyPFMet);

  int leptType;
  tree_->SetBranchAddress("leptType", &leptType);

  Float_t eLept;
  tree_->SetBranchAddress("eLept", &eLept);
  Float_t ptLept;
  tree_->SetBranchAddress("ptLept", &ptLept);
  Float_t etaLept;
  tree_->SetBranchAddress("etaLept", &etaLept);
  Float_t phiLept;
  tree_->SetBranchAddress("phiLept", &phiLept);
  Int_t chargeLept;
  tree_->SetBranchAddress("chargeLept", &chargeLept);

  Float_t eWllMC;
  tree_->SetBranchAddress("eWllMC", &eWllMC);
  Float_t ptWllMC;
  tree_->SetBranchAddress("ptWllMC", &ptWllMC);
  Float_t etaWllMC;
  tree_->SetBranchAddress("etaWllMC", &etaWllMC);
  Float_t phiWllMC;
  tree_->SetBranchAddress("phiWllMC", &phiWllMC);

  Float_t eHiggsMC;
  tree_->SetBranchAddress("eHiggsMC", &eHiggsMC);
  Float_t ptHiggsMC;
  tree_->SetBranchAddress("ptHiggsMC", &ptHiggsMC);
  Float_t etaHiggsMC;
  tree_->SetBranchAddress("etaHiggsMC", &etaHiggsMC);
  Float_t phiHiggsMC;
  tree_->SetBranchAddress("phiHiggsMC", &phiHiggsMC);

  Float_t eLeptMC;
  tree_->SetBranchAddress("eLeptMC", &eLeptMC);
  Float_t ptLeptMC;
  tree_->SetBranchAddress("ptLept1MC", &ptLeptMC);
  Float_t etaLeptMC;
  tree_->SetBranchAddress("etaLeptMC", &etaLeptMC);
  Float_t phiLeptMC;
  tree_->SetBranchAddress("phiLeptMC", &phiLeptMC);
  Float_t eNeuMC;
  tree_->SetBranchAddress("eNeuMC", &eNeuMC);
  Float_t ptNeuMC;
  tree_->SetBranchAddress("ptNeuMC", &ptNeuMC);
  Float_t etaNeuMC;
  tree_->SetBranchAddress("etaNeuMC", &etaNeuMC);
  Float_t phiNeuMC;
  tree_->SetBranchAddress("phiNeuMC", &phiNeuMC);

  Int_t nPairs;
  tree_->SetBranchAddress("nPairs", &nPairs);

  Int_t iJet1[50];
  tree_->SetBranchAddress("iJet1", iJet1);
  Float_t eJet1[50];
  tree_->SetBranchAddress("eJet1", eJet1);
  Float_t ptJet1[50];
  tree_->SetBranchAddress("ptJet1", ptJet1);
  Float_t etaJet1[50];
  tree_->SetBranchAddress("etaJet1", etaJet1);
  Float_t phiJet1[50];
  tree_->SetBranchAddress("phiJet1", phiJet1);
  Float_t eJet1Gen[50];
  tree_->SetBranchAddress("eJet1Gen", eJet1Gen);
  Float_t ptJet1Gen[50];
  tree_->SetBranchAddress("ptJet1Gen", ptJet1Gen);
  Float_t etaJet1Gen[50];
  tree_->SetBranchAddress("etaJet1Gen", etaJet1Gen);
  Float_t phiJet1Gen[50];
  tree_->SetBranchAddress("phiJet1Gen", phiJet1Gen);
  Float_t eChargedHadronsJet1[50];
  tree_->SetBranchAddress("eChargedHadronsJet1", eChargedHadronsJet1);
  Float_t rmsCandJet1[50];
  tree_->SetBranchAddress("rmsCandJet1", rmsCandJet1);
  Float_t ptDJet1[50];
  tree_->SetBranchAddress("ptDJet1", ptDJet1);
  Int_t nChargedJet1[50];
  tree_->SetBranchAddress("nChargedJet1", nChargedJet1);
  Int_t nNeutralJet1[50];
  tree_->SetBranchAddress("nNeutralJet1", nNeutralJet1);
  Float_t eMuonsJet1[50];
  tree_->SetBranchAddress("eMuonsJet1", eMuonsJet1);
  Float_t eElectronsJet1[50];
  tree_->SetBranchAddress("eElectronsJet1", eElectronsJet1);
  Float_t trackCountingHighEffBJetTagJet1[50];
  tree_->SetBranchAddress("trackCountingHighEffBJetTagJet1", trackCountingHighEffBJetTagJet1);
  Float_t trackCountingHighPurBJetTagJet1[50];
  tree_->SetBranchAddress("trackCountingHighPurBJetTagJet1", trackCountingHighPurBJetTagJet1);
  Float_t simpleSecondaryVertexHighEffBJetTagJet1[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagJet1", simpleSecondaryVertexHighEffBJetTagJet1);
  Float_t simpleSecondaryVertexHighPurBJetTagJet1[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagJet1", simpleSecondaryVertexHighPurBJetTagJet1);
  Float_t jetBProbabilityBJetTagJet1[50];
  tree_->SetBranchAddress("jetBProbabilityBJetTagJet1", jetBProbabilityBJetTagJet1);
  Float_t jetProbabilityBJetTagJet1[50];
  tree_->SetBranchAddress("jetProbabilityBJetTagJet1", jetProbabilityBJetTagJet1);

  Int_t nPFCand1;
  tree_->SetBranchAddress("nPFCand1", &nPFCand1);
  Float_t ePFCand1[100];
  tree_->SetBranchAddress("ePFCand1", ePFCand1);
  Float_t ptPFCand1[100];
  tree_->SetBranchAddress("ptPFCand1", ptPFCand1);
  Float_t etaPFCand1[100];
  tree_->SetBranchAddress("etaPFCand1", etaPFCand1);
  Float_t phiPFCand1[100];
  tree_->SetBranchAddress("phiPFCand1", phiPFCand1);
  Int_t particleTypePFCand1[100];
  tree_->SetBranchAddress("particleTypePFCand1", particleTypePFCand1);

  Int_t iJet2[50];
  tree_->SetBranchAddress("iJet2", iJet2);
  Float_t eJet2[50];
  tree_->SetBranchAddress("eJet2", eJet2);
  Float_t ptJet2[50];
  tree_->SetBranchAddress("ptJet2", ptJet2);
  Float_t etaJet2[50];
  tree_->SetBranchAddress("etaJet2", etaJet2);
  Float_t phiJet2[50];
  tree_->SetBranchAddress("phiJet2", phiJet2);
  Float_t eJet2Gen[50];
  tree_->SetBranchAddress("eJet2Gen", eJet2Gen);
  Float_t ptJet2Gen[50];
  tree_->SetBranchAddress("ptJet2Gen", ptJet2Gen);
  Float_t etaJet2Gen[50];
  tree_->SetBranchAddress("etaJet2Gen", etaJet2Gen);
  Float_t phiJet2Gen[50];
  tree_->SetBranchAddress("phiJet2Gen", phiJet2Gen);
  Float_t eChargedHadronsJet2[50];
  tree_->SetBranchAddress("eChargedHadronsJet2", eChargedHadronsJet2);
  Float_t rmsCandJet2[50];
  tree_->SetBranchAddress("rmsCandJet2", rmsCandJet2);
  Float_t ptDJet2[50];
  tree_->SetBranchAddress("ptDJet2", ptDJet2);
  Int_t nChargedJet2[50];
  tree_->SetBranchAddress("nChargedJet2", nChargedJet2);
  Int_t nNeutralJet2[50];
  tree_->SetBranchAddress("nNeutralJet2", nNeutralJet2);
  Float_t eMuonsJet2[50];
  tree_->SetBranchAddress("eMuonsJet2", eMuonsJet2);
  Float_t eElectronsJet2[50];
  tree_->SetBranchAddress("eElectronsJet2", eElectronsJet2);
  Float_t trackCountingHighEffBJetTagJet2[50];
  tree_->SetBranchAddress("trackCountingHighEffBJetTagJet2", trackCountingHighEffBJetTagJet2);
  Float_t trackCountingHighPurBJetTagJet2[50];
  tree_->SetBranchAddress("trackCountingHighPurBJetTagJet2", trackCountingHighPurBJetTagJet2);
  Float_t simpleSecondaryVertexHighEffBJetTagJet2[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagJet2", simpleSecondaryVertexHighEffBJetTagJet2);
  Float_t simpleSecondaryVertexHighPurBJetTagJet2[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagJet2", simpleSecondaryVertexHighPurBJetTagJet2);
  Float_t jetBProbabilityBJetTagJet2[50];
  tree_->SetBranchAddress("jetBProbabilityBJetTagJet2", jetBProbabilityBJetTagJet2);
  Float_t jetProbabilityBJetTagJet2[50];
  tree_->SetBranchAddress("jetProbabilityBJetTagJet2", jetProbabilityBJetTagJet2);

  Int_t nPFCand2;
  tree_->SetBranchAddress("nPFCand2", &nPFCand2);
  Float_t ePFCand2[100];
  tree_->SetBranchAddress("ePFCand2", ePFCand2);
  Float_t ptPFCand2[100];
  tree_->SetBranchAddress("ptPFCand2", ptPFCand2);
  Float_t etaPFCand2[100];
  tree_->SetBranchAddress("etaPFCand2", etaPFCand2);
  Float_t phiPFCand2[100];
  tree_->SetBranchAddress("phiPFCand2", phiPFCand2);
  Int_t particleTypePFCand2[100];
  tree_->SetBranchAddress("particleTypePFCand2", particleTypePFCand2);

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

  Float_t uncorrEnergyAK5Jet;
  tree_->SetBranchAddress("uncorrEnergyAK5Jet", &uncorrEnergyAK5Jet);//###
  Float_t SumEt;
  tree_->SetBranchAddress("SumEt", &SumEt);
  Float_t pzNeutrino;
  tree_->SetBranchAddress("pzNeutrino", &pzNeutrino);
  Float_t ptNeutrino;
  tree_->SetBranchAddress("ptNeutrino", &ptNeutrino);

  float nEventsPassed_fb_kinfit=0.;
  float nEventsPassed_fb_nokinfit=0.;
  int nEventsPassed_kinfit=0;
  // int nEventsPassed_kinfit_antiBtag=0;
  int nEventsPassed_nokinfit=0;

  int nEntries = tree_->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmshome/pernielu/CMSSW_3_8_7/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root", nPtBins);

  float nEventsTot = 0.;
  float nEvents_hiChiSquareProb = 0.;
  float nEvents_mWjj_cut = 0.;

  float nEvent_Presel=0;
  float nEvent_DileptPt=0;
  float nEvent_LeadLeptPt=0;
  float nEvent_SubleadLeptPt=0;
  float nEvent_EtaLept=0;
  float nEvent_EtaJet=0;
  float nEvent_DijetMass=0;
  float nEvent_DijetPt=0;
  float nEvent_LeadJetPt=0;
  float nEvent_SubleadJetPt=0;
  float nEvent_DrJetJet=0;
  float nEvent_DrLeptLept=0;
  float nEvent_mtW=0;
 
  for(int iEntry=0; iEntry<nEntries; ++iEntry) {
    
    if( (iEntry % 50000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;
    
    tree_->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    if( leptType_!="ALL" ) {
      if( leptType_=="ELE" && leptType==0 ) continue;
      if( leptType_=="MU" && leptType==1 ) continue;
    }



    bool isMC = (run<5);


    if( !isMC ) { 

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


    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_selected;
    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_JustPresel_selected;//###

    float cached_jetpt = 0.;

    bool ripet_Presel = false;
    bool ripet_LeadJetPt = false;//to avoid repetition in efficiency during jet selection
    bool ripet_SubleadJetPt = false;
    bool ripet_EtaJet = false;
    bool ripet_DrJetJet = false;
    bool ripet_DijetMass = false;
    bool ripet_DijetPt = false;

    for( unsigned iJetPair=0; iJetPair<nPairs; ++iJetPair ) {

      AnalysisJet jet1, jet2;
      jet1.SetPtEtaPhiE( ptJet1[iJetPair], etaJet1[iJetPair], phiJet1[iJetPair], eJet1[iJetPair]);
      jet2.SetPtEtaPhiE( ptJet2[iJetPair], etaJet2[iJetPair], phiJet2[iJetPair], eJet2[iJetPair]);

      jet1.rmsCand = rmsCandJet1[iJetPair];
      jet1.ptD = ptDJet1[iJetPair];
      jet1.nCharged = nChargedJet1[iJetPair];
      jet1.nNeutral = nNeutralJet1[iJetPair];
      jet1.muonEnergyFraction = eMuonsJet1[iJetPair]/jet1.Energy();
      jet1.electronEnergyFraction = eElectronsJet1[iJetPair]/jet1.Energy();

      jet1.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet1[iJetPair];
      jet1.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet1[iJetPair];
      jet1.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet1[iJetPair];
      jet1.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet1[iJetPair];
      jet1.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet1[iJetPair];
      jet1.jetProbabilityBJetTag               = jetProbabilityBJetTagJet1[iJetPair];


      jet1.ptGen = ptJet1Gen[iJetPair];
      jet1.etaGen = etaJet1Gen[iJetPair];
      jet1.phiGen = phiJet1Gen[iJetPair];
      jet1.eGen = eJet1Gen[iJetPair];

      jet2.rmsCand = rmsCandJet2[iJetPair];
      jet2.ptD = ptDJet2[iJetPair];
      jet2.nCharged = nChargedJet2[iJetPair];
      jet2.nNeutral = nNeutralJet2[iJetPair];
      jet2.muonEnergyFraction = eMuonsJet2[iJetPair]/jet2.Energy();
      jet2.electronEnergyFraction = eElectronsJet2[iJetPair]/jet2.Energy();

      jet2.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet2[iJetPair];
      jet2.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet2[iJetPair];
      jet2.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet2[iJetPair];
      jet2.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet2[iJetPair];
      jet2.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet2[iJetPair];
      jet2.jetProbabilityBJetTag               = jetProbabilityBJetTagJet2[iJetPair];

      jet2.ptGen =   ptJet2Gen[iJetPair];
      jet2.etaGen = etaJet2Gen[iJetPair];
      jet2.phiGen = phiJet2Gen[iJetPair];
      jet2.eGen =     eJet2Gen[iJetPair];

      TLorentzVector diJet = jet1 + jet2;

      jetPairs_JustPresel_selected.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );//#####
      if( !ripet_Presel ){ nEvent_Presel++; ripet_Presel = true;}

    if( jet1.Pt()>ptJet1_thresh_ ){
      if( !ripet_LeadJetPt ){ nEvent_LeadJetPt++; ripet_LeadJetPt = true;}
      
      if( jet2.Pt()> ptJet2_thresh_ ){
	if( !ripet_SubleadJetPt ){  nEvent_SubleadJetPt++; ripet_SubleadJetPt = true;}
	
	if(  fabs(jet1.Eta())<etaJet1_thresh_ && fabs(jet2.Eta())<etaJet1_thresh_  ){
	  if( !ripet_EtaJet ){ nEvent_EtaJet++;  ripet_EtaJet = true;}
	  
	  if( jet1.DeltaR(jet2) < deltaRjj_thresh_ ){
	    if( !ripet_DrJetJet ){  nEvent_DrJetJet++; ripet_DrJetJet = true;}
	    
	    if(  diJet.M() > mWjj_threshLo_ && diJet.M() < mWjj_threshHi_  ){ 
	      if( !ripet_DijetMass ){  nEvent_DijetMass++; ripet_DijetMass = true;}
	      
	      if(  diJet.Pt() > ptWjj_thresh_){
	      if( !ripet_DijetPt ){ nEvent_DijetPt++; ripet_DijetPt = true;}
		jetPairs_selected.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );
	      }}}}}}
      h1_mWjj_all_presel->Fill( diJet.M(), eventWeight );
      h1_ptWjj_all_presel->Fill( diJet.Pt(), eventWeight );
      h1_deltaRjj_all_presel->Fill(jet1.DeltaR(jet2), eventWeight );

      if( jet1.Pt()!=cached_jetpt ) {
        h1_ptJet_all_presel->Fill( jet1.Pt(), eventWeight );
        h1_etaJet_all_presel->Fill( jet1.Eta(), eventWeight );

        int jetPtBin = -1;
        for( unsigned iPtBin=0; iPtBin<nPtBins-1; ++iPtBin ) {
          if( jet1.Pt()>ptBins[iPtBin] && jet1.Pt()<ptBins[iPtBin+1] ) {
            jetPtBin=iPtBin;
            break;
          }
          if( iPtBin==nPtBins-2 ) jetPtBin=iPtBin; //jets with pt higher that max pt are put in last bin
        } 
        if( jetPtBin<0 ) {
          std::cout << "there must be an error this is not possible." << std::endl;
          exit(9187);
        }

        vh1_rmsCandJet_all_presel[jetPtBin]->Fill( rmsCandJet1[iJetPair], eventWeight );
        vh1_ptDJet_all_presel[jetPtBin]->Fill( ptDJet1[iJetPair], eventWeight );
        vh1_nChargedJet_all_presel[jetPtBin]->Fill( nChargedJet1[iJetPair], eventWeight );
        vh1_nNeutralJet_all_presel[jetPtBin]->Fill( nNeutralJet1[iJetPair], eventWeight );
        
        cached_jetpt = jet1.Pt();
      }

    }//for iJetPairs
      //### Scelgo i due jet migliori senza pero' tagli, ma solo preselezione.
    if( jetPairs_JustPresel_selected.size()>0. ) {
      // now look for best W mass jet pair:
      float Wmass = 80.399;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_JustPresel_selected.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_JustPresel_selected[iPair].first + jetPairs_JustPresel_selected[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Wmass) < fabs(bestMass-Wmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs
      AnalysisJet jet1_JustPresel = jetPairs_JustPresel_selected[bestPair].first;
      AnalysisJet jet2_JustPresel = jetPairs_JustPresel_selected[bestPair].second;

      TLorentzVector bestWDiJet_JustPresel = jet1_JustPresel + jet2_JustPresel;
      h1_deltaRjj_JustPresel->Fill(jet1_JustPresel.DeltaR(jet2_JustPresel), eventWeight );
      h1_mWjj_JustPresel->Fill(bestWDiJet_JustPresel.M(), eventWeight );
	  h1_ptJet1_JustPresel->Fill(jet1_JustPresel.Pt(), eventWeight );
	  h1_ptJet2_JustPresel->Fill(jet2_JustPresel.Pt(), eventWeight );
    }
 //###


    TLorentzVector jet1_presel, jet2_presel;

   //if( jetPairs_selected.size() > 1 ) continue; // Veto su altri Jet
    if( jetPairs_selected.size()>0. ) {

      // now look for best W mass jet pair:
      float Wmass = 80.399;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_selected.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_selected[iPair].first + jetPairs_selected[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Wmass) < fabs(bestMass-Wmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs

   
      AnalysisJet jet1 = jetPairs_selected[bestPair].first;
      AnalysisJet jet2 = jetPairs_selected[bestPair].second;
      TLorentzVector bestWDiJet = jet1 + jet2;
      
      // LEPTONS

    TLorentzVector lept1, lept2;
    lept1.SetPtEtaPhiE( ptLept, etaLept, phiLept, eLept );

    // Find Neutrino //@@
    TLorentzVector neuR, neuW;
    TLorentzVector NeuMC, WllMC, HiggsMC;
    NeuMC.SetPtEtaPhiE(ptNeuMC,etaNeuMC,phiNeuMC,eNeuMC);    
    WllMC.SetPtEtaPhiE(ptWllMC,etaWllMC,phiWllMC,eWllMC); 
    HiggsMC.SetPtEtaPhiE(ptHiggsMC,etaHiggsMC,phiHiggsMC,eHiggsMC); 
    
    float pn=0., pnR=0.,pnW=0., pnMH=0., pnATLAS=0.;
    //pnR=getPzRight(lept1, pxPFMet, pyPFMet, NeuMC);
    // pnW=getPzWrong(lept1, pxPFMet, pyPFMet, NeuMC);
    neuW.SetPxPyPzE(pxPFMet, pyPFMet, pnW, sqrt(pow( pxPFMet,2)+pow(pyPFMet,2)+pow(pnW,2)) );
    neuR.SetPxPyPzE(pxPFMet, pyPFMet, pnR, sqrt(pow( pxPFMet,2)+pow(pyPFMet,2)+pow(pnR,2)) );

    pn=getPz(lept1, pxPFMet, pyPFMet, jet1, jet2);
    pnMH=getPzMH(lept1, pxPFMet, pyPFMet, jet1, jet2);
    pnATLAS=getPzATLAS(lept1, pxPFMet, pyPFMet, jet1, jet2);
    lept2.SetPxPyPzE( pxPFMet, pyPFMet, pn,  sqrt(pow( pxPFMet,2)+pow(pyPFMet,2)+pow(pn,2)) );
    if( energyPFMet < 30. ) continue;

    h1_energyMet->Fill( energyPFMet );

     h1_FindPz_DRR->Fill(  lept1.Eta()-neuR.Eta() );
     h1_FindPz_DRW->Fill(  lept1.Eta()-neuW.Eta() );
     h1_FindPz_EtaR->Fill(lept1.Eta()-neuR.Eta());     h1_FindPz_EtaW->Fill(lept1.Eta()-neuW.Eta());
     h1_FindPz_EtaWnR->Fill((neuR+lept1).Eta()-neuR.Eta());  h1_FindPz_EtaWnW->Fill((neuW+lept1).Eta()-neuW.Eta()); 
     h1_FindPz_PzR->Fill(neuR.Pz());                                   h1_FindPz_PzW->Fill(neuW.Pz());

     h1_FindPz_DeltaPz->Fill(fabs(neuR.Pz()-neuW.Pz()));
     h1_FindPz_DeltaEta->Fill( fabs(neuR.Eta()-lept1.Eta())-fabs(neuW.Eta()-lept1.Eta()) );
     h1_FindPz_DeltaEtaWn->Fill(  fabs(neuR.Eta()-(lept1+neuR).Eta())-fabs(neuW.Eta()-(lept1+neuW).Eta()) );
     TLorentzVector neuR_Wstar(neuR);
     neuR_Wstar.Boost(-(neuR+lept1).BoostVector());
     h1_FindPz_WNeutWR->Fill(cos( neuR_Wstar.Angle((neuR+lept1).Vect())) );
     TLorentzVector neuW_Wstar(neuW);
     neuW_Wstar.Boost(-(neuW+lept1).BoostVector());
     h1_FindPz_WNeutWW->Fill( cos(neuW_Wstar.Angle((neuW+lept1).Vect())) );

     // Resolution on neutrino Pz
     h1_resoPz->Fill( (lept2.Pz()-NeuMC.Pz())/NeuMC.Pz() );
     h1_resoPt->Fill( (energyPFMet-NeuMC.Pt())/NeuMC.Pt() );

     h1_resoPzMH->Fill( (pnMH-NeuMC.Pz())/NeuMC.Pz() );
     h1_resoPzATLAS->Fill( (pnATLAS-NeuMC.Pz())/NeuMC.Pz() );
     h1_resoPzRight->Fill( (pnR-NeuMC.Pz())/NeuMC.Pz() );
     h1_resoPzWrong->Fill( (pnW-NeuMC.Pz())/NeuMC.Pz() );
     // Resolution with Wll
     TLorentzVector neuMH, neuATLAS;
     neuMH.SetPxPyPzE( pxPFMet, pyPFMet, pnMH, sqrt(pow( pxPFMet,2)+pow(pyPFMet,2)+pow(pnMH,2)) );
     neuATLAS.SetPxPyPzE( pxPFMet, pyPFMet, pnATLAS, sqrt(pow( pxPFMet,2)+pow(pyPFMet,2)+pow(pnATLAS,2)));
     h1_resoPzW->Fill( ((lept1+lept2).Pz()-WllMC.Pz())/WllMC.Pz() );
     h1_resoPzWMH->Fill( ((lept1+neuMH).Pz()-WllMC.Pz())/WllMC.Pz() );
     h1_resoPzWATLAS->Fill( ((lept1+neuATLAS).Pz()-WllMC.Pz())/WllMC.Pz() );
     h1_resoPtW->Fill( ((lept1+lept2).Pt()-WllMC.Pt())/WllMC.Pt() );
     h1_resomH->Fill( ((lept1+neuMH+jet1+jet2).M()-HiggsMC.M())/HiggsMC.M() );
     h1_resomHATLAS->Fill( ((lept1+neuATLAS+jet1+jet2).M()-HiggsMC.M())/HiggsMC.M() );

     // Sort leptons and them charges
     float chargeLept1=0., chargeLept2=0.;
     chargeLept1=chargeLept;
     chargeLept2=0.;
     TLorentzVector SortLeptons;
     if( lept1.Pt() < lept2.Pt() ){
       SortLeptons = lept1;
       lept1 = lept2;
       lept2 = SortLeptons;
       chargeLept1=0.;
       chargeLept2=chargeLept;
     }
      
     TLorentzVector diLepton = lept1+lept2;
     
     if( nPairs>0 ) {
       
       h1_mtW_JustPresel->Fill(sqrt(2*lept1.Pt()*lept2.Pt()*(1-cos(delta_phi(lept1.Phi(),lept2.Phi())))));
       h1_mWll_presel->Fill( diLepton.M(), eventWeight );
       if( leptType==0 )  h1_mWmumu_presel->Fill( diLepton.M(), eventWeight );
       else  h1_mWee_presel->Fill( diLepton.M(), eventWeight );

      h1_ptWll_JustPresel->Fill( diLepton.Pt(), eventWeight );

      h1_deltaRll_JustPresel->Fill(lept1.DeltaR(lept2), eventWeight );
      
      h1_ptLept1_JustPresel->Fill( lept1.Pt(), eventWeight );
      h1_ptLept2_JustPresel->Fill( lept2.Pt(), eventWeight );
      
      h1_etaLept1_presel->Fill( lept1.Eta(), eventWeight );
      h1_etaLept2_presel->Fill( lept2.Eta(), eventWeight );
      
      int nJets = getNJets(nPairs);
      h1_nJets_presel->Fill( nJets , eventWeight );
      h1_nPairs_presel->Fill( nPairs , eventWeight );

    } else {

      h1_mWll_presel_0jets->Fill( diLepton.M(), eventWeight );
      if( leptType==0 )
        h1_mWmumu_presel_0jets->Fill( diLepton.M(), eventWeight );
      else
        h1_mWee_presel_0jets->Fill( diLepton.M(), eventWeight );

    }
    
     if( jet1.trackCountingHighEffBJetTag < 1.7 && jet2.trackCountingHighEffBJetTag < 1.7 ){ //btag
      if( diLepton.Pt() > ptWll_thresh_ ){
	nEvent_DileptPt++;
	if(  lept1.Pt() > ptLept1_thresh_ ){
	  nEvent_LeadLeptPt++;
	  if(  lept2.Pt() > ptLept2_thresh_ ){
	    nEvent_SubleadLeptPt++;
	    if( fabs(lept1.Eta()) < etaLept1_thresh_ && fabs(lept2.Eta()) < etaLept2_thresh_){
	      nEvent_EtaLept++;
	      if( diLepton.M() > mtWll_threshLo_ && diLepton.M() < mtWll_threshHi_ ){
		nEvent_mtW++;
		if( lept1.DeltaR(lept2) < deltaRll_thresh_ ){
		  nEvent_DrLeptLept++;
       
		  // event has passed kinematic selection
		  
		  //  M(Higgs) Di Prova.		  
		  h1_mWW_GetPz->Fill( (bestWDiJet+ diLepton).M() );

		  /*	 
		  // ------------------------
		  //   KINEMATIC FIT (LEPTON) : BEGIN
		  // ------------------------
		  TLorentzVector OrdLept1, OrdLept2; //I take Eta from kinematic fit
		  OrdLept1.SetPtEtaPhiE(ptLept1Ord,etaLept1Ord,phiLept1Ord,eLept1Ord);
		  OrdLept2.SetPtEtaPhiE(ptLept2Ord,etaLept2Ord,phiLept2Ord,eLept2Ord);
		  
		  MissingEnergy neut; // Class defined in LeptonNeutrinoKinfitter.h
		  neut.SetuncorrEnergyAK5Jet(SumEt);
		  neut.SetNeutrino(OrdLept2);

		  LeptonNeutrinoKinFitter* fitter_lept = new LeptonNeutrinoKinFitter( "fitter_lept", "fitter_lept", Wmass );
		  std::pair<TLorentzVector,TLorentzVector> lept_kinfit = fitter_lept->fit(OrdLept1, neut);
		  TLorentzVector lept1_kinfit(lept_kinfit.first);
		  TLorentzVector lept2_kinfit(lept_kinfit.second);
		  //std::cout<<pzNeutrino<<" "<<OrdLept2.Pz()<<" "<<lept2_kinfit.Pz()<<std::endl;
		  //std::cout<<(lept1_kinfit+lept2_kinfit).M()<<endl;
		  //std::cout<<"STOP"<<std::endl;

		  h1_pzResoNeut_GetPz->Fill( (isMC) ? (OrdLept2.Pz()-pzNeutrino)/pzNeutrino : 0 );
		  h1_pzResoNeut_KinFit->Fill( (isMC) ? (lept2_kinfit.Pz()-pzNeutrino)/pzNeutrino : 0 );
		  
		  // Sort lept_kinfit
		  TLorentzVector SortLeptons;
		  if( lept1_kinfit.Pt() < lept2_kinfit.Pt() ){
		    SortLeptons = lept1_kinfit;
		    lept1_kinfit = lept2_kinfit;
		    lept2_kinfit = SortLeptons;
		  }

		  h1_mWW_kinLept->Fill( (bestWDiJet+ lept1_kinfit+ lept2_kinfit).M() );
	      */
     		     
		  // ------------------------
		  //   KINEMATIC FIT (JET) : BEGIN
		  // ------------------------
		  h1_mWW_nokinfit->Fill( (bestWDiJet+lept1+lept2).M() );
		  DiJetKinFitter* fitter_jets = new DiJetKinFitter( "fitter_jets", "fitter_jets", Wmass );
		  std::pair<TLorentzVector,TLorentzVector> jets_kinfit = fitter_jets->fit(jet1, jet2);
		  TLorentzVector jet1_kinfit(jets_kinfit.first);
		  TLorentzVector jet2_kinfit(jets_kinfit.second);
		  
		  TLorentzVector matchedPart1, matchedPart2;
		  float bestDeltaRPart1=999.;
		  float bestDeltaRPart2=999.;
		  for( unsigned iPart=0; iPart<nPart; ++iPart ) {
		    if( abs(pdgIdPart[iPart])>6 ) continue;
		    TLorentzVector thisPart;
          thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
          if( jet1.DeltaR(thisPart) < bestDeltaRPart1 ) {
            bestDeltaRPart1 = jet1.DeltaR(thisPart);
            matchedPart1 = thisPart;
          }
          if( jet2.DeltaR(thisPart) < bestDeltaRPart2 ) {
            bestDeltaRPart2 = jet2.DeltaR(thisPart);
            matchedPart2 = thisPart;
          }
        }
		  //float ptReso1_before = (isMC) ? ( jet1.Pt()-jet1.ptGen )/jet1.ptGen : 0.;
		  //float ptReso2_before = (isMC) ? ( jet2.Pt()-jet2.ptGen )/jet2.ptGen : 0.;
		  float ptReso1_before = (isMC) ? ( jet1.Pt()-matchedPart1.Pt() )/matchedPart1.Pt() : 0.;
		  float ptReso2_before = (isMC) ? ( jet2.Pt()-matchedPart2.Pt() )/matchedPart2.Pt() : 0.;
		  h1_ptResoJet1_beforeKin->Fill( ptReso1_before, eventWeight );
		  h1_ptResoJet2_beforeKin->Fill( ptReso2_before, eventWeight );
		  
		  float ptReso1_after = (isMC) ? ( jet1_kinfit.Pt()-matchedPart1.Pt() )/matchedPart1.Pt() : 0.;
		  float ptReso2_after = (isMC) ? ( jet2_kinfit.Pt()-matchedPart2.Pt() )/matchedPart2.Pt() : 0.;
		  h1_ptResoJet1_afterKin->Fill( ptReso1_after, eventWeight );
		  h1_ptResoJet2_afterKin->Fill( ptReso2_after, eventWeight );

	
		  TLorentzVector Wjj_kinfit = jet1_kinfit + jet2_kinfit;
		  
		  TLorentzVector matchedW;
		  float bestDeltaRW=999.;
        for( unsigned iPart=0; iPart<nPart; ++iPart ) {
          if( pdgIdPart[iPart]!=24 ) continue;
          TLorentzVector thisW;
          thisW.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
          if( Wjj_kinfit.DeltaR(thisW) < bestDeltaRW ) {
            bestDeltaRW = Wjj_kinfit.DeltaR(thisW);
            matchedW = thisW;
          }
        }

        float ptWreso_before = (isMC) ? ( bestWDiJet.Pt()-matchedW.Pt() )/matchedW.Pt() : 0.;
        float ptWreso_after  = (isMC) ? ( Wjj_kinfit.Pt()-matchedW.Pt() )/matchedW.Pt() : 0.;
        h1_ptWreso_beforeKin->Fill( ptWreso_before, eventWeight);
        h1_ptWreso_afterKin->Fill( ptWreso_after, eventWeight);

        TLorentzVector Wjj_constr;
        Wjj_constr.SetXYZM( bestWDiJet.Px(), bestWDiJet.Py(), bestWDiJet.Pz(), Wmass);

        TLorentzVector WW = bestWDiJet + diLepton;
        TLorentzVector WW_constr = diLepton + Wjj_constr;
        TLorentzVector WW_kinfit = diLepton + Wjj_kinfit;

        float chiSquareProb = TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF());
        h1_kinfit_chiSquare->Fill( fitter_jets->getS()/fitter_jets->getNDF(), eventWeight ); 
        h1_kinfit_chiSquareProb->Fill( chiSquareProb, eventWeight ); 

nEventsTot += eventWeight;
if( chiSquareProb>0.01 ) {
  nEvents_hiChiSquareProb += eventWeight;
  h1_mWjj_hiChiSquareProb->Fill( bestWDiJet.M(), eventWeight );
  h1_mWW_hiChiSquareProb->Fill( WW_kinfit.M(), eventWeight );
} else {
  h1_mWjj_loChiSquareProb->Fill( bestWDiJet.M(), eventWeight );
  h1_mWW_loChiSquareProb->Fill( WW_kinfit.M(), eventWeight );
}
if( bestWDiJet.M()>75. && bestWDiJet.M()<105. ) {
  nEvents_mWjj_cut += eventWeight;
  h1_mWW_mWjj_cut->Fill( WW_kinfit.M(), eventWeight );
} else {
  h1_mWW_mWjj_notcut->Fill( WW_kinfit.M(), eventWeight );
}

if( TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF())<0.2 ) h1_mWW_kinfit->Fill( WW_kinfit.M(), eventWeight );

		 
    //// match to parton:
    //int partFlavor1=0;
    //float deltaRmin1=999.;
    //for(unsigned iPart=0; iPart<nPart; ++iPart ) {
    //  TLorentzVector thisPart;
    //  thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
    //  float thisDeltaR = jet1_presel.DeltaR(thisPart);
    //  if( thisDeltaR<deltaRmin1 ) {
    //    partFlavor1 = pdgIdPart[iPart];
    //    deltaRmin1 = thisDeltaR;
    //  }
    //}

    //float deltaRmin2=999.;
    //int partFlavor2=0;
    //for(unsigned iPart=0; iPart<nPart; ++iPart ) {
    //  TLorentzVector thisPart;
    //  thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
    //  float thisDeltaR = jet2_presel.DeltaR(thisPart);
    //  if( thisDeltaR<deltaRmin2 ) {
    //    partFlavor2 = pdgIdPart[iPart];
    //    deltaRmin2 = thisDeltaR;
    //  }
    //}

    //bool bothMatched = ( deltaRmin1<0.5 && deltaRmin2<0.5 && partFlavor1!=0 && partFlavor2!=0 );

    //if( bothMatched ) {
    //  h1_mWW_MCassoc->Fill( WW.M(), eventWeight );
    //  h1_mWW_MCassoc_WjjMassConstr->Fill( WW_constr.M(), eventWeight );
    //  h1_mWW_MCassoc_kinfit->Fill( WW_kinfit.M(), eventWeight );
    //}

//    // and now full kinematic fit with PFCands:
//    TLorentzVector WW_kinfit_cands;
//    TRegexp cands_tstr("CANDS");
//    TString dataset_tstr(dataset_);
//if( dataset_tstr.Contains(cands_tstr) ) {

//    std::vector<TFitParticleEtEtaPhi*> fitCands;
//    std::vector<int> fitCandTypes;
//    TFitConstraintM *mCons_cands = new TFitConstraintM( "WMassConstraint_cands", "WMass-Constraint", 0, 0 , 91.19);


//    float testFactor = 1.;

//    // loop on PFCands of first jet
//    TLorentzVector candJet1(0., 0., 0., 0.);
//    for( unsigned iPFCand=0; iPFCand<nPFCand1; ++iPFCand ) {
//      TLorentzVector thisCand;
//      thisCand.SetPtEtaPhiE( ptPFCand1[iPFCand], etaPFCand1[iPFCand], phiPFCand1[iPFCand], ePFCand1[iPFCand] );
//      char candName[100];
//      sprintf( candName, "PFCand_1_%d", iPFCand );

//      if( particleTypePFCand1[iPFCand]==0 ) {
//        std::cout << "FOUND PARTICLE TYPE=0!!!! SKIPPING!" << std::endl;
//        continue;
//      }

//      if( particleTypePFCand1[iPFCand]==5 ) // neutral hadrons
//        testFactor=10.;
//      else
//        testFactor=1.;

//      TMatrixD m_PFCand(3,3);
//      m_PFCand(0,0) = testFactor*ErrEt ( thisCand.Et(), thisCand.Eta(), particleTypePFCand1[iPFCand] );
//      m_PFCand(1,1) = testFactor*ErrEta( thisCand.Et(), thisCand.Eta(), particleTypePFCand1[iPFCand] );
//      m_PFCand(2,2) = testFactor*ErrPhi( thisCand.Et(), thisCand.Eta(), particleTypePFCand1[iPFCand] );
//      TFitParticleEtEtaPhi* fitCand = new TFitParticleEtEtaPhi( candName, candName, &thisCand, &m_PFCand );

//      mCons_cands->addParticle1( fitCand );
//      fitCands.push_back(fitCand);
//      fitCandTypes.push_back(particleTypePFCand1[iPFCand]);

//      candJet1 += thisCand;

//    }

//    TLorentzVector cand_add1 = jet1_presel - candJet1;
//    // make it worse-resolution case:
//    int particleType_add1 = (fabs(cand_add1.Eta())<3.) ? 5 : 6; 
//    TMatrixD m_PFCand_add1(3,3);
//    m_PFCand_add1(0,0) = testFactor*ErrEt ( cand_add1.Et(), cand_add1.Eta(), particleType_add1 );
//    m_PFCand_add1(1,1) = testFactor*ErrEta( cand_add1.Et(), cand_add1.Eta(), particleType_add1 );
//    m_PFCand_add1(2,2) = testFactor*ErrPhi( cand_add1.Et(), cand_add1.Eta(), particleType_add1 );
//    TFitParticleEtEtaPhi* fitCand_add1 = new TFitParticleEtEtaPhi( "PFCand_add1", "PFCand_add1", &cand_add1, &m_PFCand_add1 );

//  //mCons_cands->addParticle1( fitCand_add1 );
//  //fitCands.push_back(fitCand_add1);
//  //fitCandTypes.push_back(particleType_add1 );


//    // loop on PFCands of second jet
//    TLorentzVector candJet2(0., 0., 0., 0.);
//    for( unsigned iPFCand=0; iPFCand<nPFCand2; ++iPFCand ) {
//      TLorentzVector thisCand;
//      thisCand.SetPtEtaPhiE( ptPFCand2[iPFCand], etaPFCand2[iPFCand], phiPFCand2[iPFCand], ePFCand2[iPFCand] );
//      char candName[100];
//      sprintf( candName, "PFCand_2_%d", iPFCand );

//      if( particleTypePFCand2[iPFCand]==0 ) {
//        std::cout << "FOUND PARTICLE TYPE=0!!!! SKIPPING!" << std::endl;
//        continue;
//      }

//      if( particleTypePFCand1[iPFCand]==5 ) // neutral hadrons
//        testFactor=10.;
//      else
//        testFactor=1.;

//      TMatrixD m_PFCand(3,3);
//      m_PFCand(0,0) = testFactor*ErrEt ( thisCand.Et(), thisCand.Eta(), particleTypePFCand2[iPFCand] );
//      m_PFCand(1,1) = testFactor*ErrEta( thisCand.Et(), thisCand.Eta(), particleTypePFCand2[iPFCand] );
//      m_PFCand(2,2) = testFactor*ErrPhi( thisCand.Et(), thisCand.Eta(), particleTypePFCand2[iPFCand] );
//      TFitParticleEtEtaPhi* fitCand = new TFitParticleEtEtaPhi( candName, candName, &thisCand, &m_PFCand );

//      mCons_cands->addParticle1( fitCand );
//      fitCands.push_back(fitCand);
//      fitCandTypes.push_back(particleTypePFCand2[iPFCand]);

//      candJet2 += thisCand;

//    }


//    TLorentzVector cand_add2 = jet2_presel - candJet2;
//    // make it worse-resolution case:
//    int particleType_add2 = (fabs(cand_add2.Eta())<3.) ? 5 : 6; 
//    TMatrixD m_PFCand_add2(3,3);
//    m_PFCand_add2(0,0) = testFactor*ErrEt ( cand_add2.Et(), cand_add2.Eta(), particleType_add2 );
//    m_PFCand_add2(1,1) = testFactor*ErrEta( cand_add2.Et(), cand_add2.Eta(), particleType_add2 );
//    m_PFCand_add2(2,2) = testFactor*ErrPhi( cand_add2.Et(), cand_add2.Eta(), particleType_add2 );
//    TFitParticleEtEtaPhi* fitCand_add2 = new TFitParticleEtEtaPhi( "PFCand_add2", "PFCand_add2", &cand_add2, &m_PFCand_add2 );

//  //mCons_cands->addParticle2( fitCand_add2 );
//  //fitCands.push_back(fitCand_add2);
//  //fitCandTypes.push_back(particleType_add2 );

//    TKinFitter* fitter_cands = new TKinFitter("fitter_cands", "fitter_cands");
//    for( unsigned iCand=0; iCand<fitCands.size(); ++iCand )
//      fitter_cands->addMeasParticle( fitCands[iCand] );
//    fitter_cands->addConstraint( mCons_cands );

//    //Set convergence criteria
//    fitter_cands->setMaxNbIter( 30 );
//    fitter_cands->setMaxDeltaS( 1e-2 );
//    fitter_cands->setMaxF( 1e-1 );
//    fitter_cands->setVerbosity(0);

//    //Perform the fit
//    fitter_cands->fit();

//    // recreate Wqq vector:
//    TLorentzVector Wqq_kinfit_cands;
//    for( unsigned iCand=0; iCand<fitCands.size(); ++iCand )
//      Wqq_kinfit_cands += *(fitCands[iCand]->getCurr4Vec());


//    for( unsigned iCand=0; iCand<fitCands.size(); ++iCand ) {
//      if( fitCands[iCand]->getCurr4Vec()->Pt()>0. && fitCands[iCand]->getIni4Vec()->Pt()>0. ) {
//        float deltaE   = (fitCands[iCand]->getCurr4Vec()->Energy() - fitCands[iCand]->getIni4Vec()->Energy())/fitCands[iCand]->getIni4Vec()->Energy();
//        float deltaEta = fitCands[iCand]->getCurr4Vec()->Eta()    - fitCands[iCand]->getIni4Vec()->Eta();
//        float deltaPhi = fitCands[iCand]->getCurr4Vec()->DeltaPhi(fitCands[iCand]->getIni4Vec()->Phi());
//        float deltaPt  = (fitCands[iCand]->getCurr4Vec()->Pt() - fitCands[iCand]->getIni4Vec()->Pt())/fitCands[iCand]->getIni4Vec()->Pt();
//        if( fitCandTypes[iCand]==1 ) { //charged hadrons
//          h1_deltaE_ch->Fill( deltaE ); 
//          h1_deltaEta_ch->Fill( deltaEta ); 
//          h1_deltaPhi_ch->Fill( deltaPhi ); 
//          h1_deltaPt_ch->Fill( deltaPt ); 
//        } else if( fitCandTypes[iCand]==4 ) { //photons
//          h1_deltaE_gamma->Fill( deltaE ); 
//          h1_deltaEta_gamma->Fill( deltaEta ); 
//          h1_deltaPhi_gamma->Fill( deltaPhi ); 
//          h1_deltaPt_gamma->Fill( deltaPt ); 
//        } else if( fitCandTypes[iCand]==5 ) { //neutral hadrons
//          h1_deltaE_nh->Fill( deltaE ); 
//          h1_deltaEta_nh->Fill( deltaEta ); 
//          h1_deltaPhi_nh->Fill( deltaPhi ); 
//          h1_deltaPt_nh->Fill( deltaPt ); 
//        }
//      }
//    } //for cands
//  

//    WW_kinfit_cands = Wll + Wqq_kinfit_cands;

//    if( bothMatched )
//      h1_mWW_MCassoc_kinfit_cands->Fill( WW_kinfit_cands.M(), eventWeight );
//    

//  } // if dataset


      // ------------------------
      //   KINEMATIC FIT: END
      // ------------------------


 		 

/*
lept1.SetPtEtaPhiE(21.71, 1.16, -0.79, 37.96);
lept2.SetPtEtaPhiE(154.6, 1.86, 0.73, 509.0);
jet1.SetPtEtaPhiE(40.5, 2.17, -2.2, 180.4);
jet2.SetPtEtaPhiE(36.2, 0.64, -0.57, 44.47);
*/


      //get helicity angles:

 HelicityLikelihoodDiscriminant::HelicityAngles hangles;
 if( chargeLept1<chargeLept2 ) hangles = computeHelicityAngles(lept1, lept2, jet1, jet2);
 else                hangles = computeHelicityAngles(lept2, lept1, jet1, jet2);
 
 HelicityLikelihoodDiscriminant::HelicityAngles hangles_kinfit;
 if( chargeLept1<chargeLept2 ) hangles_kinfit = computeHelicityAngles(lept1, lept2, jet1_kinfit, jet2_kinfit);
 else                hangles_kinfit = computeHelicityAngles(lept2, lept1, jet1_kinfit, jet2_kinfit);
 
 
      HelicityLikelihoodDiscriminant *LD = new HelicityLikelihoodDiscriminant();
      
           
      LD->setMeasurables(hangles);
      double sProb=LD->getSignalProbability();
      double bProb=LD->getBkgdProbability();
      double helicityLD=sProb/(sProb+bProb);
      h1_helicityLD->Fill(helicityLD, eventWeight);


      LD->setMeasurables(hangles_kinfit);
      double sProb_kinfit=LD->getSignalProbability();
      double bProb_kinfit=LD->getBkgdProbability();
      double helicityLD_kinfit=sProb_kinfit/(sProb_kinfit+bProb_kinfit);
      h1_helicityLD_kinfit->Fill(helicityLD_kinfit, eventWeight);
      
      
      //
      // QG LIKELIHOOD   ***BEGIN***
      //
      
      int jet1PtBin=-1;
      if( jet1.Pt() > ptBins[nPtBins] ) {
        jet1PtBin = nPtBins-1;
      } else {
        for( unsigned int iBin=0; iBin<nPtBins; ++iBin ) {
          if( jet1.Pt()>ptBins[iBin] && jet1.Pt()<ptBins[iBin+1] ) {
            jet1PtBin = iBin;
            break;
          }
        }
      }

      int jet2PtBin=-1;
      if( jet2.Pt() > ptBins[nPtBins] ) {
        jet2PtBin = nPtBins-1;
      } else {
        for( unsigned int iBin=0; iBin<nPtBins; ++iBin ) {
          if( jet2.Pt()>ptBins[iBin] && jet2.Pt()<ptBins[iBin+1] ) {
            jet2PtBin = iBin;
            break;
          }
        }
      }
      if( jet1PtBin<0 ){ std::cout<<"Jet1 not reconstructed"<<std::endl;  break;}
      if( jet2PtBin<0 ){ std::cout<<"Jet2 not reconstructed"<<std::endl;  break;} // se e' meno uno vuol dire che non ha trovato il jet

      //   bool btag_TChighPur = ( jet1.trackCountingHighPurBJetTag>5. || jet2.trackCountingHighPurBJetTag>5. );
      //  bool btag_TChighEff = ( jet1.trackCountingHighEffBJetTag>5. || jet2.trackCountingHighEffBJetTag>5. );
      //  bool btag_SSVhighPur = ( jet1.simpleSecondaryVertexHighPurBJetTag>2. || jet2.simpleSecondaryVertexHighPurBJetTag>2. );
      // bool btag_SSVhighEff = ( jet1.simpleSecondaryVertexHighEffBJetTag>2. || jet2.simpleSecondaryVertexHighEffBJetTag>2. );
      //  bool btag_SSVhighPurhighEff = ( ( jet1.simpleSecondaryVertexHighPurBJetTag>2. && jet2.simpleSecondaryVertexHighEffBJetTag>2. ) ||
      //                             ( jet1.simpleSecondaryVertexHighEffBJetTag>2. && jet2.simpleSecondaryVertexHighPurBJetTag>2. ) );

      h1_simpleSecondaryVertexHighEffBJetTagJet1->Fill(jet1.simpleSecondaryVertexHighEffBJetTag, eventWeight);
      h1_simpleSecondaryVertexHighPurBJetTagJet1->Fill(jet1.simpleSecondaryVertexHighPurBJetTag, eventWeight);
      h1_jetBProbabilityBJetTagJet1->Fill(jet1.jetBProbabilityBJetTag, eventWeight);
      h1_jetProbabilityBJetTagJet1->Fill(jet1.jetProbabilityBJetTag, eventWeight);

      h1_simpleSecondaryVertexHighEffBJetTagJet2->Fill(jet2.simpleSecondaryVertexHighEffBJetTag, eventWeight);
      h1_simpleSecondaryVertexHighPurBJetTagJet2->Fill(jet2.simpleSecondaryVertexHighPurBJetTag, eventWeight);
      h1_jetBProbabilityBJetTagJet2->Fill(jet2.jetBProbabilityBJetTag, eventWeight);
      h1_jetProbabilityBJetTagJet2->Fill(jet2.jetProbabilityBJetTag, eventWeight);



      h1_ptDJet1->Fill( jet1.ptD, eventWeight );
      h1_ptDJet2->Fill( jet2.ptD, eventWeight );

      //vh1_rmsCandJet1[jet1PtBin]->Fill( jet1.rmsCand, eventWeight );
      vh1_ptDJet1[jet1PtBin]->Fill( jet1.ptD, eventWeight );
      vh1_nChargedJet1[jet1PtBin]->Fill( jet1.nCharged, eventWeight );
      vh1_nNeutralJet1[jet1PtBin]->Fill( jet1.nNeutral, eventWeight );
      float QGLikelihoodJet1 = qglikeli->computeQGLikelihood( jet1.Pt(), jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );

      vh1_QGLikelihoodJet1[jet1PtBin]->Fill( QGLikelihoodJet1, eventWeight );
      h1_QGLikelihoodJet1->Fill( QGLikelihoodJet1, eventWeight );
     
    //  vh1_rmsCandJet2[jet2PtBin]->Fill( jet2.rmsCand, eventWeight );
       vh1_ptDJet2[jet2PtBin]->Fill( jet2.ptD, eventWeight );
      vh1_nChargedJet2[jet2PtBin]->Fill( jet2.nCharged, eventWeight ); 
      vh1_nNeutralJet2[jet2PtBin]->Fill( jet2.nNeutral, eventWeight );
      float QGLikelihoodJet2 = qglikeli->computeQGLikelihood( jet2.Pt(), jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );

      vh1_QGLikelihoodJet2[jet2PtBin]->Fill( QGLikelihoodJet2, eventWeight );
      h1_QGLikelihoodJet2->Fill( QGLikelihoodJet2, eventWeight );
 
      //  if( !btag_TChighPur ) {
      //   h1_QGLikelihoodJet1_antiBtag_TChighPur->Fill( QGLikelihoodJet1, eventWeight );
      //   h1_QGLikelihoodJet2_antiBtag_TChighPur->Fill( QGLikelihoodJet2, eventWeight );
      //  }
      //  if( !btag_TChighEff ) {
      //   h1_QGLikelihoodJet1_antiBtag_TChighEff->Fill( QGLikelihoodJet1, eventWeight );
      //  h1_QGLikelihoodJet2_antiBtag_TChighEff->Fill( QGLikelihoodJet2, eventWeight );
      //  }
      //  if( !btag_SSVhighPur ) {
      //  h1_QGLikelihoodJet1_antiBtag_SSVhighPur->Fill( QGLikelihoodJet1, eventWeight );
      //  h1_QGLikelihoodJet2_antiBtag_SSVhighPur->Fill( QGLikelihoodJet2, eventWeight );
      //  }
      //  if( !btag_SSVhighEff ) {
      //   h1_QGLikelihoodJet1_antiBtag_SSVhighEff->Fill( QGLikelihoodJet1, eventWeight );
      //   h1_QGLikelihoodJet2_antiBtag_SSVhighEff->Fill( QGLikelihoodJet2, eventWeight );
      //  }
      //  if( !btag_SSVhighPurhighEff ) {
      //   h1_QGLikelihoodJet1_antiBtag_SSVhighPurhighEff->Fill( QGLikelihoodJet1, eventWeight );
      //   h1_QGLikelihoodJet2_antiBtag_SSVhighPurhighEff->Fill( QGLikelihoodJet2, eventWeight );
      // }


      float QGLikelihoodProd = QGLikelihoodJet1*QGLikelihoodJet2;
      float QGLikelihoodRevProd = (1.-QGLikelihoodJet1)*(1.-QGLikelihoodJet2);

      h1_QGLikelihoodProd->Fill( QGLikelihoodProd, eventWeight );
      //  if( !btag_SSVhighEff ) h1_QGLikelihoodProd_antiBtag->Fill( QGLikelihoodProd, eventWeight );
      
      h1_QGLikelihoodRevProd->Fill( QGLikelihoodRevProd, eventWeight );
        
      h2_QGLikelihoodJet1_vs_Jet2->Fill( QGLikelihoodJet2, QGLikelihoodJet1, eventWeight );

      if( QGLikelihoodJet1>0.8 || QGLikelihoodJet2>0.8 ) h1_QGLikelihoodProd_hi->Fill(QGLikelihoodProd, eventWeight); 

      // last step of selection:
      // QG and helicity LD's

      if( QGLikelihoodProd < QGLikelihoodProd_thresh_ ) continue;
      if( helicityLD < helicityLD_thresh_ ) continue;


      // *****************************************
      // *****  PASSED ANALYSIS SELECTION ********
      // *****************************************
      
      if( WW_kinfit.M() > mWW_threshLo_ && WW_kinfit.M() < mWW_threshHi_ ) {
        nEventsPassed_fb_kinfit += eventWeight;
        nEventsPassed_kinfit++;
	// if( !btag_SSVhighEff ) nEventsPassed_kinfit_antiBtag++;
      }
      if( WW.M() > mWW_threshLo_ && WW.M() < mWW_threshHi_ ) {
        nEventsPassed_fb_nokinfit += eventWeight;
        nEventsPassed_nokinfit++;
      }


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
	h1_mtW->Fill(sqrt(2*lept1.Pt()*lept2.Pt()*(1-cos(delta_phi(lept1.Phi(),lept2.Phi())))));//###
        h1_deltaRjj->Fill( jet1.DeltaR(jet2), eventWeight);
        h1_ptWll->Fill( diLepton.Pt(), eventWeight);
        h1_ptWjj->Fill( bestWDiJet.Pt(), eventWeight);
        if( leptType==0 )
          h1_mWmumu->Fill( diLepton.M(), eventWeight );
        else
          h1_mWee->Fill( diLepton.M(), eventWeight );
        h1_mWll->Fill( diLepton.M(), eventWeight);
        h1_mWjj->Fill( bestWDiJet.M(), eventWeight);
        h1_mWW_UL->Fill(WW.M(), eventWeight);
        h1_mWW_hiMass->Fill(WW.M(), eventWeight);
      //  h1_mWW_highestMass->Fill(WW.M(), eventWeight);
        h1_mWW_medMass->Fill(WW.M(), eventWeight);
        h1_mWW_300Mass->Fill(WW.M(), eventWeight);
        h1_mWW_WjjMassConstr_hiMass->Fill(WW_constr.M(), eventWeight);
        h1_mWW_WjjMassConstr_medMass->Fill(WW_constr.M(), eventWeight);
        h1_mWW_WjjMassConstr_300Mass->Fill(WW_constr.M(), eventWeight);
        h1_mWW_UL_kinfit->Fill(WW_kinfit.M(), eventWeight);
        h1_mWW_kinfit_hiMass->Fill(WW_kinfit.M(), eventWeight);
	if( WW_kinfit.M()>450 && WW_kinfit.M()<550){
	h1_mWW_kinfit_cut->Fill(WW_kinfit.M(), eventWeight);
	}
      //  h1_mWW_kinfit_highestMass->Fill(WW_kinfit.M(), eventWeight);
        h1_mWW_kinfit_medMass->Fill(WW_kinfit.M(), eventWeight);
        h1_mWW_kinfit_300Mass->Fill(WW_kinfit.M(), eventWeight);
        h2_mWjj_vs_mWW->Fill( WW.M(), bestWDiJet.M() );
        h2_mWjj_vs_mWW_kinfit->Fill( WW_kinfit.M(), bestWDiJet.M() );


        h1_deltaRWW->Fill(bestWDiJet.DeltaR(diLepton), eventWeight);

        if( WW_kinfit.M()>390. && WW_kinfit.M()<460. ) {
        h1_ptWW->Fill( WW.Pt(), eventWeight );
        h1_ptWW_kinfit->Fill( WW_kinfit.Pt(), eventWeight );
        h1_etaWW->Fill( WW.Eta(), eventWeight );
        h1_etaWW_kinfit->Fill( WW_kinfit.Eta(), eventWeight );
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
	
		}}}}}  } //if you use btag
      } //if passes selection

    } //if selected jet pairs
      
   
  } //for entries

  //Histo for efficiencies.
  h1_Ev_Presel ->SetBinContent(1,nEvent_Presel);
  h1_Ev_Jet1 ->SetBinContent(1,nEvent_LeadJetPt);
  h1_Ev_Jet2  ->SetBinContent(1,nEvent_SubleadJetPt);
  h1_Ev_EtaJet ->SetBinContent(1,nEvent_EtaJet);
  h1_Ev_DrJet ->SetBinContent(1,nEvent_DrJetJet);
  h1_Ev_mWjj ->SetBinContent(1,nEvent_DijetMass);
  h1_Ev_PtWjj ->SetBinContent(1,nEvent_DijetPt);
  h1_Ev_DileptonPt ->SetBinContent(1,nEvent_DileptPt);
  h1_Ev_Lept1 ->SetBinContent(1,nEvent_LeadLeptPt);
  h1_Ev_Lept2 ->SetBinContent(1,nEvent_SubleadLeptPt);
  h1_Ev_EtaLept ->SetBinContent(1,nEvent_EtaLept);
  h1_Ev_mtW ->SetBinContent(1,nEvent_mtW);
  h1_Ev_DrLept ->SetBinContent(1,nEvent_DrLeptLept);
  h1_Ev_mWW ->SetBinContent(1,nEventsPassed_kinfit);
  h1_Ev_nCounter ->SetBinContent(1,nCounter_);


  h1_Ev_nCounterW ->SetBinContent(1,nCounterW_);
  h1_Ev_nEventsPassed_fb_kinfit ->SetBinContent(1,nEventsPassed_fb_kinfit);
  h1_Ev_nEventsPassed_kinfit ->SetBinContent(1,nEventsPassed_kinfit);
  // h1_Ev_nEventsPassed_kinfit_antiBtag ->SetBinContent(1,nEventsPassed_kinfit_antiBtag);
  h1_Ev_nEventsPassed_fb_nokinfit ->SetBinContent(1,nEventsPassed_fb_nokinfit);
  h1_Ev_nEventsPassed_nokinfit ->SetBinContent(1,nEventsPassed_nokinfit);

  std::cout<<"Ev(Presel)="<<nEvent_Presel<<"Ev(Jet1)="<<nEvent_LeadJetPt<<"Ev(Jet2)="<<nEvent_SubleadJetPt <<" Ev(EtaJet)="<<nEvent_EtaJet<<" Ev(DrJet)=" <<nEvent_DrJetJet<<" Ev(Mwjj)="
	   <<nEvent_DijetMass<<" Ev(PtWjj)="<<nEvent_DijetPt<<std::endl;
  std::cout<<"Ev(DileptonPt)="<<nEvent_DileptPt <<" Ev(Lept1)="<<nEvent_LeadLeptPt<<" Ev(Lept2)=" <<nEvent_SubleadLeptPt<<" Ev(EtaLept)="
	   <<nEvent_EtaLept<<" Ev(mtW)="<<nEvent_mtW<<" Ev(DrLept)="<<nEvent_DrLeptLept<<std::endl;

  std::cout << std::endl << std::endl;
  std::cout << "----> PASSED SELECTION: " << 1000.*nEventsPassed_fb_kinfit << " ev/fb-1(" << nEventsPassed_kinfit << " events)"<< std::endl;
  std::cout << "----> PASSED SELECTION (no kinfit): " << 1000.*nEventsPassed_fb_nokinfit << " ev/fb-1 (" << nEventsPassed_nokinfit << " events)" << std::endl;
  std::cout << std::endl;
  std::cout<<"Eff(Tot)" <<  nEventsPassed_kinfit/nCounter_ <<std::endl;
  std::cout<<"Eff(Tot)" <<  nEventsPassed_fb_kinfit/nCounterW_ <<std::endl;

  outFile_->cd();

  h1_energyMet->Write();
  h1_resoPz->Write();//@@
  h1_resoPt->Write();
  h1_resoPzMH->Write();
  h1_resoPzATLAS->Write();
  h1_resoPzRight->Write();
  h1_resoPzWrong->Write();
  h1_resoPtW->Write();
  h1_resomH->Write();
  h1_resomHATLAS->Write();
  h1_resoPzW->Write();
  h1_resoPzWMH->Write();
  h1_resoPzWATLAS->Write();

  h1_FindPz_EtaR->Write();
  h1_FindPz_EtaW->Write();
  h1_FindPz_DRW->Write();
  h1_FindPz_DRR->Write();
  h1_FindPz_EtaWnW->Write();
  h1_FindPz_EtaWnR->Write();
  h1_FindPz_PzW->Write();
  h1_FindPz_PzR->Write();
  h1_FindPz_DeltaEta->Write();
  h1_FindPz_DeltaPz->Write();
  h1_FindPz_DeltaEtaWn->Write();
  h1_FindPz_WNeutWW->Write();
  h1_FindPz_WNeutWR->Write();

  h1_Ev_Presel ->Write();
  h1_Ev_Jet1 ->Write();
  h1_Ev_Jet2 ->Write();
  h1_Ev_EtaJet ->Write();
  h1_Ev_DrJet ->Write();
  h1_Ev_mWjj ->Write();
  h1_Ev_PtWjj ->Write();
  h1_Ev_DileptonPt ->Write();
  h1_Ev_Lept1 ->Write();
  h1_Ev_Lept2 ->Write();
  h1_Ev_EtaLept ->Write();
  h1_Ev_mtW ->Write();
  h1_Ev_DrLept ->Write();
  h1_Ev_mWW -> Write();
  h1_Ev_nCounter ->Write();
  h1_Ev_nCounterW ->Write();
  h1_Ev_nEventsPassed_fb_kinfit ->Write();
  h1_Ev_nEventsPassed_kinfit ->Write();
//  h1_Ev_nEventsPassed_kinfit_antiBtag ->Write();
  h1_Ev_nEventsPassed_fb_nokinfit ->Write();
  h1_Ev_nEventsPassed_nokinfit ->Write();

  h1_run->Write();

  h1_ptJet_all_presel->Write();
  h1_etaJet_all_presel->Write();
  h1_nJets_presel->Write();
  h1_nPairs_presel->Write();

  h1_ptResoJet1_beforeKin->Write();
  h1_ptResoJet2_beforeKin->Write();
  h1_ptResoJet1_afterKin->Write();
  h1_ptResoJet2_afterKin->Write();

  h1_pzResoNeut_GetPz->Write();
  h1_pzResoNeut_KinFit->Write();

  h1_ptWreso_beforeKin->Write();
  h1_ptWreso_afterKin->Write();

  h1_deltaRll_JustPresel->Write();
  h1_deltaRjj_all_presel->Write();

  h1_ptLept1->Write();
  h1_ptLept2->Write();
  h1_mtW->Write();//###

  h1_ptLept1_JustPresel->Write();
  h1_ptLept2_JustPresel->Write();
  h1_etaLept1_presel->Write();
  h1_etaLept2_presel->Write();
  h1_mtW_JustPresel->Write();//###
  h1_deltaRjj_JustPresel->Write();//####
  h1_ptJet1_JustPresel->Write();
  h1_ptJet2_JustPresel->Write();//###
  h1_mWjj_JustPresel->Write();
  h1_mWll->Write();
  h1_mWll_presel->Write();
  h1_mWmumu->Write();
  h1_mWmumu_presel->Write();
  h1_mWmumu_presel_0jets->Write();
  h1_mWee->Write();
  h1_mWee_presel->Write();
  h1_mWee_presel_0jets->Write();

  h1_mWjj->Write();
  h1_mWjj_loChiSquareProb->Write();
  h1_mWjj_hiChiSquareProb->Write();
  h1_mWjj_all_presel->Write();

  h1_ptWll_JustPresel->Write();
  h1_ptWjj_all_presel->Write();

  h1_deltaRjj->Write();

  h1_ptWjj->Write();
  h1_ptWll->Write();

  h1_cosThetaStar->Write();
  h1_cosTheta1->Write();
  h1_cosTheta2->Write();
  h1_phi->Write();
  h1_phi1->Write();

  h1_cosThetaStar_kinfit->Write();
  h1_cosTheta1_kinfit->Write();
  h1_cosTheta2_kinfit->Write();
  h1_phi_kinfit->Write();
  h1_phi1_kinfit->Write();

  h1_kinfit_chiSquare->Write();
  h1_kinfit_chiSquareProb->Write();
  
  h1_helicityLD->Write();
  h1_helicityLD_kinfit->Write();

  //##
  h1_mWW_GetPz->Write();
  h1_mWW_kinLept->Write();

  h1_mWW_nokinfit->Write();//###
  h1_mWW_kinfit->Write();
  h1_mWW_hiChiSquareProb->Write();
  h1_mWW_loChiSquareProb->Write();
  h1_mWW_mWjj_cut->Write();
  h1_mWW_mWjj_notcut->Write();
  h1_mWW_UL->Write();
  h1_mWW_UL_kinfit->Write();
  h1_mWW_300Mass->Write();
  h1_mWW_medMass->Write();
  h1_mWW_hiMass->Write();
//  h1_mWW_highestMass->Write();
  h1_mWW_WjjMassConstr_300Mass->Write();
  h1_mWW_WjjMassConstr_medMass->Write();
  h1_mWW_WjjMassConstr_hiMass->Write();
  h1_mWW_kinfit_300Mass->Write();
  h1_mWW_kinfit_medMass->Write();
  h1_mWW_kinfit_hiMass->Write();
  h1_mWW_kinfit_cut->Write();
//  h1_mWW_kinfit_highestMass->Write();

  h1_ptWW->Write();
  h1_ptWW_kinfit->Write();
  h1_etaWW->Write();
  h1_etaWW_kinfit->Write();

  h1_deltaR_part1->Write();
  h1_ptJet1->Write();
  h1_eElectronsJet1->Write();
  h1_eMuonsJet1->Write();
  h1_partFlavorJet1->Write();

  h1_deltaR_part2->Write();
  h1_ptJet2->Write();
  h1_eElectronsJet2->Write();
  h1_eMuonsJet2->Write();
  h1_partFlavorJet2->Write();

  h1_deltaRWW->Write();

  h1_mWW_MCassoc->Write();
  h1_mWW_MCassoc_WjjMassConstr->Write();
  h1_mWW_MCassoc_kinfit->Write();
  h1_mWW_MCassoc_kinfit_cands->Write();

  h2_mWjj_vs_mWW->Write();
  h2_mWjj_vs_mWW_kinfit->Write();


  h1_deltaE_ch->Write();
  h1_deltaE_gamma->Write();
  h1_deltaE_nh->Write();
  h1_deltaEta_ch->Write();
  h1_deltaEta_gamma->Write();
  h1_deltaEta_nh->Write();
  h1_deltaPhi_ch->Write();
  h1_deltaPhi_gamma->Write();
  h1_deltaPhi_nh->Write();
  h1_deltaPt_ch->Write();
  h1_deltaPt_gamma->Write();
  h1_deltaPt_nh->Write();

  h1_ptDJet1->Write();
  h1_ptDJet2->Write();

  h1_QGLikelihoodJet1->Write();
  h1_QGLikelihoodJet2->Write();
/*  h1_QGLikelihoodJet1_antiBtag_TChighPur->Write();
  h1_QGLikelihoodJet2_antiBtag_TChighPur->Write();
  h1_QGLikelihoodJet1_antiBtag_TChighEff->Write();
  h1_QGLikelihoodJet2_antiBtag_TChighEff->Write();
  h1_QGLikelihoodJet1_antiBtag_SSVhighPur->Write();
  h1_QGLikelihoodJet2_antiBtag_SSVhighPur->Write();
  h1_QGLikelihoodJet1_antiBtag_SSVhighEff->Write();
  h1_QGLikelihoodJet2_antiBtag_SSVhighEff->Write();
  h1_QGLikelihoodJet1_antiBtag_SSVhighPurhighEff->Write();
  h1_QGLikelihoodJet2_antiBtag_SSVhighPurhighEff->Write();*/
  h1_QGLikelihoodProd->Write();
//  h1_QGLikelihoodProd_antiBtag->Write();
  h1_QGLikelihoodProd_hi->Write();
  h2_QGLikelihoodJet1_vs_Jet2->Write();
  h1_QGLikelihoodRevProd->Write();

  outFile_->mkdir("QGbins");
  outFile_->cd("QGbins");

  for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

    vh1_rmsCandJet1[iPtBin]->Write();
    vh1_ptDJet1[iPtBin]->Write();
    vh1_nChargedJet1[iPtBin]->Write();
    vh1_nNeutralJet1[iPtBin]->Write();
    vh1_QGLikelihoodJet1[iPtBin]->Write();

    vh1_rmsCandJet2[iPtBin]->Write();
    vh1_ptDJet2[iPtBin]->Write();
    vh1_nChargedJet2[iPtBin]->Write();
    vh1_nNeutralJet2[iPtBin]->Write();
    vh1_QGLikelihoodJet2[iPtBin]->Write();

  }


  outFile_->Close();


// btagFile->cd();

// h1_simpleSecondaryVertexHighEffBJetTagJet1->Write();
// h1_simpleSecondaryVertexHighPurBJetTagJet1->Write();
// h1_jetBProbabilityBJetTagJet1->Write();
//  h1_jetProbabilityBJetTagJet1->Write();

//  h1_simpleSecondaryVertexHighEffBJetTagJet2->Write();
//  h1_simpleSecondaryVertexHighPurBJetTagJet2->Write();
//  h1_jetBProbabilityBJetTagJet2->Write();
//  h1_jetProbabilityBJetTagJet2->Write();

//  btagFile->Close();


} // finalize()



void Ntp1Finalizer_HWWlvjj::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;

  if( selectionType_=="presel" ) {

    ptLept1_thresh_ = 10.;
    ptLept2_thresh_ = 10.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 40.;
    mtWll_threshHi_ = 999.;
    mWjj_threshLo_ = 60.;//###
    mWjj_threshHi_ = 100.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mWW_threshLo_ = 0.;
    mWW_threshHi_ = 10000.;

  } else if( selectionType_=="loose" ) {

    ptLept1_thresh_ = 10.;
    ptLept2_thresh_ = 10.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 40.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 60.;
    mtWll_threshHi_ = 100.;
    mWjj_threshLo_ = 60.;//## cambiata
    mWjj_threshHi_ = 100.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mWW_threshLo_ = 0.;
    mWW_threshHi_ = 10000.;
    QGLikelihoodProd_thresh_ = 0.;

  } else if( selectionType=="tight" ) {

    ptLept1_thresh_ = 100.;
    ptLept2_thresh_ = 50.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 100.;
    ptJet2_thresh_ = 50.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 86.;
    mtWll_threshHi_ = 96.;
    mWjj_threshLo_ = 80.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.5;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mWW_threshLo_ = 0.;
    mWW_threshHi_ = 10000.;
 
  } else if( selectionType=="loMass" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 55.;
    ptJet2_thresh_ = 35.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 20.; // Modificato
    mtWll_threshHi_ = 999.;
    mWjj_threshLo_ = 65.;
    mWjj_threshHi_ = 95.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mWW_threshLo_ = 0.;
    mWW_threshHi_ = 999.;

  } else if( selectionType=="opt250LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 70.;
    mtWll_threshHi_ = 110.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.4;
    QGLikelihoodProd_thresh_ = 0.11;
    mWW_threshLo_ = 237.;
    mWW_threshHi_ = 260.;

  } else if( selectionType=="opt300LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 70.;
    mtWll_threshHi_ = 110.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.58;
    QGLikelihoodProd_thresh_ = 0.13;
    mWW_threshLo_ = 280.;
    mWW_threshHi_ = 323.;

  } else if( selectionType=="opt300" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 55.;
    ptJet2_thresh_ = 35.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 40.;//Adattato
    mtWll_threshHi_ = 999.;
    mWjj_threshLo_ = 65.00;
    mWjj_threshHi_ = 95.00;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.7;
    ptWll_thresh_ = 60.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mWW_threshLo_ = 270.;
    mWW_threshHi_ = 330.;

  } else if( selectionType=="opt350LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 70.;
    mtWll_threshHi_ = 110.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.01;
    mWW_threshLo_ = 330.;
    mWW_threshHi_ = 380.;

  } else if( selectionType=="opt400" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 90.;
    ptJet2_thresh_ = 55.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 40.;
    mtWll_threshHi_ = 999.;
    mWjj_threshLo_ = 65.00;
    mWjj_threshHi_ = 95.00;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.2;
    ptWll_thresh_ = 95.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mWW_threshLo_ = 360.;
    mWW_threshHi_ = 440.;

  } else if( selectionType=="opt400noLD" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 90.;
    ptJet2_thresh_ = 55.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 81.;
    mtWll_threshHi_ = 101.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.2;
    ptWll_thresh_ = 95.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.06;
    mWW_threshLo_ = 360.;
    mWW_threshHi_ = 440.;

  } else if( selectionType=="opt400LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 70.;
    mtWll_threshHi_ = 110.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.66;
    QGLikelihoodProd_thresh_ = 0.06;
    mWW_threshLo_ = 390.;
    mWW_threshHi_ = 460.;

  } else if( selectionType=="opt400LD2" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 70.;
    mtWll_threshHi_ = 110.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.03;
    mWW_threshLo_ = 380.;
    mWW_threshHi_ = 470.;

  } else if( selectionType=="opt450LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 70.;
    mtWll_threshHi_ = 110.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.65;
    QGLikelihoodProd_thresh_ = 0.05;
    mWW_threshLo_ = 420.;
    mWW_threshHi_ = 550.;

  } else if( selectionType=="opt500" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 105.;
    ptJet2_thresh_ = 60.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 40.00;//## e' mt;
    mtWll_threshHi_ = 999.00;
    mWjj_threshLo_ = 65.00;
    mWjj_threshHi_ = 95.00;//## cambiata
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.0;
    ptWll_thresh_ = 155.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mWW_threshLo_ = 450.;
    mWW_threshHi_ = 550.;

  } else if( selectionType=="ATLAS" ) {

    ptLept1_thresh_ = 30.;
    ptLept2_thresh_ = 30.;
    etaLept1_thresh_ = 2.4;
    etaLept2_thresh_ = 2.4;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.8;
    etaJet2_thresh_ = 2.8;
    mtWll_threshLo_ = 40.00;//## e' mt;
    mtWll_threshHi_ = 999.00;
    mWjj_threshLo_ = 71.00;
    mWjj_threshHi_ = 91.00;//## cambiata
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 3.0;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mWW_threshLo_ = 350.;
    mWW_threshHi_ = 450.;

  } else if( selectionType=="opt500LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 70.;
    mtWll_threshHi_ = 110.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.75;
    QGLikelihoodProd_thresh_ = 0.15;
    mWW_threshLo_ = 470.;
    mWW_threshHi_ = 99999999.;

  } else if( selectionType=="opt500LD2" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 70.;
    mtWll_threshHi_ = 110.;
    mWjj_threshLo_ = 75.;
    mWjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptWll_thresh_ = 0.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.73;
    QGLikelihoodProd_thresh_ = 0.05;
    mWW_threshLo_ = 465.;
    mWW_threshHi_ = 99999999.;

  } else if( selectionType=="opt600" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 110.;
    ptJet2_thresh_ = 65.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mtWll_threshLo_ = 81.;
    mtWll_threshHi_ = 101.;
    mWjj_threshLo_ = 81.;
    mWjj_threshHi_ = 101.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 0.85;
    ptWll_thresh_ = 175.;
    ptWjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;

  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
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

  TLorentzVector Wll = leptPlus + leptMinus;
  TLorentzVector Wjj = jet1 + jet2;

  TLorentzVector Higgs = Wjj + Wll;

  
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
  TLorentzVector Wll_Hstar(Wll);
  Wll_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector Wjj_Hstar(Wjj);
  Wjj_Hstar.Boost(-Higgs.BoostVector());

  // boosts in Wll CoM frame:
  TLorentzVector lept1_Wllstar(lept1);
  lept1_Wllstar.Boost(-Wll.BoostVector());
  TLorentzVector H_Wllstar(Higgs);
  H_Wllstar.Boost(-Wll.BoostVector());

  // boosts in Wjj CoM frame:
  TLorentzVector jet1_Wjjstar(jet1);
  jet1_Wjjstar.Boost(-Wjj.BoostVector());
  TLorentzVector H_Wjjstar(Higgs);
  H_Wjjstar.Boost(-Wjj.BoostVector());


  returnAngles.helCosThetaStar = Wll_Hstar.CosTheta();


  TVector3 v_pbeamLAB( 0.0, 0.0, 1.0 );

  //cross prod beam x Sll
  TVector3 v_1 = (v_pbeamLAB.Cross(  (Wll_Hstar.Vect()).Unit()) ).Unit();//versor normal to z-z' plane


  //v_2 = cross prod l1 x l2 = versor normal to Wll decay plane
  // careful to the order: L1, the z-axis and W->lv make a right-handed (non-orthogonal) frame (xyz); at the end we want the angle btw x and y
  TVector3 v_2((Wll_Hstar.Vect().Cross(lept1_Hstar.Vect().Unit())).Unit());


  //v_3 = similar to v_2, BUT
  //now, if we want a right-handed set of Unit-vectors, keeping the same direction of the z-axis
  //we must swap the direction of one of the other two vectors of the W bosons. 
  //Keeping the same direction of the z-axis
  //means measuring phiWll and phiWjj w.r.t. to the same v_1 vector (i.e. w.r.t. the same z'-Wll plane)
  TVector3 v_3(((-1.0*Wjj_Hstar.Vect()).Cross(jet1_Hstar.Vect().Unit())).Unit()) ;

  //in other terms: we can define v_3 as above and then do the crss prod with v_1
  //or define v_3 in a way consistent with v_2 and then do the cross product with a newly defined
  //Unit vector v_4 =  (v_pbeamLAB.Cross(  (WjjboostedX->momentum()).Unit()) ).Unit();//versor normal to z-Wjj plane
 
  // helphiWll:
  float phiWll = fabs( acos(v_1.Dot(v_2)) );
  if(v_pbeamLAB.Dot(v_2)>0.0)phiWll=-1.0*phiWll;
  else phiWll=+1.0*phiWll;

  // helphiWjj:
  float phiWjj = fabs( acos(v_1.Dot(v_3)) );
  if(v_pbeamLAB.Dot(v_3)>0.0)phiWjj=+1.0*phiWjj; 
  else phiWjj=-1.0*phiWjj;


  float phi1 = phiWll;


  //phi
  float phi = fabs( acos(v_2.Dot(v_3)) );//two-fold ambiguity when doing the acos + pi ambiguity from sign of v_3 
  if(lept1_Hstar.Vect().Dot(v_3)>0.0)phi= +1.0 * phi;
  else phi= -1.0 * phi;

  returnAngles.helPhi1 = phi1;
  returnAngles.helPhi = phi;


  returnAngles.helCosTheta1 =  (-1.0*(lept1_Wllstar.X()* H_Wllstar.X()+
                                   lept1_Wllstar.Y()* H_Wllstar.Y()+
                                   lept1_Wllstar.Z()* H_Wllstar.Z())/
                                  (lept1_Wllstar.Vect().Mag()* H_Wllstar.Vect().Mag())  );


  returnAngles.helCosTheta2 =  fabs( (jet1_Wjjstar.X()* H_Wjjstar.X()+
                                   jet1_Wjjstar.Y()* H_Wjjstar.Y()+
                                   jet1_Wjjstar.Z()* H_Wjjstar.Z())/
                                  (jet1_Wjjstar.Vect().Mag()* H_Wjjstar.Vect().Mag())  );

  returnAngles.mzz = Higgs.M();


  return returnAngles;

}

// Delta Phi
inline double delta_phi(double phi1, double phi2) {
  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}

// Get Pz Giusto
float getPzRight( TLorentzVector lepton, float pxPFMet, float pyPFMet, TLorentzVector neuMC ) {
  float pn=0., app=0., pznp=0., pznm=0., dr1=0., dr2=0.;
  float a=0., b=0., c=0.;
  TLorentzVector neu1, neu2;

  app = pow(lepton.E(),2)+pow(pxPFMet,2)+pow(pyPFMet,2)-pow(lepton.Px()+pxPFMet,2)-pow(lepton.Py()+pyPFMet,2)-pow(lepton.Pz(),2)-pow(80.399,2);
  a= pow(lepton.E(),2)-pow(lepton.Pz(),2);
  b= lepton.Pz()*app;
  c= ( pow(pxPFMet,2)+pow(pyPFMet,2) )*pow(lepton.E(),2) - pow(app,2)/4.;

  pznp = ( -b + sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  pznm = ( -b - sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  if ( pow(b,2)-4.*a*c < 0. ){ pznp=-b/(2.*a); pznm =-b/(2.*a); }
  
  neu1.SetPxPyPzE( pxPFMet,pyPFMet,pznp,sqrt(pow(pxPFMet,2)+pow(pyPFMet,2)+pow(pznp,2)) );
  neu2.SetPxPyPzE( pxPFMet,pyPFMet,pznm,sqrt(pow(pxPFMet,2)+pow(pyPFMet,2)+pow(pznm,2)) );
  dr1=  neuMC.DeltaR(neu1);
  dr2=  neuMC.DeltaR(neu2); 
  
  if( dr1<=dr2 ){ pn=neu1.Pz(); }
  else{ pn=neu2.Pz();  }  

  return pn;
}

// Get Pz Sbagliato
float getPzWrong( TLorentzVector lepton, float pxPFMet, float pyPFMet, TLorentzVector neuMC ) {
  float pn=0., app=0., pznp=0., pznm=0., dr1=0., dr2=0.;
  float a=0., b=0., c=0.;
  TLorentzVector neu1, neu2;
  
  app = pow(lepton.E(),2)+pow(pxPFMet,2)+pow(pyPFMet,2)-pow(lepton.Px()+pxPFMet,2)-pow(lepton.Py()+pyPFMet,2)-pow(lepton.Pz(),2)-pow(80.399,2);
  a= pow(lepton.E(),2)-pow(lepton.Pz(),2);
  b= lepton.Pz()*app;
  c= ( pow(pxPFMet,2)+pow(pyPFMet,2) )*pow(lepton.E(),2) - pow(app,2)/4.;

  pznp = ( -b + sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  pznm = ( -b - sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  if ( pow(b,2)-4.*a*c < 0. ){ pznp=-b/(2.*a); pznm =-b/(2.*a); }
  
  neu1.SetPxPyPzE( pxPFMet,pyPFMet,pznp,sqrt(pow(pxPFMet,2)+pow(pyPFMet,2)+pow(pznp,2)) );
  neu2.SetPxPyPzE( pxPFMet,pyPFMet,pznm,sqrt(pow(pxPFMet,2)+pow(pyPFMet,2)+pow(pznm,2)) );
  dr1=  neuMC.DeltaR(neu1);
  dr2=  neuMC.DeltaR(neu2); 

  if( dr1>=dr2 ){ pn=neu1.Pz(); }
  else{ pn=neu2.Pz();  }

  return pn;
}

// Get Pz Scegliendo
float getPz( TLorentzVector lepton, float pxPFMet, float pyPFMet,  TLorentzVector jet1, TLorentzVector jet2) {
  float pn=0., app=0., pznp=0., pznm=0.;
  float a=0., b=0., c=0.;
  
  app = pow(lepton.E(),2)+pow(pxPFMet,2)+pow(pyPFMet,2)-pow(lepton.Px()+pxPFMet,2)-pow(lepton.Py()+pyPFMet,2)-pow(lepton.Pz(),2)-pow(80.399,2);
  a= pow(lepton.E(),2)-pow(lepton.Pz(),2);
  b= lepton.Pz()*app;
  c= ( pow(pxPFMet,2)+pow(pyPFMet,2) )*pow(lepton.E(),2) - pow(app,2)/4.;

  pznp = ( -b + sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  pznm = ( -b - sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  if ( pow(b,2)-4.*a*c < 0. ){ pznp=-b/(2.*a); pznm =-b/(2.*a); }
  
  if( fabs(pznp) < fabs(pznm) ){
    pn=pznp; }
  else{ pn=pznm;  }
  
  return pn;
}

float getPzATLAS( TLorentzVector lepton, float pxPFMet, float pyPFMet,  TLorentzVector jet1, TLorentzVector jet2) {
  float pn=0., app=0., pznp=0., pznm=0.;
  float a=0., b=0., c=0.;
  
  app = pow(lepton.E(),2)+pow(pxPFMet,2)+pow(pyPFMet,2)-pow(lepton.Px()+pxPFMet,2)-pow(lepton.Py()+pyPFMet,2)-pow(lepton.Pz(),2)-pow(80.399,2);
  a= pow(lepton.E(),2)-pow(lepton.Pz(),2);
  b= lepton.Pz()*app;
  c= ( pow(pxPFMet,2)+pow(pyPFMet,2) )*pow(lepton.E(),2) - pow(app,2)/4.;

  pznp = ( -b + sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  pznm = ( -b - sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  if ( pow(b,2)-4.*a*c < 0. ){ pznp=-b/(2.*a); pznm =-b/(2.*a); }

  if( fabs(pznp) < fabs(pznm)  ){
    pn=pznp; }
  else{ pn=pznm;  }
  
  return pn;
}

float getPzMH( TLorentzVector lepton, float pxPFMet, float pyPFMet,  TLorentzVector jet1, TLorentzVector jet2) {
  float pn=0., app=0., pznp=0., pznm=0.;
  float a=0., b=0., c=0.;
  TLorentzVector neu1, neu2;
  
  app = pow(lepton.E(),2)+pow(pxPFMet,2)+pow(pyPFMet,2)-pow(lepton.Px()+pxPFMet,2)-pow(lepton.Py()+pyPFMet,2)-pow(lepton.Pz(),2)-pow(80.399,2);
  a= pow(lepton.E(),2)-pow(lepton.Pz(),2);
  b= lepton.Pz()*app;
  c= ( pow(pxPFMet,2)+pow(pyPFMet,2) )*pow(lepton.E(),2) - pow(app,2)/4.;

  pznp = ( -b + sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  pznm = ( -b - sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  if ( pow(b,2)-4.*a*c < 0. ){ pznp=-b/(2.*a); pznm =-b/(2.*a); }
  
  neu1.SetPxPyPzE( pxPFMet,pyPFMet,pznp,sqrt(pow(pxPFMet,2)+pow(pyPFMet,2)+pow(pznp,2)) );
  neu2.SetPxPyPzE( pxPFMet,pyPFMet,pznm,sqrt(pow(pxPFMet,2)+pow(pyPFMet,2)+pow(pznm,2)) );

  if( fabs((neu1+lepton+jet1+jet2).M()-300) < fabs((neu2+lepton+jet1+jet2).M()-300)  ){
    pn=pznp; }
  else{ pn=pznm;  }
  
  return pn;
}
