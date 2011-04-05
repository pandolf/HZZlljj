#include "Ntp1Finalizer_HZZlljj.h"

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
  float simpleSecondaryVertexHighEffBJetTag;
  float simpleSecondaryVertexHighPurBJetTag;
  float jetBProbabilityBJetTag;
  float jetProbabilityBJetTag;

};





HelicityLikelihoodDiscriminant::HelicityAngles computeHelicityAngles(TLorentzVector leptMinus, TLorentzVector leptPlus, TLorentzVector jet1, TLorentzVector jet2 );




int getNJets( int nPairs );

std::vector<TH1D*> getHistoVector(int nPtBins, Double_t *ptBins, std::string histoName, int nBins, float xMin, float xMax );


// constructor:

Ntp1Finalizer_HZZlljj::Ntp1Finalizer_HZZlljj( const std::string& dataset, const std::string& selectionType, const std::string& leptType ) : Ntp1Finalizer( "HZZlljj", dataset, leptType ) {

  if( leptType!="ALL" && leptType!="MU" && leptType!="ELE" ) {
    std::cout << "Lept type '" << leptType << "' currently not supported. Exiting." << std::endl;
    exit(9176);
  }

  leptType_ = leptType;

  setSelectionType(selectionType);

  std::string fullFlags = selectionType_ + "_" + leptType_;
  this->set_flags(fullFlags); //this is for the outfile name

}





void Ntp1Finalizer_HZZlljj::finalize() {

  if( outFile_==0 ) this->createOutputFile();



  TH1F* h1_run = new TH1F("run", "", 15149, 132440, 147589);


  TH1D* h1_ptLept1_presel = new TH1D("ptLept1_presel", "", 25, 20., 300.);
  h1_ptLept1_presel->Sumw2();
  TH1D* h1_ptLept2_presel = new TH1D("ptLept2_presel", "", 25, 20., 150.);
  h1_ptLept2_presel->Sumw2();
  TH1D* h1_etaLept1_presel = new TH1D("etaLept1_presel", "", 25, -2.5, 2.5);
  h1_etaLept1_presel->Sumw2();
  TH1D* h1_etaLept2_presel = new TH1D("etaLept2_presel", "", 25, -2.5, 2.5);
  h1_etaLept2_presel->Sumw2();

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
  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 27, 30., 400.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 27, 30., 150.);
  h1_ptJet2->Sumw2();
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

  TH1D* h1_ptZreso_beforeKin = new TH1D("ptZreso_beforeKin", "", 100, -1., 1.);
  h1_ptZreso_beforeKin->Sumw2();
  TH1D* h1_ptZreso_afterKin = new TH1D("ptZreso_afterKin", "", 100, -1., 1.);
  h1_ptZreso_afterKin->Sumw2();


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

  TH1D* h1_mZjj_all_presel = new TH1D("mZjj_all_presel", "", 20, invMassMin, 400.);
  h1_mZjj_all_presel->Sumw2();

  TH1D* h1_deltaRjj_all_presel = new TH1D("deltaRjj_all_presel", "", 18, 0.5, 5.);
  h1_deltaRjj_all_presel->Sumw2();
  TH1D* h1_deltaRll_presel = new TH1D("deltaRll_presel", "", 20, 0., 5.);
  h1_deltaRll_presel->Sumw2();

  TH1D* h1_mZll = new TH1D("mZll", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZll->Sumw2();
  TH1D* h1_mZll_presel = new TH1D("mZll_presel", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZll_presel->Sumw2();
  TH1D* h1_mZll_presel_0jets = new TH1D("mZll_presel_0jets", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZll_presel_0jets->Sumw2();

  TH1D* h1_mZmumu = new TH1D("mZmumu", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZmumu->Sumw2();
  TH1D* h1_mZmumu_presel = new TH1D("mZmumu_presel", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZmumu_presel->Sumw2();
  TH1D* h1_mZmumu_presel_0jets = new TH1D("mZmumu_presel_0jets", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZmumu_presel_0jets->Sumw2();
  TH1D* h1_mZee = new TH1D("mZee", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZee->Sumw2();
  TH1D* h1_mZee_presel = new TH1D("mZee_presel", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZee_presel->Sumw2();
  TH1D* h1_mZee_presel_0jets = new TH1D("mZee_presel_0jets", "", nBins_invMass, invMassMin_ll, invMassMax);
  h1_mZee_presel_0jets->Sumw2();


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
  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 60, 0., 1.);
  h1_QGLikelihoodJet1->Sumw2();
  TH1D* h1_QGLikelihoodJet1_antiBtag_highEff = new TH1D("QGLikelihoodJet1_antiBtag_highEff", "", 60, 0., 1.);
  h1_QGLikelihoodJet1_antiBtag_highEff->Sumw2();
  TH1D* h1_QGLikelihoodJet1_antiBtag_highPur = new TH1D("QGLikelihoodJet1_antiBtag_highPur", "", 60, 0., 1.);
  h1_QGLikelihoodJet1_antiBtag_highPur->Sumw2();
  TH1D* h1_QGLikelihoodJet1_eta2 = new TH1D("QGLikelihoodJet1_eta2", "", 60, 0., 1.);
  h1_QGLikelihoodJet1_eta2->Sumw2();

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
  TH1D* h1_QGLikelihoodJet2 = new TH1D("QGLikelihoodJet2", "", 60, 0., 1.);
  h1_QGLikelihoodJet2->Sumw2();
  TH1D* h1_QGLikelihoodJet2_antiBtag_highEff = new TH1D("QGLikelihoodJet2_antiBtag_highEff", "", 60, 0., 1.);
  h1_QGLikelihoodJet2_antiBtag_highEff->Sumw2();
  TH1D* h1_QGLikelihoodJet2_antiBtag_highPur = new TH1D("QGLikelihoodJet2_antiBtag_highPur", "", 60, 0., 1.);
  h1_QGLikelihoodJet2_antiBtag_highPur->Sumw2();
  TH1D* h1_QGLikelihoodJet2_eta2 = new TH1D("QGLikelihoodJet2_eta2", "", 60, 0., 1.);
  h1_QGLikelihoodJet2_eta2->Sumw2();

  TH1D* h1_simpleSecondaryVertexHighEffBJetTagJet2 = new TH1D("simpleSecondaryVertexHighEffBJetTagJet2", "", 50, -1.5, 4.);
  h1_simpleSecondaryVertexHighEffBJetTagJet2->Sumw2();
  TH1D* h1_simpleSecondaryVertexHighPurBJetTagJet2 = new TH1D("simpleSecondaryVertexHighPurBJetTagJet2", "", 50, -1.5, 4.);
  h1_simpleSecondaryVertexHighPurBJetTagJet2->Sumw2();
  TH1D* h1_jetBProbabilityBJetTagJet2 = new TH1D("jetBProbabilityBJetTagJet2", "", 50, 0., 8.);
  h1_jetBProbabilityBJetTagJet2->Sumw2();
  TH1D* h1_jetProbabilityBJetTagJet2 = new TH1D("jetProbabilityBJetTagJet2", "", 50, 0., 2.5);
  h1_jetProbabilityBJetTagJet2->Sumw2();

  TH1D* h1_QGLikelihoodProd = new TH1D("QGLikelihoodProd", "", 60, 0., 1.);
  h1_QGLikelihoodProd->Sumw2();
  TH1D* h1_QGLikelihoodProd_hi = new TH1D("QGLikelihoodProd_hi", "", 60, 0., 1.);
  h1_QGLikelihoodProd_hi->Sumw2();
  TH1D* h1_QGLikelihoodRevProd = new TH1D("QGLikelihoodRevProd", "", 60, 0., 1.);
  h1_QGLikelihoodRevProd->Sumw2();
  TH2D* h2_QGLikelihoodJet1_vs_Jet2 = new TH2D("QGLikelihoodJet1_vs_Jet2", "", 60, 0., 1.0001, 60, 0., 1.0001);
  h2_QGLikelihoodJet1_vs_Jet2->Sumw2();


  TH1D* h1_mZjj= new TH1D("mZjj", "", nBins_invMass, 70., 120.);
  h1_mZjj->Sumw2();

  TH1D* h1_ptZll_presel = new TH1D("ptZll_presel", "", 40, 0., 160.);
  h1_ptZll_presel->Sumw2();
  TH1D* h1_ptZjj_all_presel = new TH1D("ptZjj_all_presel", "", 40, 0., 160.);
  h1_ptZjj_all_presel->Sumw2();

  TH1D* h1_ptZll = new TH1D("ptZll", "", 40, 0., 160.);
  h1_ptZll->Sumw2();
  TH1D* h1_ptZjj = new TH1D("ptZjj", "", 40, 0., 160.);
  h1_ptZjj->Sumw2();

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
  TH1D* h1_helicityLD_MW200 = new TH1D("helicityLD_MW200", "", 100, 0., 1.);
  h1_helicityLD_MW200->Sumw2();
  TH1D* h1_helicityLD_MW250 = new TH1D("helicityLD_MW250", "", 100, 0., 1.);
  h1_helicityLD_MW250->Sumw2();
  TH1D* h1_helicityLD_MW300 = new TH1D("helicityLD_MW300", "", 100, 0., 1.);
  h1_helicityLD_MW300->Sumw2();
  TH1D* h1_helicityLD_MW400 = new TH1D("helicityLD_MW400", "", 100, 0., 1.);
  h1_helicityLD_MW400->Sumw2();
  TH1D* h1_helicityLD_MW500 = new TH1D("helicityLD_MW500", "", 100, 0., 1.);
  h1_helicityLD_MW500->Sumw2();
  TH1D* h1_helicityLD_kinfit = new TH1D("helicityLD_kinfit", "", 60, 0., 1.);
  h1_helicityLD_kinfit->Sumw2();

  TH1D* h1_deltaRZZ= new TH1D("deltaRZZ", "", 60, 0., 6.);
  h1_deltaRZZ->Sumw2();

  TH1D* h1_mZZ_merda = new TH1D("mZZ_merda", "", 100, 200., 700.);
  h1_mZZ_merda->Sumw2();
  TH1D* h1_mZZ_UL = new TH1D("mZZ_UL", "", 900, 100., 1000.);
  h1_mZZ_UL->Sumw2();
  TH1D* h1_mZZ_UL_kinfit = new TH1D("mZZ_UL_kinfit", "", 900, 100., 1000.);
  h1_mZZ_UL_kinfit->Sumw2();
  TH1D* h1_mZZ_hiMass= new TH1D("mZZ_hiMass", "", 90, 250., 700.);
  h1_mZZ_hiMass->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass= new TH1D("mZZ_kinfit_hiMass", "", 90, 250., 700.);
  h1_mZZ_kinfit_hiMass->Sumw2();
//TH1D* h1_mZZ_highestMass= new TH1D("mZZ_highestMass", "", 70, 350., 700.);
//h1_mZZ_highestMass->Sumw2();
//TH1D* h1_mZZ_kinfit_highestMass= new TH1D("mZZ_kinfit_highestMass", "", 70, 350., 700.);
//h1_mZZ_kinfit_highestMass->Sumw2();
  TH1D* h1_mZZ_ZjjMassConstr_hiMass  = new TH1D("mZZ_ZjjMassConstr_hiMass", "", 70, 250., 600.);
  h1_mZZ_ZjjMassConstr_hiMass->Sumw2();
  TH1D* h1_mZZ_300Mass= new TH1D("mZZ_300Mass", "", 50, 250., 500.);
  h1_mZZ_300Mass->Sumw2();
  TH1D* h1_mZZ_ZjjMassConstr_300Mass  = new TH1D("mZZ_ZjjMassConstr_300Mass", "", 50, 250., 500.);
  h1_mZZ_ZjjMassConstr_300Mass->Sumw2();
  TH1D* h1_mZZ_kinfit_300Mass= new TH1D("mZZ_kinfit_300Mass", "", 50, 250., 500.);
  h1_mZZ_kinfit_300Mass->Sumw2();
  TH1D* h1_mZZ_medMass= new TH1D("mZZ_medMass", "", 70, 150., 350.);
  h1_mZZ_medMass->Sumw2();
  TH1D* h1_mZZ_ZjjMassConstr_medMass  = new TH1D("mZZ_ZjjMassConstr_medMass", "", 70, 150., 350.);
  h1_mZZ_ZjjMassConstr_medMass->Sumw2();
  TH1D* h1_mZZ_kinfit_medMass= new TH1D("mZZ_kinfit_medMass", "", 70, 150., 350.);
  h1_mZZ_kinfit_medMass->Sumw2();

  TH1D* h1_ptZZ  = new TH1D("ptZZ", "", 100, 0., 300.);
  h1_ptZZ->Sumw2();
  TH1D* h1_ptZZ_kinfit  = new TH1D("ptZZ_kinfit", "", 100, 0., 300.);
  h1_ptZZ_kinfit->Sumw2();
  TH1D* h1_etaZZ  = new TH1D("etaZZ", "", 100, -5.5, 5.5);
  h1_etaZZ->Sumw2();
  TH1D* h1_etaZZ_kinfit  = new TH1D("etaZZ_kinfit", "", 100, -5.5, 5.5);
  h1_etaZZ_kinfit->Sumw2();
  

  TH1D* h1_mZZ_MCassoc  = new TH1D("mZZ_MCassoc", "", 100, 200., 600.);
  h1_mZZ_MCassoc->Sumw2();
  TH1D* h1_mZZ_MCassoc_ZjjMassConstr  = new TH1D("mZZ_MCassoc_ZjjMassConstr", "", 100, 200., 600.);
  h1_mZZ_MCassoc_ZjjMassConstr->Sumw2();
  TH1D* h1_mZZ_MCassoc_kinfit  = new TH1D("mZZ_MCassoc_kinfit", "", 100, 200., 600.);
  h1_mZZ_MCassoc_kinfit->Sumw2();
  TH1D* h1_mZZ_MCassoc_kinfit_cands  = new TH1D("mZZ_MCassoc_kinfit_cands", "", 100, 200., 600.);
  h1_mZZ_MCassoc_kinfit_cands->Sumw2();

  TH1D* h1_partFlavor_tight = new TH1D("partFlavor_tight", "", 31, -8.5, 22.5);

  TH2D* h2_mZjj_vs_mZZ = new TH2D("mZjj_vs_mZZ", "", 100, 200., 600., 100, 70., 120.);
  h2_mZjj_vs_mZZ->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_kinfit = new TH2D("mZjj_vs_mZZ_kinfit", "", 100, 200., 600., 100, 70., 120.);
  h2_mZjj_vs_mZZ_kinfit->Sumw2();


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


  float nEvents_pre=0.;
  float nEvents_pre_leptPt=0.;
  float nEvents_pre_leptPt_leptMass =0.;
  float nEvents_pre_leptPt_leptMass_jetPt=0.;
  float nEvents_pre_leptPt_leptMass_jetPt_jetMass=0.;
  float nEvents_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj=0.;

  float nEventsPassed_fb_kinfit=0.;
  float nEventsPassed_fb_nokinfit=0.;
  int nEventsPassed_kinfit=0;
  int nEventsPassed_nokinfit=0;


  int nEntries = tree_->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_3_8_7/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root", nPtBins);


ofstream ofs("run_event.txt");


  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

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



    TLorentzVector lept1, lept2;
    lept1.SetPtEtaPhiE( ptLept1, etaLept1, phiLept1, eLept1 );
    lept2.SetPtEtaPhiE( ptLept2, etaLept2, phiLept2, eLept2 );

    TLorentzVector diLepton = lept1+lept2;



    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_selected;

    if( nPairs>0 ) {

      h1_mZll_presel->Fill( diLepton.M(), eventWeight );
      if( leptType==0 )
        h1_mZmumu_presel->Fill( diLepton.M(), eventWeight );
      else
        h1_mZee_presel->Fill( diLepton.M(), eventWeight );
      h1_ptZll_presel->Fill( diLepton.Pt(), eventWeight );

      h1_deltaRll_presel->Fill(lept1.DeltaR(lept2), eventWeight );

      h1_ptLept1_presel->Fill( ptLept1, eventWeight );
      h1_ptLept2_presel->Fill( ptLept2, eventWeight );

      h1_etaLept1_presel->Fill( etaLept1, eventWeight );
      h1_etaLept2_presel->Fill( etaLept2, eventWeight );

      int nJets = getNJets(nPairs);
      h1_nJets_presel->Fill( nJets , eventWeight );
      h1_nPairs_presel->Fill( nPairs , eventWeight );

    } else {

      h1_mZll_presel_0jets->Fill( diLepton.M(), eventWeight );
      if( leptType==0 )
        h1_mZmumu_presel_0jets->Fill( diLepton.M(), eventWeight );
      else
        h1_mZee_presel_0jets->Fill( diLepton.M(), eventWeight );

    }



    float cached_jetpt = 0.;

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

      jet2.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet2[iJetPair];
      jet2.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet2[iJetPair];
      jet2.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet2[iJetPair];
      jet2.jetProbabilityBJetTag               = jetProbabilityBJetTagJet2[iJetPair];

      jet2.ptGen =   ptJet2Gen[iJetPair];
      jet2.etaGen = etaJet2Gen[iJetPair];
      jet2.phiGen = phiJet2Gen[iJetPair];
      jet2.eGen =     eJet2Gen[iJetPair];

      TLorentzVector diJet = jet1 + jet2;

    //if( jet1.Pt()>ptJet1_thresh_ && jet2.Pt()>ptJet2_thresh_ && fabs(jet1.Eta())<etaJet1_thresh_ && fabs(jet2.Eta())<etaJet1_thresh_ 
    // && jet1.DeltaR(jet2) < deltaRjj_thresh_ && diJet.M() > mZjj_threshLo_ && diJet.M() < mZjj_threshHi_ && diJet.Pt() > ptZjj_thresh_ )
    //  jetPairs_selected.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );

      if( jet1.Pt()>ptJet1_thresh_ && jet2.Pt()>ptJet2_thresh_ && fabs(jet1.Eta())<etaJet1_thresh_ && fabs(jet2.Eta())<etaJet1_thresh_ 
       && jet1.DeltaR(jet2) < deltaRjj_thresh_ && diJet.M() > mZjj_threshLo_ && diJet.M() < mZjj_threshHi_ && diJet.Pt() > ptZjj_thresh_ )
//       && jet1.electronEnergyFraction < 0.5 && jet1.muonEnergyFraction < 0.4 && jet2.electronEnergyFraction < 0.5 && jet2.muonEnergyFraction < 0.4 )
        jetPairs_selected.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );


      h1_mZjj_all_presel->Fill( diJet.M(), eventWeight );
      h1_ptZjj_all_presel->Fill( diJet.Pt(), eventWeight );
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

    }


    TLorentzVector jet1_presel, jet2_presel;

    if( jetPairs_selected.size()>0. ) {

      // now look for best Z mass jet pair:
      float Zmass = 91.1876;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_selected.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_selected[iPair].first + jetPairs_selected[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs



      AnalysisJet jet1 = jetPairs_selected[bestPair].first;
      AnalysisJet jet2 = jetPairs_selected[bestPair].second;

      TLorentzVector bestZDiJet = jet1 + jet2;
  

      if( lept1.Pt() > ptLept1_thresh_ && lept2.Pt() > ptLept2_thresh_ && fabs(lept1.Eta()) < etaLept1_thresh_ && fabs(lept2.Eta()) < etaLept2_thresh_
       && diLepton.M() > mZll_threshLo_ && diLepton.M() < mZll_threshHi_ && lept1.DeltaR(lept2) < deltaRll_thresh_ && diLepton.Pt() > ptZll_thresh_ ) {

        // event has passed kinematic selection


        // ------------------------
        //   KINEMATIC FIT : BEGIN
        // ------------------------

        DiJetKinFitter* fitter_jets = new DiJetKinFitter( "fitter_jets", "fitter_jets", Zmass );
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


        TLorentzVector Zjj_kinfit = jet1_kinfit + jet2_kinfit;

        TLorentzVector matchedZ;
        float bestDeltaRZ=999.;
        for( unsigned iPart=0; iPart<nPart; ++iPart ) {
          if( pdgIdPart[iPart]!=23 ) continue;
          TLorentzVector thisZ;
          thisZ.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
          if( Zjj_kinfit.DeltaR(thisZ) < bestDeltaRZ ) {
            bestDeltaRZ = Zjj_kinfit.DeltaR(thisZ);
            matchedZ = thisZ;
          }
        }

        float ptZreso_before = (isMC) ? (bestZDiJet.Pt()-matchedZ.Pt())/matchedZ.Pt() : 0.;
        float ptZreso_after  = (isMC) ? (Zjj_kinfit.Pt()-matchedZ.Pt())/matchedZ.Pt() : 0.;
        h1_ptZreso_beforeKin->Fill( ptZreso_before, eventWeight);
        h1_ptZreso_afterKin->Fill( ptZreso_after, eventWeight);

        TLorentzVector Zjj_constr;
        Zjj_constr.SetXYZM( bestZDiJet.Px(), bestZDiJet.Py(), bestZDiJet.Pz(), Zmass);

        TLorentzVector ZZ = bestZDiJet + diLepton;
        TLorentzVector ZZ_constr = diLepton + Zjj_constr;
        TLorentzVector ZZ_kinfit = diLepton + Zjj_kinfit;

        h1_kinfit_chiSquare->Fill( fitter_jets->getS()/fitter_jets->getNDF(), eventWeight ); 
        h1_kinfit_chiSquareProb->Fill( TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF()), eventWeight ); 
if( TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF())<0.2 ) h1_mZZ_merda->Fill( ZZ_kinfit.M(), eventWeight );


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
    //  h1_mZZ_MCassoc->Fill( ZZ.M(), eventWeight );
    //  h1_mZZ_MCassoc_ZjjMassConstr->Fill( ZZ_constr.M(), eventWeight );
    //  h1_mZZ_MCassoc_kinfit->Fill( ZZ_kinfit.M(), eventWeight );
    //}

//    // and now full kinematic fit with PFCands:
//    TLorentzVector ZZ_kinfit_cands;
//    TRegexp cands_tstr("CANDS");
//    TString dataset_tstr(dataset_);
//if( dataset_tstr.Contains(cands_tstr) ) {

//    std::vector<TFitParticleEtEtaPhi*> fitCands;
//    std::vector<int> fitCandTypes;
//    TFitConstraintM *mCons_cands = new TFitConstraintM( "ZMassConstraint_cands", "ZMass-Constraint", 0, 0 , 91.19);


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

//    // recreate Zqq vector:
//    TLorentzVector Zqq_kinfit_cands;
//    for( unsigned iCand=0; iCand<fitCands.size(); ++iCand )
//      Zqq_kinfit_cands += *(fitCands[iCand]->getCurr4Vec());


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

//    ZZ_kinfit_cands = Zll + Zqq_kinfit_cands;

//    if( bothMatched )
//      h1_mZZ_MCassoc_kinfit_cands->Fill( ZZ_kinfit_cands.M(), eventWeight );
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
      if( chargeLept1<0 ) hangles = computeHelicityAngles(lept1, lept2, jet1, jet2);
      else                hangles = computeHelicityAngles(lept2, lept1, jet1, jet2);

      HelicityLikelihoodDiscriminant::HelicityAngles hangles_kinfit;
      if( chargeLept1<0 ) hangles_kinfit = computeHelicityAngles(lept1, lept2, jet1_kinfit, jet2_kinfit);
      else                hangles_kinfit = computeHelicityAngles(lept2, lept1, jet1_kinfit, jet2_kinfit);


      HelicityLikelihoodDiscriminant *LD = new HelicityLikelihoodDiscriminant();


      LD->setMeasurables(hangles);
      double sProb=LD->getSignalProbability();
      double bProb=LD->getBkgdProbability();
      double helicityLD=sProb/(sProb+bProb);


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

      bool btag_highPur = ( jet1.simpleSecondaryVertexHighPurBJetTag>2. || jet2.simpleSecondaryVertexHighPurBJetTag>2. );
      bool btag_highEff = ( jet1.simpleSecondaryVertexHighEffBJetTag>2. || jet2.simpleSecondaryVertexHighEffBJetTag>2. );

      h1_simpleSecondaryVertexHighEffBJetTagJet1->Fill(jet1.simpleSecondaryVertexHighEffBJetTag, eventWeight);
      h1_simpleSecondaryVertexHighPurBJetTagJet1->Fill(jet1.simpleSecondaryVertexHighPurBJetTag, eventWeight);
      h1_jetBProbabilityBJetTagJet1->Fill(jet1.jetBProbabilityBJetTag, eventWeight);
      h1_jetProbabilityBJetTagJet1->Fill(jet1.jetProbabilityBJetTag, eventWeight);

      h1_simpleSecondaryVertexHighEffBJetTagJet2->Fill(jet2.simpleSecondaryVertexHighEffBJetTag, eventWeight);
      h1_simpleSecondaryVertexHighPurBJetTagJet2->Fill(jet2.simpleSecondaryVertexHighPurBJetTag, eventWeight);
      h1_jetBProbabilityBJetTagJet2->Fill(jet2.jetBProbabilityBJetTag, eventWeight);
      h1_jetProbabilityBJetTagJet2->Fill(jet2.jetProbabilityBJetTag, eventWeight);

      vh1_rmsCandJet1[jet1PtBin]->Fill( jet1.rmsCand, eventWeight );
      vh1_ptDJet1[jet1PtBin]->Fill( jet1.ptD, eventWeight );
      vh1_nChargedJet1[jet1PtBin]->Fill( jet1.nCharged, eventWeight );
      vh1_nNeutralJet1[jet1PtBin]->Fill( jet1.nNeutral, eventWeight );
      float QGLikelihoodJet1 = qglikeli->computeQGLikelihood( jet1.Pt(), jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );

      vh1_QGLikelihoodJet1[jet1PtBin]->Fill( QGLikelihoodJet1, eventWeight );
      h1_QGLikelihoodJet1->Fill( QGLikelihoodJet1, eventWeight );
      if( !btag_highEff ) h1_QGLikelihoodJet1_antiBtag_highEff->Fill( QGLikelihoodJet1, eventWeight );
      if( !btag_highPur ) h1_QGLikelihoodJet1_antiBtag_highPur->Fill( QGLikelihoodJet1, eventWeight );
      if( fabs(jet1.Eta())<2. ) h1_QGLikelihoodJet1_eta2->Fill(QGLikelihoodJet1, eventWeight);
      
      vh1_rmsCandJet2[jet2PtBin]->Fill( jet2.rmsCand, eventWeight );
      vh1_ptDJet2[jet2PtBin]->Fill( jet2.ptD, eventWeight );
      vh1_nChargedJet2[jet2PtBin]->Fill( jet2.nCharged, eventWeight );
      vh1_nNeutralJet2[jet2PtBin]->Fill( jet2.nNeutral, eventWeight );
      float QGLikelihoodJet2 = qglikeli->computeQGLikelihood( jet2.Pt(), jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );

      vh1_QGLikelihoodJet2[jet2PtBin]->Fill( QGLikelihoodJet2, eventWeight );
      h1_QGLikelihoodJet2->Fill( QGLikelihoodJet2, eventWeight );
      if( !btag_highEff ) h1_QGLikelihoodJet2_antiBtag_highEff->Fill( QGLikelihoodJet2, eventWeight );
      if( !btag_highPur ) h1_QGLikelihoodJet2_antiBtag_highPur->Fill( QGLikelihoodJet2, eventWeight );
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

      if( ZZ_kinfit.M() > mZZ_threshLo_ && ZZ_kinfit.M() < mZZ_threshHi_ ) {
        nEventsPassed_fb_kinfit += eventWeight;
        nEventsPassed_kinfit++;
      }
      if( ZZ.M() > mZZ_threshLo_ && ZZ.M() < mZZ_threshHi_ ) {
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
        h1_deltaRjj->Fill( jet1.DeltaR(jet2), eventWeight);
        h1_ptZll->Fill( diLepton.Pt(), eventWeight);
        h1_ptZjj->Fill( bestZDiJet.Pt(), eventWeight);
        if( leptType==0 )
          h1_mZmumu->Fill( diLepton.M(), eventWeight );
        else
          h1_mZee->Fill( diLepton.M(), eventWeight );
        h1_mZll->Fill( diLepton.M(), eventWeight);
        h1_mZjj->Fill( bestZDiJet.M(), eventWeight);
        h1_mZZ_UL->Fill(ZZ.M(), eventWeight);
        h1_mZZ_hiMass->Fill(ZZ.M(), eventWeight);
      //  h1_mZZ_highestMass->Fill(ZZ.M(), eventWeight);
        h1_mZZ_medMass->Fill(ZZ.M(), eventWeight);
        h1_mZZ_300Mass->Fill(ZZ.M(), eventWeight);
        h1_mZZ_ZjjMassConstr_hiMass->Fill(ZZ_constr.M(), eventWeight);
        h1_mZZ_ZjjMassConstr_medMass->Fill(ZZ_constr.M(), eventWeight);
        h1_mZZ_ZjjMassConstr_300Mass->Fill(ZZ_constr.M(), eventWeight);
        h1_mZZ_UL_kinfit->Fill(ZZ_kinfit.M(), eventWeight);
        h1_mZZ_kinfit_hiMass->Fill(ZZ_kinfit.M(), eventWeight);
      //  h1_mZZ_kinfit_highestMass->Fill(ZZ_kinfit.M(), eventWeight);
        h1_mZZ_kinfit_medMass->Fill(ZZ_kinfit.M(), eventWeight);
        h1_mZZ_kinfit_300Mass->Fill(ZZ_kinfit.M(), eventWeight);
        h2_mZjj_vs_mZZ->Fill( ZZ.M(), bestZDiJet.M() );
        h2_mZjj_vs_mZZ_kinfit->Fill( ZZ_kinfit.M(), bestZDiJet.M() );


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





  } //for entries


  std::cout << std::endl << std::endl;
  std::cout << "----> PASSED SELECTION: " << 1000.*nEventsPassed_fb_kinfit << " ev/fb-1  (" << nEventsPassed_kinfit << " events)" << std::endl;
  std::cout << "----> PASSED SELECTION (no kinfit): " << 1000.*nEventsPassed_fb_nokinfit << " ev/fb-1 (" << nEventsPassed_nokinfit << " events)" << std::endl;
  std::cout << std::endl;





  outFile_->cd();

  h1_run->Write();

  h1_ptJet_all_presel->Write();
  h1_etaJet_all_presel->Write();
  h1_nJets_presel->Write();
  h1_nPairs_presel->Write();

  h1_ptResoJet1_beforeKin->Write();
  h1_ptResoJet2_beforeKin->Write();
  h1_ptResoJet1_afterKin->Write();
  h1_ptResoJet2_afterKin->Write();

  h1_ptZreso_beforeKin->Write();
  h1_ptZreso_afterKin->Write();

  h1_deltaRll_presel->Write();
  h1_deltaRjj_all_presel->Write();

  h1_ptLept1->Write();
  h1_ptLept2->Write();

  h1_ptLept1_presel->Write();
  h1_ptLept2_presel->Write();
  h1_etaLept1_presel->Write();
  h1_etaLept2_presel->Write();

  h1_mZll->Write();
  h1_mZll_presel->Write();
  h1_mZmumu->Write();
  h1_mZmumu_presel->Write();
  h1_mZmumu_presel_0jets->Write();
  h1_mZee->Write();
  h1_mZee_presel->Write();
  h1_mZee_presel_0jets->Write();

  h1_mZjj->Write();
  h1_mZjj_all_presel->Write();

  h1_ptZll_presel->Write();
  h1_ptZjj_all_presel->Write();

  h1_deltaRjj->Write();

  h1_ptZjj->Write();
  h1_ptZll->Write();

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
  h1_helicityLD_MW200->Write();
  h1_helicityLD_MW250->Write();
  h1_helicityLD_MW300->Write();
  h1_helicityLD_MW400->Write();
  h1_helicityLD_MW500->Write();
  h1_helicityLD_kinfit->Write();

  h1_mZZ_merda->Write();
  h1_mZZ_UL->Write();
  h1_mZZ_UL_kinfit->Write();
  h1_mZZ_300Mass->Write();
  h1_mZZ_medMass->Write();
  h1_mZZ_hiMass->Write();
//  h1_mZZ_highestMass->Write();
  h1_mZZ_ZjjMassConstr_300Mass->Write();
  h1_mZZ_ZjjMassConstr_medMass->Write();
  h1_mZZ_ZjjMassConstr_hiMass->Write();
  h1_mZZ_kinfit_300Mass->Write();
  h1_mZZ_kinfit_medMass->Write();
  h1_mZZ_kinfit_hiMass->Write();
//  h1_mZZ_kinfit_highestMass->Write();

  h1_ptZZ->Write();
  h1_ptZZ_kinfit->Write();
  h1_etaZZ->Write();
  h1_etaZZ_kinfit->Write();

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

  h1_deltaRZZ->Write();

  h1_mZZ_MCassoc->Write();
  h1_mZZ_MCassoc_ZjjMassConstr->Write();
  h1_mZZ_MCassoc_kinfit->Write();
  h1_mZZ_MCassoc_kinfit_cands->Write();

  h2_mZjj_vs_mZZ->Write();
  h2_mZjj_vs_mZZ_kinfit->Write();


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

  h1_QGLikelihoodJet1->Write();
  h1_QGLikelihoodJet2->Write();
  h1_QGLikelihoodJet1_antiBtag_highEff->Write();
  h1_QGLikelihoodJet2_antiBtag_highEff->Write();
  h1_QGLikelihoodJet1_antiBtag_highPur->Write();
  h1_QGLikelihoodJet2_antiBtag_highPur->Write();
  h1_QGLikelihoodProd->Write();
  h1_QGLikelihoodProd_hi->Write();
  h2_QGLikelihoodJet1_vs_Jet2->Write();
  h1_QGLikelihoodRevProd->Write();

  outFile_->mkdir("QGbins");
  outFile_->cd("QGbins");

  for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

    vh1_rmsCandJet2[iPtBin]->Write();
    vh1_ptDJet2[iPtBin]->Write();
    vh1_nChargedJet2[iPtBin]->Write();
    vh1_nNeutralJet2[iPtBin]->Write();
    vh1_QGLikelihoodJet2[iPtBin]->Write();

  }


  outFile_->Close();


  std::string btagFileName = "btagfile_" + selectionType_ + ".root";
  TFile* btagFile = TFile::Open(btagFileName.c_str(), "RECREATE");
  btagFile->cd();

  h1_simpleSecondaryVertexHighEffBJetTagJet1->Write();
  h1_simpleSecondaryVertexHighPurBJetTagJet1->Write();
  h1_jetBProbabilityBJetTagJet1->Write();
  h1_jetProbabilityBJetTagJet1->Write();

  h1_simpleSecondaryVertexHighEffBJetTagJet2->Write();
  h1_simpleSecondaryVertexHighPurBJetTagJet2->Write();
  h1_jetBProbabilityBJetTagJet2->Write();
  h1_jetProbabilityBJetTagJet2->Write();

  btagFile->Close();


} // finalize()



void Ntp1Finalizer_HZZlljj::setSelectionType( const std::string& selectionType ) {

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
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 70.;
    mZjj_threshHi_ = 110.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;

  } else if( selectionType_=="loose" ) {

    ptLept1_thresh_ = 10.;
    ptLept2_thresh_ = 10.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 40.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 70.;
    mZjj_threshHi_ = 110.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 0.;
    mZZ_threshHi_ = 10000.;

  } else if( selectionType=="tight" ) {

    ptLept1_thresh_ = 100.;
    ptLept2_thresh_ = 50.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 100.;
    ptJet2_thresh_ = 50.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 86.;
    mZll_threshHi_ = 96.;
    mZjj_threshLo_ = 80.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.5;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;

  } else if( selectionType=="loMass" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 55.;
    ptJet2_thresh_ = 35.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 40.;
    mZll_threshHi_ = 80.;
    mZjj_threshLo_ = 81.;
    mZjj_threshHi_ = 101.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;

  } else if( selectionType=="opt250LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.4;
    QGLikelihoodProd_thresh_ = 0.11;
    mZZ_threshLo_ = 237.;
    mZZ_threshHi_ = 260.;

  } else if( selectionType=="opt300LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.58;
    QGLikelihoodProd_thresh_ = 0.13;
    mZZ_threshLo_ = 280.;
    mZZ_threshHi_ = 323.;

  } else if( selectionType=="opt300" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 55.;
    ptJet2_thresh_ = 35.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 81.;
    mZjj_threshHi_ = 101.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.7;
    ptZll_thresh_ = 60.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 270.;
    mZZ_threshHi_ = 330.;

  } else if( selectionType=="opt350LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.01;
    mZZ_threshLo_ = 330.;
    mZZ_threshHi_ = 380.;

  } else if( selectionType=="opt400" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 90.;
    ptJet2_thresh_ = 55.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
  //mZll_threshLo_ = 81.;
  //mZll_threshHi_ = 101.;
    mZll_threshLo_ = 50.;
    mZll_threshHi_ = 120.;
    mZjj_threshLo_ = 81.;
    mZjj_threshHi_ = 101.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.2;
    ptZll_thresh_ = 95.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 360.;
    mZZ_threshHi_ = 440.;

  } else if( selectionType=="opt400noLD" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 90.;
    ptJet2_thresh_ = 55.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.2;
    ptZll_thresh_ = 95.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.06;
    mZZ_threshLo_ = 360.;
    mZZ_threshHi_ = 440.;

  } else if( selectionType=="opt400LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.66;
    QGLikelihoodProd_thresh_ = 0.06;
    mZZ_threshLo_ = 390.;
    mZZ_threshHi_ = 460.;

  } else if( selectionType=="opt400LD2" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.03;
    mZZ_threshLo_ = 380.;
    mZZ_threshHi_ = 470.;

  } else if( selectionType=="opt450LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.65;
    QGLikelihoodProd_thresh_ = 0.05;
    mZZ_threshLo_ = 420.;
    mZZ_threshHi_ = 550.;

  } else if( selectionType=="opt500" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 105.;
    ptJet2_thresh_ = 60.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 81.;
    mZjj_threshHi_ = 101.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.0;
    ptZll_thresh_ = 155.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 450.;
    mZZ_threshHi_ = 550.;

  } else if( selectionType=="opt500LD" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.75;
    QGLikelihoodProd_thresh_ = 0.15;
    mZZ_threshLo_ = 470.;
    mZZ_threshHi_ = 99999999.;

  } else if( selectionType=="opt500LD2" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 9999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.73;
    QGLikelihoodProd_thresh_ = 0.05;
    mZZ_threshLo_ = 465.;
    mZZ_threshHi_ = 99999999.;

  } else if( selectionType=="opt600" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 110.;
    ptJet2_thresh_ = 65.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 81.;
    mZjj_threshHi_ = 101.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 0.85;
    ptZll_thresh_ = 175.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;

  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType


// error functions for jets:

/*
// calojet resolution:
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
*/


/*

// pfjet resolutions. taken from AN-2010-371
Double_t ErrEt(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, N, S, C, m;
  if(fabs(Eta) < 0.5 ) {
    N = 3.96859;
    S = 0.18348;
    C = 0.;
    m = 0.62627;
  } else if( fabs(Eta) < 1. ) {
    N = 3.55226;
    S = 0.24026;
    C = 0.;
    m = 0.52571;
  } else if( fabs(Eta) < 1.5 ) {
    N = 4.54826;
    S = 0.22652;
    C = 0.;
    m = 0.58963;
  } else if( fabs(Eta) < 2. ) {
    N = 4.62622;
    S = 0.23664;
    C = 0.;
    m = 0.48738;
  } else if( fabs(Eta) < 2.5 ) {
    N = 2.53324;
    S = 0.34306;
    C = 0.;
    m = 0.28662;
  } else if( fabs(Eta) < 3. ) {
    N = -3.33814;
    S = 0.73360;
    C = 0.;
    m = 0.08264;
  } else if( fabs(Eta) < 5. ) {
    N = 2.95397;
    S = 0.11619;
    C = 0.;
    m = 0.96086;
  }


  // this is the absolute resolution (squared), not sigma(pt)/pt
  // so have to multiply by pt^2, thats why m+1 instead of m-1
  InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et ;
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

*/

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
  TLorentzVector Zjj = jet1 + jet2;

  TLorentzVector Higgs = Zjj + Zll;

  
  // define lept1 as negatively charged lepton:
  TLorentzVector lept1 = leptMinus;
  TLorentzVector lept2 = leptPlus;

  // no charge for jets (well, not really)
  // so choose jet with positive scalar product with Zjj 
  // in its restframe:
  TLorentzVector jet1_Zjjstar_tmp(jet1);
  jet1_Zjjstar_tmp.Boost(-Zjj.BoostVector());
  if( jet1_Zjjstar_tmp.Phi()<0. ) { //swap them
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
  TLorentzVector Zjj_Hstar(Zjj);
  Zjj_Hstar.Boost(-Higgs.BoostVector());

  // boosts in Zll CoM frame:
  TLorentzVector lept1_Zllstar(lept1);
  lept1_Zllstar.Boost(-Zll.BoostVector());
  TLorentzVector H_Zllstar(Higgs);
  H_Zllstar.Boost(-Zll.BoostVector());

  // boosts in Zjj CoM frame:
  TLorentzVector jet1_Zjjstar(jet1);
  jet1_Zjjstar.Boost(-Zjj.BoostVector());
  TLorentzVector H_Zjjstar(Higgs);
  H_Zjjstar.Boost(-Zjj.BoostVector());


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
  //means measuring phiZll and phiZjj w.r.t. to the same v_1 vector (i.e. w.r.t. the same z'-Zll plane)
  TVector3 v_3(((-1.0*Zjj_Hstar.Vect()).Cross(jet1_Hstar.Vect().Unit())).Unit()) ;

  //in other terms: we can define v_3 as above and then do the crss prod with v_1
  //or define v_3 in a way consistent with v_2 and then do the cross product with a newly defined
  //Unit vector v_4 =  (v_pbeamLAB.Cross(  (ZjjboostedX->momentum()).Unit()) ).Unit();//versor normal to z-Zjj plane
 
  // helphiZll:
  float phiZll = fabs( acos(v_1.Dot(v_2)) );
  if(v_pbeamLAB.Dot(v_2)>0.0)phiZll=-1.0*phiZll;
  else phiZll=+1.0*phiZll;

  // helphiZjj:
  float phiZjj = fabs( acos(v_1.Dot(v_3)) );
  if(v_pbeamLAB.Dot(v_3)>0.0)phiZjj=+1.0*phiZjj; 
  else phiZjj=-1.0*phiZjj;


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


  returnAngles.helCosTheta2 =  fabs( (jet1_Zjjstar.X()* H_Zjjstar.X()+
                                   jet1_Zjjstar.Y()* H_Zjjstar.Y()+
                                   jet1_Zjjstar.Z()* H_Zjjstar.Z())/
                                  (jet1_Zjjstar.Vect().Mag()* H_Zjjstar.Vect().Mag())  );

  returnAngles.mzz = Higgs.M();


  return returnAngles;

}
