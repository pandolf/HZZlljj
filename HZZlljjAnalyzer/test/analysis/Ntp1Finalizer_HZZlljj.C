#include "Ntp1Finalizer_HZZlljj.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"

#include "AnalysisJet.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "CommonTools/fitTools.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"
#include "KinematicFit/DiJetKinFitter.h"




/*

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
  float QGLikelihood;

  float muonEnergyFraction;
  float electronEnergyFraction;

  float pt_preKinFit;
  float eta_preKinFit;
  float phi_preKinFit;
  float e_preKinFit;

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
*/




//HelicityLikelihoodDiscriminant::HelicityAngles computeHelicityAngles(TLorentzVector leptMinus, TLorentzVector leptPlus, TLorentzVector jet1, TLorentzVector jet2 );




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


void Ntp1Finalizer_HZZlljj::set_useLooseBTags( bool useLooseBTags ) {
  
  useLooseBTags_ = useLooseBTags;

  std::string flags = flags_;

  flags += "_looseBTag";

  this->set_flags(flags);

}


void Ntp1Finalizer_HZZlljj::finalize() {

  if( outFile_==0 ) this->createOutputFile();


  std::string btagFileName = "btagfile_" + selectionType_ + ".root";
  TFile* btagFile = TFile::Open(btagFileName.c_str(), "RECREATE");

  TH1F* h1_run = new TH1F("run", "", 15149, 132440, 147589);


  TH1D* h1_rhoPF_presel = new TH1D("rhoPF_presel", "", 50, 0., 20.);
  h1_rhoPF_presel->Sumw2();
  TH1D* h1_rhoPF = new TH1D("rhoPF", "", 50, 0., 20.);
  h1_rhoPF->Sumw2();

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

//std::vector<TH1D*> vh1_ptDJet_all_presel = getHistoVector(nPtBins, ptBins, "ptDJet_all_presel", 60, 0., 1.);
//std::vector<TH1D*> vh1_rmsCandJet_all_presel = getHistoVector(nPtBins, ptBins, "rmsCandJet_all_presel", 60, 0., 0.07);
//std::vector<TH1D*> vh1_nChargedJet_all_presel = getHistoVector(nPtBins, ptBins, "nChargedJet_all_presel", 41, -0.5, 40.5);
//std::vector<TH1D*> vh1_nNeutralJet_all_presel = getHistoVector(nPtBins, ptBins, "nNeutralJet_all_presel", 41, -0.5, 40.5);

  TH1D* h1_nJets_presel = new TH1D("nJets_presel", "", 7, 1.5, 8.5);
  h1_nJets_presel->Sumw2();
  TH1D* h1_nPairs_presel = new TH1D("nPairs_presel", "", 21, 0.5, 21.5);
  h1_nPairs_presel->Sumw2();

  TH1D* h1_ptLept1= new TH1D("ptLept1", "", 56, 20., 300.);
  h1_ptLept1->Sumw2();
  TH1D* h1_ptLept2= new TH1D("ptLept2", "", 52, 20., 150.);
  h1_ptLept2->Sumw2();
  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 80, 0., 400.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 60, 0., 150.);
  h1_ptJet2->Sumw2();
  TH1D* h1_ptJet1_MW = new TH1D("ptJet1_MW", "", 80, 0., 400.);
  h1_ptJet1_MW->Sumw2();
  TH1D* h1_ptJet2_MW = new TH1D("ptJet2_MW", "", 60, 0., 150.);
  h1_ptJet2_MW->Sumw2();
  TH1D* h1_ptJet1_prekin = new TH1D("ptJet1_prekin", "", 54, 30., 400.);
  h1_ptJet1_prekin->Sumw2();
  TH1D* h1_ptJet2_prekin = new TH1D("ptJet2_prekin", "", 48, 30., 150.);
  h1_ptJet2_prekin->Sumw2();
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

  TH1D* h1_ptDJet1 = new TH1D("ptDJet1", "", 60, 0., 1.0001);
  h1_ptDJet1->Sumw2();
  TH1D* h1_ptDJet2 = new TH1D("ptDJet2", "", 60, 0., 1.0001);
  h1_ptDJet2->Sumw2();

  TH1D* h1_deltaR_part1 = new TH1D("deltaR_part1", "", 50, 0., 0.8);
  h1_deltaR_part1->Sumw2();
  //TH1D* h1_partFlavorJet1= new TH1D("partFlavorJet1", "", 30, -7.5, 22.5);
  TH1D* h1_partFlavorJet1= new TH1D("partFlavorJet1", "", 38, -15.5, 22.5);
  h1_partFlavorJet1->Sumw2();
//std::vector<TH1D*> vh1_ptDJet1 = getHistoVector(nPtBins, ptBins, "ptDJet1", 60, 0., 1.);
//std::vector<TH1D*> vh1_rmsCandJet1 = getHistoVector(nPtBins, ptBins, "rmsCandJet1", 60, 0., 0.1);
//std::vector<TH1D*> vh1_nChargedJet1 = getHistoVector(nPtBins, ptBins, "nChargedJet1", 51, -0.5, 50.5);
//std::vector<TH1D*> vh1_nNeutralJet1 = getHistoVector(nPtBins, ptBins, "nNeutralJet1", 51, -0.5, 50.5);
//std::vector<TH1D*> vh1_QGLikelihoodJet1 = getHistoVector(nPtBins, ptBins, "QGLikelihoodJet1", 60, 0., 1.);
  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet1->Sumw2();

  TH1D* h1_QGLikelihood_100_123 = new TH1D("QGLikelihood_100_123", "", 50, 0., 1.0001);
  h1_QGLikelihood_100_123->Sumw2();
  TH1D* h1_QGLikelihood_66_81 = new TH1D("QGLikelihood_66_81", "", 50, 0., 1.0001);
  h1_QGLikelihood_66_81->Sumw2();

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
  TH1D* h1_QGLikelihoodProd_MW = new TH1D("QGLikelihoodProd_MW", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_MW->Sumw2();


  TH1D* h1_mZjj= new TH1D("mZjj", "", nBins_invMass, 70., 120.);
  h1_mZjj->Sumw2();
  TH1D* h1_mZjj_MW = new TH1D("mZjj_MW", "", nBins_invMass, 70., 120.);
  h1_mZjj_MW->Sumw2();
  TH1D* h1_mZjj_loChiSquareProb= new TH1D("mZjj_loChiSquareProb", "", 100, 30., 200.);
  h1_mZjj_loChiSquareProb->Sumw2();
  TH1D* h1_mZjj_hiChiSquareProb= new TH1D("mZjj_hiChiSquareProb", "", 100, 30., 200.);
  h1_mZjj_hiChiSquareProb->Sumw2();

  TH1D* h1_ptZll_presel = new TH1D("ptZll_presel", "", 40, 0., 160.);
  h1_ptZll_presel->Sumw2();
  TH1D* h1_ptZjj_all_presel = new TH1D("ptZjj_all_presel", "", 40, 0., 160.);
  h1_ptZjj_all_presel->Sumw2();

  TH1D* h1_ptZll = new TH1D("ptZll", "", 50, 0., 300.);
  h1_ptZll->Sumw2();
  TH1D* h1_ptZll_MW = new TH1D("ptZll_MW", "", 50, 0., 300.);
  h1_ptZll_MW->Sumw2();
  TH1D* h1_ptZjj = new TH1D("ptZjj", "", 50, 0., 300.);
  h1_ptZjj->Sumw2();

  TH1D* h1_deltaRjj= new TH1D("deltaRjj", "", 50, 0.5, 5.);
  h1_deltaRjj->Sumw2();
  TH1D* h1_deltaRjj_prekin= new TH1D("deltaRjj_prekin", "", 50, 0.5, 5.);
  h1_deltaRjj_prekin->Sumw2();
  TH1D* h1_deltaRjj_MW= new TH1D("deltaRjj_MW", "", 50, 0.5, 5.);
  h1_deltaRjj_MW->Sumw2();

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

  TH1D* h1_mZZ_hiChiSquareProb = new TH1D("mZZ_hiChiSquareProb", "", 200, 100., 700.);
  h1_mZZ_hiChiSquareProb->Sumw2();
  TH1D* h1_mZZ_loChiSquareProb = new TH1D("mZZ_loChiSquareProb", "", 200, 100., 700.);
  h1_mZZ_loChiSquareProb->Sumw2();
  TH1D* h1_mZZ_mZjj_cut = new TH1D("mZZ_mZjj_cut", "", 200, 100., 700.);
  h1_mZZ_mZjj_cut->Sumw2();
  TH1D* h1_mZZ_mZjj_notcut = new TH1D("mZZ_mZjj_notcut", "", 200, 100., 700.);
  h1_mZZ_mZjj_notcut->Sumw2();
  TH1D* h1_mZZ_UL = new TH1D("mZZ_UL", "", 900, 100., 1000.);
  h1_mZZ_UL->Sumw2();
  TH1D* h1_mZZ_UL_kinfit = new TH1D("mZZ_UL_kinfit", "", 900, 100., 1000.);
  h1_mZZ_UL_kinfit->Sumw2();
  TH1D* h1_mZZ_hiMass= new TH1D("mZZ_hiMass", "", 90, 250., 700.);
  h1_mZZ_hiMass->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass= new TH1D("mZZ_kinfit_hiMass", "", 90, 250., 700.);
  h1_mZZ_kinfit_hiMass->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_hiQG= new TH1D("mZZ_kinfit_hiMass_hiQG", "", 90, 250., 700.);
  h1_mZZ_kinfit_hiMass_hiQG->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_loQG= new TH1D("mZZ_kinfit_hiMass_loQG", "", 90, 250., 700.);
  h1_mZZ_kinfit_hiMass_loQG->Sumw2();
//TH1D* h1_mZZ_highestMass= new TH1D("mZZ_highestMass", "", 70, 350., 700.);
//h1_mZZ_highestMass->Sumw2();
//TH1D* h1_mZZ_kinfit_highestMass= new TH1D("mZZ_kinfit_highestMass", "", 70, 350., 700.);
//h1_mZZ_kinfit_highestMass->Sumw2();
  TH1D* h1_mZZ_ZjjMassConstr_hiMass  = new TH1D("mZZ_ZjjMassConstr_hiMass", "", 70, 250., 600.);
  h1_mZZ_ZjjMassConstr_hiMass->Sumw2();
  TH1D* h1_mZZ_300Mass= new TH1D("mZZ_300Mass", "", 50, 250., 450.);
  h1_mZZ_300Mass->Sumw2();
  TH1D* h1_mZZ_ZjjMassConstr_300Mass  = new TH1D("mZZ_ZjjMassConstr_300Mass", "", 50, 250., 450.);
  h1_mZZ_ZjjMassConstr_300Mass->Sumw2();
  TH1D* h1_mZZ_kinfit_300Mass= new TH1D("mZZ_kinfit_300Mass", "", 50, 250., 450.);
  h1_mZZ_kinfit_300Mass->Sumw2();
  TH1D* h1_mZZ_kinfit_300Mass_loQG= new TH1D("mZZ_kinfit_300Mass_loQG", "", 50, 250., 450.);
  h1_mZZ_kinfit_300Mass_loQG->Sumw2();
  TH1D* h1_mZZ_kinfit_300Mass_hiQG= new TH1D("mZZ_kinfit_300Mass_hiQG", "", 50, 250., 450.);
  h1_mZZ_kinfit_300Mass_hiQG->Sumw2();
  TH1D* h1_mZZ_medMass= new TH1D("mZZ_medMass", "", 70, 150., 350.);
  h1_mZZ_medMass->Sumw2();
  TH1D* h1_mZZ_ZjjMassConstr_medMass  = new TH1D("mZZ_ZjjMassConstr_medMass", "", 70, 150., 350.);
  h1_mZZ_ZjjMassConstr_medMass->Sumw2();
  TH1D* h1_mZZ_kinfit_medMass= new TH1D("mZZ_kinfit_medMass", "", 70, 150., 350.);
  h1_mZZ_kinfit_medMass->Sumw2();
  TH1D* h1_mZZ_kinfit_medMass_hiQG= new TH1D("mZZ_kinfit_medMass_hiQG", "", 70, 150., 350.);
  h1_mZZ_kinfit_medMass_hiQG->Sumw2();
  TH1D* h1_mZZ_kinfit_medMass_loQG= new TH1D("mZZ_kinfit_medMass_loQG", "", 70, 150., 350.);
  h1_mZZ_kinfit_medMass_loQG->Sumw2();

  TH1D* h1_ptZZ  = new TH1D("ptZZ", "", 100, 0., 300.);
  h1_ptZZ->Sumw2();
  TH1D* h1_ptZZ_kinfit  = new TH1D("ptZZ_kinfit", "", 100, 0., 300.);
  h1_ptZZ_kinfit->Sumw2();
  TH1D* h1_etaZZ  = new TH1D("etaZZ", "", 100, -5.5, 5.5);
  h1_etaZZ->Sumw2();
  TH1D* h1_etaZZ_kinfit  = new TH1D("etaZZ_kinfit", "", 100, -5.5, 5.5);
  h1_etaZZ_kinfit->Sumw2();
  

//TH1D* h1_mZZ_MCassoc  = new TH1D("mZZ_MCassoc", "", 100, 200., 600.);
//h1_mZZ_MCassoc->Sumw2();
//TH1D* h1_mZZ_MCassoc_ZjjMassConstr  = new TH1D("mZZ_MCassoc_ZjjMassConstr", "", 100, 200., 600.);
//h1_mZZ_MCassoc_ZjjMassConstr->Sumw2();
//TH1D* h1_mZZ_MCassoc_kinfit  = new TH1D("mZZ_MCassoc_kinfit", "", 100, 200., 600.);
//h1_mZZ_MCassoc_kinfit->Sumw2();
//TH1D* h1_mZZ_MCassoc_kinfit_cands  = new TH1D("mZZ_MCassoc_kinfit_cands", "", 100, 200., 600.);
//h1_mZZ_MCassoc_kinfit_cands->Sumw2();

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





  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);
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


  


  float nEventsPassed_fb_kinfit=0.;
  float nEventsPassed_fb_nokinfit=0.;
  int nEventsPassed_kinfit=0;
  int nEventsPassed_kinfit_antiBtag=0;
  int nEventsPassed_nokinfit=0;


  int nEntries = tree_->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_3_8_7/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root", nPtBins);
  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_1_3/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1.root");
  float Zmass = 91.1876;
  DiJetKinFitter* fitter_jets = new DiJetKinFitter( "fitter_jets", "fitter_jets", Zmass );
  HelicityLikelihoodDiscriminant *LD = new HelicityLikelihoodDiscriminant();
  float helicityLD_selected = -1.;
  float helicityLD_kinfit_selected = -1.;


float nEventsTot = 0.;
float nEvents_hiChiSquareProb = 0.;
float nEvents_mZjj_cut = 0.;

ofstream ofs("run_event.txt");




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



    h1_rhoPF_presel->Fill( rhoPF, eventWeight);


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




    // --------------------
    // FULL EVENT VARIABLES
    // --------------------
 
    if( pfMetThresh_ > 0. && pfMet > pfMetThresh_ ) continue;

    


    // ----------------------------
    // KINEMATIC SELECTION: LEPTONS
    // ----------------------------

    if( lept1.Pt() < ptLept1_thresh_ ) continue;
    if( lept2.Pt() < ptLept2_thresh_ ) continue;
    if( fabs(lept1.Eta()) > etaLept1_thresh_ ) continue;
    if( fabs(lept2.Eta()) > etaLept2_thresh_ ) continue;
    if( diLepton.M() < mZll_threshLo_ || diLepton.M() > mZll_threshHi_ ) continue;
    if( lept1.DeltaR(lept2) > deltaRll_thresh_ ) continue;
    if( diLepton.Pt() < ptZll_thresh_ ) continue;




    float cached_jetpt = 0.;
    AnalysisJet jet1_selected, jet2_selected;
    float bestMass = 0.;
    bool foundJets = false;
    bool foundJets_ZZmass = false;

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


      // fill histos before selection
      h1_mZjj_all_presel->Fill( diJet.M(), eventWeight );
      h1_ptZjj_all_presel->Fill( diJet.Pt(), eventWeight );
      h1_deltaRjj_all_presel->Fill(jet1.DeltaR(jet2), eventWeight );

      if( jet1.Pt()!=cached_jetpt ) {
        h1_ptJet_all_presel->Fill( jet1.Pt(), eventWeight );
        h1_etaJet_all_presel->Fill( jet1.Eta(), eventWeight );
        cached_jetpt = jet1.Pt();
      }



      // -------------------------
      // KINEMATIC SELECTION: JETS
      // -------------------------

      if( jet1.Pt() < ptJet1_thresh_ ) continue;
      if( jet2.Pt() < ptJet2_thresh_ ) continue;
      if( fabs(jet1.Eta()) > etaJet1_thresh_ ) continue;
      if( fabs(jet2.Eta()) > etaJet2_thresh_ ) continue;
      if( jet1.DeltaR(jet2) > deltaRjj_thresh_ ) continue;
      if( diJet.M() < mZjj_threshLo_ || diJet.M() > mZjj_threshHi_ ) continue;
      if( diJet.Pt() < ptZjj_thresh_ ) continue;



      // -------------------
      // Q-G DISCRIMINATION:
      // -------------------

//    float QGLikelihoodJet1 = qglikeli->computeQGLikelihoodPU( jet1.Pt(), rhoPF, jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );
//    float QGLikelihoodJet2 = qglikeli->computeQGLikelihoodPU( jet2.Pt(), rhoPF, jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );
      float QGLikelihoodJet1 = qglikeli->computeQGLikelihood( jet1.Pt(), jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );
      float QGLikelihoodJet2 = qglikeli->computeQGLikelihood( jet2.Pt(), jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );
      float QGLikelihoodProd = QGLikelihoodJet1*QGLikelihoodJet2;
      if( QGLikelihoodProd < QGLikelihoodProd_thresh_ ) continue;

      jet1.QGLikelihood = QGLikelihoodJet1;
      jet2.QGLikelihood = QGLikelihoodJet2;


      // ----------
      // B-TAGGING:
      // ----------

//    bool btagJet1_loose =  ( jet1.trackCountingHighEffBJetTag>1.83 ); 
//    bool btagJet2_loose =  ( jet2.trackCountingHighEffBJetTag>1.83 ); 

//    bool btagJet1_medium =  ( jet1.trackCountingHighEffBJetTag>4. ); 
//    bool btagJet2_medium =  ( jet2.trackCountingHighEffBJetTag>4. ); 

      int nBTags;

      if( useLooseBTags_ ) {

        bool twoBTags  = ( jet1.btag_medium() && jet2.btag_loose()  )
                      || ( jet1.btag_loose()  && jet2.btag_medium() );
        bool oneBTag   = (!twoBTags) && ( jet1.btag_loose() || jet2.btag_loose() );
        //bool zeroBTags = (!twoBTags) && (!oneBTag);

        if( twoBTags ) nBTags=2;
        else if( oneBTag ) nBTags=1;
        else nBTags=0;

      } else {

        nBTags = jet1.btag_medium() + jet2.btag_medium();

      }

      if( requiredBTags_ != -1 ) if( nBTags != requiredBTags_ ) continue;



      // --------------
      // KINEMATIC FIT:
      // --------------

      std::pair<TLorentzVector,TLorentzVector> jets_kinfit = fitter_jets->fit(jet1, jet2);
      TLorentzVector jet1_kinfit(jets_kinfit.first);
      TLorentzVector jet2_kinfit(jets_kinfit.second);

      jet1.pt_preKinFit  = jet1.Pt();
      jet1.eta_preKinFit = jet1.Eta();
      jet1.phi_preKinFit = jet1.Phi();
      jet1.e_preKinFit   = jet1.Energy();

      jet2.pt_preKinFit  = jet2.Pt();
      jet2.eta_preKinFit = jet2.Eta();
      jet2.phi_preKinFit = jet2.Phi();
      jet2.e_preKinFit   = jet2.Energy();

      // update 4-vector to kinfit results:
      jet1.SetPtEtaPhiE( jet1_kinfit.Pt(), jet1_kinfit.Eta(), jet1_kinfit.Phi(), jet1_kinfit.Energy() );
      jet2.SetPtEtaPhiE( jet2_kinfit.Pt(), jet2_kinfit.Eta(), jet2_kinfit.Phi(), jet2_kinfit.Energy() );

      TLorentzVector diJet_kinfit = jet1_kinfit + jet2_kinfit;
      TLorentzVector ZZ_kinfit_tmp = diJet_kinfit + diLepton;



      // ------------
      // HELICITY LD:
      // ------------

      HelicityLikelihoodDiscriminant::HelicityAngles hangles;
      if( chargeLept1<0 ) hangles = LD->computeHelicityAngles(lept1, lept2, jet1, jet2);
      else                hangles = LD->computeHelicityAngles(lept2, lept1, jet1, jet2);
    
      LD->setMeasurables(hangles);
      double sProb=LD->getSignalProbability();
      double bProb=LD->getBkgdProbability();
      double helicityLD=sProb/(sProb+bProb);
    
      HelicityLikelihoodDiscriminant::HelicityAngles hangles_kinfit;
      if( chargeLept1<0 ) hangles_kinfit = LD->computeHelicityAngles(lept1, lept2, jet1_kinfit, jet2_kinfit);
      else                hangles_kinfit = LD->computeHelicityAngles(lept2, lept1, jet1_kinfit, jet2_kinfit);
    
      LD->setMeasurables(hangles_kinfit);
      double sProb_kinfit=LD->getSignalProbability();
      double bProb_kinfit=LD->getBkgdProbability();
      double helicityLD_kinfit=sProb_kinfit/(sProb_kinfit+bProb_kinfit);
     
      // still cutting on helicity LD pre-kinfit!! 
      if( helicityLD < helicityLD_thresh_ ) continue;



      // choose jet pair in ZZ mass window. if none found, take best Z mjj pair:

      if( ZZ_kinfit_tmp.M() > mZZ_threshLo_ && ZZ_kinfit_tmp.M() < mZZ_threshHi_ && !foundJets_ZZmass ) {
        nEventsPassed_fb_kinfit += eventWeight;
        nEventsPassed_kinfit++;
        jet1_selected = jet1;
        jet2_selected = jet2;
        helicityLD_selected = helicityLD;
        helicityLD_kinfit_selected = helicityLD_kinfit;
        foundJets = true;
        foundJets_ZZmass = true;
        // fill MW (mass window) histograms:
        h1_ptZll_MW->Fill( diLepton.Pt(), eventWeight );
        h1_ptJet1_MW->Fill( jet1_selected.Pt(), eventWeight );
        h1_ptJet2_MW->Fill( jet2_selected.Pt(), eventWeight );
        h1_deltaRjj_MW->Fill( jet1_selected.DeltaR(jet2_selected), eventWeight );
        h1_QGLikelihoodProd_MW->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
        h1_mZjj_MW->Fill( diJet.M(), eventWeight );
      }

      float invMass = diJet.M();
      if( !foundJets_ZZmass ) {
        if( !foundJets || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          jet1_selected = jet1;
          jet2_selected = jet2;
          helicityLD_selected = helicityLD;
          foundJets = true;
        }
      }

    } //for on jet pairs



    if( !foundJets ) continue;

    if( helicityLD_selected < 0. ) 
      std::cout << "helicityLD_selected is less than 0!!! THIS IS NOT POSSIBLE!!" << std::endl;

  //if( helicityLD_kinfit_selected < 0. ) 
  //  std::cout << "helicityLD_kinfit_selected is less than 0!!! THIS IS NOT POSSIBLE!!" << std::endl;



    TLorentzVector Zjj_kinfit = jet1_selected + jet2_selected;
  
    TLorentzVector jet1_nokinfit, jet2_nokinfit;
    jet1_nokinfit.SetPtEtaPhiE( jet1_selected.pt_preKinFit, jet1_selected.eta_preKinFit, jet1_selected.phi_preKinFit, jet1_selected.e_preKinFit );
    jet2_nokinfit.SetPtEtaPhiE( jet2_selected.pt_preKinFit, jet2_selected.eta_preKinFit, jet2_selected.phi_preKinFit, jet2_selected.e_preKinFit );

    TLorentzVector Zjj_nokinfit = jet1_nokinfit + jet2_nokinfit;

    TLorentzVector ZZ_nokinfit = Zjj_nokinfit + diLepton;
    TLorentzVector ZZ_kinfit = diLepton + Zjj_kinfit;

    if( ZZ_nokinfit.M() > mZZ_threshLo_ && ZZ_nokinfit.M() < mZZ_threshHi_ ) {
      nEventsPassed_fb_nokinfit += eventWeight;
      nEventsPassed_nokinfit++;
    }

    bool btag_TChighPur =  ( jet1_selected.trackCountingHighPurBJetTag>5.         || jet2_selected.trackCountingHighPurBJetTag>5. );
    bool btag_TChighEff =  ( jet1_selected.trackCountingHighEffBJetTag>5.         || jet2_selected.trackCountingHighEffBJetTag>5. );
    bool btag_SSVhighPur = ( jet1_selected.simpleSecondaryVertexHighPurBJetTag>2. || jet2_selected.simpleSecondaryVertexHighPurBJetTag>2. );
    bool btag_SSVhighEff = ( jet1_selected.simpleSecondaryVertexHighEffBJetTag>2. || jet2_selected.simpleSecondaryVertexHighEffBJetTag>2. );
    bool btag_SSVhighPurhighEff = ( ( jet1_selected.simpleSecondaryVertexHighPurBJetTag>2. && jet2_selected.simpleSecondaryVertexHighEffBJetTag>2. ) ||
                                    ( jet1_selected.simpleSecondaryVertexHighEffBJetTag>2. && jet2_selected.simpleSecondaryVertexHighPurBJetTag>2. ) );

    h1_simpleSecondaryVertexHighEffBJetTagJet1->Fill(jet1_selected.simpleSecondaryVertexHighEffBJetTag, eventWeight);
    h1_simpleSecondaryVertexHighPurBJetTagJet1->Fill(jet1_selected.simpleSecondaryVertexHighPurBJetTag, eventWeight);
    h1_jetBProbabilityBJetTagJet1->Fill(jet1_selected.jetBProbabilityBJetTag, eventWeight);
    h1_jetProbabilityBJetTagJet1->Fill(jet1_selected.jetProbabilityBJetTag, eventWeight);

    h1_simpleSecondaryVertexHighEffBJetTagJet2->Fill(jet2_selected.simpleSecondaryVertexHighEffBJetTag, eventWeight);
    h1_simpleSecondaryVertexHighPurBJetTagJet2->Fill(jet2_selected.simpleSecondaryVertexHighPurBJetTag, eventWeight);
    h1_jetBProbabilityBJetTagJet2->Fill(jet2_selected.jetBProbabilityBJetTag, eventWeight);
    h1_jetProbabilityBJetTagJet2->Fill(jet2_selected.jetProbabilityBJetTag, eventWeight);

    bool twoTags = (jet1_selected.trackCountingHighEffBJetTag>=4. && jet2_selected.trackCountingHighEffBJetTag>=4.);
    bool oneTag = (jet1_selected.trackCountingHighEffBJetTag>=4. && jet2_selected.trackCountingHighEffBJetTag<4.)
                || (jet1_selected.trackCountingHighEffBJetTag<4. && jet2_selected.trackCountingHighEffBJetTag>=4.);
    bool zeroTags = (jet1_selected.trackCountingHighEffBJetTag<4. && jet2_selected.trackCountingHighEffBJetTag<4.);


    // match to partons:
    TLorentzVector matchedPart1, matchedPart2;
    float bestDeltaRPart1=999.;
    float bestDeltaRPart2=999.;
    for( unsigned iPart=0; iPart<nPart; ++iPart ) {
      if( abs(pdgIdPart[iPart])>6 ) continue;
      TLorentzVector thisPart;
      thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
      if( jet1_selected.DeltaR(thisPart) < bestDeltaRPart1 ) {
        bestDeltaRPart1 = jet1_selected.DeltaR(thisPart);
        matchedPart1 = thisPart;
      }
      if( jet2_selected.DeltaR(thisPart) < bestDeltaRPart2 ) {
        bestDeltaRPart2 = jet2_selected.DeltaR(thisPart);
        matchedPart2 = thisPart;
      }
    }
    float ptReso1_before = (isMC) ? ( jet1_nokinfit.Pt()-matchedPart1.Pt() )/matchedPart1.Pt() : 0.;
    float ptReso2_before = (isMC) ? ( jet2_nokinfit.Pt()-matchedPart2.Pt() )/matchedPart2.Pt() : 0.;
    h1_ptResoJet1_beforeKin->Fill( ptReso1_before, eventWeight );
    h1_ptResoJet2_beforeKin->Fill( ptReso2_before, eventWeight );

    float ptReso1_after = (isMC) ? ( jet1_selected.Pt()-matchedPart1.Pt() )/matchedPart1.Pt() : 0.;
    float ptReso2_after = (isMC) ? ( jet2_selected.Pt()-matchedPart2.Pt() )/matchedPart2.Pt() : 0.;
    h1_ptResoJet1_afterKin->Fill( ptReso1_after, eventWeight );
    h1_ptResoJet2_afterKin->Fill( ptReso2_after, eventWeight );



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

    float ptZreso_before = (isMC) ? (Zjj_nokinfit.Pt()-matchedZ.Pt())/matchedZ.Pt() : 0.;
    float ptZreso_after  = (isMC) ? (Zjj_kinfit.Pt()-matchedZ.Pt())/matchedZ.Pt() : 0.;
    h1_ptZreso_beforeKin->Fill( ptZreso_before, eventWeight);
    h1_ptZreso_afterKin->Fill( ptZreso_after, eventWeight);


    //compare kinfit to Zmass constraint
    TLorentzVector Zjj_constr;
    Zjj_constr.SetXYZM( Zjj_nokinfit.Px(), Zjj_nokinfit.Py(), Zjj_nokinfit.Pz(), Zmass);

    TLorentzVector ZZ_constr = diLepton + Zjj_constr;


    float chiSquareProb = TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF());
    h1_kinfit_chiSquare->Fill( fitter_jets->getS()/fitter_jets->getNDF(), eventWeight ); 
    h1_kinfit_chiSquareProb->Fill( chiSquareProb, eventWeight ); 


 // nEventsTot += eventWeight;
 // if( chiSquareProb>0.01 ) {
 //   nEvents_hiChiSquareProb += eventWeight;
 //   h1_mZjj_hiChiSquareProb->Fill( bestZDiJet.M(), eventWeight );
 //   h1_mZZ_hiChiSquareProb->Fill( ZZ_kinfit.M(), eventWeight );
 // } else {
 //   h1_mZjj_loChiSquareProb->Fill( bestZDiJet.M(), eventWeight );
 //   h1_mZZ_loChiSquareProb->Fill( ZZ_kinfit.M(), eventWeight );
 // }
 // if( bestZDiJet.M()>75. && bestZDiJet.M()<105. ) {
 //   nEvents_mZjj_cut += eventWeight;
 //   h1_mZZ_mZjj_cut->Fill( ZZ_kinfit.M(), eventWeight );
 // } else {
 //   h1_mZZ_mZjj_notcut->Fill( ZZ_kinfit.M(), eventWeight );
 // }







 // h1_ptDJet1->Fill( jet1_selected.ptD, eventWeight );
 // h1_ptDJet2->Fill( jet2_selected.ptD, eventWeight );

 // vh1_rmsCandJet1[jet1PtBin]->Fill( jet1_selected.rmsCand, eventWeight );
 // vh1_ptDJet1[jet1PtBin]->Fill( jet1_selected.ptD, eventWeight );
 // vh1_nChargedJet1[jet1PtBin]->Fill( jet1_selected.nCharged, eventWeight );
 // vh1_nNeutralJet1[jet1PtBin]->Fill( jet1_selected.nNeutral, eventWeight );

 // //vh1_QGLikelihoodJet1[jet1PtBin]->Fill( QGLikelihoodJet1, eventWeight );
 // //h1_QGLikelihoodJet1->Fill( QGLikelihoodJet1, eventWeight );
 // //h1_QGLikelihoodJet1_4var->Fill( QGLikelihoodJet1_4var, eventWeight );
 // 
 // vh1_rmsCandJet2[jet2PtBin]->Fill( jet2_selected.rmsCand, eventWeight );
 // vh1_ptDJet2[jet2PtBin]->Fill( jet2_selected.ptD, eventWeight );
 // vh1_nChargedJet2[jet2PtBin]->Fill( jet2_selected.nCharged, eventWeight );
 // vh1_nNeutralJet2[jet2PtBin]->Fill( jet2_selected.nNeutral, eventWeight );

    //vh1_QGLikelihoodJet2[jet2PtBin]->Fill( QGLikelihoodJet2, eventWeight );
    //h1_QGLikelihoodJet2->Fill( QGLikelihoodJet2, eventWeight );
  //h1_QGLikelihoodJet2_4var->Fill( QGLikelihoodJet2_4var, eventWeight );



  //h1_QGLikelihoodProd->Fill( QGLikelihoodProd, eventWeight );
  //h1_QGLikelihoodProd_4var->Fill( QGLikelihoodProd_4var, eventWeight );
  //if( !btag_SSVhighEff ) h1_QGLikelihoodProd_antiBtag->Fill( QGLikelihoodProd, eventWeight );
  //
  //h1_QGLikelihoodRevProd->Fill( QGLikelihoodRevProd, eventWeight );
  //  
  //h2_QGLikelihoodJet1_vs_Jet2->Fill( QGLikelihoodJet2, QGLikelihoodJet1, eventWeight );

  //if( QGLikelihoodJet1>0.8 || QGLikelihoodJet2>0.8 ) h1_QGLikelihoodProd_hi->Fill(QGLikelihoodProd, eventWeight); 



    // fill histograms:

    h1_rhoPF->Fill( rhoPF, eventWeight );


    if( jet1_selected.Pt()>jet2_selected.Pt() ) {
      h1_ptJet1->Fill( jet1_selected.Pt(), eventWeight );
      h1_ptJet2->Fill( jet2_selected.Pt(), eventWeight );
      h1_ptJet1_prekin->Fill( jet1_nokinfit.Pt(), eventWeight );
      h1_ptJet2_prekin->Fill( jet2_nokinfit.Pt(), eventWeight );
    } else {
      h1_ptJet1->Fill( jet2_selected.Pt(), eventWeight );
      h1_ptJet2->Fill( jet1_selected.Pt(), eventWeight );
      h1_ptJet1_prekin->Fill( jet2_nokinfit.Pt(), eventWeight );
      h1_ptJet2_prekin->Fill( jet1_nokinfit.Pt(), eventWeight );
    }
    h1_eMuonsJet1->Fill( jet1_selected.muonEnergyFraction, eventWeight );
    h1_eMuonsJet2->Fill( jet2_selected.muonEnergyFraction, eventWeight );
    h1_eElectronsJet1->Fill( jet1_selected.electronEnergyFraction, eventWeight );
    h1_eElectronsJet2->Fill( jet2_selected.electronEnergyFraction, eventWeight );

    h1_ptLept1->Fill( lept1.Pt(), eventWeight );
    h1_ptLept2->Fill( lept2.Pt(), eventWeight );
    h1_deltaRjj->Fill( jet1_selected.DeltaR(jet2_selected), eventWeight);
    h1_deltaRjj_prekin->Fill( jet1_nokinfit.DeltaR(jet2_nokinfit), eventWeight);
    h1_ptZll->Fill( diLepton.Pt(), eventWeight);
    h1_ptZjj->Fill( Zjj_kinfit.Pt(), eventWeight);
    if( leptType==0 )
      h1_mZmumu->Fill( diLepton.M(), eventWeight );
    else
      h1_mZee->Fill( diLepton.M(), eventWeight );
    h1_mZll->Fill( diLepton.M(), eventWeight);
    h1_mZjj->Fill( Zjj_nokinfit.M(), eventWeight);

    

    h1_QGLikelihoodJet1->Fill( jet1_selected.QGLikelihood, eventWeight );
    h1_QGLikelihoodJet2->Fill( jet2_selected.QGLikelihood, eventWeight );
    h1_QGLikelihoodProd->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );

    if( jet1_selected.Pt()>100. && jet1_selected.Pt()<123. ) h1_QGLikelihood_100_123->Fill( jet1_selected.QGLikelihood, eventWeight );
    if( jet2_selected.Pt()>100. && jet2_selected.Pt()<123. ) h1_QGLikelihood_100_123->Fill( jet2_selected.QGLikelihood, eventWeight );
    if( jet1_selected.Pt()>66. && jet1_selected.Pt()<81. ) h1_QGLikelihood_66_81->Fill( jet1_selected.QGLikelihood, eventWeight );
    if( jet2_selected.Pt()>66. && jet2_selected.Pt()<81. ) h1_QGLikelihood_66_81->Fill( jet2_selected.QGLikelihood, eventWeight );


    if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood < 0.1 ) {
      h1_mZZ_kinfit_medMass_loQG->Fill(ZZ_kinfit.M(), eventWeight);
      h1_mZZ_kinfit_300Mass_loQG->Fill(ZZ_kinfit.M(), eventWeight);
      h1_mZZ_kinfit_hiMass_loQG->Fill(ZZ_kinfit.M(), eventWeight);
    } else {
      h1_mZZ_kinfit_medMass_hiQG->Fill(ZZ_kinfit.M(), eventWeight);
      h1_mZZ_kinfit_300Mass_hiQG->Fill(ZZ_kinfit.M(), eventWeight);
      h1_mZZ_kinfit_hiMass_hiQG->Fill(ZZ_kinfit.M(), eventWeight);
    }


    h1_mZZ_UL->Fill(ZZ_nokinfit.M(), eventWeight);
    h1_mZZ_hiMass->Fill(ZZ_nokinfit.M(), eventWeight);
    h1_mZZ_medMass->Fill(ZZ_nokinfit.M(), eventWeight);
    h1_mZZ_300Mass->Fill(ZZ_nokinfit.M(), eventWeight);
    h1_mZZ_ZjjMassConstr_hiMass->Fill(ZZ_constr.M(), eventWeight);
    h1_mZZ_ZjjMassConstr_medMass->Fill(ZZ_constr.M(), eventWeight);
    h1_mZZ_ZjjMassConstr_300Mass->Fill(ZZ_constr.M(), eventWeight);
    h1_mZZ_UL_kinfit->Fill(ZZ_kinfit.M(), eventWeight);
    h1_mZZ_kinfit_hiMass->Fill(ZZ_kinfit.M(), eventWeight);
    h1_mZZ_kinfit_medMass->Fill(ZZ_kinfit.M(), eventWeight);
    h1_mZZ_kinfit_300Mass->Fill(ZZ_kinfit.M(), eventWeight);
    h2_mZjj_vs_mZZ->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    h2_mZjj_vs_mZZ_kinfit->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M() );


    h1_deltaRZZ->Fill(Zjj_nokinfit.DeltaR(diLepton), eventWeight);

    h1_ptZZ->Fill( ZZ_nokinfit.Pt(), eventWeight );
    h1_ptZZ_kinfit->Fill( ZZ_kinfit.Pt(), eventWeight );
    h1_etaZZ->Fill( ZZ_nokinfit.Eta(), eventWeight );
    h1_etaZZ_kinfit->Fill( ZZ_kinfit.Eta(), eventWeight );

    h1_helicityLD->Fill( helicityLD_selected, eventWeight );

  //h1_cosThetaStar->Fill(hangles.helCosThetaStar, eventWeight);
  //h1_cosTheta1->Fill(hangles.helCosTheta1, eventWeight);
  //h1_cosTheta2->Fill(hangles.helCosTheta2, eventWeight);
  //h1_phi->Fill(hangles.helPhi, eventWeight);
  //h1_phi1->Fill(hangles.helPhi1, eventWeight);

  //h1_cosThetaStar_kinfit->Fill(hangles_kinfit.helCosThetaStar, eventWeight);
  //h1_cosTheta1_kinfit->Fill(hangles_kinfit.helCosTheta1, eventWeight);
  //h1_cosTheta2_kinfit->Fill(hangles_kinfit.helCosTheta2, eventWeight);
  //h1_phi_kinfit->Fill(hangles_kinfit.helPhi, eventWeight);
  //h1_phi1_kinfit->Fill(hangles_kinfit.helPhi1, eventWeight);



    //match to partons (any flavor):
    int partFlavor1=0;
    float deltaRmin1=999.;
    for(unsigned iPart=0; iPart<nPart; ++iPart ) {
      TLorentzVector thisPart;
      thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
      float thisDeltaR = jet1_selected.DeltaR(thisPart);
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
      float thisDeltaR = jet2_selected.DeltaR(thisPart);
      if( thisDeltaR<deltaRmin2 ) {
        partFlavor2 = pdgIdPart[iPart];
        deltaRmin2 = thisDeltaR;
      }
    }
    h1_deltaR_part2->Fill(deltaRmin2, eventWeight);
    h1_partFlavorJet2->Fill( partFlavor2, eventWeight );



  } //for entries


  std::cout << std::endl << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << "    DATASET: " << dataset_ << std::endl;
  std::cout << "----> PASSED SELECTION: " << 1000.*nEventsPassed_fb_kinfit << " ev/fb-1  (" << nEventsPassed_kinfit << " events)" << " Efficiency: " << 100.*nEventsPassed_fb_kinfit/nCounterW_ << "%" << std::endl;
  std::cout << "----> PASSED SELECTION (no kinfit): " << 1000.*nEventsPassed_fb_nokinfit << " ev/fb-1 (" << nEventsPassed_nokinfit << " events)" << std::endl;
  std::cout << std::endl;


//std::cout << "Passed chi square: " << nEvents_hiChiSquareProb << "/" << nEventsTot << " (" << 100.*nEvents_hiChiSquareProb/nEventsTot << "%)" << std::endl;
//std::cout << "Passed mZjj cut: " << nEvents_mZjj_cut << "/" << nEventsTot << " (" << 100.*nEvents_mZjj_cut/nEventsTot << "%)" << std::endl;



  outFile_->cd();

  h1_run->Write();
  h1_rhoPF_presel->Write();
  h1_rhoPF->Write();

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
  h1_mZjj_loChiSquareProb->Write();
  h1_mZjj_hiChiSquareProb->Write();
  h1_mZjj_all_presel->Write();

  h1_ptZll_presel->Write();
  h1_ptZjj_all_presel->Write();

  h1_deltaRjj->Write();
  h1_deltaRjj_prekin->Write();

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

  h1_ptZll_MW->Write();
  h1_ptJet1_MW->Write();
  h1_ptJet2_MW->Write();
  h1_deltaRjj_MW->Write();
  h1_QGLikelihoodProd_MW->Write();
  h1_mZjj_MW->Write();

  h1_mZZ_hiChiSquareProb->Write();
  h1_mZZ_loChiSquareProb->Write();
  h1_mZZ_mZjj_cut->Write();
  h1_mZZ_mZjj_notcut->Write();
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
  h1_mZZ_kinfit_300Mass_loQG->Write();
  h1_mZZ_kinfit_300Mass_hiQG->Write();
  h1_mZZ_kinfit_medMass->Write();
  h1_mZZ_kinfit_medMass_hiQG->Write();
  h1_mZZ_kinfit_medMass_loQG->Write();
  h1_mZZ_kinfit_hiMass->Write();
  h1_mZZ_kinfit_hiMass_hiQG->Write();
  h1_mZZ_kinfit_hiMass_loQG->Write();
//  h1_mZZ_kinfit_highestMass->Write();

  h1_ptZZ->Write();
  h1_ptZZ_kinfit->Write();
  h1_etaZZ->Write();
  h1_etaZZ_kinfit->Write();

  h1_deltaR_part1->Write();
  h1_ptJet1->Write();
  h1_ptJet1_prekin->Write();
  h1_eElectronsJet1->Write();
  h1_eMuonsJet1->Write();
  h1_partFlavorJet1->Write();

  h1_deltaR_part2->Write();
  h1_ptJet2->Write();
  h1_ptJet2_prekin->Write();
  h1_eElectronsJet2->Write();
  h1_eMuonsJet2->Write();
  h1_partFlavorJet2->Write();

  h1_deltaRZZ->Write();

//h1_mZZ_MCassoc->Write();
//h1_mZZ_MCassoc_ZjjMassConstr->Write();
//h1_mZZ_MCassoc_kinfit->Write();
//h1_mZZ_MCassoc_kinfit_cands->Write();

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

  h1_ptDJet1->Write();
  h1_ptDJet2->Write();

  h1_QGLikelihoodJet1->Write();
  h1_QGLikelihoodJet2->Write();
  h1_QGLikelihoodProd->Write();

  h1_QGLikelihood_100_123->Write();
  h1_QGLikelihood_66_81->Write();
//outFile_->mkdir("QGbins");
//outFile_->cd("QGbins");

//for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

//  vh1_rmsCandJet1[iPtBin]->Write();
//  vh1_ptDJet1[iPtBin]->Write();
//  vh1_nChargedJet1[iPtBin]->Write();
//  vh1_nNeutralJet1[iPtBin]->Write();
//  vh1_QGLikelihoodJet1[iPtBin]->Write();

//  vh1_rmsCandJet2[iPtBin]->Write();
//  vh1_ptDJet2[iPtBin]->Write();
//  vh1_nChargedJet2[iPtBin]->Write();
//  vh1_nNeutralJet2[iPtBin]->Write();
//  vh1_QGLikelihoodJet2[iPtBin]->Write();

//}


  outFile_->Close();


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
    mZZ_threshLo_ = 0.;
    mZZ_threshHi_ = 10000.;
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

  } else if( selectionType_=="loose_0btag" ) {

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
    requiredBTags_ = 0.;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt300" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 55.;
    ptJet2_thresh_ = 35.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 50.;
    mZll_threshHi_ = 120.;
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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt350" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 55.;
    ptJet2_thresh_ = 35.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 50.;
    mZll_threshHi_ = 120.;
    mZjj_threshLo_ = 81.;
    mZjj_threshHi_ = 101.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.7;
    ptZll_thresh_ = 60.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 315.;
    mZZ_threshHi_ = 385.;
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt250_FIXLD_0BTAGS" ) {

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
    helicityLD_thresh_ = 0.7;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 239.;
    mZZ_threshHi_ = 262.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt300_FIXLD_0BTAGS" ) {

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
    helicityLD_thresh_ = 0.7;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 290.;
    mZZ_threshHi_ = 324.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt350_FIXLD_0BTAGS" ) {

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
    helicityLD_thresh_ = 0.7;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 335.;
    mZZ_threshHi_ = 378.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt400_FIXLD_0BTAGS" ) {

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
    helicityLD_thresh_ = 0.7;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 375.;
    mZZ_threshHi_ = 456.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt450_FIXLD_0BTAGS" ) {

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
    helicityLD_thresh_ = 0.7;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 418.;
    mZZ_threshHi_ = 522.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt500_FIXLD_0BTAGS" ) {

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
    helicityLD_thresh_ = 0.7;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 467.;
    mZZ_threshHi_ = 600.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

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
    mZZ_threshLo_ = 540.;
    mZZ_threshHi_ = 660.;
    requiredBTags_ = -1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="3002l2b" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 76.;
    mZjj_threshHi_ = 106.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.8;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 270.;
    mZZ_threshHi_ = 330.;
    requiredBTags_ = 2;
    pfMetThresh_ = 35.;

  } else if( selectionType=="3502l2b" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 76.;
    mZjj_threshHi_ = 106.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.7;
    ptZll_thresh_ = 90.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 315.;
    mZZ_threshHi_ = 385.;
    requiredBTags_ = 2;
    pfMetThresh_ = 35.;

  } else if( selectionType=="4002l2b" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 76.;
    mZjj_threshHi_ = 106.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.6;
    ptZll_thresh_ = 110.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 360.;
    mZZ_threshHi_ = 440.;
    requiredBTags_ = 2;
    pfMetThresh_ = 0.;

  } else if( selectionType=="4502l2b" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 76.;
    mZjj_threshHi_ = 106.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.5;
    ptZll_thresh_ = 130.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 405.;
    mZZ_threshHi_ = 495.;
    requiredBTags_ = 2;
    pfMetThresh_ = 0.;

  } else if( selectionType=="5002l2b" ) {

    ptLept1_thresh_ = 20.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;
    mZjj_threshLo_ = 76.;
    mZjj_threshHi_ = 106.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.4;
    ptZll_thresh_ = 150.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 450.;
    mZZ_threshHi_ = 550.;
    requiredBTags_ = 2;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt250_0btag_optKin" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 59.7;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 1.5;
    ptZll_thresh_ = 56.2;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 239.;
    mZZ_threshHi_ = 262.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt300_0btag_optKin" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 65.6;
    ptJet2_thresh_ = 31.4;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 1.5;
    ptZll_thresh_ = 73.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 290.;
    mZZ_threshHi_ = 324.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt350_0btag_optKin" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 71.5;
    ptJet2_thresh_ = 38.8;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 1.5;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 89.9;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 335.;
    mZZ_threshHi_ = 378.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt400_0btag_optKin" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 77.4;
    ptJet2_thresh_ = 46.2;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 1.5;
    ptZll_thresh_ = 106.6;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 375.;
    mZZ_threshHi_ = 456.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt450_0btag_optKin" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 83.3;
    ptJet2_thresh_ = 53.6;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 1.5;
    ptZll_thresh_ = 123.4;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 418.;
    mZZ_threshHi_ = 522.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt500_0btag_optKin" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 89.2;
    ptJet2_thresh_ = 61.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    deltaRll_thresh_ = 9999.;
    deltaRjj_thresh_ = 1.5;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 140.2;
    helicityLD_thresh_ = 0.;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 467.;
    mZZ_threshHi_ = 600.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt250_0btag_optLD_QG06" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.392;
    else
      helicityLD_thresh_ = 0.459;
    QGLikelihoodProd_thresh_ = 0.06;
    mZZ_threshLo_ = 239.;
    mZZ_threshHi_ = 262.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt300_0btag_optLD_QG06" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.476;
    else
      helicityLD_thresh_ = 0.524;
    QGLikelihoodProd_thresh_ = 0.06;
    mZZ_threshLo_ = 290.;
    mZZ_threshHi_ = 324.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt350_0btag_optLD_QG06" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.56;
    else
      helicityLD_thresh_ = 0.589;
    QGLikelihoodProd_thresh_ = 0.06;
    mZZ_threshLo_ = 335.;
    mZZ_threshHi_ = 378.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt400_0btag_optLD_QG06" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.644;
    else
      helicityLD_thresh_ = 0.654;
    QGLikelihoodProd_thresh_ = 0.06;
    mZZ_threshLo_ = 375.;
    mZZ_threshHi_ = 456.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt450_0btag_optLD_QG06" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.728;
    else
      helicityLD_thresh_ = 0.719;
    QGLikelihoodProd_thresh_ = 0.06;
    mZZ_threshLo_ = 418.;
    mZZ_threshHi_ = 522.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt500_0btag_optLD_QG06" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.812;
    else
      helicityLD_thresh_ = 0.78;
    QGLikelihoodProd_thresh_ = 0.06;
    mZZ_threshLo_ = 467.;
    mZZ_threshHi_ = 600.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;
  } else if( selectionType=="opt250_0btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.453;
    else
      helicityLD_thresh_ = 0.459;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 239.;
    mZZ_threshHi_ = 262.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt300_0btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.515;
    else
      helicityLD_thresh_ = 0.524;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 290.;
    mZZ_threshHi_ = 324.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt350_0btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.577;
    else
      helicityLD_thresh_ = 0.589;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 335.;
    mZZ_threshHi_ = 378.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt400_0btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.639;
    else
      helicityLD_thresh_ = 0.654;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 375.;
    mZZ_threshHi_ = 456.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt400_0btag_optLD_noQG" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.639;
    else
      helicityLD_thresh_ = 0.654;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 375.;
    mZZ_threshHi_ = 456.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt450_0btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.701;
    else
      helicityLD_thresh_ = 0.719;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 418.;
    mZZ_threshHi_ = 522.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt500_0btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    if( useLooseBTags_ )
      helicityLD_thresh_ = 0.763;
    else
      helicityLD_thresh_ = 0.78;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 467.;
    mZZ_threshHi_ = 600.;
    requiredBTags_ = 0;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt250_1btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    helicityLD_thresh_ = 0.428;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 239.;
    mZZ_threshHi_ = 262.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt300_1btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    helicityLD_thresh_ = 0.474;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 290.;
    mZZ_threshHi_ = 324.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt350_1btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    helicityLD_thresh_ = 0.519;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 335.;
    mZZ_threshHi_ = 378.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt400_1btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    helicityLD_thresh_ = 0.565;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 375.;
    mZZ_threshHi_ = 456.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt450_1btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    helicityLD_thresh_ = 0.610;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 418.;
    mZZ_threshHi_ = 522.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="opt500_1btag_optLD" ) {

    ptLept1_thresh_ = 40.;
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
    helicityLD_thresh_ = 0.656;
    QGLikelihoodProd_thresh_ = 0.1;
    mZZ_threshLo_ = 467.;
    mZZ_threshHi_ = 600.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="250_1btag_LD055" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.55;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 225.;
    mZZ_threshHi_ = 275.;
    requiredBTags_ = 1;
    pfMetThresh_ = 35.;

  } else if( selectionType=="300_1btag_LD055" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.55;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 270.;
    mZZ_threshHi_ = 330.;
    requiredBTags_ = 1;
    pfMetThresh_ = 35.;

  } else if( selectionType=="350_1btag_LD055" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.55;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 315.;
    mZZ_threshHi_ = 385.;
    requiredBTags_ = 1;
    pfMetThresh_ = 35.;

  } else if( selectionType=="400_1btag_LD055" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.55;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 360.;
    mZZ_threshHi_ = 440.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="450_1btag_LD055" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.55;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 405.;
    mZZ_threshHi_ = 495.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="500_1btag_LD055" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.55;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 450.;
    mZZ_threshHi_ = 550.;
    requiredBTags_ = 1;
    pfMetThresh_ = 0.;

  } else if( selectionType=="250_2btag_LD05" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 225.;
    mZZ_threshHi_ = 275.;
    requiredBTags_ = 2;
    pfMetThresh_ = 35.;

  } else if( selectionType=="300_2btag_LD05" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 270.;
    mZZ_threshHi_ = 330.;
    requiredBTags_ = 2;
    pfMetThresh_ = 35.;

  } else if( selectionType=="350_2btag_LD05" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 315.;
    mZZ_threshHi_ = 385.;
    requiredBTags_ = 2;
    pfMetThresh_ = 35.;

  } else if( selectionType=="400_2btag_LD05" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 360.;
    mZZ_threshHi_ = 440.;
    requiredBTags_ = 2;
    pfMetThresh_ = 0.;

  } else if( selectionType=="450_2btag_LD05" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 405.;
    mZZ_threshHi_ = 495.;
    requiredBTags_ = 2;
    pfMetThresh_ = 0.;

  } else if( selectionType=="500_2btag_LD05" ) {

    ptLept1_thresh_ = 40.;
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
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;
    helicityLD_thresh_ = 0.5;
    QGLikelihoodProd_thresh_ = 0.;
    mZZ_threshLo_ = 450.;
    mZZ_threshHi_ = 550.;
    requiredBTags_ = 2;
    pfMetThresh_ = 0.;

  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType



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


/*
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
*/



// CODE FOR FULL PFCAND FIT:

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
