#include "Ntp1Finalizer_HZZlljj.h"

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

#include "QGLikelihoodCalculator.h"


#include "fitTools.h"


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

};


void print(TKinFitter *fitter);
Double_t ErrEt(Float_t Et, Float_t Eta);
Double_t ErrEta(Float_t Et, Float_t Eta);
Double_t ErrPhi(Float_t Et, Float_t Eta);
Double_t ErrEt(Float_t Et, Float_t Eta, int particleType);
Double_t ErrEta(Float_t Et, Float_t Eta, int particleType);
Double_t ErrPhi(Float_t Et, Float_t Eta, int particleType);

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


  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  if( dataset_=="Run2010B_runs146240_146733" )
    h1_totalLumi->SetBinContent(1, 1220000.);
  else
    h1_totalLumi->SetBinContent(1, totalLumi_);

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

  std::vector<TH1D*> vh1_ptDJet_all_presel = getHistoVector(nPtBins, ptBins, "ptDJet_all_presel", 50, 0., 1.);
  std::vector<TH1D*> vh1_rmsCandJet_all_presel = getHistoVector(nPtBins, ptBins, "rmsCandJet_all_presel", 50, 0., 0.07);
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

  TH1D* h1_ptJetBest1 = new TH1D("ptJetBest1", "", 27, 30., 400.);
  h1_ptJetBest1->Sumw2();

  TH1D* h1_ptJetRecoil = new TH1D("ptJetRecoil", "", 27, 30., 400.);
  h1_ptJetRecoil->Sumw2();
  TH1D* h1_ptDJetRecoil = new TH1D("ptDJetRecoil", "", 50, 0., 1.);
  h1_ptDJetRecoil->Sumw2();
  TH1D* h1_rmsCandJetRecoil = new TH1D("rmsCandJetRecoil", "", 50, 0., 0.07);
  h1_rmsCandJetRecoil->Sumw2();
  TH1D* h1_nChargedJetRecoil = new TH1D("nChargedJetRecoil", "", 41, -0.5, 40.5);
  h1_nChargedJetRecoil->Sumw2();
  TH1D* h1_nNeutralJetRecoil = new TH1D("nNeutralJetRecoil", "", 41, -0.5, 40.5);
  h1_nNeutralJetRecoil->Sumw2();

  int nBins_invMass = 40;
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
  TH1D* h1_partFlavorJet1= new TH1D("partFlavorJet1", "", 30, -7.5, 22.5);
  h1_partFlavorJet1->Sumw2();
  std::vector<TH1D*> vh1_ptDJet1 = getHistoVector(nPtBins, ptBins, "ptDJet1", 50, 0., 1.);
  std::vector<TH1D*> vh1_rmsCandJet1 = getHistoVector(nPtBins, ptBins, "rmsCandJet1", 50, 0., 0.1);
  std::vector<TH1D*> vh1_nChargedJet1 = getHistoVector(nPtBins, ptBins, "nChargedJet1", 51, -0.5, 50.5);
  std::vector<TH1D*> vh1_nNeutralJet1 = getHistoVector(nPtBins, ptBins, "nNeutralJet1", 51, -0.5, 50.5);
  std::vector<TH1D*> vh1_QGLikelihoodJet1 = getHistoVector(nPtBins, ptBins, "QGLikelihoodJet1", 50, 0., 1.);
  std::vector<TH1D*> vh1_QGLikelihood_normsJet1 = getHistoVector(nPtBins, ptBins, "QGLikelihood_normsJet1", 50, 0., 1.);
  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 50, 0., 1.);
  h1_QGLikelihoodJet1->Sumw2();
  TH1D* h1_QGLikelihood_normsJet1 = new TH1D("QGLikelihood_normsJet1", "", 50, 0., 1.);
  h1_QGLikelihood_normsJet1->Sumw2();
  TH1D* h1_QGLikelihoodJet1_eta2 = new TH1D("QGLikelihoodJet1_eta2", "", 50, 0., 1.);
  h1_QGLikelihoodJet1_eta2->Sumw2();


  TH1D* h1_deltaR_part2 = new TH1D("deltaR_part2", "", 50, 0., 0.8);
  h1_deltaR_part2->Sumw2();
  TH1D* h1_partFlavorJet2= new TH1D("partFlavorJet2", "", 30, -7.5, 22.5);
  h1_partFlavorJet2->Sumw2();

  std::vector<TH1D*> vh1_ptDJet2 = getHistoVector(nPtBins, ptBins, "ptDJet2", 50, 0., 1.);
  std::vector<TH1D*> vh1_rmsCandJet2 = getHistoVector(nPtBins, ptBins, "rmsCandJet2", 50, 0., 0.1);
  std::vector<TH1D*> vh1_nChargedJet2 = getHistoVector(nPtBins, ptBins, "nChargedJet2", 51, -0.5, 50.5);
  std::vector<TH1D*> vh1_nNeutralJet2 = getHistoVector(nPtBins, ptBins, "nNeutralJet2", 51, -0.5, 50.5);
  std::vector<TH1D*> vh1_QGLikelihoodJet2 = getHistoVector(nPtBins, ptBins, "QGLikelihoodJet2", 50, 0., 1.);
  std::vector<TH1D*> vh1_QGLikelihood_normsJet2 = getHistoVector(nPtBins, ptBins, "QGLikelihood_normsJet2", 50, 0., 1.);
  TH1D* h1_QGLikelihoodJet2 = new TH1D("QGLikelihoodJet2", "", 40, 0., 1.);
  h1_QGLikelihoodJet2->Sumw2();
  TH1D* h1_QGLikelihood_normsJet2 = new TH1D("QGLikelihood_normsJet2", "", 40, 0., 1.);
  h1_QGLikelihood_normsJet2->Sumw2();
  TH1D* h1_QGLikelihoodJet2_eta2 = new TH1D("QGLikelihoodJet2_eta2", "", 40, 0., 1.);
  h1_QGLikelihoodJet2_eta2->Sumw2();

  TH1D* h1_QGLikelihoodProd = new TH1D("QGLikelihoodProd", "", 50, 0., 2.);
  h1_QGLikelihoodProd->Sumw2();

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

  TH1D* h1_ZZInvMass_medMass= new TH1D("ZZInvMass_medMass", "", nBins_invMass, 150., 350.);
  h1_ZZInvMass_medMass->Sumw2();
  TH1D* h1_ZZInvMass_hiMass= new TH1D("ZZInvMass_hiMass", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass->Sumw2();
  TH1D* h1_ZZInvMass_kinfit_medMass= new TH1D("ZZInvMass_kinfit_medMass", "", nBins_invMass, 150., 350.);
  h1_ZZInvMass_kinfit_medMass->Sumw2();
  TH1D* h1_ZZInvMass_kinfit_hiMass= new TH1D("ZZInvMass_kinfit_hiMass", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_kinfit_hiMass->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_QGlikeli = new TH1D("ZZInvMass_hiMass_QGlikeli", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_QGlikeli->Sumw2();


  TH1D* h1_ZZInvMass_MCassoc  = new TH1D("ZZInvMass_MCassoc", "", 100, 200., 600.);
  h1_ZZInvMass_MCassoc->Sumw2();
  TH1D* h1_ZZInvMass_MCassoc_ZjjMassConstr  = new TH1D("ZZInvMass_MCassoc_ZjjMassConstr", "", 100, 200., 600.);
  h1_ZZInvMass_MCassoc_ZjjMassConstr->Sumw2();
  TH1D* h1_ZZInvMass_MCassoc_kinfit_jets  = new TH1D("ZZInvMass_MCassoc_kinfit_jets", "", 100, 200., 600.);
  h1_ZZInvMass_MCassoc_kinfit_jets->Sumw2();
  TH1D* h1_ZZInvMass_MCassoc_kinfit_cands  = new TH1D("ZZInvMass_MCassoc_kinfit_cands", "", 100, 200., 600.);
  h1_ZZInvMass_MCassoc_kinfit_cands->Sumw2();

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

  Float_t eLept2;
  tree_->SetBranchAddress("eLept2", &eLept2);
  Float_t ptLept2;
  tree_->SetBranchAddress("ptLept2", &ptLept2);
  Float_t etaLept2;
  tree_->SetBranchAddress("etaLept2", &etaLept2);
  Float_t phiLept2;
  tree_->SetBranchAddress("phiLept2", &phiLept2);

  Float_t eJetLead;
  tree_->SetBranchAddress("eJetLead", &eJetLead);
  Float_t ptJetLead;
  tree_->SetBranchAddress("ptJetLead", &ptJetLead);
  Float_t etaJetLead;
  tree_->SetBranchAddress("etaJetLead", &etaJetLead);
  Float_t phiJetLead;
  tree_->SetBranchAddress("phiJetLead", &phiJetLead);

  Float_t eJetLead2;
  tree_->SetBranchAddress("eJetLead2", &eJetLead2);
  Float_t ptJetLead2;
  tree_->SetBranchAddress("ptJetLead2", &ptJetLead2);
  Float_t etaJetLead2;
  tree_->SetBranchAddress("etaJetLead2", &etaJetLead2);
  Float_t phiJetLead2;
  tree_->SetBranchAddress("phiJetLead2", &phiJetLead2);

  Float_t eJetLead3;
  tree_->SetBranchAddress("eJetLead3", &eJetLead3);
  Float_t ptJetLead3;
  tree_->SetBranchAddress("ptJetLead3", &ptJetLead3);
  Float_t etaJetLead3;
  tree_->SetBranchAddress("etaJetLead3", &etaJetLead3);
  Float_t phiJetLead3;
  tree_->SetBranchAddress("phiJetLead3", &phiJetLead3);

  Float_t eJetBest1;
  tree_->SetBranchAddress("eJetBest1", &eJetBest1);
  Float_t ptJetBest1;
  tree_->SetBranchAddress("ptJetBest1", &ptJetBest1);
  Float_t etaJetBest1;
  tree_->SetBranchAddress("etaJetBest1", &etaJetBest1);
  Float_t phiJetBest1;
  tree_->SetBranchAddress("phiJetBest1", &phiJetBest1);
  Float_t rmsCandJetBest1;
  tree_->SetBranchAddress("rmsCandJetBest1", &rmsCandJetBest1);
  Float_t ptDJetBest1;
  tree_->SetBranchAddress("ptDJetBest1", &ptDJetBest1);
  Int_t nChargedJetBest1;
  tree_->SetBranchAddress("nChargedJetBest1", &nChargedJetBest1);
  Int_t nNeutralJetBest1;
  tree_->SetBranchAddress("nNeutralJetBest1", &nNeutralJetBest1);

  Float_t eJetBest2;
  tree_->SetBranchAddress("eJetBest2", &eJetBest2);
  Float_t ptJetBest2;
  tree_->SetBranchAddress("ptJetBest2", &ptJetBest2);
  Float_t etaJetBest2;
  tree_->SetBranchAddress("etaJetBest2", &etaJetBest2);
  Float_t phiJetBest2;
  tree_->SetBranchAddress("phiJetBest2", &phiJetBest2);
  Float_t rmsCandJetBest2;
  tree_->SetBranchAddress("rmsCandJetBest2", &rmsCandJetBest2);
  Float_t ptDJetBest2;
  tree_->SetBranchAddress("ptDJetBest2", &ptDJetBest2);
  Int_t nChargedJetBest2;
  tree_->SetBranchAddress("nChargedJetBest2", &nChargedJetBest2);
  Int_t nNeutralJetBest2;
  tree_->SetBranchAddress("nNeutralJetBest2", &nNeutralJetBest2);

  Float_t eJetRecoil;
  tree_->SetBranchAddress("eJetRecoil", &eJetRecoil);
  Float_t ptJetRecoil;
  tree_->SetBranchAddress("ptJetRecoil", &ptJetRecoil);
  Float_t etaJetRecoil;
  tree_->SetBranchAddress("etaJetRecoil", &etaJetRecoil);
  Float_t phiJetRecoil;
  tree_->SetBranchAddress("phiJetRecoil", &phiJetRecoil);
  Float_t rmsCandJetRecoil;
  tree_->SetBranchAddress("rmsCandJetRecoil", &rmsCandJetRecoil);
  Float_t ptDJetRecoil;
  tree_->SetBranchAddress("ptDJetRecoil", &ptDJetRecoil);
  Int_t nChargedJetRecoil;
  tree_->SetBranchAddress("nChargedJetRecoil", &nChargedJetRecoil);
  Int_t nNeutralJetRecoil;
  tree_->SetBranchAddress("nNeutralJetRecoil", &nNeutralJetRecoil);


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

    TLorentzVector jetLead;
    jetLead.SetPtEtaPhiE( ptJetLead, etaJetLead, phiJetLead, eJetLead );
    TLorentzVector jetLead2;
    jetLead2.SetPtEtaPhiE( ptJetLead2, etaJetLead2, phiJetLead2, eJetLead2 );

//  TLorentzVector jetBest1;
//  jetBest1.SetPtEtaPhiE( ptJetBest1, etaJetBest1, phiJetBest1, eJetBest1 );
//  TLorentzVector jetBest2;
//  jetBest2.SetPtEtaPhiE( ptJetBest2, etaJetBest2, phiJetBest2, eJetBest2 );

    TLorentzVector jetRecoil;
    jetRecoil.SetPtEtaPhiE( ptJetRecoil, etaJetRecoil, phiJetRecoil, eJetRecoil );

//  h1_ptJetBest1->Fill(ptJetBest1, eventWeight);
//  h1_ptJetBest2->Fill(ptJetBest2, eventWeight);

    h1_ptJetRecoil->Fill(ptJetRecoil, eventWeight);
//  if( ptJetRecoil>0. ) {
//  h1_rmsCandJetRecoil->Fill(rmsCandJetRecoil, eventWeight);
//  h1_ptDJetRecoil->Fill(ptDJetRecoil, eventWeight);
//  h1_nChargedJetRecoil->Fill(nChargedJetRecoil, eventWeight);
//  h1_nNeutralJetRecoil->Fill(nNeutralJetRecoil, eventWeight);
//  }



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

      jet2.rmsCand = rmsCandJet2[iJetPair];
      jet2.ptD = ptDJet2[iJetPair];
      jet2.nCharged = nChargedJet2[iJetPair];
      jet2.nNeutral = nNeutralJet2[iJetPair];

      TLorentzVector diJet = jet1 + jet2;

      if( jet1.Pt()>ptJet1_thresh_ && jet2.Pt()>ptJet2_thresh_ && fabs(jet1.Eta())<etaJet1_thresh_ && fabs(jet2.Eta())<etaJet1_thresh_ 
       && jet1.DeltaR(jet2) < deltaRjj_thresh_ && diJet.M() > mZjj_threshLo_ && diJet.M() < mZjj_threshHi_ && diJet.Pt() > ptZjj_thresh_ )
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
  


      // ------------------------
      //   KINEMATIC FIT: BEGIN
      // ------------------------


      TMatrixD m_jet1(3,3);
      TMatrixD m_jet2(3,3);

      m_jet1(0,0) = 0.5*ErrEt (jet1.Et(), jet1.Eta()); // et
      m_jet1(1,1) = 0.5*ErrEta(jet1.Et(), jet1.Eta()); // eta
      m_jet1(2,2) = 0.5*ErrPhi(jet1.Et(), jet1.Eta()); // phi
      m_jet2(0,0) = 0.5*ErrEt (jet2.Et(), jet2.Eta()); // et
      m_jet2(1,1) = 0.5*ErrEta(jet2.Et(), jet2.Eta()); // eta
      m_jet2(2,2) = 0.5*ErrPhi(jet2.Et(), jet2.Eta()); // phi

      TFitParticleEtEtaPhi *fitJet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &jet1, &m_jet1 );
      TFitParticleEtEtaPhi *fitJet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &jet2, &m_jet2 );
      
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


      TLorentzVector Zjj_constr;
      Zjj_constr.SetXYZM( bestZDiJet.Px(), bestZDiJet.Py(), bestZDiJet.Pz(), Zmass);


      TLorentzVector jet1_kinfit(*fitJet1->getCurr4Vec());
      TLorentzVector jet2_kinfit(*fitJet2->getCurr4Vec());
      TLorentzVector Zjj_kinfit_jets = jet1_kinfit + jet2_kinfit;

      TLorentzVector ZZ = bestZDiJet + diLepton;
      TLorentzVector ZZ_constr = diLepton + Zjj_constr;
      TLorentzVector ZZ_kinfit_jets = diLepton + Zjj_kinfit_jets;

      h2_mZjj_vs_mZZ->Fill( ZZ.M(), bestZDiJet.M() );
      h2_mZjj_vs_mZZ_kinfit->Fill( ZZ_kinfit_jets.M(), bestZDiJet.M() );


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
    //  h1_ZZInvMass_MCassoc->Fill( ZZ.M(), eventWeight );
    //  h1_ZZInvMass_MCassoc_ZjjMassConstr->Fill( ZZ_constr.M(), eventWeight );
    //  h1_ZZInvMass_MCassoc_kinfit_jets->Fill( ZZ_kinfit_jets.M(), eventWeight );
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
//      h1_ZZInvMass_MCassoc_kinfit_cands->Fill( ZZ_kinfit_cands.M(), eventWeight );
//    

//  } // if dataset


      // ------------------------
      //   KINEMATIC FIT: END
      // ------------------------




      if( lept1.Pt() > ptLept1_thresh_ && lept2.Pt() > ptLept2_thresh_ && fabs(lept1.Eta()) < etaLept1_thresh_ && fabs(lept2.Eta()) < etaLept2_thresh_
       && diLepton.M() > mZll_threshLo_ && diLepton.M() < mZll_threshHi_ && lept1.DeltaR(lept2) < deltaRll_thresh_ && diLepton.Pt() > ptZll_thresh_ ) {


        // event has passed kinematic selection: fill histograms:

        if( jet1.Pt()>jet2.Pt() ) {
          h1_ptJet1->Fill( jet1.Pt(), eventWeight );
          h1_ptJet2->Fill( jet2.Pt(), eventWeight );
        } else {
          h1_ptJet1->Fill( jet2.Pt(), eventWeight );
          h1_ptJet2->Fill( jet1.Pt(), eventWeight );
        }
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
        h1_ZZInvMass_hiMass->Fill(ZZ.M(), eventWeight);
        h1_ZZInvMass_medMass->Fill(ZZ.M(), eventWeight);
        h1_ZZInvMass_kinfit_hiMass->Fill(ZZ_kinfit_jets.M(), eventWeight);
        h1_ZZInvMass_kinfit_medMass->Fill(ZZ_kinfit_jets.M(), eventWeight);


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

        vh1_rmsCandJet1[jet1PtBin]->Fill( jet1.rmsCand, eventWeight );
        vh1_ptDJet1[jet1PtBin]->Fill( jet1.ptD, eventWeight );
        vh1_nChargedJet1[jet1PtBin]->Fill( jet1.nCharged, eventWeight );
        vh1_nNeutralJet1[jet1PtBin]->Fill( jet1.nNeutral, eventWeight );
        float QGLikelihoodJet1 = qglikeli->computeQGLikelihood( jet1.Pt(), jet1.nCharged, jet1.nNeutral, jet1.ptD, jet1.rmsCand );
        float QGLikelihood_normsJet1 = qglikeli->computeQGLikelihood( jet1.Pt(), jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );
        vh1_QGLikelihoodJet1[jet1PtBin]->Fill( QGLikelihoodJet1, eventWeight );
        h1_QGLikelihoodJet1->Fill( QGLikelihoodJet1, eventWeight );
        vh1_QGLikelihood_normsJet1[jet1PtBin]->Fill( QGLikelihood_normsJet1, eventWeight );
        h1_QGLikelihood_normsJet1->Fill( QGLikelihood_normsJet1, eventWeight );
        if( fabs(jet1.Eta())<2. ) h1_QGLikelihoodJet1_eta2->Fill(QGLikelihoodJet1, eventWeight);
        
        vh1_rmsCandJet2[jet2PtBin]->Fill( jet2.rmsCand, eventWeight );
        vh1_ptDJet2[jet2PtBin]->Fill( jet2.ptD, eventWeight );
        vh1_nChargedJet2[jet2PtBin]->Fill( jet2.nCharged, eventWeight );
        vh1_nNeutralJet2[jet2PtBin]->Fill( jet2.nNeutral, eventWeight );
        float QGLikelihoodJet2 = qglikeli->computeQGLikelihood( jet2.Pt(), jet2.nCharged, jet2.nNeutral, jet2.ptD, jet2.rmsCand );
        float QGLikelihood_normsJet2 = qglikeli->computeQGLikelihood( jet2.Pt(), jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );
        vh1_QGLikelihoodJet2[jet2PtBin]->Fill( QGLikelihoodJet2, eventWeight );
        h1_QGLikelihoodJet2->Fill( QGLikelihoodJet2, eventWeight );
        vh1_QGLikelihood_normsJet2[jet2PtBin]->Fill( QGLikelihood_normsJet2, eventWeight );
        h1_QGLikelihood_normsJet2->Fill( QGLikelihood_normsJet2, eventWeight );
        if( fabs(jet2.Eta())<2. ) h1_QGLikelihoodJet2_eta2->Fill(QGLikelihoodJet2, eventWeight);

        h1_QGLikelihoodProd->Fill( QGLikelihoodJet1*QGLikelihoodJet2, eventWeight );
        
        if( QGLikelihoodJet1<0.9 && QGLikelihoodJet2<0.9 )
          h1_ZZInvMass_hiMass_QGlikeli->Fill(ZZ_kinfit_jets.M(), eventWeight);


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


//std::string ofs400_name = "effTable400_tight_"+dataset_+".txt";
//ofstream ofs400(ofs400_name.c_str());
//ofs400 << "DATASET\tPreselection\tLepton pt\tLepton mass\tjet pt\tdijet mass\tjet deltaR" << std::endl;
//ofs400 << dataset_ << "\t"
//    << nEvents400_pre*1000. << "\t"
//    << nEvents400_pre_leptPt*1000. << "\t"
//    << nEvents400_pre_leptPt_leptMass*1000. << "\t"
//    << nEvents400_pre_leptPt_leptMass_jetPt*1000. << "\t"
//    << nEvents400_pre_leptPt_leptMass_jetPt_jetMass*1000. << "\t"
//    << nEvents400_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj*1000. << "\t"
//    << std::endl;

//std::string ofs500_name = "effTable500_tight_"+dataset_+".txt";
//ofstream ofs500(ofs500_name.c_str());
//ofs500 << "DATASET\tPreselection\tLepton pt\tLepton mass\tjet pt\tdijet mass\tjet deltaR" << std::endl;
//ofs500 << dataset_ << "\t"
//    << nEvents500_pre*1000. << "\t"
//    << nEvents500_pre_leptPt*1000. << "\t"
//    << nEvents500_pre_leptPt_leptMass*1000. << "\t"
//    << nEvents500_pre_leptPt_leptMass_jetPt*1000. << "\t"
//    << nEvents500_pre_leptPt_leptMass_jetPt_jetMass*1000. << "\t"
//    << nEvents500_pre_leptPt_leptMass_jetPt_jetMass_deltaRjj*1000. << "\t"
//    << std::endl;


  outFile_->cd();

  h1_totalLumi->Write();
  h1_run->Write();

  h1_ptJet_all_presel->Write();
  h1_etaJet_all_presel->Write();
  h1_nJets_presel->Write();
  h1_nPairs_presel->Write();

//h1_ptJetBest1->Write();
//h1_ptDJetBest1->Write();
//h1_rmsCandJetBest1->Write();
//h1_nChargedJetBest1->Write();
//h1_nNeutralJetBest1->Write();

//h1_ptJetBest2->Write();
//h1_ptDJetBest2->Write();
//h1_rmsCandJetBest2->Write();
//h1_nChargedJetBest2->Write();
//h1_nNeutralJetBest2->Write();

  h1_ptJetRecoil->Write();
  h1_ptDJetRecoil->Write();
  h1_rmsCandJetRecoil->Write();
  h1_nChargedJetRecoil->Write();
  h1_nNeutralJetRecoil->Write();

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

  h1_ZZInvMass_medMass->Write();
  h1_ZZInvMass_hiMass->Write();
  h1_ZZInvMass_kinfit_medMass->Write();
  h1_ZZInvMass_kinfit_hiMass->Write();
  h1_ZZInvMass_hiMass_QGlikeli->Write();

  h1_deltaR_part1->Write();
  h1_ptJet1->Write();
  h1_partFlavorJet1->Write();

  h1_deltaR_part2->Write();
  h1_ptJet2->Write();
  h1_partFlavorJet2->Write();


  h1_ZZInvMass_MCassoc->Write();
  h1_ZZInvMass_MCassoc_ZjjMassConstr->Write();
  h1_ZZInvMass_MCassoc_kinfit_jets->Write();
  h1_ZZInvMass_MCassoc_kinfit_cands->Write();

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
  h1_QGLikelihood_normsJet1->Write();
  h1_QGLikelihood_normsJet2->Write();
  h1_QGLikelihoodProd->Write();


  for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

    vh1_rmsCandJet_all_presel[iPtBin]->Write();
    vh1_ptDJet_all_presel[iPtBin]->Write();
    vh1_nChargedJet_all_presel[iPtBin]->Write();
    vh1_nNeutralJet_all_presel[iPtBin]->Write();

    vh1_rmsCandJet1[iPtBin]->Write();
    vh1_ptDJet1[iPtBin]->Write();
    vh1_nChargedJet1[iPtBin]->Write();
    vh1_nNeutralJet1[iPtBin]->Write();
    vh1_QGLikelihoodJet1[iPtBin]->Write();
    vh1_QGLikelihood_normsJet1[iPtBin]->Write();

    vh1_rmsCandJet2[iPtBin]->Write();
    vh1_ptDJet2[iPtBin]->Write();
    vh1_nChargedJet2[iPtBin]->Write();
    vh1_nNeutralJet2[iPtBin]->Write();
    vh1_QGLikelihoodJet2[iPtBin]->Write();
    vh1_QGLikelihood_normsJet2[iPtBin]->Write();

  }


  outFile_->Close();

} // finalize()



void Ntp1Finalizer_HZZlljj::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;

  if( selectionType_=="loose" ) {

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
    mZjj_threshHi_ = 120.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 999.;
    ptZll_thresh_ = 0.;
    ptZjj_thresh_ = 0.;

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

  } else if( selectionType=="opt400" ) {

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
    mZjj_threshLo_ = 81.;
    mZjj_threshHi_ = 101.;
    deltaRll_thresh_ = 999.;
    deltaRjj_thresh_ = 1.2;
    ptZll_thresh_ = 95.;
    ptZjj_thresh_ = 0.;

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
