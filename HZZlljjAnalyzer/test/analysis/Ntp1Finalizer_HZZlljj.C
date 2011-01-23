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


// constructor:

Ntp1Finalizer_HZZlljj::Ntp1Finalizer_HZZlljj( const std::string& dataset, const std::string& leptType ) : Ntp1Finalizer( "HZZlljj", dataset, leptType ) {

  if( leptType!="ALL" && leptType!="MU" && leptType!="ELE" ) {
    std::cout << "Lept type '" << leptType << "' currently not supported. Exiting." << std::endl;
    exit(9176);
  }

  leptType_ = leptType;

}





void Ntp1Finalizer_HZZlljj::finalize() {

  if( outFile_==0 ) this->createOutputFile();


  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  if( dataset_=="Run2010B_runs146240_146733" )
    h1_totalLumi->SetBinContent(1, 1220000.);
  else
    h1_totalLumi->SetBinContent(1, totalLumi_);

  TH1F* h1_run = new TH1F("run", "", 15149, 132440, 147589);

/*
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
*/

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
  TH1D* h1_ptDJet_all_presel = new TH1D("ptDJet_all_presel", "", 50, 0., 1.);
  h1_ptDJet_all_presel->Sumw2();
  TH1D* h1_rmsCandJet_all_presel = new TH1D("rmsCandJet_all_presel", "", 50, 0., 0.07);
  h1_rmsCandJet_all_presel->Sumw2();
  TH1D* h1_nChargedJet_all_presel = new TH1D("nChargedJet_all_presel", "", 41, -0.5, 40.5);
  h1_nChargedJet_all_presel->Sumw2();
  TH1D* h1_nNeutralJet_all_presel = new TH1D("nNeutralJet_all_presel", "", 41, -0.5, 40.5);
  h1_nNeutralJet_all_presel->Sumw2();
  TH1D* h1_nJets_presel = new TH1D("nJets_presel", "", 7, 1.5, 8.5);
  h1_nJets_presel->Sumw2();
  TH1D* h1_nPairs_presel = new TH1D("nPairs_presel", "", 21, 0.5, 21.5);
  h1_nPairs_presel->Sumw2();

  TH1D* h1_ptJetBest1 = new TH1D("ptJetBest1", "", 27, 30., 400.);
  h1_ptJetBest1->Sumw2();
  TH1D* h1_ptDJetBest1 = new TH1D("ptDJetBest1", "", 50, 0., 1.);
  h1_ptDJetBest1->Sumw2();
  TH1D* h1_rmsCandJetBest1 = new TH1D("rmsCandJetBest1", "", 50, 0., 0.07);
  h1_rmsCandJetBest1->Sumw2();
  TH1D* h1_nChargedJetBest1 = new TH1D("nChargedJetBest1", "", 41, -0.5, 40.5);
  h1_nChargedJetBest1->Sumw2();
  TH1D* h1_nNeutralJetBest1 = new TH1D("nNeutralJetBest1", "", 41, -0.5, 40.5);
  h1_nNeutralJetBest1->Sumw2();

  TH1D* h1_ptJetBest2 = new TH1D("ptJetBest2", "", 27, 30., 400.);
  h1_ptJetBest2->Sumw2();
  TH1D* h1_ptDJetBest2 = new TH1D("ptDJetBest2", "", 50, 0., 1.);
  h1_ptDJetBest2->Sumw2();
  TH1D* h1_rmsCandJetBest2 = new TH1D("rmsCandJetBest2", "", 50, 0., 0.07);
  h1_rmsCandJetBest2->Sumw2();
  TH1D* h1_nChargedJetBest2 = new TH1D("nChargedJetBest2", "", 41, -0.5, 40.5);
  h1_nChargedJetBest2->Sumw2();
  TH1D* h1_nNeutralJetBest2 = new TH1D("nNeutralJetBest2", "", 41, -0.5, 40.5);
  h1_nNeutralJetBest2->Sumw2();

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

  TH1D* h1_mZjj_presel = new TH1D("mZjj_presel", "", 20, invMassMin, 400.);
  h1_mZjj_presel->Sumw2();

  TH1D* h1_deltaRjj_presel = new TH1D("deltaRjj_presel", "", 18, 0.5, 5.);
  h1_deltaRjj_presel->Sumw2();
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

/*
  TH1D* h1_massZjj = new TH1D("massZjj", "", nBins_invMass, 20., 200.);
  h1_massZjj->Sumw2();
  TH1D* h1_massZjj_kinfit_jets = new TH1D("massZjj_kinfit_jets", "", nBins_invMass, 20., 200.);
  h1_massZjj_kinfit_jets->Sumw2();
  TH1D* h1_massZjj_MCassoc = new TH1D("massZjj_MCassoc", "", nBins_invMass, 20., 200.);
  h1_massZjj_MCassoc->Sumw2();
  TH1D* h1_massZjj_MCassoc_kinfit = new TH1D("massZjj_MCassoc_kinfit", "", nBins_invMass, 20., 200.);
  h1_massZjj_MCassoc_kinfit->Sumw2();
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
*/
/*
  TH1D* h1_ZZInvMass_loMass = new TH1D("ZZInvMass_loMass", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass->Sumw2();
  TH1D* h1_ZZInvMass_loMass_ZjjTag = new TH1D("ZZInvMass_loMass_ZjjTag", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass_ZjjTag->Sumw2();

  TH1D* h1_ZZInvMass_hiMass = new TH1D("ZZInvMass_hiMass", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_MCassoc_ZjjMassConstr = new TH1D("ZZInvMass_hiMass_MCassoc_ZjjMassConstr", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_MCassoc_ZjjMassConstr->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_kinfit_jets = new TH1D("ZZInvMass_hiMass_kinfit_jets", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_kinfit_jets->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_kinfit_cands = new TH1D("ZZInvMass_hiMass_kinfit_cands", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_kinfit_cands->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_MCassoc = new TH1D("ZZInvMass_hiMass_MCassoc", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_MCassoc->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_MCassoc_kinfit_jets = new TH1D("ZZInvMass_hiMass_MCassoc_kinfit_jets", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_MCassoc_kinfit_jets->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_MCassoc_kinfit_cands = new TH1D("ZZInvMass_hiMass_MCassoc_kinfit_cands", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_MCassoc_kinfit_cands->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_ZjjTag = new TH1D("ZZInvMass_hiMass_ZjjTag", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_ZjjTag->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_tightOLD = new TH1D("ZZInvMass_hiMass_fullSelection_tightOLD", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_fullSelection_tightOLD->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_tight = new TH1D("ZZInvMass_hiMass_fullSelection_tight", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_fullSelection_tight->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_tight_lead = new TH1D("ZZInvMass_hiMass_fullSelection_tight_lead", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_fullSelection_tight_lead->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_medium = new TH1D("ZZInvMass_hiMass_fullSelection_medium", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_fullSelection_medium->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr = new TH1D("ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_loose = new TH1D("ZZInvMass_hiMass_fullSelection_loose", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_fullSelection_loose->Sumw2();
*/


  TH1D* h1_ptJet1_loose = new TH1D("ptJet1_loose", "", 50, 30., 300.);
  h1_ptJet1_loose->Sumw2();
  TH1D* h1_ptJet1_tight = new TH1D("ptJet1_tight", "", 50, 30., 300.);
  h1_ptJet1_tight->Sumw2();

  TH1D* h1_ptJet2_loose = new TH1D("ptJet2_loose", "", 50, 30., 300.);
  h1_ptJet2_loose->Sumw2();
  TH1D* h1_ptJet2_tight = new TH1D("ptJet2_tight", "", 50, 30., 300.);
  h1_ptJet2_tight->Sumw2();

  TH1D* h1_mZll_loose = new TH1D("mZll_loose", "", nBins_invMass, 70., 120.);
  h1_mZll_loose->Sumw2();
  TH1D* h1_mZll_tight = new TH1D("mZll_tight", "", nBins_invMass, 70., 120.);
  h1_mZll_tight->Sumw2();
  TH1D* h1_mZll_opt400_LowEff = new TH1D("mZll_opt400_LowEff", "", nBins_invMass, 70., 120.);
  h1_mZll_opt400_LowEff->Sumw2();
  TH1D* h1_mZll_opt400_HighEff = new TH1D("mZll_opt400_HighEff", "", nBins_invMass, 70., 120.);
  h1_mZll_opt400_HighEff->Sumw2();
  TH1D* h1_mZll_opt500 = new TH1D("mZll_opt500", "", nBins_invMass, 70., 120.);
  h1_mZll_opt500->Sumw2();

  TH1D* h1_partFlavorJetOpt400_1 = new TH1D("partFlavorJetOpt400_1", "", 30, -7.5, 22.5);
  h1_partFlavorJetOpt400_1->Sumw2();
  TH1D* h1_rmsCandJetOpt400_1 = new TH1D("rmsCandJetOpt400_1", "", 50, 0., 0.07);
  h1_rmsCandJetOpt400_1->Sumw2();
  TH1D* h1_ptDJetOpt400_1 = new TH1D("ptDJetOpt400_1", "", 50, 0., 1.);
  h1_ptDJetOpt400_1->Sumw2();
  TH1D* h1_nChargedJetOpt400_1 = new TH1D("nChargedJetOpt400_1", "", 41, -0.5, 40.5);
  h1_nChargedJetOpt400_1->Sumw2();
  TH1D* h1_nNeutralJetOpt400_1 = new TH1D("nNeutralJetOpt400_1", "", 41, -0.5, 40.5);
  h1_nNeutralJetOpt400_1->Sumw2();

  TH1D* h1_partFlavorJetOpt400_2 = new TH1D("partFlavorJetOpt400_2", "", 30, -7.5, 22.5);
  h1_partFlavorJetOpt400_2->Sumw2();
  TH1D* h1_rmsCandJetOpt400_2 = new TH1D("rmsCandJetOpt400_2", "", 50, 0., 0.07);
  h1_rmsCandJetOpt400_2->Sumw2();
  TH1D* h1_ptDJetOpt400_2 = new TH1D("ptDJetOpt400_2", "", 50, 0., 1.);
  h1_ptDJetOpt400_2->Sumw2();
  TH1D* h1_nChargedJetOpt400_2 = new TH1D("nChargedJetOpt400_2", "", 41, -0.5, 40.5);
  h1_nChargedJetOpt400_2->Sumw2();
  TH1D* h1_nNeutralJetOpt400_2 = new TH1D("nNeutralJetOpt400_2", "", 41, -0.5, 40.5);
  h1_nNeutralJetOpt400_2->Sumw2();

  TH1D* h1_mZqq_loose = new TH1D("mZqq_loose", "", nBins_invMass, 70., 120.);
  h1_mZqq_loose->Sumw2();
  TH1D* h1_mZqq_tight = new TH1D("mZqq_tight", "", nBins_invMass, 70., 120.);
  h1_mZqq_tight->Sumw2();

  TH1D* h1_mZqq_opt400_LowEff = new TH1D("mZqq_opt400_LowEff", "", nBins_invMass, 70., 120.);
  h1_mZqq_opt400_LowEff->Sumw2();
  TH1D* h1_mZqq_opt400_HighEff = new TH1D("mZqq_opt400_HighEff", "", nBins_invMass, 70., 120.);
  h1_mZqq_opt400_HighEff->Sumw2();
  TH1D* h1_mZqq_opt500 = new TH1D("mZqq_opt500", "", nBins_invMass, 70., 120.);
  h1_mZqq_opt500->Sumw2();

  TH1D* h1_ZZInvMass_hiMass_loose = new TH1D("ZZInvMass_hiMass_loose", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_loose->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_loose_FINEBINNING = new TH1D("ZZInvMass_hiMass_loose_FINEBINNING", "", 800, 200., 1000.);
  h1_ZZInvMass_hiMass_loose_FINEBINNING->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_tight = new TH1D("ZZInvMass_hiMass_tight", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_tight->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_tight_FINEBINNING = new TH1D("ZZInvMass_hiMass_tight_FINEBINNING", "", 800, 200., 1000.);
  h1_ZZInvMass_hiMass_tight_FINEBINNING->Sumw2();

  TH1D* h1_ZZInvMass_hiMass_opt200= new TH1D("ZZInvMass_hiMass_opt200", "", nBins_invMass, 100., 300.);
  h1_ZZInvMass_hiMass_opt200->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_opt200_FINEBINNING = new TH1D("ZZInvMass_hiMass_opt200_FINEBINNING", "", 900, 100., 1000.);
  h1_ZZInvMass_hiMass_opt200_FINEBINNING->Sumw2();

  TH1D* h1_ZZInvMass_hiMass_opt300= new TH1D("ZZInvMass_hiMass_opt300", "", nBins_invMass, 200., 500.);
  h1_ZZInvMass_hiMass_opt300->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_opt300_FINEBINNING = new TH1D("ZZInvMass_hiMass_opt300_FINEBINNING", "", 900, 100., 1000.);
  h1_ZZInvMass_hiMass_opt300_FINEBINNING->Sumw2();

  TH1D* h1_ZZInvMass_hiMass_opt400= new TH1D("ZZInvMass_hiMass_opt400", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_opt400->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_opt400_FINEBINNING = new TH1D("ZZInvMass_hiMass_opt400_FINEBINNING", "", 900, 100., 1000.);
  h1_ZZInvMass_hiMass_opt400_FINEBINNING->Sumw2();

  TH1D* h1_ZZInvMass_hiMass_opt500= new TH1D("ZZInvMass_hiMass_opt500", "", nBins_invMass, 300., 700.);
  h1_ZZInvMass_hiMass_opt500->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_opt500_FINEBINNING = new TH1D("ZZInvMass_hiMass_opt500_FINEBINNING", "", 900, 100., 1000.);
  h1_ZZInvMass_hiMass_opt500_FINEBINNING->Sumw2();

  TH1D* h1_ZZInvMass_hiMass_opt600= new TH1D("ZZInvMass_hiMass_opt600", "", nBins_invMass, 400., 800.);
  h1_ZZInvMass_hiMass_opt600->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_opt600_FINEBINNING = new TH1D("ZZInvMass_hiMass_opt600_FINEBINNING", "", 900, 100., 1000.);
  h1_ZZInvMass_hiMass_opt600_FINEBINNING->Sumw2();

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

/*
  TH1D* h1_ZZInvMass_medMass = new TH1D("ZZInvMass_medMass", "", nBins_invMass, 100., 350.);
  h1_ZZInvMass_medMass->Sumw2();
  TH1D* h1_ZZInvMass_medMass_FINEBINNING = new TH1D("ZZInvMass_medMass_FINEBINNING", "", 20000, 100., 350.);
  h1_ZZInvMass_medMass_FINEBINNING->Sumw2();
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
  TH1D* h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag = new TH1D("ZZInvMass_hiMass_ZjjTag_ZllAntiTag", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag->Sumw2();

  TH1D* h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40 = new TH1D("ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40 = new TH1D("ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40", "", nBins_invMass, 200., 600.);
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40->Sumw2();
*/

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

  Double_t ptBins[16];
  fitTools::getBins_int( 16, ptBins, 20., 500.);

  TProfile* hp_ptJetGenMean = new TProfile("ptJetGenMean", "", 15, ptBins);
  std::vector<TH1F*> h1_response_vs_pt = getResponseHistos("response", 16, ptBins);
  std::vector<TH1F*> h1_response_vs_pt_Rch050 = getResponseHistos("response_Rch050", 16, ptBins);
  std::vector<TH1F*> h1_response_vs_pt_Rch5070 = getResponseHistos("response_Rch5070", 16, ptBins);
  std::vector<TH1F*> h1_response_vs_pt_Rch70100 = getResponseHistos("response_Rch70100", 16, ptBins);



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

//Int_t iJet1;
//tree_->SetBranchAddress("iJet1", &iJet1);
//Float_t eJet1;
//tree_->SetBranchAddress("eJet1", &eJet1);
//Float_t ptJet1;
//tree_->SetBranchAddress("ptJet1", &ptJet1);
//Float_t etaJet1;
//tree_->SetBranchAddress("etaJet1", &etaJet1);
//Float_t phiJet1;
//tree_->SetBranchAddress("phiJet1", &phiJet1);
//Float_t eChargedHadronsJet1;
//tree_->SetBranchAddress("eChargedHadronsJet1", &eChargedHadronsJet1);
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
//Float_t eJetGen1;
//tree_->SetBranchAddress("eJetGen1", &eJetGen1);
//Float_t ptJetGen1;
//tree_->SetBranchAddress("ptJetGen1", &ptJetGen1);
//Float_t etaJetGen1;
//tree_->SetBranchAddress("etaJetGen1", &etaJetGen1);
//Float_t phiJetGen1;
//tree_->SetBranchAddress("phiJetGen1", &phiJetGen1);

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
//Float_t eJetGen2;
//tree_->SetBranchAddress("eJetGen2", &eJetGen2);
//Float_t ptJetGen2;
//tree_->SetBranchAddress("ptJetGen2", &ptJetGen2);
//Float_t etaJetGen2;
//tree_->SetBranchAddress("etaJetGen2", &etaJetGen2);
//Float_t phiJetGen2;
//tree_->SetBranchAddress("phiJetGen2", &phiJetGen2);

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


  int nEntries = tree_->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


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

    TLorentzVector jetBest1;
    jetBest1.SetPtEtaPhiE( ptJetBest1, etaJetBest1, phiJetBest1, eJetBest1 );
    TLorentzVector jetBest2;
    jetBest2.SetPtEtaPhiE( ptJetBest2, etaJetBest2, phiJetBest2, eJetBest2 );

    TLorentzVector jetRecoil;
    jetRecoil.SetPtEtaPhiE( ptJetRecoil, etaJetRecoil, phiJetRecoil, eJetRecoil );

    h1_ptJetBest1->Fill(ptJetBest1, eventWeight);
    if( ptJetBest1>100. && ptJetBest1<150. ) {
    h1_rmsCandJetBest1->Fill(rmsCandJetBest1, eventWeight);
    h1_ptDJetBest1->Fill(ptDJetBest1, eventWeight);
    h1_nChargedJetBest1->Fill(nChargedJetBest1, eventWeight);
    h1_nNeutralJetBest1->Fill(nNeutralJetBest1, eventWeight);
    }

    h1_ptJetBest2->Fill(ptJetBest2, eventWeight);
    if( ptJetBest2>50. && ptJetBest2<80. ) {
    h1_rmsCandJetBest2->Fill(rmsCandJetBest2, eventWeight);
    h1_ptDJetBest2->Fill(ptDJetBest2, eventWeight);
    h1_nChargedJetBest2->Fill(nChargedJetBest2, eventWeight);
    h1_nNeutralJetBest2->Fill(nNeutralJetBest2, eventWeight);
    }

    h1_ptJetRecoil->Fill(ptJetRecoil, eventWeight);
    if( ptJetRecoil>0. ) {
    h1_rmsCandJetRecoil->Fill(rmsCandJetRecoil, eventWeight);
    h1_ptDJetRecoil->Fill(ptDJetRecoil, eventWeight);
    h1_nChargedJetRecoil->Fill(nChargedJetRecoil, eventWeight);
    h1_nNeutralJetRecoil->Fill(nNeutralJetRecoil, eventWeight);
    }

  //TLorentzVector jet1, jet2;
  //jet1.SetPtEtaPhiE( ptJet1, etaJet1, phiJet1, eJet1 );
  //jet2.SetPtEtaPhiE( ptJet2, etaJet2, phiJet2, eJet2 );


    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_looseSelection;
    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_tightSelection;
    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_opt200;
    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_opt300;
    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_opt400;
    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_opt500;
    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_opt600;

    if( nPairs>0 ) {

ofs << run << " " << event << " " << diLepton.M() << " " << ptLept1 << " " << etaLept1 << " " << ptLept2 << " " << etaLept2 << std::endl;
      h1_mZll_presel->Fill( diLepton.M(), eventWeight );
      if( leptType==0 )
        h1_mZmumu_presel->Fill( diLepton.M(), eventWeight );
      else
        h1_mZee_presel->Fill( diLepton.M(), eventWeight );

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

      if( jet1.Pt()>40. && jet2.Pt()>30. && fabs(jet1.Eta())<2.5 && fabs(jet2.Eta())<2.5 && diJet.M()>70. && diJet.M()<120. )
        jetPairs_looseSelection.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );

      if( jet1.Pt()>100. && jet2.Pt()>50. && jet1.DeltaR(jet2)<1.5 && diJet.M()>80. && diJet.M()<105. )
        jetPairs_tightSelection.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );


      if( jet1.Pt()>35. && jet2.Pt()>30. && jet1.DeltaR(jet2)<2.8 && diJet.M()>81. && diJet.M()<101. )
        jetPairs_opt200.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );

      if( jet1.Pt()>45. && jet2.Pt()>35. && jet1.DeltaR(jet2)<1.7 && diJet.M()>81. && diJet.M()<101. )
        jetPairs_opt300.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );

      if( jet1.Pt()>90. && jet2.Pt()>55. && jet1.DeltaR(jet2)<1.2 && diJet.M()>81. && diJet.M()<101. )
        jetPairs_opt400.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );

      if( jet1.Pt()>125. && jet2.Pt()>95. && jet1.DeltaR(jet2)<1.5 && diJet.M()>81. && diJet.M()<101. )
        jetPairs_opt500.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );

      if( jet1.Pt()>160. && jet2.Pt()>95. && jet1.DeltaR(jet2)<1.0 && diJet.M()>81. && diJet.M()<101. )
        jetPairs_opt600.push_back( std::pair<AnalysisJet,AnalysisJet>(jet1,jet2) );


      h1_mZjj_presel->Fill( diJet.M(), eventWeight );
      h1_deltaRjj_presel->Fill(jet1.DeltaR(jet2), eventWeight );

      if( jet1.Pt()!=cached_jetpt ) {
        h1_ptJet_all_presel->Fill( jet1.Pt(), eventWeight );
        h1_etaJet_all_presel->Fill( jet1.Eta(), eventWeight );
        if( jet1.Pt()>100. && jet1.Pt()<150. ) {
        h1_rmsCandJet_all_presel->Fill( rmsCandJet1[iJetPair], eventWeight );
        h1_ptDJet_all_presel->Fill( ptDJet1[iJetPair], eventWeight );
        h1_nChargedJet_all_presel->Fill( nChargedJet1[iJetPair], eventWeight );
        h1_nNeutralJet_all_presel->Fill( nNeutralJet1[iJetPair], eventWeight );
        }
        cached_jetpt = jet1.Pt();
      }

    }


    TLorentzVector jet1_presel, jet2_presel;

    if( jetPairs_looseSelection.size()>0. ) {

      // now look for best Z mass jet pair:
      float Zmass = 91.19;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_looseSelection.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_looseSelection[iPair].first + jetPairs_looseSelection[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs

      // jet1_presel and jet2_presel are used in the kinfitter
      jet1_presel = jetPairs_looseSelection[bestPair].first;
      jet2_presel = jetPairs_looseSelection[bestPair].second;

      TLorentzVector bestZDiJet = jetPairs_looseSelection[bestPair].first + jetPairs_looseSelection[bestPair].second;
      TLorentzVector ZZ_loose = diLepton + bestZDiJet; 
  


      if( lept1.Pt()>10. && lept2.Pt()>10. && diLepton.M() > 70. && diLepton.M() < 110. ) {
        if( jetPairs_looseSelection[bestPair].first.Pt()>jetPairs_looseSelection[bestPair].second.Pt() ) {
          h1_ptJet1_loose->Fill( jetPairs_looseSelection[bestPair].first.Pt(), eventWeight );
          h1_ptJet2_loose->Fill( jetPairs_looseSelection[bestPair].second.Pt(), eventWeight );
        } else {
          h1_ptJet1_loose->Fill( jetPairs_looseSelection[bestPair].second.Pt(), eventWeight );
          h1_ptJet2_loose->Fill( jetPairs_looseSelection[bestPair].first.Pt(), eventWeight );
        }
        h1_mZll_loose->Fill( diLepton.M(), eventWeight);
        h1_mZqq_loose->Fill( bestZDiJet.M(), eventWeight);
        //h2_mZjj_vs_mZZ->Fill( ZZ_loose.M(), bestZDiJet.M(), eventWeight );
        h1_ZZInvMass_hiMass_loose->Fill(ZZ_loose.M(), eventWeight);
        h1_ZZInvMass_hiMass_loose_FINEBINNING->Fill(ZZ_loose.M(), eventWeight);
      }

    }

    if( jetPairs_tightSelection.size()>0. ) {

      // now look for best Z mass jet pair:
      float Zmass = 91.19;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_tightSelection.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_tightSelection[iPair].first + jetPairs_tightSelection[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs

      TLorentzVector jet1_tight = jetPairs_tightSelection[bestPair].first;
      TLorentzVector jet2_tight = jetPairs_tightSelection[bestPair].second;

      TLorentzVector bestZDiJet = jet1_tight + jet2_tight;
      TLorentzVector ZZ_tight = diLepton + bestZDiJet; 
  
      if( lept1.Pt()>100. && lept2.Pt()>50. && diLepton.M() > 86. && diLepton.M() < 96. ) {
        if( jet1_tight.Pt()>jet2_tight.Pt() ) {
          h1_ptJet1_tight->Fill( jet1_tight.Pt(), eventWeight );
          h1_ptJet2_tight->Fill( jet2_tight.Pt(), eventWeight );
        } else {
          h1_ptJet1_tight->Fill( jet1_tight.Pt(), eventWeight );
          h1_ptJet2_tight->Fill( jet2_tight.Pt(), eventWeight );
        }
        h1_mZll_tight->Fill( diLepton.M(), eventWeight);
        h1_mZqq_tight->Fill( bestZDiJet.M(), eventWeight);
        h1_ZZInvMass_hiMass_tight->Fill(ZZ_tight.M(), eventWeight);
        h1_ZZInvMass_hiMass_tight_FINEBINNING->Fill(ZZ_tight.M(), eventWeight);

/*
        // match to parton:
        int partFlavor=0;
        float deltaRmin=999.;
        for(unsigned iPart=0; iPart<nPart; ++iPart ) {
          TLorentzVector thisPart;
          thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
          float thisDeltaR = jet1.DeltaR(thisPart);
          if( thisDeltaR<deltaRmin ) {
            partFlavor = pdgIdPart[iPart];
            deltaRmin = thisDeltaR;
          }
        }
        if( deltaRmin<0.5 && partFlavor!=0 )
          h1_partFlavor_tight->Fill(partFlavor);
if( partFlavor==21 ) std::cout << deltaRmin << std::endl;

        deltaRmin=999.;
        partFlavor=0;
        for(unsigned iPart=0; iPart<nPart; ++iPart ) {
          TLorentzVector thisPart;
          thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
          float thisDeltaR = jet2.DeltaR(thisPart);
          if( thisDeltaR<deltaRmin ) {
            partFlavor = pdgIdPart[iPart];
            deltaRmin = thisDeltaR;
          }
        }
        if( deltaRmin<0.5 && partFlavor !=0 )
          h1_partFlavor_tight->Fill(partFlavor);
if( partFlavor==21 ) std::cout << deltaRmin << std::endl;

*/
      }

    }


    if( jetPairs_opt200.size()>0. ) {

      // now look for best Z mass jet pair:
      float Zmass = 91.19;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_opt200.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_opt200[iPair].first + jetPairs_opt200[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs

      TLorentzVector jet1 = jetPairs_opt200[bestPair].first;
      TLorentzVector jet2 = jetPairs_opt200[bestPair].second;

      TLorentzVector bestZDiJet = jet1 + jet2;
      TLorentzVector ZZ_opt200= diLepton + bestZDiJet; 
      
  
      if( lept1.Pt()>40. && lept1.DeltaR(lept2)<2.8&& diLepton.M() > 86. && diLepton.M() < 96. ) {
        h1_ZZInvMass_hiMass_opt200->Fill(ZZ_opt200.M(), eventWeight);
        h1_ZZInvMass_hiMass_opt200->Fill(ZZ_opt200.M(), eventWeight);
      }

    }

    if( jetPairs_opt300.size()>0. ) {

      // now look for best Z mass jet pair:
      float Zmass = 91.19;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_opt300.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_opt300[iPair].first + jetPairs_opt300[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs

      TLorentzVector jet1 = jetPairs_opt300[bestPair].first;
      TLorentzVector jet2 = jetPairs_opt300[bestPair].second;

      TLorentzVector bestZDiJet = jet1 + jet2;
      TLorentzVector ZZ_opt300= diLepton + bestZDiJet; 
      
  
      if( lept1.Pt()>40. && lept1.DeltaR(lept2)<2.1 && diLepton.M() > 86. && diLepton.M() < 96. ) {
        h1_ZZInvMass_hiMass_opt300->Fill(ZZ_opt300.M(), eventWeight);
        h1_ZZInvMass_hiMass_opt300->Fill(ZZ_opt300.M(), eventWeight);
      }

    }


    if( jetPairs_opt400.size()>0. ) {

      // now look for best Z mass jet pair:
      float Zmass = 91.19;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_opt400.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_opt400[iPair].first + jetPairs_opt400[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs

      AnalysisJet jet1 = jetPairs_opt400[bestPair].first;
      AnalysisJet jet2 = jetPairs_opt400[bestPair].second;

      TLorentzVector bestZDiJet = jet1 + jet2;
      TLorentzVector ZZ_opt400= diLepton + bestZDiJet; 
      
  
      if( lept1.Pt()>40. && diLepton.Pt()>95. && diLepton.M() > 86. && diLepton.M() < 96. ) {
        h1_ZZInvMass_hiMass_opt400->Fill(ZZ_opt400.M(), eventWeight);
        h1_ZZInvMass_hiMass_opt400->Fill(ZZ_opt400.M(), eventWeight);
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
        if( deltaRmin1<0.5 ) {
          h1_partFlavorJetOpt400_1->Fill( partFlavor1, eventWeight );
          h1_rmsCandJetOpt400_1->Fill( jet1.rmsCand, eventWeight );
          h1_ptDJetOpt400_1->Fill( jet1.ptD, eventWeight );
          h1_nChargedJetOpt400_1->Fill( jet1.nCharged, eventWeight );
          h1_nNeutralJetOpt400_1->Fill( jet1.nNeutral, eventWeight );
        }
      
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
        if( deltaRmin2<0.5 ) {
          h1_partFlavorJetOpt400_2->Fill( partFlavor2, eventWeight );
          h1_rmsCandJetOpt400_2->Fill( jet2.rmsCand, eventWeight );
          h1_ptDJetOpt400_2->Fill( jet2.ptD, eventWeight );
          h1_nChargedJetOpt400_2->Fill( jet2.nCharged, eventWeight );
          h1_nNeutralJetOpt400_2->Fill( jet2.nNeutral, eventWeight );
        }
      }

    }


    if( jetPairs_opt500.size()>0. ) {

      // now look for best Z mass jet pair:
      float Zmass = 91.19;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_opt500.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_opt500[iPair].first + jetPairs_opt500[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs

      TLorentzVector jet1 = jetPairs_opt500[bestPair].first;
      TLorentzVector jet2 = jetPairs_opt500[bestPair].second;

      TLorentzVector bestZDiJet = jet1 + jet2;
      TLorentzVector ZZ_opt500 = diLepton + bestZDiJet; 
      
      if( lept1.Pt()>40. && lept1.DeltaR(lept2)<1.4 && diLepton.M() > 86. && diLepton.M() < 96. ) {
        h1_mZll_opt500->Fill( diLepton.M(), eventWeight);
        h1_mZqq_opt500->Fill( bestZDiJet.M(), eventWeight);
        h1_ZZInvMass_hiMass_opt500->Fill(ZZ_opt500.M(), eventWeight);
        h1_ZZInvMass_hiMass_opt500_FINEBINNING->Fill(ZZ_opt500.M(), eventWeight);
      }


    }

    

    if( jetPairs_opt600.size()>0. ) {

      // now look for best Z mass jet pair:
      float Zmass = 91.19;
      float bestMass = 0.;
      int bestPair=-1;
      
      for( unsigned iPair=0; iPair<jetPairs_opt600.size(); ++iPair ) {

        TLorentzVector dijet = jetPairs_opt600[iPair].first + jetPairs_opt600[iPair].second;
        float invMass = dijet.M();
        if( bestPair==-1 || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
          bestMass = invMass;
          bestPair = iPair;
        }

      } //for pairs

      TLorentzVector jet1 = jetPairs_opt600[bestPair].first;
      TLorentzVector jet2 = jetPairs_opt600[bestPair].second;

      TLorentzVector bestZDiJet = jet1 + jet2;
      TLorentzVector ZZ_opt600 = diLepton + bestZDiJet; 
      
      if( lept1.Pt()>40. && lept1.DeltaR(lept2)<1.2) {
        h1_ZZInvMass_hiMass_opt600->Fill(ZZ_opt600.M(), eventWeight);
        h1_ZZInvMass_hiMass_opt600_FINEBINNING->Fill(ZZ_opt600.M(), eventWeight);
      }


    }

    

/*
    h1_ptLept1->Fill( lept1.Pt(), eventWeight );
    h1_ptLept2->Fill( lept2.Pt(), eventWeight );

    h1_ptLept2OverLept1->Fill( lept2.Pt()/lept1.Pt(), eventWeight );

    float deltaRll = lept1.DeltaR(lept2);
    //h1_deltaRll->Fill( deltaRll, eventWeight );


    TLorentzVector Zll = ( lept1 + lept2 );
    h1_ptZll->Fill( Zll.Pt(), eventWeight );
    h1_pzZll->Fill( Zll.Pz(), eventWeight );



    bool jet1_is_lead = false;
    if( iJet1==0 ) jet1_is_lead=true;
    float Rch1 = eChargedHadronsJet1/eJet1;
    float Rch2 = eChargedHadronsJet2/eJet2;
  //if( fabs(jet1.Eta())<2.5 && Rch1<0.1 ) continue;
  //if( fabs(jet2.Eta())<2.5 && Rch2<0.1 ) continue;
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


    TLorentzVector Zjj_lead = jetLead + jetLead2;

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
    TLorentzVector ZZ_lead = Zjj_lead + Zll;

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
*/

    // ------------------------
    //   KINEMATIC FIT: BEGIN
    // ------------------------

    if( jet1_presel.Pt()>30. && jet2_presel.Pt()>30. ) {

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

      TMatrixD m_jet1(3,3);
      TMatrixD m_jet2(3,3);

      m_jet1(0,0) = 0.5*ErrEt (jet1_presel.Et(), jet1_presel.Eta()); // et
      m_jet1(1,1) = 0.5*ErrEta(jet1_presel.Et(), jet1_presel.Eta()); // eta
      m_jet1(2,2) = 0.5*ErrPhi(jet1_presel.Et(), jet1_presel.Eta()); // phi
      m_jet2(0,0) = 0.5*ErrEt (jet2_presel.Et(), jet2_presel.Eta()); // et
      m_jet2(1,1) = 0.5*ErrEta(jet2_presel.Et(), jet2_presel.Eta()); // eta
      m_jet2(2,2) = 0.5*ErrPhi(jet2_presel.Et(), jet2_presel.Eta()); // phi

      TFitParticleEtEtaPhi *fitJet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &jet1_presel, &m_jet1 );
      TFitParticleEtEtaPhi *fitJet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &jet2_presel, &m_jet2 );
      
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


      TLorentzVector Zjj = jet1_presel + jet2_presel;
      TLorentzVector Zjj_constr;
      Zjj_constr.SetXYZM( Zjj.Px(), Zjj.Py(), Zjj.Pz(), 91.1876);


      TLorentzVector jet1_kinfit(*fitJet1->getCurr4Vec());
      TLorentzVector jet2_kinfit(*fitJet2->getCurr4Vec());
      TLorentzVector Zjj_kinfit_jets = jet1_kinfit + jet2_kinfit;

      TLorentzVector Zll = lept1+lept2;

      TLorentzVector ZZ = Zjj + Zll;
      TLorentzVector ZZ_constr = Zll + Zjj_constr;
      TLorentzVector ZZ_kinfit_jets = Zll + Zjj_kinfit_jets;

      h2_mZjj_vs_mZZ->Fill( ZZ.M(), Zjj.M() );
      h2_mZjj_vs_mZZ_kinfit->Fill( ZZ_kinfit_jets.M(), Zjj.M() );


      // match to parton:
      int partFlavor1=0;
      float deltaRmin1=999.;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = jet1_presel.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin1 ) {
          partFlavor1 = pdgIdPart[iPart];
          deltaRmin1 = thisDeltaR;
        }
      }

      float deltaRmin2=999.;
      int partFlavor2=0;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = jet2_presel.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin2 ) {
          partFlavor2 = pdgIdPart[iPart];
          deltaRmin2 = thisDeltaR;
        }
      }

      bool bothMatched = ( deltaRmin1<0.5 && deltaRmin2<0.5 && partFlavor1!=0 && partFlavor2!=0 );

      if( bothMatched ) {
        h1_ZZInvMass_MCassoc->Fill( ZZ.M(), eventWeight );
        h1_ZZInvMass_MCassoc_ZjjMassConstr->Fill( ZZ_constr.M(), eventWeight );
        h1_ZZInvMass_MCassoc_kinfit_jets->Fill( ZZ_kinfit_jets.M(), eventWeight );
      }

      // and now full kinematic fit with PFCands:
      TLorentzVector ZZ_kinfit_cands;
      TRegexp cands_tstr("CANDS");
      TString dataset_tstr(dataset_);
  if( dataset_tstr.Contains(cands_tstr) ) {

      std::vector<TFitParticleEtEtaPhi*> fitCands;
      std::vector<int> fitCandTypes;
      TFitConstraintM *mCons_cands = new TFitConstraintM( "ZMassConstraint_cands", "ZMass-Constraint", 0, 0 , 91.19);


      float testFactor = 1.;

      // loop on PFCands of first jet
      TLorentzVector candJet1(0., 0., 0., 0.);
      for( unsigned iPFCand=0; iPFCand<nPFCand1; ++iPFCand ) {
        TLorentzVector thisCand;
        thisCand.SetPtEtaPhiE( ptPFCand1[iPFCand], etaPFCand1[iPFCand], phiPFCand1[iPFCand], ePFCand1[iPFCand] );
        char candName[100];
        sprintf( candName, "PFCand_1_%d", iPFCand );

        if( particleTypePFCand1[iPFCand]==0 ) {
          std::cout << "FOUND PARTICLE TYPE=0!!!! SKIPPING!" << std::endl;
          continue;
        }

        if( particleTypePFCand1[iPFCand]==5 ) // neutral hadrons
          testFactor=10.;
        else
          testFactor=1.;

        TMatrixD m_PFCand(3,3);
        m_PFCand(0,0) = testFactor*ErrEt ( thisCand.Et(), thisCand.Eta(), particleTypePFCand1[iPFCand] );
        m_PFCand(1,1) = testFactor*ErrEta( thisCand.Et(), thisCand.Eta(), particleTypePFCand1[iPFCand] );
        m_PFCand(2,2) = testFactor*ErrPhi( thisCand.Et(), thisCand.Eta(), particleTypePFCand1[iPFCand] );
        TFitParticleEtEtaPhi* fitCand = new TFitParticleEtEtaPhi( candName, candName, &thisCand, &m_PFCand );

        mCons_cands->addParticle1( fitCand );
        fitCands.push_back(fitCand);
        fitCandTypes.push_back(particleTypePFCand1[iPFCand]);

        candJet1 += thisCand;

      }

      TLorentzVector cand_add1 = jet1_presel - candJet1;
      // make it worse-resolution case:
      int particleType_add1 = (fabs(cand_add1.Eta())<3.) ? 5 : 6; 
      TMatrixD m_PFCand_add1(3,3);
      m_PFCand_add1(0,0) = testFactor*ErrEt ( cand_add1.Et(), cand_add1.Eta(), particleType_add1 );
      m_PFCand_add1(1,1) = testFactor*ErrEta( cand_add1.Et(), cand_add1.Eta(), particleType_add1 );
      m_PFCand_add1(2,2) = testFactor*ErrPhi( cand_add1.Et(), cand_add1.Eta(), particleType_add1 );
      TFitParticleEtEtaPhi* fitCand_add1 = new TFitParticleEtEtaPhi( "PFCand_add1", "PFCand_add1", &cand_add1, &m_PFCand_add1 );

    //mCons_cands->addParticle1( fitCand_add1 );
    //fitCands.push_back(fitCand_add1);
    //fitCandTypes.push_back(particleType_add1 );


      // loop on PFCands of second jet
      TLorentzVector candJet2(0., 0., 0., 0.);
      for( unsigned iPFCand=0; iPFCand<nPFCand2; ++iPFCand ) {
        TLorentzVector thisCand;
        thisCand.SetPtEtaPhiE( ptPFCand2[iPFCand], etaPFCand2[iPFCand], phiPFCand2[iPFCand], ePFCand2[iPFCand] );
        char candName[100];
        sprintf( candName, "PFCand_2_%d", iPFCand );

        if( particleTypePFCand2[iPFCand]==0 ) {
          std::cout << "FOUND PARTICLE TYPE=0!!!! SKIPPING!" << std::endl;
          continue;
        }

        if( particleTypePFCand1[iPFCand]==5 ) // neutral hadrons
          testFactor=10.;
        else
          testFactor=1.;

        TMatrixD m_PFCand(3,3);
        m_PFCand(0,0) = testFactor*ErrEt ( thisCand.Et(), thisCand.Eta(), particleTypePFCand2[iPFCand] );
        m_PFCand(1,1) = testFactor*ErrEta( thisCand.Et(), thisCand.Eta(), particleTypePFCand2[iPFCand] );
        m_PFCand(2,2) = testFactor*ErrPhi( thisCand.Et(), thisCand.Eta(), particleTypePFCand2[iPFCand] );
        TFitParticleEtEtaPhi* fitCand = new TFitParticleEtEtaPhi( candName, candName, &thisCand, &m_PFCand );

        mCons_cands->addParticle1( fitCand );
        fitCands.push_back(fitCand);
        fitCandTypes.push_back(particleTypePFCand2[iPFCand]);

        candJet2 += thisCand;

      }


      TLorentzVector cand_add2 = jet2_presel - candJet2;
      // make it worse-resolution case:
      int particleType_add2 = (fabs(cand_add2.Eta())<3.) ? 5 : 6; 
      TMatrixD m_PFCand_add2(3,3);
      m_PFCand_add2(0,0) = testFactor*ErrEt ( cand_add2.Et(), cand_add2.Eta(), particleType_add2 );
      m_PFCand_add2(1,1) = testFactor*ErrEta( cand_add2.Et(), cand_add2.Eta(), particleType_add2 );
      m_PFCand_add2(2,2) = testFactor*ErrPhi( cand_add2.Et(), cand_add2.Eta(), particleType_add2 );
      TFitParticleEtEtaPhi* fitCand_add2 = new TFitParticleEtEtaPhi( "PFCand_add2", "PFCand_add2", &cand_add2, &m_PFCand_add2 );

    //mCons_cands->addParticle2( fitCand_add2 );
    //fitCands.push_back(fitCand_add2);
    //fitCandTypes.push_back(particleType_add2 );

      TKinFitter* fitter_cands = new TKinFitter("fitter_cands", "fitter_cands");
      for( unsigned iCand=0; iCand<fitCands.size(); ++iCand )
        fitter_cands->addMeasParticle( fitCands[iCand] );
      fitter_cands->addConstraint( mCons_cands );

      //Set convergence criteria
      fitter_cands->setMaxNbIter( 30 );
      fitter_cands->setMaxDeltaS( 1e-2 );
      fitter_cands->setMaxF( 1e-1 );
      fitter_cands->setVerbosity(0);

      //Perform the fit
      fitter_cands->fit();

      // recreate Zqq vector:
      TLorentzVector Zqq_kinfit_cands;
      for( unsigned iCand=0; iCand<fitCands.size(); ++iCand )
        Zqq_kinfit_cands += *(fitCands[iCand]->getCurr4Vec());


      for( unsigned iCand=0; iCand<fitCands.size(); ++iCand ) {
        if( fitCands[iCand]->getCurr4Vec()->Pt()>0. && fitCands[iCand]->getIni4Vec()->Pt()>0. ) {
          float deltaE   = (fitCands[iCand]->getCurr4Vec()->Energy() - fitCands[iCand]->getIni4Vec()->Energy())/fitCands[iCand]->getIni4Vec()->Energy();
          float deltaEta = fitCands[iCand]->getCurr4Vec()->Eta()    - fitCands[iCand]->getIni4Vec()->Eta();
          float deltaPhi = fitCands[iCand]->getCurr4Vec()->DeltaPhi(fitCands[iCand]->getIni4Vec()->Phi());
          float deltaPt  = (fitCands[iCand]->getCurr4Vec()->Pt() - fitCands[iCand]->getIni4Vec()->Pt())/fitCands[iCand]->getIni4Vec()->Pt();
          if( fitCandTypes[iCand]==1 ) { //charged hadrons
            h1_deltaE_ch->Fill( deltaE ); 
            h1_deltaEta_ch->Fill( deltaEta ); 
            h1_deltaPhi_ch->Fill( deltaPhi ); 
            h1_deltaPt_ch->Fill( deltaPt ); 
          } else if( fitCandTypes[iCand]==4 ) { //photons
            h1_deltaE_gamma->Fill( deltaE ); 
            h1_deltaEta_gamma->Fill( deltaEta ); 
            h1_deltaPhi_gamma->Fill( deltaPhi ); 
            h1_deltaPt_gamma->Fill( deltaPt ); 
          } else if( fitCandTypes[iCand]==5 ) { //neutral hadrons
            h1_deltaE_nh->Fill( deltaE ); 
            h1_deltaEta_nh->Fill( deltaEta ); 
            h1_deltaPhi_nh->Fill( deltaPhi ); 
            h1_deltaPt_nh->Fill( deltaPt ); 
          }
        }
      } //for cands
    

      ZZ_kinfit_cands = Zll + Zqq_kinfit_cands;

      if( bothMatched )
        h1_ZZInvMass_MCassoc_kinfit_cands->Fill( ZZ_kinfit_cands.M(), eventWeight );
      

    } // if dataset

  } // if jet_presel
   

    // ------------------------
    //   KINEMATIC FIT: END
    // ------------------------

/*

    if( leptType==0 )
      h1_mZmumu->Fill( Zll.M(), eventWeight );
    else if( leptType==1 )
      h1_mZee->Fill( Zll.M(), eventWeight );
    else
      std::cout << "WARNING!! found incredible leptType: '" << leptType << "'." << std::endl;
    h1_massZll->Fill( Zll.M(), eventWeight );

    if( jetLead.Pt()>40. && jetLead2.Pt()>30. && fabs(jetLead.Eta())<2.5 && fabs(jetLead.Eta())<2.5 && lept1.Pt()>20. && lept2.Pt()>20. && Zll.M()>70. && Zll.M()<110. && Zjj_lead.M()>70. && Zjj_lead.M()<120. ) {
      h1_ZZInvMass_hiMass_loose->Fill( ZZ_lead.M(), eventWeight );
    }

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
    h1_ZZInvMass_hiMass_kinfit_jets->Fill( ZZ_kinfit_jets.M(), eventWeight );
    h1_ZZInvMass_hiMass_kinfit_cands->Fill( ZZ_kinfit_cands.M(), eventWeight );




    h1_massZjj->Fill( Zjj.M(), eventWeight );
    h1_massZjj_kinfit_jets->Fill( Zjj_kinfit_jets.M(), eventWeight );

    // MC association for signal only:
    bool isSignal=false;
    TString dataset_tstr(dataset_);
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

        h1_massZjj_MCassoc->Fill( Zjj.M(), eventWeight );
        h1_massZjj_MCassoc_kinfit->Fill( Zjj_kinfit_jets.M(), eventWeight );

        h1_ZZInvMass_hiMass_MCassoc->Fill( ZZ.M(), eventWeight );
        h1_ZZInvMass_hiMass_MCassoc_kinfit_jets->Fill( ZZ_kinfit_jets.M(), eventWeight );
        h1_ZZInvMass_hiMass_MCassoc_kinfit_cands->Fill( ZZ_kinfit_cands.M(), eventWeight );
        h1_ZZInvMass_hiMass_MCassoc_ZjjMassConstr->Fill( ZZ_constr.M(), eventWeight );
      }
    } //if signal


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
 // if( ZZ.M() > 190. && ZZ.M() < 200 ) {
 //   std::cout << "Run: " << run << " LS: " << LS << " event: " << event << std::endl;
 //   std::cout << " pt1: " << jet1.Pt() << " eta1: " << jet1.Eta() << " phi1: " << jet1.Phi() << " Rch1: " << Rch1 << std::endl;
 //   std::cout << " pt2: " << jet2.Pt() << " eta2: " << jet2.Eta() << " phi2: " << jet2.Phi() << " Rch2: " << Rch2 << std::endl;
 //   std::cout << " M(jj): " << Zjj.M() << std::endl;
 // }

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

    if( ptLept1 > 100. &&
        ptLept2 > 50. ) {

      if( Zll.M() > 86. &&
          Zll.M() < 96. ) {

        if( ptJetLead > 100. &&
            ptJetLead2 > 60. ) {
      
          if( Zjj_lead.M() > 80. &&
              Zjj_lead.M() < 105. ) {
        
              if( deltaRjj < 1.5 ) {

                  h1_ZZInvMass_hiMass_fullSelection_tight_lead->Fill( ZZ_lead.M(), eventWeight );

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


*/


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
  h1_rmsCandJet_all_presel->Write();
  h1_ptDJet_all_presel->Write();
  h1_nChargedJet_all_presel->Write();
  h1_nNeutralJet_all_presel->Write();
  h1_nJets_presel->Write();
  h1_nPairs_presel->Write();

  h1_ptJetBest1->Write();
  h1_ptDJetBest1->Write();
  h1_rmsCandJetBest1->Write();
  h1_nChargedJetBest1->Write();
  h1_nNeutralJetBest1->Write();

  h1_ptJetBest2->Write();
  h1_ptDJetBest2->Write();
  h1_rmsCandJetBest2->Write();
  h1_nChargedJetBest2->Write();
  h1_nNeutralJetBest2->Write();

  h1_ptJetRecoil->Write();
  h1_ptDJetRecoil->Write();
  h1_rmsCandJetRecoil->Write();
  h1_nChargedJetRecoil->Write();
  h1_nNeutralJetRecoil->Write();

  h1_deltaRll_presel->Write();
  h1_deltaRjj_presel->Write();

  h1_ptLept1_presel->Write();
  h1_ptLept2_presel->Write();

  h1_etaLept1_presel->Write();
  h1_etaLept2_presel->Write();

  h1_ptJet1_loose->Write();
  h1_ptJet1_tight->Write();

  h1_ptJet2_loose->Write();
  h1_ptJet2_tight->Write();

  h1_mZjj_presel->Write();

  h1_mZll_presel->Write();
  h1_mZmumu->Write();
  h1_mZmumu_presel->Write();
  h1_mZmumu_presel_0jets->Write();
  h1_mZee->Write();
  h1_mZee_presel->Write();
  h1_mZee_presel_0jets->Write();

  h1_mZll_loose->Write();
  h1_mZll_tight->Write();
  h1_mZll_opt400_LowEff->Write();
  h1_mZll_opt400_HighEff->Write();
  h1_mZll_opt500->Write();

  h1_mZqq_loose->Write();
  h1_mZqq_tight->Write();
  h1_mZqq_opt400_LowEff->Write();
  h1_mZqq_opt400_HighEff->Write();
  h1_mZqq_opt500->Write();

  h1_ZZInvMass_hiMass_loose->Write();
  h1_ZZInvMass_hiMass_loose_FINEBINNING->Write();
  h1_ZZInvMass_hiMass_tight->Write();
  h1_ZZInvMass_hiMass_tight_FINEBINNING->Write();

  h1_ZZInvMass_hiMass_opt200->Write();
  h1_ZZInvMass_hiMass_opt200_FINEBINNING->Write();
  h1_ZZInvMass_hiMass_opt300->Write();
  h1_ZZInvMass_hiMass_opt300_FINEBINNING->Write();
  h1_ZZInvMass_hiMass_opt400->Write();
  h1_ZZInvMass_hiMass_opt400_FINEBINNING->Write();
  h1_ZZInvMass_hiMass_opt500->Write();
  h1_ZZInvMass_hiMass_opt500_FINEBINNING->Write();
  h1_ZZInvMass_hiMass_opt600->Write();
  h1_ZZInvMass_hiMass_opt600_FINEBINNING->Write();

  h1_partFlavorJetOpt400_1->Write();
  h1_rmsCandJetOpt400_1->Write();
  h1_ptDJetOpt400_1->Write();
  h1_nChargedJetOpt400_1->Write();
  h1_nNeutralJetOpt400_1->Write();

  h1_partFlavorJetOpt400_2->Write();
  h1_rmsCandJetOpt400_2->Write();
  h1_ptDJetOpt400_2->Write();
  h1_nChargedJetOpt400_2->Write();
  h1_nNeutralJetOpt400_2->Write();

  h1_partFlavor_tight->Write();

  h1_ZZInvMass_MCassoc->Write();
  h1_ZZInvMass_MCassoc_ZjjMassConstr->Write();
  h1_ZZInvMass_MCassoc_kinfit_jets->Write();
  h1_ZZInvMass_MCassoc_kinfit_cands->Write();

  h2_mZjj_vs_mZZ->Write();
  h2_mZjj_vs_mZZ_kinfit->Write();

/*
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
  h1_mZmumu->Write();
  h1_mZee->Write();

  h1_massZjj->Write();
  h1_massZjj_kinfit_jets->Write();
  h1_massZjj_MCassoc->Write();
  h1_massZjj_MCassoc_kinfit->Write();
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
  h1_ZZInvMass_hiMass_kinfit_jets->Write();
  h1_ZZInvMass_hiMass_kinfit_cands->Write();
  h1_ZZInvMass_hiMass_MCassoc->Write();
  h1_ZZInvMass_hiMass_MCassoc_kinfit_jets->Write();
  h1_ZZInvMass_hiMass_MCassoc_kinfit_cands->Write();
  h1_ZZInvMass_hiMass_MCassoc_ZjjMassConstr->Write();


  h1_ZZInvMass_loMass_ZjjTag->Write();
  h1_ZZInvMass_medMass->Write();
  h1_ZZInvMass_hiMass_ZjjTag->Write();
  h1_ZZInvMass_hiMass_fullSelection_tightOLD->Write();
  h1_ZZInvMass_hiMass_fullSelection_tight->Write();
  h1_ZZInvMass_hiMass_fullSelection_tight_lead->Write();
  h1_ZZInvMass_hiMass_fullSelection_medium->Write();
  h1_ZZInvMass_hiMass_fullSelection_medium_ZjjMassConstr->Write();
  h1_ZZInvMass_hiMass_fullSelection_loose->Write();
  h1_ZZInvMass_hiMass_loose->Write();

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

  this->writeResponseHistos( outFile_, h1_response_vs_pt, "response" );
  this->writeResponseHistos( outFile_, h1_response_vs_pt_Rch050, "response_Rch050" );
  this->writeResponseHistos( outFile_, h1_response_vs_pt_Rch5070, "response_Rch5070" );
  this->writeResponseHistos( outFile_, h1_response_vs_pt_Rch70100, "response_Rch70100" );
*/

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

  outFile_->Close();

}



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
