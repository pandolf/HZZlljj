#include <TH2F.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>


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


void finalize(const std::string& dataset) {

  tree = new TChain("reducedTree");


  std::string infileName, treeName;


  if( dataset=="DATA_EG_37X" ) {

    addFile( "EG_Run2010A_Jul15thReReco_v1" );
    addFile( "EG_Run2010A_Jul26thReReco_v1" );

  } else if( dataset=="all" ) {

    finalize( "HZZ_qqll_gluonfusion_M130" );
    finalize( "HZZ_qqll_gluonfusion_M150" );
    finalize( "HZZ_qqll_gluonfusion_M200" );
    finalize( "HZZ_qqll_gluonfusion_M300" );
    finalize( "HZZ_qqll_gluonfusion_M400" );
    finalize( "HZZ_qqll_gluonfusion_M500" );
    finalize( "ZJets_madgraph" );
    return;

  } else {
  
    addFile( dataset );

  }



  std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);


  TH1D* h1_ptLept1 = new TH1D("ptLept1", "", 50, 10., 220.);
  h1_ptLept1->Sumw2();
  TH1D* h1_ptLept2 = new TH1D("ptLept2", "", 50, 10., 220.);
  h1_ptLept2->Sumw2();

  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 50, 20., 220.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 50, 20., 140.);
  h1_ptJet2->Sumw2();
  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 50, -5., 5.);
  h1_etaJet1->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 50, -5., 5.);
  h1_etaJet2->Sumw2();
  TH1D* h1_RchJet1 = new TH1D("RchJet1", "", 50, 0., 1.001);
  h1_RchJet1->Sumw2();
  TH1D* h1_RchJet2 = new TH1D("RchJet2", "", 50, 0., 1.001);
  h1_RchJet2->Sumw2();

  TH1D* h1_LeptLeptPt = new TH1D("LeptLeptPt", "", 50, 0., 200.);
  h1_LeptLeptPt->Sumw2();
  TH1D* h1_JetJetPt = new TH1D("JetJetPt", "", 50, 0., 200.);
  h1_JetJetPt->Sumw2();
  TH1D* h1_ptHardestZ = new TH1D("ptHardestZ", "", 50, 0., 240.);
  h1_ptHardestZ->Sumw2();

  TH1D* h1_deltaRjj = new TH1D("deltaRjj", "", 50, 0.5, 5.);
  h1_deltaRjj->Sumw2();
  TH1D* h1_deltaRZZ = new TH1D("deltaRZZ", "", 50, 0.5, 5.);
  h1_deltaRZZ->Sumw2();
  TH1D* h1_ptZZ = new TH1D("ptZZ", "", 50, 0., 100.);
  h1_ptZZ->Sumw2();

  int nBins_invMass = 50;
  float invMassMin = 0.;
  float invMassMax = 120.;

  TH1D* h1_LeptLeptInvMass = new TH1D("LeptLeptInvMass", "", nBins_invMass, invMassMin, invMassMax);
  h1_LeptLeptInvMass->Sumw2();
  TH1D* h1_MuMuInvMass = new TH1D("MuMuInvMass", "", nBins_invMass, invMassMin, invMassMax);
  h1_MuMuInvMass->Sumw2();
  TH1D* h1_EleEleInvMass = new TH1D("EleEleInvMass", "", nBins_invMass, invMassMin, invMassMax);
  h1_EleEleInvMass->Sumw2();
  TH1D* h1_JetJetInvMass = new TH1D("JetJetInvMass", "", nBins_invMass, 20., 200.);
  h1_JetJetInvMass->Sumw2();

  TH1D* h1_ZZInvMass_loMass = new TH1D("ZZInvMass_loMass", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass->Sumw2();
  TH1D* h1_ZZInvMass_medMass = new TH1D("ZZInvMass_medMass", "", nBins_invMass, 120., 280.);
  h1_ZZInvMass_medMass->Sumw2();
  TH1D* h1_ZZInvMass_hiMass = new TH1D("ZZInvMass_hiMass", "", nBins_invMass, 400., 600.);
  h1_ZZInvMass_hiMass->Sumw2();

  TH1D* h1_ZZInvMass_loMass_ZjjTag = new TH1D("ZZInvMass_loMass_ZjjTag", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass_ZjjTag->Sumw2();
  TH1D* h1_ZZInvMass_medMass_ZjjTag = new TH1D("ZZInvMass_medMass_ZjjTag", "", nBins_invMass, 120., 280.);
  h1_ZZInvMass_medMass_ZjjTag->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_ZjjTag = new TH1D("ZZInvMass_hiMass_ZjjTag", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_ZjjTag->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_tight = new TH1D("ZZInvMass_hiMass_fullSelection_tight", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_fullSelection_tight->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_medium = new TH1D("ZZInvMass_hiMass_fullSelection_medium", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_fullSelection_medium->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_fullSelection_loose = new TH1D("ZZInvMass_hiMass_fullSelection_loose", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_fullSelection_loose->Sumw2();

  TH1D* h1_ZZInvMass_medMass_ZjjTag_kinem = new TH1D("ZZInvMass_medMass_ZjjTag_kinem", "", nBins_invMass, 120., 280.);
  h1_ZZInvMass_medMass_ZjjTag_kinem->Sumw2();

  TH1D* h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag = new TH1D("ZZInvMass_loMass_ZjjTag_ZllAntiTag", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag = new TH1D("ZZInvMass_hiMass_ZjjTag_ZllAntiTag", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag->Sumw2();

  TH1D* h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40 = new TH1D("ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40", "", nBins_invMass, 90., 190.);
  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40->Sumw2();
  TH1D* h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40 = new TH1D("ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40", "", nBins_invMass, 250., 600.);
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40->Sumw2();


  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

  Float_t eMet;
  tree->SetBranchAddress("epfMet", &eMet);
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

  Float_t eJet1;
  tree->SetBranchAddress("eJet1", &eJet1);
  Float_t ptJet1;
  tree->SetBranchAddress("ptJet1", &ptJet1);
  Float_t etaJet1;
  tree->SetBranchAddress("etaJet1", &etaJet1);
  Float_t phiJet1;
  tree->SetBranchAddress("phiJet1", &phiJet1);
  Float_t eChargedHadronsJet1;
  tree->SetBranchAddress("eChargedHadronsJet1", &eChargedHadronsJet1);

  Float_t eJet2;
  tree->SetBranchAddress("eJet2", &eJet2);
  Float_t ptJet2;
  tree->SetBranchAddress("ptJet2", &ptJet2);
  Float_t etaJet2;
  tree->SetBranchAddress("etaJet2", &etaJet2);
  Float_t phiJet2;
  tree->SetBranchAddress("phiJet2", &phiJet2);
  Float_t eChargedHadronsJet2;
  tree->SetBranchAddress("eChargedHadronsJet2", &eChargedHadronsJet2);



  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    h1_ptLept1->Fill( ptLept1, eventWeight );
    h1_ptLept2->Fill( ptLept2, eventWeight );

    TLorentzVector lept1, lept2;
    lept1.SetPtEtaPhiE( ptLept1, etaLept1, phiLept1, eLept1 );
    lept2.SetPtEtaPhiE( ptLept2, etaLept2, phiLept2, eLept2 );

    TLorentzVector Zll = ( lept1 + lept2 );
    h1_LeptLeptPt->Fill( Zll.Pt(), eventWeight );
    h1_LeptLeptInvMass->Fill( Zll.M(), eventWeight );
    if( leptType==0 )
      h1_MuMuInvMass->Fill( Zll.M(), eventWeight );
    else if( leptType==1 )
      h1_EleEleInvMass->Fill( Zll.M(), eventWeight );
    else
      std::cout << "WARNING!! found incredible leptType: '" << leptType << "'." << std::endl;

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiE( ptJet1, etaJet1, phiJet1, eJet1 );
    jet2.SetPtEtaPhiE( ptJet2, etaJet2, phiJet2, eJet2 );

    h1_ptJet1->Fill( ptJet1, eventWeight );
    h1_ptJet2->Fill( ptJet2, eventWeight );

    h1_etaJet1->Fill( etaJet1, eventWeight );
    h1_etaJet2->Fill( etaJet2, eventWeight );

    float Rch1 = eChargedHadronsJet1/eJet1;
    float Rch2 = eChargedHadronsJet2/eJet2;
    h1_RchJet1->Fill( Rch1, eventWeight );
    h1_RchJet2->Fill( Rch2, eventWeight );

    float deltaRjj = jet1.DeltaR(jet2);
    h1_deltaRjj->Fill( deltaRjj, eventWeight );

    TLorentzVector Zjj = ( jet1 + jet2 );
    h1_JetJetInvMass->Fill( Zjj.M(), eventWeight );
    h1_JetJetPt->Fill( Zjj.Pt(), eventWeight );
 
    float ptHardestZ;
    if( Zjj.Pt() > Zll.Pt() )
      ptHardestZ = Zjj.Pt();
    else
      ptHardestZ = Zll.Pt();

    h1_ptHardestZ->Fill( ptHardestZ, eventWeight );

    float deltaRZZ = Zll.DeltaR(Zjj);
    h1_deltaRZZ->Fill( deltaRZZ, eventWeight );

    TLorentzVector ZZ = Zll + Zjj;

    float ptZZ = ZZ.Pt();
    h1_ptZZ->Fill( ptZZ, eventWeight );

    // ----------------------------------
    //   LOW MASS ANALYSIS (mH < 180)
    // ----------------------------------

    h1_ZZInvMass_loMass->Fill( ZZ.M(), eventWeight );

    if( Zjj.M() > 80. && ptJet1>35.) {

      h1_ZZInvMass_loMass_ZjjTag->Fill( ZZ.M(), eventWeight );

      if( Zll.M() < 80. ) {

        h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag->Fill( ZZ.M(), eventWeight );

        if( Rch1>0.4 && Rch2>0.4 ) {

          h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40->Fill( ZZ.M(), eventWeight );
 
        }
      }
    }


    // ----------------------------------
    //   HIGH MASS ANALYSIS (mH > 180)
    // ----------------------------------

    if( Zll.M()<80. || Zll.M()>100. ) continue;

    h1_ZZInvMass_medMass->Fill( ZZ.M(), eventWeight );
    h1_ZZInvMass_hiMass->Fill( ZZ.M(), eventWeight );

    if( Zjj.M() > 80. && Zjj.M()<110. && ptJet1>115. && ptJet2>55. && deltaRjj<1.5 && ptZZ>30.) {

      h1_ZZInvMass_hiMass_ZjjTag->Fill( ZZ.M(), eventWeight );

    }

    if( ptJet1 > 100. &&
        ptJet2 > 50. &&
        ptLept1 > 60. &&
        ptLept2 > 60. &&
        Zjj.M() > 80. &&
        Zjj.M() < 110. &&
        Zll.M() > 86. &&
        Zll.M() < 96. &&
        ptHardestZ > 120. &&
        deltaRjj < 1.5 &&
        ptZZ > 20. &&
        deltaRZZ > 2.8 ) {

      h1_ZZInvMass_hiMass_fullSelection_tight->Fill( ZZ.M(), eventWeight );

    }

    if( ptJet1 > 80. &&
        ptJet2 > 40. &&
        ptLept1 > 60. &&
        ptLept2 > 60. &&
        Zjj.M() > 80. &&
        Zjj.M() < 110. &&
        Zll.M() > 86. &&
        Zll.M() < 96. &&
        ptHardestZ > 100. &&
        deltaRjj < 2. &&
        ptZZ > 20. &&
        deltaRZZ > 2.8 ) {

      h1_ZZInvMass_hiMass_fullSelection_medium->Fill( ZZ.M(), eventWeight );

    }

    if( ptJet1 > 80. &&
        ptJet2 > 30. &&
        ptLept1 > 50. &&
        ptLept2 > 50. &&
        Zjj.M() > 80. &&
        Zjj.M() < 110. &&
        Zll.M() > 86. &&
        Zll.M() < 96. &&
        ptHardestZ > 80. &&
        deltaRjj < 2. &&
        ptZZ > 20. &&
        deltaRZZ > 2.8 ) {

      h1_ZZInvMass_hiMass_fullSelection_loose->Fill( ZZ.M(), eventWeight );

    }

    if( Zjj.M() > 80. && Zjj.M()<110. && ptJet1>45. && ptJet2>30. && deltaRjj>2.2 ) {

      h1_ZZInvMass_medMass_ZjjTag->Fill( ZZ.M(), eventWeight );

      if( ptZZ>15. && deltaRZZ < 3.5 ) {

        h1_ZZInvMass_medMass_ZjjTag_kinem->Fill( ZZ.M(), eventWeight );
   
      }

    }

  } //for entries



  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaHZZlljj_"+dataset;
  else {
   if(dataset!="") outfileName = "HZZlljj_"+dataset;
   else outfileName = "HZZlljj";
  }

  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();

  h1_totalLumi->Write();

  h1_ptLept1->Write();
  h1_ptLept2->Write();

  h1_ptJet1->Write();
  h1_ptJet2->Write();
  h1_etaJet1->Write();
  h1_etaJet2->Write();
  h1_RchJet1->Write();
  h1_RchJet2->Write();

  h1_deltaRjj->Write();
  h1_ptZZ->Write();
  h1_deltaRZZ->Write();

  h1_LeptLeptPt->Write();
  h1_JetJetPt->Write();
  h1_ptHardestZ->Write();

  h1_LeptLeptInvMass->Write();
  h1_MuMuInvMass->Write();
  h1_EleEleInvMass->Write();
  h1_JetJetInvMass->Write();

  h1_ZZInvMass_loMass->Write();
  h1_ZZInvMass_medMass->Write();
  h1_ZZInvMass_hiMass->Write();

  h1_ZZInvMass_loMass_ZjjTag->Write();
  h1_ZZInvMass_medMass_ZjjTag->Write();
  h1_ZZInvMass_hiMass_ZjjTag->Write();
  h1_ZZInvMass_hiMass_fullSelection_tight->Write();
  h1_ZZInvMass_hiMass_fullSelection_medium->Write();
  h1_ZZInvMass_hiMass_fullSelection_loose->Write();

  h1_ZZInvMass_medMass_ZjjTag_kinem->Write();

  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag->Write();
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag->Write();

  h1_ZZInvMass_loMass_ZjjTag_ZllAntiTag_Rch40->Write();
  h1_ZZInvMass_hiMass_ZjjTag_ZllAntiTag_Rch40->Write();

  outFile->Close();

  delete h1_ptLept1;
  h1_ptLept1 = 0;
  delete h1_ptLept2;
  h1_ptLept2 = 0;

  delete h1_ptJet1;
  h1_ptJet1 = 0;
  delete h1_ptJet2;
  h1_ptJet2 = 0;
  delete h1_etaJet1;
  h1_etaJet1 = 0;
  delete h1_etaJet2;
  h1_etaJet2 = 0;
  delete h1_RchJet1;
  h1_RchJet1 = 0;
  delete h1_RchJet2;
  h1_RchJet2 = 0;

  delete h1_deltaRjj;
  h1_deltaRjj = 0;
  delete h1_ptZZ;
  h1_ptZZ = 0;
  delete h1_deltaRZZ;
  h1_deltaRZZ = 0;

  delete h1_totalLumi;
  h1_totalLumi = 0;

  delete h1_JetJetPt;
  h1_JetJetPt = 0;
  delete h1_LeptLeptPt;
  h1_LeptLeptPt = 0;

  delete h1_LeptLeptInvMass;
  h1_LeptLeptInvMass = 0;
  delete h1_MuMuInvMass;
  h1_MuMuInvMass = 0;
  delete h1_EleEleInvMass;
  h1_EleEleInvMass = 0;
  delete h1_JetJetInvMass;
  h1_JetJetInvMass = 0;
  delete h1_ptHardestZ;
  h1_ptHardestZ = 0;

  delete h1_ZZInvMass_loMass;
  h1_ZZInvMass_loMass = 0;
  delete h1_ZZInvMass_medMass;
  h1_ZZInvMass_medMass = 0;
  delete h1_ZZInvMass_hiMass;
  h1_ZZInvMass_hiMass = 0;

  delete h1_ZZInvMass_loMass_ZjjTag;
  h1_ZZInvMass_loMass_ZjjTag = 0;
  delete h1_ZZInvMass_medMass_ZjjTag;
  h1_ZZInvMass_medMass_ZjjTag = 0;
  delete h1_ZZInvMass_hiMass_ZjjTag;
  h1_ZZInvMass_hiMass_ZjjTag = 0;
  delete h1_ZZInvMass_hiMass_fullSelection_tight;
  h1_ZZInvMass_hiMass_fullSelection_tight = 0;
  delete h1_ZZInvMass_hiMass_fullSelection_medium;
  h1_ZZInvMass_hiMass_fullSelection_medium = 0;
  delete h1_ZZInvMass_hiMass_fullSelection_loose;
  h1_ZZInvMass_hiMass_fullSelection_loose = 0;

  delete h1_ZZInvMass_medMass_ZjjTag_kinem;
  h1_ZZInvMass_medMass_ZjjTag_kinem = 0;

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


  
