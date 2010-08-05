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


void addInput(const std::string& dataset);


void finalize(const std::string& dataset) {

  tree = new TChain("reducedTree");


  std::string infileName, treeName;


  if( dataset=="DATA_EG_37X" ) {

    addInput( "EG_Run2010A_Jul15thReReco_v1" );
    addInput( "EG_Run2010A_Jul26thReReco_v1" );

  } else {
  
    addInput( dataset );

  }



  std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);


  int nBins_invMass = 50;

  TH1D* h1_LeptLeptInvMass = new TH1D("LeptLeptInvMass", "", nBins_invMass, 40., 140.);
  h1_LeptLeptInvMass->Sumw2();
  TH1D* h1_JetJetInvMass = new TH1D("JetJetInvMass", "", nBins_invMass, 40., 140.);
  h1_JetJetInvMass->Sumw2();
  TH1D* h1_ZZInvMass = new TH1D("ZZInvMass", "", nBins_invMass, 80., 230.);
  h1_ZZInvMass->Sumw2();



  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Int_t event;
  tree->SetBranchAddress("event", &event);

  Float_t eMet;
  tree->SetBranchAddress("epfMet", &eMet);
  Float_t phiMet;
  tree->SetBranchAddress("phipfMet", &phiMet);


  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

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


  Float_t eventWeight = 1.;


  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    TLorentzVector lept1, lept2;
    lept1.SetPtEtaPhiE( ptLept1, etaLept1, phiLept1, eLept1 );
    lept2.SetPtEtaPhiE( ptLept2, etaLept2, phiLept2, eLept2 );

    TLorentzVector Zll = ( lept1 + lept2 );
    h1_LeptLeptInvMass->Fill( Zll.M(), eventWeight );

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiE( ptJet1, etaJet1, phiJet1, eJet1 );
    jet2.SetPtEtaPhiE( ptJet2, etaJet2, phiJet2, eJet2 );

    TLorentzVector Zjj = ( jet1 + jet2 );
    h1_JetJetInvMass->Fill( Zjj.M(), eventWeight );

    TLorentzVector ZZ = Zll + Zjj;
    h1_ZZInvMass->Fill( ZZ.M(), eventWeight );


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


  h1_LeptLeptInvMass->Write();
  h1_JetJetInvMass->Write();
  h1_ZZInvMass->Write();


  outFile->Close();

  delete h1_LeptLeptInvMass;
  h1_LeptLeptInvMass = 0;
  delete h1_JetJetInvMass;
  h1_JetJetInvMass = 0;
  delete h1_ZZInvMass;
  h1_ZZInvMass = 0;

  delete tree;
  tree = 0;

  totalLumi = 0.;

}


void addInput( const std::string& dataset ) {

  std::string infileName = "files_HZZlljj_2ndLevel_" + dataset+"_" +".txt";
  TH1F* h1_lumi;


  //open from file.txt:
  FILE* iff = fopen(infileName.c_str(),"r");
  if(iff == 0) {
    std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    infileName = "HZZlljj_2ndLevelTree_" + dataset + ".root";
    std::string treeName = infileName +"/reducedTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    TFile* infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
      std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    infile->Close();

  } else {

    char singleLine[500];

    while( fscanf(iff, "%s", singleLine) !=EOF ) {

      std::string rootfilename(singleLine);
      std::string treename = rootfilename + "/reducedTree";
      std::cout << "-> Added " << treename;
      tree->Add(treename.c_str());
      TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
      h1_lumi = (TH1F*)infile->Get("lumi");
      if( h1_lumi!=0 ) {
        totalLumi += h1_lumi->GetBinContent(1);
        std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
      } else {
        std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
      }
      infile->Close();

    }
    fclose(iff);

  }

} //addinput


  
