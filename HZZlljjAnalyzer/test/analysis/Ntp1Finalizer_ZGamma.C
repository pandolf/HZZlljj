#include "Ntp1Finalizer_ZGamma.h"

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





// constructor:

Ntp1Finalizer_ZGamma::Ntp1Finalizer_ZGamma( const std::string& dataset ) : Ntp1Finalizer( "ZGamma", dataset ) {

}





void Ntp1Finalizer_ZGamma::finalize() {

  if( outFile_==0 ) this->createOutputFile();


  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  if( dataset_=="Run2010B_runs146240_146733" )
    h1_totalLumi->SetBinContent(1, 1220000.);
  else
    h1_totalLumi->SetBinContent(1, totalLumi_);


  TH1D* h1_nJets = new TH1D("nJets", "", 10, -0.5, 9.5);
  h1_nJets->Sumw2();

  Double_t ptBins[101];
  fitTools::getBins( 101, ptBins, 20., 1000. );

  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 100, ptBins);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 100, ptBins);
  h1_ptJet2->Sumw2();
  TH1D* h1_ptJet3 = new TH1D("ptJet3", "", 100, ptBins);
  h1_ptJet3->Sumw2();
  TH1D* h1_ptJet4 = new TH1D("ptJet4", "", 100, ptBins);
  h1_ptJet4->Sumw2();

  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 100, -6., 6.);
  h1_etaJet1->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 100, -6., 6.);
  h1_etaJet2->Sumw2();
  TH1D* h1_etaJet3 = new TH1D("etaJet3", "", 100, -6., 6.);
  h1_etaJet3->Sumw2();
  TH1D* h1_etaJet4 = new TH1D("etaJet4", "", 100, -6., 6.);
  h1_etaJet4->Sumw2();

  TH1D* h1_deltaR12 = new TH1D("deltaR12", "", 100, 0., 7.);
  h1_deltaR12->Sumw2();



  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);

  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Float_t ptZGamma;
  tree_->SetBranchAddress("ptZGamma", &ptZGamma);
  Float_t etaZGamma;
  tree_->SetBranchAddress("etaZGamma", &etaZGamma);

  Int_t nJet;
  tree_->SetBranchAddress("nJet", &nJet);
  Float_t eJet[10];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[10];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t etaJet[10];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[10];
  tree_->SetBranchAddress("phiJet", phiJet);



  int nEntries = tree_->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);

//std::cout << "n: " << nJet <<  std::endl;

//for( unsigned i=0; i<nJet; ++i ) std::cout << "pt: " << ptJet[i] << std::endl;

    if( eventWeight <= 0. ) eventWeight = 1.;


    if( fabs(etaZGamma)>2.5 || ptZGamma<50. ) continue;

    if( nJet<2 ) continue;

    if( ptJet[1]<20. ) continue;

    h1_nJets->Fill( nJet, eventWeight );

    if( nJet>0 ) {
      if( ptJet[0]>20. ) {
        h1_ptJet1->Fill( ptJet[0], eventWeight );
        h1_etaJet1->Fill( etaJet[0], eventWeight );
      }
    }

    if( nJet>1 ) {
      if( ptJet[1]>20. ) {
        h1_ptJet2->Fill( ptJet[1], eventWeight );
        h1_etaJet2->Fill( etaJet[1], eventWeight );
        TLorentzVector jet1, jet2;
        jet1.SetPtEtaPhiE( ptJet[0], etaJet[0], phiJet[0], eJet[0]); 
        jet2.SetPtEtaPhiE( ptJet[1], etaJet[1], phiJet[1], eJet[1]); 
        h1_deltaR12->Fill( jet1.DeltaR(jet2), eventWeight);
      }
    }

    if( nJet>2 ) {
      if( ptJet[2]>20. ) {
        h1_ptJet3->Fill( ptJet[2], eventWeight );
        h1_etaJet3->Fill( etaJet[2], eventWeight );
      }
    }

    if( nJet>3 ) {
      if( ptJet[3]>20. ) {
        h1_ptJet4->Fill( ptJet[3], eventWeight );
        h1_etaJet4->Fill( etaJet[3], eventWeight );
      }
    }


  } //for entries



  outFile_->cd();

  h1_totalLumi->Write();

  h1_nJets->Write();

  h1_ptJet1->Write();
  h1_ptJet2->Write();
  h1_ptJet3->Write();
  h1_ptJet4->Write();

  h1_etaJet1->Write();
  h1_etaJet2->Write();
  h1_etaJet3->Write();
  h1_etaJet4->Write();

  h1_deltaR12->Write();

  outFile_->Close();

}


