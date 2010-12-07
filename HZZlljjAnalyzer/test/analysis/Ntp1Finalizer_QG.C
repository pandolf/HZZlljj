#include "Ntp1Finalizer_QG.h"

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

Ntp1Finalizer_QG::Ntp1Finalizer_QG( const std::string& dataset ) : Ntp1Finalizer( "QG", dataset ) {

}





void Ntp1Finalizer_QG::finalize() {

  if( outFile_==0 ) this->createOutputFile();


  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  if( dataset_=="Run2010B_runs146240_146733" )
    h1_totalLumi->SetBinContent(1, 1220000.);
  else
    h1_totalLumi->SetBinContent(1, totalLumi_);

  Double_t ptBins_ZG[11];
  fitTools::getBins( 11, ptBins_ZG, 100., 500. );

  Double_t ptBins[11];
  fitTools::getBins( 11, ptBins, 30., 500. );


  TH1D* h1_ptJet_gluon = new TH1D("ptJet_gluon", "", 30, 30., 300.);
  h1_ptJet_gluon->Sumw2();
  TH1D* h1_ptJet_quark = new TH1D("ptJet_quark", "", 30, 30., 300.);
  h1_ptJet_quark->Sumw2();

  TH1D* h1_nCharged_gluon = new TH1D("nCharged_gluon", "", 51, -0.5, 50.5);
  h1_nCharged_gluon->Sumw2();
  TH1D* h1_nCharged_quark = new TH1D("nCharged_quark", "", 51, -0.5, 50.5);
  h1_nCharged_quark->Sumw2();

  TH1D* h1_nNeutral_gluon = new TH1D("nNeutral_gluon", "", 51, -0.5, 50.5);
  h1_nNeutral_gluon->Sumw2();
  TH1D* h1_nNeutral_quark = new TH1D("nNeutral_quark", "", 51, -0.5, 50.5);
  h1_nNeutral_quark->Sumw2();

  TH1D* h1_rmsCand_gluon = new TH1D("rmsCand_gluon", "", 50, 0., 0.1);
  h1_rmsCand_gluon->Sumw2();
  TH1D* h1_rmsCand_quark = new TH1D("rmsCand_quark", "", 50, 0., 0.1);
  h1_rmsCand_quark->Sumw2();

  TH1D* h1_ptD_gluon = new TH1D("ptD_gluon", "", 50, 0., 1.);
  h1_ptD_gluon->Sumw2();
  TH1D* h1_ptD_quark = new TH1D("ptD_quark", "", 50, 0., 1.);
  h1_ptD_quark->Sumw2();



  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);

  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Int_t nJet;
  tree_->SetBranchAddress("nJet", &nJet);
  Float_t eJet[20];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[20];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t etaJet[20];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[20];
  tree_->SetBranchAddress("phiJet", phiJet);
  Int_t nChargedJet[20];
  tree_->SetBranchAddress("nChargedJet", nChargedJet);
  Int_t nNeutralJet[20];
  tree_->SetBranchAddress("nNeutralJet", nNeutralJet);
  Float_t rmsCandJet[20];
  tree_->SetBranchAddress("rmsCandJet", rmsCandJet);
  Float_t ptDJet[20];
  tree_->SetBranchAddress("ptDJet", ptDJet);

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

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;


    for( unsigned iJet=0; iJet<nJet; ++iJet ) {

      TLorentzVector thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      if( fabs(thisJet.Eta())>2. ) continue;

      float deltaRmin=999.;
      int partFlavor=-1;

      for( unsigned iPart=0; iPart<nPart; iPart++ ) {

        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart]);

        float thisDeltaR = thisJet.DeltaR(thisPart);

        if( thisDeltaR < deltaRmin ) {
          deltaRmin = thisDeltaR;
          partFlavor = pdgIdPart[iPart];
        }

      } //for partons

      //if( deltaRmin > 0.5 ) continue;
      if( deltaRmin > 0.5 ) partFlavor=21; //lets try this

      if( abs(partFlavor)< 7 ) { //quark
        h1_ptJet_quark->Fill( ptJet[iJet], eventWeight );
        h1_nCharged_quark->Fill( nChargedJet[iJet], eventWeight );
        h1_nNeutral_quark->Fill( nNeutralJet[iJet], eventWeight );
        h1_ptD_quark->Fill( ptDJet[iJet], eventWeight );
        h1_rmsCand_quark->Fill( rmsCandJet[iJet], eventWeight );
      } else if( partFlavor==21 ) { //gluon
        h1_ptJet_gluon->Fill( ptJet[iJet], eventWeight );
        h1_nCharged_gluon->Fill( nChargedJet[iJet], eventWeight );
        h1_nNeutral_gluon->Fill( nNeutralJet[iJet], eventWeight );
        h1_ptD_gluon->Fill( ptDJet[iJet], eventWeight );
        h1_rmsCand_gluon->Fill( rmsCandJet[iJet], eventWeight );
      }

    } // for jets

  } //for entries



  outFile_->cd();


  h1_ptJet_gluon->Write();
  h1_nCharged_gluon->Write();
  h1_nNeutral_gluon->Write();
  h1_ptD_gluon->Write();
  h1_rmsCand_gluon->Write();

  h1_ptJet_quark->Write();
  h1_nCharged_quark->Write();
  h1_nNeutral_quark->Write();
  h1_ptD_quark->Write();
  h1_rmsCand_quark->Write();


  outFile_->Close();

}


