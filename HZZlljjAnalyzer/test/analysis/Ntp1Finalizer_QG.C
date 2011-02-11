#include "Ntp1Finalizer_QG.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"


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

//Double_t ptBins_ZG[11];
//fitTools::getBins( 11, ptBins_ZG, 100., 500. );


  const int nBins = 20;
  Double_t ptBins[nBins+1];
  fitTools::getBins_int( nBins+1, ptBins, 15., 1000. );


  std::vector<TH1D*> vh1_nCharged_gluon;
  std::vector<TH1D*> vh1_nNeutral_gluon;
  std::vector<TH1D*> vh1_ptD_gluon;
  std::vector<TH1D*> vh1_rmsCand_gluon;

  std::vector<TH1D*> vh1_nCharged_quark;
  std::vector<TH1D*> vh1_nNeutral_quark;
  std::vector<TH1D*> vh1_ptD_quark;
  std::vector<TH1D*> vh1_rmsCand_quark;

  std::vector<TH2D*> vh2_ptD_vs_nCharged_gluon;
  std::vector<TH2D*> vh2_ptD_vs_rmsCand_gluon;
  std::vector<TH2D*> vh2_rmsCand_vs_nCharged_gluon;
  std::vector<TH2D*> vh2_nCharged_vs_nNeutral_gluon;

  std::vector<TH2D*> vh2_ptD_vs_nCharged_quark;
  std::vector<TH2D*> vh2_ptD_vs_rmsCand_quark;
  std::vector<TH2D*> vh2_rmsCand_vs_nCharged_quark;
  std::vector<TH2D*> vh2_nCharged_vs_nNeutral_quark;

//TH1D* h1_ptJet_gluon = new TH1D("ptJet_gluon", "", 30, 30., 300.);
//h1_ptJet_gluon->Sumw2();
//TH1D* h1_ptJet_quark = new TH1D("ptJet_quark", "", 30, 30., 300.);
//h1_ptJet_quark->Sumw2();

  for( unsigned iBin=0; iBin<nBins; iBin++ ) {

  //float ptMin = 100.+20.*iBin;
  //float ptMax = ptMin + 20.;

    float ptMin = ptBins[iBin];
    float ptMax = ptBins[iBin+1];

    char histoname[200];
    
    sprintf( histoname, "nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_gluon_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_gluon_new->Sumw2();
    vh1_nCharged_gluon.push_back(h1_nCharged_gluon_new);
    sprintf( histoname, "nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_quark_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_quark_new->Sumw2();
    vh1_nCharged_quark.push_back(h1_nCharged_quark_new);

    sprintf( histoname, "nNeutral_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_gluon_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_gluon_new->Sumw2();
    vh1_nNeutral_gluon.push_back(h1_nNeutral_gluon_new);
    sprintf( histoname, "nNeutral_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_quark_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_quark_new->Sumw2();
    vh1_nNeutral_quark.push_back(h1_nNeutral_quark_new);

    sprintf( histoname, "rmsCand_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_gluon_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_gluon_new->Sumw2();
    vh1_rmsCand_gluon.push_back(h1_rmsCand_gluon_new);
    sprintf( histoname, "rmsCand_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_quark_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_quark_new->Sumw2();
    vh1_rmsCand_quark.push_back(h1_rmsCand_quark_new);

    sprintf( histoname, "ptD_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_gluon_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_ptD_gluon_new->Sumw2();
    vh1_ptD_gluon.push_back(h1_ptD_gluon_new);
    sprintf( histoname, "ptD_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_quark_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_ptD_quark_new->Sumw2();
    vh1_ptD_quark.push_back(h1_ptD_quark_new);

    sprintf( histoname, "ptD_vs_rmsCand_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_vs_rmsCand_gluon_new = new TH2D(histoname, "", 50, 0., 0.1, 50, 0., 1.);
    h2_ptD_vs_rmsCand_gluon_new->Sumw2();
    vh2_ptD_vs_rmsCand_gluon.push_back(h2_ptD_vs_rmsCand_gluon_new);
    sprintf( histoname, "ptD_vs_rmsCand_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_vs_rmsCand_quark_new = new TH2D(histoname, "", 50, 0., 0.1, 50, 0., 1.);
    h2_ptD_vs_rmsCand_quark_new->Sumw2();
    vh2_ptD_vs_rmsCand_quark.push_back(h2_ptD_vs_rmsCand_quark_new);

    sprintf( histoname, "ptD_vs_nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_vs_nCharged_gluon_new = new TH2D(histoname, "", 51, -0.5, 50.5, 50, 0., 1.);
    h2_ptD_vs_nCharged_gluon_new->Sumw2();
    vh2_ptD_vs_nCharged_gluon.push_back(h2_ptD_vs_nCharged_gluon_new);
    sprintf( histoname, "ptD_vs_nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_vs_nCharged_quark_new = new TH2D(histoname, "", 51, -0.5, 50.5, 50, 0., 1.);
    h2_ptD_vs_nCharged_quark_new->Sumw2();
    vh2_ptD_vs_nCharged_quark.push_back(h2_ptD_vs_nCharged_quark_new);

    sprintf( histoname, "nCharged_vs_nNeutral_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_nCharged_vs_nNeutral_gluon_new = new TH2D(histoname, "", 51, -0.5, 50.5, 51, -0.5, 50.5);
    h2_nCharged_vs_nNeutral_gluon_new->Sumw2();
    vh2_nCharged_vs_nNeutral_gluon.push_back(h2_nCharged_vs_nNeutral_gluon_new);
    sprintf( histoname, "nCharged_vs_nNeutral_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_nCharged_vs_nNeutral_quark_new = new TH2D(histoname, "", 51, -0.5, 50.5, 51, -0.5, 50.5);
    h2_nCharged_vs_nNeutral_quark_new->Sumw2();
    vh2_nCharged_vs_nNeutral_quark.push_back(h2_nCharged_vs_nNeutral_quark_new);

    sprintf( histoname, "rmsCand_vs_nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_rmsCand_vs_nCharged_gluon_new = new TH2D(histoname, "", 51, -0.5, 50.5, 50, 0., 0.1);
    h2_rmsCand_vs_nCharged_gluon_new->Sumw2();
    vh2_rmsCand_vs_nCharged_gluon.push_back(h2_rmsCand_vs_nCharged_gluon_new);
    sprintf( histoname, "rmsCand_vs_nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_rmsCand_vs_nCharged_quark_new = new TH2D(histoname, "", 51, -0.5, 50.5, 50, 0., 0.1);
    h2_rmsCand_vs_nCharged_quark_new->Sumw2();
    vh2_rmsCand_vs_nCharged_quark.push_back(h2_rmsCand_vs_nCharged_quark_new);

  } //for bins


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


    //for( unsigned iJet=0; iJet<nJet; ++iJet ) {
    for( unsigned iJet=0; (iJet<nJet && iJet<3); ++iJet ) { //only 3 leading jets considered

      TLorentzVector thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      if( fabs(thisJet.Eta())>2. ) continue;

      float deltaRmin=999.;
      int partFlavor=-1;
      TLorentzVector foundPart;

      for( unsigned iPart=0; iPart<nPart; iPart++ ) {

        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart]);

        float thisDeltaR = thisJet.DeltaR(thisPart);

        if( thisDeltaR < deltaRmin ) {
          deltaRmin = thisDeltaR;
          partFlavor = pdgIdPart[iPart];
          foundPart = thisPart;
        }

      } //for partons

      if( deltaRmin > 0.5 ) continue;
      //if( deltaRmin > 0.5 ) partFlavor=21; //lets try this


      int thisBin=-1;
      if( thisJet.Pt() > ptBins[nBins] ) {
        thisBin = nBins-1;
      } else {
        for( unsigned int iBin=0; iBin<nBins; ++iBin ) {
          if( thisJet.Pt()>ptBins[iBin] && thisJet.Pt()<ptBins[iBin+1] ) {
            thisBin = iBin;
            break;
          }
        }
      }

/*
if( iEntry>500000 ) std::cout << "c" << std::endl;
      if( thisBin<0 ) continue; //shouldnt be possible because of the preselection pt cut
      if( thisBin>nBins ) std::cout << "!!!!!!!!!!!!!        thisBin: " << thisBin << std::endl;

if( iEntry>500000 ) {
  std::cout << "iJet: " << iJet << " this Bin: " << thisBin << std::endl;
  std::cout << vh1_nCharged_quark.size() << std::endl;
  std::cout << vh1_nNeutral_quark.size() << std::endl;
  std::cout << vh1_ptD_quark.size() << std::endl;
  std::cout << vh1_rmsCand_quark.size() << std::endl;

  std::cout << vh2_ptD_vs_nCharged_quark.size() << std::endl;
  std::cout << vh2_ptD_vs_rmsCand_quark.size() << std::endl;
  std::cout << vh2_rmsCand_vs_nCharged_quark.size() << std::endl;
  std::cout << vh2_nCharged_vs_nNeutral_quark.size() << std::endl;
  std::cout << "ptD: " << ptDJet[iJet] << std::endl;
  std::cout << "rmsCand: " << rmsCandJet[iJet] << std::endl;
  std::cout << "nCharged: " << nChargedJet[iJet] << std::endl;
  std::cout << "nNeutral: " << nNeutralJet[iJet] << std::endl;
  std::cout << "eventWeight: " << eventWeight << std::endl;
  std::cout << "ptD: " << ptDJet[iJet] << std::endl;
}
*/
      if( abs(partFlavor)< 7 ) { //quark
        //h1_ptJet_quark[thisBin]->Fill( ptJet[iJet], eventWeight );
        vh1_nCharged_quark[thisBin]->Fill( nChargedJet[iJet], eventWeight );
        vh1_nNeutral_quark[thisBin]->Fill( nNeutralJet[iJet], eventWeight );
        vh1_ptD_quark[thisBin]->Fill( ptDJet[iJet], eventWeight );
        vh1_rmsCand_quark[thisBin]->Fill( rmsCandJet[iJet], eventWeight );

        vh2_ptD_vs_nCharged_quark[thisBin]->Fill( nChargedJet[iJet], ptDJet[iJet], eventWeight );
        vh2_ptD_vs_rmsCand_quark[thisBin]->Fill( rmsCandJet[iJet], ptDJet[iJet], eventWeight );
        vh2_rmsCand_vs_nCharged_quark[thisBin]->Fill( nChargedJet[iJet], rmsCandJet[iJet], eventWeight );
        vh2_nCharged_vs_nNeutral_quark[thisBin]->Fill( nNeutralJet[iJet], nChargedJet[iJet], eventWeight );

      } else if( partFlavor==21 ) { //gluon
        //h1_ptJet_gluon[thisBin]->Fill( ptJet[iJet], eventWeight );
        vh1_nCharged_gluon[thisBin]->Fill( nChargedJet[iJet], eventWeight );
        vh1_nNeutral_gluon[thisBin]->Fill( nNeutralJet[iJet], eventWeight );
        vh1_ptD_gluon[thisBin]->Fill( ptDJet[iJet], eventWeight );
        vh1_rmsCand_gluon[thisBin]->Fill( rmsCandJet[iJet], eventWeight );

        vh2_ptD_vs_nCharged_gluon[thisBin]->Fill( nChargedJet[iJet], ptDJet[iJet], eventWeight );
        vh2_ptD_vs_rmsCand_gluon[thisBin]->Fill( rmsCandJet[iJet], ptDJet[iJet], eventWeight );
        vh2_rmsCand_vs_nCharged_gluon[thisBin]->Fill( nChargedJet[iJet], rmsCandJet[iJet], eventWeight );
        vh2_nCharged_vs_nNeutral_gluon[thisBin]->Fill( nNeutralJet[iJet], nChargedJet[iJet], eventWeight );

      }

    } // for jets

  } //for entries




  outFile_->cd();


  for( unsigned iBin=0; iBin<nBins; ++iBin ) {
    //h1_ptJet_gluon[iBin]->Write();
    vh1_nCharged_gluon[iBin]->Write();
    vh1_nNeutral_gluon[iBin]->Write();
    vh1_ptD_gluon[iBin]->Write();
    vh1_rmsCand_gluon[iBin]->Write();

    //h1_ptJet_quark[iBin]->Write();
    vh1_nCharged_quark[iBin]->Write();
    vh1_nNeutral_quark[iBin]->Write();
    vh1_ptD_quark[iBin]->Write();
    vh1_rmsCand_quark[iBin]->Write();

    vh2_ptD_vs_nCharged_gluon[iBin]->Write();
    vh2_ptD_vs_rmsCand_gluon[iBin]->Write();
    vh2_rmsCand_vs_nCharged_gluon[iBin]->Write();
    vh2_nCharged_vs_nNeutral_gluon[iBin]->Write();

    vh2_ptD_vs_nCharged_quark[iBin]->Write();
    vh2_ptD_vs_rmsCand_quark[iBin]->Write();
    vh2_rmsCand_vs_nCharged_quark[iBin]->Write();
    vh2_nCharged_vs_nNeutral_quark[iBin]->Write();

  }


  outFile_->Close();

}

