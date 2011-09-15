#include "Ntp1Finalizer_ComputeQGLikelihood.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "CommonTools/fitTools.h"





// constructor:

Ntp1Finalizer_ComputeQGLikelihood::Ntp1Finalizer_ComputeQGLikelihood( const std::string& dataset ) : Ntp1Finalizer( "ComputeQGLikelihood", dataset ) {

}


void Ntp1Finalizer_ComputeQGLikelihood::addFile(const std::string& dataset) {

  std::string infileName = "QG_2ndLevelTreeW_" + dataset + ".root"; //the W is important: means that files have passed treatment (merging and weights)
  std::string treeName = infileName +"/reducedTree";
  tree_->Add(treeName.c_str());
  std::cout << "-> Added " << treeName << ". Tree has " << tree_->GetEntries() << " entries." << std::endl;
  TFile* infile = TFile::Open(infileName.c_str(), "READ");
  TH1F* h1_lumi = (TH1F*)infile->Get("lumi");
  infile->Close();


}




void Ntp1Finalizer_ComputeQGLikelihood::finalize() {

  if( outFile_==0 ) this->createOutputFile();



  const int nBins = 20;
  Double_t ptBins[nBins+1];
  fitTools::getBins_int( nBins+1, ptBins, 15., 1000. );


  std::vector<TH1D*> vh1_QGLikelihood_gluon;
  std::vector<TH1D*> vh1_QGLikelihood_quark;

  std::vector<TH1D*> vh1_QGLikelihood_norms_gluon;
  std::vector<TH1D*> vh1_QGLikelihood_norms_quark;

  std::vector<TH1D*> vh1_QGLikelihood_noptD_gluon;
  std::vector<TH1D*> vh1_QGLikelihood_noptD_quark;

  std::vector<TH1D*> vh1_QGLikelihood_norms_noptD_gluon;
  std::vector<TH1D*> vh1_QGLikelihood_norms_noptD_quark;

  std::vector<TH1D*> vh1_QGLikelihood_onlyNch_gluon;
  std::vector<TH1D*> vh1_QGLikelihood_onlyNch_quark;

  std::vector<TH1D*> vh1_rmsCand_cutOnQGLikelihood_norms_gluon;
  std::vector<TH1D*> vh1_rmsCand_cutOnQGLikelihood_norms_quark;


  for( unsigned iBin=0; iBin<nBins; iBin++ ) {

    float ptMin = ptBins[iBin];
    float ptMax = ptBins[iBin+1];

    char histoname[200];
    
    sprintf( histoname, "QGLikelihood_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_gluon_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_gluon_new->Sumw2();
    vh1_QGLikelihood_gluon.push_back(h1_QGLikelihood_gluon_new);

    sprintf( histoname, "QGLikelihood_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_quark_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_quark_new->Sumw2();
    vh1_QGLikelihood_quark.push_back(h1_QGLikelihood_quark_new);


    sprintf( histoname, "QGLikelihood_norms_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_norms_gluon_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_norms_gluon_new->Sumw2();
    vh1_QGLikelihood_norms_gluon.push_back(h1_QGLikelihood_norms_gluon_new);

    sprintf( histoname, "QGLikelihood_norms_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_norms_quark_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_norms_quark_new->Sumw2();
    vh1_QGLikelihood_norms_quark.push_back(h1_QGLikelihood_norms_quark_new);


    sprintf( histoname, "QGLikelihood_noptD_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_noptD_gluon_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_noptD_gluon_new->Sumw2();
    vh1_QGLikelihood_noptD_gluon.push_back(h1_QGLikelihood_noptD_gluon_new);

    sprintf( histoname, "QGLikelihood_noptD_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_noptD_quark_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_noptD_quark_new->Sumw2();
    vh1_QGLikelihood_noptD_quark.push_back(h1_QGLikelihood_noptD_quark_new);


    sprintf( histoname, "QGLikelihood_norms_noptD_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_norms_noptD_gluon_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_norms_noptD_gluon_new->Sumw2();
    vh1_QGLikelihood_norms_noptD_gluon.push_back(h1_QGLikelihood_norms_noptD_gluon_new);

    sprintf( histoname, "QGLikelihood_norms_noptD_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_norms_noptD_quark_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_norms_noptD_quark_new->Sumw2();
    vh1_QGLikelihood_norms_noptD_quark.push_back(h1_QGLikelihood_norms_noptD_quark_new);


    sprintf( histoname, "QGLikelihood_onlyNch_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_onlyNch_gluon_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_onlyNch_gluon_new->Sumw2();
    vh1_QGLikelihood_onlyNch_gluon.push_back(h1_QGLikelihood_onlyNch_gluon_new);

    sprintf( histoname, "QGLikelihood_onlyNch_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_QGLikelihood_onlyNch_quark_new = new TH1D(histoname, "", 50, 0., 1.);
    h1_QGLikelihood_onlyNch_quark_new->Sumw2();
    vh1_QGLikelihood_onlyNch_quark.push_back(h1_QGLikelihood_onlyNch_quark_new);


    sprintf( histoname, "rmsCand_cutOnQGLikelihood_norms_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_cutOnQGLikelihood_norms_gluon_new = new TH1D(histoname, "", 100, 0., 0.1 );
    h1_rmsCand_cutOnQGLikelihood_norms_gluon_new->Sumw2();
    vh1_rmsCand_cutOnQGLikelihood_norms_gluon.push_back(h1_rmsCand_cutOnQGLikelihood_norms_gluon_new);

    sprintf( histoname, "rmsCand_cutOnQGLikelihood_norms_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_cutOnQGLikelihood_norms_quark_new = new TH1D(histoname, "", 100, 0., 0.1 );
    h1_rmsCand_cutOnQGLikelihood_norms_quark_new->Sumw2();
    vh1_rmsCand_cutOnQGLikelihood_norms_quark.push_back(h1_rmsCand_cutOnQGLikelihood_norms_quark_new);

  } //for bins


  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);

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


  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root", nBins);
  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1.root", nBins);


  int nEntries = tree_->GetEntries();
  //nEntries = 10000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {
  //for(int iEntry=0; iEntry<500000; ++iEntry) {

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

      if( thisBin==-1 ) continue;

    //float QGLikelihood = qglikeli->computeQGLikelihood( thisJet.Pt(), nChargedJet[iJet], nNeutralJet[iJet], ptDJet[iJet], rmsCandJet[iJet] );
    //float QGLikelihood_norms = qglikeli->computeQGLikelihood( thisJet.Pt(), nChargedJet[iJet], nNeutralJet[iJet], ptDJet[iJet], -1. );
    //float QGLikelihood_noptD = qglikeli->computeQGLikelihood( thisJet.Pt(), nChargedJet[iJet], nNeutralJet[iJet], -1., rmsCandJet[iJet] );
    //float QGLikelihood_norms_noptD = qglikeli->computeQGLikelihood( thisJet.Pt(), nChargedJet[iJet], nNeutralJet[iJet], -1., -1. );
    //float QGLikelihood_onlyNch = qglikeli->computeQGLikelihood( thisJet.Pt(), nChargedJet[iJet], -1., -1., -1. );

      float QGLikelihood = qglikeli->computeQGLikelihoodPU( thisJet.Pt(), rhoPF, nChargedJet[iJet], nNeutralJet[iJet], ptDJet[iJet], rmsCandJet[iJet] );
      float QGLikelihood_norms = qglikeli->computeQGLikelihoodPU( thisJet.Pt(), rhoPF, nChargedJet[iJet], nNeutralJet[iJet], ptDJet[iJet], -1. );
      float QGLikelihood_noptD = qglikeli->computeQGLikelihoodPU( thisJet.Pt(), rhoPF, nChargedJet[iJet], nNeutralJet[iJet], -1., rmsCandJet[iJet] );
      float QGLikelihood_norms_noptD = qglikeli->computeQGLikelihoodPU( thisJet.Pt(), rhoPF, nChargedJet[iJet], nNeutralJet[iJet], -1., -1. );
      float QGLikelihood_onlyNch = qglikeli->computeQGLikelihoodPU( thisJet.Pt(), rhoPF, nChargedJet[iJet], -1., -1., -1. );

      if( abs(partFlavor)< 7 ) { //quark
        
        vh1_QGLikelihood_quark[thisBin]->Fill( QGLikelihood, eventWeight );
        vh1_QGLikelihood_norms_quark[thisBin]->Fill( QGLikelihood_norms, eventWeight );
        vh1_QGLikelihood_noptD_quark[thisBin]->Fill( QGLikelihood_noptD, eventWeight );
        vh1_QGLikelihood_norms_noptD_quark[thisBin]->Fill( QGLikelihood_norms_noptD, eventWeight );
        if( QGLikelihood_norms<0.9 ) vh1_rmsCand_cutOnQGLikelihood_norms_quark[thisBin]->Fill( rmsCandJet[iJet], eventWeight );
        vh1_QGLikelihood_onlyNch_quark[thisBin]->Fill( QGLikelihood_onlyNch, eventWeight );

      } else if( partFlavor==21 ) { //gluon

        vh1_QGLikelihood_gluon[thisBin]->Fill( QGLikelihood, eventWeight );
        vh1_QGLikelihood_norms_gluon[thisBin]->Fill( QGLikelihood_norms, eventWeight );
        vh1_QGLikelihood_noptD_gluon[thisBin]->Fill( QGLikelihood_noptD, eventWeight );
        vh1_QGLikelihood_norms_noptD_gluon[thisBin]->Fill( QGLikelihood_norms_noptD, eventWeight );
        if( QGLikelihood_norms<0.9 ) vh1_rmsCand_cutOnQGLikelihood_norms_gluon[thisBin]->Fill( rmsCandJet[iJet], eventWeight );
        vh1_QGLikelihood_onlyNch_gluon[thisBin]->Fill( QGLikelihood_onlyNch, eventWeight );

      }

    } // for jets


  } //for entries


  std::vector<TGraphErrors*> vgr_eff_vs_rej;
  std::vector<TGraphErrors*> vgr_eff_vs_rej_norms;
  std::vector<TGraphErrors*> vgr_eff_vs_rej_noptD;
  std::vector<TGraphErrors*> vgr_eff_vs_rej_norms_noptD;
  std::vector<TGraphErrors*> vgr_eff_vs_rej_onlyNch;

  for( unsigned iBin=0; iBin<nBins; ++iBin ) {

    float ptMin = ptBins[iBin];
    float ptMax = ptBins[iBin+1];

    char graphName[200];
    sprintf( graphName, "eff_vs_rej_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* newGraph = new TGraphErrors(0);
    newGraph->SetName(graphName);
    vgr_eff_vs_rej.push_back(newGraph);

    sprintf( graphName, "eff_vs_rej_norms_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* newGraph_norms = new TGraphErrors(0);
    newGraph_norms->SetName(graphName);
    vgr_eff_vs_rej_norms.push_back(newGraph_norms);

    sprintf( graphName, "eff_vs_rej_noptD_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* newGraph_noptD = new TGraphErrors(0);
    newGraph_noptD->SetName(graphName);
    vgr_eff_vs_rej_noptD.push_back(newGraph_noptD);

    sprintf( graphName, "eff_vs_rej_norms_noptD_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* newGraph_norms_noptD = new TGraphErrors(0);
    newGraph_norms_noptD->SetName(graphName);
    vgr_eff_vs_rej_norms_noptD.push_back(newGraph_norms_noptD);

    sprintf( graphName, "eff_vs_rej_onlyNch_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* newGraph_onlyNch = new TGraphErrors(0);
    newGraph_onlyNch->SetName(graphName);
    vgr_eff_vs_rej_onlyNch.push_back(newGraph_onlyNch);

  }


  // compute eff-rej curves:
  for( unsigned iBin=0; iBin<nBins; ++iBin ) {

    for( unsigned iHistoBin=1; iHistoBin<vh1_QGLikelihood_quark[iBin]->GetNbinsX(); ++iHistoBin ) {

      float eff_q = vh1_QGLikelihood_quark[iBin]->Integral(iHistoBin, vh1_QGLikelihood_quark[iBin]->GetNbinsX())/vh1_QGLikelihood_quark[iBin]->Integral(1, vh1_QGLikelihood_quark[iBin]->GetNbinsX());
      float eff_g = vh1_QGLikelihood_gluon[iBin]->Integral(iHistoBin, vh1_QGLikelihood_gluon[iBin]->GetNbinsX())/vh1_QGLikelihood_gluon[iBin]->Integral(1, vh1_QGLikelihood_gluon[iBin]->GetNbinsX());

      float eff_q_norms = vh1_QGLikelihood_norms_quark[iBin]->Integral(iHistoBin, vh1_QGLikelihood_norms_quark[iBin]->GetNbinsX())/vh1_QGLikelihood_norms_quark[iBin]->Integral(1, vh1_QGLikelihood_norms_quark[iBin]->GetNbinsX());
      float eff_g_norms = vh1_QGLikelihood_norms_gluon[iBin]->Integral(iHistoBin, vh1_QGLikelihood_norms_gluon[iBin]->GetNbinsX())/vh1_QGLikelihood_norms_gluon[iBin]->Integral(1, vh1_QGLikelihood_norms_gluon[iBin]->GetNbinsX());

      float eff_q_noptD = vh1_QGLikelihood_noptD_quark[iBin]->Integral(iHistoBin, vh1_QGLikelihood_noptD_quark[iBin]->GetNbinsX())/vh1_QGLikelihood_noptD_quark[iBin]->Integral(1, vh1_QGLikelihood_noptD_quark[iBin]->GetNbinsX());
      float eff_g_noptD = vh1_QGLikelihood_noptD_gluon[iBin]->Integral(iHistoBin, vh1_QGLikelihood_noptD_gluon[iBin]->GetNbinsX())/vh1_QGLikelihood_noptD_gluon[iBin]->Integral(1, vh1_QGLikelihood_noptD_gluon[iBin]->GetNbinsX());

      float eff_q_norms_noptD = vh1_QGLikelihood_norms_noptD_quark[iBin]->Integral(iHistoBin, vh1_QGLikelihood_norms_noptD_quark[iBin]->GetNbinsX())/vh1_QGLikelihood_norms_noptD_quark[iBin]->Integral(1, vh1_QGLikelihood_norms_noptD_quark[iBin]->GetNbinsX());
      float eff_g_norms_noptD = vh1_QGLikelihood_norms_noptD_gluon[iBin]->Integral(iHistoBin, vh1_QGLikelihood_norms_noptD_gluon[iBin]->GetNbinsX())/vh1_QGLikelihood_norms_noptD_gluon[iBin]->Integral(1, vh1_QGLikelihood_norms_noptD_gluon[iBin]->GetNbinsX());

      float eff_q_onlyNch = vh1_QGLikelihood_onlyNch_quark[iBin]->Integral(iHistoBin, vh1_QGLikelihood_onlyNch_quark[iBin]->GetNbinsX())/vh1_QGLikelihood_onlyNch_quark[iBin]->Integral(1, vh1_QGLikelihood_onlyNch_quark[iBin]->GetNbinsX());
      float eff_g_onlyNch = vh1_QGLikelihood_onlyNch_gluon[iBin]->Integral(iHistoBin, vh1_QGLikelihood_onlyNch_gluon[iBin]->GetNbinsX())/vh1_QGLikelihood_onlyNch_gluon[iBin]->Integral(1, vh1_QGLikelihood_onlyNch_gluon[iBin]->GetNbinsX());

      vgr_eff_vs_rej[iBin]->SetPoint( vgr_eff_vs_rej[iBin]->GetN(), eff_q, 1.-eff_g ); 
      vgr_eff_vs_rej_norms[iBin]->SetPoint( vgr_eff_vs_rej_norms[iBin]->GetN(), eff_q_norms, 1.-eff_g_norms ); 
      vgr_eff_vs_rej_noptD[iBin]->SetPoint( vgr_eff_vs_rej_noptD[iBin]->GetN(), eff_q_noptD, 1.-eff_g_noptD ); 
      vgr_eff_vs_rej_norms_noptD[iBin]->SetPoint( vgr_eff_vs_rej_norms_noptD[iBin]->GetN(), eff_q_norms_noptD, 1.-eff_g_norms_noptD ); 
      vgr_eff_vs_rej_onlyNch[iBin]->SetPoint( vgr_eff_vs_rej_onlyNch[iBin]->GetN(), eff_q_onlyNch, 1.-eff_g_onlyNch ); 

    }

  }


  outFile_->cd();


  for( unsigned iBin=0; iBin<nBins; ++iBin ) {

    vh1_QGLikelihood_gluon[iBin]->Write();
    vh1_QGLikelihood_quark[iBin]->Write();

    vh1_QGLikelihood_norms_gluon[iBin]->Write();
    vh1_QGLikelihood_norms_quark[iBin]->Write();

    vh1_QGLikelihood_noptD_gluon[iBin]->Write();
    vh1_QGLikelihood_noptD_quark[iBin]->Write();

    vh1_QGLikelihood_norms_noptD_gluon[iBin]->Write();
    vh1_QGLikelihood_norms_noptD_quark[iBin]->Write();

    vh1_QGLikelihood_onlyNch_gluon[iBin]->Write();
    vh1_QGLikelihood_onlyNch_quark[iBin]->Write();

    vh1_rmsCand_cutOnQGLikelihood_norms_gluon[iBin]->Write();
    vh1_rmsCand_cutOnQGLikelihood_norms_quark[iBin]->Write();

    vgr_eff_vs_rej[iBin]->Write();
    vgr_eff_vs_rej_norms[iBin]->Write();
    vgr_eff_vs_rej_noptD[iBin]->Write();
    vgr_eff_vs_rej_norms_noptD[iBin]->Write();
    vgr_eff_vs_rej_onlyNch[iBin]->Write();

  }

  outFile_->Close();

}

