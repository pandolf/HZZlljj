#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <cstring>
#include <string>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>

void OptimizingHelicity_HWWlvjj( ){

  TCanvas* myc1 = new TCanvas("myc1", "myc1", 600, 600);
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptFit(111110);
  gStyle->SetOptFile(1);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);

  TH1F* Sensitivity = new TH1F("Sensitivity","",100,0.,1.);

  double Ev_WW=0., Ev_WJets=0., Ev_tt=0.;
  double Sens=0., Best_Sens=0., Best_cut=0.;

  TFile* file_WW = TFile::Open("HWWlvjj_WW500_helicity_ALL.root");
  TTree* tree_WW = (TTree*)file_WW->Get("Tree_optim_hely");
  Float_t eventWeight_WW; 
  tree_WW->SetBranchAddress("eventWeight",&eventWeight_WW);
  tree_WW->GetEntry(1);

  TFile* file_tt = TFile::Open("HWWlvjj_TT_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_3_helicity_ALL.root");
  TTree* tree_tt = (TTree*)file_tt->Get("Tree_optim_hely");
  Float_t eventWeight_tt;
  tree_tt->SetBranchAddress("eventWeight",&eventWeight_tt);
  tree_tt->GetEntry(1);

  TFile* file_WJets = TFile::Open("HWWlvjj_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  TTree* tree_WJets = (TTree*)file_WJets->Get("Tree_optim_hely");
  Float_t eventWeight_WJets;
  tree_WJets->SetBranchAddress("eventWeight",&eventWeight_WJets);
  tree_WJets->GetEntry(1);
  std::cout<<"Pesi: WW, tt, WJets: "<<eventWeight_WW<<" "<<eventWeight_tt<<" "<<eventWeight_WJets<<std::endl;
  char cut[400];

	for(int i=1; i<100; i++ ){
          
          float thresh = (float)i/100.;
          sprintf(cut,"(helicityLD_kinfit>%f)", thresh);

     	  Ev_WW = (tree_WW->GetEntries( cut ))*eventWeight_WW;
          Ev_tt = (tree_tt->GetEntries( cut ))*eventWeight_tt;
          Ev_WJets = (tree_WJets->GetEntries( cut ))*eventWeight_WJets;

          Sens = ( Ev_WW>0. && (Ev_tt>0. && Ev_WJets>0.) ) ? (Ev_WW*1000.)/sqrt((Ev_WJets*1000.)+(Ev_tt*1000.)) : 0.;

          if(Sens != 0.) std::cout<<Ev_WW*1000.<<"  "<<Ev_WJets*1000.<<"  "<<Ev_tt*1000.<<" Sens "<<Sens<<std::endl;

          Sensitivity->SetBinContent(i,Sens);

          if( Sens > Best_Sens ){
             Best_cut = (float)i/100.;
             Best_Sens = Sens;
          }
      }

  std::cout<<"The best cut in helicity is:  "<<Best_cut<<"  Which give us a S/sqrt(B): "<<Best_Sens<<std::endl;
  myc1->cd();
  Sensitivity->Draw();
  }
