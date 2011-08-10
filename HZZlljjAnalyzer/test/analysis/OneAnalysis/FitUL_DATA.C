#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <cstring>
#include <string>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TChain.h>
#include <TMath.h>
#include <TLegend.h>
void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog);


TCanvas* FitUL_DATA(){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptFit(111110);
  gStyle->SetOptFile(1);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);
  gStyle->SetFillColor(0);

  TCanvas* myc1 = new TCanvas("myc1", "myc1", 600, 600);
 
  TH1F* mWW_peak_DATA = new TH1F("mWW_peak_DATA","Before Correction",56, 140., 700.);
  mWW_peak_DATA->Sumw2();
  mWW_peak_DATA->SetXTitle("WW Invariant Mass [GeV]");
  TH1F* mWW_side_DATA = new TH1F("mWW_side_DATA","Before Correction",56, 140., 700.);
  mWW_side_DATA->Sumw2();
  mWW_side_DATA->SetXTitle("WW Invariant Mass [GeV]");
  TLegend *leg = new TLegend( 0.5,0.7,0.8,0.8 );
  leg->AddEntry(mWW_side_DATA,"Side-bands region","F");
  leg->AddEntry(mWW_peak_DATA,"Peak region","F");

  TH1F* mJJ_tot_DATA = new TH1F("mJJ_tot_DATA","",100, 0., 250.);
  mJJ_tot_DATA->Sumw2();


  // HISTO with Xaxis log
  double Bin[30]={1.};
  double* BinP;
  BinP=Bin;
  getBins_int( 10, BinP, 140., 700., true);
  TH1F* mWW_peak_log_DATA = new TH1F("mWW_peak_log_DATA","", 9, BinP);
  mWW_peak_log_DATA->Sumw2();
  TH1F* mWW_side_log_DATA = new TH1F("mWW_side_log_DATA","", 9, BinP);
  mWW_side_log_DATA->Sumw2();

  TChain ch("Tree_FITUL");
  ch.Add("HWWlvjj_DATA_6july_helicity_ALL.root");
  ch.Draw("mWW");

  float mJJ, mWW, eventWeight;
  ch.SetBranchAddress("mJJ", &mJJ);
  ch.SetBranchAddress("mWW", &mWW);
  ch.SetBranchAddress("eventWeight", &eventWeight);

  for( int iEntry=0; iEntry<ch.GetEntries() ; iEntry++ ){
  ch.GetEntry(iEntry);
   if( mJJ>60. && mJJ<100 )                          {  mWW_peak_DATA->Fill(mWW,eventWeight); mWW_peak_log_DATA->Fill(mWW,eventWeight); }
   if( (mJJ>40 && mJJ<60 ) || (mJJ>100 && mJJ<160) ) {  mWW_side_DATA->Fill(mWW,eventWeight); mWW_side_log_DATA->Fill(mWW,eventWeight); } 
   mJJ_tot_DATA->Fill(mJJ,eventWeight);
  }

  // M(JJ)
  myc1->cd();
  mJJ_tot_DATA->Draw();
  myc1->SaveAs("mJJ_DATA.eps");

  // M(WW) COMPARISON
  mWW_side_DATA->SetLineColor(2);
  mWW_side_DATA->Scale( 1/mWW_side_DATA->Integral() );
  mWW_peak_DATA->Scale( 1/mWW_peak_DATA->Integral() );
  mWW_peak_DATA->Draw("HISTO");
  mWW_side_DATA->Draw("sameHISTO");
  leg->Draw("same");
  myc1->SaveAs("mWW_DATA.eps");

  // M(WW)_LOG peak & side
  mWW_peak_log_DATA->Draw();
  myc1->SaveAs("mWW_peak_log_DATA.eps");
  mWW_side_log_DATA->Draw();
  myc1->SaveAs("mWW_side_log_DATA.eps");

  // Now I correct my distribution
  TH1F* mWW_peak_log2_DATA = new TH1F("mWW_peak_log2_DATA","After Correction", 56, 140., 700);
  mWW_peak_log2_DATA->Sumw2();
  mWW_peak_log2_DATA->SetXTitle("WW Invariant Mass [GeV]");
  TH1F* mWW_side_log2_DATA = new TH1F("mWW_side_log2_DATA","After Correction", 56, 140., 700) ;
  mWW_side_log2_DATA->Sumw2();
  mWW_side_log2_DATA->SetXTitle("WW Invariant Mass [GeV]");
  TLegend *leg_CORR = new TLegend( 0.5,0.7,0.8,0.8 );
  leg_CORR->AddEntry(mWW_side_log2_DATA,"Side-bands region","F");
  leg_CORR->AddEntry(mWW_peak_log2_DATA,"Peak region","F");

  for( int iEntry=0; iEntry<ch.GetEntries() ; iEntry++ ){
   ch.GetEntry(iEntry);
   if( mJJ>60. && mJJ<100 )                          {  mWW_peak_log2_DATA->Fill( mWW,eventWeight ); }
   if( (mJJ>40 && mJJ<60 ) || (mJJ>100 && mJJ<160) ) {  mWW_side_log2_DATA->Fill( mWW,eventWeight*( 1.26421 + (-1.79074*pow(10,-3))*mWW +
                                                                           (-1.51523*pow(10,-6))*pow(mWW,2) + (4.98491*pow(10,-9))*pow(mWW,3) ) ); }
  }
  // M(WW) COMPARISON_RATIO
  mWW_side_log2_DATA->SetLineColor(2);
  mWW_side_log2_DATA->Scale( 1/mWW_side_log2_DATA->Integral() );
  mWW_peak_log2_DATA->Scale( 1/mWW_peak_log2_DATA->Integral() );
  mWW_side_log2_DATA->Draw("HISTO");
  mWW_peak_log2_DATA->Draw("HISTOsame");
  leg_CORR->Draw("same");
  myc1->SaveAs("mWW_CORR_DATA.eps");
  return myc1;

}


void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog) {

  Double_t Lower_exact;
  int nBins = nBins_total-1;
  const double dx = (plotLog) ? pow((xmax / xmin), (1. / (double)nBins)) : ((xmax - xmin) / (double)nBins);
  Lower[0] = xmin;
  Lower_exact = Lower[0];
  for (int i = 1; i != nBins; ++i) {

    if (plotLog) {
      Lower_exact *= dx;
      Lower[i] = TMath::Ceil(Lower_exact);
    } else {
      Lower[i] = TMath::Ceil(Lower[i-1] + dx);
    }

  }

  Lower[nBins] = xmax;

}
