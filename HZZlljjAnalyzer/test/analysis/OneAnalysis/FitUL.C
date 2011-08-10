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


TCanvas* FitUL(){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptFit(111110);
  gStyle->SetOptFile(1);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);
  gStyle->SetFillColor(0);

  TCanvas* myc1 = new TCanvas("myc1", "myc1", 600, 600);
 
  TH1F* mWW_peak = new TH1F("mWW_peak","Before Correction",56, 140., 700.);
  mWW_peak->Sumw2();
  mWW_peak->SetXTitle("WW Invariant Mass [GeV]");
  TH1F* mWW_side = new TH1F("mWW_side","Before Correction",56, 140., 700.);
  mWW_side->Sumw2();
  mWW_side->SetXTitle("WW Invariant Mass [GeV]");
  TLegend *leg = new TLegend( 0.5,0.7,0.8,0.8 );
  leg->AddEntry(mWW_side,"Side-bands region","F");
  leg->AddEntry(mWW_peak,"Peak region","F");

  TH1F* mJJ_tot = new TH1F("mWW_side","",100, 0., 250.);
  mJJ_tot->Sumw2();


  // HISTO with Xaxis log
  double Bin[30]={1.};
  double* BinP;
  BinP=Bin;
  getBins_int( 10, BinP, 140., 700., true);
  TH1F* mWW_peak_log = new TH1F("mWW_peak_log","", 9, BinP);
  mWW_peak_log->Sumw2();
  TH1F* mWW_side_log = new TH1F("mWW_side_log","", 9, BinP);
  mWW_side_log->Sumw2();

  TChain ch("Tree_FITUL");
  ch.Add("HWWlvjj_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  ch.Add("HWWlvjj_DY_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  ch.Add("HWWlvjj_TToBLNu_TuneZ2_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  ch.Add("HWWlvjj_TT_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_3_helicity_ALL.root");
  ch.Add("HWWlvjj_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2_helicity_ALL.root");
  ch.Add("HWWlvjj_QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6_helicity_ALL.root");
  ch.Add("HWWlvjj_QCD_EMEnriched_TuneZ2_7TeV-pythia6_3_helicity_ALL.root");
  ch.Add("HWWlvjj_QCD_BCtoE_TuneZ2_7TeV-pythia6_helicity_ALL.root");
  ch.Add("HWWlvjj_GJet_TuneZ2_7TeV-alpgen_helicity_ALL.root");
  ch.Draw("mWW");

  float mJJ, mWW, eventWeight;
  ch.SetBranchAddress("mJJ", &mJJ);
  ch.SetBranchAddress("mWW", &mWW);
  ch.SetBranchAddress("eventWeight", &eventWeight);

  for( int iEntry=0; iEntry<ch.GetEntries() ; iEntry++ ){
  ch.GetEntry(iEntry);
  //std::cout<<iEntry<<"/"<<ch.GetEntries()<<":  "<< mJJ <<", "<< mWW <<std::endl;
   if( mJJ>60. && mJJ<100 )                          {  mWW_peak->Fill(mWW,eventWeight); mWW_peak_log->Fill(mWW,eventWeight); }
   if( (mJJ>40 && mJJ<60 ) || (mJJ>100 && mJJ<160) ) {  mWW_side->Fill(mWW,eventWeight); mWW_side_log->Fill(mWW,eventWeight); } 
   mJJ_tot->Fill(mJJ,eventWeight);
   //if(mWW>600 && mWW<650 && ( mJJ>60. && mJJ<100 ) )                          { std::cout<<"IN "<<mJJ<<"  "<<eventWeight<<std::endl; }
   //if(mWW>600 && mWW<650 && ( (mJJ>40 && mJJ<60 ) || (mJJ>100 && mJJ<160) ) ) { std::cout<<"OUT"<<mJJ<<"  "<<eventWeight<<std::endl; }
  }

  // M(JJ)
  myc1->cd();
  mJJ_tot->Draw();
  myc1->SaveAs("mJJ.eps");
  // M(WW) COMPARISON
  mWW_side->SetLineColor(2);
  mWW_side->Scale( 1/mWW_side->Integral() );
  mWW_peak->Scale( 1/mWW_peak->Integral() );
  mWW_peak->Draw("HISTO");
  mWW_side->Draw("sameHISTO");
  leg->Draw("same");
  myc1->SaveAs("mWW.eps");

  // M(WW) RATIO
  TH1F* peak_side = new TH1F( *mWW_peak );
  peak_side->Divide( mWW_side );
  peak_side->Draw();
  myc1->SaveAs("peak_side.eps");
  // M(WW)_LOG peak & side
  mWW_peak_log->Draw();
  myc1->SaveAs("mWW_peak_log.eps");
  mWW_side_log->Draw();
  myc1->SaveAs("mWW_side_log.eps");
  // M(WW) RATIO xLOG

  myc1->cd();
  TH1F* peak_side_log = new TH1F( *mWW_peak_log );
  peak_side_log->Sumw2();
  peak_side_log->SetXTitle("WW Invariant Mass [GeV]");
  peak_side_log->SetYTitle("Ratio peak/side");
  peak_side_log->Divide( mWW_side_log );
  peak_side_log->Draw();
  myc1->SaveAs("peak_side_log.eps");

  // NOW FIT
  TF1 *myfit = new TF1("myfit","[0]+[1]*x+[2]*x^2+[3]*x^3", 140, 700);
  myfit->SetParName(0,"c0");
  myfit->SetParName(1,"c1");
  myfit->SetParName(2,"c2");
  myfit->SetParName(3,"c3");
  peak_side_log->Fit("myfit");
  peak_side_log->Draw();
  myc1->SaveAs("Ratio_fitted.eps");
    
  // Now I correct my distribution
  TH1F* mWW_peak_log2 = new TH1F("mWW_peak_log2","After Correction", 56, 140., 700);
  mWW_peak_log2->Sumw2();
  mWW_peak_log2->SetXTitle("WW Invariant Mass [GeV]");
  TH1F* mWW_side_log2 = new TH1F("mWW_side_log2","After Correction", 56, 140., 700) ;
  mWW_side_log2->Sumw2();
  mWW_side_log2->SetXTitle("WW Invariant Mass [GeV]");
  TLegend *leg_CORR = new TLegend( 0.5,0.7,0.8,0.8 );
  leg_CORR->AddEntry(mWW_side_log2,"Side-bands region","F");
  leg_CORR->AddEntry(mWW_peak_log2,"Peak region","F");

  for( int iEntry=0; iEntry<ch.GetEntries() ; iEntry++ ){
   ch.GetEntry(iEntry);
   if( mJJ>60. && mJJ<100 )                          {  mWW_peak_log2->Fill( mWW,eventWeight ); }
   if( (mJJ>40 && mJJ<60 ) || (mJJ>100 && mJJ<160) ) {  mWW_side_log2->Fill( mWW,eventWeight*( 1.26421 + (-1.79074*pow(10,-3))*mWW +
                                                                           (-1.51523*pow(10,-6))*pow(mWW,2) + (4.98491*pow(10,-9))*pow(mWW,3) ) ); }
  }
  // M(WW) COMPARISON_RATIO
  mWW_side_log2->SetLineColor(2);
  mWW_side_log2->Scale( 1/mWW_side_log2->Integral() );
  mWW_peak_log2->Scale( 1/mWW_peak_log2->Integral() );
  mWW_side_log2->Draw("HISTO");
  mWW_peak_log2->Draw("HISTOsame");
  leg_CORR->Draw("same");
  myc1->SaveAs("mWW_CORR.eps");
  return myc1;  

  //RATIO CORRECTED
  TH1F* peak_side_log2 = new TH1F( *mWW_peak_log2 );
  peak_side_log2->Sumw2();
  peak_side_log2->Divide( mWW_side_log2 );
  peak_side_log2->Draw();
  myc1->SaveAs("peak_side_log2.eps");


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
