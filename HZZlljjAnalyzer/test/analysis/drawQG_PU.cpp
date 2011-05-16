#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"



void drawCompareQCD( DrawBase* db, const std::string& ptRange, TFile* fileSignal_NEW, TFile* fileSignal_OLD );
void drawEffvsRej( DrawBase* db, TFile* qcdFile_OLD, TFile* qcdFile_NEW, const std::string& varName, const std::string& ptRange );


int main(int argc, char* argv[]) {

  if(  argc != 1 ) {
    std::cout << "USAGE: ./drawQG_PU" << std::endl;
    exit(23);
  }



  std::string lept_type = "ALL";

  DrawBase* db = new DrawBase("QG");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "QGPUPlots";
  db->set_outputdir(outputdir_str);



  TFile* signalFile_NEW = TFile::Open("HZZlljj_SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_2_opt400_0btag_optLD_ALL_looseBTag_NEW.root");
  db->add_mcFile( signalFile_NEW, "signalNEW", "Flat10 PU PDFs", 38, 3004);

  TFile* signalFile_OLD = TFile::Open("HZZlljj_SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6_Spring11-PU_S1_START311_V1G1-v1_2_opt400_0btag_optLD_ALL_looseBTag_OLD.root");
  db->add_mcFile( signalFile_OLD, "signalOLD", "No PU PDFs", 46, 3005);


  db->set_shapeNormalization();

  db->drawHisto( "QGLikelihoodJet1", "Leading Jet Q-G Likelihood");
  db->drawHisto( "QGLikelihoodJet2", "Subleading Jet Q-G Likelihood");
  db->drawHisto( "QGLikelihoodProd", "Q-G Likelihood Product");


//drawCompareQCD( db, "66_81", signalFile_NEW, signalFile_OLD );
//drawCompareQCD( db, "100_123", signalFile_NEW, signalFile_OLD );



  DrawBase* db2 = new DrawBase("db2");
  db2->set_outputdir(outputdir_str);

  db2->add_mcFile( signalFile_NEW, "signalNEW", "HZZ (400)", 38, 3004);

  TFile* qcdFile= TFile::Open("QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1.root");
  db2->add_mcFile( qcdFile, "qcd", "QCD Flat 15-1000", 46, 3005);

  db2->set_shapeNormalization();

  db2->drawHisto( "rhoPF", "Particle Flow Energy Density", "GeV", "Events", true );

  // draw th2 rho vs nvertex :: BEGIN

  TH2D* axes = new TH2D("axes", "", 10, 0.5, 15.5, 10, 0., 20.);
  axes->SetXTitle("Number of Reconstructed Vertexes");
  axes->SetYTitle("Particle Flow Energy Density [GeV]");
  
  TH2D* h2_rho_vs_nvertex = (TH2D*)qcdFile->Get("rhoPF_vs_nvertex");
  TH1D* h1_profY = h2_rho_vs_nvertex->ProfileY();

  h1_profY->SetMarkerStyle(20);
  h1_profY->SetMarkerColor(kRed);


  TPaveText* labelCMS = db2->get_labelCMS();
  TPaveText* labelSqrt = db2->get_labelSqrt();

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  c2->cd();
  
  axes->Draw();
  h2_rho_vs_nvertex->Draw("same");
  h1_profY->Draw("psame");

  labelCMS->Draw("same");
  labelSqrt->Draw("same");

  gPad->RedrawAxis();
 
  std::string canvasName = db2->get_outputdir() + "/rhoPF_vs_nvertex";
  std::string canvasName_eps = canvasName + ".eps";
  std::string canvasName_png = canvasName + ".png";
  c2->SaveAs( canvasName_eps.c_str());
  c2->SaveAs( canvasName_png.c_str());

  // draw th2 rho vs nvertex :: END

  DrawBase* db3 = new DrawBase("db3");
  db3->set_outputdir(outputdir_str);

  TFile* qcdFile_NEW  = TFile::Open("ComputeQGLikelihood_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1.root");
  db3->add_mcFile( qcdFile_NEW, "qcdNEW", "Flat10 PU MC", kRed+3, 3004);

  TFile* qcdFile_OLD = TFile::Open("ComputeQGLikelihood_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root");
  db3->add_mcFile( qcdFile_OLD, "qcdOLD", "No PU MC", kOrange, 3005);

  db3->set_shapeNormalization();
  db3->set_logx(false);

  db3->set_legendTitle( "35 < p_{T} < 43 GeV" );
  db3->drawHisto( "QGLikelihood_norms_quark_pt35_43", "Expected Likelihood for Quark Jets");
  db3->drawHisto( "QGLikelihood_norms_gluon_pt35_43", "Expected Likelihood for Gluon Jets");

  db3->set_legendTitle( "81 < p_{T} < 100 GeV" );
  db3->drawHisto( "QGLikelihood_norms_quark_pt81_100", "Expected Likelihood for Quark Jets");
  db3->drawHisto( "QGLikelihood_norms_gluon_pt81_100", "Expected Likelihood for Gluon Jets");

  db3->set_legendTitle( "230 < p_{T} < 284 GeV" );
  db3->drawHisto( "QGLikelihood_norms_quark_pt230_284", "Expected Likelihood for Quark Jets");
  db3->drawHisto( "QGLikelihood_norms_gluon_pt230_284", "Expected Likelihood for Gluon Jets");

  db3->set_legendTitle( "533 < p_{T} < 658 GeV" );
  db3->drawHisto( "QGLikelihood_norms_quark_pt533_658", "Expected Likelihood for Quark Jets");
  db3->drawHisto( "QGLikelihood_norms_gluon_pt533_658", "Expected Likelihood for Gluon Jets");

  drawEffvsRej( db3, qcdFile_OLD, qcdFile_NEW, "QGLikelihood_norms", "35_43" );
  drawEffvsRej( db3, qcdFile_OLD, qcdFile_NEW, "QGLikelihood_norms", "81_100" );
  drawEffvsRej( db3, qcdFile_OLD, qcdFile_NEW, "QGLikelihood_norms", "230_284" );
  drawEffvsRej( db3, qcdFile_OLD, qcdFile_NEW, "QGLikelihood_norms", "533_658" );


  return 0;

}  




void drawCompareQCD( DrawBase* db, const std::string& ptRange, TFile* fileSignal_NEW, TFile* fileSignal_OLD ) {


  //TFile* file_QCD = TFile::Open("ComputeQGLikelihood_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root");
  TFile* file_QCD = TFile::Open("ComputeQGLikelihood_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1.root");

  std::string histoNameQCD = "QGLikelihood_norms_quark_pt"+ptRange;

  TH1D* h1_QCD = (TH1D*)file_QCD->Get(histoNameQCD.c_str());
  //TH1D* h1_QCD_PU = (TH1D*)file_QCD_PU->Get(histoNameQCD.c_str());

  h1_QCD->SetFillColor( kOrange );
  h1_QCD->Rebin( 2 );

  std::string histoNameSignal = "QGLikelihood_"+ptRange;

  TH1D* h1_signalOLD = (TH1D*)fileSignal_OLD->Get( histoNameSignal.c_str() );
  TH1D* h1_signalNEW = (TH1D*)fileSignal_NEW->Get( histoNameSignal.c_str() );

  h1_signalOLD->SetFillColor( kRed+1 ); 
  h1_signalOLD->SetLineColor( kRed+2 ); 
  h1_signalOLD->SetLineWidth( 2 );
  h1_signalOLD->SetFillStyle( 3004 );
  h1_signalOLD->Rebin( 2 );

  h1_signalNEW->SetMarkerStyle( 20 );
  h1_signalNEW->SetMarkerSize( 1.3 );
  h1_signalNEW->Rebin( 2 );


  std::string legendTitle;
  if( ptRange=="100_123" ) legendTitle = "100 < p_{T} < 123 GeV";
  if( ptRange=="66_81" ) legendTitle = "66 < p_{T} < 81 GeV";

  TLegend* legend = new TLegend(0.22, 0.65, 0.55, 0.9, legendTitle.c_str());
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_QCD, "Expected uds with PU", "F");
  legend->AddEntry( h1_signalOLD, "Signal, No PU PDFs", "F");
  legend->AddEntry( h1_signalNEW, "Signal, #rho-Binned PDFs", "P");

  TPaveText* labelCMS = db->get_labelCMS();
  TPaveText* labelSqrt = db->get_labelSqrt();

  TH2D* axes = new TH2D("axes", "", 10, 0., 1.00001, 10, 0., 0.15);
  axes->SetXTitle("Q-G Likelihood");
  axes->SetYTitle("Normalized to Unity");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  axes->Draw("");
  legend->Draw("same");
  labelCMS->Draw("same");
  labelSqrt->Draw("same");
  h1_QCD->DrawNormalized("histosame");
  h1_signalOLD->DrawNormalized("histosame");
  h1_signalNEW->DrawNormalized("Psame");

  gPad->RedrawAxis();

  std::string canvasName = db->get_outputdir() + "/compareQCD_" + ptRange;
  std::string canvasName_eps = canvasName + ".eps";
  std::string canvasName_png = canvasName + ".png";
  c1->SaveAs(canvasName_eps.c_str());
  c1->SaveAs(canvasName_png.c_str());

}


void drawEffvsRej( DrawBase* db, TFile* qcdFile_OLD, TFile* qcdFile_NEW, const std::string& varName, const std::string& ptRange ) {


  std::string histoName_quark = varName + "_quark_pt" + ptRange;
  std::string histoName_gluon = varName + "_gluon_pt" + ptRange;

  TH1D* h1_gluon_OLD = (TH1D*)qcdFile_OLD->Get(histoName_gluon.c_str());
  TH1D* h1_quark_OLD = (TH1D*)qcdFile_OLD->Get(histoName_quark.c_str());

  TH1D* h1_gluon_NEW = (TH1D*)qcdFile_NEW->Get(histoName_gluon.c_str());
  TH1D* h1_quark_NEW = (TH1D*)qcdFile_NEW->Get(histoName_quark.c_str());

  TGraph* effRej_OLD = new TGraph(0);
  TGraph* effRej_NEW = new TGraph(0);

  int nBins = h1_gluon_OLD->GetXaxis()->GetNbins();

  for( unsigned iBin=0; iBin<nBins; ++iBin ) {

    float effQ_OLD = h1_quark_OLD->Integral( iBin+1, nBins ) / h1_quark_OLD->Integral( 1, nBins );
    float effG_OLD = h1_gluon_OLD->Integral( iBin+1, nBins ) / h1_gluon_OLD->Integral( 1, nBins );

    float effQ_NEW = h1_quark_NEW->Integral( iBin+1, nBins ) / h1_quark_NEW->Integral( 1, nBins );
    float effG_NEW = h1_gluon_NEW->Integral( iBin+1, nBins ) / h1_gluon_NEW->Integral( 1, nBins );

    effRej_OLD->SetPoint( effRej_OLD->GetN(), effQ_OLD, 1.-effG_OLD );
    effRej_NEW->SetPoint( effRej_NEW->GetN(), effQ_NEW, 1.-effG_NEW );

  } // for bins


  effRej_OLD->SetMarkerStyle(21);
  effRej_OLD->SetMarkerSize(2.);
  effRej_OLD->SetMarkerColor(kOrange);

  effRej_NEW->SetMarkerStyle(29);
  effRej_NEW->SetMarkerSize(2.);
  effRej_NEW->SetMarkerColor(kRed+3);


  TPaveText* labelCMS = db->get_labelCMS();
  TPaveText* labelSqrt = db->get_labelSqrt();

  TLine* diagonal = new TLine(0., 1., 1., 0.);

  std::string legendTitle;
  if( ptRange=="35_43" ) legendTitle = "35 < p_{T} < 43 GeV";
  else if( ptRange=="81_100" ) legendTitle = "81 < p_{T} < 100 GeV";
  else if( ptRange=="230_284" ) legendTitle = "230 < p_{T} < 284 GeV";
  else if( ptRange=="533_658" ) legendTitle = "533 < p_{T} < 658 GeV";

  TLegend* legend = new TLegend( 0.2, 0.2, 0.44, 0.45, legendTitle.c_str() );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( effRej_OLD, "No PU MC", "P" );
  legend->AddEntry( effRej_NEW, "PU MC, #rho Binning", "P" );

  TH2D* axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., 1.0001 );
  axes->SetXTitle( "Quark Efficiency" );
  axes->SetYTitle( "Gluon Rejection" );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();

  axes->Draw();

  diagonal->Draw("same");

  effRej_OLD->Draw("psame");
  effRej_NEW->Draw("psame");

  legend->Draw("same");

  labelCMS->Draw("same");
  labelSqrt->Draw("same");

  gPad->RedrawAxis();


  std::string canvasName = db->get_outputdir() + "/rejEff_"+varName+"_"+ptRange;
  std::string canvasName_eps = canvasName + ".eps";
  std::string canvasName_png = canvasName + ".png";

  c1->SaveAs(canvasName_eps.c_str());
  c1->SaveAs(canvasName_png.c_str());

}

