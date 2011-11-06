#include <stdlib.h>
#include <iostream>
#include <string>
#include "DrawBase.h"
#include "fitTools.h"

#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#include "TString.h"


bool withSignal_=true;



/*
struct BGFitParameters {

  float fermi_cutoff;
  float fermi_beta;
  float CB_m;
  float CB_wdth;
  float CB_alpha;
  float CB_n;
  float CB_theta;

  float fermi_cutoff_err;
  float fermi_beta_err;
  float CB_m_err;
  float CB_wdth_err;
  float CB_alpha_err;
  float CB_n_err;
  float CB_theta_err;

};
*/



void drawHistoWithCurve( DrawBase* db, const std::string& data_dataset, const std::string& PUType, int nbtags, std::string flags="" );




int main(int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./drawMZZ_with_curve [(string)data_dataset] [(int)signalScaleFactor] [(string)selType=\"optLD_looseBTags_v2\"]" << std::endl;
    exit(23);
  }

  std::string leptType = "ALL";

  std::string data_prefix(argv[1]);
  std::string data_dataset = "DATA_" + data_prefix;

  float signalScaleFactor=3.;
  if( argc>2 ) {
    std::string signalScaleFactor_str(argv[2]);
    signalScaleFactor = (float)atoi(signalScaleFactor_str.c_str());;
  }

  std::string selType = "optLD_looseBTags_v2";
  if( argc>3 ) {
    std::string selType_tmp(argv[3]);
    selType = selType_tmp;
  }

  //std::string ZJetsMC = "madgraph";


  TString dataset_tstr(data_prefix);
  std::string PUType = "Run2011A";
  if( data_dataset=="HR11" )
    PUType = "HR11";
  if( dataset_tstr.BeginsWith("Run2011B") )
    PUType = "2011B";


  DrawBase* db = new DrawBase("HZZlljjRM");
  db->set_pdf_aussi((bool)false);



  std::string outputdir_str = "HZZlljjRMPlots_" + data_dataset;
  if( withSignal_ ) {
    outputdir_str += "_plusSignal";
    if( signalScaleFactor!=1. ) {
      char scaleFactorText[100];
      sprintf( scaleFactorText, "_times%.0f", signalScaleFactor );
      std::string scaleFactorText_str(scaleFactorText);
      outputdir_str += scaleFactorText_str;
    }
  }
  outputdir_str += "_" + selType + "_PU" + PUType + "_" + leptType;
  db->set_outputdir(outputdir_str);


  std::string dataFileName = "HZZlljjRM_" + data_dataset + "_"+selType+"_"+leptType+".root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db->add_dataFile( dataFile, "THEDATA" );

  std::string signalFileName = "HZZlljjRM_GluGluToHToZZTo2L2Q_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1";
  signalFileName += "_" + selType;
  signalFileName += "_PU" + PUType;
  signalFileName += "_" + leptType;
  signalFileName += ".root";
  TFile* signalFile = TFile::Open(signalFileName.c_str());
  if( withSignal_ ) {
    char signalLegendText[400];
    if( signalScaleFactor==1. ) 
      sprintf( signalLegendText, "H(400)" );
    else
      sprintf( signalLegendText, "H(400) #times %.0f", signalScaleFactor);
    std::string signalLegendText_str(signalLegendText);
    db->add_mcFile( signalFile, signalScaleFactor, "H400", signalLegendText_str, kRed+2, 3004);
  }

  std::string mcZJetsFileName;
  mcZJetsFileName = "HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcZJetsFileName += "_" + selType;
  signalFileName += "_PU" + PUType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + Jets", 30, 3001);



  std::string mcVVFileName = "HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1";
  mcVVFileName += "_" + selType;
  signalFileName += "_PU" + PUType;
  mcVVFileName += "_" + leptType;
  mcVVFileName += ".root";
  TFile* mcVVFile = TFile::Open(mcVVFileName.c_str());
  db->add_mcFile( mcVVFile, "VVtoAnything_TuneZ2", "ZZ/WZ/WW", 38, 3003);


  std::string mcTTbarFileName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1";
  mcTTbarFileName += "_" + selType;
  signalFileName += "_PU" + PUType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  db->add_mcFile( mcTTbarFile, "TTtW", "tt/tW", 39, 3002);






  if( data_dataset=="DATA_Run2011A_v2_Sub2" )
    db->set_lumiNormalization(175.);
  else if( data_dataset=="DATA_1fb" )
    db->set_lumiNormalization(859.);
  else if( data_dataset=="DATA_EPS" )
    db->set_lumiNormalization(960.); 
  else if( data_dataset=="DATA_EPS_FINAL" )
    db->set_lumiNormalization(1000.); 
  else if( data_dataset=="DATA_EPS_FINAL_FULL" )
    db->set_lumiNormalization(1143.); 
  else if( data_dataset=="DATA_EPS_FINAL_FULL_plusSingleMu" )
    db->set_lumiNormalization(1143.); 
  else if( data_dataset=="DoubleElectron_Aug05ReReco" )
    db->set_lumiNormalization(227.); 
  else if( data_dataset=="DoubleMu_Aug05ReReco" )
    db->set_lumiNormalization(285.); 
  else if( data_dataset=="DATA_EPS_FINAL_plusSingleMu" )
    db->set_lumiNormalization(1143.);
  else if( data_dataset=="DATA_LP11" )
    db->set_lumiNormalization(1580.);
  else if( data_dataset=="DATA_Run2011A_FULL" )
    db->set_lumiNormalization(2100.);
  else if( data_dataset=="DATA_HR11" )
    db->set_lumiNormalization(4200.);




  bool log = true;

  db->set_rebin(20);
  db->set_xAxisMax(750.);

  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag", "m_{lljj}", "GeV", "Events", log);
  drawHistoWithCurve( db, data_prefix, PUType, 0);

  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag", "m_{lljj}", "GeV", "Events", log);
  drawHistoWithCurve( db, data_prefix, PUType, 1);

  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag", "m_{lljj}", "GeV", "Events", log);
  drawHistoWithCurve( db, data_prefix, PUType, 2);

  db->set_xAxisMax();
  db->set_rebin(1);

  db->set_legendTitle("Gluon- and 0 b-tag");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags<=0)", 30, 150., 750., "mZZ_g0btag", "m_{ZZ}", "GeV");
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0)", 30, 150., 750., "mZZ_0btag", "m_{ZZ}", "GeV");
  db->set_legendTitle("Gluon-tag Category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==-1)", 30, 150., 750., "mZZ_gtag", "m_{ZZ}", "GeV");


  // long range (up to 1300 gev):
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0)", 60, 150., 1350., "mZZ_0btag_longRange", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, 0, "longRange");

  db->set_legendTitle("1 b-tag Category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==1)", 60, 150., 1350., "mZZ_1btag_longRange", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, 1, "longRange");

  db->set_legendTitle("2 b-tag Category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==2)", 60, 150., 1350., "mZZ_2btag_longRange", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, 2, "longRange");


//db->set_legendTitle("0 b-tag Sidebands");
//db->drawHisto("mZZ_kinfit_hiMass_sidebands_0btag", "m_{lljj}", "GeV", "Events", log);

//db->set_legendTitle("1 b-tag Sidebands");
//db->drawHisto("mZZ_kinfit_hiMass_sidebands_1btag", "m_{lljj}", "GeV", "Events", log);

//db->set_legendTitle("2 b-tag Sidebands");
//db->drawHisto("mZZ_kinfit_hiMass_sidebands_2btag", "m_{lljj}", "GeV", "Events", log);



  delete db;
  db = 0;

  return 0;

}  



void drawHistoWithCurve( DrawBase* db, const std::string& data_dataset, const std::string& PUType, int nbtags, std::string flags ) {

  if( flags!="" ) flags = "_" + flags;

  TH1F::AddDirectory(kTRUE);

  // get histograms:

  std::vector< TH1D* > lastHistos_data = db->get_lastHistos_data();
  std::vector< TH1D* > lastHistos_mc   = db->get_lastHistos_mc();


  TH1D* h1_data = new TH1D(*(lastHistos_data[0]));
  float xMin = (db->get_xAxisMin()!=9999.) ? db->get_xAxisMin() : h1_data->GetXaxis()->GetXmin();
  float xMax = (db->get_xAxisMax()!=9999.) ? db->get_xAxisMax() : h1_data->GetXaxis()->GetXmax();

  // create data graph (poisson asymm errors):
  TGraphAsymmErrors* graph_data_poisson = new TGraphAsymmErrors(0);
  graph_data_poisson = fitTools::getGraphPoissonErrors(h1_data);
  graph_data_poisson->SetMarkerStyle(20);

  THStack* mc_stack = new THStack();
  for( unsigned ihisto=0; ihisto<lastHistos_mc.size(); ++ihisto ) 
    mc_stack->Add(lastHistos_mc[lastHistos_mc.size()-ihisto-1]);





  // open fit results file:
  char fitResultsFileName[200];
  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_ALL_PU%s.root", data_dataset.c_str(), nbtags, PUType.c_str());
  TFile* fitResultsFile = TFile::Open(fitResultsFileName);

  // get bg workspace:
  char workspaceName[200];
  sprintf( workspaceName, "fitWorkspace_%dbtag", nbtags );
  RooWorkspace* bgws = (RooWorkspace*)fitResultsFile->Get(workspaceName);

  // get mZZ variable:
  RooRealVar* CMS_hzz2l2q_mZZ = (RooRealVar*)bgws->var("CMS_hzz2l2q_mZZ");

  // get bg shape:
  RooAbsPdf* background = (RooAbsPdf*)bgws->pdf("background_decorr");


  // get bg normalization:
  TTree* treeSidebandsDATA_alphaCorr = (TTree*)fitResultsFile->Get("sidebandsDATA_alpha");
  TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, xMin, xMax);
  char sidebandsCut_alpha[500];
  sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d)", nbtags);
  treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
  float expBkg = h1_mZZ_sidebands_alpha->Integral();

  RooPlot *plot_MCbkg = CMS_hzz2l2q_mZZ->frame(xMin,xMax,(int)(xMax-xMin)/h1_data->GetXaxis()->GetBinWidth(1));
  background->plotOn(plot_MCbkg,RooFit::Normalization(expBkg));

  TF1* f1_bgForLegend = new TF1("bgForLegend", "[0]");
  f1_bgForLegend->SetLineColor(kBlue);
  f1_bgForLegend->SetLineWidth(3);
  

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.3*h1_data->GetMaximum());
  char yTitle[200];
  sprintf( yTitle, "Events / (%.0f GeV)", h1_data->GetXaxis()->GetBinWidth(1) );
  h2_axes->SetYTitle(yTitle);
  h2_axes->SetXTitle("m_{ZZ} [GeV]");

  float legend_xMin = 0.42;
  float legend_yMax = 0.91;
  float legend_yMin = legend_yMax - 0.07*6.;
  float legend_xMax = 0.92;

  TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax, (db->get_legendTitle()).c_str());
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( graph_data_poisson, "Data", "P");
  legend->AddEntry( f1_bgForLegend, "Expected Background", "L");
  for( unsigned imc=0; imc<lastHistos_mc.size(); ++imc ) 
    legend->AddEntry( lastHistos_mc[imc], (db->get_mcFile(imc).legendName).c_str(), "F");

  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600);
  c1->cd();

  h2_axes->Draw();
  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  legend->Draw("same");
  mc_stack->Draw("histo same");
  plot_MCbkg->Draw("same");
  graph_data_poisson->Draw("P same");

  gPad->RedrawAxis();

  char canvasName[1000];
  sprintf( canvasName, "%s/mZZ_%dbtag_withCurve%s", (db->get_outputdir()).c_str(), nbtags, flags.c_str() );
  std::string canvasName_str(canvasName);
  std::string canvasName_eps = canvasName_str+".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->Clear();
  c1->SetLogy();

  TLegend* legend_log = new TLegend(0.6, legend_yMin, legend_xMax, legend_yMax, (db->get_legendTitle()).c_str());
  legend_log->SetTextSize(0.04);
  legend_log->SetFillColor(0);
  legend_log->AddEntry( graph_data_poisson, "Data", "P");
  legend_log->AddEntry( f1_bgForLegend, "Exp. BG", "L");
  for( unsigned imc=0; imc<lastHistos_mc.size(); ++imc ) 
    legend_log->AddEntry( lastHistos_mc[imc], (db->get_mcFile(imc).legendName).c_str(), "F");



  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.05, 200.*h1_data->GetMaximum());
  h2_axes_log->SetYTitle(yTitle);
  h2_axes_log->SetXTitle("m_{ZZ} [GeV]");
  h2_axes_log->GetYaxis()->SetNoExponent();

  h2_axes_log->Draw();
  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  legend_log->Draw("same");
  mc_stack->Draw("histo same");
  plot_MCbkg->Draw("same");
  graph_data_poisson->Draw("P same");

  gPad->RedrawAxis();

  std::string canvasName_log_eps = canvasName_str+"_log.eps";
  c1->SaveAs(canvasName_log_eps.c_str());

}



