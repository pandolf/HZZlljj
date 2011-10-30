#include <stdlib.h>
#include <iostream>
#include <string>
#include "DrawBase.h"
#include "fitTools.h"

#include "RooRealVar.h"
#include "RooFermi.h"
#include "RooCB.h"
#include "RooProdPdf.h"
#include "RooPlot.h"


bool withSignal_=true;



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



void drawHistoWithCurve( DrawBase* db, const std::string& data_dataset, int nbtags, std::string flags="" );
BGFitParameters get_BGFitParameters( const std::string& dataset, int nbtags );




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
  outputdir_str += "_" + selType + "_" + leptType;
  db->set_outputdir(outputdir_str);


  std::string dataFileName = "HZZlljjRM_" + data_dataset + "_"+selType+"_"+leptType+".root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db->add_dataFile( dataFile, "THEDATA" );

  std::string signalFileName = "HZZlljjRM_GluGluToHToZZTo2L2Q_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1";
  signalFileName += "_" + selType;
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
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + Jets", 30, 3001);



  std::string mcVVFileName = "HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1";
  mcVVFileName += "_" + selType;
  mcVVFileName += "_" + leptType;
  mcVVFileName += ".root";
  TFile* mcVVFile = TFile::Open(mcVVFileName.c_str());
  db->add_mcFile( mcVVFile, "VVtoAnything_TuneZ2", "ZZ/WZ/WW", 38, 3003);


  std::string mcTTbarFileName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1";
  mcTTbarFileName += "_" + selType;
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
  drawHistoWithCurve( db, data_prefix, 0);

  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag", "m_{lljj}", "GeV", "Events", log);
  drawHistoWithCurve( db, data_prefix, 1);

  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag", "m_{lljj}", "GeV", "Events", log);
  drawHistoWithCurve( db, data_prefix, 2);

  db->set_xAxisMax();
  db->set_rebin(1);

  // long range (up to 1300 gev):
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto_fromTree("tree_passedEvents", "mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0)", 60, 150., 1350., "mZZ_0btag_longRange", "m_{ZZ} [GeV]", "GeV");
  drawHistoWithCurve( db, data_prefix, 0, "longRange");

  db->set_legendTitle("1 b-tag Category");
  db->drawHisto_fromTree("tree_passedEvents", "mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==1)", 60, 150., 1350., "mZZ_1btag_longRange", "m_{ZZ} [GeV]", "GeV");
  drawHistoWithCurve( db, data_prefix, 1, "longRange");

  db->set_legendTitle("2 b-tag Category");
  db->drawHisto_fromTree("tree_passedEvents", "mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==2)", 60, 150., 1350., "mZZ_2btag_longRange", "m_{ZZ} [GeV]", "GeV");
  drawHistoWithCurve( db, data_prefix, 2, "longRange");


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



void drawHistoWithCurve( DrawBase* db, const std::string& data_dataset, int nbtags, std::string flags ) {

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




  // define mZZ variable
  RooRealVar CMS_hzz2l2q_mZZ("CMS_hzz2l2q_mZZfull", "zz inv mass", xMin, xMax );


  // define background PDF:
  BGFitParameters bgfp = get_BGFitParameters( data_dataset, nbtags );

  RooRealVar fermi_cutoff("fermi_cutoff", "position of fermi", bgfp.fermi_cutoff, 0., 1000.);
  fermi_cutoff.setConstant(kTRUE);
  RooRealVar fermi_beta("fermi_beta", "width of fermi", bgfp.fermi_beta, 0., 50.);
  fermi_beta.setConstant(kTRUE);

  RooFermi fermi_BKG("fermi_BKG", "fermi function", CMS_hzz2l2q_mZZ, fermi_cutoff, fermi_beta);


  RooRealVar m("m", "m", bgfp.CB_m, 100., 1000.);
  m.setConstant(kTRUE);
  RooRealVar wdth("wdth", "wdth", bgfp.CB_wdth, 0., 1000.);
  wdth.setConstant(kTRUE);
  RooRealVar n("n", "n", bgfp.CB_n, 0., 1001.);
  n.setConstant(kTRUE);
  RooRealVar alpha("alpha", "alpha", bgfp.CB_alpha, -100., 100.);
  alpha.setConstant(kTRUE);
  RooRealVar theta("theta", "theta", bgfp.CB_theta, -3.14159, 3.14159); 
  theta.setConstant(kTRUE);
  

  RooCB CB_BKG("CB_BKG", "Crystal ball", CMS_hzz2l2q_mZZ, m, wdth, alpha, n, theta);
  RooProdPdf background("background", "background", RooArgSet(fermi_BKG,CB_BKG));


  //get expected bg normalization:
  char alphaFileName[200];
  sprintf( alphaFileName, "alphaFile_%s_%dbtag_ALL.root", data_dataset.c_str(), nbtags);
  TFile* alphaFile = TFile::Open(alphaFileName);
  TTree* treeSidebandsDATA_alphaCorr = (TTree*)alphaFile->Get("sidebandsDATA_alpha");
  TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, 150., 800.);
  char sidebandsCut_alpha[500];
  sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d)", nbtags);
  treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "mZZ", sidebandsCut_alpha);
  float expBkg = h1_mZZ_sidebands_alpha->Integral();

  RooPlot *plot_MCbkg = CMS_hzz2l2q_mZZ.frame(xMin,xMax,(int)(xMax-xMin)/h1_data->GetXaxis()->GetBinWidth(1));
  background.plotOn(plot_MCbkg,RooFit::Normalization(expBkg));

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.3*h1_data->GetMaximum());
  char yTitle[200];
  sprintf( yTitle, "Events / (%.0f GeV)", h1_data->GetXaxis()->GetBinWidth(1) );
  h2_axes->SetYTitle(yTitle);
  h2_axes->SetXTitle("m_{ZZ} [GeV]");

  float legend_xMin = 0.9*0.63;
  float legend_yMax = 0.91;
  float legend_yMin = legend_yMax - 0.07*5.;
  float legend_xMax = 0.92;

  TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax, (db->get_legendTitle()).c_str());
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( graph_data_poisson, "Data", "P");
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


  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.05, 50.*h1_data->GetMaximum());
  h2_axes_log->SetYTitle(yTitle);
  h2_axes_log->SetXTitle("m_{ZZ} [GeV]");
  h2_axes_log->GetYaxis()->SetNoExponent();

  h2_axes_log->Draw();
  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  legend->Draw("same");
  mc_stack->Draw("histo same");
  plot_MCbkg->Draw("same");
  graph_data_poisson->Draw("P same");

  gPad->RedrawAxis();

  std::string canvasName_log_eps = canvasName_str+"_log.eps";
  c1->SaveAs(canvasName_log_eps.c_str());

}



BGFitParameters get_BGFitParameters( const std::string& dataset, int nbtags ) {

  // read background parametrizations from fit results file
  char fitResultsFile[900];
  sprintf( fitResultsFile, "FitSidebands_%s/fitresultsDATA_%dbtag.txt", dataset.c_str(), nbtags);
  
  ifstream ifs(fitResultsFile);
  ifs.clear();
  ifs.seekg(0);

  BGFitParameters bgfp;

  while( ifs.good() ) {

    std::string varName;
    float value, error;
    ifs >> varName >> value >> error;

    if( varName=="beta" ) {
      bgfp.fermi_beta = value;
      bgfp.fermi_beta_err = error;
    }
    if( varName=="cutOff" ) {
      bgfp.fermi_cutoff = value;
      bgfp.fermi_cutoff_err = error;
    }
    if( varName=="m" ) {
      bgfp.CB_m = value;
      bgfp.CB_m_err = error;
    }
    if( varName=="n" ) {
      bgfp.CB_n = value;
      bgfp.CB_n_err = error;
    }
    if( varName=="alpha_rot" ) {
      bgfp.CB_alpha = value;
      bgfp.CB_alpha_err = error;
    }
    if( varName=="wdth_rot" ) {
      bgfp.CB_wdth = value;
      bgfp.CB_wdth_err = error;
    }
    if( varName=="theta_best" ) {
      bgfp.CB_theta = value;
      bgfp.CB_theta_err = error;
    }

  } // while ifs fitresults

  return bgfp;

}


