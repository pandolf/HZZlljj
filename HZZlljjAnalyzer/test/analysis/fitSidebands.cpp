#include <cstdlib>
#include <fstream>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "CommonTools/DrawBase.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooFermi.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

using namespace RooFit;



TH1D* drawSingleAlphaHisto( DrawBase* db, int nBTags, const std::string leptType, TFile* file_ZJets_alpgen, TFile* file_ZJets_madgraph, TFile* file_TT_TW, TFile* file_Diboson, TFile* file_Zbb );
//void fitSidebands( DrawBase* db, TH1D* h1_alpha, TFile* file_DATA, int btagCategory, const std::string& leptType );
void fitSidebands( DrawBase* db, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType );
void plotOnFrame(RooPlot* rooPlot, RooDataSet* dataset, int nBins, RooRealVar* mZZ, RooRealVar* a_exp, RooExponential* exp);


int main() {


  DrawBase* db = new DrawBase("FitSidebands");

  TFile* file_ZJets_alpgen = TFile::Open("HZZlljjRM_ZJets_alpgen_TuneZ2_Spring11_v2_optLD_looseBTags_v2_ALL.root");
  TFile* file_ZJets_madgraph = TFile::Open("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");

  db->add_mcFile( file_ZJets_alpgen, "ZJets_alpgen", "Alpgen", kOrange+1);
  db->add_mcFile( file_ZJets_madgraph, "ZJets_madgraph", "Madgraph", kRed+3);

  db->set_outputdir("FitSidebands");
  db->set_shapeNormalization();

  bool log = true;

  db->set_rebin(5);
  db->drawHisto("mZjj", "M(jj)", "GeV", "Events", log);

  db->set_rebin(20);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag", "m_(lljj)", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag", "m_(lljj)", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag", "m_(lljj)", "GeV", "Events", log);

  TFile* file_TT_TW = TFile::Open("HZZlljjRM_TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_optLD_looseBTags_v2_ALL.root");
  TFile* file_Diboson = TFile::Open("HZZlljjRM_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_optLD_looseBTags_v2_ALL.root");
  TFile* file_Zbb = TFile::Open("HZZlljjRM_ZBB_alpgen_TuneZ2_Spring11_v2_optLD_looseBTags_v2_ALL.root");

  TH1D* alpha_0btag = drawSingleAlphaHisto( db, 0, "ALL", file_ZJets_alpgen, file_ZJets_madgraph, file_TT_TW, file_Diboson, file_Zbb );
  TH1D* alpha_1btag = drawSingleAlphaHisto( db, 1, "ALL", file_ZJets_alpgen, file_ZJets_madgraph, file_TT_TW, file_Diboson, file_Zbb );
  TH1D* alpha_2btag = drawSingleAlphaHisto( db, 2, "ALL", file_ZJets_alpgen, file_ZJets_madgraph, file_TT_TW, file_Diboson, file_Zbb );

  TFile* file_DATA = TFile::Open("HZZlljjRM_DATA_EPS_optLD_looseBTags_v2_ALL.root");
  TTree* treeDATA = (TTree*)file_DATA->Get("tree_passedEvents");

  TChain* chainMC = new TChain("tree_passedEvents");
  chainMC->Add("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root/tree_passedEvents");
  chainMC->Add("HZZlljjRM_TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_optLD_looseBTags_v2_ALL.root/tree_passedEvents");
  chainMC->Add("HZZlljjRM_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_optLD_looseBTags_v2_ALL.root/tree_passedEvents");

  fitSidebands( db, (TTree*)chainMC, treeDATA, 0, "ALL" );
  fitSidebands( db, (TTree*)chainMC, treeDATA, 1, "ALL" );
  fitSidebands( db, (TTree*)chainMC, treeDATA, 2, "ALL" );

  return 0;

}



TH1D* drawSingleAlphaHisto( DrawBase* db, int nBTags, const std::string leptType, TFile* file_ZJets_alpgen, TFile* file_ZJets_madgraph, TFile* file_TT_TW, TFile* file_Diboson, TFile* file_Zbb ) {


  std::string leptType_text;
  if( leptType=="ELE" ) leptType_text = "_ELE";
  else if( leptType=="MU" ) leptType_text = "_MU";
  else if( leptType=="ALL" ) leptType_text = "";
  else {
    std::cout << "UNKNOWN LEPT TYPE: " << leptType << ". EXITING." << std::endl;
    exit(1111);
  }


  char histoName[300];
  sprintf( histoName, "mZZ_kinfit_hiMass_%dbtag%s", nBTags, leptType_text.c_str());
  char histoName_sidebands[300];
  sprintf( histoName_sidebands, "mZZ_kinfit_hiMass_sidebands_%dbtag%s", nBTags, leptType_text.c_str());

  TH1D* h1_ZJets_alpgen = (TH1D*)file_ZJets_alpgen->Get(histoName);
  TH1D* h1_ZJets_sidebands_alpgen = (TH1D*)file_ZJets_alpgen->Get(histoName_sidebands);

  TH1D* h1_ZJets_madgraph = (TH1D*)file_ZJets_madgraph->Get(histoName);
  TH1D* h1_ZJets_sidebands_madgraph = (TH1D*)file_ZJets_madgraph->Get(histoName_sidebands);

  TH1D* h1_TT_TW= (TH1D*)file_TT_TW->Get(histoName);
  TH1D* h1_TT_TW_sidebands= (TH1D*)file_TT_TW->Get(histoName_sidebands);

  TH1D* h1_Diboson= (TH1D*)file_Diboson->Get(histoName);
  TH1D* h1_Diboson_sidebands= (TH1D*)file_Diboson->Get(histoName_sidebands);

  TH1D* h1_Zbb= (TH1D*)file_Zbb->Get(histoName);
  TH1D* h1_Zbb_sidebands= (TH1D*)file_Zbb->Get(histoName_sidebands);


  int rebin = 10;

  h1_ZJets_alpgen->Rebin(rebin);
  h1_ZJets_sidebands_alpgen->Rebin(rebin);
  h1_ZJets_madgraph->Rebin(rebin);
  h1_ZJets_sidebands_madgraph->Rebin(rebin);
  h1_TT_TW->Rebin(rebin);
  h1_TT_TW_sidebands->Rebin(rebin);
  h1_Diboson->Rebin(rebin);
  h1_Diboson_sidebands->Rebin(rebin);
  h1_Zbb->Rebin(rebin);
  h1_Zbb_sidebands->Rebin(rebin);

  TH1D* h1_alpgen = new TH1D(*h1_ZJets_alpgen);
  h1_alpgen->SetName("alpgen");
  TH1D* h1_sidebands_alpgen = new TH1D(*h1_ZJets_sidebands_alpgen);
  h1_sidebands_alpgen->SetName("sidebands_alpgen");
  TH1D* h1_madgraph = new TH1D(*h1_ZJets_madgraph);
  h1_madgraph->SetName("madgraph");
  TH1D* h1_sidebands_madgraph = new TH1D(*h1_ZJets_sidebands_madgraph);
  h1_sidebands_madgraph->SetName("sidebands_madgraph");

  h1_alpgen->Add(h1_TT_TW);
  h1_alpgen->Add(h1_Diboson);
  h1_alpgen->Add(h1_Zbb);

  h1_sidebands_alpgen->Add(h1_TT_TW_sidebands);
  h1_sidebands_alpgen->Add(h1_Diboson_sidebands);
  h1_sidebands_alpgen->Add(h1_Zbb_sidebands);

  h1_madgraph->Add(h1_TT_TW);
  h1_madgraph->Add(h1_Diboson);

  h1_sidebands_madgraph->Add(h1_TT_TW_sidebands);
  h1_sidebands_madgraph->Add(h1_Diboson_sidebands);
              

  // normalize:
  h1_alpgen->Scale( 1./h1_ZJets_alpgen->Integral(1, h1_ZJets_alpgen->GetNbinsX()) );
  h1_sidebands_alpgen->Scale( 1./h1_ZJets_sidebands_alpgen->Integral(1, h1_ZJets_sidebands_alpgen->GetNbinsX()) );
  h1_madgraph->Scale( 1./h1_ZJets_madgraph->Integral(1, h1_ZJets_madgraph->GetNbinsX()) );
  h1_sidebands_madgraph->Scale( 1./h1_ZJets_sidebands_madgraph->Integral(1, h1_ZJets_sidebands_madgraph->GetNbinsX()) );

  TH1D* h1_alpha_alpgen = new TH1D(*h1_ZJets_alpgen);
  h1_alpha_alpgen->Divide(h1_ZJets_sidebands_alpgen);

  TH1D* h1_alpha_madgraph = new TH1D(*h1_ZJets_madgraph);
  h1_alpha_madgraph->Divide(h1_ZJets_sidebands_madgraph);

  h1_alpha_alpgen->SetLineColor(kOrange+1);
  h1_alpha_alpgen->SetLineWidth(2);

  h1_alpha_madgraph->SetLineColor(kRed+3);
  h1_alpha_madgraph->SetLineWidth(2);


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  char yTitle[200];
  sprintf( yTitle, "Alpha Factor / (%d GeV)", rebin );

  TH2D* axes = new TH2D("axes", "", 10, 150., 800., 10, 0., 3.);
  axes->SetXTitle("M(lljj) [GeV]");
  axes->SetYTitle(yTitle);

  axes->Draw();

  //h1_alpha_alpgen->Draw("histo");
  //h1_alpha_madgraph->Draw("histo same");
    
  h1_alpha_alpgen->Draw("same");
  h1_alpha_madgraph->Draw("same");
    
  std::string leptType_legendText;
  if( leptType=="ALL" ) leptType_legendText = "";
  else if( leptType=="ELE" ) leptType_legendText = ", Electron Channel";
  else if( leptType=="MU" ) leptType_legendText = ", Muon Channel";

  char legendTitle[100];
  sprintf( legendTitle, "%d-btag Category%s", nBTags, leptType_legendText.c_str() );

  TLegend* legend = new TLegend(0.2, 0.7, 0.5, 0.9, legendTitle);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(h1_alpha_alpgen, "Alpgen", "L");
  legend->AddEntry(h1_alpha_madgraph, "Madgraph", "L");
  legend->Draw("same");

  TPaveText* labelCMS = db->get_labelCMS();
  TPaveText* labelSqrt = db->get_labelSqrt();

  labelCMS->Draw("same");
  labelSqrt->Draw("same");



  char canvasName[400];
  sprintf( canvasName, "FitSidebands/alpha_%dbtag_%s", nBTags, leptType.c_str());

  std::string canvasName_str(canvasName);
  std::string canvasName_eps = canvasName_str + ".eps";
  std::string canvasName_png = canvasName_str + ".png";

  c1->SaveAs(canvasName_eps.c_str()); 
  c1->SaveAs(canvasName_png.c_str()); 

  delete c1;

if( nBTags==0 ) {

  TFile* file_PROVA = TFile::Open("PROVA.root", "recreate");
  file_PROVA->cd();
  h1_alpgen->Write();
  h1_sidebands_alpgen->Write();
  h1_madgraph->Write();
  h1_sidebands_madgraph->Write();
  h1_alpha_alpgen->Write();
  h1_alpha_madgraph->Write();
  file_PROVA->Close();

}
    
 

  return h1_alpha_alpgen;

 
}




//void fitSidebands( DrawBase* db, TH1D* h1_alpha, TFile* file_DATA, int btagCategory, const std::string& leptType ) {
void fitSidebands( DrawBase* db, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType ) {



  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    return;
  }
  


  char ofs_name[400];
  sprintf( ofs_name, "FitSidebands/fitresults_%dbtag.txt", btagCategory);

  ofstream ofs(ofs_name);

  //TTree* tree_passedEvents = (TTree*)file_DATA->Get("tree_passedEvents");



  char cut_base[500];
  sprintf( cut_base, "nBTags==%d %s", btagCategory, leptType_cut.c_str());
  char cut_sidebands[500];
  sprintf( cut_sidebands, "%s && ( (mZjj>60. && mZjj<75.)||(mZjj>105. && mZjj<130.) )", cut_base);
  char cut_signal[500];
  sprintf( cut_signal, "%s && ( mZjj>75. && mZjj<105. )", cut_base);
  

  float mZZ_min = 230.;
  float mZZ_max = 810.;
  int nBins = (int)(mZZ_max-mZZ_min)/20.;

  RooRealVar* eventWeight = new RooRealVar("eventWeight", "event weight", 0., 2., "");
  //RooRealVar* mZZ = new RooRealVar("mZZ", "m_{lljj}", 150., 800., "GeV");
  RooRealVar* mZZ = new RooRealVar("mZZ", "m_{lljj}", mZZ_min, mZZ_max, "GeV");
  RooRealVar* nBTags = new RooRealVar("nBTags", "number of BTags", -1., 2., "");
  RooRealVar* mZjj = new RooRealVar("mZjj", "mZjj", 60., 130., "GeV");

  RooFormulaVar* weight_lumi = new RooFormulaVar("weight_lumi", "@0*1000.", RooArgList(*eventWeight));

  //RooDataSet sidebandsMC("sidebandsMC","sidebandsMC",treeMC,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj,*weight_lumi),cut_sidebands,"weight_lumi");
  //RooDataSet signalMC("signalMC","signalMC",treeMC,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj,*weight_lumi),cut_signal,"weight_lumi");
  RooDataSet sidebandsMC("sidebandsMC","sidebandsMC",treeMC,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight");
  RooDataSet signalMC("signalMC","signalMC",treeMC,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_signal,"eventWeight");

  RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",treeDATA,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_sidebands);
  RooDataSet signalDATA("signalDATA","signalDATA",treeDATA,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_signal);



//  // ------------------------ fermi ------------------------------
//  RooRealVar cutOff("cutOff","position of fermi",191.12,0,1000);
//  //cutOff.setConstant(kTRUE);
//  RooRealVar beta("beta","width of fermi",4.698,0,50);
//  beta.setConstant(kTRUE);
//  RooFermi fermi("fermi","fermi function",*mZZ,cutOff,beta);
//
//  // -------------------- crystal ball ---------------------------
//  RooRealVar m("m","m",200.17,200,1000);
//  //m.setConstant(kTRUE);
//  RooRealVar wdth("wdth","wdth",85.73,0,1000);
//  RooRealVar n("n","n",13.067,0,100);
//  RooRealVar alpha("alpha","alpha",-1.395,-100,100); 
//  RooCBShape CB("CB","Crystal ball",*mZZ,m,wdth,alpha,n);

//  RooProdPdf background("background","background",RooArgSet(fermi,CB));
 

  // -------------------- exponential ---------------------------
  RooRealVar a_exp("a_exp","a_exp",-0.001, -1., 0.);
  RooExponential exp("exp","exp",*mZZ,a_exp);
  
  //RooFitResult *r = background.fitTo(sidebandsDATA,SumW2Error(kTRUE),InitialHesse(kTRUE),Save());
  //RooFitResult *r = background.fitTo(signalDATA,SumW2Error(kTRUE));



  // FIRST: fit MC sidebands:

  RooFitResult *r_sidebandsMC = exp.fitTo(sidebandsMC,SumW2Error(kTRUE));

  // save value of parameter:
  float a_exp_sidebandsMC = a_exp.getVal();
  ofs << "a_sidebandsMC: \t" << a_exp_sidebandsMC << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_sidebandsMC = mZZ->frame();

  //plotOnFrame(plot_sidebandsMC, &sidebandsMC, nBins, mZZ, &a_exp, &exp);

  sidebandsMC.plotOn(plot_sidebandsMC, Binning(nBins));

  RooRealVar* a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  RooExponential* exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mZZ,*a_exp_plusSigma);
  exp_plusSigma->plotOn(plot_sidebandsMC,LineColor(38),LineStyle(2));

  RooRealVar* a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  RooExponential* exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mZZ,*a_exp_minusSigma);
  exp_minusSigma->plotOn(plot_sidebandsMC,LineColor(38),LineStyle(2));

  exp.plotOn(plot_sidebandsMC, LineColor(kRed));
  sidebandsMC.plotOn(plot_sidebandsMC, Binning(nBins));

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  plot_sidebandsMC->Draw();

  char canvasName[500];
  sprintf( canvasName, "FitSidebands/mZZ_sidebandsMC_%dbtag_%s", btagCategory, leptType.c_str());
  std::string* canvasName_str = new std::string(canvasName);
  std::string canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());


  // SECOND: fit MC signal region:

  c1->Clear();
  c1->SetLogy(false);


  RooFitResult *r_signalMC = exp.fitTo(signalMC,SumW2Error(kTRUE));

  // save value of parameter:
  float a_exp_signalMC = a_exp.getVal();
  ofs << "a_signalMC: \t" << a_exp_signalMC << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_signalMC = mZZ->frame();

  //plotOnFrame(plot_signalMC, &signalMC, nBins, mZZ, &a_exp, &exp);

  signalMC.plotOn(plot_signalMC, Binning(nBins));

  a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mZZ,*a_exp_plusSigma);
  exp_plusSigma->plotOn(plot_signalMC,LineColor(38),LineStyle(2));

  a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mZZ,*a_exp_minusSigma);
  exp_minusSigma->plotOn(plot_signalMC,LineColor(38),LineStyle(2));

  exp.plotOn(plot_signalMC, LineColor(kRed));
  signalMC.plotOn(plot_signalMC, Binning(nBins));


  plot_signalMC->Draw();

  sprintf( canvasName, "FitSidebands/mZZ_signalMC_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());


  // THIRD: define alpha:

  float alpha = a_exp_sidebandsMC - a_exp_signalMC;


  // FOURTH: fit DATA sidebands:

  c1->Clear();
  c1->SetLogy(false);

  RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA);
  //RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());

  // save value of parameter:
  float a_exp_sidebandsDATA = a_exp.getVal();
  ofs << "a_sidebandsDATA: \t" << a_exp_sidebandsDATA << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_sidebandsDATA = mZZ->frame();


  sidebandsDATA.plotOn(plot_sidebandsDATA, Binning(nBins));

  a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mZZ,*a_exp_plusSigma);
  exp_plusSigma->plotOn(plot_sidebandsDATA,LineColor(38),LineStyle(2));

  a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mZZ,*a_exp_minusSigma);
  exp_minusSigma->plotOn(plot_sidebandsDATA,LineColor(38),LineStyle(2));

  exp.plotOn(plot_sidebandsDATA, LineColor(kRed));
  sidebandsDATA.plotOn(plot_sidebandsDATA, Binning(nBins));

  plot_sidebandsDATA->Draw();

  sprintf( canvasName, "FitSidebands/mZZ_sidebandsDATA_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());


  // FIFTH: scale data sidebands fit with alpha and superimpose to signal region:

  c1->Clear();
  c1->SetLogy(false);

  RooRealVar a_exp_data("a_exp_data","a_exp_data", -alpha+a_exp_sidebandsDATA, -1., 0.);
  a_exp_data.setConstant(kTRUE);
  RooExponential exp_data("exp_data","exp_data",*mZZ,a_exp_data);

  RooFitResult *r_signalDATA = exp_data.fitTo(signalDATA);

  RooPlot *plot_signalDATA = mZZ->frame();
  signalDATA.plotOn(plot_signalDATA, Binning(nBins));
  exp_data.plotOn(plot_signalDATA, LineColor(kRed));

  plot_signalDATA->Draw();

  sprintf( canvasName, "FitSidebands/mZZ_signalDATA_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  std::cout << std::endl << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << " -- " << btagCategory << " btags" << std::endl;
  std::cout << " a_exp_sidebandsMC: " << a_exp_sidebandsMC << std::endl;
  std::cout << " a_exp_signalMC: " << a_exp_signalMC << std::endl;
  std::cout << " a_exp_sidebandsDATA: " << a_exp_sidebandsDATA << std::endl;
  std::cout << " Alpha: " << alpha << std::endl;
  std::cout << " alpha/a_sidebandsDATA: " << alpha/a_exp_sidebandsDATA << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl << std::endl;

  ofs.close();

  delete eventWeight;
  delete mZZ;
  delete nBTags;
  delete mZjj;
  delete c1;
  delete plot_signalMC;
  delete plot_sidebandsMC;
  delete plot_signalDATA;
  delete plot_sidebandsDATA;
  delete r_sidebandsMC;
  delete r_signalMC;
  delete r_sidebandsDATA;
  delete r_signalDATA;
  delete canvasName_str;

}



void plotOnFrame(RooPlot* rooPlot, RooDataSet* dataset, int nBins, RooRealVar* mZZ, RooRealVar* a_exp, RooExponential* exp) {

  dataset->plotOn(rooPlot, Binning(nBins));

  RooRealVar* a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp->getVal()+a_exp->getError(), -1., 0.);
  RooExponential* exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mZZ,*a_exp_plusSigma);
  exp_plusSigma->plotOn(rooPlot,LineColor(38),LineStyle(2));

  RooRealVar* a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp->getVal()-a_exp->getError(), -1., 0.);
  RooExponential* exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mZZ,*a_exp_minusSigma);
  exp_minusSigma->plotOn(rooPlot,LineColor(38),LineStyle(2));

  exp->plotOn(rooPlot, LineColor(kRed));
  dataset->plotOn(rooPlot, Binning(nBins));

  delete a_exp_plusSigma;
  delete exp_plusSigma;

  delete a_exp_minusSigma;
  delete exp_minusSigma;

}
