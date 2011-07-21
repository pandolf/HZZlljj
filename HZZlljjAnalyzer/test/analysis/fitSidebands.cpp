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
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooCBShape.h"
#include "RooCruijffPdf.h"
#include "RooExponential.h"
#include "RooArgusBG.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

using namespace RooFit;



TH1D* drawSingleAlphaHisto( DrawBase* db, int nBTags, const std::string leptType, TFile* file_ZJets_alpgen, TFile* file_ZJets_madgraph, TFile* file_TT_TW, TFile* file_Diboson, TFile* file_Zbb );
//void fitSidebands( DrawBase* db, TH1D* h1_alpha, TFile* file_DATA, int btagCategory, const std::string& leptType );
void fitSidebands( DrawBase* db, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha=0 );
TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name );


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
  db->drawHisto("mZjj", "m_{jj}", "GeV", "Events", log);

  db->set_rebin(20);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag", "m_{lljj}", "GeV", "Events", log);

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

  TTree* treeDATA_0btag = treeDATA->CopyTree("nBTags==0");
  TTree* treeDATA_1btag = treeDATA->CopyTree("nBTags==1");
  TTree* treeDATA_2btag = treeDATA->CopyTree("nBTags==2");

  TTree* treeMC_0btag = chainMC->CopyTree("nBTags==0");
  TTree* treeMC_1btag = chainMC->CopyTree("nBTags==1");
  TTree* treeMC_2btag = chainMC->CopyTree("nBTags==2");


  fitSidebands( db, treeMC_0btag, treeDATA_0btag, 0, "ALL", alpha_0btag );
  fitSidebands( db, treeMC_1btag, treeDATA_1btag, 1, "ALL", alpha_1btag );
  fitSidebands( db, treeMC_2btag, treeDATA_2btag, 2, "ALL", alpha_2btag );

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


  // default: return madgraph
  return h1_alpha_madgraph;

 
}




//void fitSidebands( DrawBase* db, TH1D* h1_alpha, TFile* file_DATA, int btagCategory, const std::string& leptType ) {
void fitSidebands( DrawBase* db, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha ) {



  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    return;
  }
  
//TFile* file_prova = TFile::Open("PROVA.root", "recreate");
//file_prova->cd();
  std::cout << "Correcting signal (MC): " << std::endl;
  TTree* tree_sidebandsMC_alpha = correctTreeWithAlpha( treeMC, h1_alpha, btagCategory, "sidebandsMC_alpha");
  std::cout << "Correcting signal (DATA): " << std::endl;
  TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( treeDATA, h1_alpha, btagCategory, "sidebandsDATA_alpha");

//tree_sidebandsMC_alpha->Write();
//tree_sidebandsDATA_alpha->Write();
//treeMC->SetName("tree_passedEventsMC");
//treeDATA->SetName("tree_passedEventsDATA");
//treeMC->Write();
//treeDATA->Write();
//file_prova->Write();
//exit(11);

  char ofs_name[400];
  sprintf( ofs_name, "FitSidebands/fitresults_%dbtag.txt", btagCategory);

  ofstream ofs(ofs_name);



  char cut_base[500];
  sprintf( cut_base, "nBTags==%d %s", btagCategory, leptType_cut.c_str());
  char cut_sidebands[500];
  sprintf( cut_sidebands, "%s && ( (mZjj>60. && mZjj<75.)||(mZjj>105. && mZjj<130.) )", cut_base);
  char cut_signal[500];
  sprintf( cut_signal, "%s && ( mZjj>75. && mZjj<105. )", cut_base);
  

  //float mZZ_min = 230.;
  float mZZ_min = 190.;
  //float mZZ_max = 300.;
  float mZZ_max = 810.;
  float binWidth = 20.;
  int nBins = (int)(mZZ_max-mZZ_min)/binWidth;

  RooRealVar* eventWeight = new RooRealVar("eventWeight", "event weight", 0., 2., "");
  RooRealVar* eventWeight_alpha = new RooRealVar("eventWeight_alpha", "event weight (alpha corrected)", 0., 2., "");
  //RooRealVar* mZZ = new RooRealVar("mZZ", "m_{lljj}", 150., 800., "GeV");
  RooRealVar* mZZ = new RooRealVar("mZZ", "m_{lljj}", mZZ_min, mZZ_max, "GeV");
  RooRealVar* nBTags = new RooRealVar("nBTags", "number of BTags", -1., 2., "");
  RooRealVar* mZjj = new RooRealVar("mZjj", "mZjj", 60., 130., "GeV");

  RooFormulaVar* weight_lumi = new RooFormulaVar("weight_lumi", "@0*1000.", RooArgList(*eventWeight));

  //RooDataSet sidebandsMC("sidebandsMC","sidebandsMC",treeMC,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj,*weight_lumi),cut_sidebands,"weight_lumi");
  //RooDataSet signalMC("signalMC","signalMC",treeMC,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj,*weight_lumi),cut_signal,"weight_lumi");
  RooDataSet sidebandsMC("sidebandsMC","sidebandsMC",treeMC,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight");
  RooDataSet signalMC("signalMC","signalMC",treeMC,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_signal,"eventWeight");
  RooDataSet sidebandsMC_alpha("sidebandsMC_alpha","sidebandsMC_alpha",tree_sidebandsMC_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");

  RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",treeDATA,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_sidebands);
  RooDataSet signalDATA("signalDATA","signalDATA",treeDATA,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_signal);
  RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");





  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",191.12,0,1000);
  //cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",4.698,0,50);
  //beta.setConstant(kTRUE);
  RooFermi fermi("fermi","fermi function",*mZZ,cutOff,beta);
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
  RooRealVar a_exp("a_exp","a_exp",-0.001, -2., 1.);
  RooExponential exp("exp","exp",*mZZ,a_exp);
  

  // -------------------- gauss ---------------------------
  RooRealVar m_gaus("m_gaus","m_gaus",220., 150., 300.);
  RooRealVar s_gaus("s_gaus","s_gaus",50., 0., 100.);
  RooGaussian gaus("gaus","gaus",*mZZ,m_gaus,s_gaus);
  

  // -------------------- landau ---------------------------
  RooRealVar m_land("m_land","m_land",220., 190., 300.);
  RooRealVar s_land("s_land","s_land",50., 0., 100.);
  RooLandau landau("landau","landau",*mZZ,m_land,s_land);
  

  // -------------------- ARGUS ---------------------------
  RooRealVar m0_arg("m0_arg","m0_arg", 200., 150., 400.);
  RooRealVar c_arg("c_arg","c_arg", 0.1, -100., 100.);
  RooArgusBG argus("argus","argus",*mZZ,m0_arg,c_arg);
  


  // -------------------- CRUIJFF ---------------------------
  RooRealVar m0_cru("m0_cru","m0_cru", 250., 150., 300.);
  RooRealVar sL_cru("sL_cru","sL_cru", 100., 0., 200.);
  sL_cru.setConstant(kTRUE);
  RooRealVar aL_cru("aL_cru","aL_cru", 0., 0., 1.);
  aL_cru.setConstant(kTRUE);
  RooRealVar sR_cru("sR_cru","sR_cru", 10., 0., 100.);
  RooRealVar aR_cru("aR_cru","aR_cru", 0.001, -2., 2.);
  RooCruijffPdf cruijff("cruijff","cruijff",*mZZ,m0_cru,sL_cru,sR_cru,aL_cru,aR_cru);
  
  RooProdPdf argus_exp("argus_exp","argus_exp",RooArgSet(argus,exp));
  RooProdPdf fermi_exp("fermi_exp","fermi_exp",RooArgSet(fermi,exp));
  RooProdPdf landau_exp("landau_exp","landau_exp",RooArgSet(landau,exp));
//  RooProdPdf gaus_exp("gaus_exp","gaus_exp",RooArgSet(gaus,exp));




/*
  // FIRST: fit MC sidebands:

  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIRST STEP: FIT MC SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;

  RooFitResult *r_sidebandsMC = exp.fitTo(sidebandsMC,SumW2Error(kTRUE));
  RooFitResult *r_sidebandsMC2 = landau_exp.fitTo(sidebandsMC,SumW2Error(kTRUE));
//RooFitResult *r_sidebandsMC2 = cruijff.fitTo(sidebandsMC,SumW2Error(kTRUE));
//RooFitResult *r_sidebandsMC3 = gaus.fitTo(sidebandsMC,SumW2Error(kTRUE));
//RooFitResult *r_sidebandsMC4 = landau.fitTo(sidebandsMC,SumW2Error(kTRUE));
  TF1* f1_landau_exp = landau_exp.asTF(*mZZ);
  f1_landau_exp->SetLineColor(kPink);
  f1_landau_exp->SetLineStyle(2);
  f1_landau_exp->SetLineWidth(2);

  // save value/error of parameter:
  float a_exp_sidebandsMC = a_exp.getVal();
  float a_exp_sidebandsMC_error = a_exp.getError();
  ofs << "a_sidebandsMC: \t" << a_exp_sidebandsMC << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_sidebandsMC = mZZ->frame();

  //plotOnFrame(plot_sidebandsMC, &sidebandsMC, nBins, mZZ, &a_exp, &exp);

  sidebandsMC.plotOn(plot_sidebandsMC, Binning(nBins));

  RooRealVar* a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  RooExponential* exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mZZ,*a_exp_plusSigma);
  //exp_plusSigma->plotOn(plot_sidebandsMC,LineColor(38),LineStyle(2));

  RooRealVar* a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  RooExponential* exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mZZ,*a_exp_minusSigma);
  //exp_minusSigma->plotOn(plot_sidebandsMC,LineColor(38),LineStyle(2));

  //exp.plotOn(plot_sidebandsMC, LineColor(kRed));
  //cruijff.plotOn(plot_sidebandsMC, LineColor(kBlue));
  //gaus.plotOn(plot_sidebandsMC, LineColor(kBlue));
  //landau.plotOn(plot_sidebandsMC, LineColor(kBlack));
  landau_exp.plotOn(plot_sidebandsMC, LineColor(kRed));
  sidebandsMC.plotOn(plot_sidebandsMC, Binning(nBins));

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  plot_sidebandsMC->Draw();
  f1_landau_exp->Draw("same");

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

  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  SECOND STEP: FIT MC SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;

  c1->Clear();
  c1->SetLogy(false);


  RooFitResult *r_signalMC = exp.fitTo(signalMC,SumW2Error(kTRUE));
  RooFitResult *r_signalMC2 = landau_exp.fitTo(signalMC,SumW2Error(kTRUE));


  // save value/error of parameter:
  float a_exp_signalMC = a_exp.getVal();
  float a_exp_signalMC_error = a_exp.getError();
  ofs << "a_signalMC: \t" << a_exp_signalMC << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_signalMC = mZZ->frame();

  //plotOnFrame(plot_signalMC, &signalMC, nBins, mZZ, &a_exp, &exp);

  signalMC.plotOn(plot_signalMC, Binning(nBins));

  a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mZZ,*a_exp_plusSigma);
  //exp_plusSigma->plotOn(plot_signalMC,LineColor(38),LineStyle(2));

  a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mZZ,*a_exp_minusSigma);
  //exp_minusSigma->plotOn(plot_signalMC,LineColor(38),LineStyle(2));

  //exp.plotOn(plot_signalMC, LineColor(kRed));
  landau_exp.plotOn(plot_signalMC, LineColor(kGreen));
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

  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  THIRD STEP: FIT DATA SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);

  RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA);
  RooFitResult *r_sidebandsDATA2 = landau_exp.fitTo(sidebandsDATA);
  //RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());

  // save value/error of parameter:
  float a_exp_sidebandsDATA = a_exp.getVal();
  float a_exp_sidebandsDATA_error = a_exp.getError();
  ofs << "a_sidebandsDATA: \t" << a_exp_sidebandsDATA << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_sidebandsDATA = mZZ->frame();


  sidebandsDATA.plotOn(plot_sidebandsDATA, Binning(nBins));

  a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mZZ,*a_exp_plusSigma);
  //exp_plusSigma->plotOn(plot_sidebandsDATA,LineColor(38),LineStyle(2));

  a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mZZ,*a_exp_minusSigma);
  //exp_minusSigma->plotOn(plot_sidebandsDATA,LineColor(38),LineStyle(2));

  //exp.plotOn(plot_sidebandsDATA, LineColor(kRed));
  landau_exp.plotOn(plot_sidebandsDATA, LineColor(kGreen));
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


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FOURTH STEP: FIT DATA SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);

  RooRealVar a_exp_data("a_exp_data","a_exp_data", -alpha+a_exp_sidebandsDATA, -1., 0.);
  a_exp_data.setConstant(kTRUE);
  RooExponential exp_data("exp_data","exp_data",*mZZ,a_exp_data);

  float a_exp_data_error = sqrt( a_exp_sidebandsMC_error*a_exp_sidebandsMC_error + 
                                 a_exp_signalMC_error*a_exp_signalMC_error +
                                 a_exp_sidebandsDATA_error*a_exp_sidebandsDATA_error );


  ofs << "a_signalDATA: \t" << a_exp_data.getVal() << "\t+-" << a_exp_data_error << std::endl;

  RooFitResult *r_signalDATA = exp_data.fitTo(signalDATA);

  RooPlot *plot_signalDATA = mZZ->frame();
  signalDATA.plotOn(plot_signalDATA, Binning(nBins));

  a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp_data.getVal()+a_exp_data_error, -1., 0.);
  exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mZZ,*a_exp_plusSigma);
  //exp_plusSigma->plotOn(plot_signalDATA,LineColor(38),LineStyle(2));

  a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp_data.getVal()-a_exp_data_error, -1., 0.);
  exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mZZ,*a_exp_minusSigma);
  //exp_minusSigma->plotOn(plot_signalDATA,LineColor(38),LineStyle(2));

  //exp_data.plotOn(plot_signalDATA, LineColor(kRed));
  signalDATA.plotOn(plot_signalDATA, Binning(nBins));

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

*/

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT ALPHA-CORRECTED MC SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);

  RooFitResult *r_sidebandsMC_alpha = landau_exp.fitTo(sidebandsMC_alpha);
  //RooFitResult *r_sidebandsMC = exp.fitTo(sidebandsMC,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());


  RooPlot *plot_sidebandsMC_alpha = mZZ->frame();

  sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins));

  landau_exp.plotOn(plot_sidebandsMC_alpha, LineColor(kRed));
  sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins));

  plot_sidebandsMC_alpha->Draw();

  char canvasName[400];
  sprintf( canvasName, "FitSidebands/mZZ_sidebandsMC_alpha_%dbtag_%s", btagCategory, leptType.c_str());
  std::string* canvasName_str = new std::string(canvasName);
  std::string canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());




  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT MC SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;



  // -------------------- const landau ---------------------------
  RooRealVar m_land_const("m_land_const","m_land_const",m_land.getVal(), 190., 300.);
  m_land_const.setConstant(kTRUE);
  RooRealVar s_land_const("s_land_const","s_land_const",s_land.getVal(), 0., 100.);
  s_land_const.setConstant(kTRUE);
  RooLandau landau_const("landau_const","landau_const",*mZZ,m_land_const,s_land_const);

  // -------------------- const exp ---------------------------
  RooRealVar a_exp_const("a_exp_const","a_exp_const",a_exp.getVal(), -2., 1.);
  a_exp_const.setConstant(kTRUE);
  RooExponential exp_const("exp_const","exp_const",*mZZ,a_exp_const);

  RooProdPdf landau_exp_const("landau_exp_const","landau_exp_const",RooArgSet(landau_const,exp_const));



  c1->Clear();
  c1->SetLogy(false);

  RooFitResult *r_signalMC_alpha = landau_exp_const.fitTo(signalMC);
  //RooFitResult *r_signalMC = exp.fitTo(signalMC,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());


  RooPlot *plot_signalMC_alpha = mZZ->frame();

  signalMC.plotOn(plot_signalMC_alpha, Binning(nBins));

  landau_exp_const.plotOn(plot_signalMC_alpha, LineColor(kRed));
  signalMC.plotOn(plot_signalMC_alpha, Binning(nBins));

  plot_signalMC_alpha->Draw();

  sprintf( canvasName, "FitSidebands/mZZ_signalMC_alpha_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());




  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT ALPHA-CORRECTED DATA SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);

  RooFitResult *r_sidebandsDATA_alpha = landau_exp.fitTo(sidebandsDATA_alpha);
  //RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());


  RooPlot *plot_sidebandsDATA_alpha = mZZ->frame();

  sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

  landau_exp.plotOn(plot_sidebandsDATA_alpha, LineColor(kRed));
  sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

  plot_sidebandsDATA_alpha->Draw();

  sprintf( canvasName, "FitSidebands/mZZ_sidebandsDATA_alpha_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());





  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT DATA SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  // -------------------- const landau ---------------------------
  RooRealVar m_land_const_data("m_land_const_data","m_land_const_data",m_land.getVal(), 190., 300.);
  m_land_const_data.setConstant(kTRUE);
  RooRealVar s_land_const_data("s_land_const_data","s_land_const_data",s_land.getVal(), 0., 100.);
  s_land_const_data.setConstant(kTRUE);
  RooLandau landau_const_data("landau_const_data","landau_const_data",*mZZ,m_land_const_data,s_land_const_data);

  // -------------------- const exp ---------------------------
  RooRealVar a_exp_const_data("a_exp_const","a_exp_const",a_exp.getVal(), -2., 1.);
  a_exp_const_data.setConstant(kTRUE);
  RooExponential exp_const_data("exp_const_data","exp_const_data",*mZZ,a_exp_const_data);

  RooProdPdf landau_exp_const_data("landau_exp_const_data","landau_exp_const_data",RooArgSet(landau_const_data,exp_const_data));


  c1->Clear();
  c1->SetLogy(false);

  RooFitResult *r_signalDATA_alpha = landau_exp_const_data.fitTo(signalDATA);
  //RooFitResult *r_signalDATA = exp.fitTo(signalDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());


  RooPlot *plot_signalDATA_alpha = mZZ->frame();

  signalDATA.plotOn(plot_signalDATA_alpha, Binning(nBins));

  landau_exp_const_data.plotOn(plot_signalDATA_alpha, LineColor(kRed));
  signalDATA.plotOn(plot_signalDATA_alpha, Binning(nBins));

  plot_signalDATA_alpha->Draw();

  sprintf( canvasName, "FitSidebands/mZZ_signalDATA_alpha_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());








  delete eventWeight;
  delete eventWeight_alpha;
  delete mZZ;
  delete nBTags;
  delete mZjj;
  delete c1;
 // delete plot_signalMC;
 // delete plot_sidebandsMC;
  delete plot_sidebandsMC_alpha;
 // delete plot_signalDATA;
 // delete plot_sidebandsDATA;
  delete plot_sidebandsDATA_alpha;
 // delete r_sidebandsMC;
 // delete r_signalMC;
  delete r_sidebandsMC_alpha;
//  delete r_sidebandsDATA;
  delete r_sidebandsDATA_alpha;
 // delete r_signalDATA;
  delete canvasName_str;
  //delete treeMC;
  //delete treeDATA;
  //delete tree_sidebandsMC_alpha;
  //delete tree_sidebandsDATA_alpha;

}



TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name ) {

  Int_t leptType;
  tree->SetBranchAddress( "leptType", &leptType );
  Int_t nBTags;
  tree->SetBranchAddress( "nBTags", &nBTags );
  Float_t mZZ;
  tree->SetBranchAddress( "mZZ", &mZZ );
  Float_t mZjj;
  tree->SetBranchAddress( "mZjj", &mZjj );
  Float_t eventWeight;
  tree->SetBranchAddress( "eventWeight", &eventWeight );
  Bool_t isSidebands;
  tree->SetBranchAddress( "isSidebands", &isSidebands );


  TTree* newTree = tree->CloneTree(0);
  newTree->SetName(name.c_str());

  Float_t newWeight;
  newTree->Branch( "eventWeight_alpha", &newWeight, "newWeight/F" );

  
  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry( iEntry );
    if( (iEntry % 10000)==0 ) std::cout << "Entry: " << iEntry << "/" << nentries << std::endl;

    if( !isSidebands ) continue;
    if( nBTags!=btagCategory ) continue;
    if( mZZ>800. ) continue;

    int alphabin = h1_alpha->FindBin( mZZ );
    float alpha = h1_alpha->GetBinContent( alphabin );

    // alpha correction
    newWeight = eventWeight*alpha;

    newTree->Fill();

  }

  return newTree;

}


