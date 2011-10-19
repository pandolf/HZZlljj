#include <cstdlib>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
//#include "CommonTools/DrawBase.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "PDFs/RooFermi.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooArgusBG.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

using namespace RooFit;



TH1D* getAlphaHisto( int nBTags, const std::string leptType, TFile* file_ZJets_madgraph, TFile* file_TT_TW, TFile* file_Diboson );
void fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha=0 );
TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name );


int main() {


  TFile* file_ZJets_madgraph = TFile::Open("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  TFile* file_TT_TW = TFile::Open("HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  TFile* file_Diboson = TFile::Open("HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");

  TH1D* alpha_0btag = getAlphaHisto( 0, "ALL", file_ZJets_madgraph, file_TT_TW, file_Diboson );
  TH1D* alpha_1btag = getAlphaHisto( 1, "ALL", file_ZJets_madgraph, file_TT_TW, file_Diboson );
  TH1D* alpha_2btag = getAlphaHisto( 2, "ALL", file_ZJets_madgraph, file_TT_TW, file_Diboson );

  TFile* file_DATA = TFile::Open("HZZlljjRM_DATA_LP11_optLD_looseBTags_v2_ALL.root");
  TTree* treeDATA = (TTree*)file_DATA->Get("tree_passedEvents");

  TChain* chainMC = new TChain("tree_passedEvents");
  chainMC->Add("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root/tree_passedEvents");
  chainMC->Add("HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root/tree_passedEvents");
  chainMC->Add("HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root/tree_passedEvents");

  TTree* treeDATA_0btag = treeDATA->CopyTree("nBTags==0");
  TTree* treeDATA_1btag = treeDATA->CopyTree("nBTags==1");
  TTree* treeDATA_2btag = treeDATA->CopyTree("nBTags==2");

  TTree* treeMC_0btag = chainMC->CopyTree("nBTags==0");
  TTree* treeMC_1btag = chainMC->CopyTree("nBTags==1");
  TTree* treeMC_2btag = chainMC->CopyTree("nBTags==2");


  fitSidebands( treeMC_0btag, treeDATA_0btag, 0, "ALL", alpha_0btag );
  fitSidebands( treeMC_1btag, treeDATA_1btag, 1, "ALL", alpha_1btag );
  fitSidebands( treeMC_2btag, treeDATA_2btag, 2, "ALL", alpha_2btag );

  return 0;

}



TH1D* getAlphaHisto( int nBTags, const std::string leptType, TFile* file_ZJets_madgraph, TFile* file_TT_TW, TFile* file_Diboson ) {


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


  TH1D* h1_ZJets_madgraph = (TH1D*)file_ZJets_madgraph->Get(histoName);
  TH1D* h1_ZJets_sidebands_madgraph = (TH1D*)file_ZJets_madgraph->Get(histoName_sidebands);

  TH1D* h1_TT_TW= (TH1D*)file_TT_TW->Get(histoName);
  TH1D* h1_TT_TW_sidebands= (TH1D*)file_TT_TW->Get(histoName_sidebands);

  TH1D* h1_Diboson= (TH1D*)file_Diboson->Get(histoName);
  TH1D* h1_Diboson_sidebands= (TH1D*)file_Diboson->Get(histoName_sidebands);


  int rebin = 10;

  h1_ZJets_madgraph->Rebin(rebin);
  h1_ZJets_sidebands_madgraph->Rebin(rebin);
  h1_TT_TW->Rebin(rebin);
  h1_TT_TW_sidebands->Rebin(rebin);
  h1_Diboson->Rebin(rebin);
  h1_Diboson_sidebands->Rebin(rebin);

  TH1D* h1_madgraph = new TH1D(*h1_ZJets_madgraph);
  h1_madgraph->SetName("madgraph");
  TH1D* h1_sidebands_madgraph = new TH1D(*h1_ZJets_sidebands_madgraph);
  h1_sidebands_madgraph->SetName("sidebands_madgraph");


  h1_madgraph->Add(h1_TT_TW);
  h1_madgraph->Add(h1_Diboson);

  h1_sidebands_madgraph->Add(h1_TT_TW_sidebands);
  h1_sidebands_madgraph->Add(h1_Diboson_sidebands);
              

  // normalize:
  h1_madgraph->Scale( 1./h1_ZJets_madgraph->Integral(1, h1_ZJets_madgraph->GetNbinsX()) );
  h1_sidebands_madgraph->Scale( 1./h1_ZJets_sidebands_madgraph->Integral(1, h1_ZJets_sidebands_madgraph->GetNbinsX()) );

  TH1D* h1_alpha_madgraph = new TH1D(*h1_ZJets_madgraph);
  h1_alpha_madgraph->SetName("alpha_madgraph");
  h1_alpha_madgraph->Divide(h1_ZJets_sidebands_madgraph);
  h1_alpha_madgraph->GetXaxis()->SetRangeUser(150.,800.);


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
  legend->AddEntry(h1_alpha_madgraph, "Madgraph", "L");
  legend->Draw("same");


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




void fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha ) {



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

//tree_sidebandsMC_alpha->Write();
//tree_sidebandsDATA_alpha->Write();
//treeMC->SetName("tree_passedEventsMC");
//treeDATA->SetName("tree_passedEventsDATA");
//treeMC->Write();
//treeDATA->Write();
//file_prova->Write();
//exit(11);

  char ofsMC_name[400];
  sprintf( ofsMC_name, "FitSidebands/fitresultsMC_%dbtag.txt", btagCategory);
  ofstream ofsMC(ofsMC_name);

  char ofsDATA_name[400];
  sprintf( ofsDATA_name, "FitSidebands/fitresultsDATA_%dbtag.txt", btagCategory);
  ofstream ofsDATA(ofsDATA_name);


  char cut_base[500];
  sprintf( cut_base, "nBTags==%d %s", btagCategory, leptType_cut.c_str());
  char cut_sidebands[500];
  sprintf( cut_sidebands, "%s && ( (mZjj>60. && mZjj<75.)||(mZjj>105. && mZjj<130.) )", cut_base);
  char cut_signal[500];
  sprintf( cut_signal, "%s && ( mZjj>75. && mZjj<105. )", cut_base);
  

  //float mZZ_min = 230.;
  //float mZZ_min = (btagCategory==1) ? 150 : 170.;
  float mZZ_min = 150.;
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
  std::cout << "Correcting signal (MC): " << std::endl;
  TTree* tree_sidebandsMC_alpha = correctTreeWithAlpha( treeMC, h1_alpha, btagCategory, "sidebandsMC_alpha");
  RooDataSet sidebandsMC_alpha("sidebandsMC_alpha","sidebandsMC_alpha",tree_sidebandsMC_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");

  RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",treeDATA,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_sidebands);
  RooDataSet signalDATA("signalDATA","signalDATA",treeDATA,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_signal);
  std::cout << "Correcting signal (DATA): " << std::endl;
  TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( treeDATA, h1_alpha, btagCategory, "sidebandsDATA_alpha");
  RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");


  char alphaFileName[500];
  sprintf( alphaFileName, "alphaFile_%dbtag_%s.root", btagCategory, leptType.c_str());
  TFile* file_alpha = TFile::Open(alphaFileName, "recreate");
  file_alpha->cd();
  h1_alpha->Write();
  tree_sidebandsDATA_alpha->Write();
  tree_sidebandsMC_alpha->Write();
  file_alpha->Close();





  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",191.12,175.,220.);
  //cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",4.698,0.,30.);
  //beta.setConstant(kTRUE);
  RooFermi fermi("fermi","fermi function",*mZZ,cutOff,beta);

  // -------------------- crystal ball ---------------------------
  RooRealVar m("m","m",200.17,190.,300.);
  //m.setConstant(kTRUE);
  RooRealVar wdth("wdth","wdth",85.73,0.,200.);
  RooRealVar n("n","n",13.067,0.,100.);
  RooRealVar alpha("alpha","alpha",-1.395,-10.,10.); 
  RooCBShape CB("CB","Crystal ball",*mZZ,m,wdth,alpha,n);

  RooProdPdf background("background","background",RooArgSet(fermi,CB));
 





  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT ALPHA-CORRECTED MC SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);

  //RooFitResult *r_sidebandsMC_alpha = background.fitTo(sidebandsMC,SumW2Error(kTRUE));
  RooFitResult *r_sidebandsMC_alpha = background.fitTo(sidebandsMC_alpha,SumW2Error(kTRUE));
  //RooFitResult *r_sidebandsMC = exp.fitTo(sidebandsMC,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());

  ofsMC << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
  ofsMC << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
  ofsMC << "m " << m.getVal() << " " << m.getError() << std::endl;
  ofsMC << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
  ofsMC << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
  ofsMC << "n " << n.getVal() << " " << n.getError() << std::endl;

  ofsMC.close();

  RooPlot *plot_sidebandsMC_alpha = mZZ->frame();

  //sidebandsMC.plotOn(plot_sidebandsMC, Binning(nBins));
  sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins));

  background.plotOn(plot_sidebandsMC_alpha, LineColor(kRed));
  //sidebandsMC.plotOn(plot_sidebandsMC_alpha, Binning(nBins));
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


  //fix shape:
  cutOff.setConstant(kTRUE);
  beta.setConstant(kTRUE);
  m.setConstant(kTRUE);
  wdth.setConstant(kTRUE);
  n.setConstant(kTRUE);
  alpha.setConstant(kTRUE);


  c1->Clear();
  c1->SetLogy(false);

  RooFitResult *r_signalMC_alpha = background.fitTo(signalMC);
  //RooFitResult *r_signalMC = exp.fitTo(signalMC,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());


  RooPlot *plot_signalMC_alpha = mZZ->frame();

  signalMC.plotOn(plot_signalMC_alpha, Binning(nBins));

  background.plotOn(plot_signalMC_alpha, LineColor(kRed));
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

  // unset const-ness of alpha and wdth:
  wdth.setConstant(kFALSE);
  alpha.setConstant(kFALSE);

  RooFitResult *r_sidebandsDATA_alpha = background.fitTo(sidebandsDATA_alpha);
  //RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());

  ofsDATA << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
  ofsDATA << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
  ofsDATA << "m " << m.getVal() << " " << m.getError() << std::endl;
  ofsDATA << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
  ofsDATA << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
  ofsDATA << "n " << n.getVal() << " " << n.getError() << std::endl;

  ofsDATA.close();

  RooPlot *plot_sidebandsDATA_alpha = mZZ->frame();

  sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

  background.plotOn(plot_sidebandsDATA_alpha, LineColor(kRed));
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


/*


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT DATA SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


//// -------------------- const landau ---------------------------
//RooRealVar m_land_const_data("m_land_const_data","m_land_const_data",m_land.getVal(), 190., 300.);
//m_land_const_data.setConstant(kTRUE);
//RooRealVar s_land_const_data("s_land_const_data","s_land_const_data",s_land.getVal(), 0., 100.);
//s_land_const_data.setConstant(kTRUE);
//RooLandau landau_const_data("landau_const_data","landau_const_data",*mZZ,m_land_const_data,s_land_const_data);

//// -------------------- const exp ---------------------------
//RooRealVar a_exp_const_data("a_exp_const","a_exp_const",a_exp.getVal(), -2., 1.);
//a_exp_const_data.setConstant(kTRUE);
//RooExponential exp_const_data("exp_const_data","exp_const_data",*mZZ,a_exp_const_data);

//RooProdPdf landau_exp_const_data("landau_exp_const_data","landau_exp_const_data",RooArgSet(landau_const_data,exp_const_data));


  c1->Clear();
  c1->SetLogy(false);

  //RooFitResult *r_signalDATA_alpha = landau_exp_const_data.fitTo(signalDATA);
  //RooFitResult *r_signalDATA = exp.fitTo(signalDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());


  RooPlot *plot_signalDATA_alpha = mZZ->frame();

  //signalDATA.plotOn(plot_signalDATA_alpha, Binning(nBins));

//landau_exp_const_data.plotOn(plot_signalDATA_alpha, LineColor(kRed));
  landau_exp.plotOn(plot_signalDATA_alpha, LineColor(kRed));
  //signalDATA.plotOn(plot_signalDATA_alpha, Binning(nBins));

  plot_signalDATA_alpha->Draw();

  sprintf( canvasName, "FitSidebands/mZZ_signalDATA_alpha_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());



  char alphaFileName[500];
  sprintf( alphaFileName, "alphaFile_%dbtag_%s.root", nBTags, leptType.c_str());
  TFile* file_alpha = TFile::Open(alphaFileName, "recreate");
  file_alpha->cd();
  h1_alpha_madgraph->Write();
  file_alpha->Close();



  delete eventWeight;
  delete eventWeight_alpha;
  delete mZZ;
  delete nBTags;
  delete mZjj;
  delete c1;
 // delete plot_signalMC;
 // delete plot_sidebandsMC;
  delete plot_sidebandsMC_alpha;
  delete plot_signalMC_alpha;
 // delete plot_signalDATA;
 // delete plot_sidebandsDATA;
  delete plot_sidebandsDATA_alpha;
  delete plot_signalDATA_alpha;
 // delete r_sidebandsMC;
 // delete r_signalMC;
  delete r_sidebandsMC_alpha;
  delete r_signalMC_alpha;
//  delete r_sidebandsDATA;
  delete r_sidebandsDATA_alpha;
  //delete r_signalDATA_alpha;
 // delete r_signalDATA;
  delete canvasName_str;
  //delete treeMC;
  //delete treeDATA;
  //delete tree_sidebandsMC_alpha;
  //delete tree_sidebandsDATA_alpha;
*/


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

    if( nBTags!=btagCategory ) continue;

    int alphabin = h1_alpha->FindBin( mZZ );
    float alpha = h1_alpha->GetBinContent( alphabin );

    // alpha correction
    newWeight = eventWeight;
    if( isSidebands && mZZ>183. && mZZ<800. ) newWeight *= alpha;

    newTree->Fill();

  }

  return newTree;

}


