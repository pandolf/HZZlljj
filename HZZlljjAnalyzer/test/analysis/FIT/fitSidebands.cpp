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
#include "TMatrixDSym.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooFermi.h"
#include "RooGaussian.h"
#include "RooCB.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

using namespace RooFit;



TH1D* getAlphaHisto( int nBTags, const std::string leptType_str, TTree* treeMC );
//TH1D* getAlphaHisto( int nBTags, const std::string leptType, TFile* file_ZJets_madgraph, TFile* file_TT_TW, TFile* file_Diboson );
void fitSidebands( const std::string& dataset, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha=0 );
TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name );


int main( int argc, char* argv[] ) {

  std::string dataset = "LP11";
  if( argc==2 ) {
    std::string dataset_str(argv[1]);
    dataset = dataset_str;
  } 

  TFile* file_ZJets_madgraph = TFile::Open("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  TFile* file_TT_TW = TFile::Open("HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  TFile* file_Diboson = TFile::Open("HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");

  std::string datafileName = "HZZlljjRM_DATA_" + dataset + "_optLD_looseBTags_v2_ALL.root";
  TFile* file_DATA = TFile::Open(datafileName.c_str());
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

  //TH1D* alpha_0btag = getAlphaHisto( 0, "ALL", file_ZJets_madgraph, file_TT_TW, file_Diboson );
  //TH1D* alpha_1btag = getAlphaHisto( 1, "ALL", file_ZJets_madgraph, file_TT_TW, file_Diboson );
  //TH1D* alpha_2btag = getAlphaHisto( 2, "ALL", file_ZJets_madgraph, file_TT_TW, file_Diboson );

  TH1D* alpha_0btag = getAlphaHisto( 0, "ALL", treeMC_0btag );
  TH1D* alpha_1btag = getAlphaHisto( 1, "ALL", treeMC_1btag );
  TH1D* alpha_2btag = getAlphaHisto( 2, "ALL", treeMC_2btag );

  fitSidebands( dataset, treeMC_0btag, treeDATA_0btag, 0, "ALL", alpha_0btag );
  fitSidebands( dataset, treeMC_1btag, treeDATA_1btag, 1, "ALL", alpha_1btag );
  fitSidebands( dataset, treeMC_2btag, treeDATA_2btag, 2, "ALL", alpha_2btag );

  return 0;

}



TH1D* getAlphaHisto( int btagCategory, const std::string leptType_str, TTree* treeMC ) {


  std::string leptType_text;
  if( leptType_str=="ELE" ) leptType_text = "_ELE";
  else if( leptType_str=="MU" ) leptType_text = "_MU";
  else if( leptType_str=="ALL" ) leptType_text = "";
  else {
    std::cout << "UNKNOWN LEPT TYPE: " << leptType_str << ". EXITING." << std::endl;
    exit(1111);
  }

  float mZZ;
  float eventWeight;
  int nBTags;
  float mZjj;
  int leptType;

  treeMC->SetBranchAddress("mZZ",&mZZ);
  treeMC->SetBranchAddress("eventWeight",&eventWeight);
  treeMC->SetBranchAddress("nBTags",&nBTags);
  treeMC->SetBranchAddress("mZjj",&mZjj);
  treeMC->SetBranchAddress("leptType",&leptType);

  
  float bins0[26]={150,165,180,195,210,225,240,255,270,285,300,320,340,360,380,400,430,460,490,520,550,600,650,700,750,800};
   
  TH1D* h1_mZZ_signalRegion = new TH1D("mZZ_signalRegion", "", 25, bins0);
  h1_mZZ_signalRegion->Sumw2();
  TH1D* h1_mZZ_sidebands = new TH1D("mZZ_sidebands", "", 25, bins0);
  h1_mZZ_sidebands->Sumw2();

  for( unsigned iEntry=0; iEntry<treeMC->GetEntries(); ++iEntry ) {

    treeMC->GetEntry(iEntry);
    if( iEntry%10000 == 0 ) std::cout << "Entry: " << iEntry << "/" << treeMC->GetEntries() << std::endl;

    if( leptType_str=="MU" && leptType!=0 ) continue;
    if( leptType_str=="ELE" && leptType!=1 ) continue;
    if( nBTags!=btagCategory ) continue;
    if( mZZ>800. || mZZ < 183. ) continue;
 
    bool isSignalRegion = (mZjj>75. && mZjj<105.);
    if( isSignalRegion ) h1_mZZ_signalRegion->Fill(mZZ, eventWeight);
    if( !isSignalRegion && mZjj>60. && mZjj<130.) h1_mZZ_sidebands->Fill(mZZ, eventWeight);

  }

  TH1D* h1_alpha = new TH1D(*h1_mZZ_signalRegion);
  h1_alpha->SetName("alpha");
  h1_alpha->Sumw2();
  h1_alpha->Divide(h1_mZZ_sidebands);

  // smooth it:
  double BinContent=0;
  double SmoothingThreshold=3.0;
  for(int iBin=1; iBin<h1_alpha->GetNbinsX(); iBin++){
    if(h1_alpha->GetBinContent(iBin)>SmoothingThreshold){
      //if(iBin!=h1_alpha->GetNbinsX()){
        if( h1_alpha->GetBinContent(iBin+1)<SmoothingThreshold && h1_alpha->GetBinContent(iBin-1)<SmoothingThreshold )
          h1_alpha->SetBinContent(iBin,(h1_alpha->GetBinContent(iBin+1)+h1_alpha->GetBinContent(iBin-1))/2.);     
        else if( h1_alpha->GetBinContent(iBin+1)<SmoothingThreshold )
          h1_alpha->SetBinContent(iBin,h1_alpha->GetBinContent(iBin+1));     
        else if( h1_alpha->GetBinContent(iBin-1)<SmoothingThreshold )
          h1_alpha->SetBinContent(iBin,h1_alpha->GetBinContent(iBin-1));
        else
          h1_alpha->SetBinContent(iBin,1.);
      //} else if(iBin==h1_alpha->GetNbinsX()){
      //  h1_alpha->SetBinContent(iBin,h1_alpha->GetBinContent(iBin-1));
      //} else if(iBin==1){
      //  h1_alpha->SetBinContent(iBin,(h1_alpha->GetBinContent(iBin+1)+h1_alpha->GetBinContent(iBin))/2);
      //}
    }
  }

  return h1_alpha;
  
}




void fitSidebands( const std::string& dataset, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha ) {

  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    return;
  }
  

  std::string outdir = "FitSidebands_" + dataset;

  char ofsMC_name[400];
  sprintf( ofsMC_name, "%s/fitresultsMC_%dbtag.txt", outdir.c_str(), btagCategory);
  ofstream ofsMC(ofsMC_name);


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
  sprintf( alphaFileName, "alphaFile_%s_%dbtag_%s.root", dataset.c_str(), btagCategory, leptType.c_str());
  TFile* file_alpha = TFile::Open(alphaFileName, "recreate");
  file_alpha->cd();
  h1_alpha->Write();
  tree_sidebandsDATA_alpha->Write();
  tree_sidebandsMC_alpha->Write();
  file_alpha->Close();



//float theta_val=0.;
//if( btagCategory==0 ) theta_val = -1.5545;
//if( btagCategory==1 ) theta_val = -1.552;
//if( btagCategory==2 ) theta_val = -1.5575;

  double a0 = -1.395;
  double w0 = 85.73;
  //double a=cos(-theta_val)*a0 - sin(-theta_val)*w0;
  //double w=sin(-theta_val)*a0 + cos(-theta_val)*w0;

  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",191.12,175.,220.);
  RooRealVar cutOff2("cutOff2","position of fermi",191.12,175.,220.);
  //cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",4.698,0.,30.);
  RooRealVar beta2("beta2","width of fermi",4.698,0.,30.);
  //beta.setConstant(kTRUE);
  RooFermi fermi("fermi","fermi function",*mZZ,cutOff2,beta2);
  RooFermi fermi2("fermi2","fermi function",*mZZ,cutOff2,beta2);

  // -------------------- crystal ball ---------------------------
  RooRealVar m("m","m",200.17,190.,300.);
  RooRealVar m2("m2","m2",200.17,190.,300.);
  //m.setConstant(kTRUE);
  RooRealVar wdth("wdth","wdth",w0,-200.,200.);
  RooRealVar wdth0("wdth0","wdth0",w0,-200.,200.);
  RooRealVar n("n","n",13.067,0.,100.);
  RooRealVar n2("n2","n2",13.067,0.,100.);
  RooRealVar alpha("alpha","alpha",a0,-200.,200.); 
  RooRealVar alpha0("alpha0","alpha0",a0,-200.,200.); 

  RooRealVar theta("theta","theta",0.,-3.1416,3.1416); 
  theta.setConstant(kTRUE);


  RooCB CB("CB","Crystal ball",*mZZ,m,wdth,alpha,n, theta);
  RooCBShape CBShape("CB","Crystal ball",*mZZ,m2,wdth0,alpha0,n2);

  RooProdPdf background("background","background",RooArgSet(fermi,CB));
  RooProdPdf background2("background","background",RooArgSet(fermi2,CBShape));
 





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
  RooFitResult *r_sidebandsMC_alpha = background.fitTo(sidebandsMC_alpha,SumW2Error(kTRUE), Save());
  RooFitResult *r_sidebandsMC_alpha_2 = background2.fitTo(sidebandsMC_alpha,SumW2Error(kTRUE), Save());
  //RooFitResult *r_sidebandsMC = exp.fitTo(sidebandsMC,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());

  ofsMC << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
  ofsMC << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
  ofsMC << "m " << m.getVal() << " " << m.getError() << std::endl;
  ofsMC << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
  ofsMC << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
  ofsMC << "n " << n.getVal() << " " << n.getError() << std::endl;
  ofsMC << "theta " << theta.getVal() << " " << theta.getError() << std::endl;

  Double_t rhoMC_sidebands_alpha = r_sidebandsMC_alpha->correlation("alpha", "wdth");
  TMatrixDSym corrMatrixMC_sidebands_alpha = r_sidebandsMC_alpha->correlationMatrix();
  TMatrixDSym covMatrixMC_sidebands_alpha = r_sidebandsMC_alpha->covarianceMatrix();
  
  //ofsMC << std::endl;
  //ofsMC << "Correlation matrix: " << std::endl;
  //ofsMC << corrMatrixMC_sidebands_alpha[0][0] << " " << corrMatrixMC_sidebands_alpha[0][1] << std::endl;
  //ofsMC << corrMatrixMC_sidebands_alpha[1][0] << " " << corrMatrixMC_sidebands_alpha[1][1] << std::endl;

  ofsMC.close();

  

  RooPlot *plot_sidebandsMC_alpha = mZZ->frame();

  //sidebandsMC.plotOn(plot_sidebandsMC, Binning(nBins));
  sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins));

  background.plotOn(plot_sidebandsMC_alpha, LineColor(kRed));
  background2.plotOn(plot_sidebandsMC_alpha, LineColor(38), LineStyle(2));
  //sidebandsMC.plotOn(plot_sidebandsMC_alpha, Binning(nBins));
  sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins));

  plot_sidebandsMC_alpha->Draw();

  char canvasName[400];
  sprintf( canvasName, "%s/mZZ_sidebandsMC_alpha_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
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
  theta.setConstant(kTRUE);


  c1->Clear();
  c1->SetLogy(false);

  //RooFitResult *r_signalMC  = background.fitTo(signalMC);
  //RooFitResult *r_signalMC = exp.fitTo(signalMC,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());

//char ofsMCsig_name[400];
//sprintf( ofsMCsig_name, "FitSidebands/fitresultsMCsig_%dbtag.txt", btagCategory);
//ofstream ofsMCsig(ofsMCsig_name);


//ofsMCsig << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
//ofsMCsig << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
//ofsMCsig << "m " << m.getVal() << " " << m.getError() << std::endl;
//ofsMCsig << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
//ofsMCsig << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
//ofsMCsig << "n " << n.getVal() << " " << n.getError() << std::endl;
//ofsMCsig << "theta " << theta.getVal() << " " << theta.getError() << std::endl;

//ofsMCsig.close();


  RooPlot *plot_signalMC  = mZZ->frame();

  signalMC.plotOn(plot_signalMC , Binning(nBins));

  background.plotOn(plot_signalMC, LineColor(kRed));
  signalMC.plotOn(plot_signalMC, Binning(nBins));

  plot_signalMC->Draw();

  sprintf( canvasName, "%s/mZZ_signalMC_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
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

  RooFitResult *r_sidebandsDATA_alpha = background.fitTo(sidebandsDATA_alpha, SumW2Error(kFALSE), Save());
  //RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());

  char ofsDATA_name[600];
  sprintf( ofsDATA_name, "%s/fitresultsDATA_%dbtag.txt", outdir.c_str(), btagCategory);
  ofstream ofsDATA(ofsDATA_name);


  ofsDATA << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
  ofsDATA << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
  ofsDATA << "m " << m.getVal() << " " << m.getError() << std::endl;
  ofsDATA << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
  ofsDATA << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
  ofsDATA << "n " << n.getVal() << " " << n.getError() << std::endl;
  ofsDATA << "theta " << theta.getVal() << " " << theta.getError() << std::endl;

  Double_t rhoDATA_sidebands_alpha = r_sidebandsDATA_alpha->correlation("alpha", "wdth");
  TMatrixDSym corrMatrixDATA_sidebands_alpha = r_sidebandsDATA_alpha->correlationMatrix();
  TMatrixDSym covMatrixDATA_sidebands_alpha = r_sidebandsDATA_alpha->covarianceMatrix();
  
  //ofsDATA << std::endl;
  //ofsDATA << "Correlation matrix: " << std::endl;
  //ofsDATA << corrMatrixDATA_sidebands_alpha[0][0] << " " << corrMatrixDATA_sidebands_alpha[0][1] << std::endl;
  //ofsDATA << corrMatrixDATA_sidebands_alpha[1][0] << " " << corrMatrixDATA_sidebands_alpha[1][1] << std::endl;


  RooPlot *plot_sidebandsDATA_alpha = mZZ->frame();

  sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

  background.plotOn(plot_sidebandsDATA_alpha, LineColor(kRed));
  sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

  plot_sidebandsDATA_alpha->Draw();

  sprintf( canvasName, "%s/mZZ_sidebandsDATA_alpha_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  Trying to find decorrelation " << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);


  RooPlot *plot_rot = mZZ->frame();

  background.plotOn(plot_rot, LineColor(kRed));

  
  double precision=0.05;
  double lowerBound = -2.;
  double upperBound = 0.;
  double bestValue = r_sidebandsDATA_alpha->correlation("alpha", "wdth");
  double bestTheta = theta.getVal();
  double alpha_fit = alpha.getVal();
  double width_fit = wdth.getVal();

  while(fabs(bestValue) > precision){

    double last = 0.;

    for(int i =0; i < 30; i++){

      theta.setVal(lowerBound+i*(upperBound-lowerBound)/30.);
      double a=cos(-theta.getVal())*alpha_fit - sin(-theta.getVal())*width_fit;
      double w=sin(-theta.getVal())*alpha_fit + cos(-theta.getVal())*width_fit;
      alpha.setVal(a);
      wdth.setVal(w);
      RooFitResult *r_sidebandsDATA_alpha = background.fitTo(sidebandsDATA_alpha, SumW2Error(kFALSE), Save());
      double newCor = r_sidebandsDATA_alpha->correlation("alpha", "wdth");
      if(fabs(newCor)<fabs(bestValue)){
        bestValue=newCor;
        bestTheta=theta.getVal();
      }
      if(newCor * last < 0. ){// found a zero-crossing
        double oldstep = (upperBound-lowerBound)/30.;
        lowerBound = theta.getVal()-5.*oldstep;
        upperBound = theta.getVal()+5.*oldstep;
        break;
      } else{
        last = newCor;
      }

    } //for i 0-30

  } //while precision


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  found best angle "<< bestTheta << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;
  
  double a_rot = cos(-bestTheta)*alpha_fit - sin(-bestTheta)*width_fit;
  double w_rot = sin(-bestTheta)*alpha_fit + cos(-bestTheta)*width_fit;

  RooRealVar wdth_rot("wdth_rot","wdth_rot",w_rot,-200.,200.);
  RooRealVar alpha_rot("alpha_rot","alpha_rot",a_rot,-200.,200.);

  RooRealVar theta_best("theta_best","theta_best",bestTheta,-3.1416,3.1416);


  RooCB CB_rot("CB","Crystal ball",*mZZ,m,wdth_rot,alpha_rot,n, theta_best);

  RooProdPdf background_rot("background_rot","background_rot",RooArgSet(fermi,CB_rot));
  background_rot.plotOn(plot_rot, LineColor(38), LineStyle(2));

  plot_rot->Draw();

  char canvasName_rot[400];
  sprintf( canvasName_rot, "%s/check_rot_%dbtag.eps", outdir.c_str(), btagCategory);
  c1->SaveAs(canvasName_rot);
  

  c1->SetLogy();
  sprintf( canvasName_rot, "%s/check_rot_%dbtag_log.eps", outdir.c_str(), btagCategory);
  c1->SaveAs(canvasName_rot);

  ofsDATA <<  "alpha_rot " << a_rot << std::endl;
  ofsDATA <<  "wdth_rot " << w_rot << std::endl;
  ofsDATA << "theta_best " << bestTheta << std::endl;
  ofsDATA.close();




  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT DATA SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;

  //fix shape:
  cutOff.setConstant(kTRUE);
  beta.setConstant(kTRUE);
  m.setConstant(kTRUE);
  wdth.setConstant(kTRUE);
  n.setConstant(kTRUE);
  alpha.setConstant(kTRUE);

//char ofsDATAsig_name[400];
//sprintf( ofsDATAsig_name, "FitSidebands/fitresultsDATAsig_%dbtag.txt", btagCategory);
//ofstream ofsDATAsig(ofsDATAsig_name);


//ofsDATAsig << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
//ofsDATAsig << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
//ofsDATAsig << "m " << m.getVal() << " " << m.getError() << std::endl;
//ofsDATAsig << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
//ofsDATAsig << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
//ofsDATAsig << "n " << n.getVal() << " " << n.getError() << std::endl;
//ofsDATAsig << "theta " << theta.getVal() << " " << theta.getError() << std::endl;

//ofsDATAsig.close();



  c1->Clear();
  c1->SetLogy(false);

  //RooFitResult *r_signalDATA = background.fitTo(signalDATA);
  //RooFitResult *r_signalDATA = exp.fitTo(signalDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());


  RooPlot *plot_signalDATA = mZZ->frame();

  signalDATA.plotOn(plot_signalDATA, Binning(nBins));

  background.plotOn(plot_signalDATA, LineColor(kRed));
  signalDATA.plotOn(plot_signalDATA, Binning(nBins));

  plot_signalDATA->Draw();

  sprintf( canvasName, "%s/mZZ_signalDATA_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
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
  delete plot_signalMC;
 // delete plot_signalDATA;
 // delete plot_sidebandsDATA;
  delete plot_sidebandsDATA_alpha;
  delete plot_signalDATA;
 // delete r_sidebandsMC;
 // delete r_signalMC;
  delete r_sidebandsMC_alpha;
  //delete r_signalMC;
//  delete r_sidebandsDATA;
  delete r_sidebandsDATA_alpha;
  //delete r_signalDATA_alpha;
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


