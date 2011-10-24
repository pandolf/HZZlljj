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
#include "TRandom3.h"

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



struct FitResults {

  float fermi_beta;
  float fermi_cutoff;
  float CB_m;
  float CB_wdth;
  float CB_alpha;
  float CB_n;
  float CB_theta;

  float fermi_beta_err;
  float fermi_cutoff_err;
  float CB_m_err;
  float CB_wdth_err;
  float CB_alpha_err;
  float CB_n_err;
  float CB_theta_err;

  float CB_alpha_rot;
  float CB_wdth_rot;
  float CB_theta_best;

};



TH1D* getAlphaHisto( int nBTags, const std::string leptType_str, TTree* treeMC );
//TH1D* getAlphaHisto( int nBTags, const std::string leptType, TFile* file_ZJets_madgraph, TFile* file_TT_TW, TFile* file_Diboson );
//FitResults fitSidebands( const std::string& dataset, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha=0, bool writeFile=true );
FitResults fitSidebands( const std::string& dataset, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha=0, int seed=-1 );
TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name );
//TH1D* shuffle(TH1D* inhist, TRandom3 random);
TH1D* shuffle( TH1D* inhist, TRandom3* random, char *histname );


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

  TH1D* alpha_0btag = getAlphaHisto( 0, "ALL", treeMC_0btag );
  TH1D* alpha_1btag = getAlphaHisto( 1, "ALL", treeMC_1btag );
  TH1D* alpha_2btag = getAlphaHisto( 2, "ALL", treeMC_2btag );

  FitResults fr_0btag = fitSidebands( dataset, treeMC_0btag, treeDATA_0btag, 0, "ALL", alpha_0btag );
  FitResults fr_1btag = fitSidebands( dataset, treeMC_1btag, treeDATA_1btag, 1, "ALL", alpha_1btag );
  FitResults fr_2btag = fitSidebands( dataset, treeMC_2btag, treeDATA_2btag, 2, "ALL", alpha_2btag );


  // and now generate toys to get the error on the relevant fit parameters:
  TH1D* h1_alphaFit_0btag = new TH1D("alphaFit_0btag", "", 100, -200., 50.); 
  TH1D* h1_wdthFit_0btag = new TH1D("wdthFit_0btag", "", 100, 0., 1.); 

  TH1D* h1_alphaFit_1btag = new TH1D("alphaFit_1btag", "", 100, -200., 50.); 
  TH1D* h1_wdthFit_1btag = new TH1D("wdthFit_1btag", "", 100, 0., 1.); 

  TH1D* h1_alphaFit_2btag = new TH1D("alphaFit_2btag", "", 100, -200., 50.); 
  TH1D* h1_wdthFit_2btag = new TH1D("wdthFit_2btag", "", 100, 0., 1.); 

  TRandom3* rand = new TRandom3(13);
  int nToys=2;
  TFile* prova = TFile::Open("PROVA.root", "recreate");
  prova->cd();

  for( unsigned iToy=0; iToy<nToys; ++iToy ) {

    std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;;
    std::cout << " +++ TOY:   " << iToy+1 << " / " << nToys << std::endl;
    std::cout << std::endl << std::endl;


    char histName[100];
    sprintf( histName, "alphaHist_0btag_%d", iToy );
    TH1D* randomAlpha_0btag = shuffle( alpha_0btag, rand, histName );
    sprintf( histName, "alphaHist_1btag_%d", iToy );
    TH1D* randomAlpha_1btag = shuffle( alpha_1btag, rand, histName );
    sprintf( histName, "alphaHist_2btag_%d", iToy );
    TH1D* randomAlpha_2btag = shuffle( alpha_2btag, rand, histName );

    randomAlpha_0btag->Write();
    randomAlpha_1btag->Write();
    randomAlpha_2btag->Write();

    FitResults fr_0btag_iToy = fitSidebands( dataset, treeMC_0btag, treeDATA_0btag, 0, "ALL", randomAlpha_0btag, iToy );
    FitResults fr_1btag_iToy = fitSidebands( dataset, treeMC_1btag, treeDATA_1btag, 1, "ALL", randomAlpha_1btag, iToy );
    FitResults fr_2btag_iToy = fitSidebands( dataset, treeMC_2btag, treeDATA_2btag, 2, "ALL", randomAlpha_2btag, iToy );

    if( fr_0btag_iToy.CB_theta!=0. ) {
      h1_alphaFit_0btag->Fill( fr_0btag_iToy.CB_alpha_rot );
      h1_wdthFit_0btag->Fill( fr_0btag_iToy.CB_wdth_rot );
    }

    if( fr_1btag_iToy.CB_theta!=0. ) {
      h1_alphaFit_1btag->Fill( fr_1btag_iToy.CB_alpha_rot );
      h1_wdthFit_1btag->Fill( fr_1btag_iToy.CB_wdth_rot );
    }

    if( fr_2btag_iToy.CB_theta!=0. ) {
      h1_alphaFit_2btag->Fill( fr_2btag_iToy.CB_alpha_rot );
      h1_wdthFit_2btag->Fill( fr_2btag_iToy.CB_wdth_rot );
    }

  }

  h1_alphaFit_0btag->Write();
  h1_wdthFit_0btag->Write();
  prova->Close();

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
  for(int iBin=1; iBin<h1_alpha->GetNbinsX()+1; iBin++) {
    if(h1_alpha->GetBinContent(iBin)>SmoothingThreshold) {
        if(iBin!=h1_alpha->GetNbinsX()) {
          if( h1_alpha->GetBinContent(iBin+1)<SmoothingThreshold && h1_alpha->GetBinContent(iBin-1)<SmoothingThreshold )
            h1_alpha->SetBinContent(iBin,(h1_alpha->GetBinContent(iBin+1)+h1_alpha->GetBinContent(iBin-1))/2.);     
          else if( h1_alpha->GetBinContent(iBin+1)<SmoothingThreshold )
            h1_alpha->SetBinContent(iBin,h1_alpha->GetBinContent(iBin+1));     
          else if( h1_alpha->GetBinContent(iBin-1)<SmoothingThreshold )
            h1_alpha->SetBinContent(iBin,h1_alpha->GetBinContent(iBin-1));
          else
            h1_alpha->SetBinContent(iBin,1.);
        } else if(iBin==h1_alpha->GetNbinsX()){
          h1_alpha->SetBinContent(iBin,h1_alpha->GetBinContent(iBin-1));
        } else if(iBin==1){
          h1_alpha->SetBinContent(iBin,h1_alpha->GetBinContent(iBin+1));
        }
     } //if over thresh
  } //for bins

  return h1_alpha;
  
}




//FitResults fitSidebands( const std::string& dataset, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, bool writeFile ) {
FitResults fitSidebands( const std::string& dataset, TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed ) {

  bool writeFile = (seed==-1);

  bool warnings;
  int warningLevel;
  if( writeFile ) {
    warnings = true;
    warningLevel = 1;
  } else {
    warnings = false;
    warningLevel = -1;
  }


  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    exit(109);
  }
  

  std::string outdir = "FitSidebands_" + dataset;
  std::string mkdir_command = "mkdir -p " + outdir;
  system(mkdir_command.c_str());


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

  char suffix[20];
  if( !writeFile )
    sprintf( suffix, "_%d", seed );
  char treeName_MC[200];
  sprintf( treeName_MC, "sidebandsMC_alpha%s", suffix );
  std::string treeName_MC_str(treeName_MC);
  std::cout << "Correcting signal (MC): " << std::endl;
  TTree* tree_sidebandsMC_alpha = correctTreeWithAlpha( treeMC, h1_alpha, btagCategory, treeName_MC_str );
  RooDataSet sidebandsMC_alpha("sidebandsMC_alpha","sidebandsMC_alpha",tree_sidebandsMC_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");

  RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",treeDATA,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_sidebands);
  RooDataSet signalDATA("signalDATA","signalDATA",treeDATA,RooArgSet(*eventWeight,*mZZ,*nBTags,*mZjj),cut_signal);
  char treeName_DATA[200];
  sprintf( treeName_DATA, "sidebandsDATA_alpha%s", suffix );
  std::string treeName_DATA_str(treeName_DATA);
  std::cout << "Correcting signal (DATA): " << std::endl;
  TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( treeDATA, h1_alpha, btagCategory, treeName_DATA_str );
  RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");


  TFile* file_alpha = 0;

  if( writeFile ) {

    char alphaFileName[500];
    sprintf( alphaFileName, "alphaFile_%s_%dbtag_%s.root", dataset.c_str(), btagCategory, leptType.c_str());
    file_alpha = TFile::Open(alphaFileName, "recreate");
    file_alpha->cd();
    h1_alpha->Write();
    tree_sidebandsDATA_alpha->Write();
    tree_sidebandsMC_alpha->Write();
    file_alpha->Close();

  }



  double a0 = -1.395;
  double w0 = 85.73;

  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",191.12,175.,220.);
  RooRealVar cutOff2("cutOff2","position of fermi",191.12,175.,220.);
  RooRealVar beta("beta","width of fermi",4.698,0.,30.);
  RooRealVar beta2("beta2","width of fermi",4.698,0.,30.);
  RooFermi fermi("fermi","fermi function",*mZZ,cutOff,beta);
  RooFermi fermi2("fermi2","fermi function",*mZZ,cutOff2,beta2);

  // -------------------- crystal ball ---------------------------
  RooRealVar m("m","m",200.17,190.,300.);
  RooRealVar m2("m2","m2",200.17,190.,300.);
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

  RooFitResult *r_sidebandsMC_alpha = background.fitTo(sidebandsMC_alpha,SumW2Error(kTRUE), Save(), Warnings(warnings), PrintLevel(warningLevel));
  RooFitResult *r_sidebandsMC_alpha_2 = background2.fitTo(sidebandsMC_alpha,SumW2Error(kTRUE), Save(), Warnings(warnings), PrintLevel(warningLevel));

  if( writeFile ) {

    char ofsMC_name[400];
    sprintf( ofsMC_name, "%s/fitresultsMC_%dbtag.txt", outdir.c_str(), btagCategory);
    ofstream ofsMC(ofsMC_name);

    ofsMC << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
    ofsMC << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
    ofsMC << "m " << m.getVal() << " " << m.getError() << std::endl;
    ofsMC << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
    ofsMC << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
    ofsMC << "n " << n.getVal() << " " << n.getError() << std::endl;
    ofsMC << "theta " << theta.getVal() << " " << theta.getError() << std::endl;

    ofsMC.close();

    RooPlot *plot_sidebandsMC_alpha = mZZ->frame(mZZ_min, mZZ_max, nBins);

    sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins));

    background.plotOn(plot_sidebandsMC_alpha, LineColor(kRed));
    background2.plotOn(plot_sidebandsMC_alpha, LineColor(38), LineStyle(2));
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

    delete plot_sidebandsMC_alpha;

  } //if writeFile



  if( writeFile ) {

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


    RooPlot *plot_signalMC  = mZZ->frame(mZZ_min, mZZ_max, nBins);

    background.plotOn(plot_signalMC, LineColor(kRed), Normalization(sidebandsMC_alpha.sumEntries()));
    signalMC.plotOn(plot_signalMC, Binning(nBins));

    plot_signalMC->Draw();

    char canvasName[400];
    sprintf( canvasName, "%s/mZZ_signalMC_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
    std::string* canvasName_str = new std::string(canvasName);
    std::string canvasName_eps = *canvasName_str + ".eps";
    c1->SaveAs(canvasName_eps.c_str());

    c1->SetLogy();
    *canvasName_str += "_log";
    canvasName_eps = *canvasName_str + ".eps";
    c1->SaveAs(canvasName_eps.c_str());

    delete plot_signalMC;

  } //if writeFile



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

  RooFitResult *r_sidebandsDATA_alpha = background.fitTo(sidebandsDATA_alpha, SumW2Error(kFALSE), Save(), Warnings(warnings), PrintLevel(warningLevel));


  RooPlot *plot_sidebandsDATA_alpha = mZZ->frame(mZZ_min, mZZ_max, nBins);


  if( writeFile ) {

    sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

    background.plotOn(plot_sidebandsDATA_alpha, LineColor(kRed));
    sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

    plot_sidebandsDATA_alpha->Draw();

    char canvasName[400];
    sprintf( canvasName, "%s/mZZ_sidebandsDATA_alpha_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
    std::string* canvasName_str = new std::string(canvasName);
    std::string canvasName_eps = *canvasName_str + ".eps";
    c1->SaveAs(canvasName_eps.c_str());

    c1->SetLogy();
    *canvasName_str += "_log";
    canvasName_eps = *canvasName_str + ".eps";
    c1->SaveAs(canvasName_eps.c_str());

  }


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  Trying to find decorrelation (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);


  RooPlot *plot_rot = mZZ->frame(mZZ_min, mZZ_max, nBins);

  background.plotOn(plot_rot, LineColor(kRed));

  
  double precision=0.05;
  double lowerBound = -2.;
  double upperBound = 0.;
  double bestValue = r_sidebandsDATA_alpha->correlation("alpha", "wdth");
  double bestTheta = theta.getVal();
  double alpha_fit = alpha.getVal();
  double width_fit = wdth.getVal();
  int iTry = 0;

  while( fabs(bestValue) > precision && (iTry<=100 || writeFile) ) { //(no more than 100 tries if not writing)

    double last = 0.;

    for(int i =0; i < 30; i++){

      theta.setVal(lowerBound+i*(upperBound-lowerBound)/30.);
      double a=cos(-theta.getVal())*alpha_fit - sin(-theta.getVal())*width_fit;
      double w=sin(-theta.getVal())*alpha_fit + cos(-theta.getVal())*width_fit;
      alpha.setVal(a);
      wdth.setVal(w);
      // now fit but please shut up:
      RooFitResult *r_sidebandsDATA_alpha_rot = background.fitTo(sidebandsDATA_alpha, SumW2Error(kFALSE), Save(), Verbose(false), Warnings(false), PrintEvalErrors(false), PrintLevel(-1));
      double newCor = r_sidebandsDATA_alpha_rot->correlation("alpha", "wdth");
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

      delete r_sidebandsDATA_alpha_rot;

    } //for i 0-30

    iTry++;

  } //while precision


  if( iTry==200 && !writeFile ) bestTheta=0.;

  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  found best angle " << bestTheta << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;
  
  double a_rot = cos(-bestTheta)*alpha_fit - sin(-bestTheta)*width_fit;
  double w_rot = sin(-bestTheta)*alpha_fit + cos(-bestTheta)*width_fit;



  if( writeFile ) {

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

    ofsDATA << "alpha_rot " << a_rot << " 0" << std::endl;
    ofsDATA << "wdth_rot " << w_rot << " 0" << std::endl;
    ofsDATA << "theta_best " << bestTheta << " 0" << std::endl;
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



    c1->Clear();
    c1->SetLogy(false);


    RooPlot *plot_signalDATA = mZZ->frame(mZZ_min, mZZ_max, nBins);

    background.plotOn(plot_signalDATA, LineColor(kRed), Normalization(sidebandsDATA_alpha.sumEntries()));
    signalDATA.plotOn(plot_signalDATA, Binning(nBins));

    plot_signalDATA->Draw();

    char canvasName[400];
    sprintf( canvasName, "%s/mZZ_signalDATA_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
    std::string* canvasName_str = new std::string(canvasName);
    std::string canvasName_eps = *canvasName_str + ".eps";
    c1->SaveAs(canvasName_eps.c_str());

    c1->SetLogy();
    *canvasName_str += "_log";
    canvasName_eps = *canvasName_str + ".eps";
    c1->SaveAs(canvasName_eps.c_str());
  
    delete plot_signalDATA;

  } //if writeFile



  FitResults fr;

  fr.fermi_beta   = beta.getVal();
  fr.fermi_cutoff = cutOff.getVal();
  fr.CB_m         = m.getVal();
  fr.CB_wdth      = wdth.getVal();
  fr.CB_alpha     = alpha.getVal();
  fr.CB_n         = n.getVal();
  fr.CB_theta     = theta.getVal();

  fr.fermi_beta_err   = beta.getError();
  fr.fermi_cutoff_err = cutOff.getError();
  fr.CB_m_err         = m.getError();
  fr.CB_wdth_err      = wdth.getError();
  fr.CB_alpha_err     = alpha.getError();
  fr.CB_n_err         = n.getError();
  fr.CB_theta_err     = theta.getError();

  fr.CB_alpha_rot  = a_rot;
  fr.CB_wdth_rot   = w_rot;
  fr.CB_theta_best = bestTheta;



  delete eventWeight;
  delete eventWeight_alpha;
  delete mZZ;
  delete nBTags;
  delete mZjj;
  delete c1;
  delete r_sidebandsMC_alpha;
  delete r_sidebandsMC_alpha_2;
  delete r_sidebandsDATA_alpha;

  return fr; 

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



TH1D* shuffle( TH1D* inhist, TRandom3* random, char *histName ) {

  TH1D* outhist = (TH1D*) inhist->Clone();
  outhist->SetName(histName);

  for(int i=1 ; i < outhist->GetNbinsX() ; i++) {

    float val = outhist->GetBinContent(i);
    float err = outhist->GetBinError(i);
    if(val==0. || err==0.)
      continue;

    outhist->SetBinContent(i,random->Gaus(val,err));

  }

  return outhist;

}
