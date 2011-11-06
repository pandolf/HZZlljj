#include "SidebandFitter.h"

#include <cstdlib>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
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
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"

#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

#include "PdfDiagonalizer.h"



using namespace RooFit;


SidebandFitter::SidebandFitter( const std::string& dataset, const std::string PUType ) {

  dataset_ = dataset;
  PUType_ = PUType;

  mZZmin_ = 160.;
  mZZmax_ = 800.;

}



TH1D* SidebandFitter::getAlphaHisto( int btagCategory, const std::string leptType_str, TTree* treeMC ) {


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
    if( mZZ>mZZmax_ || mZZ < mZZmin_ ) continue;
 
    bool isSignalRegion = (mZjj>75. && mZjj<105.);
    if( isSignalRegion ) h1_mZZ_signalRegion->Fill(mZZ, eventWeight);
    if( !isSignalRegion && mZjj>60. && mZjj<130.) h1_mZZ_sidebands->Fill(mZZ, eventWeight);

  }

  TH1D* h1_alpha = new TH1D(*h1_mZZ_signalRegion);
  h1_alpha->SetName("alpha");
  h1_alpha->Sumw2();
  h1_alpha->Divide(h1_mZZ_sidebands);

  // smooth it:
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



RooFitResult* SidebandFitter::fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed ) {

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
  

  std::string outdir = get_outdir();
  std::string mkdir_command = "mkdir -p " + outdir;
  system(mkdir_command.c_str());


  char cut_base[500];
  sprintf( cut_base, "nBTags==%d %s", btagCategory, leptType_cut.c_str());
  char cut_sidebands[500];
  sprintf( cut_sidebands, "%s && ( (mZjj>60. && mZjj<75.)||(mZjj>105. && mZjj<130.) )", cut_base);
  char cut_signal[500];
  sprintf( cut_signal, "%s && ( mZjj>75. && mZjj<105. )", cut_base);
  

  float binWidth = 20.;
  int nBins = (int)(mZZmax_-mZZmin_)/binWidth;
  RooRealVar* CMS_hzz2l2q_mZZ = new RooRealVar("CMS_hzz2l2q_mZZ", "m_{lljj}", mZZmin_, mZZmax_, "GeV");

  RooRealVar* eventWeight = new RooRealVar("eventWeight", "event weight", 0., 2., "");
  RooRealVar* eventWeight_alpha = new RooRealVar("eventWeight_alpha", "event weight (alpha corrected)", 0., 2., "");
  RooRealVar* nBTags = new RooRealVar("nBTags", "number of BTags", -1., 2., "");
  RooRealVar* mZjj = new RooRealVar("mZjj", "mZjj", 60., 130., "GeV");


  RooDataSet sidebandsMC("sidebandsMC","sidebandsMC",treeMC,RooArgSet(*eventWeight,*CMS_hzz2l2q_mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight");
  RooDataSet signalMC("signalMC","signalMC",treeMC,RooArgSet(*eventWeight,*CMS_hzz2l2q_mZZ,*nBTags,*mZjj),cut_signal,"eventWeight");

  char suffix[20];
  if( !writeFile )
    sprintf( suffix, "_%d", seed );
  else
    sprintf( suffix, "" );
  char treeName_MC[200];
  sprintf( treeName_MC, "sidebandsMC_alpha%s", suffix );
  std::string treeName_MC_str(treeName_MC);
  std::cout << "Correcting signal (MC): " << std::endl;
  TTree* tree_sidebandsMC_alpha = correctTreeWithAlpha( treeMC, h1_alpha, btagCategory, treeName_MC_str );
  tree_sidebandsMC_alpha->GetBranch("mZZ")->SetName("CMS_hzz2l2q_mZZ"); 
  RooDataSet sidebandsMC_alpha("sidebandsMC_alpha","sidebandsMC_alpha",tree_sidebandsMC_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*CMS_hzz2l2q_mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");

  RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",treeDATA,RooArgSet(*eventWeight,*CMS_hzz2l2q_mZZ,*nBTags,*mZjj),cut_sidebands);
  RooDataSet signalDATA("signalDATA","signalDATA",treeDATA,RooArgSet(*eventWeight,*CMS_hzz2l2q_mZZ,*nBTags,*mZjj),cut_signal);
  char treeName_DATA[200];
  sprintf( treeName_DATA, "sidebandsDATA_alpha%s", suffix );
  std::string treeName_DATA_str(treeName_DATA);
  std::cout << "Correcting signal (DATA): " << std::endl;
  TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( treeDATA, h1_alpha, btagCategory, treeName_DATA_str );
  tree_sidebandsDATA_alpha->GetBranch("mZZ")->SetName("CMS_hzz2l2q_mZZ"); 
  RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*CMS_hzz2l2q_mZZ,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");


  TFile* file_alpha = 0;


  double a0 = -1.395;
  double w0 = 85.73;

  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",191.12,175.,220.);
  RooRealVar beta("beta","width of fermi",4.698,0.,30.);
  RooFermi fermi("fermi","fermi function",*CMS_hzz2l2q_mZZ,cutOff,beta);

  // -------------------- crystal ball ---------------------------
  RooRealVar m("m","m",200.17,190.,300.);
  RooRealVar wdth("wdth","wdth",w0,-200.,200.);
  RooRealVar n("n","n",13.067,0.,100.);
  RooRealVar alpha("alpha","alpha",a0,-200.,200.); 


  RooCBShape CBShape("CB","Crystal ball",*CMS_hzz2l2q_mZZ,m,wdth,alpha,n);

  RooProdPdf background("background","background",RooArgSet(fermi,CBShape));
 





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

  if( writeFile ) {

    std::string ofsMCName = get_fitResultsName( btagCategory, "MC" );
    ofstream ofsMC(ofsMCName.c_str());

    ofsMC << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
    ofsMC << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
    ofsMC << "m " << m.getVal() << " " << m.getError() << std::endl;
    ofsMC << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
    ofsMC << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
    ofsMC << "n " << n.getVal() << " " << n.getError() << std::endl;

    ofsMC.close();

    RooPlot *plot_sidebandsMC_alpha = CMS_hzz2l2q_mZZ->frame(mZZmin_, mZZmax_, nBins);

    sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins));

    background.plotOn(plot_sidebandsMC_alpha, LineColor(kRed));
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


    c1->Clear();
    c1->SetLogy(false);


    RooPlot *plot_signalMC  = CMS_hzz2l2q_mZZ->frame(mZZmin_, mZZmax_, nBins);

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
  char fitResultName[200];
  if( leptType!="ALL" )
    sprintf( fitResultName, "fitResults_%dbtag_%s", btagCategory, leptType.c_str() );
  else 
    sprintf( fitResultName, "fitResults_%dbtag", btagCategory );
  r_sidebandsDATA_alpha->SetName(fitResultName);


  RooPlot *plot_sidebandsDATA_alpha = CMS_hzz2l2q_mZZ->frame(mZZmin_, mZZmax_, nBins);


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


  char fitWorkspaceName[200];
  sprintf( fitWorkspaceName, "fitWorkspace_%dbtag", btagCategory );
  RooWorkspace* fitWorkspace = new RooWorkspace(fitWorkspaceName, fitWorkspaceName);
  //fitWorkspace->addClassDeclImportDir("/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms3/include/");

  //fitWorkspace->importClassCode(RooFermi::Class(),kTRUE);
  //fitWorkspace->importClassCode("RooFermi",kTRUE);



  //now decorrelate parameters:
  char diagonalizerName[200];
  sprintf( diagonalizerName, "CMS_hzz2l2q_bkg_%db", btagCategory);
  PdfDiagonalizer diago(diagonalizerName, fitWorkspace, *r_sidebandsDATA_alpha);
  RooAbsPdf *background_decorr = diago.diagonalize(background);
  background_decorr->SetName("background_decorr");
  RooFitResult *r_sidebandsDATA_alpha_decorr = background_decorr->fitTo(sidebandsDATA_alpha, SumW2Error(kFALSE), Save(), Warnings(warnings), PrintLevel(warningLevel));
  char fitResultName_eig[200];
  if( leptType!="ALL" )
    sprintf( fitResultName_eig, "%s_%s_decorr", fitResultName, leptType.c_str() );
  else 
    sprintf( fitResultName_eig, "%s_decorr", fitResultName );
  r_sidebandsDATA_alpha_decorr->SetName(fitResultName_eig);

  //import both pdfs in the workspace:
  fitWorkspace->import(*CMS_hzz2l2q_mZZ);
  fitWorkspace->import(background);
  fitWorkspace->import(*background_decorr, RooFit::RecycleConflictNodes());


  RooPlot *plot_rot = CMS_hzz2l2q_mZZ->frame(mZZmin_, mZZmax_, nBins);
  sidebandsDATA_alpha.plotOn(plot_rot);
  background.plotOn(plot_rot,RooFit::LineColor(kRed));
  background_decorr->plotOn(plot_rot,RooFit::LineColor(kBlue), RooFit::LineStyle(2));
  TCanvas* c_rot = new TCanvas("rot", "", 600, 600);
  c_rot->cd();
  plot_rot->Draw();
  char canvasName_rot[200];
  sprintf( canvasName_rot, "checkrot_%dbtag.eps", btagCategory );
  c_rot->SaveAs(canvasName_rot);





  if( writeFile ) {

  //RooRealVar wdth_rot("wdth_rot","wdth_rot",w_rot,-200.,200.);
  //RooRealVar alpha_rot("alpha_rot","alpha_rot",a_rot,-200.,200.);

  //RooRealVar theta_best("theta_best","theta_best",bestTheta,-3.1416,3.1416);


  //RooCB CB_rot("CB","Crystal ball",*mZZ,m,wdth_rot,alpha_rot,n, theta_best);

  //RooProdPdf background_rot("background_rot","background_rot",RooArgSet(fermi,CB_rot));
  //background_rot.plotOn(plot_rot, LineColor(38), LineStyle(2));

  //plot_rot->Draw();

  //char canvasName_rot[400];
  //sprintf( canvasName_rot, "%s/check_rot_%dbtag.eps", outdir.c_str(), btagCategory);
  //c1->SaveAs(canvasName_rot);
  //

  //c1->SetLogy();
  //sprintf( canvasName_rot, "%s/check_rot_%dbtag_log.eps", outdir.c_str(), btagCategory);
  //c1->SaveAs(canvasName_rot);

    std::string ofsDATAName = get_fitResultsName( btagCategory, "DATA" );
    ofstream ofsDATA(ofsDATAName.c_str());


    ofsDATA << "beta " << beta.getVal() << " " << beta.getError() << std::endl;
    ofsDATA << "cutOff " << cutOff.getVal() << " " << cutOff.getError() << std::endl;
    ofsDATA << "m " << m.getVal() << " " << m.getError() << std::endl;
    ofsDATA << "wdth " << wdth.getVal() << " " << wdth.getError() << std::endl;
    ofsDATA << "alpha " << alpha.getVal() << " " << alpha.getError() << std::endl;
    ofsDATA << "n " << n.getVal() << " " << n.getError() << std::endl;

  //ofsDATA << "alpha_rot " << a_rot << " 0" << std::endl;
  //ofsDATA << "wdth_rot " << w_rot << " 0" << std::endl;
  //ofsDATA << "theta_best " << bestTheta << " 0" << std::endl;
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


    RooPlot *plot_signalDATA = CMS_hzz2l2q_mZZ->frame(mZZmin_, mZZmax_, nBins);

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




  if( writeFile ) {

    char fitResultsFileName[500];
    sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_%s_PU%s.root", dataset_.c_str(), btagCategory, leptType.c_str(), PUType_.c_str());
    file_alpha = TFile::Open(fitResultsFileName, "recreate");
    file_alpha->cd();
    h1_alpha->Write();
    tree_sidebandsDATA_alpha->Write();
    tree_sidebandsMC_alpha->Write();
    r_sidebandsDATA_alpha->Write();
    r_sidebandsDATA_alpha_decorr->Write();
    fitWorkspace->Write();
    file_alpha->Close();

  }




  delete eventWeight;
  delete eventWeight_alpha;
  delete CMS_hzz2l2q_mZZ;
  delete nBTags;
  delete mZjj;
  delete c1;
  delete r_sidebandsMC_alpha;


  return r_sidebandsDATA_alpha;

}



std::string SidebandFitter::get_fitResultsName( int nbtags, const std::string& data_mc ) {

  std::string outdir = get_outdir();

  char fitResultsName[600];
  sprintf( fitResultsName, "%s/fitresults%s_%dbtag.txt", outdir.c_str(), data_mc.c_str(), nbtags);
  std::string returnString(fitResultsName);

  return returnString;

}



std::string SidebandFitter::get_outdir() {

  std::string returnString = "FitSidebands_" + dataset_;

  return returnString;

}



TTree* SidebandFitter::correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name ) {

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
    if( isSidebands && mZZ>mZZmin_ && mZZ<mZZmax_ ) newWeight *= alpha;

    newTree->Fill();

  }

  return newTree;

}



TH1D* SidebandFitter::shuffle( TH1D* inhist, TRandom3* random, char *histName ) {

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


