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
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooProdPdf.h"
#include "RooEllipse.h"

#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

#include "fitTools.h"

#include "PdfDiagonalizer.h"



using namespace RooFit;


SidebandFitter::SidebandFitter( const std::string& dataset, const std::string PUType ) {

  dataset_ = dataset;
  PUType_ = PUType;

  mZZmin_ = 160.;
  mZZmax_ = 800.;

  workspace_ = new RooWorkspace("SideBandWS","Workspace for Sideband Fits");
  workspace_->addClassDeclImportDir("/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms3/include/");
  workspace_->importClassCode(RooFermi::Class(),kTRUE);
  workspace_->importClassCode("RooFermi",kTRUE);

  //add tree variables;
  char range[200];
  sprintf(range,"CMS_hzz2l2q_mZZ[%f,%f]",mZZmin_, mZZmax_);
  workspace_->factory(range);
  float binWidth = 20.;
  int nBins = (int)(mZZmax_-mZZmin_)/binWidth;
  workspace_->var("CMS_hzz2l2q_mZZ")->setBins(nBins);
  workspace_->factory("eventWeight[0.,2.]");
  workspace_->factory("eventWeight_alpha[0.,2.]");
  workspace_->factory("nBTags[-1., 2.]");
  workspace_->factory("mZjj[60., 130.]");
  workspace_->defineSet("tree", "CMS_hzz2l2q_mZZ,eventWeight,nBTags,mZjj");
  workspace_->defineSet("treeAlpha", "CMS_hzz2l2q_mZZ,eventWeight,nBTags,mZjj,eventWeight_alpha");
  
  // ------------------------ fermi ------------------------------
  workspace_->factory("cutOff[191.12,175.,220]");
  workspace_->factory("beta[4.698,0.,30.]");
  RooFermi fermi("fermi","fermi function",*(workspace_->var("CMS_hzz2l2q_mZZ")),*(workspace_->var("cutOff")),*(workspace_->var("beta")));
  workspace_->import(fermi);

  // -------------------- crystal ball ---------------------------
  workspace_->factory("m[200.17,190.,300.]");
  workspace_->factory("wdth[85.7,-200.,200.]");
  workspace_->factory("n[13.067,0.,100.]");
  workspace_->factory("alpha[-1.395,-200.,200.]"); 
  
  workspace_->factory("RooCBShape::CB(CMS_hzz2l2q_mZZ,m,wdth,alpha,n)");
  workspace_->factory("PROD::background({fermi,CB})");

  workspace_->defineSet("vars", "cutOff,beta,m,wdth,n,alpha");

  
  std::string outdir = get_outdir();
  std::string mkdir_command = "mkdir -p " + outdir;
  system(mkdir_command.c_str());
  
  std::string ofsMCName = get_fitResultsName( 0 , "FreeInitial" );
  workspace_->set("vars")->writeToFile(ofsMCName.c_str());

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

  treeMC->SetBranchAddress("CMS_hzz2l2q_mZZ",&mZZ);
  treeMC->SetBranchAddress("eventWeight",&eventWeight);
  treeMC->SetBranchAddress("nBTags",&nBTags);
  treeMC->SetBranchAddress("mZjj",&mZjj);
  treeMC->SetBranchAddress("leptType",&leptType);

  

  int nbins=0;  Double_t* bins;
  if(btagCategory==0){
    nbins=25;
  }
  if(btagCategory==1){
    nbins=25;
  }
  if(btagCategory==2){
    nbins=15;
  }
  bins = new Double_t[nbins+1];
  
  fitTools::getBins(nbins+1 ,bins, 150., 800., true);
   
  TH1D* h1_mZZ_signalRegion = new TH1D("mZZ_signalRegion", "", nbins, bins);
  h1_mZZ_signalRegion->Sumw2();
  TH1D* h1_mZZ_sidebands = new TH1D("mZZ_sidebands", "", nbins, bins);
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





// generate plots and text-files for MC signal and extrapolated sideband
// the textfiles can be used to fix the parameters for the data-fit
void SidebandFitter::generateFixedPars(TTree* treeMC,int btagCategory, const std::string& leptType, TH1D* h1_alpha){

  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    exit(109);
  }

  char cut_base[500];
  sprintf( cut_base, "nBTags==%d %s", btagCategory, leptType_cut.c_str());
  char cut_sidebands[500];
  sprintf( cut_sidebands, "%s && ( (mZjj>60. && mZjj<75.)||(mZjj>105. && mZjj<130.) )", cut_base);
  char cut_signal[500];
  sprintf( cut_signal, "%s && ( mZjj>75. && mZjj<105. )", cut_base);
  
  
  float binWidth = 20.;
  int nBins = (int)(mZZmax_-mZZmin_)/binWidth;
  
  RooDataSet signalMC("signalMC","signalMC",treeMC,*(workspace_->set("tree")),cut_signal,"eventWeight");

  char treeName_MC[200];
  sprintf( treeName_MC, "sidebandsMC_alpha");
  std::string treeName_MC_str(treeName_MC);
  std::cout << "Correcting sidebands (MC): " << std::endl;
  TTree* tree_sidebandsMC_alpha = correctTreeWithAlpha( treeMC, h1_alpha, btagCategory, treeName_MC_str );
  RooDataSet sidebandsMC_alpha("sidebandsMC_alpha","sidebandsMC_alpha",tree_sidebandsMC_alpha,*(workspace_->set("treeAlpha")),cut_sidebands,"eventWeight_alpha");

  
  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT ALPHA-CORRECTED MC SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);

  std::string ofsMCName = get_fitResultsName( 0 , "FreeInitial" );
  workspace_->argSet("cutOff,beta,m,wdth,n,alpha").readFromFile(ofsMCName.c_str());

  workspace_->Print("v");
  workspace_->pdf("background")->fitTo(sidebandsMC_alpha,SumW2Error(kTRUE));
  
  
  
  ofsMCName = get_fitResultsName( btagCategory, "MC" );
  workspace_->var("cutOff")->setConstant(kTRUE);
  workspace_->var("beta")->setConstant(kTRUE);
  workspace_->var("m")->setConstant(kTRUE);
  workspace_->var("n")->setConstant(kTRUE);
  workspace_->set("vars")->writeToFile(ofsMCName.c_str());
  workspace_->var("cutOff")->setConstant(kFALSE);
  workspace_->var("beta")->setConstant(kFALSE);
  workspace_->var("m")->setConstant(kFALSE);
  workspace_->var("n")->setConstant(kFALSE);

  
  RooPlot *plot_sidebandsMC_alpha = workspace_->var("CMS_hzz2l2q_mZZ")->frame();
  
  sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins));
  workspace_->pdf("background")->plotOn(plot_sidebandsMC_alpha, LineColor(kRed));
  
  plot_sidebandsMC_alpha->Draw();
  
  std::string outdir = get_outdir();
  std::string mkdir_command = "mkdir -p " + outdir;
  char canvasName[400];
  sprintf( canvasName, "%s/mZZ_sidebandsMC_alpha_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  std::string canvasName_str(canvasName);
  std::string canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
  
  c1->SetLogy();
  canvasName_str += "_log";
  canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
  c1->SetLogy(false);
  
  delete plot_sidebandsMC_alpha;
  
  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT MC Signal (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;
  ofsMCName = get_fitResultsName( 0, "FreeInitial" );
  workspace_->argSet("cutOff,beta,m,wdth,n,alpha").readFromFile(ofsMCName.c_str());
  workspace_->pdf("background")->fitTo(signalMC,SumW2Error(kTRUE));
    
  //write variables with apppropriate fixed values
  ofsMCName = get_fitResultsName( btagCategory, "MCSignal" );
  workspace_->var("cutOff")->setConstant(kTRUE);
  workspace_->var("beta")->setConstant(kTRUE);
  workspace_->var("m")->setConstant(kTRUE);
  workspace_->var("n")->setConstant(kTRUE);
  workspace_->set("vars")->writeToFile(ofsMCName.c_str());
  
  plot_sidebandsMC_alpha = workspace_->var("CMS_hzz2l2q_mZZ")->frame(mZZmin_, mZZmax_, nBins);
  
  sidebandsMC_alpha.plotOn(plot_sidebandsMC_alpha, Binning(nBins),MarkerColor(kRed));
  signalMC.plotOn(plot_sidebandsMC_alpha, Binning(nBins));
  
  workspace_->pdf("background")->plotOn(plot_sidebandsMC_alpha, LineColor(kBlack));
  ofsMCName = get_fitResultsName( btagCategory, "MC" );
  workspace_->argSet("cutOff,beta,m,wdth,n,alpha").readFromFile(ofsMCName.c_str());
  workspace_->pdf("background")->plotOn(plot_sidebandsMC_alpha, LineColor(kRed));
  
  plot_sidebandsMC_alpha->Draw();
  
  sprintf( canvasName, "%s/mZZ_sidebandsMCvsSignal_alpha_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  canvasName_str = canvasName;
  canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
  
  c1->SetLogy();
  canvasName_str += "_log";
  canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
  c1->Clear();
  c1->SetLogy(false);

  TFile* file_alpha = 0;
  char fitResultsFileName[500];
  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_%s_PU%s.root", dataset_.c_str(), btagCategory, leptType.c_str(), PUType_.c_str());
  file_alpha = TFile::Open(fitResultsFileName, "recreate");
  file_alpha->cd();
  h1_alpha->Write();
  tree_sidebandsMC_alpha->Write();
  file_alpha->Close();
  
  delete plot_sidebandsMC_alpha;
  delete c1;
} 





RooFitResult* SidebandFitter::fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed , std::string init) {

  std::string outdir = get_outdir();
  

  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    exit(109);
  }

  char cut_base[500];
  sprintf( cut_base, "nBTags==%d %s", btagCategory, leptType_cut.c_str());
  char cut_sidebands[500];
  sprintf( cut_sidebands, "%s && ( (mZjj>60. && mZjj<75.)||(mZjj>105. && mZjj<130.) )", cut_base);
  char cut_signal[500];
  sprintf( cut_signal, "%s && ( mZjj>75. && mZjj<105. )", cut_base);
  

  //RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",treeDATA,RooArgSet(*eventWeight,*CMS_hzz2l2q_mZZ,*nBTags,*mZjj),cut_sidebands);
  RooDataSet signalDATA("signalDATA","signalDATA",treeDATA,*(workspace_->set("tree")),cut_signal);
  char treeName_DATA[200];
  sprintf( treeName_DATA, "sidebandsDATA_alpha");
  std::string treeName_DATA_str(treeName_DATA);
  std::cout << "Correcting sidebands (DATA): " << std::endl;
  TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( treeDATA, h1_alpha, btagCategory, treeName_DATA_str );
  //tree_sidebandsDATA_alpha->GetBranch("mZZ")->SetName("CMS_hzz2l2q_mZZ"); 
  RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,*(workspace_->set("treeAlpha")),cut_sidebands,"eventWeight_alpha");


  TFile* file_alpha = 0;
  
  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT ALPHA-CORRECTED DATA SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;
  
  std::string ofsMCName = get_fitResultsName( btagCategory, init );
  workspace_->argSet("cutOff,beta,m,wdth,n,alpha").readFromFile(ofsMCName.c_str());

  RooFitResult *r_sidebandsDATA_alpha;
  r_sidebandsDATA_alpha = workspace_->pdf("background")->fitTo(sidebandsDATA_alpha, SumW2Error(kTRUE), Save());
  char fitResultName[200];
  if( leptType!="ALL" )
    sprintf( fitResultName, "fitResults_%dbtag_%s", btagCategory, leptType.c_str() );
  else 
    sprintf( fitResultName, "fitResults_%dbtag", btagCategory );
  r_sidebandsDATA_alpha->SetName(fitResultName);
  
  
  RooPlot *plot_sidebandsDATA_alpha = workspace_->var("CMS_hzz2l2q_mZZ")->frame();
  
  sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha);
  workspace_->pdf("background")->plotOn(plot_sidebandsDATA_alpha);
  sidebandsDATA_alpha.plotOn(plot_sidebandsDATA_alpha);
    
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
  
  //now decorrelate parameters:
  char diagonalizerName[200];
  sprintf( diagonalizerName, "CMS_hzz2l2q_bkg_%db", btagCategory);
  PdfDiagonalizer diago(diagonalizerName, workspace_, *r_sidebandsDATA_alpha);
  RooAbsPdf *background_decorr = diago.diagonalize(*(workspace_->pdf("background")));
  RooFitResult *r_sidebandsDATA_alpha_decorr = background_decorr->fitTo(sidebandsDATA_alpha, SumW2Error(kTRUE), Save());
  char fitResultName_eig[200];
  if( leptType!="ALL" )
    sprintf( fitResultName_eig, "%s_%s_decorr", fitResultName, leptType.c_str() );
  else 
    sprintf( fitResultName_eig, "%s_decorr", fitResultName );
  r_sidebandsDATA_alpha_decorr->SetName(fitResultName_eig);
  background_decorr->SetName("background_decorr");
  workspace_->import(*background_decorr, RooFit::RecycleConflictNodes());
  

  RooPlot *plot_rot = workspace_->var("CMS_hzz2l2q_mZZ")->frame();
  sidebandsDATA_alpha.plotOn(plot_rot);
  workspace_->pdf("background")->plotOn(plot_rot,RooFit::LineColor(kRed));
  background_decorr->plotOn(plot_rot,RooFit::LineColor(kBlue), RooFit::LineStyle(2));
  TCanvas* c_rot = new TCanvas("rot", "", 600, 600);
  c_rot->cd();
  plot_rot->Draw();
  char canvasName_rot[200];
  sprintf( canvasName_rot, "%s/checkrot_%dbtag.eps", outdir.c_str(),btagCategory );
  c_rot->SaveAs(canvasName_rot);

  ContourPlot("alpha","wdth" ,r_sidebandsDATA_alpha);
  sprintf( canvasName_rot, "%s/checkrot_Elli_%dbtag.eps", outdir.c_str(),btagCategory );
  c_rot->SaveAs(canvasName_rot);
  char var1[50];
  char var2[50];
  sprintf(var1,"CMS_hzz2l2q_bkg_%db_eig0",btagCategory);
  sprintf(var2,"CMS_hzz2l2q_bkg_%db_eig1",btagCategory);
  ContourPlot(var1,var2 ,r_sidebandsDATA_alpha_decorr);
  
  sprintf( canvasName_rot, "%s/checkrot_ElliDecorr_%dbtag.eps", outdir.c_str(),btagCategory );
  c_rot->SaveAs(canvasName_rot);
  std::string ofsDATAName = get_fitResultsName( btagCategory, "DATA" );
  workspace_->set("vars")->writeToFile(ofsDATAName.c_str());
  ofsDATAName = get_fitResultsName( btagCategory, "DATADCORR" );
  RooArgSet tmpset(r_sidebandsDATA_alpha_decorr->floatParsFinal());
  tmpset.writeToFile(ofsDATAName.c_str());
  
  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  PLOT DATA SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;
  
  c1->Clear();
  c1->SetLogy(false);


  RooPlot *plot_signalDATA = workspace_->var("CMS_hzz2l2q_mZZ")->frame();

  background_decorr->plotOn(plot_signalDATA,VisualizeError(*r_sidebandsDATA_alpha_decorr,2.0,kFALSE),FillColor(kYellow),Normalization(sidebandsDATA_alpha.sumEntries()));
  background_decorr->plotOn(plot_signalDATA,VisualizeError(*r_sidebandsDATA_alpha_decorr,1.0,kFALSE),FillColor(kGreen),Normalization(sidebandsDATA_alpha.sumEntries()));
  background_decorr->plotOn(plot_signalDATA, LineColor(kRed),Normalization(sidebandsDATA_alpha.sumEntries()));
  signalDATA.plotOn(plot_signalDATA);


  plot_signalDATA->Draw();
  
  sprintf( canvasName, "%s/mZZ_signalDATA_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  plot_signalDATA->SetMinimum(0.1);
  plot_signalDATA->SetMaximum(1500);
  plot_signalDATA->Draw();  
  c1->SetLogy();

  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
  
  delete plot_signalDATA;


  //dynamic_cast<RooRealVar*>(r_sidebandsDATA_alpha->floatParsFinal().find("wdth"))->setError(10.);
  
  
  char fitResultsFileName[500];
  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_%s_PU%s.root", dataset_.c_str(), btagCategory, leptType.c_str(), PUType_.c_str());
  file_alpha = TFile::Open(fitResultsFileName, "UPDATE");
  file_alpha->cd();
  tree_sidebandsDATA_alpha->Write();
  r_sidebandsDATA_alpha->Write(); // << this one is final
  //r_sidebandsDATA_alpha_decorr->Write(); // this one needs to get adjusted errors
  char wnam[50];
  sprintf(wnam,"fitWorkspace_%dbtag",btagCategory);
  workspace_->SetName(wnam);
  workspace_->Write();
  file_alpha->Close();


  
  return r_sidebandsDATA_alpha_decorr; // we still need this to adjust the errors.

}


void SidebandFitter::fitPseudo( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed , std::string init) {

  std::string outdir = get_outdir();
  

  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    exit(109);
  }
  
  char cut_base[500];
  sprintf( cut_base, "nBTags==%d %s", btagCategory, leptType_cut.c_str());
  char cut_sidebands[500];
  sprintf( cut_sidebands, "%s && ( (mZjj>60. && mZjj<75.)||(mZjj>105. && mZjj<130.) )", cut_base);
  char cut_signal[500];
  sprintf( cut_signal, "%s && ( mZjj>75. && mZjj<105. )", cut_base);
  
  char treeName_DATA[200];
  sprintf( treeName_DATA, "sidebandsDATA_alpha");
  std::string treeName_DATA_str(treeName_DATA);
  std::cout << "Correcting sidebands (DATA): " << std::endl;
  TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( treeDATA, h1_alpha, btagCategory, treeName_DATA_str );
  RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,*(workspace_->set("treeAlpha")),cut_sidebands,"eventWeight_alpha");

  std::string ofsMCName = get_fitResultsName( btagCategory, init );
  workspace_->argSet("cutOff,beta,m,n").readFromFile(ofsMCName.c_str());// read fixed parameters
  
  char var1[50];
  char var2[50];
  sprintf(var1,"CMS_hzz2l2q_bkg_%db_eig0",btagCategory);
  sprintf(var2,"CMS_hzz2l2q_bkg_%db_eig1",btagCategory);
  char both[100];
  sprintf(both,"%s,%s",var1,var2);
  ofsMCName = get_fitResultsName( btagCategory, "DATADCORR" );
  workspace_->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value

  RooFitResult* r_pseudo = workspace_->pdf("background_decorr")->fitTo(sidebandsDATA_alpha, SumW2Error(kTRUE),Save(), PrintLevel(-1));
  
  char indexstring[200];
  sprintf(indexstring,"DATADCORR%d",seed);
  ofsMCName = get_fitResultsName( btagCategory, indexstring );
  RooArgSet tmpset(r_pseudo->floatParsFinal());
  tmpset.writeToFile(ofsMCName.c_str());

  delete  r_pseudo;
  delete tree_sidebandsDATA_alpha;

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
  //tree->SetBranchAddress( "mZZ", &mZZ );
  tree->SetBranchAddress( "CMS_hzz2l2q_mZZ", &mZZ );
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
    //float alpha = 1.;

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



void SidebandFitter::pseudoMassge(int btagCategory , const std::string& leptType, std::string init, RooFitResult* r_nominal){

  std::string outdir = get_outdir();
  

  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    exit(109);
  }
  
  std::string ofsMCName = get_fitResultsName( btagCategory, init );
  workspace_->argSet("cutOff,beta,m,n").readFromFile(ofsMCName.c_str());// read fixed parameters
  char var1[50];
  char var2[50];
  sprintf(var1,"CMS_hzz2l2q_bkg_%db_eig0",btagCategory);
  sprintf(var2,"CMS_hzz2l2q_bkg_%db_eig1",btagCategory);
  char both[100];
  sprintf(both,"%s,%s",var1,var2);
  ofsMCName = get_fitResultsName( btagCategory, "DATADCORR" );
  workspace_->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value
  
  RooPlot *plot_MCbkg = workspace_->var("CMS_hzz2l2q_mZZ")->frame();
  

  std::vector<float> vals;
  std::vector<float> vals1;
  std::vector<float> vals2;
  vals.reserve(500);
  vals1.reserve(500);
  vals2.reserve(500);

  char indexstring[200];
  for(int i =0 ; i < 500 ; i++){
    sprintf(indexstring,"DATADCORR%d",i);
    ofsMCName = get_fitResultsName( btagCategory, indexstring );
    workspace_->argSet(both).readFromFile(ofsMCName.c_str());
    vals1.push_back(workspace_->var(var1)->getVal());
    vals2.push_back(workspace_->var(var2)->getVal());
    workspace_->pdf("background_decorr")->plotOn(plot_MCbkg,LineWidth(1),LineColor(1));
    std::string mkdir_command = "rm " + ofsMCName;
    system(mkdir_command.c_str());
  }
  RooCurve* upper = new RooCurve();
  RooCurve* lower = new RooCurve();
  float min = plot_MCbkg->GetXaxis()->GetXmin();
  float max = plot_MCbkg->GetXaxis()->GetXmax();
  float range=max-min;
  for(unsigned int x =0 ; x < 200 ; x++  ){
    float xval = min +x*range/200.;
    for(int i =0 ; i < 500 ; i++){
      vals.push_back(((RooCurve*)(plot_MCbkg->getObject(i)))->interpolate(xval));
    }

    std::sort(vals.begin(),vals.end());
    lower->addPoint(xval,vals[83]);
    upper->addPoint(xval,vals[416]);    
    vals.clear();
  }

  ofsMCName = get_fitResultsName( btagCategory, "DATADCORR" );
  workspace_->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value
  double x1= workspace_->var(var1)->getVal();
  double x2= workspace_->var(var2)->getVal();
  double e1= workspace_->var(var1)->getError();
  double e2= workspace_->var(var2)->getError();
  
  lower->SetLineWidth(2);
  lower->SetLineColor(2);
  upper->SetLineWidth(2);
  upper->SetLineColor(2);


  workspace_->pdf("background_decorr")->plotOn(plot_MCbkg,VisualizeError(*r_nominal,2.0,kFALSE),FillColor(kYellow));
  workspace_->pdf("background_decorr")->plotOn(plot_MCbkg,VisualizeError(*r_nominal,1.0,kFALSE),FillColor(kGreen));
  workspace_->pdf("background_decorr")->plotOn(plot_MCbkg);

  plot_MCbkg->addPlotable(lower);
  plot_MCbkg->addPlotable(upper);
  
  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  plot_MCbkg->Draw();
  char canvasName[400];
  sprintf( canvasName, "%s/mZZ_sidenbandData_alphaVar_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  std::string canvasName_str(canvasName);
  std::string canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  plot_MCbkg->SetMinimum(0.0001);
  plot_MCbkg->SetMaximum(.3);
  plot_MCbkg->Draw();  
  c1->SetLogy();

  canvasName_str += "_log";
  canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->Clear();
  c1->SetLogy(false);
  delete plot_MCbkg;


  // plotting value distributions
  std::sort(vals1.begin(),vals1.end());
  std::sort(vals2.begin(),vals2.end());
  TH1F* histo1= new TH1F("test1","test1",100,vals1[0],vals1[499]);
  char tit[20];
  sprintf(tit,"#alpha_{%d}",btagCategory);
  histo1->GetXaxis()->SetTitle(tit);
  histo1->GetYaxis()->SetTitle("Nr trials");
  TH1F* histo2= new TH1F("test2","test2",100,vals2[0],vals2[499]);
  sprintf(tit,"#beta_{%d}",btagCategory);
  histo2->GetYaxis()->SetTitle("Nr trials");
  histo2->GetXaxis()->SetTitle(tit);
  TLine* line = new TLine();
  line->SetLineColor(2);
  line->SetLineWidth(2);

  for(int i =0 ; i < 500 ; i++){
    histo1->Fill(vals1[i]);
    histo2->Fill(vals2[i]);
  }
  histo1->Fit("gaus");
  double s1 = histo1->GetFunction("gaus")->GetParameter("Sigma"); // we need this to update the errors
  histo2->Fit("gaus");
  double s2 = histo2->GetFunction("gaus")->GetParameter("Sigma"); // we need this to update the errors
  
  histo1->Draw();
  line->SetLineColor(2);
  line->DrawLine(x1,histo1->GetMinimum(),x1,histo1->GetMaximum());
  line->SetLineColor(4);
  line->DrawLine(x1+e1,histo1->GetMinimum(),x1+e1,histo1->GetMaximum());
  line->DrawLine(x1-e1,histo1->GetMinimum(),x1-e1,histo1->GetMaximum());

  sprintf( canvasName, "%s/alphaVar_par1_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  canvasName_str = canvasName;
  canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  histo2->Draw();
  line->SetLineColor(2);
  line->DrawLine(x2,histo2->GetMinimum(),x2,histo2->GetMaximum());
  line->SetLineColor(4);
  line->DrawLine(x2+e2,histo2->GetMinimum(),x2+e2,histo2->GetMaximum());
  line->DrawLine(x2-e2,histo2->GetMinimum(),x2-e2,histo2->GetMaximum());

  sprintf( canvasName, "%s/alphaVar_par2_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  canvasName_str = canvasName;
  canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());



  std::cout << "increasing errors from " << e1 << " : " << e2  << std::endl;
  double newErr1 = sqrt(e1*e1+s1*s1);
  double newErr2 = sqrt(e2*e2+s2*s2);  
  dynamic_cast<RooRealVar*>(r_nominal->floatParsFinal().find(var1))->setError(newErr1);
  dynamic_cast<RooRealVar*>(r_nominal->floatParsFinal().find(var2))->setError(newErr2);
  std::cout << "to " << dynamic_cast<RooRealVar*>(r_nominal->floatParsFinal().find(var1))->getError() << " : " << dynamic_cast<RooRealVar*>(r_nominal->floatParsFinal().find(var2))->getError() << std::endl;

  std::cout << " writing adjusted fitresult "<< std::endl;
  TFile* file_alpha = 0;
  char fitResultsFileName[500];
  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_%s_PU%s.root", dataset_.c_str(), btagCategory, leptType.c_str(), PUType_.c_str());
  file_alpha = TFile::Open(fitResultsFileName, "UPDATE");
  file_alpha->cd();
  r_nominal->Write();
  file_alpha->Close();


  // control plot for new error
  plot_MCbkg = workspace_->var("CMS_hzz2l2q_mZZ")->frame();
  TRandom3 random;
  for(int i =0 ; i < 500 ; i++){
    ofsMCName = get_fitResultsName( btagCategory,  "DATADCORR" );
    workspace_->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value
    workspace_->var(var1)->setVal(random.Gaus(x1,newErr1));
    workspace_->var(var2)->setVal(random.Gaus(x2,newErr2));
    workspace_->pdf("background_decorr")->plotOn(plot_MCbkg,LineWidth(1),LineColor(1));
  }
  vals.clear();
  upper = new RooCurve();
  lower = new RooCurve();
  for(unsigned int x =0 ; x < 200 ; x++  ){
    float xval = min +x*range/200.;
    for(int i =0 ; i < 500 ; i++){
      vals.push_back(((RooCurve*)(plot_MCbkg->getObject(i)))->interpolate(xval));
    }
    std::sort(vals.begin(),vals.end());
    lower->addPoint(xval,vals[83]);
    upper->addPoint(xval,vals[416]);
    vals.clear();
  }
  
  ofsMCName = get_fitResultsName( btagCategory,  "DATADCORR" );
  workspace_->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value
  workspace_->pdf("background_decorr")->plotOn(plot_MCbkg,VisualizeError(*r_nominal,2.0,kFALSE),FillColor(kYellow));
  workspace_->pdf("background_decorr")->plotOn(plot_MCbkg,VisualizeError(*r_nominal,1.0,kFALSE),FillColor(kGreen));
  workspace_->pdf("background_decorr")->plotOn(plot_MCbkg);

  
  lower->SetLineWidth(2);
  lower->SetLineColor(2);
  upper->SetLineWidth(2);
  upper->SetLineColor(2);

  plot_MCbkg->addPlotable(lower);
  plot_MCbkg->addPlotable(upper);

  plot_MCbkg->Draw();  
  sprintf( canvasName, "%s/mZZ_sidenbandData_alphaVarAdded_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  canvasName_str = canvasName;
  canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  plot_MCbkg->SetMinimum(0.0001);
  plot_MCbkg->SetMaximum(.3);
  plot_MCbkg->Draw();  
  c1->SetLogy();

  canvasName_str += "_log";
  canvasName_eps = canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
 
  delete plot_MCbkg;
  

}


RooPlot* SidebandFitter::ContourPlot(std::string var1, std::string var2, RooFitResult* r){

  Double_t x1= workspace_->var(var1.c_str())->getVal();
  Double_t x2= workspace_->var(var2.c_str())->getVal();
  Double_t s1= workspace_->var(var1.c_str())->getError();
  Double_t s2= workspace_->var(var2.c_str())->getError();
  Double_t rho= r->correlation(var1.c_str(),var2.c_str());

  RooEllipse *contour= new RooEllipse("contour",x1,x2,s1,s2,rho,500);
  contour->SetLineWidth(2) ;


  RooPlot *plot = new RooPlot(*(workspace_->var(var1.c_str())),*(workspace_->var(var2.c_str())),x1-2*s1,x1+2*s1,
                              x2-2*s2,x2+2*s2);
  //RooPlot *plot = new RooPlot(*alpha_0,*beta_0,40,100,-1.5,.4);
  plot->addPlotable(contour);

  r->plotOn(plot,*(workspace_->var(var1.c_str())),*(workspace_->var(var2.c_str())),"ME12");

  plot->Draw();
  return plot;

}


  
// this method returns only rate:
Double_t SidebandFitter::get_backgroundNormalization( int nbtags, const std::string& leptType, const std::string& data_mc ) {

  std::pair<float,float> rate_and_error = this->get_backgroundNormalizationAndError( nbtags, leptType, data_mc );

  return rate_and_error.first;

}


// this method return both rate (first) and error on rate (second):
std::pair<Double_t,Double_t> SidebandFitter::get_backgroundNormalizationAndError( int nbtags, const std::string& leptType, const std::string& data_mc ) {

  
  // open fit results file:
  char fitResultsFileName[200];
  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_ALL_PU%s.root", dataset_.c_str(), nbtags, PUType_.c_str());
  TFile* fitResultsFile = TFile::Open(fitResultsFileName);


  // get alpha-corrected tree:
  std::string treeName = "sidebands" + data_mc + "_alpha"; 
  TTree* treeSidebandsDATA_alphaCorr = (TTree*)fitResultsFile->Get(treeName.c_str());
  

  // compute expected BG yield from observed sideband events
  Double_t rate_background;
  Double_t rate_background_error;

  // special treatment for 2 btag category:
  // fix relative ele/mu normalization by taking MC ratio
  // in order to minimize sideband fluctuations in data

  if( nbtags==2 && leptType!="ALL" ) { 

    TTree* treeMC = (TTree*)fitResultsFile->Get("sidebandsMC_alpha");

    TH1D* h1_mZZ_signalMC_ELE = new TH1D("mZZ_signalMC_ELE", "", 65, mZZmin_, mZZmax_ );
    TH1D* h1_mZZ_signalMC_MU = new TH1D("mZZ_signalMC_MU", "", 65, mZZmin_, mZZmax_ );
    h1_mZZ_signalMC_ELE->Sumw2();
    h1_mZZ_signalMC_MU->Sumw2();

    char signalCutMC[500];
    sprintf( signalCutMC, "eventWeight*(mZjj>75. && mZjj<105. && leptType==0 && nBTags==%d)", nbtags );
    treeMC->Project("mZZ_signalMC_MU", "CMS_hzz2l2q_mZZ", signalCutMC);
    sprintf( signalCutMC, "eventWeight*(mZjj>75. && mZjj<105. && leptType==1 && nBTags==%d)", nbtags );
    treeMC->Project("mZZ_signalMC_ELE", "CMS_hzz2l2q_mZZ", signalCutMC);

    float eleMC = h1_mZZ_signalMC_ELE->Integral();
    float muMC = h1_mZZ_signalMC_MU->Integral();
    float ratioMC = (leptType=="MU") ? eleMC/muMC : muMC/eleMC;

    TH1D* h1_mZZ_sidebandsDATA = new TH1D("mZZ_sidebandsDATA", "", 65, mZZmin_, mZZmax_ );
    h1_mZZ_sidebandsDATA->Sumw2();
    char sidebandsCut_alpha[500];
    sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d)", nbtags ); //electrons+muons
    treeSidebandsDATA_alphaCorr->Project("mZZ_sidebandsDATA", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
    double sumDATA = h1_mZZ_sidebandsDATA->IntegralAndError( h1_mZZ_sidebandsDATA->GetXaxis()->GetFirst(), h1_mZZ_sidebandsDATA->GetXaxis()->GetLast(), rate_background_error );

    rate_background = sumDATA / ( ratioMC+1.);
    rate_background_error /= ( ratioMC+1.);

  } else { //nbtags =0,1 or 2-tag but ele+mu

    TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, mZZmin_, mZZmax_ );
    h1_mZZ_sidebands_alpha->Sumw2();
    char sidebandsCut_alpha[500];
    if( leptType=="ALL" )
      sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d)", nbtags );
    else
      sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d && leptType==%d)", nbtags, SidebandFitter::convert_leptType(leptType) );
    treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
    rate_background = h1_mZZ_sidebands_alpha->IntegralAndError( h1_mZZ_sidebands_alpha->GetXaxis()->GetFirst(), h1_mZZ_sidebands_alpha->GetXaxis()->GetLast(), rate_background_error );

  }


  fitResultsFile->Close();

  std::pair<Double_t,Double_t> rate_and_error;
  rate_and_error.first = rate_background;
  rate_and_error.second = rate_background_error;

  return rate_and_error;

}




RooDataSet* SidebandFitter::get_observedDataset( RooRealVar* CMS_hzz2l2q_mZZ, const std::string& leptType_str, int nbtags ) {


  std::string dataFileName = "HZZlljjRM_DATA_" + dataset_ + "_optLD_looseBTags_v2_ALL.root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  TTree* tree_data = (TTree*)dataFile->Get("tree_passedEvents");

  

  RooRealVar nBTags("nBTags","nBTags",-1.,3.);
  RooRealVar eventWeight("eventWeight","eventWeight",0,100.);
  RooRealVar mZjj("mZjj","mZjj",0,150.);
  RooRealVar leptType("leptType","lepton type",-1,2);

  char selection[900];
  if( leptType_str=="ALL" )
    sprintf( selection, "mZjj>75. && mZjj<105. && nBTags==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f", nbtags, mZZmin_, mZZmax_ );
  else {
    int leptType_int = SidebandFitter::convert_leptType(leptType_str);
    sprintf( selection, "mZjj>75. && mZjj<105. && nBTags==%d && leptType==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f", nbtags, leptType_int, mZZmin_, mZZmax_ );
  }


  RooFormulaVar rooselection("selection", selection, RooArgList(*CMS_hzz2l2q_mZZ,nBTags,mZjj,leptType));
  RooDataSet *dataset_obs = new RooDataSet("dataset_obs", "dataset_obs", tree_data,
                                     RooArgSet(*CMS_hzz2l2q_mZZ, nBTags, mZjj, leptType, eventWeight),
                                     rooselection, "eventWeight");


  return dataset_obs;

}



 


int SidebandFitter::convert_leptType( const std::string& leptType ) {
  
  if( leptType!="ELE" && leptType!="MU" ) {
    std::cout << "WARNING!!! LeptType '" << leptType << "' is NOT supported!!! Returning -1." << std::endl;
    return -1;
  }
  
  int leptType_int = (leptType=="MU" ) ? 0 : 1;
  
  return leptType_int;
    
} 
