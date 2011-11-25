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
#include "RooFitResult.h"
#include "RooWorkspace.h"

#include "RooEllipse.h"

//#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

#include "PdfDiagonalizer.h"
#include <algorithm>


using namespace RooFit;


SidebandFitter::SidebandFitter( const std::string& dataset, const std::string& PUType, const std::string& init, const std::string& flags ) {

  dataset_ = dataset;
  PUType_ = PUType;
  init_ = init;
  flags_ = flags;

  mZZmin_ = 160.;
  mZZmax_ = 800.;


  double a0 = -1.395;
  double w0 = 85.73;


  CMS_hzz2l2q_mZZ_ = new RooRealVar("CMS_hzz2l2q_mZZ", "m_{lljj}", mZZmin_, mZZmax_, "GeV");


  // ------------------------ fermi ------------------------------
  cutOff_ = new RooRealVar("cutOff","position of fermi",191.12,175.,220.);
  beta_ = new RooRealVar("beta","width of fermi",4.698,0.,30.);
  fermi_ = new RooFermi("fermi","fermi function",*CMS_hzz2l2q_mZZ_,*cutOff_,*beta_);

  // -------------------- crystal ball ---------------------------
  m_ = new RooRealVar("m","m",200.17,190.,300.);
  wdth_ = new RooRealVar("wdth","wdth",w0,-200.,200.);
  n_ = new RooRealVar("n","n",13.067,0.,100.);
  alpha_ = new RooRealVar("alpha","alpha",a0,-200.,200.); 


  CBShape_ = new RooCBShape("CB","Crystal ball",*CMS_hzz2l2q_mZZ_,*m_,*wdth_,*alpha_,*n_);

  background_ = new RooProdPdf("background","background",RooArgSet(*fermi_,*CBShape_));
 
}



TH1D* SidebandFitter::getAlphaHisto( int btagCategory, const std::string& leptType_str, TTree* treeMC ) {


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

  //treeMC->SetBranchAddress("mZZ",&mZZ);
  treeMC->SetBranchAddress("CMS_hzz2l2q_mZZ",&mZZ);
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


  // put bin values in pathological region equal to sherpa values:
  if( btagCategory==0 ) {
    double err4 = h1_alpha->GetBinError(4);
    double err5 = h1_alpha->GetBinError(5);
    h1_alpha->SetBinContent(4, 1.09979);
    h1_alpha->SetBinContent(5, 0.9638);
    h1_alpha->SetBinError(4, err4+0.0574);
    h1_alpha->SetBinError(5, err5+0.0412);
  }

  return (TH1D*)(h1_alpha->Clone());
  
}



RooFitResult* SidebandFitter::fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha ) {


  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    exit(109);
  }
  

  if( init_!="MCSignal" && init_!="MC" && init_!="DATA" ) {
    std::cout << "Can't initialize on '" << init_ << "'. (Don't know what it is.) Exiting." << std::endl;
    exit(311);
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

  RooRealVar* eventWeight = new RooRealVar("eventWeight", "event weight", 0., 2., "");
  RooRealVar* eventWeight_alpha = new RooRealVar("eventWeight_alpha", "event weight (alpha corrected)", 0., 2., "");
  RooRealVar* nBTags = new RooRealVar("nBTags", "number of BTags", -1., 2., "");
  RooRealVar* mZjj = new RooRealVar("mZjj", "mZjj", 60., 130., "GeV");


  std::cout << "Correcting sidebands (MC): " << std::endl;
  TTree* tree_sidebandsMC_alpha = correctTreeWithAlpha( treeMC, h1_alpha, btagCategory, "sidebandsMC_alpha" );
  std::cout << "Correcting sidebands (DATA): " << std::endl;
  TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( treeDATA, h1_alpha, btagCategory, "sidebandsDATA_alpha" );


  RooDataSet *sidebandsDATA_alpha = new RooDataSet("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*CMS_hzz2l2q_mZZ_,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");

  RooDataSet *fitDataset;
 
  if( init_=="MCSignal" ) {

    fitDataset = new RooDataSet("signalMC","signalMC",treeMC,RooArgSet(*eventWeight,*CMS_hzz2l2q_mZZ_,*nBTags,*mZjj),cut_signal,"eventWeight");

  } else if( init_=="MC" ) {

    fitDataset = new RooDataSet("sidebandsMC_alpha","sidebandsMC_alpha",tree_sidebandsMC_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*CMS_hzz2l2q_mZZ_,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");
  } else if( init_=="DATA" ) {
    
    fitDataset = sidebandsDATA_alpha;

  }





  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();



  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIXING PARAMETERS WITH A FIT TO " << init_ << " (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);
  fitDataset->write("test.txt");
  
  // lets try this:
  if( btagCategory==2 ) {
    // set "good" inital values:
    cutOff_->setVal(187.383);
    beta_->setVal(3.43718);
    m_->setVal(218.161);
  }
  RooFitResult *r_sidebandsMC_alpha = background_->fitTo(*fitDataset, SumW2Error(kTRUE), Save());


  std::string ofsMCName = get_fitResultsName( btagCategory, init_ );
  ofstream ofsMC(ofsMCName.c_str());

  ofsMC << "beta " << beta_->getVal() << " " << beta_->getError() << std::endl;
  ofsMC << "cutOff " << cutOff_->getVal() << " " << cutOff_->getError() << std::endl;
  ofsMC << "m " << m_->getVal() << " " << m_->getError() << std::endl;
  ofsMC << "wdth " << wdth_->getVal() << " " << wdth_->getError() << std::endl;
  ofsMC << "alpha " << alpha_->getVal() << " " << alpha_->getError() << std::endl;
  ofsMC << "n " << n_->getVal() << " " << n_->getError() << std::endl;

  ofsMC.close();

  RooPlot *fitPlot = CMS_hzz2l2q_mZZ_->frame(mZZmin_, mZZmax_, nBins);

  fitDataset->plotOn(fitPlot, Binning(nBins));

  background_->plotOn(fitPlot, LineColor(kRed));
  fitDataset->plotOn(fitPlot, Binning(nBins));

  fitPlot->Draw();

  char canvasName[400];
  sprintf( canvasName, "%s/mZZ_%s_alpha_%dbtag_%s", outdir.c_str(), init_.c_str(), btagCategory, leptType.c_str());
  std::string* canvasName_str = new std::string(canvasName);
  std::string canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  delete fitPlot;


  // fix parameters
  cutOff_->setConstant(kTRUE);
  beta_->setConstant(kTRUE);
  m_->setConstant(kTRUE);
  n_->setConstant(kTRUE);


  // and now (re)fit data sidebands to fix alpha and wdth


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT ALPHA-CORRECTED DATA SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);



  // just to be sure:
  wdth_->setConstant(kFALSE);
  alpha_->setConstant(kFALSE);


  RooFitResult *r_sidebandsDATA_alpha = background_->fitTo(*sidebandsDATA_alpha, SumW2Error(kTRUE), Save());
  char fitResultName[200];
  if( leptType!="ALL" )
    sprintf( fitResultName, "fitResults_%dbtag_%s", btagCategory, leptType.c_str() );
  else 
    sprintf( fitResultName, "fitResults_%dbtag", btagCategory );
  r_sidebandsDATA_alpha->SetName(fitResultName);


  RooPlot *plot_sidebandsDATA_alpha = CMS_hzz2l2q_mZZ_->frame(mZZmin_, mZZmax_, nBins);

  sidebandsDATA_alpha->plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

  //background_->plotOn(plot_sidebandsDATA_alpha,VisualizeError(*r_sidebandsDATA_alpha,2.0,kFALSE),FillColor(kYellow), Normalization(sidebandsDATA_alpha->sumEntries()));
  //background_->plotOn(plot_sidebandsDATA_alpha,VisualizeError(*r_sidebandsDATA_alpha,1.0,kFALSE),FillColor(kGreen), Normalization(sidebandsDATA_alpha->sumEntries()));
  background_->plotOn(plot_sidebandsDATA_alpha, LineColor(kRed));
  sidebandsDATA_alpha->plotOn(plot_sidebandsDATA_alpha, Binning(nBins));

  plot_sidebandsDATA_alpha->Draw();

  sprintf( canvasName, "%s/mZZ_sidebandsDATA_alpha_%dbtag_%s", outdir.c_str(), btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());



  char fitWorkspaceName[200];
  sprintf( fitWorkspaceName, "fitWorkspace_%dbtag", btagCategory );
  //RooWorkspace* fitWorkspace = new RooWorkspace(fitWorkspaceName, fitWorkspaceName);
  fitWorkspace = new RooWorkspace(fitWorkspaceName, fitWorkspaceName);



  //now decorrelate parameters:
  char diagonalizerName[200];
  sprintf( diagonalizerName, "CMS_hzz2l2q_bkg_%db", btagCategory);
  PdfDiagonalizer diago(diagonalizerName, fitWorkspace, *r_sidebandsDATA_alpha);
  background_decorr_ = diago.diagonalize(*background_);
  background_decorr_->SetName("background_decorr");
  RooFitResult *r_sidebandsDATA_alpha_decorr = background_decorr_->fitTo(*sidebandsDATA_alpha, SumW2Error(kFALSE), Save());
  char fitResultName_eig[200];
  if( leptType!="ALL" )
    sprintf( fitResultName_eig, "%s_%s_decorr", fitResultName, leptType.c_str() );
  else 
    sprintf( fitResultName_eig, "%s_decorr", fitResultName );
  r_sidebandsDATA_alpha_decorr->SetName(fitResultName_eig);

  //import both pdfs in the workspace:
  fitWorkspace->import(*CMS_hzz2l2q_mZZ_);
  fitWorkspace->import(*background_);
  fitWorkspace->import(*background_decorr_, RooFit::RecycleConflictNodes());
  
  //plot the correlated/decorrelated uncertainties

  TCanvas* c_rot = new TCanvas("rot", "", 600, 600);
  c_rot->cd();
  ContourPlot(alpha_,wdth_ ,r_sidebandsDATA_alpha);
  char canvasName_rot[200];
  sprintf( canvasName_rot, "%s/checkrot_Elli_%dbtag.eps", outdir.c_str(),btagCategory );
  c_rot->SaveAs(canvasName_rot);
  
  char var1[50];
  char var2[50];
  sprintf(var1,"CMS_hzz2l2q_bkg_%db_eig0",btagCategory);
  sprintf(var2,"CMS_hzz2l2q_bkg_%db_eig1",btagCategory);
  ContourPlot(fitWorkspace->var(var1),fitWorkspace->var(var2),r_sidebandsDATA_alpha_decorr);
  sprintf( canvasName_rot, "%s/checkrot_ElliDecorr_%dbtag.eps", outdir.c_str(),btagCategory );
  c_rot->SaveAs(canvasName_rot);
  delete c_rot;
  std::string ofsDATAName = get_fitResultsName( btagCategory, "DATADCORR" );
  RooArgSet tmpset(r_sidebandsDATA_alpha_decorr->floatParsFinal());
  tmpset.writeToFile(ofsDATAName.c_str());


  //control plot  to check functions agree before/after rotation

  RooPlot *plot_rot = CMS_hzz2l2q_mZZ_->frame(mZZmin_, mZZmax_, nBins);
  sidebandsDATA_alpha->plotOn(plot_rot);
  background_->plotOn(plot_rot,RooFit::LineColor(kRed));
  background_decorr_->plotOn(plot_rot,RooFit::LineColor(kBlue), RooFit::LineStyle(2));
  c_rot = new TCanvas("rot", "", 600, 600);
  c_rot->cd();
  plot_rot->Draw();
  sprintf( canvasName, "checkrot_%dbtag.eps", btagCategory );
  c_rot->SaveAs(canvasName);





  ofsDATAName = get_fitResultsName( btagCategory, "DATA" );
  ofstream ofsDATA(ofsDATAName.c_str());


  ofsDATA << "beta " << beta_->getVal() << " " << beta_->getError() << std::endl;
  ofsDATA << "cutOff " << cutOff_->getVal() << " " << cutOff_->getError() << std::endl;
  ofsDATA << "m " << m_->getVal() << " " << m_->getError() << std::endl;
  ofsDATA << "wdth " << wdth_->getVal() << " " << wdth_->getError() << std::endl;
  ofsDATA << "alpha " << alpha_->getVal() << " " << alpha_->getError() << std::endl;
  ofsDATA << "n " << n_->getVal() << " " << n_->getError() << std::endl;

  ofsDATA.close();


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIT DATA SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  //fix shape:
  cutOff_->setConstant(kTRUE);
  beta_->setConstant(kTRUE);
  m_->setConstant(kTRUE);
  wdth_->setConstant(kTRUE);
  n_->setConstant(kTRUE);
  alpha_->setConstant(kTRUE);



  c1->Clear();
  c1->SetLogy(false);


  RooDataSet signalDATA("signalDATA","signalDATA",treeDATA,RooArgSet(*eventWeight,*CMS_hzz2l2q_mZZ_,*nBTags,*mZjj),cut_signal);
  RooPlot *plot_signalDATA = CMS_hzz2l2q_mZZ_->frame(mZZmin_, mZZmax_, nBins);

  background_->plotOn(plot_signalDATA, LineColor(kRed), Normalization(sidebandsDATA_alpha->sumEntries()));
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

  delete plot_signalDATA;






  std::string fitResultsFileName = get_fitResultsRootFileName( btagCategory, leptType );
  TFile* file_alpha = TFile::Open(fitResultsFileName.c_str(), "recreate");
  file_alpha->cd();
  h1_alpha->Write();
  tree_sidebandsDATA_alpha->Write();
  tree_sidebandsMC_alpha->Write();
  r_sidebandsDATA_alpha->Write();
  r_sidebandsDATA_alpha_decorr->Write();
  fitWorkspace->Write();
  file_alpha->Close();



//   delete  fitDataset; 
//   delete sidebandsDATA_alpha;


  delete eventWeight;
  delete eventWeight_alpha;
  delete nBTags;
  delete mZjj;
  delete c1;
  //delete r_sidebandsMC_alpha;


  return r_sidebandsDATA_alpha_decorr;

}



std::string SidebandFitter::get_fitResultsName( int nbtags, const std::string& init ) {

  std::string outdir = get_outdir();

  char fitResultsName[600];
  sprintf( fitResultsName, "%s/fitresults%s_%dbtag.txt", outdir.c_str(), init.c_str(), nbtags);
  std::string returnString(fitResultsName);

  return returnString;

}



std::string SidebandFitter::get_outdir() {

  std::string returnString = "FitSidebands_" + dataset_ + "_fit" + init_ + flags_;

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


  
// this method returns only rate:
Double_t SidebandFitter::get_backgroundNormalization( int nbtags, const std::string& leptType, const std::string& data_mc, float mZZmin, float mZZmax ) {

  std::pair<float,float> rate_and_error = this->get_backgroundNormalizationAndError( nbtags, leptType, data_mc, mZZmin, mZZmax );

  return rate_and_error.first;

}



// this method return both rate (first) and error on rate (second):
std::pair<Double_t,Double_t> SidebandFitter::get_backgroundNormalizationAndError( int nbtags, const std::string& leptType, const std::string& data_mc, float mZZmin, float mZZmax ) {

  // use fit ranges as default
  if( mZZmin < 0. ) mZZmin = mZZmin_;
  if( mZZmax < 0. ) mZZmax = mZZmax_;
  
  // open fit results file:
  std::string fitResultsFileName = get_fitResultsRootFileName( nbtags, "ALL" );
  TFile* fitResultsFile = TFile::Open(fitResultsFileName.c_str());


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

    TH1D* h1_mZZ_signalMC_ELE = new TH1D("mZZ_signalMC_ELE", "", 65, mZZmin, mZZmax );
    TH1D* h1_mZZ_signalMC_MU = new TH1D("mZZ_signalMC_MU", "", 65, mZZmin, mZZmax );
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

    TH1D* h1_mZZ_sidebandsDATA = new TH1D("mZZ_sidebandsDATA", "", 65, mZZmin, mZZmax );
    h1_mZZ_sidebandsDATA->Sumw2();
    char sidebandsCut_alpha[500];
    sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f)", nbtags, mZZmin, mZZmax ); //electrons+muons
    treeSidebandsDATA_alphaCorr->Project("mZZ_sidebandsDATA", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
    double sumDATA = h1_mZZ_sidebandsDATA->IntegralAndError( h1_mZZ_sidebandsDATA->GetXaxis()->GetFirst(), h1_mZZ_sidebandsDATA->GetXaxis()->GetLast(), rate_background_error );

    rate_background = sumDATA / ( ratioMC+1.);
    rate_background_error /= ( ratioMC+1.);

  } else { //nbtags =0,1 or 2-tag but ele+mu

    TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, mZZmin, mZZmax );
    h1_mZZ_sidebands_alpha->Sumw2();
    char sidebandsCut_alpha[500];
    if( leptType=="ALL" )
      sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f)", nbtags, mZZmin, mZZmax );
    else
      sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d && leptType==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f)", nbtags, SidebandFitter::convert_leptType(leptType), mZZmin, mZZmax );
    treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
    rate_background = h1_mZZ_sidebands_alpha->IntegralAndError( h1_mZZ_sidebands_alpha->GetXaxis()->GetFirst(), h1_mZZ_sidebands_alpha->GetXaxis()->GetLast(), rate_background_error );

  }


  fitResultsFile->Close();

  std::pair<Double_t,Double_t> rate_and_error;
  rate_and_error.first = rate_background;
  rate_and_error.second = rate_background_error;

  return rate_and_error;

}



//// this method return both rate (first) and error on rate (second):
//std::pair<Double_t,Double_t> SidebandFitter::get_backgroundNormalizationAndError( int nbtags, const std::string& leptType, const std::string& data_mc ) {
//
//  
//  // open fit results file:
//  char fitResultsFileName[200];
//  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_ALL_PU%s.root", dataset_.c_str(), nbtags, PUType_.c_str());
//  TFile* fitResultsFile = TFile::Open(fitResultsFileName);
//
//
//  // get alpha-corrected tree:
//  std::string treeName = "sidebands" + data_mc + "_alpha"; 
//  TTree* treeSidebandsDATA_alphaCorr = (TTree*)fitResultsFile->Get(treeName.c_str());
//  
//
//  // compute expected BG yield from observed sideband events
//  Double_t rate_background;
//  Double_t rate_background_error;
//
//  // special treatment for 2 btag category:
//  // fix relative ele/mu normalization by taking MC ratio
//  // in order to minimize sideband fluctuations in data
//
//  if( nbtags==2 && leptType!="ALL" ) { 
//
//    TTree* treeMC = (TTree*)fitResultsFile->Get("sidebandsMC_alpha");
//
//    TH1D* h1_mZZ_signalMC_ELE = new TH1D("mZZ_signalMC_ELE", "", 65, 150., 800.);
//    TH1D* h1_mZZ_signalMC_MU = new TH1D("mZZ_signalMC_MU", "", 65, 150., 800.);
//    h1_mZZ_signalMC_ELE->Sumw2();
//    h1_mZZ_signalMC_MU->Sumw2();
//
//    char signalCutMC[500];
//    sprintf( signalCutMC, "eventWeight*(mZjj>75. && mZjj<105. && leptType==0 && nBTags==%d)", nbtags );
//    treeMC->Project("mZZ_signalMC_MU", "CMS_hzz2l2q_mZZ", signalCutMC);
//    sprintf( signalCutMC, "eventWeight*(mZjj>75. && mZjj<105. && leptType==1 && nBTags==%d)", nbtags );
//    treeMC->Project("mZZ_signalMC_ELE", "CMS_hzz2l2q_mZZ", signalCutMC);
//
//    float eleMC = h1_mZZ_signalMC_ELE->Integral();
//    float muMC = h1_mZZ_signalMC_MU->Integral();
//    float ratioMC = (leptType=="MU") ? eleMC/muMC : muMC/eleMC;
//
//    TH1D* h1_mZZ_sidebandsDATA = new TH1D("mZZ_sidebandsDATA", "", 65, 150., 800.);
//    char sidebandsCut_alpha[500];
//    sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d)", nbtags ); //electrons+muons
//    treeSidebandsDATA_alphaCorr->Project("mZZ_sidebandsDATA", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
//    double sumDATA = h1_mZZ_sidebandsDATA->Integral();
//
//    rate_background = sumDATA / ( ratioMC+1.);
//
//  } else { //nbtags =0,1 or 2-tag but ele+mu
//
//    TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, 150., 800.);
//    h1_mZZ_sidebands_alpha->Sumw2();
//    char sidebandsCut_alpha[500];
//    if( leptType=="ALL" )
//      sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d)", nbtags );
//    else
//      sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d && leptType==%d)", nbtags, SidebandFitter::convert_leptType(leptType) );
//    treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
//    rate_background = h1_mZZ_sidebands_alpha->IntegralAndError( h1_mZZ_sidebands_alpha->GetXaxis()->GetFirst(), h1_mZZ_sidebands_alpha->GetXaxis()->GetLast(), rate_background_error );
//
//  }
//
//
//  fitResultsFile->Close();
//
//  std::pair<Double_t,Double_t> rate_and_error;
//  rate_and_error.first = rate_background;
//  rate_and_error.second = rate_background_error;
//
//  return rate_and_error;
//
//}



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
    sprintf( selection, "mZjj>75. && mZjj<105. && nBTags==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f", nbtags, CMS_hzz2l2q_mZZ->getMin(), CMS_hzz2l2q_mZZ->getMax() );
  else {
    int leptType_int = SidebandFitter::convert_leptType(leptType_str);
    sprintf( selection, "mZjj>75. && mZjj<105. && nBTags==%d && leptType==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f", nbtags, leptType_int, CMS_hzz2l2q_mZZ->getMin(), CMS_hzz2l2q_mZZ->getMax() );
  }


  RooFormulaVar rooselection("selection", selection, RooArgList(*CMS_hzz2l2q_mZZ,nBTags,mZjj,leptType));
  RooDataSet *dataset_obs = new RooDataSet("dataset_obs", "dataset_obs", tree_data,
                                     RooArgSet(*CMS_hzz2l2q_mZZ, nBTags, mZjj, leptType, eventWeight),
                                     rooselection, "eventWeight");


  return dataset_obs;

}


std::string SidebandFitter::get_fitResultsRootFileName( int btagCategory, const std::string& leptType ) {

  char fitResultsFileName[500];
  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_%s_PU%s_fit%s%s.root", dataset_.c_str(), btagCategory, leptType.c_str(), PUType_.c_str(), init_.c_str(), flags_.c_str() );

  std::string fitResultsFileName_str(fitResultsFileName);

  return fitResultsFileName_str;

}



 


int SidebandFitter::convert_leptType( const std::string& leptType ) {
  
  if( leptType!="ELE" && leptType!="MU" ) {
    std::cout << "WARNING!!! LeptType '" << leptType << "' is NOT supported!!! Returning -1." << std::endl;
    return -1;
  }
  
  int leptType_int = (leptType=="MU" ) ? 0 : 1;
  
  return leptType_int;
    
} 

RooPlot* SidebandFitter::ContourPlot(RooRealVar* var1,RooRealVar* var2, RooFitResult* r){

  Double_t x1= var1->getVal();
  Double_t x2= var2->getVal();
  Double_t s1= var1->getError();
  Double_t s2= var2->getError();
  Double_t rho= r->correlation(var1->GetName(),var2->GetName());

  RooEllipse *contour= new RooEllipse("contour",x1,x2,s1,s2,rho,500);
  contour->SetLineWidth(2) ;


  RooPlot *plot = new RooPlot(*var1,*var2,x1-2*s1,x1+2*s1,x2-2*s2,x2+2*s2);
  plot->addPlotable(contour);

  r->plotOn(plot,*var1,*var2,"M12");

  plot->Draw();
  return plot;

}
void SidebandFitter::fitPseudo( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed) {

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
  RooRealVar* eventWeight = new RooRealVar("eventWeight", "event weight", 0., 2., "");
  RooRealVar* eventWeight_alpha = new RooRealVar("eventWeight_alpha", "event weight (alpha corrected)", 0., 2., "");
  RooRealVar* nBTags = new RooRealVar("nBTags", "number of BTags", -1., 2., "");
  RooRealVar* mZjj = new RooRealVar("mZjj", "mZjj", 60., 130., "GeV");

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
  RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,RooArgSet(*eventWeight,*eventWeight_alpha,*CMS_hzz2l2q_mZZ_,*nBTags,*mZjj),cut_sidebands,"eventWeight_alpha");

  std::string ofsMCName = get_fitResultsName( btagCategory, init_ );
  //workspace_->argSet("cutOff,beta,m,n").readFromFile(ofsMCName.c_str());// we assume that these are still correct

  char var1[50];
  char var2[50];
  sprintf(var1,"CMS_hzz2l2q_bkg_%db_eig0",btagCategory);
  sprintf(var2,"CMS_hzz2l2q_bkg_%db_eig1",btagCategory);
  char both[100];
  sprintf(both,"%s,%s",var1,var2);
  ofsMCName = get_fitResultsName( btagCategory, "DATADCORR" );
  fitWorkspace->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value

  RooFitResult* r_pseudo = fitWorkspace->pdf("background_decorr")->fitTo(sidebandsDATA_alpha, SumW2Error(kTRUE), Save(), PrintLevel(-1));

  char indexstring[200];
  sprintf(indexstring,"DATADCORR%d",seed);
  ofsMCName = get_fitResultsName( btagCategory, indexstring );
  RooArgSet tmpset(r_pseudo->floatParsFinal());
  tmpset.writeToFile(ofsMCName.c_str());

  delete  r_pseudo;
  delete tree_sidebandsDATA_alpha;
  delete eventWeight;
  delete eventWeight_alpha;
  delete nBTags;
  delete mZjj;

}

void SidebandFitter::pseudoMassge(int ntoys, int btagCategory , const std::string& leptType, RooFitResult* r_nominal){

  std::string outdir = get_outdir();
  std::string ofsMCName;

  std::string leptType_cut="";
  if( leptType=="MU" ) {
    leptType_cut=" && leptType==0";
  } else if( leptType=="ELE" ) {
    leptType_cut=" && leptType==1";
  } else if( leptType!="ALL" ) {
    std::cout << "Unknown leptType: '" << leptType << "'. Exiting." << std::endl;
    exit(109);
  }

  char var1[50];
  char var2[50];
  sprintf(var1,"CMS_hzz2l2q_bkg_%db_eig0",btagCategory);
  sprintf(var2,"CMS_hzz2l2q_bkg_%db_eig1",btagCategory);
  char both[100];
  sprintf(both,"%s,%s",var1,var2);
  ofsMCName = get_fitResultsName( btagCategory, "DATADCORR" );
  fitWorkspace->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value
  
  RooPlot *plot_MCbkg = CMS_hzz2l2q_mZZ_->frame();

  std::vector<float> vals;
  std::vector<float> vals1;
  std::vector<float> vals2;
  vals.reserve(500);
  vals1.reserve(500);
  vals2.reserve(500);

  char indexstring[200];
  for(int i =0 ; i < ntoys ; i++){
    sprintf(indexstring,"DATADCORR%d",i);
    ofsMCName = get_fitResultsName( btagCategory, indexstring );
    fitWorkspace->argSet(both).readFromFile(ofsMCName.c_str());
    vals1.push_back(fitWorkspace->var(var1)->getVal());
    vals2.push_back(fitWorkspace->var(var2)->getVal());
    fitWorkspace->pdf("background_decorr")->plotOn(plot_MCbkg,LineWidth(1),LineColor(1));
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
  fitWorkspace->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value
  double x1= fitWorkspace->var(var1)->getVal();
  double x2= fitWorkspace->var(var2)->getVal();
  double e1= fitWorkspace->var(var1)->getError();
  double e2= fitWorkspace->var(var2)->getError();

  
  lower->SetLineWidth(2);
  lower->SetLineColor(2);
  upper->SetLineWidth(2);
  upper->SetLineColor(2);


  fitWorkspace->pdf("background_decorr")->plotOn(plot_MCbkg,VisualizeError(*r_nominal,2.0,kFALSE),FillColor(kYellow));
  fitWorkspace->pdf("background_decorr")->plotOn(plot_MCbkg,VisualizeError(*r_nominal,1.0,kFALSE),FillColor(kGreen));
  fitWorkspace->pdf("background_decorr")->plotOn(plot_MCbkg);

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
  r_nominal->floatParsFinal().Print("v");
  dynamic_cast<RooRealVar*>(r_nominal->floatParsFinal().find(var1))->setError(newErr1);
  dynamic_cast<RooRealVar*>(r_nominal->floatParsFinal().find(var2))->setError(newErr2);
  std::cout << "to " << dynamic_cast<RooRealVar*>(r_nominal->floatParsFinal().find(var1))->getError() << " : " << dynamic_cast<RooRealVar*>(r_nominal->floatParsFinal().find(var2))->getError() << std::endl;

  std::cout << " writing adjusted fitresult "<< std::endl;
  std::string fitResultsFileName = get_fitResultsRootFileName( btagCategory, leptType );
  TFile* file_alpha = TFile::Open(fitResultsFileName.c_str(), "UPDATE");
  file_alpha->cd();
  r_nominal->Write(); //will over-write the old one
  file_alpha->Close();

  // control plot for new error
  plot_MCbkg = CMS_hzz2l2q_mZZ_->frame();
  TRandom3 random;
  for(int i =0 ; i < 500 ; i++){
    ofsMCName = get_fitResultsName( btagCategory,  "DATADCORR" );
    fitWorkspace->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value
    fitWorkspace->var(var1)->setVal(random.Gaus(x1,newErr1));
    fitWorkspace->var(var2)->setVal(random.Gaus(x2,newErr2));
    fitWorkspace->pdf("background_decorr")->plotOn(plot_MCbkg,LineWidth(1),LineColor(1));
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
  fitWorkspace->argSet(both).readFromFile(ofsMCName.c_str());// read nominal best fit value
  fitWorkspace->pdf("background_decorr")->plotOn(plot_MCbkg,VisualizeError(*r_nominal,2.0,kFALSE),FillColor(kYellow));
  fitWorkspace->pdf("background_decorr")->plotOn(plot_MCbkg,VisualizeError(*r_nominal,1.0,kFALSE),FillColor(kGreen));
  fitWorkspace->pdf("background_decorr")->plotOn(plot_MCbkg);

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

