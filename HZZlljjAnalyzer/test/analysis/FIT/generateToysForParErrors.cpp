#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"

#include "SidebandFitter.h"




int main( int argc, char* argv[] ) {


  std::string dataset = "LP11";
  if( argc>1 ) {
    std::string dataset_str(argv[1]);
    dataset = dataset_str;
  }

  int nToys=10;
  if( argc>2 ) {
    std::string nToys_str(argv[2]);
    nToys = atoi(nToys_str.c_str());
  }

  int seed=13;
  if( argc>3 ) {
    std::string seed_str(argv[3]);
    seed = atoi(seed_str.c_str());
  }



  SidebandFitter *sf = new SidebandFitter(dataset);



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

  TH1D* alpha_0btag = sf->getAlphaHisto( 0, "ALL", treeMC_0btag );
  TH1D* alpha_1btag = sf->getAlphaHisto( 1, "ALL", treeMC_1btag );
  TH1D* alpha_2btag = sf->getAlphaHisto( 2, "ALL", treeMC_2btag );

  
  TH1D* h1_alphaToys_0btag = new TH1D("alphaToys_0btag", "", 100, -200., 50.); 
  TH1D* h1_wdthToys_0btag = new TH1D("wdthToys_0btag", "", 100, -1., 2.); 

  TH1D* h1_alphaToys_1btag = new TH1D("alphaToys_1btag", "", 100, -200., 50.); 
  TH1D* h1_wdthToys_1btag = new TH1D("wdthToys_1btag", "", 100, -1., 2.); 

  TH1D* h1_alphaToys_2btag = new TH1D("alphaToys_2btag", "", 100, -200., 50.); 
  TH1D* h1_wdthToys_2btag = new TH1D("wdthToys_2btag", "", 100, -1., 2.); 

  TRandom3* rand = new TRandom3(seed);

  for( unsigned iToy=0; iToy<nToys; ++iToy ) {

    std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;;
    std::cout << " +++ TOY:   " << iToy+1 << " / " << nToys << std::endl;
    std::cout << std::endl << std::endl;


    char histName[100];
    sprintf( histName, "alphaHist_0btag_%d", iToy );
    TH1D* randomAlpha_0btag = sf->shuffle( alpha_0btag, rand, histName );
    sprintf( histName, "alphaHist_1btag_%d", iToy );
    TH1D* randomAlpha_1btag = sf->shuffle( alpha_1btag, rand, histName );
    sprintf( histName, "alphaHist_2btag_%d", iToy );
    TH1D* randomAlpha_2btag = sf->shuffle( alpha_2btag, rand, histName );


    FitResults fr_0btag_iToy = sf->fitSidebands( treeMC_0btag, treeDATA_0btag, 0, "ALL", randomAlpha_0btag, iToy );
    FitResults fr_1btag_iToy = sf->fitSidebands( treeMC_1btag, treeDATA_1btag, 1, "ALL", randomAlpha_1btag, iToy );
    FitResults fr_2btag_iToy = sf->fitSidebands( treeMC_2btag, treeDATA_2btag, 2, "ALL", randomAlpha_2btag, iToy );

    if( fr_0btag_iToy.CB_theta!=0. ) {
      h1_alphaToys_0btag->Fill( fr_0btag_iToy.CB_alpha_rot );
      h1_wdthToys_0btag->Fill( fr_0btag_iToy.CB_wdth_rot );
    }

    if( fr_1btag_iToy.CB_theta!=0. ) {
      h1_alphaToys_1btag->Fill( fr_1btag_iToy.CB_alpha_rot );
      h1_wdthToys_1btag->Fill( fr_1btag_iToy.CB_wdth_rot );
    }

    if( fr_2btag_iToy.CB_theta!=0. ) {
      h1_alphaToys_2btag->Fill( fr_2btag_iToy.CB_alpha_rot );
      h1_wdthToys_2btag->Fill( fr_2btag_iToy.CB_wdth_rot );
    }

  }


//std::cout << std::endl << std::endl << std::endl;
//std::cout << "->  Now modifying fitresults files." << std::endl;

//TF1* f1_gaus = new TF1("gauss", "gaus");

//h1_alphaToys_0btag->Fit(f1_gaus, "Q");
//modifyFitResultError( "alpha_rot", f1_gaus->GetParameter(2), dataset, 0 );
//h1_alphaToys_1btag->Fit(f1_gaus, "Q");
//modifyFitResultError( "alpha_rot", f1_gaus->GetParameter(2), dataset, 1 );
//h1_alphaToys_2btag->Fit(f1_gaus, "Q");
//modifyFitResultError( "alpha_rot", f1_gaus->GetParameter(2), dataset, 2 );

//h1_wdthToys_0btag->Fit(f1_gaus, "Q");
//modifyFitResultError( "wdth_rot", f1_gaus->GetParameter(2), dataset, 0 );
//h1_wdthToys_1btag->Fit(f1_gaus, "Q");
//modifyFitResultError( "wdth_rot", f1_gaus->GetParameter(2), dataset, 1 );
//h1_wdthToys_2btag->Fit(f1_gaus, "Q");
//modifyFitResultError( "wdth_rot", f1_gaus->GetParameter(2), dataset, 2 );



  char outfileName[500];
  sprintf( outfileName, "fitParamErrors_%s_%d.root", dataset.c_str(), seed );
  TFile* outfile = TFile::Open(outfileName, "recreate");
  outfile->cd();

  h1_alphaToys_0btag->Write();
  h1_alphaToys_1btag->Write();
  h1_alphaToys_2btag->Write();

  h1_wdthToys_0btag->Write();
  h1_wdthToys_1btag->Write();
  h1_wdthToys_2btag->Write();

  outfile->Close();
  
  
  return 0;

}
