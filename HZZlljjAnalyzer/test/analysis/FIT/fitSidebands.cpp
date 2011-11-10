#include <cstdlib>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TString.h"
#include "TROOT.h"

#include "SidebandFitter.h"




int main( int argc, char* argv[] ) {


  std::string dataset = "LP11";
  if( argc==2 ) {
    std::string dataset_str(argv[1]);
    dataset = dataset_str;
  }


  TString dataset_tstr(dataset);
  std::string PUReweighing = "Run2011A";
  if( dataset=="HR11" ) PUReweighing = "HR11";
  if( dataset=="HR11_v2" ) PUReweighing = "HR11_73pb";
  if( dataset_tstr.BeginsWith("Run2011B") ) PUReweighing = "Run2011B";



  SidebandFitter *sf = new SidebandFitter(dataset, PUReweighing);


  std::string datafileName = "HZZlljjRM_DATA_" + dataset + "_optLD_looseBTags_v2_ALL.root";
  TFile* file_DATA = TFile::Open(datafileName.c_str());
  TTree* treeDATA = (TTree*)file_DATA->Get("tree_passedEvents");




  TChain* chainMC = new TChain("tree_passedEvents");
  std::string bgTreeName;
  bgTreeName = "HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_PU" + PUReweighing + "_ALL.root/tree_passedEvents";
  chainMC->Add(bgTreeName.c_str());
  bgTreeName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_PU" + PUReweighing + "_ALL.root/tree_passedEvents";
  chainMC->Add(bgTreeName.c_str());
  bgTreeName = "HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_PU" + PUReweighing + "_ALL.root/tree_passedEvents";
  chainMC->Add(bgTreeName.c_str());

  gROOT->cd(); //magic!

  TTree* treeDATA_0btag = treeDATA->CopyTree("nBTags==0");
  TTree* treeDATA_1btag = treeDATA->CopyTree("nBTags==1");
  TTree* treeDATA_2btag = treeDATA->CopyTree("nBTags==2");

  TTree* treeMC_0btag = chainMC->CopyTree("nBTags==0");
  TTree* treeMC_1btag = chainMC->CopyTree("nBTags==1");
  TTree* treeMC_2btag = chainMC->CopyTree("nBTags==2");

  TH1D* alpha_0btag = sf->getAlphaHisto( 0, "ALL", treeMC_0btag );
  TH1D* alpha_1btag = sf->getAlphaHisto( 1, "ALL", treeMC_1btag );
  TH1D* alpha_2btag = sf->getAlphaHisto( 2, "ALL", treeMC_2btag );

  RooFitResult* fr_0btag = sf->fitSidebands( treeMC_0btag, treeDATA_0btag, 0, "ALL", alpha_0btag );
  RooFitResult* fr_1btag = sf->fitSidebands( treeMC_1btag, treeDATA_1btag, 1, "ALL", alpha_1btag );
  RooFitResult* fr_2btag = sf->fitSidebands( treeMC_2btag, treeDATA_2btag, 2, "ALL", alpha_2btag );

  delete fr_0btag;
  delete fr_1btag;
  delete fr_2btag;

  return 0;

}
