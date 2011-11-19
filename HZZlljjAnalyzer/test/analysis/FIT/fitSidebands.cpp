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


  if( argc!=1 && argc!=2 && argc!=3 ) {
    std::cout << "USAGE: ./fitSidebands [data_dataset] [init]" << std::endl;
    exit(1111);
  }
 

  std::string dataset = "HR11_v2";
  if( argc==2 ) {
    std::string dataset_str(argv[1]);
    dataset = dataset_str;
  }

  std::string init="MC";
  if( argc==3 ) {
    std::string init_str(argv[2]);
    init = init_str;
  }

  std::cout << "-> Going to fix fit parameters on : " << init << std::endl;


  int nToys = 500;
  TRandom3* random = new TRandom3(0);


  TString dataset_tstr(dataset);
  std::string PUReweighing = "Run2011A";
  if( dataset=="HR11" ) PUReweighing = "HR11";
  if( dataset=="HR11_v2" ) PUReweighing = "HR11_73pb";
  if( dataset_tstr.BeginsWith("Run2011B") ) PUReweighing = "Run2011B";





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



  for( unsigned ibtag=0; ibtag<3; ++ibtag ) {

    SidebandFitter *sf = new SidebandFitter(dataset, PUReweighing, init);

    char btagCut[100];
    sprintf( btagCut, "nBTags==%d", ibtag );

    TTree* treeDATA_Xbtag = treeDATA->CopyTree(btagCut);

    TTree* treeMC_Xbtag = chainMC->CopyTree(btagCut);

    TH1D* alpha_Xbtag = sf->getAlphaHisto( ibtag, "ALL", treeMC_Xbtag );

    RooFitResult* fr = sf->fitSidebands( treeMC_Xbtag, treeDATA_Xbtag, ibtag, "ALL", alpha_Xbtag );

    for(int i = 0 ; i <nToys ; i++) {
      std::cout << std::endl << "[ " << ibtag << " b-tags ]  Toy: " << i << "/" << nToys << std::endl;
      TH1D* variedHisto = sf->shuffle(alpha_Xbtag, random ,"tmp");
      sf->fitPseudo( treeMC_Xbtag, treeDATA_Xbtag, ibtag, "ALL", variedHisto,i);
      delete variedHisto;
    }
    sf->pseudoMassge(nToys, ibtag,"ALL",fr);

    delete fr;
    delete sf;

  } //for btags

  return 0;

}
