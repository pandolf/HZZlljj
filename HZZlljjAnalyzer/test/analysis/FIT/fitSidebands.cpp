
//////////////////////////////////////////////////////////////////////
//
//    USAGE: ./fitSidebands [data_dataset] doPseudo
//
//    ( omit doPseudo in command line if you dont want to compute
//      alpha-related errors. done for quick debug. doPseudo is
//      needed in full analysis chain. )
//
//////////////////////////////////////////////////////////////////////


#include <cstdlib>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TString.h"
#include "TROOT.h"
#include <iostream>

#include "SidebandFitter.h"

TTree* getTreeData(std::string whichPeriod, bool useCMG);// useCMG=false => Franceescos trees
TTree* getTreeBckg(std::string whichPeriod, bool useCMG);// useCMG=true  => CMG         trees
TTree* conditionTree(std::string whichPeriod,TTree* inTree);

int main( int argc, char* argv[] ) {

  //get dataset from CommandLine
  std::string dataset = "LP11";
  if( argc>=2 ) {
    std::string dataset_str(argv[1]);
    dataset = dataset_str;
  }
  TString dataset_tstr(dataset);

  int nToys = 500;

  bool useCMG=false;
  std::string init("MCSignal");
  bool doPseudo=false;
  bool fitData=false;
  for(int i = 1 ; i< argc ; i++){
    std::string par(argv[i]);
    if(par == "CMG") {
      std::cout << "-> Will use CMG ntuples." << std::endl;
      useCMG=true;
    }
    if(par == "FixToSide") {
      std::cout << "-> Will fix fit params on MC sidebands." << std::endl;
      init="MC";
    }
    if(par=="fixtodata") { 
      std::cout << "-> Will fix fit params on DATA sidebands." << std::endl;
      fitData=true;
    }
    if(par=="DoPseudo"||par=="doPseudo"||par=="dopseudo") {
      std::cout << "-> Will compute alpha-error with " << nToys << " toys." << std::endl;
      doPseudo=true;
    }
  }
  
  // try to deduce the proper weight
  std::string PUReweighing = "Run2011A";
  if( dataset=="HR11" ) PUReweighing = "HR11";
  if( dataset=="HR11_v2" ) PUReweighing = "HR11_73pb";
  if( dataset_tstr.BeginsWith("Run2011B") ) PUReweighing = "Run2011B";
  if( dataset_tstr.BeginsWith("All11") ) PUReweighing = "All11";




  TTree* treeDATA = getTreeData(dataset,useCMG);
  TTree* chainMC = getTreeBckg(PUReweighing,useCMG);

  if(useCMG){
    treeDATA=conditionTree(PUReweighing,treeDATA);
    chainMC=conditionTree(PUReweighing,chainMC);
  }

  TRandom3* random = new TRandom3(0);
  
  for( int b=0 ; b < 3 ; b++){
    SidebandFitter *sf = new SidebandFitter(dataset,PUReweighing );
    char cutstring[200];
    sprintf(cutstring,"nBTags==%d",b);
    TTree* treeDATA_Xbtag = treeDATA->CopyTree(cutstring);

    TTree* treeMC_Xbtag = chainMC->CopyTree(cutstring);
    
    TH1D* alpha_Xbtag = sf->getAlphaHisto( b, "ALL", treeMC_Xbtag );

    if( !fitData ) sf->generateFixedPars(treeMC_Xbtag, b , "ALL", alpha_Xbtag );

    RooFitResult* fr_Xbtag = sf->fitSidebands( treeMC_Xbtag, treeDATA_Xbtag, b, "ALL", alpha_Xbtag,-1,init);

    if(doPseudo) {
      for(int i = 0 ; i <nToys ; i++) {
        std::cout << std::endl << "[ " << b << " b-tags ]  Toy: " << i << "/" << nToys << std::endl;
        TH1D* variedHisto = sf->shuffle(alpha_Xbtag, random ,"tmp");
        if( fitData )
          sf->fitPseudo( treeMC_Xbtag, treeDATA_Xbtag, b, "ALL", variedHisto,i,"DATA");
        else
          sf->fitPseudo( treeMC_Xbtag, treeDATA_Xbtag, b, "ALL", variedHisto,i,init);
        delete variedHisto;
      }
      sf->pseudoMassge(nToys, b,"ALL",init,fr_Xbtag);
    }

    delete treeDATA_Xbtag;
    delete treeMC_Xbtag;
    delete alpha_Xbtag;
    delete fr_Xbtag;
    delete sf;
  }


  delete random;
  return 0;

}

TTree* getTreeData(std::string whichPeriod, bool useCMG){
  std::vector<std::string> filenames;

  if(useCMG){ //CMG files
    if(whichPeriod=="Run2011A_FULL"){
      filenames.push_back(std::string("summer11_data_highmass.root/AngularInfo"));
    }
    if(whichPeriod=="Run2011B_FULL"){
      filenames.push_back(std::string("summer11_data_2011B.root/AngularInfo"));
    }
    if(whichPeriod=="All11"){
      filenames.push_back(std::string("summer11_data_highmass.root/AngularInfo"));
      filenames.push_back(std::string("summer11_data_2011B.root/AngularInfo"));
    }
  }
  else{ //Roma Files
    std::string datafileName = "HZZlljjRM_DATA_" + whichPeriod + "_optLD_looseBTags_v2_ALL.root/tree_passedEvents";    
    filenames.push_back(datafileName);
  }


  TChain chain;
  for(unsigned int i = 0; i < filenames.size(); i++){
    chain.Add(filenames[i].c_str());
  }

  gROOT->cd();
  
  return chain.CopyTree("");

}

TTree* getTreeBckg(std::string whichPeriod, bool useCMG){
  std::vector<std::string> filenames;

  if(useCMG){ //CMG files
    filenames.push_back(std::string("summer11_TTbarIncl_highmass.root/AngularInfo"));
    filenames.push_back(std::string("summer11_WZ_highmass.root/AngularInfo"));
    filenames.push_back(std::string("summer11_ZZ_highmass.root/AngularInfo"));
    filenames.push_back(std::string("summer11_ZJets_madgraph_highmass.root/AngularInfo"));
  }
  else{ //Roma Files
    std::string bgTreeName;
    bgTreeName = "HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_PU" + whichPeriod + "_ALL.root/tree_passedEvents";
    filenames.push_back(bgTreeName);
    bgTreeName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_PU" + whichPeriod + "_ALL.root/tree_passedEvents";
    filenames.push_back(bgTreeName);
    bgTreeName = "HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_PU" + whichPeriod + "_ALL.root/tree_passedEvents";
    filenames.push_back(bgTreeName);
  }


  TChain chain;
  for(unsigned int i = 0; i < filenames.size(); i++){
    chain.Add(filenames[i].c_str());
  }

  gROOT->cd();
  
  return chain.CopyTree("");
}

TTree* conditionTree(std::string whichPeriod,TTree* inTree){

  Int_t    leptTypeNew;
  Double_t leptTypeOld;
  inTree->SetBranchAddress( "lep", &leptTypeOld);
  Int_t    nBTagsNew;
  Double_t nBTagsOld;
  inTree->SetBranchAddress( "nBTags", &nBTagsOld );
  Float_t  mZZNew;
  Double_t mZZOld;
  inTree->SetBranchAddress( "mZZ", &mZZOld );
  Float_t  mZjjNew;
  Double_t mZjjOld;
  inTree->SetBranchAddress( "mJJ", &mZjjOld );
  Float_t  eventWeightNew;
  Double_t eventWeightOld;
  inTree->SetBranchAddress( "weight", &eventWeightOld );
  Bool_t   isSidebandsNew;
  
  TTree* newTree = new TTree("tmptree","tmptree");

  newTree->Branch( "eventWeight", &eventWeightNew, "newWeight/F" );
  newTree->Branch( "mZjj", &mZjjNew, "mZjj/F" );
  newTree->Branch( "CMS_hzz2l2q_mZZ", &mZZNew, "CMS_hzz2l2q_mZZ/F" ); 
  newTree->Branch( "nBTags", &nBTagsNew, "nBTags/I" );
  newTree->Branch( "leptType", &leptTypeNew, "leptTypeNew/I" );
  newTree->Branch( "isSidebands", &isSidebandsNew, "isSidebands/O" );
 
  int nentries = inTree->GetEntries();

  std::cout << " conditioning... " << std::endl;
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    inTree->GetEntry( iEntry );
    if( (iEntry % 10000)==0 ) std::cout << "Entry: " << iEntry << "/" << nentries << std::endl;

    leptTypeNew= leptTypeOld ? 0 : 1 ;// Francesco has opposite convention
    nBTagsNew= nBTagsOld;
    mZZNew= mZZOld;
    mZjjNew= mZjjOld;
    eventWeightNew= eventWeightOld;
    isSidebandsNew= ( (mZjjOld>60. &&mZjjOld < 75.) || (mZjjOld >105. && mZjjOld<130.));
    
    newTree->Fill();

  }

  TTree* tmp = newTree->CloneTree();
  delete newTree;
  delete inTree;
  return tmp;

}


