#include "Ntp1Finalizer.h"
#include <iostream>
#include "TROOT.h"



Ntp1Finalizer::Ntp1Finalizer( const std::string& analyzerType, const std::string& dataset, const std::string& flags) {

  DEBUG_ = false;


  tree_ = new TChain("reducedTree");

  analyzerType_ = analyzerType;
  inputAnalyzerType_ = analyzerType;
  dataset_ = dataset;
  flags_ = flags;

  outFile_ = 0;

  nCounter_ = 0.;
  nCounterW_ = 0.;
  nCounterPU_ = 0.;

} //constructor



Ntp1Finalizer::~Ntp1Finalizer() {

  std::cout << std::endl << "-> Histograms saved in: " << outFile_->GetName() << std::endl;

  if( tree_!=0 ) {
    delete tree_;
    tree_=0;
  }

} //destructor



void Ntp1Finalizer::createOutputFile( const std::string& additionalFlags ) {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova_"+dataset_;
   else {
    if(dataset_!="") outfileName = analyzerType_ + "_" + dataset_;
    else outfileName = analyzerType_;
   }


   if( flags_!="" )
     outfileName = outfileName + "_" + flags_;
   if( additionalFlags!="" )
     outfileName = outfileName + "_" + additionalFlags;

   outfileName = outfileName + ".root";

   outFile_ = TFile::Open(outfileName.c_str(), "RECREATE");
   
   outFile_->cd();


}






void Ntp1Finalizer::addFile(const std::string& dataset, const std::string& selection) {

  std::string infileName = inputAnalyzerType_ + "_2ndLevelTreeW_" + dataset + ".root"; //the W is important: means that files have passed treatment (merging and weights)
  TFile* infile = TFile::Open(infileName.c_str(), "READ");
  if( infile==0 ) {
    std::cout << "---> Didn't find file: " << infileName << std::endl;
    std::cout << "---> Exiting!!" << std::endl;
    exit(11);
  }
  infile->cd();
  std::string treeName = infileName +"/reducedTree";
  tree_->Add(treeName.c_str());
  std::cout << "-> Added " << treeName << ". Tree now has " << tree_->GetEntries() << " entries." << std::endl;
  TH1F* h1_nCounter = (TH1F*)infile->Get("nCounter");
  TH1F* h1_nCounterW = (TH1F*)infile->Get("nCounterW");
  TH1F* h1_nCounterPU = (TH1F*)infile->Get("nCounterPU");
  if( h1_nCounter!= 0 && h1_nCounterW != 0 && h1_nCounterPU ) {
    nCounter_ += h1_nCounter->GetBinContent(1);
    nCounterW_ += h1_nCounterW->GetBinContent(1);
    nCounterPU_ += h1_nCounterPU->GetBinContent(1);
  } else {
    std::cout << std::endl << std::endl << "WARNING!! Dataset '" << dataset << "' has no nCounter information!!!" << std::endl;
  }
  infile->Close();


}



std::vector<TH1F*> Ntp1Finalizer::getResponseHistos(const std::string& name, unsigned binArraySize, Double_t* ptBins) {

  std::vector<TH1F*> returnVector;

  for( unsigned i=0; i<(binArraySize-1); ++i ) {
    char histoName[100];
    sprintf( histoName, "%s_ptBin_%.0f_%.0f", name.c_str(), ptBins[i], ptBins[i+1]);
    int nbins = 50;
    float xmin = 0.5;
    float xmax = 2.;
    TH1F* newHisto = new TH1F(histoName, "", nbins, xmin, xmax);
    newHisto->Sumw2();
    returnVector.push_back(newHisto);
  }

  return returnVector;

}


void Ntp1Finalizer::writeResponseHistos( TFile* file, std::vector<TH1F*> h1_response, std::string dirName ) {

  file->mkdir( dirName.c_str() );
  file->cd( dirName.c_str() );

  for( unsigned iHisto=0; iHisto<h1_response.size(); ++iHisto ) h1_response[iHisto]->Write();

  file->cd();

}



int Ntp1Finalizer::get_nBTags( const AnalysisJet& jet1, const AnalysisJet& jet2, BTagSFUtil* btsfutil, bool loosebtags ) {

  int nBTags;

  bool jet1_tagged_medium = jet1.btag_medium();
  bool jet1_tagged_loose  = jet1.btag_loose();
  bool jet2_tagged_medium = jet2.btag_medium();
  bool jet2_tagged_loose  = jet2.btag_loose();

//btsfutil->modifyBTagsWithSF( jet1_tagged_loose, jet1_tagged_medium, jet1.Pt(), jet1.Eta(), jet1.pdgIdPart );
//btsfutil->modifyBTagsWithSF( jet2_tagged_loose, jet2_tagged_medium, jet2.Pt(), jet2.Eta(), jet2.pdgIdPart );

  if( loosebtags ) {

    bool twoBTags  = ( jet1_tagged_medium && jet2_tagged_loose  )
                  || ( jet1_tagged_loose  && jet2_tagged_medium );
    bool oneBTag   = (!twoBTags) && ( jet1_tagged_loose || jet2_tagged_loose );

    if( twoBTags ) nBTags=2;
    else if( oneBTag ) nBTags=1;
    else nBTags=0;

  } else {

    bool twoBTags  = ( jet1_tagged_medium && jet2_tagged_medium );
    bool oneBTag   = (!twoBTags) && ( jet1_tagged_medium || jet2_tagged_medium );

    if( twoBTags ) nBTags=2;
    else if( oneBTag ) nBTags=1;
    else nBTags=0;

  }

  return nBTags;

}

