#include "Ntp1Finalizer.h"
#include <iostream>



Ntp1Finalizer::Ntp1Finalizer( const std::string& analyzerType, const std::string& dataset, const std::string& flags) {

  DEBUG_ = false;

  tree_ = new TChain("reducedTree");

  totalLumi_ = 0.;

  analyzerType_ = analyzerType;
  dataset_ = dataset;
  flags_ = flags;

  outFile_ = 0;

} //constructor



Ntp1Finalizer::~Ntp1Finalizer() {

  if( tree_!=0 ) {
    delete tree_;
    tree_=0;
  }

} //destructor



void Ntp1Finalizer::createOutputFile() {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova_"+dataset_;
   else {
    if(dataset_!="") outfileName = analyzerType_ + "_" + dataset_;
    else outfileName = analyzerType_;
   }


   if( flags_!="" )
     outfileName = outfileName + "_" + flags_ + ".root";

   outFile_ = TFile::Open(outfileName.c_str(), "RECREATE");
   
   outFile_->cd();

}


/*
void Ntp1Finalizer::set_outFile( const std::string& fileName, const std::string& suffix ) {

  std::string outfileName;

  if( fileName!="" ) {

    outfileName = fileName;

  } else {

    if( DEBUG_ ) outfileName = "provaHZZlljj_"+dataset_;
    else {
     if(dataset_!="") outfileName = "HZZlljj_"+dataset_;
     else outfileName = "HZZlljj";
    }

    if( suffix!="" ) outfileName += "_"+suffix;

    outfileName += ".root";

  }

  outFile_ = new TFile(outfileName.c_str(), "RECREATE");
  outFile_->cd();

}
*/




void Ntp1Finalizer::addFile(const std::string& dataset) {

  std::string infileName = "HZZlljj_2ndLevelTreeW_" + dataset + ".root"; //the W is important: means that files have passed treatment (merging and weights)
  std::string treeName = infileName +"/reducedTree";
  tree_->Add(treeName.c_str());
  std::cout << "-> Added " << treeName << ". Tree has " << tree_->GetEntries() << " entries." << std::endl;
  TFile* infile = TFile::Open(infileName.c_str(), "READ");
  TH1F* h1_lumi = (TH1F*)infile->Get("lumi");
  if( h1_lumi!=0 ) {
    totalLumi_ += h1_lumi->GetBinContent(1);
    std::cout << "\tTotal lumi: " << totalLumi_ << " ub-1" << std::endl;
  } else {
    std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
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
