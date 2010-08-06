#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <cstdlib>


TChain* tree = 0;

struct EventsAndLumi {
  int nTotalEvents;
  float totalLumi;
};



EventsAndLumi addInput( const std::string& dataset );
float getWeight( const std::string& dataset, int nEvents );


int main( int argc, char* argv[] ) {

  if( argc!=2 ) {
    std::cout << "USAGE: ./merge_and_setWeights [dataset]" << std::endl;
    exit(917);
  }

  std::string dataset = argv[1];

  tree = new TChain("reducedTree");

  EventsAndLumi evlu;
  evlu = addInput( dataset );

  float weight = getWeight( dataset, evlu.nTotalEvents );

  // and now set the weights
  tree->SetBranchStatus( "eventWeight", 0 );
  
  std::string outfilename = "HZZlljj_2ndLevelTreeW_"+dataset+".root";
  TFile* outfile = new TFile(outfilename.c_str(), "recreate");
  outfile->cd();

  TH1F* h1_lumi = new TH1F("lumi", "", 1, 0., 1.);
  h1_lumi->SetBinContent(1, evlu.totalLumi);

  TTree* newTree = tree->CloneTree(0);
  Float_t newWeight;
  newTree->Branch( "eventWeight", &newWeight, "newWeight/F" );

  int nentries = tree->GetEntries();
  for( unsigned ientry = 0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( (ientry % 10000) ==0 ) std::cout << "Entry: " << ientry << " /" << nentries << std::endl;

    newWeight = weight;

    newTree->Fill();

  } //for entries

  h1_lumi->Write();
  newTree->Write();
  outfile->Write();
  outfile->Close();

  return 0;

}


EventsAndLumi addInput( const std::string& dataset ) {

  std::string infileName = "files_HZZlljj_2ndLevel_" + dataset + ".txt";
  TH1F* h1_lumi;
  TH1F* h1_nCounter;

  int totalEvents = 0;
  float totalLumi = 0.;

  //open from file.txt:
  FILE* iff = fopen(infileName.c_str(),"r");
  if(iff == 0) {
    std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9172);
//  infileName = "HZZlljj_2ndLevelTree_" + dataset + ".root";
//  std::string treeName = infileName +"/reducedTree";
//  tree->Add(treeName.c_str());
//  std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
//  TFile* infile = TFile::Open(infileName.c_str(), "READ");
//  h1_nCounter = (TH1F*)infile->Get("nCounter");
//  if( h1_nCounter!=0 ) {
//    totalEvents += h1_nCounter->GetBinContent(1);
//  } else {
//    std::cout << " WARNING! File '" << infileName << "' has no nCounter information. Skipping." << std::endl;
//  }
//  h1_lumi = (TH1F*)infile->Get("lumi");
//  if( h1_lumi!=0 ) {
//    totalLumi += h1_lumi->GetBinContent(1);
//    std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
//  } else {
//    std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
//  }
//  infile->Close();

  } else {

    char singleLine[500];

    while( fscanf(iff, "%s", singleLine) !=EOF ) {

      std::string rootfilename(singleLine);
      std::string treename = rootfilename + "/reducedTree";
      std::cout << "-> Added " << treename;
      tree->Add(treename.c_str());
      TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
      h1_nCounter = (TH1F*)infile->Get("nCounter");
      if( h1_nCounter!=0 ) {
        totalEvents += h1_nCounter->GetBinContent(1);
      } else {
        std::cout << " WARNING! File '" << infileName << "' has no nCounter information. Skipping." << std::endl;
      }
      h1_lumi = (TH1F*)infile->Get("lumi");
      if( h1_lumi!=0 ) {
        totalLumi += h1_lumi->GetBinContent(1);
        std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
      } else {
        std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
      }
      infile->Close();

    }
    fclose(iff);

  }

  EventsAndLumi evlu;
  evlu.nTotalEvents = totalEvents;
  evlu.totalLumi = totalLumi;

  return evlu;

} //addinput


float getWeight( const std::string& dataset, int nEvents ) {

  float xSection = -1.;

  if( dataset=="Zjets-madgraph" ) {
    xSection = 2800.;
  } else if( dataset=="HZZ_qqll_gluonfusion_M130" ) {
    xSection = 25.560*0.03913; //sigma x BR
  } else if( dataset=="HZZ_qqll_gluonfusion_M150" ) {
    xSection = 19.568*0.08234; //sigma x BR
  } else if( dataset=="HZZ_qqll_gluonfusion_M500" ) {
    xSection = 2.1914*0.2602; //sigma x BR
  } else {
    std::cout << std::endl;
    std::cout << "-> WARNING!! Dataset: '" << dataset << "' not present in database. Cross section unknown." << std::endl;
    std::cout << "-> Will set unitary weights." << std::endl;
    return 1.;
  }

  float weight = xSection/((float)nEvents);

  return weight;

}
