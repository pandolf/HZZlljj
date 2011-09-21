#include "TFile.h"
#include <cstdlib>
#include <string>
#include "TTree.h"
#include "TH1F.h"

// the following script converts the bool, isSidebands, into an int, isSB, and 
//rewrites the tree.  It also, in this version, adds a wght according to alpha(mZZ)
// NOTE: this was meant to use on DATA sidebands only!!!

void readWriteTree(string fileName, string treeName){

  TFile* f = new TFile(fileName.c_str());
  TTree* t = (TTree*) f->Get(treeName.c_str());

  TFile* f_alpha0 = new TFile("alpha_0btag.root");
  TFile* f_alpha1 = new TFile("alpha_1btag.root");
  TFile* f_alpha2 = new TFile("alpha_2btag.root");
  TH1F *alpha0 = (TH1F*) f_alpha0->Get("alpha_MADGRAPH");
  TH1F *alpha1 = (TH1F*) f_alpha1->Get("alpha_MADGRAPH");
  TH1F *alpha2 = (TH1F*) f_alpha2->Get("alpha_MADGRAPH");

  float mZZ;
  float mZjj;
  float wght;
  int nBTags;
  bool isSidebands;
  int isSB;

  t->SetBranchAddress("mZZ",&mZZ);
  t->SetBranchAddress("mZjj",&mZjj);
  t->SetBranchAddress("nBTags",&nBTags);
  t->SetBranchAddress("isSidebands",&isSidebands);
  t->SetBranchAddress("eventWeight",&wght);
  
  string tempFile = "NEW_"+fileName;
  TFile *outFile  = new TFile(tempFile.c_str(),"RECREATE");
  TTree* selectedEvents = new TTree("selectedEvents","selectedEvents");

  selectedEvents->Branch("mZZ",&mZZ);
  selectedEvents->Branch("mZjj",&mZjj);
  selectedEvents->Branch("nBTags",&nBTags);
  selectedEvents->Branch("isSB",&isSB);
  selectedEvents->Branch("wght",&wght);
  
  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);

    if(isSidebands) isSB=1;
    else isSB=0;

    if(!(mZjj>75 && mZjj<105)){

      // 150 and 10 should not be hard coded... get from histograms!

      if(nBTags==0)
	wght=alpha0->GetBinContent(alpha0->FindBin(mZZ));
      else if(nBTags==1)
	wght=alpha1->GetBinContent(alpha1->FindBin(mZZ));
      else if(nBTags==2)
	wght=alpha2->GetBinContent(alpha2->FindBin(mZZ));
    }

    selectedEvents->Fill();

  }
  outFile->cd();
  selectedEvents->Write();

}
