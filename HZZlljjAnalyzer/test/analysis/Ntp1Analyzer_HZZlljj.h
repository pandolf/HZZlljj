//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_HZZlljj_h
#define Ntp1Analyzer_HZZlljj_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_HZZlljj : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_HZZlljj( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_HZZlljj();

   virtual void CreateOutputFile();
   virtual void Loop();



 private:

   int leptType_; //0: muon; 1: electron

   Float_t eLept1_;
   Float_t ptLept1_;
   Float_t etaLept1_;
   Float_t phiLept1_;

   Float_t eLept1Gen_;
   Float_t ptLept1Gen_;
   Float_t etaLept1Gen_;
   Float_t phiLept1Gen_;

   Float_t eLept2_;
   Float_t ptLept2_;
   Float_t etaLept2_;
   Float_t phiLept2_;

   Float_t eLept2Gen_;
   Float_t ptLept2Gen_;
   Float_t etaLept2Gen_;
   Float_t phiLept2Gen_;

   Float_t  ptJet1_;
   Float_t   eJet1_;
   Float_t  ptCorrJet1_;
   Float_t phiJet1_;
   Float_t etaJet1_;

   Float_t   ptJet1Gen_;
   Float_t    eJet1Gen_;
   Float_t  phiJet1Gen_;
   Float_t  etaJet1Gen_;

   Float_t  eChargedHadronsJet1_;
   Float_t  ePhotonsJet1_;
   Float_t  eNeutralHadronsJet1_;
   Float_t  eMuonsJet1_;
   Float_t  eElectronsJet1_;
   Float_t  eHFHadronsJet1_;
   Float_t  eHFEMJet1_;

   Int_t  nChargedHadronsJet1_;
   Int_t  nPhotonsJet1_;
   Int_t  nNeutralHadronsJet1_;
   Int_t  nMuonsJet1_;
   Int_t  nElectronsJet1_;
   Int_t  nHFHadronsJet1_;
   Int_t  nHFEMJet1_;

   Float_t  ptJet2_;
   Float_t   eJet2_;
   Float_t  ptCorrJet2_;
   Float_t phiJet2_;
   Float_t etaJet2_;

   Float_t   ptJet2Gen_;
   Float_t    eJet2Gen_;
   Float_t  phiJet2Gen_;
   Float_t  etaJet2Gen_;

   Float_t  eChargedHadronsJet2_;
   Float_t  ePhotonsJet2_;
   Float_t  eNeutralHadronsJet2_;
   Float_t  eMuonsJet2_;
   Float_t  eElectronsJet2_;
   Float_t  eHFHadronsJet2_;
   Float_t  eHFEMJet2_;

   Int_t  nChargedHadronsJet2_;
   Int_t  nPhotonsJet2_;
   Int_t  nNeutralHadronsJet2_;
   Int_t  nMuonsJet2_;
   Int_t  nElectronsJet2_;
   Int_t  nHFHadronsJet2_;
   Int_t  nHFEMJet2_;

   Int_t   pdgIdPartJet1_;
   Int_t   pdgIdPartJet2_;

   Float_t epfMet_;
   Float_t phipfMet_;

   TH1F* h1_nEvents_vs_ptEle; 
   TH1F* h1_nEvents_vs_ptMuon; 
   TH1F* h1_passed_vs_ptEle; 
   TH1F* h1_passed_vs_ptMuon; 
   TH1F* h1_deltaRmatching_muons; 
   TH1F* h1_deltaRmatching_electrons; 
   TH1F* h1_ptHadronicZ; 
   TH1F* h1_deltaRqq; 

   bool DEBUG_VERBOSE_;

};




#endif
