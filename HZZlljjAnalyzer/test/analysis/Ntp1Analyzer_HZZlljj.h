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
   //Double_t computeQGLikelihood(const Double_t jtpt, Int_t ncharged, Int_t nneutral, Double_t PtD, Double_t r);




 private:

   int leptType_; //0: muon; 1: electron

   Float_t eZqqMC_;
   Float_t ptZqqMC_;
   Float_t etaZqqMC_;
   Float_t phiZqqMC_;

   Float_t eZllMC_;
   Float_t ptZllMC_;
   Float_t etaZllMC_;
   Float_t phiZllMC_;

   Float_t eHiggsMC_;
   Float_t ptHiggsMC_;
   Float_t etaHiggsMC_;
   Float_t phiHiggsMC_;

   Float_t eLept1_;
   Float_t ptLept1_;
   Float_t etaLept1_;
   Float_t phiLept1_;
   Int_t   chargeLept1_;

   Float_t eLept1Gen_;
   Float_t ptLept1Gen_;
   Float_t etaLept1Gen_;
   Float_t phiLept1Gen_;

   Float_t eLept2_;
   Float_t ptLept2_;
   Float_t etaLept2_;
   Float_t phiLept2_;
   Int_t   chargeLept2_;

   Float_t eLept2Gen_;
   Float_t ptLept2Gen_;
   Float_t etaLept2Gen_;
   Float_t phiLept2Gen_;

   Float_t  ptJetBest1_;
   Float_t   eJetBest1_;
   Float_t phiJetBest1_;
   Float_t etaJetBest1_;
   Float_t rmsCandJetBest1_;
   Float_t ptDJetBest1_;
   Int_t nChargedJetBest1_;
   Int_t nNeutralJetBest1_;

   Float_t  ptJetBest2_;
   Float_t   eJetBest2_;
   Float_t phiJetBest2_;
   Float_t etaJetBest2_;
   Float_t rmsCandJetBest2_;
   Float_t ptDJetBest2_;
   Int_t nChargedJetBest2_;
   Int_t nNeutralJetBest2_;

   Float_t  ptJetRecoil_;
   Float_t   eJetRecoil_;
   Float_t phiJetRecoil_;
   Float_t etaJetRecoil_;
   Float_t rmsCandJetRecoil_;
   Float_t ptDJetRecoil_;
   Int_t nChargedJetRecoil_;
   Int_t nNeutralJetRecoil_;

   Float_t  ptJetLead_;
   Float_t   eJetLead_;
   Float_t phiJetLead_;
   Float_t etaJetLead_;

   Float_t  ptJetLead2_;
   Float_t   eJetLead2_;
   Float_t phiJetLead2_;
   Float_t etaJetLead2_;

   Float_t  ptJetLead3_;
   Float_t   eJetLead3_;
   Float_t phiJetLead3_;
   Float_t etaJetLead3_;

   Int_t nPairs_;

   Int_t  iJet1_[50];
   Float_t  ptJet1_[50];
   Float_t   eJet1_[50];
   Float_t phiJet1_[50];
   Float_t etaJet1_[50];

   Float_t ptDJet1_[50];
   Float_t rmsCandJet1_[50];
   Int_t nChargedJet1_[50];
   Int_t nNeutralJet1_[50];
   Float_t QGlikelihoodJet1_[50];

   Float_t  eChargedHadronsJet1_[50];
   Float_t  ePhotonsJet1_[50];
   Float_t  eNeutralEmJet1_[50];
   Float_t  eNeutralHadronsJet1_[50];
   Float_t  eMuonsJet1_[50];
   Float_t  eElectronsJet1_[50];
   Float_t  eHFHadronsJet1_[50];
   Float_t  eHFEMJet1_[50];

   Int_t  nChargedHadronsJet1_[50];
   Int_t  nPhotonsJet1_[50];
   Int_t  nNeutralHadronsJet1_[50];
   Int_t  nMuonsJet1_[50];
   Int_t  nElectronsJet1_[50];
   Int_t  nHFHadronsJet1_[50];
   Int_t  nHFEMJet1_[50];

   Float_t  ptJet1Gen_[50];
   Float_t   eJet1Gen_[50];
   Float_t phiJet1Gen_[50];
   Float_t etaJet1Gen_[50];

   Int_t  nPFCand1_;
   Float_t  ePFCand1_[100];
   Float_t  ptPFCand1_[100];
   Float_t  etaPFCand1_[100];
   Float_t  phiPFCand1_[100];
   Int_t  particleTypePFCand1_[100];

   Int_t  iJet2_[50];
   Float_t  ptJet2_[50];
   Float_t   eJet2_[50];
   Float_t phiJet2_[50];
   Float_t etaJet2_[50];

   Float_t ptDJet2_[50];
   Float_t rmsCandJet2_[50];
   Int_t nChargedJet2_[50];
   Int_t nNeutralJet2_[50];
   Float_t QGlikelihoodJet2_[50];

   Float_t  eChargedHadronsJet2_[50];
   Float_t  ePhotonsJet2_[50];
   Float_t  eNeutralEmJet2_[50];
   Float_t  eNeutralHadronsJet2_[50];
   Float_t  eMuonsJet2_[50];
   Float_t  eElectronsJet2_[50];
   Float_t  eHFHadronsJet2_[50];
   Float_t  eHFEMJet2_[50];

   Int_t  nChargedHadronsJet2_[50];
   Int_t  nPhotonsJet2_[50];
   Int_t  nNeutralHadronsJet2_[50];
   Int_t  nMuonsJet2_[50];
   Int_t  nElectronsJet2_[50];
   Int_t  nHFHadronsJet2_[50];
   Int_t  nHFEMJet2_[50];

   Float_t  ptJet2Gen_[50];
   Float_t   eJet2Gen_[50];
   Float_t phiJet2Gen_[50];
   Float_t etaJet2Gen_[50];


   Int_t  nPFCand2_;
   Float_t  ePFCand2_[100];
   Float_t  ptPFCand2_[100];
   Float_t  etaPFCand2_[100];
   Float_t  phiPFCand2_[100];
   Int_t  particleTypePFCand2_[100];

   Int_t nPart_;
   Float_t  ptPart_[20];
   Float_t   ePart_[20];
   Float_t phiPart_[20];
   Float_t etaPart_[20];
   Int_t pdgIdPart_[20];


   Float_t epfMet_;
   Float_t phipfMet_;

   TH1F* h1_nEvents_vs_ptEle; 
   TH1F* h1_nEvents_vs_ptMuon; 
   TH1F* h1_passed_vs_ptEle; 
   TH1F* h1_passed_vs_ptMuon; 
   TH1F* h1_deltaRmatching_muons; 
   TH1F* h1_deltaRmatching_electrons; 
   TH1F* h1_deltaRmatching_jet_parton; 
   TH1F* h1_deltaRmatching_genjet_parton; 
   TH1F* h1_deltaRmatching_jet_genjet; 
   TH1F* h1_deltaRmatching_jet_leptonParton;
   TH1F* h1_nJets30;
// TH1F* h1_indexMatchedJet;
// TH1F* h1_indexMatched05Jet;
// TH1F* h1_nMatched_per_event;
// TH1F* h1_nMatched05_per_event;
// TH1F* h1_pdgIdParton2;
// TH1F* h1_ptHadronicZ; 
// TH1F* h1_deltaRqq; 

   bool DEBUG_VERBOSE_;

};




#endif
