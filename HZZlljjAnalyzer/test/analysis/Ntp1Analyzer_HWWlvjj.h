//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->WW->lvjj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_HWWlvjj_h
#define Ntp1Analyzer_HWWlvjj_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_HWWlvjj : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_HWWlvjj( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_HWWlvjj();

   virtual void CreateOutputFile();
   virtual void Loop();
   //Double_t computeQGLikelihood(const Double_t jtpt, Int_t ncharged, Int_t nneutral, Double_t PtD, Double_t r);




 private:

   int leptType_; //0: muon; 1: electron

   Int_t TotEvent_;//###


   Float_t rhoPF_;
   Float_t eWqqMC_;
   Float_t ptWqqMC_;
   Float_t etaWqqMC_;
   Float_t phiWqqMC_;

   Float_t eWllMC_;
   Float_t ptWllMC_;
   Float_t etaWllMC_;
   Float_t phiWllMC_;

   Float_t eHiggsMC_;
   Float_t ptHiggsMC_;
   Float_t etaHiggsMC_;
   Float_t phiHiggsMC_;

   Float_t eQuark1_;
   Float_t ptQuark1_;
   Float_t etaQuark1_;
   Float_t phiQuark1_;
  
   Float_t eQuark2_;
   Float_t ptQuark2_;
   Float_t etaQuark2_;
   Float_t phiQuark2_;

   Float_t eLept_;
   Float_t ptLept_;
   Float_t etaLept_;
   Float_t phiLept_;
   Int_t   chargeLept_;

   Float_t eLeptMC_;
   Float_t ptLeptMC_;
   Float_t etaLeptMC_;
   Float_t phiLeptMC_;

   Float_t eNeuMC_;
   Float_t ptNeuMC_;
   Float_t etaNeuMC_;
   Float_t phiNeuMC_;


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

   Float_t trackCountingHighEffBJetTagJet1_[50];
   Float_t trackCountingHighEffBJetTagJet2_[50];

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
   Int_t motherPart_[20];


   Float_t energyPFMet_;
   Float_t phiPFMet_;
   Float_t pxPFMet_;
   Float_t pyPFMet_;

   Float_t uncorrEnergyAK5Jet_;//###
   Float_t SumEt_;

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
 TH1F* h1_Cont_inclusive;
 TH1F* h1_Cont_PV;
 TH1F* h1_Cont_TightMu;
 TH1F* h1_Cont_TightEle;
 TH1F* h1_Cont_VetoMU;
 TH1F* h1_Cont_VetoELE;
 TH1F* h1_Cont_JetsELE;
 TH1F* h1_Cont_JetsMU;

   TH1F *h1_resoPz;
   TH1F *h1_resoPt;

   bool DEBUG_VERBOSE_;

};




#endif
