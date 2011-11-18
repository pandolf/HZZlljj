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

   Bool_t passed_HLT_DoubleMu6_;
   Bool_t passed_HLT_DoubleMu7_;
   Bool_t passed_HLT_Mu13_Mu8_;
   Bool_t passed_HLT_Mu17_Mu8_;
   Bool_t passed_HLT_IsoMu17_;
   Bool_t passed_HLT_IsoMu24_;
   Bool_t passed_HLT_IsoMu24_eta2p1_;
   Bool_t passed_HLT_Mu8_Jet40_;
   Bool_t passed_HLT_L2DoubleMu23_NoVertex_;
   Bool_t passed_HLT_L2DoubleMu30_NoVertex_;
   Bool_t passed_HLT_TripleMu5_;

   Bool_t passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_;
   Bool_t passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_;


   int leptType_; //0: muon; 1: electron
   int leptTypeMC_; //0: muon; 1: electron

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

   Float_t trackCountingHighEffBJetTagJet1_[50];
   Float_t trackCountingHighPurBJetTagJet1_[50];
   Float_t simpleSecondaryVertexHighEffBJetTagJet1_[50];
   Float_t simpleSecondaryVertexHighPurBJetTagJet1_[50];
   Float_t jetBProbabilityBJetTagJet1_[50];
   Float_t jetProbabilityBJetTagJet1_[50];

   Float_t SFTCHEJet1_[50];
   Float_t SFerrTCHEJet1_[50];

   Float_t  ptGenJet1_[50];
   Float_t   eGenJet1_[50];
   Float_t phiGenJet1_[50];
   Float_t etaGenJet1_[50];

   Float_t  ptPartJet1_[50];
   Float_t   ePartJet1_[50];
   Float_t phiPartJet1_[50];
   Float_t etaPartJet1_[50];
   Int_t pdgIdPartJet1_[50];

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

   Float_t trackCountingHighEffBJetTagJet2_[50];
   Float_t trackCountingHighPurBJetTagJet2_[50];
   Float_t simpleSecondaryVertexHighEffBJetTagJet2_[50];
   Float_t simpleSecondaryVertexHighPurBJetTagJet2_[50];
   Float_t jetBProbabilityBJetTagJet2_[50];
   Float_t jetProbabilityBJetTagJet2_[50];

   Float_t SFTCHEJet2_[50];
   Float_t SFerrTCHEJet2_[50];

   Float_t  ptGenJet2_[50];
   Float_t   eGenJet2_[50];
   Float_t phiGenJet2_[50];
   Float_t etaGenJet2_[50];

   Float_t  ptPartJet2_[50];
   Float_t   ePartJet2_[50];
   Float_t phiPartJet2_[50];
   Float_t etaPartJet2_[50];
   Int_t pdgIdPartJet2_[50];


   Int_t  iJetRecoil_[50];
   Float_t  ptJetRecoil_[50];
   Float_t   eJetRecoil_[50];
   Float_t phiJetRecoil_[50];
   Float_t etaJetRecoil_[50];

   Float_t ptDJetRecoil_[50];
   Float_t rmsCandJetRecoil_[50];
   Int_t nChargedJetRecoil_[50];
   Int_t nNeutralJetRecoil_[50];
   Float_t QGlikelihoodJetRecoil_[50];

   Float_t  eChargedHadronsJetRecoil_[50];
   Float_t  ePhotonsJetRecoil_[50];
   Float_t  eNeutralEmJetRecoil_[50];
   Float_t  eNeutralHadronsJetRecoil_[50];
   Float_t  eMuonsJetRecoil_[50];
   Float_t  eElectronsJetRecoil_[50];
   Float_t  eHFHadronsJetRecoil_[50];
   Float_t  eHFEMJetRecoil_[50];

   Int_t  nChargedHadronsJetRecoil_[50];
   Int_t  nPhotonsJetRecoil_[50];
   Int_t  nNeutralHadronsJetRecoil_[50];
   Int_t  nMuonsJetRecoil_[50];
   Int_t  nElectronsJetRecoil_[50];
   Int_t  nHFHadronsJetRecoil_[50];
   Int_t  nHFEMJetRecoil_[50];

   Float_t trackCountingHighEffBJetTagJetRecoil_[50];
   Float_t trackCountingHighPurBJetTagJetRecoil_[50];
   Float_t simpleSecondaryVertexHighEffBJetTagJetRecoil_[50];
   Float_t simpleSecondaryVertexHighPurBJetTagJetRecoil_[50];
   Float_t jetBProbabilityBJetTagJetRecoil_[50];
   Float_t jetProbabilityBJetTagJetRecoil_[50];

   Float_t SFTCHEJetRecoil_[50];
   Float_t SFerrTCHEJetRecoil_[50];

   Float_t  ptGenJetRecoil_[50];
   Float_t   eGenJetRecoil_[50];
   Float_t phiGenJetRecoil_[50];
   Float_t etaGenJetRecoil_[50];

   Float_t  ptPartJetRecoil_[50];
   Float_t   ePartJetRecoil_[50];
   Float_t phiPartJetRecoil_[50];
   Float_t etaPartJetRecoil_[50];
   Int_t pdgIdPartJetRecoil_[50];


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


   Float_t epfMet_;
   Float_t sumEtpfMet_;
   Float_t metSignificance_;
   Float_t mEtSig_;
   Float_t phipfMet_;

   TH1D* h1_nCounter_Zee_;
   TH1D* h1_nCounter_Zmumu_;

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
