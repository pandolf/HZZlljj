//------------------------------------------------------------------

//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_WZ_h
#define Ntp1Analyzer_WZ_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_WZ : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_WZ( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_WZ();

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


   Int_t nJet_;

   Int_t  iJet_[50];
   Float_t  ptJet_[50];
   Float_t   eJet_[50];
   Float_t phiJet_[50];
   Float_t etaJet_[50];

   Float_t ptDJet_[50];
   Float_t rmsCandJet_[50];
   Int_t nChargedJet_[50];
   Int_t nNeutralJet_[50];
   Float_t QGlikelihoodJet_[50];

   Float_t  eChargedHadronsJet_[50];
   Float_t  ePhotonsJet_[50];
   Float_t  eNeutralEmJet_[50];
   Float_t  eNeutralHadronsJet_[50];
   Float_t  eMuonsJet_[50];
   Float_t  eElectronsJet_[50];
   Float_t  eHFHadronsJet_[50];
   Float_t  eHFEMJet_[50];

   Int_t  nChargedHadronsJet_[50];
   Int_t  nPhotonsJet_[50];
   Int_t  nNeutralHadronsJet_[50];
   Int_t  nMuonsJet_[50];
   Int_t  nElectronsJet_[50];
   Int_t  nHFHadronsJet_[50];
   Int_t  nHFEMJet_[50];

   Int_t nPart_;
   Float_t  ptPart_[20];
   Float_t   ePart_[20];
   Float_t phiPart_[20];
   Float_t etaPart_[20];
   Int_t pdgIdPart_[20];


   Float_t epfMet_;
   Float_t phipfMet_;


   bool DEBUG_VERBOSE_;

};




#endif
