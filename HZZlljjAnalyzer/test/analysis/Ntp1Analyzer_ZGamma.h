//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_ZGamma_h
#define Ntp1Analyzer_ZGamma_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_ZGamma : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_ZGamma( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_ZGamma();

   virtual void CreateOutputFile();
   virtual void Loop();



 private:

   Int_t leptType_;

   Float_t eZGamma_;
   Float_t ptZGamma_;
   Float_t etaZGamma_;
   Float_t phiZGamma_;

   Float_t eZGammaGen_;
   Float_t ptZGammaGen_;
   Float_t etaZGammaGen_;
   Float_t phiZGammaGen_;

   Int_t nJets_total_;

   Int_t nJet_;
   Float_t eJet_[10];
   Float_t ptJet_[10];
   Float_t etaJet_[10];
   Float_t phiJet_[10];

   Int_t nJetsGen_total_;

   Int_t nJetGen_;
   Float_t eJetGen_[10];
   Float_t ptJetGen_[10];
   Float_t etaJetGen_[10];
   Float_t phiJetGen_[10];


   bool DEBUG_VERBOSE_;

};




#endif
