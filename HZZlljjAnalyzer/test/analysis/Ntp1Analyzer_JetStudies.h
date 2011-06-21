//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_JetStudies_h
#define Ntp1Analyzer_JetStudies_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_JetStudies : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_JetStudies( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_JetStudies();

   virtual void CreateOutputFile();
   virtual void Loop();



 private:

   Int_t nvertex;
   Float_t rhoPF;

   Int_t nJet_;
   Float_t eJet_[5];
   Float_t ptJet_[5];
   Float_t etaJet_[5];
   Float_t phiJet_[5];

   Float_t ptDJet_[5];
   Float_t rmsCandJet_[5];
   Int_t nChargedJet_[5];
   Int_t nNeutralJet_[5];

   bool DEBUG_VERBOSE_;

};




#endif
