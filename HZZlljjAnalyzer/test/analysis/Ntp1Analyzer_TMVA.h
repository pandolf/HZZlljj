//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_TMVA_h
#define Ntp1Analyzer_TMVA_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_TMVA : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_TMVA( const std::string& dataset, const std::string& jetChoice="BESTZ", const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_TMVA();

   virtual void CreateOutputFile();
   virtual void Loop();



 private:

   Int_t leptType_;

   Float_t epfMet_;

   Float_t ptLept1_;
   Float_t absetaLept1_;

   Float_t ptLept2_;
   Float_t absetaLept2_;

   Float_t  ptJet1_;
   Float_t absetaJet1_;

   Float_t  ptJet2_;
   Float_t absetaJet2_;

   Float_t ptJetRecoil_;
   Float_t absetaJetRecoil_;
   Float_t deltaR_recoil_jet1_;
   Float_t deltaR_recoil_Zjj_;
   Float_t deltaR_recoil_Higgs_;

   Float_t mZjj_;
   Float_t mZll_;
   Float_t deltaRjj_;
   Float_t deltaRll_;

   Float_t ptZjj_;
   Float_t ptZll_;

   Float_t ptZZ_;
   Float_t mZZ_;
   Float_t absEtaZZ_;
   Float_t deltaRZZ_;
   Float_t deltaAbsEtaZZ_;
   Float_t absDeltaEtaZZ_;
   Float_t absDeltaPhiZZ_;

   TH1F* h1_mZjj; 
   TH1F* h1_mZjj_matched; 
   TH1F* h1_mZjj_bestZ; 
   TH1F* h1_mZjj_bestH; 
   TH1F* h1_mZjj_closestPair; 


   std::string jetChoice_;

   bool DEBUG_VERBOSE_;

};




#endif