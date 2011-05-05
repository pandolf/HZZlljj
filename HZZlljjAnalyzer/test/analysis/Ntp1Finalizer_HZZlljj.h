// ------------------------------------------------------------
//  
//    Ntp1Finalizer_HZZlljj - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"



class Ntp1Finalizer_HZZlljj : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_HZZlljj( const std::string& dataset, const std::string& selectionType, const std::string& leptType="ALL" );
  virtual ~Ntp1Finalizer_HZZlljj() {};

  virtual void finalize();
  void setSelectionType( const std::string& selectionType );



 private:

   std::string leptType_;
   std::string selectionType_;

   float ptLept1_thresh_;
   float ptLept2_thresh_;
   float etaLept1_thresh_;
   float etaLept2_thresh_;
   float ptJet1_thresh_;
   float ptJet2_thresh_;
   float etaJet1_thresh_;
   float etaJet2_thresh_;
   float mZll_threshLo_;
   float mZll_threshHi_;
   float mZjj_threshLo_;
   float mZjj_threshHi_;
   float deltaRll_thresh_;
   float deltaRjj_thresh_;
   float ptZll_thresh_;
   float ptZjj_thresh_;
   float QGLikelihoodProd_thresh_;
   float helicityLD_thresh_;
   float mZZ_threshLo_;
   float mZZ_threshHi_;
   int requiredBTags_;

};

