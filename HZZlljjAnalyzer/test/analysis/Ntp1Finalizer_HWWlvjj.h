// ------------------------------------------------------------
//  
//    Ntp1Finalizer_HWWlvjj - Derived class 
//    for the finalization of the H->WW->lvjj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"



class Ntp1Finalizer_HWWlvjj : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_HWWlvjj( const std::string& dataset, const std::string& selectionType, const std::string& leptType="ALL" );
  virtual ~Ntp1Finalizer_HWWlvjj() {};

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
   float mtWll_threshLo_;
   float mtWll_threshHi_;
   float mWjj_threshLo_;
   float mWjj_threshHi_;
   float deltaRll_thresh_;
   float deltaRjj_thresh_;
   float ptWll_thresh_;
   float ptWjj_thresh_;
   float QGLikelihoodProd_thresh_;
   float helicityLD_thresh_;
   float mWW_threshLo_;
   float mWW_threshHi_;

};

