// ------------------------------------------------------------
//  
//    Ntp1Finalizer_WZ - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"



class Ntp1Finalizer_WZ : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_WZ( const std::string& dataset, const std::string& selectionType, const std::string& leptType="ALL" );
  virtual ~Ntp1Finalizer_WZ() {};

  virtual void finalize();
  void setSelectionType( const std::string& selectionType );



 private:

   std::string leptType_;
   std::string selectionType_;

   float ptLept1_thresh_;
   float ptLept2_thresh_;
   float etaLept1_thresh_;
   float etaLept2_thresh_;
   float ptJet_thresh_;
   float etaJet_thresh_;
   float mZll_threshLo_;
   float mZll_threshHi_;
   float etaZll_thresh_;

};

