// ------------------------------------------------------------
//  
//    Ntp1Finalizer_ZGamma - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"



class Ntp1Finalizer_ZGamma : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_ZGamma( const std::string& dataset );
  virtual ~Ntp1Finalizer_ZGamma() {};

  void finalize();



 private:

   std::string leptType_;

};

