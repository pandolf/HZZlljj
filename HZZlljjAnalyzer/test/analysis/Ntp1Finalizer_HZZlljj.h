// ------------------------------------------------------------
//  
//    Ntp1Finalizer_HZZlljj - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"



class Ntp1Finalizer_HZZlljj : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_HZZlljj( const std::string& dataset, const std::string& leptType="ALL" );
  virtual ~Ntp1Finalizer_HZZlljj() {};

  void finalize();



 private:

   std::string leptType_;

};

