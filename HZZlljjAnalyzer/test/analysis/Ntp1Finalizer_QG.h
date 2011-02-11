// ------------------------------------------------------------
//  
//    Ntp1Finalizer_QG - Derived class 
//    for the computation of the basic quark-gluon
//    discrimination variables.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"



class Ntp1Finalizer_QG : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_QG( const std::string& dataset );
  virtual ~Ntp1Finalizer_QG() {};

  void finalize();



 private:

};

