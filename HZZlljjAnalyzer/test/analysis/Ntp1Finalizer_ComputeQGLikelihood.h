// ------------------------------------------------------------
//  
//    Ntp1Finalizer_ComputeQGLikelihood - Derived class 
//    for the computation of the basic quark-gluon
//    discrimination variables.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"



class Ntp1Finalizer_ComputeQGLikelihood : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_ComputeQGLikelihood( const std::string& dataset );
  virtual ~Ntp1Finalizer_ComputeQGLikelihood() {};

  virtual void addFile(const std::string& dataset);

  virtual void finalize();



 private:

};

