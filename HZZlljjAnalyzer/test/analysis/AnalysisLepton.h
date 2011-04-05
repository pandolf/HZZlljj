// -------------------------------------------------------------
//
//  AnalysisLepton - lepton class used in the HZZlljj analysis
//
// -------------------------------------------------------------

#ifndef AnalysisLepton_h
#define AnalysisLepton_h

#include "TLorentzVector.h"



class AnalysisLepton : public TLorentzVector {

 public:

  AnalysisLepton( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {
    charge=0;
  }

  AnalysisLepton( const TLorentzVector &v) : TLorentzVector( v ) {
    charge=0;
  }

  int charge;

};

#endif
