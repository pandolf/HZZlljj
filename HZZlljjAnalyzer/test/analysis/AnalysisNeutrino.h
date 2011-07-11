// ---------------------------------------------------------------
//
//  AnalysisNeutrino - neutrino class used in the HWWlnujj analysis
//
// ---------------------------------------------------------------


#ifndef AnalysisNeutrino_h
#define AnalysisNeutrino_h

#include "AnalysisLepton.h"




class AnalysisNeutrino : public AnalysisLepton {
  
 public:
  
  AnalysisNeutrino( float x=0., float y=0., float z=0., float t=0.) : AnalysisLepton( x, y, z, t ) {};
    AnalysisNeutrino( double x=0., double y=0., double z=0., double t=0.) : AnalysisLepton( x, y, z, t ) {};
      
      AnalysisNeutrino( const TLorentzVector &v) : AnalysisLepton( v ) {};

};

#endif
