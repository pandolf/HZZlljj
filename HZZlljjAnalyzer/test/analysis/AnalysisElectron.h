// ---------------------------------------------------------------
//
//  AnalysisElectron - electron class used in the HZZlljj analysis
//
// ---------------------------------------------------------------


#ifndef AnalysisElectron_h
#define AnalysisElectron_h

#include "AnalysisLepton.h"




class AnalysisElectron : public AnalysisLepton {

 public:

  AnalysisElectron( float x=0., float y=0., float z=0., float t=0.) : AnalysisLepton( x, y, z, t ) {};
  AnalysisElectron( double x=0., double y=0., double z=0., double t=0.) : AnalysisLepton( x, y, z, t ) {};

  AnalysisElectron( const TLorentzVector &v) : AnalysisLepton( v ) {};

  bool isIsolatedVBTF80();
  bool isIsolatedVBTF95();

  bool electronIDVBTF80();
  bool electronIDVBTF95();

  bool conversionRejectionVBTF80();
  bool conversionRejectionVBTF95();

  bool passedVBTF80();
  bool passedVBTF95();

  double combinedIsoRel();




  // public data members:

  // isolation:
  double dr03TkSumPt;
  double dr03EcalRecHitSumEt;
  double dr03HcalTowerSumEt;

  // electron ID:
  double sigmaIetaIeta; 
  double deltaPhiAtVtx; 
  double deltaEtaAtVtx; 
  double hOverE; 

  // conversion rejection:
  int expInnerLayersGsfTrack;
  double convDist;
  double convDcot;

};

#endif
