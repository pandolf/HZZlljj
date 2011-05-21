#include "AnalysisElectron.h"


bool AnalysisElectron::isIsolatedVBTF80() {

  double combinedIsoRel = this->combinedIsoRel();

  double combinedIsoRelThresh;
  if( fabs(this->Eta())<1.4442 )
    combinedIsoRelThresh = 0.07;
  else
    combinedIsoRelThresh = 0.06;


  bool isIsolated = (combinedIsoRel < combinedIsoRelThresh);

  return isIsolated;

}


bool AnalysisElectron::isIsolatedVBTF95() {

  double combinedIsoRel = this->combinedIsoRel();

  double combinedIsoRelThresh;
  if( fabs(this->Eta())<1.4442 )
    combinedIsoRelThresh = 0.15;
  else
    combinedIsoRelThresh = 0.1;

  bool isIsolated = (combinedIsoRel < combinedIsoRelThresh);

  return isIsolated;

}

bool AnalysisElectron::electronIDVBTF80() {

  double sigmaIetaIeta_thresh80;
  double deltaPhiAtVtx_thresh80;
  double deltaEtaAtVtx_thresh80;
  double hOverE_thresh80;

  if( fabs(this->Eta())<1.4442 ) {
    sigmaIetaIeta_thresh80 = 0.01;
    deltaPhiAtVtx_thresh80 = 0.06;
    deltaEtaAtVtx_thresh80 = 0.004;
    hOverE_thresh80 = 0.04;
  } else {
    sigmaIetaIeta_thresh80 = 0.03;
    deltaPhiAtVtx_thresh80 = 0.7;
    deltaEtaAtVtx_thresh80 = 0.007;
    hOverE_thresh80 = 0.15;
  }


  bool eleID_VBTF80 = (sigmaIetaIeta < sigmaIetaIeta_thresh80) &&
                      (fabs(deltaPhiAtVtx) < deltaPhiAtVtx_thresh80) &&
                      (fabs(deltaEtaAtVtx) < deltaEtaAtVtx_thresh80) &&
                      (hOverE < hOverE_thresh80);

  return eleID_VBTF80;

}

bool AnalysisElectron::electronIDVBTF95() {

  double sigmaIetaIeta_thresh95;
  double deltaPhiAtVtx_thresh95;
  double deltaEtaAtVtx_thresh95;
  double hOverE_thresh95;

  if( fabs(this->Eta())<1.4442 ) {
    sigmaIetaIeta_thresh95 = 0.01;
    deltaPhiAtVtx_thresh95 = 0.8;
    deltaEtaAtVtx_thresh95 = 0.007;
    hOverE_thresh95 = 0.15;
  } else {
    sigmaIetaIeta_thresh95 = 0.03;
    deltaPhiAtVtx_thresh95 = 0.7;
    deltaEtaAtVtx_thresh95 = 0.01; 
    hOverE_thresh95 = 0.07;
  }


  bool eleID_VBTF95 = (sigmaIetaIeta < sigmaIetaIeta_thresh95) &&
                      (fabs(deltaPhiAtVtx) < deltaPhiAtVtx_thresh95) &&
                      (fabs(deltaEtaAtVtx) < deltaEtaAtVtx_thresh95) &&
                      (hOverE < hOverE_thresh95);

  return eleID_VBTF95;

}


bool AnalysisElectron::conversionRejectionVBTF80() {

   bool convRej_VBTF80 = (expInnerLayersGsfTrack<=0) && (fabs(convDist)>0.02 || fabs(convDcot)>0.02);
   
   return convRej_VBTF80;

}

bool AnalysisElectron::conversionRejectionVBTF95() {

   bool convRej_VBTF95 = (expInnerLayersGsfTrack<=1) && (fabs(convDist)>-1. || fabs(convDcot)>-1.);
   
   return convRej_VBTF95;

}

bool AnalysisElectron::passedVBTF80() {

  bool isIsolated = this->isIsolatedVBTF80();
  bool electronID = this->electronIDVBTF80();
  bool conversionRejection = this->conversionRejectionVBTF80();

  bool passed = ( isIsolated && electronID && conversionRejection );

  return passed;

}

bool AnalysisElectron::passedVBTF95() {

  bool isIsolated = this->isIsolatedVBTF95();
  bool electronID = this->electronIDVBTF95();
  bool conversionRejection = this->conversionRejectionVBTF95();

  bool passed = ( isIsolated && electronID && conversionRejection );

  return passed;

}



double AnalysisElectron::combinedIsoRel() {

  double combinedIsoRel;

  if( fabs(this->Eta())<1.4442 )
    combinedIsoRel = ( dr03TkSumPt + TMath::Max(0., dr03EcalRecHitSumEt - 1.) + dr03HcalTowerSumEt ) / this->Pt();
  else
    combinedIsoRel = ( dr03TkSumPt + dr03EcalRecHitSumEt + dr03HcalTowerSumEt ) / this->Pt();

  return combinedIsoRel;

}
