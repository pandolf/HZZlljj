#include "AnalysisMuon.h"



bool AnalysisMuon::passedMuonID() {

  bool passed = ( isGlobalMuonPromptTight && isAllTrackerMuon && (dxy<0.02) && (dz<1.) );

  return passed;

}



bool AnalysisMuon::isIsolated() {

  bool isIsolated = (sumPt03 + emEt03 + hadEt03) < 0.15*this->Pt();

  return isIsolated;

}




bool AnalysisMuon::passedVBTF() {

  bool isIsolated = this->isIsolated();
  bool passedMuonID = this->passedMuonID();

  return ( isIsolated && passedMuonID );

}
