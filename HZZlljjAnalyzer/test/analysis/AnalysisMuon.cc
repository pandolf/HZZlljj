#include "AnalysisMuon.h"



bool AnalysisMuon::passedMuonID() {

  bool passed = ( isGlobalMuonPromptTight && isAllTrackerMuon 
                && pixelHits>0 && trackerHits>10 
                && nMatchedStations>=2
                && (dxy<0.02) && (dz<1.) );

  return passed;

}



bool AnalysisMuon::isIsolated() {

  bool isIsolated = this->combinedIsoRel() < 0.15;

  return isIsolated;

}




bool AnalysisMuon::passedVBTF() {

  bool isIsolated = this->isIsolated();
  bool passedMuonID = this->passedMuonID();

  return ( isIsolated && passedMuonID );

}


float AnalysisMuon::combinedIsoRel() {

 return (sumPt03 + emEt03 + hadEt03)/this->Pt();

}
