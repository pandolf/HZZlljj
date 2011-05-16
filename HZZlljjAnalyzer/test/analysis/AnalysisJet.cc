#include "AnalysisJet.h"



bool AnalysisJet::btag_loose() {

  return trackCountingHighEffBJetTag>1.83;

}


bool AnalysisJet::btag_medium() {

  return trackCountingHighEffBJetTag>4.;

}

