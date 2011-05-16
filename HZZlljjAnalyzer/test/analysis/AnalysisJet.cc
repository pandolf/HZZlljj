#include "AnalysisJet.h"



bool AnalysisJet::btag_loose() const {

  return trackCountingHighEffBJetTag>1.83;

}


bool AnalysisJet::btag_medium() const {

  return trackCountingHighEffBJetTag>4.;

}

