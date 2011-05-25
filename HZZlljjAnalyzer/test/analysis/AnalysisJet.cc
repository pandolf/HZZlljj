#include "AnalysisJet.h"



bool AnalysisJet::btag_loose() const {

  return trackCountingHighEffBJetTag>1.7;

}


bool AnalysisJet::btag_medium() const {

  return trackCountingHighEffBJetTag>3.3;

}

