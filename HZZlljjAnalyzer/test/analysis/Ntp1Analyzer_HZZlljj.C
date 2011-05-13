#include "Ntp1Analyzer_HZZlljj.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"
#include "TMVA/Reader.h"

#include "AnalysisElectron.h"
#include "AnalysisMuon.h"


//#include "fitTools.h"


double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);


class AnalysisJet : public TLorentzVector {

 public:

  AnalysisJet( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {
    eChargedHadrons=0.;
    ePhotons=0.;
    eNeutralEm=0.;
    eNeutralHadrons=0.;
    eElectrons=0.;
    nChargedHadrons=0;
    nPhotons=0;
    nNeutralHadrons=0;
  }

  float eChargedHadrons;
  float ePhotons;
  float eNeutralEm;
  float eNeutralHadrons;
  float eMuons;
  float eElectrons;
//float eHFHadrons;
//float eHFEM;

  int nChargedHadrons;
  int nPhotons;
  int nNeutralHadrons;
  int nMuons;
  int nElectrons;
//int nHFHadrons;
//int nHFEM;

  float ptD;
  float rmsCand;
  int nCharged;
  int nNeutral;
  float QGlikelihood;

  float ptGen;
  float etaGen;
  float phiGen;
  float eGen;

  //btags:
  float trackCountingHighEffBJetTag;
  float trackCountingHighPurBJetTag;
  float simpleSecondaryVertexHighEffBJetTag;
  float simpleSecondaryVertexHighPurBJetTag;
  float jetBProbabilityBJetTag;
  float jetProbabilityBJetTag;

};



/*
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

*/


Ntp1Analyzer_HZZlljj::Ntp1Analyzer_HZZlljj( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "HZZlljj", dataset, flags, tree ) {


  h1_nCounter_Zee_ = new TH1D("nCounter_Zee", "", 1, 0., 1.);
  h1_nCounter_Zmumu_ = new TH1D("nCounter_Zmumu", "", 1, 0., 1.);


} //constructor



void Ntp1Analyzer_HZZlljj::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");
  reducedTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");
  reducedTree_->Branch("leptTypeMC",&leptTypeMC_,"leptTypeMC_/I");

//// triggers:
//reducedTree_->Branch("HLT_Mu11", &HLT_Mu11_, "HLT_Mu11_/O");
//reducedTree_->Branch("HLT_Ele17_SW_EleId_L1R", &HLT_Ele17_SW_EleId_L1R_, "HLT_Ele17_SW_EleId_L1R_/O");
//reducedTree_->Branch("HLT_DoubleMu3", &HLT_DoubleMu3_, "HLT_DoubleMu3_/O");


  reducedTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  reducedTree_->Branch("leptType",  &leptType_,  "leptType_/I");
  
  reducedTree_->Branch("eZqqMC",  &eZqqMC_,  "eZqqMC_/F");
  reducedTree_->Branch("ptZqqMC",  &ptZqqMC_,  "ptZqqMC_/F");
  reducedTree_->Branch("etaZqqMC",  &etaZqqMC_,  "etaZqqMC_/F");
  reducedTree_->Branch("phiZqqMC",  &phiZqqMC_,  "phiZqqMC_/F");

  reducedTree_->Branch("eZllMC",  &eZllMC_,  "eZllMC_/F");
  reducedTree_->Branch("ptZllMC",  &ptZllMC_,  "ptZllMC_/F");
  reducedTree_->Branch("etaZllMC",  &etaZllMC_,  "etaZllMC_/F");
  reducedTree_->Branch("phiZllMC",  &phiZllMC_,  "phiZllMC_/F");

  reducedTree_->Branch("eHiggsMC",  &eHiggsMC_,  "eHiggsMC_/F");
  reducedTree_->Branch("ptHiggsMC",  &ptHiggsMC_,  "ptHiggsMC_/F");
  reducedTree_->Branch("etaHiggsMC",  &etaHiggsMC_,  "etaHiggsMC_/F");
  reducedTree_->Branch("phiHiggsMC",  &phiHiggsMC_,  "phiHiggsMC_/F");

  reducedTree_->Branch("eLept1",  &eLept1_,  "eLept1_/F");
  reducedTree_->Branch("ptLept1",  &ptLept1_,  "ptLept1_/F");
  reducedTree_->Branch("etaLept1",  &etaLept1_,  "etaLept1_/F");
  reducedTree_->Branch("phiLept1",  &phiLept1_,  "phiLept1_/F");
  reducedTree_->Branch("chargeLept1",  &chargeLept1_,  "chargeLept1_/I");

  reducedTree_->Branch("eLept1Gen",  &eLept1Gen_,  "eLept1Gen_/F");
  reducedTree_->Branch("ptLept1Gen",  &ptLept1Gen_,  "ptLept1Gen_/F");
  reducedTree_->Branch("etaLept1Gen",  &etaLept1Gen_,  "etaLept1Gen_/F");
  reducedTree_->Branch("phiLept1Gen",  &phiLept1Gen_,  "phiLept1Gen_/F");

  reducedTree_->Branch("eLept2",  &eLept2_,  "eLept2_/F");
  reducedTree_->Branch("ptLept2",  &ptLept2_,  "ptLept2_/F");
  reducedTree_->Branch("etaLept2",  &etaLept2_,  "etaLept2_/F");
  reducedTree_->Branch("phiLept2",  &phiLept2_,  "phiLept2_/F");
  reducedTree_->Branch("chargeLept2",  &chargeLept2_,  "chargeLept2_/I");

  reducedTree_->Branch("eLept2Gen",  &eLept2Gen_,  "eLept2Gen_/F");
  reducedTree_->Branch("ptLept2Gen",  &ptLept2Gen_,  "ptLept2Gen_/F");
  reducedTree_->Branch("etaLept2Gen",  &etaLept2Gen_,  "etaLept2Gen_/F");
  reducedTree_->Branch("phiLept2Gen",  &phiLept2Gen_,  "phiLept2Gen_/F");

  reducedTree_->Branch("nPairs", &nPairs_, "nPairs_/I");

  reducedTree_->Branch("iJet1",  iJet1_,  "iJet1_[nPairs_]/I");
  reducedTree_->Branch("eJet1",  eJet1_,  "eJet1_[nPairs_]/F");
  reducedTree_->Branch( "ptJet1",  ptJet1_,  "ptJet1_[nPairs_]/F");
  reducedTree_->Branch("etaJet1", etaJet1_, "etaJet1_[nPairs_]/F");
  reducedTree_->Branch("phiJet1", phiJet1_, "phiJet1_[nPairs_]/F");

  reducedTree_->Branch("ptDJet1", ptDJet1_, "ptDJet1_[nPairs_]/F");
  reducedTree_->Branch("rmsCandJet1", rmsCandJet1_, "rmsCandJet1_[nPairs_]/F");
  reducedTree_->Branch("nChargedJet1", nChargedJet1_, "nChargedJet1_[nPairs_]/F");
  reducedTree_->Branch("nNeutralJet1", nNeutralJet1_, "nNeutralJet1_[nPairs_]/F");
  reducedTree_->Branch("QGlikelihoodJet1", QGlikelihoodJet1_, "QGlikelihoodJet1_[nPairs_]/F");

  reducedTree_->Branch("eChargedHadronsJet1", eChargedHadronsJet1_, "eChargedHadronsJet1_[nPairs_]/F");
  reducedTree_->Branch("ePhotonsJet1", ePhotonsJet1_, "ePhotonsJet1_[nPairs_]/F");
  reducedTree_->Branch("eNeutralEmJet1", eNeutralEmJet1_, "eNeutralEmJet1_[nPairs_]/F");
  reducedTree_->Branch("eNeutralHadronsJet1", eNeutralHadronsJet1_, "eNeutralHadronsJet1_[nPairs_]/F");
  reducedTree_->Branch("eMuonsJet1", eMuonsJet1_, "eMuonsJet1_[nPairs_]/F");
  reducedTree_->Branch("eElectronsJet1", eElectronsJet1_, "eElectronsJet1_[nPairs_]/F");
  reducedTree_->Branch("eHFHadronsJet1", eHFHadronsJet1_, "eHFHadronsJet1_[nPairs_]/F");
  reducedTree_->Branch("eHFEMJet1", eHFEMJet1_, "eHFEMJet1_[nPairs_]/F");

  reducedTree_->Branch("nChargedHadronsJet1", nChargedHadronsJet1_, "nChargedHadronsJet1_[nPairs_]/I");
  reducedTree_->Branch("nPhotonsJet1", nPhotonsJet1_, "nPhotonsJet1_[nPairs_]/I");
  reducedTree_->Branch("nNeutralHadronsJet1", nNeutralHadronsJet1_, "nNeutralHadronsJet1_[nPairs_]/I");
  reducedTree_->Branch("nMuonsJet1", nMuonsJet1_, "nMuonsJet1_[nPairs_]/I");
  reducedTree_->Branch("nElectronsJet1", nElectronsJet1_, "nElectronsJet1_[nPairs_]/I");
  reducedTree_->Branch("nHFHadronsJet1", nHFHadronsJet1_, "nHFHadronsJet1_[nPairs_]/I");
  reducedTree_->Branch("nHFEMJet1", nHFEMJet1_, "nHFEMJet1_[nPairs_]/I");

  reducedTree_->Branch("trackCountingHighEffBJetTagJet1", trackCountingHighEffBJetTagJet1_, "trackCountingHighEffBJetTagJet1_[nPairs_]/F");
  reducedTree_->Branch("trackCountingHighPurBJetTagJet1", trackCountingHighPurBJetTagJet1_, "trackCountingHighPurBJetTagJet1_[nPairs_]/F");
  reducedTree_->Branch("simpleSecondaryVertexHighEffBJetTagJet1", simpleSecondaryVertexHighEffBJetTagJet1_, "simpleSecondaryVertexHighEffBJetTagJet1_[nPairs_]/F");
  reducedTree_->Branch("simpleSecondaryVertexHighPurBJetTagJet1", simpleSecondaryVertexHighPurBJetTagJet1_, "simpleSecondaryVertexHighPurBJetTagJet1_[nPairs_]/F");
  reducedTree_->Branch("jetBProbabilityBJetTagJet1", jetBProbabilityBJetTagJet1_, "jetBProbabilityBJetTagJet1_[nPairs_]/F");
  reducedTree_->Branch("jetProbabilityBJetTagJet1", jetProbabilityBJetTagJet1_, "jetProbabilityBJetTagJet1_[nPairs_]/F");

  reducedTree_->Branch("eJet1Gen",  eJet1Gen_,  "eJet1Gen_[nPairs_]/F");
  reducedTree_->Branch( "ptJet1Gen",  ptJet1Gen_,  "ptJet1Gen_[nPairs_]/F");
  reducedTree_->Branch("etaJet1Gen", etaJet1Gen_, "etaJet1Gen_[nPairs_]/F");
  reducedTree_->Branch("phiJet1Gen", phiJet1Gen_, "phiJet1Gen_[nPairs_]/F");

  reducedTree_->Branch("nPFCand1",  &nPFCand1_,  "nPFCand1_/I");
  reducedTree_->Branch("ePFCand1",  &ePFCand1_,  "ePFCand1_[nPFCand1_]/F");
  reducedTree_->Branch("ptPFCand1",  &ptPFCand1_,  "ptPFCand1_[nPFCand1_]/F");
  reducedTree_->Branch("etaPFCand1",  &etaPFCand1_,  "etaPFCand1_[nPFCand1_]/F");
  reducedTree_->Branch("phiPFCand1",  &phiPFCand1_,  "phiPFCand1_[nPFCand1_]/F");
  reducedTree_->Branch("particleTypePFCand1",  &particleTypePFCand1_,  "particleTypePFCand1_[nPFCand1_]/I");

  reducedTree_->Branch("iJet2",  iJet2_,  "iJet2_[nPairs_]/I");
  reducedTree_->Branch("eJet2",  eJet2_,  "eJet2_[nPairs_]/F");
  reducedTree_->Branch( "ptJet2",  ptJet2_,  "ptJet2_[nPairs_]/F");
  reducedTree_->Branch("etaJet2", etaJet2_, "etaJet2_[nPairs_]/F");
  reducedTree_->Branch("phiJet2", phiJet2_, "phiJet2_[nPairs_]/F");

  reducedTree_->Branch("ptDJet2", ptDJet2_, "ptDJet2_[nPairs_]/F");
  reducedTree_->Branch("rmsCandJet2", rmsCandJet2_, "rmsCandJet2_[nPairs_]/F");
  reducedTree_->Branch("nChargedJet2", nChargedJet2_, "nChargedJet2_[nPairs_]/F");
  reducedTree_->Branch("nNeutralJet2", nNeutralJet2_, "nNeutralJet2_[nPairs_]/F");
  reducedTree_->Branch("QGlikelihoodJet2", QGlikelihoodJet2_, "QGlikelihoodJet2_[nPairs_]/F");

  reducedTree_->Branch("eChargedHadronsJet2", eChargedHadronsJet2_, "eChargedHadronsJet2_[nPairs_]/F");
  reducedTree_->Branch("ePhotonsJet2", ePhotonsJet2_, "ePhotonsJet2_[nPairs_]/F");
  reducedTree_->Branch("eNeutralEmJet2", eNeutralEmJet2_, "eNeutralEmJet2_[nPairs_]/F");
  reducedTree_->Branch("eNeutralHadronsJet2", eNeutralHadronsJet2_, "eNeutralHadronsJet2_[nPairs_]/F");
  reducedTree_->Branch("eMuonsJet2", eMuonsJet2_, "eMuonsJet2_[nPairs_]/F");
  reducedTree_->Branch("eElectronsJet2", eElectronsJet2_, "eElectronsJet2_[nPairs_]/F");
  reducedTree_->Branch("eHFHadronsJet2", eHFHadronsJet2_, "eHFHadronsJet2_[nPairs_]/F");
  reducedTree_->Branch("eHFEMJet2", eHFEMJet2_, "eHFEMJet2_[nPairs_]/F");

  reducedTree_->Branch("nChargedHadronsJet2", nChargedHadronsJet2_, "nChargedHadronsJet2_[nPairs_]/I");
  reducedTree_->Branch("nPhotonsJet2", nPhotonsJet2_, "nPhotonsJet2_[nPairs_]/I");
  reducedTree_->Branch("nNeutralHadronsJet2", nNeutralHadronsJet2_, "nNeutralHadronsJet2_[nPairs_]/I");
  reducedTree_->Branch("nMuonsJet2", nMuonsJet2_, "nMuonsJet2_[nPairs_]/I");
  reducedTree_->Branch("nElectronsJet2", nElectronsJet2_, "nElectronsJet2_[nPairs_]/I");
  reducedTree_->Branch("nHFHadronsJet2", nHFHadronsJet2_, "nHFHadronsJet2_[nPairs_]/I");
  reducedTree_->Branch("nHFEMJet2", nHFEMJet2_, "nHFEMJet2_[nPairs_]/I");

  reducedTree_->Branch("trackCountingHighEffBJetTagJet2", trackCountingHighEffBJetTagJet2_, "trackCountingHighEffBJetTagJet2_[nPairs_]/F");
  reducedTree_->Branch("trackCountingHighPurBJetTagJet2", trackCountingHighPurBJetTagJet2_, "trackCountingHighPurBJetTagJet2_[nPairs_]/F");
  reducedTree_->Branch("simpleSecondaryVertexHighEffBJetTagJet2", simpleSecondaryVertexHighEffBJetTagJet2_, "simpleSecondaryVertexHighEffBJetTagJet2_[nPairs_]/F");
  reducedTree_->Branch("simpleSecondaryVertexHighPurBJetTagJet2", simpleSecondaryVertexHighPurBJetTagJet2_, "simpleSecondaryVertexHighPurBJetTagJet2_[nPairs_]/F");
  reducedTree_->Branch("jetBProbabilityBJetTagJet2", jetBProbabilityBJetTagJet2_, "jetBProbabilityBJetTagJet2_[nPairs_]/F");
  reducedTree_->Branch("jetProbabilityBJetTagJet2", jetProbabilityBJetTagJet2_, "jetProbabilityBJetTagJet2_[nPairs_]/F");

  reducedTree_->Branch("eJet2Gen",  eJet2Gen_,  "eJet2Gen_[nPairs_]/F");
  reducedTree_->Branch( "ptJet2Gen",  ptJet2Gen_,  "ptJet2Gen_[nPairs_]/F");
  reducedTree_->Branch("etaJet2Gen", etaJet2Gen_, "etaJet2Gen_[nPairs_]/F");
  reducedTree_->Branch("phiJet2Gen", phiJet2Gen_, "phiJet2Gen_[nPairs_]/F");


  reducedTree_->Branch("nPFCand2",  &nPFCand2_,  "nPFCand2_/I");
  reducedTree_->Branch("ePFCand2",  &ePFCand2_,  "ePFCand2_[nPFCand2_]/F");
  reducedTree_->Branch("ptPFCand2",  &ptPFCand2_,  "ptPFCand2_[nPFCand2_]/F");
  reducedTree_->Branch("etaPFCand2",  &etaPFCand2_,  "etaPFCand2_[nPFCand2_]/F");
  reducedTree_->Branch("phiPFCand2",  &phiPFCand2_,  "phiPFCand2_[nPFCand2_]/F");
  reducedTree_->Branch("particleTypePFCand2",  &particleTypePFCand2_,  "particleTypePFCand2_[nPFCand2_]/I");

  reducedTree_->Branch("nPart", &nPart_, "nPart_/I");
  reducedTree_->Branch("ePart",  ePart_,  "ePart_[nPart_]/F");
  reducedTree_->Branch( "ptPart",  ptPart_,  "ptPart_[nPart_]/F");
  reducedTree_->Branch("etaPart", etaPart_, "etaPart_[nPart_]/F");
  reducedTree_->Branch("phiPart", phiPart_, "phiPart_[nPart_]/F");
  reducedTree_->Branch("pdgIdPart", pdgIdPart_, "pdgIdPart_[nPart_]/I");
  reducedTree_->Branch("motherPart", motherPart_, "motherPart_[nPart_]/I");


  reducedTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  reducedTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");

  
  int nBins_eff = 20;
  float ptMin_eff = 10.;
  float ptMax_eff = 150.;
  h1_nEvents_vs_ptEle = new TH1F("nEvents_vs_ptEle", "", nBins_eff, ptMin_eff, ptMax_eff);
  h1_nEvents_vs_ptMuon = new TH1F("nEvents_vs_ptMuon", "", nBins_eff, ptMin_eff, ptMax_eff);
  h1_passed_vs_ptEle = new TH1F("passed_vs_ptEle", "", nBins_eff, ptMin_eff, ptMax_eff);
  h1_passed_vs_ptMuon = new TH1F("passed_vs_ptMuon", "", nBins_eff, ptMin_eff, ptMax_eff);
  h1_deltaRmatching_muons = new TH1F("deltaRmatching_muons", "", 100, 0., 0.01);
  h1_deltaRmatching_electrons = new TH1F("deltaRmatching_electrons", "", 100, 0., 0.01);
  h1_deltaRmatching_jet_parton = new TH1F("deltaRmatching_jet_parton", "", 100, 0., 0.6);
  h1_deltaRmatching_genjet_parton = new TH1F("deltaRmatching_genjet_parton", "", 100, 0., 0.6);
  h1_deltaRmatching_jet_genjet = new TH1F("deltaRmatching_jet_genjet", "", 100, 0., 0.6);
  h1_deltaRmatching_jet_leptonParton = new TH1F("deltaRmatching_leptonParton", "", 100, 0., 4.);
  h1_nJets30 = new TH1F("nJets30", "", 31, -0.5, 30.5);
//h1_indexMatchedJet = new TH1F("indexMatchedJet", "", 6, -0.5, 5.5);
//h1_indexMatched05Jet = new TH1F("indexMatched05Jet", "", 6, -0.5, 5.5);
//h1_nMatched_per_event = new TH1F("nMatched_per_event", "", 6, -0.5, 5.5);
//h1_nMatched05_per_event = new TH1F("nMatched05_per_event", "", 6, -0.5, 5.5);
//h1_pdgIdParton1 = new TH1F("pdgIdParton1", "", 36, -10.5, 25.5);
//h1_pdgIdParton2 = new TH1F("pdgIdParton2", "", 36, -10.5, 25.5);
//h1_ptHadronicZ = new TH1F("ptHadronicZ", "", 50, 0., 400.);
//h1_deltaRqq = new TH1F("deltaRqq", "", 50, 0., 3.);

} 



Ntp1Analyzer_HZZlljj::~Ntp1Analyzer_HZZlljj() {

  outfile_->cd();
  h1_nCounter_Zee_->Write();
  h1_nCounter_Zmumu_->Write();
  h1_nEvents_vs_ptEle->Write();
  h1_nEvents_vs_ptMuon->Write();
  h1_passed_vs_ptEle->Write();
  h1_passed_vs_ptMuon->Write();
  h1_deltaRmatching_muons->Write();
  h1_deltaRmatching_electrons->Write();
  h1_deltaRmatching_jet_parton->Write();
  h1_deltaRmatching_genjet_parton->Write();
  h1_deltaRmatching_jet_genjet->Write();
  h1_deltaRmatching_jet_leptonParton->Write();
  h1_nJets30->Write();
//h1_indexMatchedJet->Write();
//h1_indexMatched05Jet->Write();
//h1_nMatched_per_event->Write();
//h1_nMatched05_per_event->Write();
//h1_pdgIdParton1->Write();
//h1_pdgIdParton2->Write();
//h1_ptHadronicZ->Write();
//h1_deltaRqq->Write();
  

}



void Ntp1Analyzer_HZZlljj::Loop()
{


   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;


   // to fix Z->ee BR bug, compute number of events having Z->ee and Z->mumu
   int nCounterZee = (isMC_) ? fChain->GetEntries("statusMc==3 && idMc==11") : 0;
   int nCounterZmumu = (isMC_) ? fChain->GetEntries("statusMc==3 && idMc==13") : 0;
   h1_nCounter_Zee_->SetBinContent( 1, nCounterZee );
   h1_nCounter_Zmumu_->SetBinContent( 1, nCounterZmumu );


   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();


   Long64_t nbytes = 0, nb = 0;

   TRandom3 rand;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%100000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;



     run_ = runNumber;
     LS_ = lumiBlock;
     event_ = eventNumber;
     eventWeight_ = -1.; //default
     leptTypeMC_ = -1;


     if( !isGoodEvent() ) continue; //this takes care also of trigger

     if( nPV==0 ) continue;
     bool goodVertex = (ndofPV[0] >= 4.0 && sqrt(PVxPV[0]*PVxPV[0]+PVyPV[0]*PVyPV[0]) < 2. && fabs(PVzPV[0]) < 24. );
     if( !goodVertex ) continue;
  
     nvertex_ = nPV;
     rhoPF_ = rhoFastjet;


     //bool isMC = ( runNumber < 5 );


     ptHat_ = (isMC_) ? genPtHat : ptHat_;



     bool noLeptons = false;
     TLorentzVector lept1MC, lept2MC;
     int zIndexqq=-1;
     int zIndexll=-1;

     if( isMC_ ) {


       // first look for Z->qq
       std::vector<TLorentzVector> quarksMC;

       for( unsigned iMc=0; iMc<nMc && quarksMC.size()<2; ++iMc ) {

         // quarks have status 3
         if( statusMc[iMc] != 3 ) continue;

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

         if( fabs(idMc[iMc])<7 && idMc[mothMc[iMc]]==23 ) {
           zIndexqq = mothMc[iMc];
           quarksMC.push_back( *thisParticle );
         }

       }

       // (checked that always 2 quarks are found)
       if( quarksMC.size()==2 && zIndexqq!=-1 ) {

         TLorentzVector ZqqMC;
         ZqqMC.SetPtEtaPhiE( pMc[zIndexqq]*sin(thetaMc[zIndexqq]), etaMc[zIndexqq], phiMc[zIndexqq], energyMc[zIndexqq] );

         ptZqqMC_  = ZqqMC.Pt();
         eZqqMC_   = ZqqMC.Energy();
         etaZqqMC_ = ZqqMC.Eta();
         phiZqqMC_ = ZqqMC.Phi();

      // float ptZqq = pMc[zIndexqq]*sin(thetaMc[zIndexqq]);
      // h1_ptHadronicZ->Fill( ptZqq );

      // float deltaRqq = quarksMC[0].DeltaR(quarksMC[1]);
      // h1_deltaRqq->Fill(deltaRqq);

       }

       // now look for Z->ll

       std::vector<TLorentzVector> electronsMC;
       std::vector<TLorentzVector> muonsMC;

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         // partons only
         if( statusMc[iMc] != 3 ) continue;

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

       
         if( idMc[mothMc[iMc]]==23 ) {
           zIndexll = mothMc[iMc]; 
           if( fabs(idMc[iMc])==11 && idMc[mothMc[iMc]]==23 ) electronsMC.push_back( *thisParticle );
           if( fabs(idMc[iMc])==13 && idMc[mothMc[iMc]]==23 ) muonsMC.push_back( *thisParticle );
         }

         delete thisParticle;
         thisParticle = 0;

       }

       if( electronsMC.size()==2 ) {
         if( electronsMC[0].Pt() > electronsMC[1].Pt() ) {
           lept1MC = electronsMC[0];
           lept2MC = electronsMC[1];
         } else {
           lept1MC = electronsMC[1];
           lept2MC = electronsMC[0];
         }
         if( (fabs(lept1MC.Eta()) < 2.5) && ( fabs(lept1MC.Eta())<1.4442 || fabs(lept1MC.Eta())>1.566) ) h1_nEvents_vs_ptEle->Fill( lept1MC.Pt() );
         if( (fabs(lept2MC.Eta()) < 2.5) && ( fabs(lept2MC.Eta())<1.4442 || fabs(lept2MC.Eta())>1.566) ) h1_nEvents_vs_ptEle->Fill( lept2MC.Pt() );
       } else if( muonsMC.size()==2 ) {
         if( muonsMC[0].Pt() > muonsMC[1].Pt() ) {
           lept1MC = muonsMC[0];
           lept2MC = muonsMC[1];
         } else {
           lept1MC = muonsMC[1];
           lept2MC = muonsMC[0];
         }
         if( fabs(lept1MC.Eta()) < 2.4 ) h1_nEvents_vs_ptMuon->Fill( lept1MC.Pt() );
         if( fabs(lept1MC.Eta()) < 2.1 || fabs(lept2MC.Eta()) < 2.1 ) h1_nEvents_vs_ptMuon->Fill( lept2MC.Pt() );
       } else {
         //taus
         noLeptons = true;
       }



       if( !noLeptons ) {

         TLorentzVector ZllMC;
         ZllMC.SetPtEtaPhiE( pMc[zIndexll]*sin(thetaMc[zIndexll]), etaMc[zIndexll], phiMc[zIndexll], energyMc[zIndexll] );

         ptZllMC_  = ZllMC.Pt();
         eZllMC_   = ZllMC.Energy();
         etaZllMC_ = ZllMC.Eta();
         phiZllMC_ = ZllMC.Phi();

         if( muonsMC.size() > 0 ) leptTypeMC_ = 0;
         else if( electronsMC.size() > 0 ) leptTypeMC_ = 1;

       }


       // now look for the higgs:
       if( zIndexll!=-1 && zIndexqq!=-1 ) {

         int higgsIndex = mothMc[zIndexll];

         if( idMc[higgsIndex] == 25 ) {

           TLorentzVector HiggsMC;
           HiggsMC.SetPtEtaPhiE( pMc[higgsIndex]*sin(thetaMc[higgsIndex]), etaMc[higgsIndex], phiMc[higgsIndex], energyMc[higgsIndex] );

           eHiggsMC_   = HiggsMC.Energy(); 
           ptHiggsMC_  = HiggsMC.Pt(); 
           etaHiggsMC_ = HiggsMC.Eta(); 
           phiHiggsMC_ = HiggsMC.Phi(); 

         } // if higgs

       } //if found two Z's

     } //if isMC

     if( !noLeptons )
       if( lept1MC.Pt() < lept2MC.Pt() ) std::cout << "WARNING MC leptons not ordered in pt!!" << std::endl;



     // -----------------------------
     //      FROM NOW ON RECO
     // -----------------------------


     epfMet_ = energyPFMet[0];
     phipfMet_ = phiPFMet[0];


     // ------------------
     // MUONS
     // ------------------

     std::vector<AnalysisMuon> muons;
     int chargeFirstMuon;


     for( unsigned int iMuon=0; iMuon<nMuon && (muons.size()<2); ++iMuon ) {

       AnalysisMuon thisMuon( pxMuon[iMuon], pyMuon[iMuon], pzMuon[iMuon], energyMuon[iMuon] );
       thisMuon.charge = chargeMuon[iMuon];

       // --------------
       // kinematics:
       // --------------
       if( thisMuon.Pt() < 20. ) continue;

       thisMuon.isGlobalMuonPromptTight = (muonIdMuon[iMuon]>>8)&1;
       thisMuon.isAllTrackerMuon = (muonIdMuon[iMuon]>>11)&1;


       // to compute dxy, look for primary vertex:
       int hardestPV = -1;
       float sumPtMax = 0.0;
       for(int v=0; v<nPV; v++) {
         if(SumPtPV[v] > sumPtMax) {
           sumPtMax = SumPtPV[v];
           hardestPV = v;
         }
       }  
   
       float dxy;
       if( hardestPV==-1 ) {
         dxy = 0.;
       } else {
         dxy = fabs(trackDxyPV(PVxPV[hardestPV], PVyPV[hardestPV], PVzPV[hardestPV],
                              trackVxTrack[trackIndexMuon[iMuon]], trackVyTrack[trackIndexMuon[iMuon]], trackVzTrack[trackIndexMuon[iMuon]],
                              pxTrack[trackIndexMuon[iMuon]], pyTrack[trackIndexMuon[iMuon]], pzTrack[trackIndexMuon[iMuon]]));
       }


       float dz = fabs(trackVzTrack[trackIndexMuon[iMuon]]-PVzPV[hardestPV]);

       thisMuon.dxy = dxy;
       thisMuon.dz = dz;

       thisMuon.sumPt03 = sumPt03Muon[iMuon];
       thisMuon.emEt03  = emEt03Muon[iMuon];
       thisMuon.hadEt03 = hadEt03Muon[iMuon];

       if( !thisMuon.passedVBTF() ) continue;


       // for now simple selection, will have to optimize this (T&P?)
       if( muons.size()==0 ) {
         muons.push_back( thisMuon );
         chargeFirstMuon = chargeMuon[iMuon];
       } else {
         if( chargeMuon[iMuon]==chargeFirstMuon ) continue;
         //if( fabs(muons[0].Eta())>2.1 && fabs(thisMuon.Eta())>2.1 ) continue;
         muons.push_back(thisMuon);
       }

     } //for muons



     // ------------------
     // ELECTRONS
     // ------------------

     std::vector<AnalysisElectron> electrons;
     int chargeFirstEle = 0;
     bool firstPassedVBTF80 = false;

     for( unsigned int iEle=0; (iEle<nEle) && (electrons.size()<2); ++iEle ) {

       AnalysisElectron thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );
       thisEle.charge = chargeEle[iEle];


       // --------------
       // kinematics:
       // --------------
       if( thisEle.Pt() < 20. ) continue;
       if( (fabs(thisEle.Eta()) > 2.5) || ( fabs(thisEle.Eta())>1.4442 && fabs(thisEle.Eta())<1.566) ) continue;

       // isolation
       thisEle.dr03TkSumPt = dr03TkSumPtEle[iEle];
       thisEle.dr03EcalRecHitSumEt = dr03EcalRecHitSumEtEle[iEle];
       thisEle.dr03HcalTowerSumEt = dr03HcalTowerSumEtEle[iEle];

       // electron ID
       thisEle.sigmaIetaIeta = (superClusterIndexEle[iEle]>=0) ? covIEtaIEtaSC[superClusterIndexEle[iEle]] : covIEtaIEtaSC[PFsuperClusterIndexEle[iEle]];
       thisEle.deltaPhiAtVtx = deltaPhiAtVtxEle[iEle];
       thisEle.deltaEtaAtVtx = deltaEtaAtVtxEle[iEle];
       thisEle.hOverE = hOverEEle[iEle];

       // conversion rejection
       thisEle.expInnerLayersGsfTrack = expInnerLayersGsfTrack[gsfTrackIndexEle[iEle]];
       thisEle.convDist = convDistEle[iEle];
       thisEle.convDcot = convDcotEle[iEle];


       bool passed_VBTF95 = thisEle.passedVBTF95();
       bool passed_VBTF80 = thisEle.passedVBTF80();

       if( !passed_VBTF95 ) continue;

       // check that not matched to muon (clean electrons faked by muon MIP):
       bool matchedtomuon=false;
       for( std::vector<AnalysisMuon>::iterator iMu=muons.begin(); iMu!=muons.end(); ++iMu )
         if( iMu->DeltaR(thisEle)<0.1 ) matchedtomuon=true;

       if( matchedtomuon ) continue;


       // for now simple selection, will have to optimize this (T&P?)
       // one electron required to pass VBTF80, the other VBTF95
       if( electrons.size()==0 ) {
         electrons.push_back( thisEle );
         chargeFirstEle = chargeEle[iEle];
         if( passed_VBTF80 ) firstPassedVBTF80 = true;
       } else if( chargeEle[iEle] != chargeFirstEle && ( firstPassedVBTF80||passed_VBTF80 ) ) {
         electrons.push_back( thisEle );
       }


     } //for electrons


     if( electrons.size() < 2 && muons.size() < 2 ) continue;



     std::vector< AnalysisLepton > leptons;

     if( electrons.size() == 2 && muons.size() == 2 ) { //veto H->ZZ->4l

       continue;

     } else if( electrons.size() == 2 ) {

       leptType_ = 1;

       if( electrons[0].Pt() > electrons[1].Pt() ) {

         leptons.push_back( electrons[0] );
         leptons.push_back( electrons[1] );

       } else {

         leptons.push_back( electrons[1] );
         leptons.push_back( electrons[0] );

       }

     } else if( muons.size() == 2 ) {

       leptType_ = 0;

       if( muons[0].Pt() > muons[1].Pt() ) {

         leptons.push_back( muons[0] );
         leptons.push_back( muons[1] );

       } else {

         leptons.push_back( muons[1] );
         leptons.push_back( muons[0] );

       }

     } else {

       std::cout << "There must be an error this is not possible." << std::endl;
       exit(9101);

     }

     eLept1_ = leptons[0].Energy();
     ptLept1_ = leptons[0].Pt();
     etaLept1_ = leptons[0].Eta();
     phiLept1_ = leptons[0].Phi();
     chargeLept1_ = leptons[0].charge;
     
     eLept2_ = leptons[1].Energy();
     ptLept2_ = leptons[1].Pt();
     etaLept2_ = leptons[1].Eta();
     phiLept2_ = leptons[1].Phi();
     chargeLept2_ = leptons[1].charge;


     // --------------------
     // match leptons to MC:
     // --------------------
     int correctIdMc = (leptType_==0 ) ? 13 : 11;

     for( unsigned iLept=0; iLept<leptons.size(); ++iLept ) {

       float deltaRmin = 100.;
       TLorentzVector matchedLeptonMC;

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         if( statusMc[iMc]==1 && fabs(idMc[iMc])==correctIdMc && idMc[mothMc[mothMc[iMc]]]==23 ) {

           TLorentzVector* thisParticle = new TLorentzVector();
           thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );
           float thisDeltaR = leptons[iLept].DeltaR( *thisParticle );
           if( thisDeltaR < deltaRmin ) {
             deltaRmin = thisDeltaR;
             matchedLeptonMC = *thisParticle;
           }

           delete thisParticle;
           thisParticle = 0;

         } //if correct id mc

       } // for i mc

       if( !noLeptons ) {
         if( leptType_==0 ) {
           h1_deltaRmatching_muons->Fill( deltaRmin );
           if( deltaRmin<0.1 ) {
             h1_passed_vs_ptMuon->Fill( matchedLeptonMC.Pt() );
           }
         } else if( leptType_==1 ) { 
           h1_deltaRmatching_electrons->Fill( deltaRmin );
           if( deltaRmin<0.1 ) {
             h1_passed_vs_ptEle->Fill( matchedLeptonMC.Pt() );
           }
         }  //if lept type
       } //if yes leptons


     } //for i leptons



     // ------------------
     // JETS
     // ------------------

     float jetPt_thresh = 30.;

     // first save leading jets in event:
     std::vector<AnalysisJet> leadJets;
     std::vector<int> leadJetsIndex; //index in the event collection (needed afterwards for PFCandidates)
     int nJets30=0;

     for( unsigned int iJet=0; iJet<nAK5PFPUcorrJet; ++iJet ) {

       AnalysisJet thisJet( pxAK5PFPUcorrJet[iJet], pyAK5PFPUcorrJet[iJet], pzAK5PFPUcorrJet[iJet], energyAK5PFPUcorrJet[iJet] );

       thisJet.eChargedHadrons = chargedHadronEnergyAK5PFPUcorrJet[iJet];
       thisJet.ePhotons        = photonEnergyAK5PFPUcorrJet[iJet];
       thisJet.eNeutralEm      = neutralEmEnergyAK5PFPUcorrJet[iJet];
       thisJet.eNeutralHadrons = neutralHadronEnergyAK5PFPUcorrJet[iJet];
       thisJet.eElectrons      = electronEnergyAK5PFPUcorrJet[iJet];
       thisJet.eMuons          = muonEnergyAK5PFPUcorrJet[iJet];

       thisJet.nChargedHadrons = chargedHadronMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nPhotons        = photonMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nNeutralHadrons = neutralHadronMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nElectrons      = electronMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nMuons          = muonMultiplicityAK5PFPUcorrJet[iJet];

       thisJet.nCharged = chargedHadronMultiplicityAK5PFPUcorrJet[iJet]+electronMultiplicityAK5PFPUcorrJet[iJet]+muonMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nNeutral = neutralHadronMultiplicityAK5PFPUcorrJet[iJet]+photonMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.rmsCand =  rmsCandAK5PFPUcorrJet[iJet];
       thisJet.ptD =  ptDAK5PFPUcorrJet[iJet];

       thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.jetBProbabilityBJetTag = jetBProbabilityBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.jetProbabilityBJetTag = jetProbabilityBJetTagsAK5PFPUcorrJet[iJet];

       //thisJet.QGlikelihood = qglc.ComputeLikelihood( thisJet.Pt(), thisJet.nCharged, thisJet.nNeutral, thisJet.ptD, thisJet.rmsCand );

       if( thisJet.Pt()>jetPt_thresh ) nJets30++;

       // save at least 3 lead jets (if event has them) and all jets with pt>thresh:
       if( leadJets.size()>=3 && thisJet.Pt()<jetPt_thresh ) break;

       // far away from leptons:
       if( thisJet.DeltaR( leptons[0] ) <= 0.5 ) continue;
       if( thisJet.DeltaR( leptons[1] ) <= 0.5 ) continue;

       // jet ID:
       int multiplicity = thisJet.nCharged +  thisJet.nNeutral + HFEMMultiplicityAK5PFPUcorrJet[iJet] + HFHadronMultiplicityAK5PFPUcorrJet[iJet];
       if( multiplicity < 2 ) continue;
       if( fabs(thisJet.Eta())<2.4 && thisJet.nChargedHadrons == 0 ) continue;
       if( thisJet.eNeutralHadrons >= 0.99*thisJet.Energy() ) continue;
       if( thisJet.ePhotons >= 0.99*thisJet.Energy() ) continue;

       // match to genjet:
       float bestDeltaR=999.;
       TLorentzVector matchedGenJet;
       for( unsigned iGenJet=0; iGenJet<nAK5GenJet; ++iGenJet ) {
         TLorentzVector thisGenJet(pxAK5GenJet[iGenJet], pyAK5GenJet[iGenJet], pzAK5GenJet[iGenJet], energyAK5GenJet[iGenJet]);
         if( thisGenJet.DeltaR(thisJet) < bestDeltaR ) {
           bestDeltaR=thisGenJet.DeltaR(thisJet);
           matchedGenJet=thisGenJet;
         }
       }
   
       
       thisJet.ptGen  = (isMC_) ? matchedGenJet.Pt() : 0.;
       thisJet.etaGen = (isMC_) ? matchedGenJet.Eta() : 20.;
       thisJet.phiGen = (isMC_) ? matchedGenJet.Phi() : 0.;
       thisJet.eGen   = (isMC_) ? matchedGenJet.Energy() : 0.;

       leadJets.push_back(thisJet);
       leadJetsIndex.push_back(iJet);

     }


     h1_nJets30->Fill(nJets30);
     if( leadJets.size()<2 ) continue;
     if( leadJets[1].Pt()<jetPt_thresh ) continue; //at least 2 jets over thresh



   //// now look for best invariant mass jet pair 
   //float Zmass = 91.19;
   //float bestMass = 0.;
   //int best_i=-1;
   //int best_j=-1;
   //int best_i_eventIndex=-1;
   //int best_j_eventIndex=-1;

     nPairs_ = 0;
     nPart_ = 0;


     for( unsigned iJet=0; iJet<leadJets.size(); ++iJet ) {
   
       AnalysisJet thisJet = leadJets[iJet];

       // --------------
       // kinematics:
       // --------------
       if( thisJet.Pt() < jetPt_thresh ) continue;
       if( fabs(thisJet.Eta()) > 2.4 ) continue;


       for( unsigned int jJet=iJet+1; jJet<leadJets.size(); ++jJet ) {

         AnalysisJet otherJet = leadJets[jJet];

         // --------------
         // kinematics:
         // --------------
         if( otherJet.Pt() < jetPt_thresh ) continue;
         if( fabs(otherJet.Eta()) > 2.4 ) continue;


         if( nPairs_>=50 ) {
        
           std::cout << "MORE than 50 jet pairs found. SKIPPING!!" << std::endl;

         } else {

           eJet1_[nPairs_] = leadJets[iJet].Energy();
           ptJet1_[nPairs_] = leadJets[iJet].Pt();
           etaJet1_[nPairs_] = leadJets[iJet].Eta();
           phiJet1_[nPairs_] = leadJets[iJet].Phi();
           eChargedHadronsJet1_[nPairs_] = leadJets[iJet].eChargedHadrons;
           ePhotonsJet1_[nPairs_]        = leadJets[iJet].ePhotons;
           eNeutralEmJet1_[nPairs_]      = leadJets[iJet].eNeutralEm;
           eNeutralHadronsJet1_[nPairs_] = leadJets[iJet].eNeutralHadrons;
           eElectronsJet1_[nPairs_]      = leadJets[iJet].eElectrons;
           eMuonsJet1_[nPairs_]          = leadJets[iJet].eMuons;
           nChargedHadronsJet1_[nPairs_] = leadJets[iJet].nChargedHadrons;
           nPhotonsJet1_[nPairs_]        = leadJets[iJet].nPhotons;
           nNeutralHadronsJet1_[nPairs_] = leadJets[iJet].nNeutralHadrons;
           nElectronsJet1_[nPairs_]      = leadJets[iJet].nElectrons;
           nMuonsJet1_[nPairs_]          = leadJets[iJet].nMuons;

           ptDJet1_[nPairs_] = leadJets[iJet].ptD;
           rmsCandJet1_[nPairs_] = leadJets[iJet].rmsCand;
           nChargedJet1_[nPairs_] = leadJets[iJet].nCharged;
           nNeutralJet1_[nPairs_] = leadJets[iJet].nNeutral;
           QGlikelihoodJet1_[nPairs_] = leadJets[iJet].QGlikelihood;

           trackCountingHighEffBJetTagJet1_[nPairs_] = leadJets[iJet].trackCountingHighEffBJetTag;
           trackCountingHighPurBJetTagJet1_[nPairs_] = leadJets[iJet].trackCountingHighPurBJetTag;
           simpleSecondaryVertexHighEffBJetTagJet1_[nPairs_] = leadJets[iJet].simpleSecondaryVertexHighEffBJetTag;
           simpleSecondaryVertexHighPurBJetTagJet1_[nPairs_] = leadJets[iJet].simpleSecondaryVertexHighPurBJetTag;
           jetBProbabilityBJetTagJet1_[nPairs_] = leadJets[iJet].jetBProbabilityBJetTag;
           jetProbabilityBJetTagJet1_[nPairs_] = leadJets[iJet].jetProbabilityBJetTag;


           eJet1Gen_[nPairs_] = leadJets[iJet].eGen;
           ptJet1Gen_[nPairs_] = leadJets[iJet].ptGen;
           etaJet1Gen_[nPairs_] = leadJets[iJet].etaGen;
           phiJet1Gen_[nPairs_] = leadJets[iJet].phiGen;
            
           eJet2_[nPairs_] = leadJets[jJet].Energy();
           ptJet2_[nPairs_] = leadJets[jJet].Pt();
           etaJet2_[nPairs_] = leadJets[jJet].Eta();
           phiJet2_[nPairs_] = leadJets[jJet].Phi();
           eChargedHadronsJet2_[nPairs_] = leadJets[jJet].eChargedHadrons;
           ePhotonsJet2_[nPairs_]        = leadJets[jJet].ePhotons;
           eNeutralEmJet2_[nPairs_]      = leadJets[jJet].eNeutralEm;
           eNeutralHadronsJet2_[nPairs_] = leadJets[jJet].eNeutralHadrons;
           eElectronsJet2_[nPairs_]      = leadJets[jJet].eElectrons;
           eMuonsJet2_[nPairs_]          = leadJets[jJet].eMuons;
           nChargedHadronsJet2_[nPairs_] = leadJets[jJet].nChargedHadrons;
           nPhotonsJet2_[nPairs_]        = leadJets[jJet].nPhotons;
           nNeutralHadronsJet2_[nPairs_] = leadJets[jJet].nNeutralHadrons;
           nElectronsJet2_[nPairs_]      = leadJets[jJet].nElectrons;
           nMuonsJet2_[nPairs_]          = leadJets[jJet].nMuons;

           ptDJet2_[nPairs_] = leadJets[jJet].ptD;
           rmsCandJet2_[nPairs_] = leadJets[jJet].rmsCand;
           nChargedJet2_[nPairs_] = leadJets[jJet].nCharged;
           nNeutralJet2_[nPairs_] = leadJets[jJet].nNeutral;
           QGlikelihoodJet2_[nPairs_] = leadJets[jJet].QGlikelihood;

           trackCountingHighEffBJetTagJet2_[nPairs_] = leadJets[jJet].trackCountingHighEffBJetTag;
           trackCountingHighPurBJetTagJet2_[nPairs_] = leadJets[jJet].trackCountingHighPurBJetTag;
           simpleSecondaryVertexHighEffBJetTagJet2_[nPairs_] = leadJets[jJet].simpleSecondaryVertexHighEffBJetTag;
           simpleSecondaryVertexHighPurBJetTagJet2_[nPairs_] = leadJets[jJet].simpleSecondaryVertexHighPurBJetTag;
           jetBProbabilityBJetTagJet2_[nPairs_] = leadJets[jJet].jetBProbabilityBJetTag;
           jetProbabilityBJetTagJet2_[nPairs_] = leadJets[jJet].jetProbabilityBJetTag;

           eJet2Gen_[nPairs_] = leadJets[jJet].eGen;
           ptJet2Gen_[nPairs_] = leadJets[jJet].ptGen;
           etaJet2Gen_[nPairs_] = leadJets[jJet].etaGen;
           phiJet2Gen_[nPairs_] = leadJets[jJet].phiGen;
            

           nPairs_++;
          
         }


       } //for j
     } //for i

     
     //if( jets.size() < 2 ) continue;
     //if( best_i==-1 || best_j==-1 ) continue; //means that less than 2 jets were found




/*

    // and now full kinematic fit with PFCands:
    TLorentzVector ZZ_kinfit_cands;
    TRegexp cands_tstr("CANDS");
    TString dataset_tstr(dataset_);
    if( dataset_tstr.Contains(cands_tstr) ) {

      nPFCand1_=0;
      nPFCand2_=0;
    
      // get the candidates for 2 candidate jets:
      for( unsigned iCand=0; iCand<nPFCand; ++iCand) {

        TLorentzVector thisCand( pxPFCand[iCand], pyPFCand[iCand], pzPFCand[iCand], energyPFCand[iCand] );

        if( iPFJetPFCand[iCand]==best_i_eventIndex ) {

          if( nPFCand1_>=100 ) {

            std::cout << "More than 100 candidates found. Skipping" << std::endl;

          } else {

            ePFCand1_[nPFCand1_] = thisCand.Energy();
            ptPFCand1_[nPFCand1_] = thisCand.Pt();
            etaPFCand1_[nPFCand1_] = thisCand.Eta();
            phiPFCand1_[nPFCand1_] = thisCand.Phi();
            particleTypePFCand1_[nPFCand1_] = particleTypePFCand[iCand];

            nPFCand1_++;

          }


        } // if best_i

        if( iPFJetPFCand[iCand]==best_j_eventIndex ) {

          if( nPFCand2_>=100 ) {

            std::cout << "More than 100 candidates found. Skipping" << std::endl;

          } else {

            ePFCand2_[nPFCand2_] = thisCand.Energy();
            ptPFCand2_[nPFCand2_] = thisCand.Pt();
            etaPFCand2_[nPFCand2_] = thisCand.Eta();
            phiPFCand2_[nPFCand2_] = thisCand.Phi();
            particleTypePFCand2_[nPFCand2_] = particleTypePFCand[iCand];

            nPFCand2_++;

          }

        } // if best_j

      }  // for candidates
     
     }

*/


     if( isMC_ ) {

       // store event partons in tree:
       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         //if( statusMc[iMc]==3 && (fabs(idMc[iMc])<=6 || idMc[iMc]==21) ) {
         if( statusMc[iMc]==3 && pMc[iMc]*sin(thetaMc[iMc])>0.1 ) {

           TLorentzVector* thisParticle = new TLorentzVector();
           thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

           if( nPart_<20 ) {

             ptPart_[nPart_] = thisParticle->Pt();
             etaPart_[nPart_] = thisParticle->Eta();
             phiPart_[nPart_] = thisParticle->Phi();
             ePart_[nPart_] = thisParticle->Energy();
             pdgIdPart_[nPart_] = idMc[iMc];
             motherPart_[nPart_] = idMc[mothMc[iMc]];

             nPart_++;

           } else {
      
             std::cout << "Found more than 20 partons, skipping." << std::endl;

           }

           delete thisParticle;
           thisParticle = 0;

         } //if correct id mc

       } // for i mc

     } // if is mc


     reducedTree_->Fill(); 


   } //for entries


} //loop



double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}



/*
Double_t Ntp1Analyzer_HZZlljj::computeQGLikelihood(const Double_t jtpt, Int_t ncharged, Int_t nneutral, Double_t PtD, Double_t r) {

//std::vector<variable*> vars;
//  std::vector<variable*> spects;
  TMVA::Reader *reader=new TMVA::Reader("Reader");
    Double_t result;
//AddVariables(vars);
//  AddSpectators(spects);
//  AddVariableAndSpectators(reader,vars,spects);
  
std::cout << "jtpt: " << jtpt << std::endl;
  if(15.<jtpt&&jtpt<30.)reader->BookMVA("PDEFoam","Bins/15-30/weights/TMVA Analysis_PDEFoam.weights.xml"); else
  if(30.<jtpt&&jtpt<50.)reader->BookMVA("PDEFoam","Bins/30-50/weights/TMVA Analysis_PDEFoam.weights.xml"); else
  if(50.<jtpt&&jtpt<80.)reader->BookMVA("PDEFoam","Bins/50-80/weights/TMVA Analysis_PDEFoam.weights.xml"); else
  if(80.<jtpt&&jtpt<120.){
    std::cout << "in here" << std::endl;
    reader->BookMVA("PDEFoam","Bins/80-120/weights/TMVA Analysis_PDEFoam.weights.xml"); 
    std::cout << "done" << std::endl;
    } else
  if(120.<jtpt&&jtpt<170.)reader->BookMVA("PDEFoam","Bins/120-170/weights/TMVA Analysis_PDEFoam.weights.xml"); else
  if(170.<jtpt&&jtpt<300.)reader->BookMVA("PDEFoam","Bins/170-300/weights/TMVA Analysis_PDEFoam.weights.xml"); else
  if(300.<jtpt&&jtpt<470.)reader->BookMVA("PDEFoam","Bins/300-470/weights/TMVA Analysis_PDEFoam.weights.xml"); else
  if(470.<jtpt&&jtpt<600.)reader->BookMVA("PDEFoam","Bins/470-6000/weights/TMVA Analysis_PDEFoam.weights.xml");else
  return -1.0; //dentro l'else
  
//std::vector<variable*>::iterator it;
//  for(it=vars.begin();it!=vars.end();it++)
//          {
//          if((*it)->name=="ncharged") (*it)->SetVar(&ncharged); else
//          if((*it)->name=="nneutral") (*it)->SetVar(&nneutral); else
//          if((*it)->name=="PtD") (*it)->SetVar(&PtD);else
//          if((*it)->name== "r") (*it)->SetVar(&r);else
//          {return -2.0;}
//          
//          }
std::cout << "jaja" << std::endl;
  result=reader->EvaluateMVA("PDEFoam");
  return result;
  
}
*/
