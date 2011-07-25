#include "Ntp1Analyzer_HWWlvjj.h"

#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"
#include "TMVA/Reader.h"

#include "AnalysisElectron.h"
#include "AnalysisMuon.h"
#include "AnalysisNeutrino.h"

#include "Utils.hh"

//#include "fitTools.h"

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
};

Ntp1Analyzer_HWWlvjj::Ntp1Analyzer_HWWlvjj( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "HWWlvjj", dataset, flags, tree ) {

     std::string dataset_=dataset;
  //nothing to do here

} //constructor

/// specific for HWW that has multiple channels with different HLT requirements
void Ntp1Analyzer_HWWlvjj::reloadTriggerMask(int runN){
  std::vector<int> triggerMask;

  // load the triggers required for E
  for (std::vector< std::string >::const_iterator fIter=requiredTriggerElectron.begin();fIter!=requiredTriggerElectron.end();++fIter)
    {   
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          //if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
          // nameHLT[i] has ..._vXXX
          if(nameHLT->at(i).find(pathName) != std::string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggersElectron = triggerMask;

  // load the triggers required for M
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=requiredTriggerMuon.begin();fIter!=requiredTriggerMuon.end();++fIter)
    {   
      //      std::cout << "For MM required: " << *fIter << std::endl;
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(pathName) != std::string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggersMuon = triggerMask;
}

bool Ntp1Analyzer_HWWlvjj::hasPassedHLT(int channel) {
  Utils anaUtils;
  if(channel==0) return anaUtils.getTriggersOR(m_requiredTriggersElectron, firedTrg);
  else if(channel==1) {
    bool required = anaUtils.getTriggersOR(m_requiredTriggersMuon, firedTrg);
    //bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersMuon, firedTrg);
    return (required );//&& !notRequired);
  //} else if(channel==em) {
  //  bool required = anaUtils.getTriggersOR(m_requiredTriggersEM, firedTrg);
   // bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersEM, firedTrg);
   // return (required && !notRequired);
  }
  return true;
}

std::string Ntp1Analyzer_HWWlvjj::getHLTPathForRun(int runN, std::string fullname) {
  TString fullName = TString(fullname.c_str());
  TObjArray* selectionTokens = fullName.Tokenize(":");
  if (selectionTokens->GetEntries()!=2) {
    std::cout << "Wrong trigger strings " << selectionTokens->GetEntries() << std::endl;
    return std::string("NOPATH");
  }
  TString RunRange =((TObjString*)(*selectionTokens)[0])->GetString();
  TString HLTPathName =((TObjString*)(*selectionTokens)[1])->GetString();
  
  TObjArray* runs = RunRange.Tokenize("-");
  if (runs->GetEntries()!=2) {
    std::cout << "Wrong trigger run range strings " << runs->GetEntries() << std::endl;
    return std::string("NOPATH");    
  }
  
  const char *minStr = (((TObjString*)(*runs)[0])->GetString()).Data();
  const char *maxStr = (((TObjString*)(*runs)[1])->GetString()).Data();

  int min = atoi(minStr);
  int max = atoi(maxStr);

  if(runN>=min && runN<=max) return std::string(HLTPathName.Data());
  else return std::string("NOPATH");
}


void Ntp1Analyzer_HWWlvjj::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();
  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");
  reducedTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");

//// triggers:
//reducedTree_->Branch("HLT_Mu11", &HLT_Mu11_, "HLT_Mu11_/O");
//reducedTree_->Branch("HLT_Ele17_SW_EleId_L1R", &HLT_Ele17_SW_EleId_L1R_, "HLT_Ele17_SW_EleId_L1R_/O");
//reducedTree_->Branch("HLT_DoubleMu3", &HLT_DoubleMu3_, "HLT_DoubleMu3_/O");


  reducedTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  reducedTree_->Branch("leptType",  &leptType_,  "leptType_/I");
  
  reducedTree_->Branch("eWqqMC",  &eWqqMC_,  "eWqqMC_/F");
  reducedTree_->Branch("ptWqqMC",  &ptWqqMC_,  "ptWqqMC_/F");
  reducedTree_->Branch("etaWqqMC",  &etaWqqMC_,  "etaWqqMC_/F");
  reducedTree_->Branch("phiWqqMC",  &phiWqqMC_,  "phiWqqMC_/F");

  reducedTree_->Branch("eWllMC",  &eWllMC_,  "eWllMC_/F");
  reducedTree_->Branch("ptWllMC",  &ptWllMC_,  "ptWllMC_/F");
  reducedTree_->Branch("etaWllMC",  &etaWllMC_,  "etaWllMC_/F");
  reducedTree_->Branch("phiWllMC",  &phiWllMC_,  "phiWllMC_/F");

  reducedTree_->Branch("eHiggsMC",  &eHiggsMC_,  "eHiggsMC_/F");
  reducedTree_->Branch("ptHiggsMC",  &ptHiggsMC_,  "ptHiggsMC_/F");
  reducedTree_->Branch("etaHiggsMC",  &etaHiggsMC_,  "etaHiggsMC_/F");
  reducedTree_->Branch("phiHiggsMC",  &phiHiggsMC_,  "phiHiggsMC_/F");

  reducedTree_->Branch("eQuark1",  &eQuark1_,  "eQuark1_/F");
  reducedTree_->Branch("ptQuark1",  &ptQuark1_,  "ptQuark1_/F");
  reducedTree_->Branch("etaQuark1",  &etaQuark1_,  "etaQuark1_/F");
  reducedTree_->Branch("phiQuark1",  &phiQuark1_,  "phiQuark1_/F");

  reducedTree_->Branch("eQuark2",  &eQuark2_,  "eQuark2_/F");
  reducedTree_->Branch("ptQuark2",  &ptQuark2_,  "ptQuark2_/F");
  reducedTree_->Branch("etaQuark2",  &etaQuark2_,  "etaQuark2_/F");
  reducedTree_->Branch("phiQuark2",  &phiQuark2_,  "phiQuark2_/F");

  reducedTree_->Branch("eLept",  &eLept_,  "eLept_/F");
  reducedTree_->Branch("ptLept",  &ptLept_,  "ptLept_/F");
  reducedTree_->Branch("etaLept",  &etaLept_,  "etaLept_/F");
  reducedTree_->Branch("phiLept",  &phiLept_,  "phiLept_/F");
  reducedTree_->Branch("chargeLept",  &chargeLept_,  "chargeLept_/I");

  reducedTree_->Branch("eLeptMC",  &eLeptMC_,  "eLeptMC_/F");
  reducedTree_->Branch("ptLeptMC",  &ptLeptMC_,  "ptLeptMC_/F");
  reducedTree_->Branch("etaLeptMC",  &etaLeptMC_,  "etaLeptMC_/F");
  reducedTree_->Branch("phiLeptMC",  &phiLeptMC_,  "phiLeptMC_/F");
  reducedTree_->Branch("eNeuMC",  &eNeuMC_,  "eNeuMC_/F");
  reducedTree_->Branch("ptNeuMC",  &ptNeuMC_,  "ptNeuMC_/F");
  reducedTree_->Branch("etaNeuMC",  &etaNeuMC_,  "etaNeuMC_/F");
  reducedTree_->Branch("phiNeuMC",  &phiNeuMC_,  "phiNeuMC_/F");

  reducedTree_->Branch("nPairs", &nPairs_, "nPairs_/I");

  reducedTree_->Branch("trackCountingHighEffBJetTagJet1", trackCountingHighEffBJetTagJet1_, "trackCountingHighEffBJetTagJet1_[nPairs_]/F");
  reducedTree_->Branch("trackCountingHighEffBJetTagJet2", trackCountingHighEffBJetTagJet2_, "trackCountingHighEffBJetTagJet2[nPairs_]/F");

  reducedTree_->Branch("iJet1",  iJet1_,  "iJet1_[nPairs_]/I");
  reducedTree_->Branch("eJet1",  eJet1_,  "eJet1_[nPairs_]/F");
  reducedTree_->Branch("ptJet1",  ptJet1_,  "ptJet1_[nPairs_]/F");
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

  reducedTree_->Branch("eJet1Gen",  eJet1Gen_,  "eJet1Gen_[nPairs_]/F");
  reducedTree_->Branch("ptJet1Gen",  ptJet1Gen_,  "ptJet1Gen_[nPairs_]/F");
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
  reducedTree_->Branch("ptPart",  ptPart_,  "ptPart_[nPart_]/F");
  reducedTree_->Branch("etaPart", etaPart_, "etaPart_[nPart_]/F");
  reducedTree_->Branch("phiPart", phiPart_, "phiPart_[nPart_]/F");
  reducedTree_->Branch("pdgIdPart", pdgIdPart_, "pdgIdPart_[nPart_]/I");
  reducedTree_->Branch("motherPart", motherPart_, "motherPart_[nPart_]/I");

  reducedTree_->Branch("energyPFMet",&energyPFMet_,"energyPFMet_/F");
  reducedTree_->Branch("phiPFMet",&phiPFMet_,"phiPFMet_/F");
  reducedTree_->Branch("pxPFMet",&pxPFMet_,"pxPFMet_/F");
  reducedTree_->Branch("pyPFMet",&pyPFMet_,"pyPFMet_/F");

  reducedTree_->Branch("uncorrEnergyAK5Jet",&uncorrEnergyAK5Jet_,"uncorrEnergyAK5Jet_/F");
  reducedTree_->Branch("SumEt",&SumEt_,"SumEt_/F");
 
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
  h1_deltaRmatching_jet_leptonParton = new TH1F("deltaRmatching_jet_leptonParton", "", 100, 0., 4.);
  h1_nJets30 = new TH1F("nJets30", "", 31, -0.5, 30.5);
//h1_indexMatchedJet = new TH1F("indexMatchedJet", "", 6, -0.5, 5.5);
//h1_indexMatched05Jet = new TH1F("indexMatched05Jet", "", 6, -0.5, 5.5);
//h1_nMatched_per_event = new TH1F("nMatched_per_event", "", 6, -0.5, 5.5);
//h1_nMatched05_per_event = new TH1F("nMatched05_per_event", "", 6, -0.5, 5.5);
//h1_pdgIdParton1 = new TH1F("pdgIdParton1", "", 36, -10.5, 25.5);
//h1_pdgIdParton2 = new TH1F("pdgIdParton2", "", 36, -10.5, 25.5);
//h1_ptHadronicW = new TH1F("ptHadronicW", "", 50, 0., 400.);
//h1_deltaRqq = new TH1F("deltaRqq", "", 50, 0., 3.);
//h1_deltaRqq = new TH1F("deltaRqq", "", 50, 0., 3.);

  h1_Cont_inclusive = new TH1F("Cont_inclusive", "", 1, 0., 1.);
  h1_Cont_PV = new TH1F("Cont_PV", "", 1, 0., 1.);
  h1_Cont_TightMu = new TH1F("Cont_TightMu", "", 1, 0., 1.);
  h1_Cont_TightEle = new TH1F("Cont_TightEle", "", 1, 0., 1.);
  h1_Cont_VetoMU = new TH1F("Cont_VetoMU", "", 1, 0., 1.);
  h1_Cont_VetoELE = new TH1F("Cont_VetoELE", "", 1, 0., 1.);
  h1_Cont_JetsELE = new TH1F("Cont_JetsELE", "", 1, 0., 1.);
  h1_Cont_JetsMU = new TH1F("Cont_JetsMU", "", 1, 0., 1.);

} 
 
Ntp1Analyzer_HWWlvjj::~Ntp1Analyzer_HWWlvjj() {

  outfile_->cd();

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
//h1_ptHadronicW->Write();
//h1_deltaRqq->Write(); 

  h1_Cont_inclusive->Write();
  h1_Cont_PV->Write();
  h1_Cont_TightMu->Write();
  h1_Cont_TightEle->Write();
  h1_Cont_VetoMU->Write();
  h1_Cont_VetoELE->Write();
  h1_Cont_JetsELE->Write();
  h1_Cont_JetsMU->Write();

}

void Ntp1Analyzer_HWWlvjj::Loop(){

   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;

   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   TRandom3 rand;
   float Cont_inclusive=0., Cont_PV=0., Cont_MU=0., Cont_ELE=0., Cont_VetoMU=0., Cont_VetoELE=0., Cont_JetsELE=0., Cont_JetsMU=0.;
   Long64_t Jentry;

// QUI ci metti il 18 di HiggsApp
std::vector< std::string > requiredTriggerElectron;//sarebbe maskEE
std::vector< std::string > requiredTriggerMuon;
requiredTriggerElectron.push_back("1-164237:HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v");
requiredTriggerElectron.push_back("165085-166967:HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v");
requiredTriggerElectron.push_back("166968-999999:HLT_Ele52_CaloIdVT_TrkIdT_v");
requiredTriggerMuon.push_back("1-163261:HLT_Mu15_v");
requiredTriggerMuon.push_back("163262-164237:HLT_Mu24_v");
requiredTriggerMuon.push_back("165085-166967:HLT_Mu30_v");
requiredTriggerMuon.push_back("163262-166967:HLT_IsoMu17_v");
requiredTriggerMuon.push_back("167039-999999:HLT_IsoMu20_eta2p1_v");

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

     Jentry=jentry;

     if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;
     

     if( (jentry%100000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;
  
   //HLT_Mu11_ = this->PassedHLT("HLT_Mu11");
   //HLT_Ele17_SW_EleId_L1R_ = this->PassedHLT("HLT_Ele17_SW_EleId_L1R");
   //HLT_DoubleMu3_ = this->PassedHLT("HLT_DoubleMu3");

     run_ = runNumber;
     LS_ = lumiBlock;
     event_ = eventNumber;
     eventWeight_ = -1.; //default

     Cont_inclusive++;//Other;

     if( !isGoodEvent( jentry ) ) continue; //this takes care also of integrated luminosity and trigger

       if( nPV==0 ) continue;
     bool goodVertex = (ndofPV[0] >= 4.0 && sqrt(PVxPV[0]*PVxPV[0]+PVyPV[0]*PVyPV[0]) < 2. && fabs(PVzPV[0]) < 24. );
     if( !goodVertex ) continue;
     nvertex_ = nPV;
     rhoPF_ = rhoFastjet;

     Cont_PV++; //Other

     //trigger:
     // not yet
     
     if( !isMC_ ){
     reloadTriggerMask(runNumber);
     bool passedElectronTrigger=true, passedMuonTrigger=true;
     std::string SingEle("SingleElectron");
     std::string SingMuo("SingleMu");
     int elect=dataset_.compare(SingEle), muo=dataset_.compare(SingMuo);
     if( elect==0 ){ passedElectronTrigger = hasPassedHLT(0);}
     if( muo==0 ){ passedMuonTrigger = hasPassedHLT(1);}
     if( elect!=0 && muo!=0  ){ std::cout<<"Trigger doesn't works; dataset wrong"<<std::endl; }
     if( !passedElectronTrigger || !passedMuonTrigger ) continue;
    }

     ptHat_ = (isMC_) ? genPtHat : ptHat_;

     //if( isMC_ ) 
     //  if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;

     bool noLeptons = false;
     TLorentzVector lept1MC, lept2MC;
     int wIndexqq=-1;
     int wIndexll=-1;
    
     if( isMC_ ) {
       
       // first look for W->qq'
       std::vector<TLorentzVector> quarksMC;

       for( unsigned iMc=0; iMc<nMc && quarksMC.size()<2; ++iMc ) {

         // quarks have status 3
         if( statusMc[iMc] != 3 ) continue;

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

         if( fabs(idMc[iMc])<7 && fabs(idMc[mothMc[iMc]])==24 ) {
           wIndexqq = mothMc[iMc];
           quarksMC.push_back( *thisParticle );
         }
       }

       // (checked that always 2 quarks are found)
       if( quarksMC.size()==2 && wIndexqq!=-1 ) {
	   eQuark1_=quarksMC[0].E();
	   ptQuark1_=quarksMC[0].Pt();
	   etaQuark1_=quarksMC[0].Eta();
	   phiQuark1_=quarksMC[0].Phi();
	   eQuark2_=quarksMC[1].E();
	   ptQuark2_=quarksMC[1].Pt();
	   etaQuark2_=quarksMC[1].Eta();
	   phiQuark2_=quarksMC[1].Phi();
 
         TLorentzVector WqqMC;
         WqqMC.SetPtEtaPhiE( pMc[wIndexqq]*sin(thetaMc[wIndexqq]), etaMc[wIndexqq], phiMc[wIndexqq], energyMc[wIndexqq] );
         ptWqqMC_  = WqqMC.Pt();
         eWqqMC_   = WqqMC.Energy();
         etaWqqMC_ = WqqMC.Eta();
         phiWqqMC_ = WqqMC.Phi();
	 //float deltaRqq = quarksMC[0].DeltaR(quarksMC[1]);
	 // h1_deltaRqq->Fill(deltaRqq);
       }
       
       // now look for W->lv

       std::vector<TLorentzVector> electronMC;
       std::vector<TLorentzVector> muonMC;
       std::vector<TLorentzVector> neutrinoMC;

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {
	 
	 if( statusMc[iMc] != 3 ) continue;
	 
         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

           if( fabs(idMc[iMc])==11 && fabs(idMc[mothMc[iMc]])==24 ) {electronMC.push_back( *thisParticle ); wIndexll = mothMc[iMc]; }
	   if( fabs(idMc[iMc])==13 && fabs(idMc[mothMc[iMc]])==24 ) {muonMC.push_back( *thisParticle ); wIndexll = mothMc[iMc]; }
	   if( (fabs(idMc[iMc])==12 || fabs(idMc[iMc])==14 ) && fabs(idMc[mothMc[iMc]])==24 ) {
	     neutrinoMC.push_back( *thisParticle );  wIndexll = mothMc[iMc];
	       }//for comparison with fit
 
         delete thisParticle;
         thisParticle = 0;
       }
       
       if( electronMC.size()==1 && neutrinoMC.size()==1 ) {
	 lept1MC = electronMC[0];
	 lept2MC = neutrinoMC[0];
	
       	 if( (fabs(lept1MC.Eta()) < 2.5) && ( fabs(lept1MC.Eta())<1.4442 || fabs(lept1MC.Eta())>1.566) ){
       	 h1_nEvents_vs_ptEle->Fill( lept1MC.Pt() );}
       	 h1_nEvents_vs_ptEle->Fill( lept2MC.Pt() );
       	 }
       else if( muonMC.size()==1 && neutrinoMC.size()==1 ) {
	 lept1MC = muonMC[0];
	 lept2MC = neutrinoMC[0];

	 if( fabs(lept1MC.Eta()) < 2.4 ){ 
         h1_nEvents_vs_ptMuon->Fill( lept1MC.Pt() );} h1_nEvents_vs_ptMuon->Fill( lept2MC.Pt() );
       }
       else {
	 //taus
	 noLeptons=true;
       }

       if( !noLeptons ) {
       eLeptMC_=lept1MC.E();
       ptLeptMC_=lept1MC.Pt();
       etaLeptMC_=lept1MC.Eta();
       phiLeptMC_=lept1MC.Phi();       
       eNeuMC_=lept2MC.E();
       ptNeuMC_=lept2MC.Pt();
       etaNeuMC_=lept2MC.Eta();
       phiNeuMC_=lept2MC.Phi();
       }
       //if is a tau or a not recognized e/mu default values are setted
       if( noLeptons ) {
       eLeptMC_=1.;
       ptLeptMC_=1.;
       etaLeptMC_=10.;
       phiLeptMC_=1.;       
       eNeuMC_=1.;
       ptNeuMC_=1.;
       etaNeuMC_=10.;
       phiNeuMC_=1.;    
       }

       if( !noLeptons ) {
       TLorentzVector WllMC;
       WllMC.SetPtEtaPhiE( pMc[wIndexll]*sin(thetaMc[wIndexll]), etaMc[wIndexll], phiMc[wIndexll], energyMc[wIndexll] );
       
       ptWllMC_  = WllMC.Pt();
       eWllMC_   = WllMC.Energy();
       etaWllMC_ = WllMC.Eta();
       phiWllMC_ = WllMC.Phi();
       }
       
       // now look for the higgs:
       if( wIndexll!=-1 && wIndexqq!=-1 ) {
	 
         int higgsIndex = mothMc[wIndexll];
         if( idMc[higgsIndex] == 25 ) {
           TLorentzVector HiggsMC;
           HiggsMC.SetPtEtaPhiE( pMc[higgsIndex]*sin(thetaMc[higgsIndex]), etaMc[higgsIndex], phiMc[higgsIndex], energyMc[higgsIndex] );

           eHiggsMC_   = HiggsMC.Energy(); 
           ptHiggsMC_  = HiggsMC.Pt(); 
           etaHiggsMC_ = HiggsMC.Eta(); 
           phiHiggsMC_ = HiggsMC.Phi(); 
         } // if higgs

       } //if found two W's
     
     } //if isMC
     
     // -----------------------------
     //      FROM NOW ON RECO
     // -----------------------------

     energyPFMet_ = energyPFMet[0];
     phiPFMet_ = phiPFMet[0];
     pxPFMet_ = pxPFMet[0];
     pyPFMet_ = pyPFMet[0];     
     
     // -----------------------------
     //###      uncorrEnergyJet (to use in Neutrino's fit if you don't have SumEt)
     // -----------------------------
     
       for( int i=0;i<nAK5Jet;i++ ){
       uncorrEnergyAK5Jet_ += uncorrEnergyAK5Jet[i];
       }

     // -----------------------------
     //###      SumEt
     // -----------------------------
       
     SumEt_ = sumEtPFMet[0]; 
 
     // ------------------
     // MUON
     // ------------------

     std::vector<AnalysisMuon> muon;
     bool looseMuon=false;

     for( unsigned int iMuon=0; iMuon<nMuon; ++iMuon ) {

         AnalysisMuon thisMuon( pxMuon[iMuon], pyMuon[iMuon], pzMuon[iMuon], energyMuon[iMuon] );

       // --------------
       // kinematics:
       // --------------

       if( thisMuon.Pt() < 10. ) continue;
       if( fabs(thisMuon.Eta()) > 2.4 /*2.1 Other*/ ) continue;

	 thisMuon.isGlobalMuonPromptTight = (muonIdMuon[iMuon]>>8)&1;
	 thisMuon.isAllTrackerMuon = (muonIdMuon[iMuon]>>11)&1;

       //     // --------------
       //     // ID:
       //     // --------------
       //     if( !( (muonIdMuon[iMuon]>>8)&1 ) ) continue; //GlobalMuonPromptTight
       //     if( !( (muonIdMuon[iMuon]>>11)&1 ) ) continue; //AllTrackerMuon
       //if( numberOfValidPixelBarrelHitsTrack[trackIndexMuon[iMuon]]==0 && numberOfValidPixelEndcapHitsTrack[trackIndexMuon[iMuon]]==0 ) continue;      
         
	thisMuon.pixelHits = numberOfValidPixelBarrelHitsTrack[trackIndexMuon[iMuon]]+numberOfValidPixelEndcapHitsTrack[trackIndexMuon[iMuon]];
        thisMuon.trackerHits = trackValidHitsTrack[trackIndexMuon[iMuon]];
        thisMuon.nMatchedStations = numberOfMatchesMuon[iMuon]; // This branch not exists yet in WW500

       // to compute dxy, look for leading primary vertex:
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
         dxy = fabs(this->trackDxyPV(PVxPV[hardestPV], PVyPV[hardestPV], PVzPV[hardestPV],
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
     
       // --------------
       // isolation:
       // --------------
       // (this is sum pt tracks)
       //if( sumPt03Muon[iMuon] >= 3. ) continue;
       // combined isolation < 15%:
       //if( (sumPt03Muon[iMuon] + emEt03Muon[iMuon] + hadEt03Muon[iMuon]) >= 0.15*thisMuon.Pt() ) continue;

       // looking at loose muon
       looseMuon=true;
       if( thisMuon.Pt()<20. ) continue;

       // for now simple selection, will have to optimize this (T&P?)
       chargeLept_=chargeMuon[iMuon];
     muon.push_back( thisMuon );
       looseMuon=false;

   } //for muon

     // ------------------
     // ELECTRON
     // ------------------

     std::vector<AnalysisElectron> electron;    
     bool looseEle = false; //Other
     
     for( unsigned int iEle=0; iEle<nEle ; ++iEle ) {

       AnalysisElectron thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );

       // --------------
       // kinematics:
       // --------------
       if( thisEle.Pt() < 15. ) continue;
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

       if( !passed_VBTF80 && passed_VBTF95 ) looseEle=true; //Other
       if( !passed_VBTF80 ) continue; 
       if( thisEle.Pt() < 20. ) continue; // if is not a loose ele I require a higher Pt     

       // check that not matched to muon (clean electron faked by muon MIP):
       bool matchedtomuon=false;
       for( std::vector<AnalysisMuon>::iterator iMu=muon.begin(); iMu!=muon.end(); ++iMu )
         if( iMu->DeltaR(thisEle)<0.1 ) matchedtomuon=true;

       if( matchedtomuon ) continue;

       // for now simple selection, will have to optimize this (T&P?)
       // one electron required to pass VBTF80, the other VBTF95
       chargeLept_=chargeEle[iEle];
       electron.push_back( thisEle );
     } //for electron
                                         //  FILL LEPTONS
// Fill Efficiency for Others
bool findJetELE=false;
bool findJetMU=false;

if( electron.size() >= 1 ) { Cont_ELE++;
  if( electron.size() == 1 && muon.size()==0 && !looseEle && !looseMuon ) { Cont_VetoELE++; findJetELE=true; }
  }

if( muon.size() >= 1 ) { Cont_MU++;
  if( muon.size() == 1 && electron.size()==0 && !looseEle && !looseMuon ) { Cont_VetoMU++; findJetMU=true; }
  }
// Fill lepton
    if( electron.size() > 1 || muon.size() > 1 ) continue;

     if( electron.size() < 1 &&  muon.size() < 1 ) continue;
     
     //   // clean electron faked by muon MIP in ECAL
     //   for( std::vector<TLorentzVector>::iterator iEle=electron.begin(); iEle!=electron.end(); ++iEle ) {
     //     for( std::vector<TLorentzVector>::iterator iMu=muon.begin(); iMu!=muon.end(); ++iMu ) {
     //       if( iMu->DeltaR(*iEle)<0.1 ) {
     //         std::cout << "lll" << std::endl;
     //         electron.erase(iEle);
     //       }
     //     } //for ele
     //   } //for mu
     
     //   if( electron.size() < 2 && muon.size() < 1 ) continue;
    
     std::vector< AnalysisLepton > leptons;

     if( electron.size() == 1 && muon.size() == 1 ) continue;
     if( electron.size() == 1 && (looseEle || looseMuon) ) continue;
     if( muon.size() == 1 && (looseEle || looseMuon) ) continue;

     if( electron.size() == 1 && (!looseEle || !looseMuon) ) { //Other 
       leptType_ = 1;
       leptons.push_back( electron[0] );

     } else if( muon.size() == 1 && (!looseEle || !looseMuon) ) {
       leptType_ = 0;
       leptons.push_back( muon[0] );

     } else {
       std::cout << "There must be an error this is not possible." << std::endl;
       exit(9101);
     }

     eLept_ = leptons[0].Energy();
     ptLept_ = leptons[0].Pt();
     etaLept_ = leptons[0].Eta();
     phiLept_ = leptons[0].Phi(); //chargeLept_ filled before
    
     // --------------------
     // match leptons to MC:
     // --------------------
     int correctIdMc = (leptType_==0 ) ? 13 : 11;

	float deltaRmin = 100.;
	TLorentzVector matchedLeptonMC;
	for( unsigned iMc=0; iMc<nMc; ++iMc ){
	   if( statusMc[iMc]==1 && fabs(idMc[iMc])==correctIdMc && idMc[mothMc[mothMc[iMc]]]==24 ) {
	     
	     TLorentzVector* thisParticle = new TLorentzVector();
	     thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );
	     float thisDeltaR = leptons[0].DeltaR( *thisParticle );
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
   
      // ------------------
      // JETS
      // ------------------
	     
	float jetPt_thresh = 20.;
	bool matched=false;      
	
	// first save leading jets in event:
	std::vector<AnalysisJet> leadJets;
	std::vector<int> leadJetsIndex; //index in the event collection (needed afterwards for PFCandidates)
	int nJets30=0;
	
	//QGLikelihoodCalculator qglc;
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

       if( thisJet.Pt()>jetPt_thresh ) nJets30++;

       // save at least 3 lead jets (if event has them) and all jets with pt>thresh:
       if( leadJets.size()>=3 && thisJet.Pt()<jetPt_thresh ) break;

       // far away from leptons:
      if( leptType_==0){ if(thisJet.DeltaR( leptons[0] ) < 0.5 ) continue;}    
      if( leptType_==1){ if(thisJet.DeltaR( leptons[0] ) < 0.5 ) continue;} 

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

       // match to parton:
       float bestDeltaR_part=999.;
       TLorentzVector matchedPart;
       int pdgIdPart=0;
       for( unsigned iPart=0; iPart<nMc; ++iPart ) {
         if( statusMc[iPart]!=3 ) continue; //partons
         if( idMc[iPart]!=21 && abs(idMc[iPart])>6 ) continue; //quarks or gluons
         TLorentzVector thisPart;
         thisPart.SetPtEtaPhiE(pMc[iPart]*sin(thetaMc[iPart]), etaMc[iPart], phiMc[iPart], energyMc[iPart]);
         if( thisPart.DeltaR(thisJet) < bestDeltaR_part ) {
           bestDeltaR_part=thisPart.DeltaR(thisJet);
           matchedPart=thisPart;
           pdgIdPart=idMc[iPart];
         }
       }

       leadJets.push_back(thisJet);
       leadJetsIndex.push_back(iJet);
     }//iJet

     h1_nJets30->Fill(nJets30);

     if( leadJets.size()<2 ) continue;
     if( leadJets[1].Pt()<jetPt_thresh ) continue; //at least 2 jets over thresh
   
     // now look for best invariant mass jet pair 
     float Wmass = 80.399;
     float bestMass = 0.;
     int best_i=-1;
     int best_j=-1;
     int best_i_eventIndex=-1;
     int best_j_eventIndex=-1;

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

           trackCountingHighEffBJetTagJet1_[nPairs_] = leadJets[iJet].trackCountingHighEffBJetTag;
           trackCountingHighEffBJetTagJet2_[nPairs_] = leadJets[jJet].trackCountingHighEffBJetTag;

           ptDJet1_[nPairs_] = leadJets[iJet].ptD;
           rmsCandJet1_[nPairs_] = leadJets[iJet].rmsCand;
           nChargedJet1_[nPairs_] = leadJets[iJet].nCharged;
           nNeutralJet1_[nPairs_] = leadJets[iJet].nNeutral;
           QGlikelihoodJet1_[nPairs_] = leadJets[iJet].QGlikelihood;

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

           eJet2Gen_[nPairs_] = leadJets[jJet].eGen;
           ptJet2Gen_[nPairs_] = leadJets[jJet].ptGen;
           etaJet2Gen_[nPairs_] = leadJets[jJet].etaGen;
           phiJet2Gen_[nPairs_] = leadJets[jJet].phiGen;
            

           nPairs_++;
          
         }
	
       } //for j
     } //for i
   if( findJetELE && nPairs_>=1 ) {Cont_JetsELE++;} //Other 
   if( findJetMU && nPairs_>=1 ) {Cont_JetsMU++;} //Other 

 
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
                       // Other
float proportion = 2061760. /*109989.*/ /Cont_inclusive;
	if( Jentry == (nentries-1) ){
  std::cout << std::fixed << std::setprecision(6) << "Inclusive: " << Cont_inclusive*proportion << " PV: " << Cont_PV*proportion << " MU: " << Cont_MU*proportion << " ELE: " << Cont_ELE*proportion << " VetoMU: " << Cont_VetoMU*proportion 
  << " VetoELE: " << Cont_VetoELE*proportion << " JetsELE: " << Cont_JetsELE*proportion << " JetsMU: " << Cont_JetsMU*proportion << std::endl;

h1_Cont_inclusive->SetBinContent(1,Cont_inclusive*proportion);
h1_Cont_PV->SetBinContent(1,Cont_PV*proportion);
h1_Cont_TightMu->SetBinContent(1,Cont_MU*proportion);
h1_Cont_TightEle->SetBinContent(1,Cont_ELE*proportion);
h1_Cont_VetoMU->SetBinContent(1,Cont_VetoMU*proportion);
h1_Cont_VetoELE->SetBinContent(1,Cont_VetoELE*proportion);
h1_Cont_JetsELE->SetBinContent(1,Cont_JetsELE*proportion);
h1_Cont_JetsMU->SetBinContent(1,Cont_JetsMU*proportion);
}

} //loop



