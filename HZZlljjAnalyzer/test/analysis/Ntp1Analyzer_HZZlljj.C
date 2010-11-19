#include "Ntp1Analyzer_HZZlljj.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include "fitTools.h"


double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);


class AnalysisJet : public TLorentzVector {

 public:

  AnalysisJet( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {
    eChargedHadrons=0.;
    ePhotons=0.;
    eNeutralHadrons=0.;
    eElectrons=0.;
  }

  float eChargedHadrons;
  float ePhotons;
  float eNeutralHadrons;
//float eMuons;
  float eElectrons;
//float eHFHadrons;
//float eHFEM;

//int nChargedHadrons;
//int nPhotons;
//int nNeutralHadrons;
//int nMuons;
//int nElectrons;
//int nHFHadrons;
//int nHFEM;

};



Ntp1Analyzer_HZZlljj::Ntp1Analyzer_HZZlljj( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "HZZlljj", dataset, flags, tree ) {


  //nothing to do here


} //constructor



void Ntp1Analyzer_HZZlljj::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");

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

  reducedTree_->Branch("eLept1Gen",  &eLept1Gen_,  "eLept1Gen_/F");
  reducedTree_->Branch("ptLept1Gen",  &ptLept1Gen_,  "ptLept1Gen_/F");
  reducedTree_->Branch("etaLept1Gen",  &etaLept1Gen_,  "etaLept1Gen_/F");
  reducedTree_->Branch("phiLept1Gen",  &phiLept1Gen_,  "phiLept1Gen_/F");

  reducedTree_->Branch("eLept2",  &eLept2_,  "eLept2_/F");
  reducedTree_->Branch("ptLept2",  &ptLept2_,  "ptLept2_/F");
  reducedTree_->Branch("etaLept2",  &etaLept2_,  "etaLept2_/F");
  reducedTree_->Branch("phiLept2",  &phiLept2_,  "phiLept2_/F");

  reducedTree_->Branch("eLept2Gen",  &eLept2Gen_,  "eLept2Gen_/F");
  reducedTree_->Branch("ptLept2Gen",  &ptLept2Gen_,  "ptLept2Gen_/F");
  reducedTree_->Branch("etaLept2Gen",  &etaLept2Gen_,  "etaLept2Gen_/F");
  reducedTree_->Branch("phiLept2Gen",  &phiLept2Gen_,  "phiLept2Gen_/F");

  reducedTree_->Branch("eJetLead",  &eJetLead_,  "eJetLead_/F");
  reducedTree_->Branch( "ptJetLead",  &ptJetLead_,  "ptJetLead_/F");
  reducedTree_->Branch("etaJetLead", &etaJetLead_, "etaJetLead_/F");
  reducedTree_->Branch("phiJetLead", &phiJetLead_, "phiJetLead_/F");

  reducedTree_->Branch("eJetLead2",  &eJetLead2_,  "eJetLead2_/F");
  reducedTree_->Branch( "ptJetLead2",  &ptJetLead2_,  "ptJetLead2_/F");
  reducedTree_->Branch("etaJetLead2", &etaJetLead2_, "etaJetLead2_/F");
  reducedTree_->Branch("phiJetLead2", &phiJetLead2_, "phiJetLead2_/F");

  reducedTree_->Branch("eJetLead3",  &eJetLead3_,  "eJetLead3_/F");
  reducedTree_->Branch( "ptJetLead3",  &ptJetLead3_,  "ptJetLead3_/F");
  reducedTree_->Branch("etaJetLead3", &etaJetLead3_, "etaJetLead3_/F");
  reducedTree_->Branch("phiJetLead3", &phiJetLead3_, "phiJetLead3_/F");

  reducedTree_->Branch("eJetRecoil",  &eJetRecoil_,  "eJetRecoil_/F");
  reducedTree_->Branch( "ptJetRecoil",  &ptJetRecoil_,  "ptJetRecoil_/F");
  reducedTree_->Branch("etaJetRecoil", &etaJetRecoil_, "etaJetRecoil_/F");
  reducedTree_->Branch("phiJetRecoil", &phiJetRecoil_, "phiJetRecoil_/F");

  reducedTree_->Branch("nPairs", &nPairs_, "nPairs_/I");

  reducedTree_->Branch("iJet1",  iJet1_,  "iJet1_[nPairs_]/I");
  reducedTree_->Branch("eJet1",  eJet1_,  "eJet1_[nPairs_]/F");
  reducedTree_->Branch( "ptJet1",  ptJet1_,  "ptJet1_[nPairs_]/F");
  reducedTree_->Branch("etaJet1", etaJet1_, "etaJet1_[nPairs_]/F");
  reducedTree_->Branch("phiJet1", phiJet1_, "phiJet1_[nPairs_]/F");

  reducedTree_->Branch("eChargedHadronsJet1", eChargedHadronsJet1_, "eChargedHadronsJet1_/F");
  reducedTree_->Branch("ePhotonsJet1", &ePhotonsJet1_, "ePhotonsJet1_/F");
  reducedTree_->Branch("eNeutralHadronsJet1", &eNeutralHadronsJet1_, "eNeutralHadronsJet1_/F");
  reducedTree_->Branch("eMuonsJet1", &eMuonsJet1_, "eMuonsJet1_/F");
  reducedTree_->Branch("eElectronsJet1", &eElectronsJet1_, "eElectronsJet1_/F");
  reducedTree_->Branch("eHFHadronsJet1", &eHFHadronsJet1_, "eHFHadronsJet1_/F");
  reducedTree_->Branch("eHFEMJet1", &eHFEMJet1_, "eHFEMJet1_/F");

  reducedTree_->Branch("nChargedHadronsJet1", &nChargedHadronsJet1_, "nChargedHadronsJet1_/I");
  reducedTree_->Branch("nPhotonsJet1", &nPhotonsJet1_, "nPhotonsJet1_/I");
  reducedTree_->Branch("nNeutralHadronsJet1", &nNeutralHadronsJet1_, "nNeutralHadronsJet1_/I");
  reducedTree_->Branch("nMuonsJet1", &nMuonsJet1_, "nMuonsJet1_/I");
  reducedTree_->Branch("nElectronsJet1", &nElectronsJet1_, "nElectronsJet1_/I");
  reducedTree_->Branch("nHFHadronsJet1", &nHFHadronsJet1_, "nHFHadronsJet1_/I");
  reducedTree_->Branch("nHFEMJet1", &nHFEMJet1_, "nHFEMJet1_/I");

  reducedTree_->Branch("nPFCand1",  &nPFCand1_,  "nPFCand1_/I");
  reducedTree_->Branch("ePFCand1",  &ePFCand1_,  "ePFCand1_[nPFCand1_]/F");
  reducedTree_->Branch("ptPFCand1",  &ptPFCand1_,  "ptPFCand1_[nPFCand1_]/F");
  reducedTree_->Branch("etaPFCand1",  &etaPFCand1_,  "etaPFCand1_[nPFCand1_]/F");
  reducedTree_->Branch("phiPFCand1",  &phiPFCand1_,  "phiPFCand1_[nPFCand1_]/F");
  reducedTree_->Branch("particleTypePFCand1",  &particleTypePFCand1_,  "particleTypePFCand1_[nPFCand1_]/I");

//reducedTree_->Branch(   "eJetGen1",    &eJetGen1_,    "eJetGen1_/F");
//reducedTree_->Branch(  "ptJetGen1",   &ptJetGen1_,   "ptJetGen1_/F");
//reducedTree_->Branch( "etaJetGen1",  &etaJetGen1_,  "etaJetGen1_/F");
//reducedTree_->Branch( "phiJetGen1",  &phiJetGen1_,  "phiJetGen1_/F");
//reducedTree_->Branch("partIdJetGen1", &partIdJetGen1_, "partIdJetGen1_/I");

//reducedTree_->Branch(   "ePart1",    &ePart1_,    "ePart1_/F");
//reducedTree_->Branch(  "ptPart1",   &ptPart1_,   "ptPart1_/F");
//reducedTree_->Branch( "etaPart1",  &etaPart1_,  "etaPart1_/F");
//reducedTree_->Branch( "phiPart1",  &phiPart1_,  "phiPart1_/F");

  reducedTree_->Branch("iJet2",  iJet2_,  "iJet2_[nPairs_]/I");
  reducedTree_->Branch("eJet2",  eJet2_,  "eJet2_[nPairs_]/F");
  reducedTree_->Branch( "ptJet2",  ptJet2_,  "ptJet2_[nPairs_]/F");
  reducedTree_->Branch("etaJet2", etaJet2_, "etaJet2_[nPairs_]/F");
  reducedTree_->Branch("phiJet2", phiJet2_, "phiJet2_[nPairs_]/F");

  reducedTree_->Branch("eChargedHadronsJet2", &eChargedHadronsJet2_, "eChargedHadronsJet2_/F");
  reducedTree_->Branch("ePhotonsJet2", &ePhotonsJet2_, "ePhotonsJet2_/F");
  reducedTree_->Branch("eNeutralHadronsJet2", &eNeutralHadronsJet2_, "eNeutralHadronsJet2_/F");
  reducedTree_->Branch("eMuonsJet2", &eMuonsJet2_, "eMuonsJet2_/F");
  reducedTree_->Branch("eElectronsJet2", &eElectronsJet2_, "eElectronsJet2_/F");
  reducedTree_->Branch("eHFHadronsJet2", &eHFHadronsJet2_, "eHFHadronsJet2_/F");
  reducedTree_->Branch("eHFEMJet2", &eHFEMJet2_, "eHFEMJet2_/F");

  reducedTree_->Branch("nChargedHadronsJet2", &nChargedHadronsJet2_, "nChargedHadronsJet2_/I");
  reducedTree_->Branch("nPhotonsJet2", &nPhotonsJet2_, "nPhotonsJet2_/I");
  reducedTree_->Branch("nNeutralHadronsJet2", &nNeutralHadronsJet2_, "nNeutralHadronsJet2_/I");
  reducedTree_->Branch("nMuonsJet2", &nMuonsJet2_, "nMuonsJet2_/I");
  reducedTree_->Branch("nElectronsJet2", &nElectronsJet2_, "nElectronsJet2_/I");
  reducedTree_->Branch("nHFHadronsJet2", &nHFHadronsJet2_, "nHFHadronsJet2_/I");
  reducedTree_->Branch("nHFEMJet2", &nHFEMJet2_, "nHFEMJet2_/I");

  reducedTree_->Branch("nPFCand2",  &nPFCand2_,  "nPFCand2_/I");
  reducedTree_->Branch("ePFCand2",  &ePFCand2_,  "ePFCand2_[nPFCand2_]/F");
  reducedTree_->Branch("ptPFCand2",  &ptPFCand2_,  "ptPFCand2_[nPFCand2_]/F");
  reducedTree_->Branch("etaPFCand2",  &etaPFCand2_,  "etaPFCand2_[nPFCand2_]/F");
  reducedTree_->Branch("phiPFCand2",  &phiPFCand2_,  "phiPFCand2_[nPFCand2_]/F");
  reducedTree_->Branch("particleTypePFCand2",  &particleTypePFCand2_,  "particleTypePFCand2_[nPFCand2_]/I");

//reducedTree_->Branch(   "eJetGen2",    &eJetGen2_,    "eJetGen2_/F");
//reducedTree_->Branch(  "ptJetGen2",   &ptJetGen2_,   "ptJetGen2_/F");
//reducedTree_->Branch( "etaJetGen2",  &etaJetGen2_,  "etaJetGen2_/F");
//reducedTree_->Branch( "phiJetGen2",  &phiJetGen2_,  "phiJetGen2_/F");
//reducedTree_->Branch("partIdJetGen2", &partIdJetGen2_, "partIdJetGen2_/I");

//reducedTree_->Branch(   "ePart2",    &ePart2_,    "ePart2_/F");
//reducedTree_->Branch(  "ptPart2",   &ptPart2_,   "ptPart2_/F");
//reducedTree_->Branch( "etaPart2",  &etaPart2_,  "etaPart2_/F");
//reducedTree_->Branch( "phiPart2",  &phiPart2_,  "phiPart2_/F");

  reducedTree_->Branch("nPart", &nPart_, "nPart_/I");
  reducedTree_->Branch("ePart",  ePart_,  "ePart_[nPart_]/F");
  reducedTree_->Branch( "ptPart",  ptPart_,  "ptPart_[nPart_]/F");
  reducedTree_->Branch("etaPart", etaPart_, "etaPart_[nPart_]/F");
  reducedTree_->Branch("phiPart", phiPart_, "phiPart_[nPart_]/F");
  reducedTree_->Branch("pdgIdPart", pdgIdPart_, "pdgIdPart_[nPart_]/I");


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


   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();


   Long64_t nbytes = 0, nb = 0;

   TRandom3 rand;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%100000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;

     run_ = runNumber;
     LS_ = lumiBlock;
     event_ = eventNumber;
     eventWeight_ = -1.; //default

     if( !isGoodEvent() ) continue; //this takes care also of integrated luminosity and trigger

     //trigger:
     // not yet

     bool isMC = ( runNumber < 5 );


     ptHat_ = (isMC) ? genPtHat : ptHat_;

     if( isMC ) 
       if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;


     bool noLeptons = false;
     TLorentzVector lept1MC, lept2MC;
     int zIndexqq=-1;
     int zIndexll=-1;

     if( isMC ) {


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

         if( statusMc[iMc] != 1 ) continue;

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

         // remember: a stable lepton is daughter of a parton lepton, which is daughter of the Z:
         if( idMc[mothMc[mothMc[iMc]]]==23 ) {
           zIndexll = mothMc[mothMc[iMc]]; 
           if( fabs(idMc[iMc])==11 && idMc[mothMc[mothMc[iMc]]]==23 ) electronsMC.push_back( *thisParticle );
           if( fabs(idMc[iMc])==13 && idMc[mothMc[mothMc[iMc]]]==23 ) muonsMC.push_back( *thisParticle );
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
     // ELECTRONS
     // ------------------

     std::vector<TLorentzVector> electrons;
     int chargeFirstEle = 0;
     bool firstPassedVBTF80 = false;

     for( unsigned int iEle=0; (iEle<nEle) && (electrons.size()<2); ++iEle ) {

       TLorentzVector thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );

       // --------------
       // kinematics:
       // --------------
       if( thisEle.Pt() < 20. ) continue;
       if( (fabs(thisEle.Eta()) > 2.5) || ( fabs(thisEle.Eta())>1.4442 && fabs(thisEle.Eta())<1.566) ) continue;


       Float_t dr03TkSumPt_thresh95;
       Float_t dr03EcalRecHitSumEt_thresh95;
       Float_t dr03HcalTowerSumEt_thresh95;
       Float_t combinedIsoRel_thresh95;
       Float_t sigmaIetaIeta_thresh95;
       Float_t deltaPhiAtVtx_thresh95;
       Float_t deltaEtaAtVtx_thresh95;
       Float_t hOverE_thresh95;

       Float_t dr03TkSumPt_thresh80;
       Float_t dr03EcalRecHitSumEt_thresh80;
       Float_t dr03HcalTowerSumEt_thresh80;
       Float_t combinedIsoRel_thresh80;
       Float_t sigmaIetaIeta_thresh80;
       Float_t deltaPhiAtVtx_thresh80;
       Float_t deltaEtaAtVtx_thresh80;
       Float_t hOverE_thresh80;

       if( fabs(thisEle.Eta())<1.4442 ) {
         dr03TkSumPt_thresh95 = 0.15;
         dr03EcalRecHitSumEt_thresh95 = 2.;
         dr03HcalTowerSumEt_thresh95 = 0.12;
         combinedIsoRel_thresh95 = 0.15;

         dr03TkSumPt_thresh80 = 0.09;
         dr03EcalRecHitSumEt_thresh80 = 0.07;
         dr03HcalTowerSumEt_thresh80 = 0.10;
         combinedIsoRel_thresh80 = 0.07;

         sigmaIetaIeta_thresh95 = 0.01;
         deltaPhiAtVtx_thresh95 = 0.8;
         deltaEtaAtVtx_thresh95 = 0.007;
         hOverE_thresh95 = 0.15;

         sigmaIetaIeta_thresh80 = 0.01;
         deltaPhiAtVtx_thresh80 = 0.06;
         deltaEtaAtVtx_thresh80 = 0.004;
         hOverE_thresh80 = 0.04;
       } else {
         dr03TkSumPt_thresh95 = 0.08;
         dr03EcalRecHitSumEt_thresh95 = 0.06;
         dr03HcalTowerSumEt_thresh95 = 0.05;
         combinedIsoRel_thresh95 = 0.1;

         dr03TkSumPt_thresh80 = 0.04;
         dr03EcalRecHitSumEt_thresh80 = 0.05;
         dr03HcalTowerSumEt_thresh80 = 0.025;
         combinedIsoRel_thresh80 = 0.06;

         sigmaIetaIeta_thresh80 = 0.03;
         deltaPhiAtVtx_thresh80 = 0.7;
         deltaEtaAtVtx_thresh80 = 0.007; //no cut
         hOverE_thresh80 = 0.025;

         sigmaIetaIeta_thresh95 = 0.03;
         deltaPhiAtVtx_thresh95 = 0.7;
         deltaEtaAtVtx_thresh95 = 0.01; //no cut
         hOverE_thresh95 = 0.07;
       }


       // --------------
       // isolation:
       // --------------
     //// no relative iso, using combined
     //if( dr03TkSumPtEle[iEle]/thisEle.Pt() > dr03TkSumPt_thresh ) continue;
     //if( dr03EcalRecHitSumEtEle[iEle]/thisEle.Pt() > dr03EcalRecHitSumEt_thresh ) continue;
     //if( dr03HcalTowerSumEtEle[iEle]/thisEle.Pt() > dr03HcalTowerSumEt_thresh ) continue;

       Float_t combinedIsoRel;
       if( fabs(thisEle.Eta())<1.4442 )
         combinedIsoRel = ( dr03TkSumPtEle[iEle] + TMath::Max(0., dr03EcalRecHitSumEtEle[iEle] - 1.) + dr03HcalTowerSumEtEle[iEle] ) / thisEle.Pt();
       else
         combinedIsoRel = ( dr03TkSumPtEle[iEle] + dr03EcalRecHitSumEtEle[iEle] + dr03HcalTowerSumEtEle[iEle] ) / thisEle.Pt();

       bool iso_VBTF95 = (combinedIsoRel < combinedIsoRel_thresh95);
       bool iso_VBTF80 = (combinedIsoRel < combinedIsoRel_thresh80);

       
       // --------------
       // electron ID:
       // --------------
       bool eleID_VBTF95 = (covIEtaIEtaSC[iEle] < sigmaIetaIeta_thresh95) &&
                           (fabs(deltaPhiAtVtxEle[iEle]) < deltaPhiAtVtx_thresh95) &&
                           (fabs(deltaEtaAtVtxEle[iEle]) < deltaEtaAtVtx_thresh95) &&
                           (hOverEEle[iEle] < hOverE_thresh95);

       bool eleID_VBTF80 = (covIEtaIEtaSC[iEle] < sigmaIetaIeta_thresh80) &&
                           (fabs(deltaPhiAtVtxEle[iEle]) < deltaPhiAtVtx_thresh80) &&
                           (fabs(deltaEtaAtVtxEle[iEle]) < deltaEtaAtVtx_thresh80) &&
                           (hOverEEle[iEle] < hOverE_thresh80);

       bool passed_VBTF95 = (iso_VBTF95 && eleID_VBTF95);
       bool passed_VBTF80 = (iso_VBTF80 && eleID_VBTF80);



       // for now simple selection, will have to optimize this (T&P?)
       // one electron required to pass VBTF80, the other VBTF95
       if( electrons.size()==0 && passed_VBTF95 ) {
         electrons.push_back( thisEle );
         chargeFirstEle = chargeEle[iEle];
         if( passed_VBTF80 ) firstPassedVBTF80 = true;
       } else if( chargeEle[iEle] != chargeFirstEle && ( (firstPassedVBTF80&&passed_VBTF95)||passed_VBTF80 ) ) {
         electrons.push_back( thisEle );
       }


     } //for electrons


     // ------------------
     // MUONS
     // ------------------

     std::vector<TLorentzVector> muons;
     int chargeFirstMuon;

     for( unsigned int iMuon=0; iMuon<nMuon && (muons.size()<2); ++iMuon ) {

       TLorentzVector thisMuon( pxMuon[iMuon], pyMuon[iMuon], pzMuon[iMuon], energyMuon[iMuon] );

       // --------------
       // kinematics:
       // --------------
       if( thisMuon.Pt() < 20. ) continue;


       // --------------
       // ID:
       // --------------
       if( !( (muonIdMuon[iMuon]>>8)&1 ) ) continue; //GlobalMuonPromptTight
       if( !( (muonIdMuon[iMuon]>>11)&1 ) ) continue; //AllTrackerMuons
       if( pixelHitsTrack[trackIndexMuon[iMuon]]==0 ) continue;


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

       if( dxy > 0.02 ) continue;


       float dz = fabs(trackVzTrack[trackIndexMuon[iMuon]]-PVzPV[hardestPV]);
       if(dz > 1.0) continue;



       // --------------
       // isolation:
       // --------------
       // (this is sum pt tracks)
       //if( sumPt03Muon[iMuon] >= 3. ) continue;
       // combined isolation < 15%:
       if( (sumPt03Muon[iMuon] + emEt03Muon[iMuon] + hadEt03Muon[iMuon]) >= 0.15*thisMuon.Pt() ) continue;



       // for now simple selection, will have to optimize this (T&P?)
       if( muons.size()==0 ) {
         muons.push_back( thisMuon );
         chargeFirstMuon = chargeMuon[iMuon];
       } else {
         if( chargeMuon[iMuon]==chargeFirstMuon ) continue;
         if( fabs(muons[0].Eta())>2.1 && fabs(thisMuon.Eta())>2.1 ) continue;
         muons.push_back(thisMuon);
       }

     } //for muons


     if( electrons.size() < 2 && muons.size() < 2 ) continue;


     std::vector< TLorentzVector > leptons;

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
     
     eLept2_ = leptons[1].Energy();
     ptLept2_ = leptons[1].Pt();
     etaLept2_ = leptons[1].Eta();
     phiLept2_ = leptons[1].Phi();


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

     for( unsigned int iJet=0; iJet<nAK5PFJet; ++iJet ) {

       AnalysisJet thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );

       // save at least 3 lead jets (if event has them) and all jets with pt>thresh:
       if( leadJets.size()>=3 && thisJet.Pt()<jetPt_thresh ) break;

       // far away from leptons:
       if( thisJet.DeltaR( leptons[0] ) <= 0.5 ) continue;
       if( thisJet.DeltaR( leptons[1] ) <= 0.5 ) continue;

       thisJet.eChargedHadrons = chargedHadronEnergyAK5PFJet[iJet];

       leadJets.push_back(thisJet);
       leadJetsIndex.push_back(iJet);

     }

     if( leadJets.size()<2 ) continue;
     if( leadJets[1].Pt()<jetPt_thresh ) continue; //at least 2 jets over thresh



     // now look for best invariant mass jet pair 
     float Zmass = 91.19;
     float bestMass = 0.;
     int best_i=-1;
     int best_j=-1;
     int best_i_eventIndex=-1;
     int best_j_eventIndex=-1;

     nPairs_ = 0;
     nPart_ = 0;


//std::cout << std::endl;
     for( unsigned iJet=0; iJet<leadJets.size(); ++iJet ) {
   
       AnalysisJet thisJet = leadJets[iJet];

       // --------------
       // kinematics:
       // --------------
       if( thisJet.Pt() < jetPt_thresh ) continue;

//std::cout << "First jet: pt: " << thisJet.Pt() << " eta: " << thisJet.Eta() << " phi: " << thisJet.Phi();

       for( unsigned int jJet=iJet+1; jJet<leadJets.size(); ++jJet ) {

         AnalysisJet otherJet = leadJets[jJet];

         // --------------
         // kinematics:
         // --------------
         if( otherJet.Pt() < jetPt_thresh ) continue;

//std::cout << "   Second jet: pt: " << otherJet.Pt() << " eta: " << otherJet.Eta() << " phi: " << otherJet.Phi() << std::endl;

         if( nPairs_>=50 ) {
        
           std::cout << "MORE than 50 jet pairs found. SKIPPING!!" << std::endl;

         } else {

           eJet1_[nPairs_] = leadJets[iJet].Energy();
           ptJet1_[nPairs_] = leadJets[iJet].Pt();
           etaJet1_[nPairs_] = leadJets[iJet].Eta();
           phiJet1_[nPairs_] = leadJets[iJet].Phi();
           eChargedHadronsJet1_[nPairs_] = leadJets[iJet].eChargedHadrons;
            
           eJet2_[nPairs_] = leadJets[jJet].Energy();
           ptJet2_[nPairs_] = leadJets[jJet].Pt();
           etaJet2_[nPairs_] = leadJets[jJet].Eta();
           phiJet2_[nPairs_] = leadJets[jJet].Phi();
           eChargedHadronsJet2_[nPairs_] = leadJets[jJet].eChargedHadrons;

           nPairs_++;
          
         }


         TLorentzVector dijet = thisJet + otherJet;
         float invMass = dijet.M();
         if( (best_i==-1 && best_j==-1 ) || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
           bestMass = invMass;
           best_i = iJet;
           best_j = jJet;
           best_i_eventIndex = leadJetsIndex[iJet];
           best_j_eventIndex = leadJetsIndex[jJet];

         }
       } //for j
     } //for i

     
     //if( jets.size() < 2 ) continue;
     if( best_i==-1 || best_j==-1 ) continue; //means that less than 2 jets were found


     
     // look for hardest jet in event which is not already picked as best-Z pair (but no pt cut):
     TLorentzVector recoilJet(0., 0., 0., 0.);
     for( unsigned int iJet=0; iJet<leadJets.size() && recoilJet.Energy()==0.; ++iJet ) {

       if( iJet==best_i || iJet==best_j ) {

         continue;

       } else {

         recoilJet = leadJets[iJet];

       }

     }

         

     eJetRecoil_ = (recoilJet.Energy()==0.) ? 0. : recoilJet.Energy();
     ptJetRecoil_ = (recoilJet.Energy()==0.) ? 0. : recoilJet.Pt();
     etaJetRecoil_ = (recoilJet.Energy()==0.) ? 20. : recoilJet.Eta();
     phiJetRecoil_ = (recoilJet.Energy()==0.) ? 0. : recoilJet.Phi();
   //eChargedHadronsJetRecoil_ = recoilJet.eChargedHadrons;
   //eNeutralHadronsJetRecoil_ = recoilJet.eNeutralHadrons;
   //ePhotonsJetRecoil_ = recoilJet.ePhotons;
   //eElectronsJetRecoil_ = recoilJet.eElectrons;
   
//     TLorentzVector leadJet(pxAK5PFJet[iLeadJet], pyAK5PFJet[iLeadJet], pzAK5PFJet[iLeadJet], energyAK5PFJet[iLeadJet] );

     eJetLead_ = leadJets[0].Energy();
     ptJetLead_ = leadJets[0].Pt();
     etaJetLead_ = leadJets[0].Eta();
     phiJetLead_ = leadJets[0].Phi();

     eJetLead2_ = leadJets[1].Energy();
     ptJetLead2_ = leadJets[1].Pt();
     etaJetLead2_ = leadJets[1].Eta();
     phiJetLead2_ = leadJets[1].Phi();

     eJetLead3_ = (leadJets.size()>2) ? leadJets[2].Energy() : 0.;
     ptJetLead3_ = (leadJets.size()>2) ? leadJets[2].Pt() : 0.;
     etaJetLead3_ = (leadJets.size()>2) ? leadJets[2].Eta() : 0.;
     phiJetLead3_ = (leadJets.size()>2) ? leadJets[2].Phi() : 0.;

//   iJet1_ = best_i;
//   eJet1_ = leadJets[best_i].Energy();
//   ptJet1_ = leadJets[best_i].Pt();
//   etaJet1_ = leadJets[best_i].Eta();
//   phiJet1_ = leadJets[best_i].Phi();
//   eChargedHadronsJet1_ = leadJets[best_i].eChargedHadrons;
//  // ePhotonsJet1_ = leadJets[best_i];


//   iJet2_ = best_j;
//   eJet2_ = leadJets[best_j].Energy();
//   ptJet2_ = leadJets[best_j].Pt();
//   etaJet2_ = leadJets[best_j].Eta();
//   phiJet2_ = leadJets[best_j].Phi();
//   eChargedHadronsJet2_ = leadJets[best_j].eChargedHadrons;
//  // ePhotonsJet2_ = leadJets[best_j];


/*
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
     
*/



     // store event partons in tree:
     for( unsigned iMc=0; iMc<nMc; ++iMc ) {

       if( statusMc[iMc]==3 && (fabs(idMc[iMc])<=6 || idMc[iMc]==21) ) {

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

         if( nPart_<20 ) {

           ptPart_[nPart_] = thisParticle->Pt();
           etaPart_[nPart_] = thisParticle->Eta();
           phiPart_[nPart_] = thisParticle->Phi();
           ePart_[nPart_] = thisParticle->Energy();
           pdgIdPart_[nPart_] = idMc[iMc];

           nPart_++;

         } else {
    
           std::cout << "Found more than 20 partons, skipping." << std::endl;

         }

         delete thisParticle;
         thisParticle = 0;

       } //if correct id mc

     } // for i mc


     reducedTree_->Fill(); 


   } //for entries

} //loop



double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}

