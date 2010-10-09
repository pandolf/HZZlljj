#include "Ntp1Analyzer_HZZlljj.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"


double delta_phi(double phi1, double phi2) {

  double dphi = fabs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}


float delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}




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

  reducedTree_->Branch("eJetRecoil",  &eJetRecoil_,  "eJetRecoil_/F");
  reducedTree_->Branch( "ptJetRecoil",  &ptJetRecoil_,  "ptJetRecoil_/F");
  reducedTree_->Branch("etaJetRecoil", &etaJetRecoil_, "etaJetRecoil_/F");
  reducedTree_->Branch("phiJetRecoil", &phiJetRecoil_, "phiJetRecoil_/F");

  reducedTree_->Branch("iJet1",  &iJet1_,  "iJet1_/I");
  reducedTree_->Branch("eJet1",  &eJet1_,  "eJet1_/F");
  reducedTree_->Branch( "ptJet1",  &ptJet1_,  "ptJet1_/F");
  reducedTree_->Branch("etaJet1", &etaJet1_, "etaJet1_/F");
  reducedTree_->Branch("phiJet1", &phiJet1_, "phiJet1_/F");

  reducedTree_->Branch("eChargedHadronsJet1", &eChargedHadronsJet1_, "eChargedHadronsJet1_/F");
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

  reducedTree_->Branch(   "eJetGen1",    &eJetGen1_,    "eJetGen1_/F");
  reducedTree_->Branch(  "ptJetGen1",   &ptJetGen1_,   "ptJetGen1_/F");
  reducedTree_->Branch( "etaJetGen1",  &etaJetGen1_,  "etaJetGen1_/F");
  reducedTree_->Branch( "phiJetGen1",  &phiJetGen1_,  "phiJetGen1_/F");
  reducedTree_->Branch("partIdJetGen1", &partIdJetGen1_, "partIdJetGen1_/I");

  reducedTree_->Branch(   "ePart1",    &ePart1_,    "ePart1_/F");
  reducedTree_->Branch(  "ptPart1",   &ptPart1_,   "ptPart1_/F");
  reducedTree_->Branch( "etaPart1",  &etaPart1_,  "etaPart1_/F");
  reducedTree_->Branch( "phiPart1",  &phiPart1_,  "phiPart1_/F");

  reducedTree_->Branch("iJet2",  &iJet2_,  "iJet2_/I");
  reducedTree_->Branch("eJet2",  &eJet2_,  "eJet2_/F");
  reducedTree_->Branch( "ptJet2",  &ptJet2_,  "ptJet2_/F");
  reducedTree_->Branch("etaJet2", &etaJet2_, "etaJet2_/F");
  reducedTree_->Branch("phiJet2", &phiJet2_, "phiJet2_/F");

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

  reducedTree_->Branch(   "eJetGen2",    &eJetGen2_,    "eJetGen2_/F");
  reducedTree_->Branch(  "ptJetGen2",   &ptJetGen2_,   "ptJetGen2_/F");
  reducedTree_->Branch( "etaJetGen2",  &etaJetGen2_,  "etaJetGen2_/F");
  reducedTree_->Branch( "phiJetGen2",  &phiJetGen2_,  "phiJetGen2_/F");
  reducedTree_->Branch("partIdJetGen2", &partIdJetGen2_, "partIdJetGen2_/I");

  reducedTree_->Branch(   "ePart2",    &ePart2_,    "ePart2_/F");
  reducedTree_->Branch(  "ptPart2",   &ptPart2_,   "ptPart2_/F");
  reducedTree_->Branch( "etaPart2",  &etaPart2_,  "etaPart2_/F");
  reducedTree_->Branch( "phiPart2",  &phiPart2_,  "phiPart2_/F");


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
  h1_indexMatchedJet = new TH1F("indexMatchedJet", "", 6, -0.5, 5.5);
  h1_indexMatched05Jet = new TH1F("indexMatched05Jet", "", 6, -0.5, 5.5);
  h1_nMatched_per_event = new TH1F("nMatched_per_event", "", 6, -0.5, 5.5);
  h1_nMatched05_per_event = new TH1F("nMatched05_per_event", "", 6, -0.5, 5.5);
  h1_pdgIdParton1 = new TH1F("pdgIdParton1", "", 36, -10.5, 25.5);
  h1_pdgIdParton2 = new TH1F("pdgIdParton2", "", 36, -10.5, 25.5);
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
  h1_indexMatchedJet->Write();
  h1_indexMatched05Jet->Write();
  h1_nMatched_per_event->Write();
  h1_nMatched05_per_event->Write();
  h1_pdgIdParton1->Write();
  h1_pdgIdParton2->Write();
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
     ptHat_ = genPtHat;
     eventWeight_ = -1.; //default

     if( !isGoodEvent() ) continue; //this takes care also of integrated luminosity and trigger

     //trigger:
     // not yet

     bool isMC = ( runNumber < 5 );


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

     if( !noLeptons && lept1MC.Pt() < lept2MC.Pt() ) std::cout << "WARNING MC leptons not ordered in pt!!" << std::endl;



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

     for( unsigned int iEle=0; (iEle<nEle) && (electrons.size()<2); ++iEle ) {

       TLorentzVector thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );

       // --------------
       // kinematics:
       // --------------
       if( thisEle.Pt() < 10. ) continue;
       if( (fabs(thisEle.Eta()) > 2.5) || ( fabs(thisEle.Eta())>1.4442 && fabs(thisEle.Eta())<1.566) ) continue;


       Float_t dr03TkSumPt_thresh;
       Float_t dr03EcalRecHitSumEt_thresh;
       Float_t dr03HcalTowerSumEt_thresh;
       Float_t combinedIsoRel_thresh;
       Float_t sigmaIetaIeta_thresh;
       Float_t deltaPhiAtVtx_thresh;
       Float_t deltaEtaAtVtx_thresh;
       Float_t hOverE_thresh;

       if( fabs(thisEle.Eta())<1.4442 ) {
         dr03TkSumPt_thresh = 0.15;
         dr03EcalRecHitSumEt_thresh = 2.;
         dr03HcalTowerSumEt_thresh = 0.12;
         combinedIsoRel_thresh = 0.15;

         sigmaIetaIeta_thresh = 0.01;
         deltaPhiAtVtx_thresh = 0.8;
         deltaEtaAtVtx_thresh = 0.007;
         hOverE_thresh = 0.15;
       } else {
         dr03TkSumPt_thresh = 0.08;
         dr03EcalRecHitSumEt_thresh = 0.06;
         dr03HcalTowerSumEt_thresh = 0.05;
         combinedIsoRel_thresh = 0.1;

         sigmaIetaIeta_thresh = 0.03;
         deltaPhiAtVtx_thresh = 0.7;
         deltaEtaAtVtx_thresh = 10.; //no cut
         hOverE_thresh = 0.07;
       }


       // --------------
       // isolation:
       // --------------
       if( dr03TkSumPtEle[iEle]/thisEle.Pt() > dr03TkSumPt_thresh ) continue;
       if( dr03EcalRecHitSumEtEle[iEle]/thisEle.Pt() > dr03EcalRecHitSumEt_thresh ) continue;
       if( dr03HcalTowerSumEtEle[iEle]/thisEle.Pt() > dr03HcalTowerSumEt_thresh ) continue;

       Float_t combinedIsoRel;
       if( fabs(thisEle.Eta())<1.4442 )
         combinedIsoRel = ( dr03TkSumPtEle[iEle] + TMath::Max(0., dr03EcalRecHitSumEtEle[iEle] - 1.) + dr03HcalTowerSumEtEle[iEle] ) / thisEle.Pt();
       else
         combinedIsoRel = ( dr03TkSumPtEle[iEle] + dr03EcalRecHitSumEtEle[iEle] + dr03HcalTowerSumEtEle[iEle] ) / thisEle.Pt();

       if( combinedIsoRel > combinedIsoRel_thresh ) continue;

       
       // --------------
       // cluster shape:
       // --------------
       if( covIEtaIEtaSC[iEle] > sigmaIetaIeta_thresh ) continue;
       if( fabs(deltaPhiAtVtxEle[iEle]) > deltaPhiAtVtx_thresh ) continue;
       if( fabs(deltaEtaAtVtxEle[iEle]) > deltaEtaAtVtx_thresh ) continue;
       if( hOverEEle[iEle] > hOverE_thresh ) continue;

       // for now simple selection, will have to optimize this (T&P?)
       if( electrons.size()==0 ) {
         electrons.push_back( thisEle );
         chargeFirstEle = chargeEle[iEle];
       } else if( chargeEle[iEle] != chargeFirstEle ) {
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
       if( thisMuon.Pt() < 10. ) continue;


       // --------------
       // ID:
       // --------------
       if( !( (muonIdMuon[iMuon]>>8)&1 ) ) continue; //GlobalMuonPromptTight
       if( !( (muonIdMuon[iMuon]>>11)&1 ) ) continue; //AllTrackerMuons
       if( pixelHitsTrack[trackIndexMuon[iMuon]]==0 ) continue;
       if( transvImpactParTrack[trackIndexMuon[iMuon]] > 0.2 ) continue;    


       // --------------
       // isolation:
       // --------------
       // (this is sum pt tracks)
       if( sumPt03Muon[iMuon] >= 3. ) continue;



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


     // look for jet pair which has invariant mass closest to Z:
     float Zmass = 91.19;
     float bestMass = 0.;
     int best_i=-1;
     int best_j=-1;
     int iLeadJet=-1;
     bool firstJet=true;

     for( unsigned int iJet=0; iJet<nAK5PFJet && iJet<6; ++iJet ) {

       AnalysisJet thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );

       // --------------
       // kinematics:
       // --------------
       if( thisJet.Pt() < 30. ) continue;

       // far away from leptons:
       Float_t deltaEta1 = thisJet.Eta() - etaLept1_;
       Float_t deltaPhi1 = delta_phi((Float_t)thisJet.Phi(), phiLept1_);
       Float_t deltaR1 = sqrt( deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1 );
       if( deltaR1 <= 0.25 ) continue;

       Float_t deltaEta2 = thisJet.Eta() - etaLept2_;
       Float_t deltaPhi2 = delta_phi((Float_t)thisJet.Phi(), phiLept2_);
       Float_t deltaR2 = sqrt( deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2 );
       if( deltaR2 <= 0.25 ) continue;

       if( firstJet ) {
         iLeadJet = iJet;
         firstJet=false;
       }


       for( unsigned int jJet=iJet+1; jJet<nAK5PFJet && jJet<6; ++jJet ) {

         AnalysisJet otherJet( pxAK5PFJet[jJet], pyAK5PFJet[jJet], pzAK5PFJet[jJet], energyAK5PFJet[jJet] );

         // --------------
         // kinematics:
         // --------------
         if( otherJet.Pt() < 30. ) continue;

         // far away from leptons:
         Float_t deltaEta1 = otherJet.Eta() - etaLept1_;
         Float_t deltaPhi1 = delta_phi((Float_t)otherJet.Phi(), phiLept1_);
         Float_t deltaR1 = sqrt( deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1 );
         if( deltaR1 <= 0.25 ) continue;

         Float_t deltaEta2 = otherJet.Eta() - etaLept2_;
         Float_t deltaPhi2 = delta_phi((Float_t)otherJet.Phi(), phiLept2_);
         Float_t deltaR2 = sqrt( deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2 );
         if( deltaR2 <= 0.25 ) continue;


         TLorentzVector dijet = thisJet + otherJet;
         float invMass = dijet.M();
         if( (best_i==-1 && best_j==-1 ) || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
           bestMass = invMass;
           best_i = iJet;
           best_j = jJet;
         }
       } //for j
     } //for i

     
     //if( jets.size() < 2 ) continue;
     if( best_i==-1 || best_j==-1 || iLeadJet==-1 ) continue; //means that less than 2 jets were found


//   for( unsigned i=0; i<jets.size(); ++i ) {
//     for( unsigned j=i+1; j<jets.size(); ++j ) {
//       TLorentzVector dijet = jets[i]+jets[j];
//       float invMass = dijet.M();
//       if( (best_i==-1 && best_j==-1 ) || ( fabs(invMass-Zmass) < fabs(bestMass-Zmass) ) ) {
//         bestMass = invMass;
//         best_i = i;
//         best_j = j;
//       }
//     } //jets j
//   } //jets i

     
     // look for hardest jet in event which is not already picked as best-Z pair (but no pt cut):
     TLorentzVector recoilJet(0., 0., 0., 0.);
     for( unsigned int iJet=0; iJet<nAK5PFJet && recoilJet.Energy()==0.; ++iJet ) {

       if( iJet==best_i || iJet==best_j ) {

         continue;

       } else {

         AnalysisJet thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );
       
         // far away from leptons:
         Float_t deltaEta1 = thisJet.Eta() - etaLept1_;
         Float_t deltaPhi1 = delta_phi((Float_t)thisJet.Phi(), phiLept1_);
         Float_t deltaR1 = sqrt( deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1 );
         if( deltaR1 <= 0.25 ) continue;
       
         Float_t deltaEta2 = thisJet.Eta() - etaLept2_;
         Float_t deltaPhi2 = delta_phi((Float_t)thisJet.Phi(), phiLept2_);
         Float_t deltaR2 = sqrt( deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2 );
         if( deltaR2 <= 0.25 ) continue;
       
         recoilJet = thisJet;

       }

     }

         

     eJetRecoil_ = recoilJet.Energy();
     ptJetRecoil_ = recoilJet.Pt();
     etaJetRecoil_ = recoilJet.Eta();
     phiJetRecoil_ = recoilJet.Phi();
   //eChargedHadronsJetRecoil_ = recoilJet.eChargedHadrons;
   //eNeutralHadronsJetRecoil_ = recoilJet.eNeutralHadrons;
   //ePhotonsJetRecoil_ = recoilJet.ePhotons;
   //eElectronsJetRecoil_ = recoilJet.eElectrons;
   
     TLorentzVector leadJet(pxAK5PFJet[iLeadJet], pyAK5PFJet[iLeadJet], pzAK5PFJet[iLeadJet], energyAK5PFJet[iLeadJet] );

     eJetLead_ = leadJet.Energy();
     ptJetLead_ = leadJet.Pt();
     etaJetLead_ = leadJet.Eta();
     phiJetLead_ = leadJet.Phi();
   //eChargedHadronsJetLead_ = jets[0].eChargedHadrons;
   //eNeutralHadronsJetLead_ = jets[0].eNeutralHadrons;
   //ePhotonsJetLead_ = jets[0].ePhotons;
   //eElectronsJetLead_ = jets[0].eElectrons;

     TLorentzVector jet1(pxAK5PFJet[best_i], pyAK5PFJet[best_i], pzAK5PFJet[best_i], energyAK5PFJet[best_i] );

     iJet1_ = best_i;
     eJet1_ = jet1.Energy();
     ptJet1_ = jet1.Pt();
     etaJet1_ = jet1.Eta();
     phiJet1_ = jet1.Phi();
     eChargedHadronsJet1_ = chargedHadronEnergyAK5PFJet[best_i];
     eChargedHadronsJet1_ = neutralHadronEnergyAK5PFJet[best_i];
     ePhotonsJet1_ = neutralEmEnergyAK5PFJet[best_i];
     ePhotonsJet1_ = chargedEmEnergyAK5PFJet[best_i];


     TLorentzVector jet2(pxAK5PFJet[best_j], pyAK5PFJet[best_j], pzAK5PFJet[best_j], energyAK5PFJet[best_j] );

     iJet2_ = best_j;
     eJet2_ = jet2.Energy();
     ptJet2_ = jet2.Pt();
     etaJet2_ = jet2.Eta();
     phiJet2_ = jet2.Phi();
     eChargedHadronsJet2_ = chargedHadronEnergyAK5PFJet[best_j];
     eChargedHadronsJet2_ = neutralHadronEnergyAK5PFJet[best_j];
     ePhotonsJet2_ = neutralEmEnergyAK5PFJet[best_j];
     ePhotonsJet2_ = chargedEmEnergyAK5PFJet[best_j];

      
     
     std::vector<TLorentzVector> jets;
     jets.push_back(jet1);
     jets.push_back(jet2);


     // --------------------
     // match jets to MC:
     // --------------------
     std::vector<TLorentzVector> matchedGenJets;
     std::vector<TLorentzVector> matchedPartons;
     std::vector<int> pdgIdPartonGenJets;

     int nMatched_per_event = 0;
     int nMatched05_per_event = 0;

     for( unsigned iJet=0; iJet<jets.size(); ++iJet ) {

       float deltaRmin = 100.;
       TLorentzVector matchedPartonMC;

       //if( iJet!=0 && iJet!=1 ) continue;
       //if( iJet!=best_i && iJet!=best_j ) continue;
       
       // first match to partons (only for signal):
       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         if( statusMc[iMc]==3 && fabs(idMc[iMc])<=6 && idMc[mothMc[iMc]]==23 ) { //partons from Z

           TLorentzVector* thisParticle = new TLorentzVector();
           thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );
           float thisDeltaR = jets[iJet].DeltaR( *thisParticle );
           if( thisDeltaR < deltaRmin ) {
             deltaRmin = thisDeltaR;
             matchedPartonMC = *thisParticle;
           }

           delete thisParticle;
           thisParticle = 0;

         } //if correct id mc

       } // for i mc

       h1_deltaRmatching_jet_parton->Fill( deltaRmin );
       if( deltaRmin < 0.25 )  { //matched
         h1_indexMatchedJet->Fill( iJet );
         nMatched_per_event++;
       }
       if( deltaRmin < 0.5 )  { //matched
         h1_indexMatched05Jet->Fill( iJet );
         nMatched05_per_event++;
       }
       matchedPartons.push_back(matchedPartonMC);





       // now match to genjet (general mc case):

       float deltaRmin_genjet = 100.;
       TLorentzVector matchedGenJet;

       for( unsigned iGenJet=0; iGenJet<nAK5GenJet; ++iGenJet ) {

         TLorentzVector* thisGenJet = new TLorentzVector( pxAK5GenJet[iGenJet], pyAK5GenJet[iGenJet], pzAK5GenJet[iGenJet], energyAK5GenJet[iGenJet]);
         float thisDeltaR = jets[iJet].DeltaR( *thisGenJet );
         if( thisDeltaR < deltaRmin_genjet ) {
           deltaRmin_genjet = thisDeltaR;
           matchedGenJet = *thisGenJet;
         }

         delete thisGenJet;
         thisGenJet = 0;

       } // for i GenJet

       h1_deltaRmatching_jet_genjet->Fill( deltaRmin_genjet );
       matchedGenJets.push_back(matchedGenJet);



       // match genjet to parton:

       float deltaRmin_parton = 100.;
       int foundPdgId;

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         if( statusMc[iMc]==3 && (fabs(idMc[iMc])<=6 || idMc[iMc]==21) ) {

           TLorentzVector* thisParticle = new TLorentzVector();
           thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );
           float thisDeltaR = matchedGenJet.DeltaR( *thisParticle );
           if( thisDeltaR < deltaRmin_parton ) {
             deltaRmin_parton = thisDeltaR;
             foundPdgId = idMc[iMc];
           }

           delete thisParticle;
           thisParticle = 0;

         } //if correct id mc

       } // for i mc

       h1_deltaRmatching_genjet_parton->Fill( deltaRmin_parton );
       pdgIdPartonGenJets.push_back(foundPdgId); 


     } //for i jets

     h1_nMatched_per_event->Fill( nMatched_per_event );
     h1_nMatched05_per_event->Fill( nMatched05_per_event );

     if( nMatched_per_event==0 ) {

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         if( statusMc[mothMc[iMc]]==3 && (fabs(idMc[iMc])==11 || fabs(idMc[iMc])==13) ) {
           TLorentzVector lepton;
           lepton.SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );
           h1_deltaRmatching_jet_leptonParton->Fill( lepton.DeltaR(jets[0]) );
           h1_deltaRmatching_jet_leptonParton->Fill( lepton.DeltaR(jets[1]) );
         }

       }
     }


     int iPart1 = -1;
     int iPart2 = -1;
     float ptPartMax=0.;

     for( unsigned iMc=0; iMc<nMc; ++iMc ) {

       if( statusMc[iMc]==3 && (fabs(idMc[iMc])<=6 || idMc[iMc]==21) ) {

         if( pMc[iMc]*sin(thetaMc[iMc])>ptPartMax ) {
           iPart2 = iPart1;
           iPart1 = iMc;
         }

       }          

     }

     if( iPart1>=0 ) h1_pdgIdParton1->Fill( idMc[iPart1] );
     if( iPart2>=0 ) h1_pdgIdParton2->Fill( idMc[iPart2] );

     eJetGen1_ = matchedGenJets[0].Energy();
     ptJetGen1_ = matchedGenJets[0].Pt();
     etaJetGen1_ = matchedGenJets[0].Eta();
     phiJetGen1_ = matchedGenJets[0].Phi();
     partIdJetGen1_ = pdgIdPartonGenJets[0];

     eJetGen2_ = matchedGenJets[1].Energy();
     ptJetGen2_ = matchedGenJets[1].Pt();
     etaJetGen2_ = matchedGenJets[1].Eta();
     phiJetGen2_ = matchedGenJets[1].Phi();
     partIdJetGen2_ = pdgIdPartonGenJets[1];

     ePart1_ = matchedPartons[0].Energy();
     ptPart1_ = matchedPartons[0].Pt();
     etaPart1_ = matchedPartons[0].Eta();
     phiPart1_ = matchedPartons[0].Phi();

     ePart2_ = matchedPartons[1].Energy();
     ptPart2_ = matchedPartons[1].Pt();
     etaPart2_ = matchedPartons[1].Eta();
     phiPart2_ = matchedPartons[1].Phi();


     reducedTree_->Fill(); 


   } //for entries

} //loop



