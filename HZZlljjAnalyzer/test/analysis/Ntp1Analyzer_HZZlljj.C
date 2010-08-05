#include "Ntp1Analyzer_HZZlljj.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include "fitTools.h"



class AnalysisJet : public TLorentzVector {

 public:

  AnalysisJet( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {
    eChargedHadrons=0.;
  }

  float eChargedHadrons;
//float ePhotons;
//float eNeutralHadrons;
//float eMuons;
//float eElectrons;
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
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");

  reducedTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  reducedTree_->Branch("leptType",  &leptType_,  "leptType_/S");
  
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

  reducedTree_->Branch(  "eJet1Gen",   &eJet1Gen_,   "eJet1Gen_/F");
  reducedTree_->Branch(  "ptJet1Gen",   &ptJet1Gen_,   "ptJet1Gen_/F");
  reducedTree_->Branch( "etaJet1Gen",  &etaJet1Gen_,  "etaJet1Gen_/F");
  reducedTree_->Branch( "phiJet1Gen",  &phiJet1Gen_,  "phiJet1Gen_/F");
  reducedTree_->Branch("pdgIdPartJet1", &pdgIdPartJet1_, "pdgIdPartJet1_/I");

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

  reducedTree_->Branch(  "eJet2Gen",   &eJet2Gen_,   "eJet2Gen_/F");
  reducedTree_->Branch(  "ptJet2Gen",   &ptJet2Gen_,   "ptJet2Gen_/F");
  reducedTree_->Branch( "etaJet2Gen",  &etaJet2Gen_,  "etaJet2Gen_/F");
  reducedTree_->Branch( "phiJet2Gen",  &phiJet2Gen_,  "phiJet2Gen_/F");
  reducedTree_->Branch("pdgIdPartJet2", &pdgIdPartJet2_, "pdgIdPartJet2_/I");

  reducedTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  reducedTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");

} 



Ntp1Analyzer_HZZlljj::~Ntp1Analyzer_HZZlljj() {

  outfile_->cd();

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
     event_ = eventNumber;
     ptHat_ = genPtHat;

     if( !isGoodLS() ) continue; //this takes care also of integrated luminosity

     //trigger:
     // not yet

     bool isMC = ( runNumber < 5 );


     if( isMC ) 
       if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;


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
       if( (fabs(thisEle.Eta() > 2.5) || ( fabs(thisEle.Eta())>1.4442 && fabs(thisEle.Eta())<1.566) ) ) continue;


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

//   std::vector<TLorentzVector> muons;

//   for( unsigned int iMuon=0; iMuon<nMuon; ++iMuon ) {

//     TLorentzVector thisMuon( pxMuon[iMuon], pyMuon[iMuon], pzMuon[iMuon], energyMuon[iMuon] );

//     // --------------
//     // kinematics:
//     // --------------
//     if( thisMuon.Pt() < 10. ) continue;
//     //if( (fabs(thisMu.Eta() > 2.5) || ( fabs(thisMu.Eta())>1.4442 && fabs(thisMu.Eta())<1.566) ) ) continue;

//     // --------------
//     // isolation:
//     // --------------

//   } //for muons


     if( electrons.size() < 2 ) continue;

     eLept1_ = electrons[0].Energy();
     ptLept1_ = electrons[0].Pt();
     etaLept1_ = electrons[0].Eta();
     phiLept1_ = electrons[0].Phi();

     eLept2_ = electrons[1].Energy();
     ptLept2_ = electrons[1].Pt();
     etaLept2_ = electrons[1].Eta();
     phiLept2_ = electrons[1].Phi();


     // ------------------
     // JETS
     // ------------------

     std::vector<AnalysisJet> jets;

     for( unsigned int iJet=0; iJet<nAK5PFJet && jets.size()<2; ++iJet ) {

       AnalysisJet thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );

       // --------------
       // kinematics:
       // --------------
       if( thisJet.Pt() < 20. ) continue;

       // far away from leptons:
       Float_t deltaEta1 = thisJet.Eta() - etaLept1_;
       Float_t deltaPhi1 = fitTools::delta_phi((Float_t)thisJet.Phi(), phiLept1_);
       Float_t deltaR1 = sqrt( deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1 );
       if( deltaR1 <= 0.25 ) continue;

       Float_t deltaEta2 = thisJet.Eta() - etaLept2_;
       Float_t deltaPhi2 = fitTools::delta_phi((Float_t)thisJet.Phi(), phiLept2_);
       Float_t deltaR2 = sqrt( deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2 );
       if( deltaR2 <= 0.25 ) continue;

       thisJet.eChargedHadrons = chargedHadronEnergyAK5PFJet[iJet];

       jets.push_back( thisJet );

     } //for jets

     
     if( jets.size() < 2 ) continue;


     eJet1_ = jets[0].Energy();
     ptJet1_ = jets[0].Pt();
     etaJet1_ = jets[0].Eta();
     phiJet1_ = jets[0].Phi();
     eChargedHadronsJet1_ = jets[0].eChargedHadrons;

     eJet2_ = jets[1].Energy();
     ptJet2_ = jets[1].Pt();
     etaJet2_ = jets[1].Eta();
     phiJet2_ = jets[1].Phi();
     eChargedHadronsJet2_ = jets[1].eChargedHadrons;
     

     reducedTree_->Fill(); 


   } //for entries

} //loop



