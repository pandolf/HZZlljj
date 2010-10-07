#include "Ntp1Analyzer_JetStudies.h"


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



Ntp1Analyzer_JetStudies::Ntp1Analyzer_JetStudies( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "JetStudies", dataset, flags, tree ) {


  //nothing to do here


} //constructor



void Ntp1Analyzer_JetStudies::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("nJet",  &nJet_,  "nJet_/I");

  reducedTree_->Branch("eJet",  &eJet_,  "eJet_[nJet_]/F");
  reducedTree_->Branch("ptJet",  &ptJet_,  "ptJet_[nJet_]/F");
  reducedTree_->Branch("etaJet",  &etaJet_,  "etaJet_[nJet_]/F");
  reducedTree_->Branch("phiJet",  &phiJet_,  "phiJet_[nJet_]/F");


} 



Ntp1Analyzer_JetStudies::~Ntp1Analyzer_JetStudies() {

  outfile_->cd();
  

}



void Ntp1Analyzer_JetStudies::Loop()
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

     if( !isGoodEvent() ) continue; //this takes care also of integrated luminosity and trigger

     //trigger:
     // not yet

     bool isMC = ( runNumber < 5 );


     if( isMC ) 
       if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;


     nJet_=0;

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

       if( electrons[0].Pt() > electrons[1].Pt() ) {

         leptons.push_back( electrons[0] );
         leptons.push_back( electrons[1] );

       } else {

         leptons.push_back( electrons[1] );
         leptons.push_back( electrons[0] );

       }

     } else if( muons.size() == 2 ) {

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


     float eLept1_ = leptons[0].Energy();
     float ptLept1_ = leptons[0].Pt();
     float etaLept1_ = leptons[0].Eta();
     float phiLept1_ = leptons[0].Phi();
     
     float eLept2_ = leptons[1].Energy();
     float ptLept2_ = leptons[1].Pt();
     float etaLept2_ = leptons[1].Eta();
     float phiLept2_ = leptons[1].Phi();




     // ------------------
     // JETS
     // ------------------

     std::vector<AnalysisJet> jets;

     for( unsigned int iJet=0; iJet<nAK5PFJet && iJet<5; ++iJet ) {

       AnalysisJet thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );

       // --------------
       // kinematics:
       // --------------
       if( thisJet.Pt() < 20. ) continue;

       // far away from leptons:
       Float_t deltaEta1 = thisJet.Eta() - etaLept1_;
       Float_t deltaPhi1 = delta_phi((Float_t)thisJet.Phi(), phiLept1_);
       Float_t deltaR1 = sqrt( deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1 );
       if( deltaR1 <= 0.25 ) continue;

       Float_t deltaEta2 = thisJet.Eta() - etaLept2_;
       Float_t deltaPhi2 = delta_phi((Float_t)thisJet.Phi(), phiLept2_);
       Float_t deltaR2 = sqrt( deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2 );
       if( deltaR2 <= 0.25 ) continue;

       thisJet.eChargedHadrons = chargedHadronEnergyAK5PFJet[iJet];

       ptJet_[nJet_] = thisJet.Pt();
       etaJet_[nJet_] = thisJet.Eta();
       phiJet_[nJet_] = thisJet.Phi();
       eJet_[nJet_] = thisJet.Energy();
       nJet_++;

     } //for jets

    
     reducedTree_->Fill(); 


   } //for entries

} //loop



