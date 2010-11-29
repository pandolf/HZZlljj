#include "Ntp1Analyzer_ZGamma.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"



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



Ntp1Analyzer_ZGamma::Ntp1Analyzer_ZGamma( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "ZGamma", dataset, flags, tree ) {


  //nothing to do here


} //constructor



void Ntp1Analyzer_ZGamma::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");

  reducedTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  reducedTree_->Branch("eZGamma",  &eZGamma_,  "eZGamma_/F");
  reducedTree_->Branch("ptZGamma",  &ptZGamma_,  "ptZGamma_/F");
  reducedTree_->Branch("etaZGamma",  &etaZGamma_,  "etaZGamma_/F");
  reducedTree_->Branch("phiZGamma",  &phiZGamma_,  "phiZGamma_/F");

  reducedTree_->Branch(  "eZGammaGen",    &eZGammaGen_,    "eZGammaGen_/F");
  reducedTree_->Branch( "ptZGammaGen",   &ptZGammaGen_,   "ptZGammaGen_/F");
  reducedTree_->Branch("etaZGammaGen",  &etaZGammaGen_,  "etaZGammaGen_/F");
  reducedTree_->Branch("phiZGammaGen",  &phiZGammaGen_,  "phiZGammaGen_/F");

  reducedTree_->Branch("nJets_total",  &nJets_total_,  "nJets_total_/I");

  reducedTree_->Branch("nJet",  &nJet_,  "nJet_/I");
  reducedTree_->Branch("eJet",  &eJet_,  "eJet_[nJet_]/F");
  reducedTree_->Branch("ptJet",  &ptJet_,  "ptJet_[nJet_]/F");
  reducedTree_->Branch("etaJet",  &etaJet_,  "etaJet_[nJet_]/F");
  reducedTree_->Branch("phiJet",  &phiJet_,  "phiJet_[nJet_]/F");

  reducedTree_->Branch("nJetsGen_total",  &nJetsGen_total_,  "nJetsGen_total_/I");

  reducedTree_->Branch(  "nJetGen",    &nJetGen_,    "nJetGen_/I");
  reducedTree_->Branch(  "eJetGen",    &eJetGen_,    "eJetGen_[nJetGen_]/F");
  reducedTree_->Branch( "ptJetGen",   &ptJetGen_,   "ptJetGen_[nJetGen_]/F");
  reducedTree_->Branch("etaJetGen",  &etaJetGen_,  "etaJetGen_[nJetGen_]/F");
  reducedTree_->Branch("phiJetGen",  &phiJetGen_,  "phiJetGen_[nJetGen_]/F");


} 



Ntp1Analyzer_ZGamma::~Ntp1Analyzer_ZGamma() {

  outfile_->cd();

}



void Ntp1Analyzer_ZGamma::Loop()
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




     if( isMC ) {

       bool noLeptonsGen = false;
       TLorentzVector lept1MC, lept2MC;
       int zIndexll=-1;


       // look for Z->ll

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
       } else if( muonsMC.size()==2 ) {
         if( muonsMC[0].Pt() > muonsMC[1].Pt() ) {
           lept1MC = muonsMC[0];
           lept2MC = muonsMC[1];
         } else {
           lept1MC = muonsMC[1];
           lept2MC = muonsMC[0];
         }
       } else {
         //taus
         noLeptonsGen = true;
       }


       if( !noLeptonsGen ) {

         TLorentzVector ZGammaGen;
         ZGammaGen.SetPtEtaPhiE( pMc[zIndexll]*sin(thetaMc[zIndexll]), etaMc[zIndexll], phiMc[zIndexll], energyMc[zIndexll] );

         ptZGammaGen_  = ZGammaGen.Pt();
         eZGammaGen_   = ZGammaGen.Energy();
         etaZGammaGen_ = ZGammaGen.Eta();
         phiZGammaGen_ = ZGammaGen.Phi();

         nJetGen_ = nJetsGen_total_ = 0;

         std::vector<AnalysisJet> jets;

         for( unsigned int iJet=0; iJet<nAK5GenJet; ++iJet ) {

           AnalysisJet thisJet( pxAK5GenJet[iJet], pyAK5GenJet[iJet], pzAK5GenJet[iJet], energyAK5GenJet[iJet] );

           // far away from leptons:
           if( thisJet.DeltaR(lept1MC) <= 0.5 ) continue;
           if( thisJet.DeltaR(lept2MC) <= 0.5 ) continue;

           nJetsGen_total_++;

           if( nJetGen_<10 && thisJet.Pt()>30. ) {
             ptJetGen_[nJetGen_] = thisJet.Pt();
             etaJetGen_[nJetGen_] = thisJet.Eta();
             phiJetGen_[nJetGen_] = thisJet.Phi();
             eJetGen_[nJetGen_] = thisJet.Energy();
             nJetGen_++;
           }

         } //for jets

       } // if !noleptonsgen

     } // if is MC



     // and now on RECO:

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
         deltaEtaAtVtx_thresh80 = 0.007;
         hOverE_thresh80 = 0.025;

         sigmaIetaIeta_thresh95 = 0.03;
         deltaPhiAtVtx_thresh95 = 0.7;
         deltaEtaAtVtx_thresh95 = 0.01; 
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


       if( !passed_VBTF95 ) continue;

       // check that not matched to muon (clean electrons faked by muon MIP):
       bool matchedtomuon=false;
       for( std::vector<TLorentzVector>::iterator iMu=muons.begin(); iMu!=muons.end(); ++iMu )
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



     bool noLeptons=false;
     if( electrons.size() < 2 && muons.size() < 2 ) noLeptons=true;


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

     // do nothing (will look for photon)
     //std::cout << "There must be an error this is not possible." << std::endl;
     //exit(9101);

     }


     TLorentzVector diLepton;
     if(leptons.size()==2) diLepton = leptons[0] + leptons[1];

     bool foundZ = false;
     if( leptons.size()==2 )
       if( diLepton.M() > 70. && diLepton.M() < 110. ) foundZ=true;

     // for now, consider only reco Z's (photon done in GammaJet framework)
     if( foundZ ) {

       eZGamma_ = diLepton.Energy();
       ptZGamma_ = diLepton.Pt();
       etaZGamma_ = diLepton.Eta();
       phiZGamma_ = diLepton.Phi();

       nJet_ = nJets_total_ = 0;

       std::vector<AnalysisJet> jets;

       for( unsigned int iJet=0; iJet<nAK5PFJet; ++iJet ) {

         AnalysisJet thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );

         // far away from leptons:
         if( thisJet.DeltaR(leptons[0]) <= 0.5 ) continue;
         if( thisJet.DeltaR(leptons[1]) <= 0.5 ) continue;

         nJets_total_++;

         if( nJet_<10 && thisJet.Pt()>30. ) {
           ptJet_[nJet_] = thisJet.Pt();
           etaJet_[nJet_] = thisJet.Eta();
           phiJet_[nJet_] = thisJet.Phi();
           eJet_[nJet_] = thisJet.Energy();
           nJet_++;
         }

       } //for jets

     } else { // if not found Z


       // look for isolated photon
       TLorentzVector foundPhoton;
       float maxPtPhot=0.;
       for( unsigned iSC=0; iSC<nSC; ++iSC ) {

         TLorentzVector thisPhoton;
         thisPhoton.SetPtEtaPhiE( energySC[iSC]*sin(thetaSC[iSC]), etaSC[iSC], phiSC[iSC], energySC[iSC] );
if( thisPhoton.M()>0. ) std::cout << "mass: " << thisPhoton.M() << std::endl;

         if( thisPhoton.Pt()<50. ) continue;

         // cluster shape:
         if( sMajSC[iSC] > 0.35 ) continue;
         if( sMinSC[iSC] > 0.3 || sMinSC[iSC]<0.15 ) continue;

         // isolation:
         if( scBasedEcalSum04SC[iSC]/thisPhoton.Energy() > 0.05 && scBasedEcalSum04SC[iSC]>3. ) continue;
         if( hcalTowerSumEtConeDR04SC[iSC]/thisPhoton.Energy() > 0.05 && hcalTowerSumEtConeDR04SC[iSC]>2.4 ) continue;
         // compute track isolation with tracks:
         int nTracks035=0;
         int ptTracks035=0.;
         for( unsigned iTrack=0; iTrack<nTrack; ++iTrack ) {
           TLorentzVector thisTrack;
           thisTrack.SetXYZM( pxTrack[iTrack], pyTrack[iTrack], pzTrack[iTrack], 0.140 );
           if( thisTrack.DeltaR(thisPhoton)<0.35 ) {
             nTracks035++;
             ptTracks035+=thisTrack.Pt();
           } // if deltaR
         } // for tracks
         
         if( nTracks035>=3 ) continue;
         if( ptTracks035/thisPhoton.Pt()>0.1 ) continue;

         if( thisPhoton.Pt() > maxPtPhot ) {
           foundPhoton = thisPhoton;
           maxPtPhot = thisPhoton.Pt();
         }

       } // for SC's


       eZGamma_ = (maxPtPhot>0.) ? foundPhoton.Energy() : 0.;
       ptZGamma_ = (maxPtPhot>0.) ? foundPhoton.Pt() : 0.;
       etaZGamma_ = (maxPtPhot>0.) ? foundPhoton.Eta() : 10.;
       phiZGamma_ = (maxPtPhot>0.) ? foundPhoton.Phi() : 0.;

       nJet_ = nJets_total_ = 0;

       std::vector<AnalysisJet> jets;

       for( unsigned int iJet=0; iJet<nAK5PFJet; ++iJet ) {

         AnalysisJet thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );

         // far away from leptons:
         if( thisJet.DeltaR(foundPhoton) <= 0.5 ) continue;

         nJets_total_++;

         if( nJet_<10 && thisJet.Pt()>30. ) {
           ptJet_[nJet_] = thisJet.Pt();
           etaJet_[nJet_] = thisJet.Eta();
           phiJet_[nJet_] = thisJet.Phi();
           eJet_[nJet_] = thisJet.Energy();
           nJet_++;
         }

       } //for jets
     }


     reducedTree_->Fill(); 


   } //for entries

} //loop


double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}


