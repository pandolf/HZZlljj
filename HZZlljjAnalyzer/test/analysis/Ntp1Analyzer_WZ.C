#include "Ntp1Analyzer_WZ.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"
#include "TMVA/Reader.h"



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

};



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




Ntp1Analyzer_WZ::Ntp1Analyzer_WZ( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "WZ", dataset, flags, tree ) {


  //nothing to do here


} //constructor



void Ntp1Analyzer_WZ::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");

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

  reducedTree_->Branch("nJet", &nJet_, "nJet_/I");

  reducedTree_->Branch("iJet",  iJet_,  "iJet_[nJet_]/I");
  reducedTree_->Branch("eJet",  eJet_,  "eJet_[nJet_]/F");
  reducedTree_->Branch( "ptJet",  ptJet_,  "ptJet_[nJet_]/F");
  reducedTree_->Branch("etaJet", etaJet_, "etaJet_[nJet_]/F");
  reducedTree_->Branch("phiJet", phiJet_, "phiJet_[nJet_]/F");

  reducedTree_->Branch("ptDJet", ptDJet_, "ptDJet_[nJet_]/F");
  reducedTree_->Branch("rmsCandJet", rmsCandJet_, "rmsCandJet_[nJet_]/F");
  reducedTree_->Branch("nChargedJet", nChargedJet_, "nChargedJet_[nJet_]/F");
  reducedTree_->Branch("nNeutralJet", nNeutralJet_, "nNeutralJet_[nJet_]/F");
  reducedTree_->Branch("QGlikelihoodJet", QGlikelihoodJet_, "QGlikelihoodJet_[nJet_]/F");

  reducedTree_->Branch("eChargedHadronsJet", eChargedHadronsJet_, "eChargedHadronsJet_[nJet_]/F");
  reducedTree_->Branch("ePhotonsJet", ePhotonsJet_, "ePhotonsJet_[nJet_]/F");
  reducedTree_->Branch("eNeutralEmJet", eNeutralEmJet_, "eNeutralEmJet_[nJet_]/F");
  reducedTree_->Branch("eNeutralHadronsJet", eNeutralHadronsJet_, "eNeutralHadronsJet_[nJet_]/F");
  reducedTree_->Branch("eMuonsJet", eMuonsJet_, "eMuonsJet_[nJet_]/F");
  reducedTree_->Branch("eElectronsJet", eElectronsJet_, "eElectronsJet_[nJet_]/F");
  reducedTree_->Branch("eHFHadronsJet", eHFHadronsJet_, "eHFHadronsJet_[nJet_]/F");
  reducedTree_->Branch("eHFEMJet", eHFEMJet_, "eHFEMJet_[nJet_]/F");

  reducedTree_->Branch("nChargedHadronsJet", nChargedHadronsJet_, "nChargedHadronsJet_[nJet_]/I");
  reducedTree_->Branch("nPhotonsJet", nPhotonsJet_, "nPhotonsJet_[nJet_]/I");
  reducedTree_->Branch("nNeutralHadronsJet", nNeutralHadronsJet_, "nNeutralHadronsJet_[nJet_]/I");
  reducedTree_->Branch("nMuonsJet", nMuonsJet_, "nMuonsJet_[nJet_]/I");
  reducedTree_->Branch("nElectronsJet", nElectronsJet_, "nElectronsJet_[nJet_]/I");
  reducedTree_->Branch("nHFHadronsJet", nHFHadronsJet_, "nHFHadronsJet_[nJet_]/I");
  reducedTree_->Branch("nHFEMJet", nHFEMJet_, "nHFEMJet_[nJet_]/I");


  reducedTree_->Branch("nPart", &nPart_, "nPart_/I");
  reducedTree_->Branch("ePart",  ePart_,  "ePart_[nPart_]/F");
  reducedTree_->Branch( "ptPart",  ptPart_,  "ptPart_[nPart_]/F");
  reducedTree_->Branch("etaPart", etaPart_, "etaPart_[nPart_]/F");
  reducedTree_->Branch("phiPart", phiPart_, "phiPart_[nPart_]/F");
  reducedTree_->Branch("pdgIdPart", pdgIdPart_, "pdgIdPart_[nPart_]/I");


  reducedTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  reducedTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");

  
} 



Ntp1Analyzer_WZ::~Ntp1Analyzer_WZ() {

  outfile_->cd();

}



void Ntp1Analyzer_WZ::Loop()
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

     // good primary vertex requirement:
     if( nPV==0 ) continue;
     bool goodVertex = (ndofPV[0] >= 4.0 && sqrt(PVxPV[0]*PVxPV[0]+PVyPV[0]*PVyPV[0]) < 2. && fabs(PVzPV[0]) < 24. );
     if( !goodVertex ) continue;


     bool isMC = ( runNumber < 5 );


     ptHat_ = (isMC_) ? genPtHat : ptHat_;

   //if( isMC_ ) 
   //  if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;


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
     // MUONS
     // ------------------

     std::vector<AnalysisLepton> muons;
     int chargeFirstMuon;


     for( unsigned int iMuon=0; iMuon<nMuon && (muons.size()<2); ++iMuon ) {

       AnalysisLepton thisMuon( pxMuon[iMuon], pyMuon[iMuon], pzMuon[iMuon], energyMuon[iMuon] );
       thisMuon.charge = chargeMuon[iMuon];

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
         //if( fabs(muons[0].Eta())>2.1 && fabs(thisMuon.Eta())>2.1 ) continue;
         muons.push_back(thisMuon);
       }

     } //for muons



     // ------------------
     // ELECTRONS
     // ------------------

     std::vector<AnalysisLepton> electrons;
     int chargeFirstEle = 0;
     bool firstPassedVBTF80 = false;

     for( unsigned int iEle=0; (iEle<nEle) && (electrons.size()<2); ++iEle ) {

       AnalysisLepton thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );
       thisEle.charge = chargeEle[iEle];


       // --------------
       // kinematics:
       // --------------
       if( thisEle.Pt() < 20. ) continue;
       if( (fabs(thisEle.Eta()) > 2.5) || ( fabs(thisEle.Eta())>1.4442 && fabs(thisEle.Eta())<1.566) ) continue;


       // ELE ID vars:
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

       // CONVERSION REJECTION VARS:
       Int_t nMissingHits_thresh95 = 1;
       Float_t deltaCotTheta_thresh95 = -1.;
       Float_t dist_thresh95 = -1.;

       Int_t nMissingHits_thresh80 = 0;
       Float_t deltaCotTheta_thresh80 = 0.02;
       Float_t dist_thresh80 = 0.02;


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
       
       // ---------------------
       // conversion rejection:
       // ---------------------
       int nMissingHits = expInnerLayersGsfTrack[gsfTrackIndexEle[iEle]];
       bool convRej_VBTF95 = (nMissingHits<=nMissingHits_thresh95) && (fabs(convDistEle[iEle])>dist_thresh95 || fabs(convDcotEle[iEle])>deltaCotTheta_thresh95);
       bool convRej_VBTF80 = (nMissingHits<=nMissingHits_thresh80) && (fabs(convDistEle[iEle])>dist_thresh80 || fabs(convDcotEle[iEle])>deltaCotTheta_thresh80);



       bool passed_VBTF95 = (iso_VBTF95 && eleID_VBTF95 && convRej_VBTF95);
       bool passed_VBTF80 = (iso_VBTF80 && eleID_VBTF80 && convRej_VBTF80);


       if( !passed_VBTF95 ) continue;

       // check that not matched to muon (clean electrons faked by muon MIP):
       bool matchedtomuon=false;
       for( std::vector<AnalysisLepton>::iterator iMu=muons.begin(); iMu!=muons.end(); ++iMu )
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

//   // clean electrons faked by muon MIP in ECAL
//   for( std::vector<TLorentzVector>::iterator iEle=electrons.begin(); iEle!=electrons.end(); ++iEle ) {
//     for( std::vector<TLorentzVector>::iterator iMu=muons.begin(); iMu!=muons.end(); ++iMu ) {
//       if( iMu->DeltaR(*iEle)<0.1 ) {
//         std::cout << "lll" << std::endl;
//         electrons.erase(iEle);
//       }
//     } //for ele
//   } //for mu

//   if( electrons.size() < 2 && muons.size() < 2 ) continue;


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


     } //for i leptons



     // ------------------
     // JETS
     // ------------------

     float jetPt_thresh = 20.;

     // first save leading jets in event:
     std::vector<AnalysisJet> leadJets;
     std::vector<int> leadJetsIndex; //index in the event collection (needed afterwards for PFCandidates)

     //QGLikelihoodCalculator qglc;
     for( unsigned int iJet=0; iJet<nAK5PFJet; ++iJet ) {

       AnalysisJet thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );

       thisJet.eChargedHadrons = chargedHadronEnergyAK5PFJet[iJet];
       thisJet.ePhotons        = photonEnergyAK5PFJet[iJet];
       thisJet.eNeutralEm      = neutralEmEnergyAK5PFJet[iJet];
       thisJet.eNeutralHadrons = neutralHadronEnergyAK5PFJet[iJet];
       thisJet.eElectrons      = electronEnergyAK5PFJet[iJet];
       thisJet.eMuons          = muonEnergyAK5PFJet[iJet];

       thisJet.nChargedHadrons = chargedHadronMultiplicityAK5PFJet[iJet];
       thisJet.nPhotons        = photonMultiplicityAK5PFJet[iJet];
       thisJet.nNeutralHadrons = neutralHadronMultiplicityAK5PFJet[iJet];
       thisJet.nElectrons      = electronMultiplicityAK5PFJet[iJet];
       thisJet.nMuons          = muonMultiplicityAK5PFJet[iJet];

       thisJet.nCharged = chargedHadronMultiplicityAK5PFJet[iJet]+electronMultiplicityAK5PFJet[iJet]+muonMultiplicityAK5PFJet[iJet];
       thisJet.nNeutral = neutralHadronMultiplicityAK5PFJet[iJet]+photonMultiplicityAK5PFJet[iJet];
       thisJet.rmsCand =  rmsCandAK5PFJet[iJet];
       thisJet.ptD =  ptDAK5PFJet[iJet];

       //thisJet.QGlikelihood = qglc.ComputeLikelihood( thisJet.Pt(), thisJet.nCharged, thisJet.nNeutral, thisJet.ptD, thisJet.rmsCand );


       // save at least 3 lead jets (if event has them) and all jets with pt>thresh:
       if( leadJets.size()>=3 && thisJet.Pt()<jetPt_thresh ) break;

       // far away from leptons:
       if( thisJet.DeltaR( leptons[0] ) <= 0.5 ) continue;
       if( thisJet.DeltaR( leptons[1] ) <= 0.5 ) continue;

       // jet ID:
       int multiplicity = thisJet.nCharged +  thisJet.nNeutral + HFEMMultiplicityAK5PFJet[iJet] + HFHadronMultiplicityAK5PFJet[iJet];
       if( multiplicity < 2 ) continue;
       if( fabs(thisJet.Eta())<2.4 && thisJet.nChargedHadrons == 0 ) continue;
       if( thisJet.eNeutralHadrons >= 0.99*thisJet.Energy() ) continue;
       if( thisJet.ePhotons >= 0.99*thisJet.Energy() ) continue;

       leadJets.push_back(thisJet);
       leadJetsIndex.push_back(iJet);

     }




     nJet_ = 0;
     nPart_ = 0;


     for( unsigned iJet=0; iJet<leadJets.size(); ++iJet ) {
   
       AnalysisJet thisJet = leadJets[iJet];

       // --------------
       // kinematics:
       // --------------
       if( thisJet.Pt() < jetPt_thresh ) continue;

       if( nJet_>=10 ) {
      
         std::cout << "MORE than 10 jet pairs found. SKIPPING!!" << std::endl;

       } else {

         eJet_[nJet_] = leadJets[iJet].Energy();
         ptJet_[nJet_] = leadJets[iJet].Pt();
         etaJet_[nJet_] = leadJets[iJet].Eta();
         phiJet_[nJet_] = leadJets[iJet].Phi();
         eChargedHadronsJet_[nJet_] = leadJets[iJet].eChargedHadrons;
         ePhotonsJet_[nJet_]        = leadJets[iJet].ePhotons;
         eNeutralEmJet_[nJet_]      = leadJets[iJet].eNeutralEm;
         eNeutralHadronsJet_[nJet_] = leadJets[iJet].eNeutralHadrons;
         eElectronsJet_[nJet_]      = leadJets[iJet].eElectrons;
         eMuonsJet_[nJet_]          = leadJets[iJet].eMuons;
         nChargedHadronsJet_[nJet_] = leadJets[iJet].nChargedHadrons;
         nPhotonsJet_[nJet_]        = leadJets[iJet].nPhotons;
         nNeutralHadronsJet_[nJet_] = leadJets[iJet].nNeutralHadrons;
         nElectronsJet_[nJet_]      = leadJets[iJet].nElectrons;
         nMuonsJet_[nJet_]          = leadJets[iJet].nMuons;

         ptDJet_[nJet_] = leadJets[iJet].ptD;
         rmsCandJet_[nJet_] = leadJets[iJet].rmsCand;
         nChargedJet_[nJet_] = leadJets[iJet].nCharged;
         nNeutralJet_[nJet_] = leadJets[iJet].nNeutral;
         QGlikelihoodJet_[nJet_] = leadJets[iJet].QGlikelihood;
          
         nJet_++;
        
       }

     } //for i

     



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



