#include "Ntp1Analyzer_TMVA.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"

#include "KinematicFit/DiJetKinFitter.h"
#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"



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


HelicityLikelihoodDiscriminant::HelicityAngles computeHelicityAngles(TLorentzVector leptMinus, TLorentzVector leptPlus, TLorentzVector jet1, TLorentzVector jet2 );


double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);

Double_t ErrEt(Float_t Et, Float_t Eta);
Double_t ErrEta(Float_t Et, Float_t Eta);
Double_t ErrPhi(Float_t Et, Float_t Eta);



Ntp1Analyzer_TMVA::Ntp1Analyzer_TMVA( const std::string& dataset, const std::string& jetChoice, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "TMVA", dataset, flags, tree ) {


  jetChoice_ = jetChoice;
  PRESEL_ = 0;


} //constructor



void Ntp1Analyzer_TMVA::CreateOutputFile() {

  if( PRESEL_!=0 ) {

    char presel[50];
    sprintf( presel, "%dPRESEL", PRESEL_);
    std::string preselstr(presel);

    if( flags_=="" ) flags_ = preselstr;
    else flags_ += ("_" + preselstr);
  }

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");

  reducedTree_->Branch("ptLept1",  &ptLept1_,  "ptLept1_/F");
  reducedTree_->Branch("absEtaLept1",  &absEtaLept1_,  "absEtaLept1_/F");

  reducedTree_->Branch("ptLept2",  &ptLept2_,  "ptLept2_/F");
  reducedTree_->Branch("absEtaLept2",  &absEtaLept2_,  "absEtaLept2_/F");

  reducedTree_->Branch("mZll",  &mZll_,  "mZll_/F");
  reducedTree_->Branch("ptZll",  &ptZll_,  "ptZll_/F");
  reducedTree_->Branch("deltaRll",  &deltaRll_,  "deltaRll_/F");

  reducedTree_->Branch( "ptJet1",  &ptJet1_,  "ptJet1_/F");
  reducedTree_->Branch( "ptJet1_preKin",  &ptJet1_preKin_,  "ptJet1_preKin_/F");
  reducedTree_->Branch("absEtaJet1", &absEtaJet1_, "absEtaJet1_/F");
  reducedTree_->Branch("nChargedJet1", &nChargedJet1_, "nChargedJet1_/I");
  reducedTree_->Branch("nNeutralJet1", &nNeutralJet1_, "nNeutralJet1_/I");
  reducedTree_->Branch("rmsCandJet1", &rmsCandJet1_, "rmsCandJet1_/F");
  reducedTree_->Branch("ptDJet1", &ptDJet1_, "ptDJet1_/F");

  reducedTree_->Branch( "ptJet2",  &ptJet2_,  "ptJet2_/F");
  reducedTree_->Branch( "ptJet2_preKin",  &ptJet2_preKin_,  "ptJet2_preKin_/F");
  reducedTree_->Branch("absEtaJet2", &absEtaJet2_, "absEtaJet2_/F");
  reducedTree_->Branch("nChargedJet2", &nChargedJet2_, "nChargedJet2_/I");
  reducedTree_->Branch("nNeutralJet2", &nNeutralJet2_, "nNeutralJet2_/I");
  reducedTree_->Branch("rmsCandJet2", &rmsCandJet2_, "rmsCandJet2_/F");
  reducedTree_->Branch("ptDJet2", &ptDJet2_, "ptDJet2_/F");

  reducedTree_->Branch( "ptJetRecoil",  &ptJetRecoil_,  "ptJetRecoil_/F");
  reducedTree_->Branch("absEtaJetRecoil", &absEtaJetRecoil_, "absEtaJetRecoil_/F");
  reducedTree_->Branch("deltaR_recoil_jet1", &deltaR_recoil_jet1_, "deltaR_recoil_jet1_/F");
  reducedTree_->Branch("deltaR_recoil_Zjj", &deltaR_recoil_Zjj_, "deltaR_recoil_Zjj_/F");
  reducedTree_->Branch("deltaR_recoil_Higgs", &deltaR_recoil_Higgs_, "deltaR_recoil_Higgs_/F");

  reducedTree_->Branch("QGLikelihoodJet1", &QGLikelihoodJet1_, "QGLikelihoodJet1_/F");
  reducedTree_->Branch("QGLikelihoodJet2", &QGLikelihoodJet2_, "QGLikelihoodJet2_/F");
  reducedTree_->Branch("QGLikelihoodJetRecoil", &QGLikelihoodJetRecoil_, "QGLikelihoodJetRecoil_/F");
  reducedTree_->Branch("QGLikelihoodJet1Jet2", &QGLikelihoodJet1Jet2_, "QGLikelihoodJet1Jet2_/F");
  reducedTree_->Branch("QGLikelihoodJet1Jet2Recoil", &QGLikelihoodJet1Jet2Recoil_, "QGLikelihoodJet1Jet2Recoil_/F");

  reducedTree_->Branch("mZjj",  &mZjj_,  "mZjj_/F");
  reducedTree_->Branch("ptZjj",  &ptZjj_,  "ptZjj_/F");
  reducedTree_->Branch("ptZjj_preKin",  &ptZjj_preKin_,  "ptZjj_preKin_/F");
  reducedTree_->Branch("deltaRjj",  &deltaRjj_,  "deltaRjj_/F");
  reducedTree_->Branch("deltaRjj_preKin",  &deltaRjj_preKin_,  "deltaRjj_preKin_/F");

  reducedTree_->Branch("deltaRZZ",  &deltaRZZ_,  "deltaRZZ_/F");
  reducedTree_->Branch("deltaAbsEtaZZ",  &deltaAbsEtaZZ_,  "deltaAbsEtaZZ_/F");
  reducedTree_->Branch("absDeltaEtaZZ",  &absDeltaEtaZZ_,  "absDeltaEtaZZ_/F");
  reducedTree_->Branch("absDeltaPhiZZ",  &absDeltaPhiZZ_,  "absDeltaPhiZZ_/F");
  reducedTree_->Branch("ptZZ",  &ptZZ_,  "ptZZ_/F");
  reducedTree_->Branch("mZZ",  &mZZ_,  "mZZ_/F");
  reducedTree_->Branch("absEtaZZ",  &absEtaZZ_,  "absEtaZZ_/F");
  reducedTree_->Branch("mZZ_preKin",  &mZZ_preKin_,  "mZZ_preKin_/F");

  reducedTree_->Branch("pfMet",  &epfMet_,  "epfMet_/F");

  reducedTree_->Branch("helicityLD", &helicityLD_, "helicityLD_/F");
  reducedTree_->Branch("helicityLD_kinFit", &helicityLD_kinFit_, "helicityLD_kinFit_/F");

  h1_mZjj = new TH1F("mZjj", "", 50, 0., 300.);
  h1_mZjj_matched = new TH1F("mZjj_matched", "", 50, 0., 300.);

  h1_mZjj_bestZ = new TH1F("mZjj_bestZ", "", 50, 0., 300.);
  h1_mZjj_bestH = new TH1F("mZjj_bestH", "", 50, 0., 300.);
  h1_mZjj_closestPair = new TH1F("mZjj_closestPair", "", 50, 0., 300.);

} 



Ntp1Analyzer_TMVA::~Ntp1Analyzer_TMVA() {

  outfile_->cd();
  h1_mZjj->Write();
  h1_mZjj_matched->Write();
  h1_mZjj_bestZ->Write();
  h1_mZjj_bestH->Write();
  h1_mZjj_closestPair->Write();

}



void Ntp1Analyzer_TMVA::Loop()
{


   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;


   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();


   Long64_t nbytes = 0, nb = 0;


   int nCorrectPairs_bestZ = 0;
   int nIncorrectPairs_bestZ = 0;

   int nCorrectPairs_bestH = 0;
   int nIncorrectPairs_bestH = 0;

   int nCorrectPairs_closestPair = 0;
   int nIncorrectPairs_closestPair = 0;

   QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator();


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%100000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;

     eventWeight_ = -1.; //default

     if( !isGoodEvent() ) continue; //this takes care also of integrated luminosity and trigger

     //trigger:
     // not yet

     TRegexp re("Higgs");
     TString dataset_str(dataset_);
     bool isSignal = dataset_str.Contains(re);


     ptHat_ = genPtHat;

     if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;


     std::vector<TLorentzVector> quarksMC;


     if( isSignal ) { //only signal has Z->qq

       // first look for Z->qq

       for( unsigned iMc=0; iMc<nMc && quarksMC.size()<2; ++iMc ) {

         // quarks have status 3
         if( statusMc[iMc] != 3 ) continue;

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

         if( fabs(idMc[iMc])<7 && idMc[mothMc[iMc]]==23 ) {
           quarksMC.push_back( *thisParticle );
         }

       }

       // (checked that always 2 quarks are found)
       if( quarksMC.size()!=2 ) {

         std::cout << "Found a number of quarks different from 2 (" << quarksMC.size() << ")!! PROBLEM!!" << std::endl;
         continue; // skip event

       }

     }



     // now look for Z->ll

     bool noLeptons = false;
     TLorentzVector lept1MC, lept2MC;
     std::vector<TLorentzVector> electronsMC;
     std::vector<TLorentzVector> muonsMC;

     for( unsigned iMc=0; iMc<nMc; ++iMc ) {

       if( statusMc[iMc] != 1 ) continue;

       TLorentzVector* thisParticle = new TLorentzVector();
       thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

       // remember: a stable lepton is daughter of a parton lepton, which is daughter of the Z:
       if( idMc[mothMc[mothMc[iMc]]]==23 ) {
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


     if( !noLeptons )
       if( lept1MC.Pt() < lept2MC.Pt() ) std::cout << "WARNING MC leptons not ordered in pt!!" << std::endl;



     // -----------------------------
     //      FROM NOW ON RECO
     // -----------------------------


     epfMet_ = energyPFMet[0];
     //phipfMet_ = phiPFMet[0];


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




     ptLept1_ = leptons[0].Pt();
     absEtaLept1_ = fabs(leptons[0].Eta());
     int chargeLept1 = leptons[0].charge;
     
     ptLept2_ = leptons[1].Pt();
     absEtaLept2_ = fabs(leptons[1].Eta());
     int chargeLept2 = leptons[1].charge;

     TLorentzVector diLepton = leptons[0] + leptons[1];

     mZll_ = diLepton.M();
     ptZll_ = diLepton.Pt();
     deltaRll_ = leptons[0].DeltaR(leptons[1]);



     // ------------------
     // JETS
     // ------------------

     float jetPt_thresh = 30.;


     // first save leading jets in event:
     std::vector<AnalysisJet> leadJets;
     std::vector<int> leadJetsIndex; //index in the event collection (needed afterwards for PFCandidates)
     int nJets30=0;

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

       if( thisJet.Pt()>jetPt_thresh ) nJets30++;

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


     if( leadJets.size()<2 ) continue;
     if( leadJets[1].Pt()<jetPt_thresh ) continue; //at least 2 jets over thresh


     int matchedIndex_quark0 = -1;     
     int matchedIndex_quark1 = -1;     
     if( isSignal ) {
 
       // look for jets matched to partons:
       float deltaRmax0 = 999.;
       for( unsigned i=0; i<leadJets.size(); ++i ) {
         float deltaR = quarksMC[0].DeltaR(leadJets[i]);
         if( deltaR<deltaRmax0 ) {
           deltaRmax0 = deltaR;
           matchedIndex_quark0 = i;
         }
       }

       float deltaRmax1 = 999.;
       for( unsigned i=0; i<leadJets.size(); ++i ) {
         float deltaR = quarksMC[1].DeltaR(leadJets[i]);
         if( deltaR<deltaRmax1 && i!=matchedIndex_quark0 ) {
           deltaRmax1 = deltaR;
           matchedIndex_quark1 = i;
         }
       }

     } //is signal



     // now look for best  jet pair 
     float Zmass = 91.19;
     float bestMass_Z = 0.;
     float bestMass_H = 0.;
     float deltaRMin_jetPair = 999.;
     int best_i_bestZ=-1;
     int best_j_bestZ=-1;
     int recoil_i_bestZ=-1;
     int best_i_bestH=-1;
     int best_j_bestH=-1;
     int recoil_i_bestH=-1;
     int best_i_closestPair=-1;
     int best_j_closestPair=-1;
     int recoil_i_closestPair=-1;

     int nPairs_ = 0;


     for( unsigned iJet=0; iJet<leadJets.size(); ++iJet ) {
   
       TLorentzVector thisJet = leadJets[iJet];

       // --------------
       // kinematics:
       // --------------
       if( thisJet.Pt() < jetPt_thresh ) continue;
       if( fabs(thisJet.Eta()) > 2.4 ) continue;



       for( unsigned int jJet=iJet+1; jJet<leadJets.size(); ++jJet ) {

         TLorentzVector otherJet = leadJets[jJet];

         // --------------
         // kinematics:
         // --------------
         if( otherJet.Pt() < jetPt_thresh ) continue;
         if( fabs(otherJet.Eta()) > 2.4 ) continue;


         if( nPairs_>=50 ) {
        
           std::cout << "MORE than 50 jet pairs found. SKIPPING!!" << std::endl;

         } else {

   //      ptJet1_[nPairs_] = leadJets[iJet].Pt();
   //      absEtaJet1_[nPairs_] = fabs(leadJets[iJet].Eta());
   //       
   //      ptJet2_[nPairs_] = leadJets[jJet].Pt();
   //      etaJet2_[nPairs_] = leadJets[jJet].Eta();
   //      phiJet2_[nPairs_] = leadJets[jJet].Phi();

           nPairs_++;
          
         }


         // once jet 1 and 2 are defined, define recoil jet as hardest jet in event excluding 1 and 2
         int i_recoil = -1;
         for( unsigned ii =0; ii<leadJets.size(); ++ii ) {
           if( ii!=iJet && ii!=jJet ) {
             i_recoil = ii;
             break;
           }
         }
         

         TLorentzVector dijet = thisJet + otherJet;
         float invMass_Z = dijet.M();
         
         h1_mZjj->Fill( invMass_Z );
         if( (iJet==matchedIndex_quark0 && jJet==matchedIndex_quark1) || (iJet==matchedIndex_quark1 && jJet==matchedIndex_quark0) ) {
           h1_mZjj_matched->Fill( invMass_Z );
         }

         if( (best_i_bestZ==-1 && best_j_bestZ==-1 ) || ( fabs(invMass_Z-Zmass) < fabs(bestMass_Z-Zmass) ) ) {
           bestMass_Z = invMass_Z;
           best_i_bestZ = iJet;
           best_j_bestZ = jJet;
           recoil_i_bestZ = i_recoil;
         }

         TLorentzVector higgs = dijet + leptons[0] + leptons[1];
         float invMass_H = higgs.M();

         if( (best_i_bestH==-1 && best_j_bestH==-1 ) || ( fabs(invMass_H-400.) < fabs(bestMass_H-400.) ) ) {
           bestMass_H = invMass_H;
           best_i_bestH = iJet;
           best_j_bestH = jJet;
           recoil_i_bestH = i_recoil;
         }

         float thisDeltaR = thisJet.DeltaR(otherJet);
         if( (best_i_closestPair==-1 && best_j_closestPair==-1 ) || ( thisDeltaR<deltaRMin_jetPair ) ) {
           deltaRMin_jetPair = thisDeltaR;
           best_i_closestPair = iJet;
           best_j_closestPair = jJet;
           recoil_i_closestPair = i_recoil;
         }

       } //for j
     } //for i

     
     //if( jets.size() < 2 ) continue;
     if( best_i_bestZ==-1 || best_j_bestZ==-1 ) {
       //std::cout << "leadjets size: " << leadJets.size() << " best_i: " << best_i_bestZ << " best_j: " << best_j_bestZ << std::endl;
       continue; //means that less than 2 jets were found
     }

     AnalysisJet jet1, jet2, jetRecoil;

     if( jetChoice_ == "BESTZ" ) {
       jet1 = leadJets[best_i_bestZ];
       jet2 = leadJets[best_j_bestZ];
       if( recoil_i_bestZ >=0 ) jetRecoil = leadJets[recoil_i_bestZ];
     } else if( jetChoice_ == "BESTH" ) {
       jet1 = leadJets[best_i_bestH];
       jet2 = leadJets[best_j_bestH];
       if( recoil_i_bestH >=0 ) jetRecoil = leadJets[recoil_i_bestH];
     } else if( jetChoice_ == "CLOSESTPAIR" ) {
       jet1 = leadJets[best_i_closestPair];
       jet2 = leadJets[best_j_closestPair];
       if( recoil_i_closestPair >=0 ) jetRecoil = leadJets[recoil_i_closestPair];
     } else {
       std::cout << "Unknown jet choice option '" << jetChoice_ << "'. Exiting." << std::endl;
       exit(919);
     }
     

     if( (best_i_bestZ==matchedIndex_quark0 && best_j_bestZ==matchedIndex_quark1) || (best_i_bestZ==matchedIndex_quark1 && best_j_bestZ==matchedIndex_quark0) )
       nCorrectPairs_bestZ++;
     else
       nIncorrectPairs_bestZ++;

     if( (best_i_bestH==matchedIndex_quark0 && best_j_bestH==matchedIndex_quark1) || (best_i_bestH==matchedIndex_quark1 && best_j_bestH==matchedIndex_quark0) )
       nCorrectPairs_bestH++;
     else
       nIncorrectPairs_bestH++;

     if( (best_i_closestPair==matchedIndex_quark0 && best_j_closestPair==matchedIndex_quark1) || (best_i_closestPair==matchedIndex_quark1 && best_j_closestPair==matchedIndex_quark0) )
       nCorrectPairs_closestPair++;
     else
       nIncorrectPairs_closestPair++;


     // compute QG likelihood before kinfit:
     QGLikelihoodJet1_ = qglikeli->computeQGLikelihood( jet1.Pt(), jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );
     QGLikelihoodJet2_ = qglikeli->computeQGLikelihood( jet2.Pt(), jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );
   //float QGLikeli_recoil = -1.;
   //if( jetRecoil.Energy()>0. jetRecoil.ptD>=0. && jetRecoil.nCharged>) {
   //  if( fabs(jetRecoil.Eta())<2. ) {
   //    QGLikeli_recoil = qglikeli->computeQGLikelihood( jetRecoil.Pt(), jetRecoil.nCharged, jetRecoil.nNeutral, jetRecoil.ptD, -1. );
     QGLikelihoodJetRecoil_ = (jetRecoil.Pt()>30.) ? qglikeli->computeQGLikelihood( jetRecoil.Pt(), jetRecoil.nCharged, jetRecoil.nNeutral, jetRecoil.ptD, -1. ) : -1.;

     QGLikelihoodJet1Jet2_ = QGLikelihoodJet1_*QGLikelihoodJet2_;
     QGLikelihoodJet1Jet2Recoil_ = (QGLikelihoodJetRecoil_>=0.) ? QGLikelihoodJet1_*QGLikelihoodJet2_*QGLikelihoodJetRecoil_ : QGLikelihoodJet1_*QGLikelihoodJet2_;

     // ------------------------
     //   KINEMATIC FIT: BEGIN
     // ------------------------
     
     TMatrixD m_jet1(3,3);
     TMatrixD m_jet2(3,3);
     
     m_jet1(0,0) = 0.5*ErrEt (jet1.Et(), jet1.Eta()); // et
     m_jet1(1,1) = 0.5*ErrEta(jet1.Et(), jet1.Eta()); // eta
     m_jet1(2,2) = 0.5*ErrPhi(jet1.Et(), jet1.Eta()); // phi
     m_jet2(0,0) = 0.5*ErrEt (jet2.Et(), jet2.Eta()); // et
     m_jet2(1,1) = 0.5*ErrEta(jet2.Et(), jet2.Eta()); // eta
     m_jet2(2,2) = 0.5*ErrPhi(jet2.Et(), jet2.Eta()); // phi
     
     TFitParticleEtEtaPhi *fitJet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &jet1, &m_jet1 );
     TFitParticleEtEtaPhi *fitJet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &jet2, &m_jet2 );
     
     TFitConstraintM *mCons_jets = new TFitConstraintM( "ZMassConstraint_jets", "ZMass-Constraint", 0, 0 , 91.19);
     mCons_jets->addParticles1( fitJet1, fitJet2 );
     
     TKinFitter* fitter_jets = new TKinFitter("fitter_jets", "fitter_jets");
     fitter_jets->addMeasParticle( fitJet1 );
     fitter_jets->addMeasParticle( fitJet2 );
     fitter_jets->addConstraint( mCons_jets );
     
     //Set convergence criteria
     fitter_jets->setMaxNbIter( 30 );
     fitter_jets->setMaxDeltaS( 1e-2 );
     fitter_jets->setMaxF( 1e-1 );
     fitter_jets->setVerbosity(0);
     
     //Perform the fit
     fitter_jets->fit();
     

     TLorentzVector jet1_kinfit(*fitJet1->getCurr4Vec());
     TLorentzVector jet2_kinfit(*fitJet2->getCurr4Vec());

     // fill mZjj before the kinfit:
     TLorentzVector diJet_pre = jet1 + jet2;
     TLorentzVector ZZ_pre = diJet_pre + diLepton;
     mZjj_ = diJet_pre.M();
     ptJet1_preKin_ = jet1.Pt();
     ptJet2_preKin_ = jet2.Pt();
     ptZjj_preKin_ = diJet_pre.Pt();
     deltaRjj_preKin_ = jet1.DeltaR(jet2);
     mZZ_preKin_ = ZZ_pre.M();


     //helicity angles: (also before the kinfit)
     HelicityLikelihoodDiscriminant::HelicityAngles hangles;
     if( chargeLept1<0 ) hangles = computeHelicityAngles(leptons[0], leptons[1], jet1, jet2);
     else                hangles = computeHelicityAngles(leptons[1], leptons[0], jet1, jet2);

     HelicityLikelihoodDiscriminant *LD = new HelicityLikelihoodDiscriminant();
     LD->setMeasurables(hangles);
     double sProb = LD->getSignalProbability();
     double bProb = LD->getBkgdProbability();
     helicityLD_ = sProb/(sProb+bProb);


if( jet2.Energy()<0. ) std::cout << "BLABLA: " << jentry << std::endl;

     jet1.SetPtEtaPhiE( jet1_kinfit.Pt(), jet1_kinfit.Eta(), jet1_kinfit.Phi(), jet1_kinfit.E() );
     jet2.SetPtEtaPhiE( jet2_kinfit.Pt(), jet2_kinfit.Eta(), jet2_kinfit.Phi(), jet2_kinfit.E() );

if( jet2.Energy()<0. ) std::cout << "jet2_kinfit.E(): " << jet2_kinfit.E() << std::endl;
     //helicity angles after kinfit:
     HelicityLikelihoodDiscriminant::HelicityAngles hangles_kinFit;
     if( chargeLept1<0 ) hangles_kinFit = computeHelicityAngles(leptons[0], leptons[1], jet1, jet2);
     else                hangles_kinFit = computeHelicityAngles(leptons[1], leptons[0], jet1, jet2);

     LD->setMeasurables(hangles_kinFit);
     double sProb_kinFit = LD->getSignalProbability();
     double bProb_kinFit = LD->getBkgdProbability();
     if( TMath::IsNaN(sProb_kinFit) || TMath::IsNaN(bProb_kinFit) ) helicityLD_kinFit_ = 0.;
     else  helicityLD_kinFit_ = sProb_kinFit/(sProb_kinFit+bProb_kinFit);




     TLorentzVector diJet_post = jet1 + jet2;
     TLorentzVector ZZ = diJet_post + diLepton;


     TLorentzVector diJet_bestZ = leadJets[best_i_bestZ] + leadJets[best_j_bestZ];
     TLorentzVector diJet_bestH = leadJets[best_i_bestH] + leadJets[best_j_bestH];
     TLorentzVector diJet_closestPair = leadJets[best_i_closestPair] + leadJets[best_j_closestPair];

     h1_mZjj_bestZ->Fill( diJet_bestZ.M() );
     h1_mZjj_bestH->Fill( diJet_bestH.M() );
     h1_mZjj_closestPair->Fill( diJet_closestPair.M() );

     ptJet1_ = jet1.Pt();
     absEtaJet1_ = fabs(jet1.Eta());
     nChargedJet1_ = jet1.nCharged;
     nNeutralJet1_ = jet1.nNeutral;
     ptDJet1_ = jet1.ptD;
     rmsCandJet1_ = jet1.rmsCand;
     
     ptJet2_ = jet2.Pt();
     absEtaJet2_ = fabs(jet2.Eta());
     nChargedJet2_ = jet2.nCharged;
     nNeutralJet2_ = jet2.nNeutral;
     ptDJet2_ = jet2.ptD;
     rmsCandJet2_ = jet2.rmsCand;


     ptJetRecoil_ = (jetRecoil.E()>0.) ? jetRecoil.Pt() : 0.;
     absEtaJetRecoil_ = (jetRecoil.E()>0.) ? fabs(jetRecoil.Eta()) : -1.;
     deltaR_recoil_jet1_ = (jetRecoil.E()>0.) ? jetRecoil.DeltaR(jet1) : -1.;


     deltaRjj_ = jet1.DeltaR(jet2);
     ptZjj_ = diJet_post.Pt();

     deltaR_recoil_Zjj_ = (jetRecoil.E()>0.) ? jetRecoil.DeltaR(diJet_post) : -1.;
     deltaR_recoil_Higgs_ = (jetRecoil.E()>0.) ? jetRecoil.DeltaR(ZZ) : -1.;

     deltaRZZ_ = diLepton.DeltaR(diJet_post);
     deltaAbsEtaZZ_ = fabs(diLepton.Eta()) - fabs(diJet_post.Eta());
     absDeltaEtaZZ_ = fabs(diLepton.Eta() - diJet_post.Eta());
     absDeltaPhiZZ_ = fabs(diJet_post.DeltaPhi(diLepton));
     ptZZ_ = ZZ.Pt();
     mZZ_ = ZZ.M();
     absEtaZZ_ = fabs(ZZ.Eta());

/*
if( !(sProb>0. && sProb<1.) ) {
  std::cout << std::endl << "jentry: " << jentry << std::endl;
  std::cout << "sProb: " << sProb << std::endl;
  std::cout << "bProb: " << bProb << std::endl;
  std::cout << "helicityLD_: " << helicityLD_ << std::endl;
  std::cout << "hangles.helPhi: " << hangles.helPhi << std::endl;
  std::cout << "hangles.helPhi1: " << hangles.helPhi1 << std::endl;
  std::cout << "hangles.helCosThetaStar: " << hangles.helCosThetaStar << std::endl;
  std::cout << "hangles.helCosTheta1: " << hangles.helCosTheta1 << std::endl;
  std::cout << "hangles.helCosTheta2: " << hangles.helCosTheta2 << std::endl;
  std::cout << "hangles.mzz: " << hangles.mzz << std::endl;
exit(1);
}*/
     
     bool eventOK = true;

  // if( PRESEL_!=0 ) {
  //   if( ptLept1_ < 50. ) eventOK = false;
  //   if( absEtaLept1_ > 2.1 ) eventOK = false;
  //   if( deltaRll_ > 2. ) eventOK = false;
  //   if( ptJet1_ < 70. ) eventOK = false;
  //   if( mZjj_ < 60. ) eventOK = false;
  //   if( mZjj_ > 200. ) eventOK = false;
  //   if( deltaRjj_ > 2. ) eventOK = false;
  //   if( PRESEL_==400 ) {
  //     if( mZZ_ < 320. ) eventOK = false;
  //     if( mZZ_ > 480. ) eventOK = false;
  //   } else if( PRESEL_==500 ) {
  //     if( mZZ_ < 420. ) eventOK = false;
  //     if( mZZ_ > 580. ) eventOK = false;
  //   }
  // }
     
     if( eventOK )
       reducedTree_->Fill(); 


   } //for entries


   std::cout << "BEST Z: " << std::endl;
   std::cout << "Correct Pairs: " << nCorrectPairs_bestZ << " (" << (float)nCorrectPairs_bestZ/((float)nCorrectPairs_bestZ+(float)nIncorrectPairs_bestZ)*100. << "%)" << std::endl;
   std::cout << "Incorrect Pairs: " << nIncorrectPairs_bestZ << " (" << (float)nIncorrectPairs_bestZ/((float)nCorrectPairs_bestZ+(float)nIncorrectPairs_bestZ)*100. << "%)" << std::endl;

   std::cout << std::endl << "BEST H: " << std::endl;
   std::cout << "Correct Pairs: " << nCorrectPairs_bestH << " (" << (float)nCorrectPairs_bestH/((float)nCorrectPairs_bestH+(float)nIncorrectPairs_bestH)*100. << "%)" << std::endl;
   std::cout << "Incorrect Pairs: " << nIncorrectPairs_bestH << " (" << (float)nIncorrectPairs_bestH/((float)nCorrectPairs_bestH+(float)nIncorrectPairs_bestH)*100. << "%)" << std::endl;

   std::cout << std::endl << "CLOSEST PAIR:" << std::endl;
   std::cout << "Correct Pairs: " << nCorrectPairs_closestPair << " (" << (float)nCorrectPairs_closestPair/((float)nCorrectPairs_closestPair+(float)nIncorrectPairs_closestPair)*100. << "%)" << std::endl;
   std::cout << "Incorrect Pairs: " << nIncorrectPairs_closestPair << " (" << (float)nIncorrectPairs_closestPair/((float)nCorrectPairs_closestPair+(float)nIncorrectPairs_closestPair)*100. << "%)" << std::endl;

} //loop


double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}

// error functions for jets:

Double_t ErrEt(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 5.6;
    b = 1.25;
    c = 0.033;
  }
  else{
    a = 4.8;
    b = 0.89;
    c = 0.043;
  }
  InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
  return InvPerr2;
}



Double_t ErrEta(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 1.215;
    b = 0.037;
    c = 7.941 * 0.0001;
  }
  else{
    a = 1.773;
    b = 0.034;
    c = 3.56 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

Double_t ErrPhi(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 6.65;
    b = 0.04;
    c = 8.49 * 0.00001;
  }
  else{
    a = 2.908;
    b = 0.021;
    c = 2.59 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}
               

HelicityLikelihoodDiscriminant::HelicityAngles computeHelicityAngles(TLorentzVector leptMinus, TLorentzVector leptPlus, TLorentzVector jet1, TLorentzVector jet2 ) {

  HelicityLikelihoodDiscriminant::HelicityAngles returnAngles;

  TLorentzVector Zll = leptPlus + leptMinus;
  TLorentzVector Zjj = jet1 + jet2;

  TLorentzVector Higgs = Zjj + Zll;

  
  // define lept1 as negatively charged lepton:
  TLorentzVector lept1 = leptMinus;
  TLorentzVector lept2 = leptPlus;

  // no charge for jets (well, not really)
  // so choose jet with positive scalar product with Zjj 
  // in its restframe:
  TLorentzVector jet1_Zjjstar_tmp(jet1);
  jet1_Zjjstar_tmp.Boost(-Zjj.BoostVector());
  if( jet1_Zjjstar_tmp.Phi()<0. ) { //swap them
    TLorentzVector jet1_tmp(jet1);
    TLorentzVector jet2_tmp(jet2);
    jet1 = jet2_tmp;
    jet2 = jet1_tmp;
  }


  //     BOOSTS:

  // boosts in Higgs CoM frame:
  TLorentzVector lept1_Hstar(lept1);
  lept1_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector lept2_Hstar(lept2);
  lept2_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector jet1_Hstar(jet1);
  jet1_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector jet2_Hstar(jet2);
  jet2_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector Zll_Hstar(Zll);
  Zll_Hstar.Boost(-Higgs.BoostVector());
  TLorentzVector Zjj_Hstar(Zjj);
  Zjj_Hstar.Boost(-Higgs.BoostVector());

  // boosts in Zll CoM frame:
  TLorentzVector lept1_Zllstar(lept1);
  lept1_Zllstar.Boost(-Zll.BoostVector());
  TLorentzVector H_Zllstar(Higgs);
  H_Zllstar.Boost(-Zll.BoostVector());

  // boosts in Zjj CoM frame:
  TLorentzVector jet1_Zjjstar(jet1);
  jet1_Zjjstar.Boost(-Zjj.BoostVector());
  TLorentzVector H_Zjjstar(Higgs);
  H_Zjjstar.Boost(-Zjj.BoostVector());


  returnAngles.helCosThetaStar = Zll_Hstar.CosTheta();

if( jet1.Energy() < 0. ) std::cout << "jet1 energy < 0!!!!!!" << std::endl;
if( jet2.Energy() < 0. ) std::cout << "jet2 energy < 0!!!!!!" << std::endl;
if( jet1_Zjjstar.Energy() < 0. ) std::cout << "jet1_Zjjstar energy < 0!!!!!!" << std::endl;

  TVector3 v_pbeamLAB( 0.0, 0.0, 1.0 );

  //cross prod beam x Zll
  TVector3 v_1 = (v_pbeamLAB.Cross(  (Zll_Hstar.Vect()).Unit()) ).Unit();//versor normal to z-z' plane


  //v_2 = cross prod l1 x l2 = versor normal to Zll decay plane
  // careful to the order: L1, the z-axis and Z->ll make a right-handed (non-orthogonal) frame (xyz); at the end we want the angle btw x and y
  TVector3 v_2((Zll_Hstar.Vect().Cross(lept1_Hstar.Vect().Unit())).Unit());


  //v_3 = similar to v_2, BUT
  //now, if we want a right-handed set of Unit-vectors, keeping the same direction of the z-axis
  //we must swap the direction of one of the other two vectors of the Z bosons. 
  //Keeping the same direction of the z-axis
  //means measuring phiZll and phiZjj w.r.t. to the same v_1 vector (i.e. w.r.t. the same z'-Zll plane)
  TVector3 v_3(((-1.0*Zjj_Hstar.Vect()).Cross(jet1_Hstar.Vect().Unit())).Unit()) ;

  //in other terms: we can define v_3 as above and then do the crss prod with v_1
  //or define v_3 in a way consistent with v_2 and then do the cross product with a newly defined
  //Unit vector v_4 =  (v_pbeamLAB.Cross(  (ZjjboostedX->momentum()).Unit()) ).Unit();//versor normal to z-Zjj plane
 
  // helphiZll:
  float phiZll = fabs( acos(v_1.Dot(v_2)) );
  if(v_pbeamLAB.Dot(v_2)>0.0)phiZll=-1.0*phiZll;
  else phiZll=+1.0*phiZll;

  // helphiZjj:
  float phiZjj = fabs( acos(v_1.Dot(v_3)) );
  if(v_pbeamLAB.Dot(v_3)>0.0)phiZjj=+1.0*phiZjj; 
  else phiZjj=-1.0*phiZjj;


  float phi1 = phiZll;


  //phi
  float phi = fabs( acos(v_2.Dot(v_3)) );//two-fold ambiguity when doing the acos + pi ambiguity from sign of v_3 
  if(lept1_Hstar.Vect().Dot(v_3)>0.0)phi= +1.0 * phi;
  else phi= -1.0 * phi;

  returnAngles.helPhi1 = phi1;
  returnAngles.helPhi = phi;


  returnAngles.helCosTheta1 =  (-1.0*(lept1_Zllstar.X()* H_Zllstar.X()+
                                   lept1_Zllstar.Y()* H_Zllstar.Y()+
                                   lept1_Zllstar.Z()* H_Zllstar.Z())/
                                  (lept1_Zllstar.Vect().Mag()* H_Zllstar.Vect().Mag())  );


  returnAngles.helCosTheta2 =  fabs( (jet1_Zjjstar.X()* H_Zjjstar.X()+
                                   jet1_Zjjstar.Y()* H_Zjjstar.Y()+
                                   jet1_Zjjstar.Z()* H_Zjjstar.Z())/
                                  (jet1_Zjjstar.Vect().Mag()* H_Zjjstar.Vect().Mag())  );

  returnAngles.mzz = Higgs.M();


  return returnAngles;

}
