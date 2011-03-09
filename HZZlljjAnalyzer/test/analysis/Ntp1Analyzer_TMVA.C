#include "Ntp1Analyzer_TMVA.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"

#include "TFitConstraintM.h"
#include "TFitParticleEtEtaPhi.h"
#include "TKinFitter.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"



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



double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);
void calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22,
double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2, bool &swappedZ1);

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
  reducedTree_->Branch("absEtaJet1", &absEtaJet1_, "absEtaJet1_/F");

  reducedTree_->Branch( "ptJet2",  &ptJet2_,  "ptJet2_/F");
  reducedTree_->Branch("absEtaJet2", &absEtaJet2_, "absEtaJet2_/F");

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
  reducedTree_->Branch("mZZ",  &mZZ_,  "mZZ_/F");

  reducedTree_->Branch("pfMet",  &epfMet_,  "epfMet_/F");

  reducedTree_->Branch("costheta1", &costheta1_, "costheta1_/F");
  reducedTree_->Branch("costheta2", &costheta2_, "costheta2_/F");
  reducedTree_->Branch("phi", &phi_, "phi_/F");
  reducedTree_->Branch("costhetastar", &costhetastar_, "costhetastar_/F");
  reducedTree_->Branch("phistar1", &phistar1_, "phistar1_/F");
  reducedTree_->Branch("phistar2", &phistar2_, "phistar2_/F");
  reducedTree_->Branch("phistar12", &phistar12_, "phistar12_/F");
  reducedTree_->Branch("phi1", &phi1_, "phi1_/F");
  reducedTree_->Branch("phi2", &phi2_, "phi2_/F");
  reducedTree_->Branch("swappedZ1", &swappedZ1_, "swappedZ1_/O");

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
     
     ptLept2_ = leptons[1].Pt();
     absEtaLept2_ = fabs(leptons[1].Eta());

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
     mZjj_ = diJet_pre.M();
     ptZjj_preKin_ = diJet_pre.Pt();
     deltaRjj_preKin_ = jet1.DeltaR(jet2);

     jet1.SetPtEtaPhiE( jet1_kinfit.Pt(), jet1_kinfit.Eta(), jet1_kinfit.Phi(), jet1_kinfit.E() );
     jet2.SetPtEtaPhiE( jet2_kinfit.Pt(), jet2_kinfit.Eta(), jet2_kinfit.Phi(), jet2_kinfit.E() );

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
     
     ptJet2_ = jet2.Pt();
     absEtaJet2_ = fabs(jet2.Eta());


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
               


/*
void calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22,
double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2, bool &swappedZ1){
        //std::cout << "In calculate angles..." << std::endl;

  double norm;
  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame( thep4Z1 );
  TLorentzVector thep4Z2inXFrame( thep4Z2 );       thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );
  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
     // calculate phi1, phi2, costhetastar
  phi1 = theZ1X_p3.Phi();
  phi2 = theZ2X_p3.Phi();
     ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4H = thep4H;
  if ((phi1 < 0)&&(phi1 >= -TMath::Pi())){
    std::cout<<">>> Swapping the Z in calculateAngles (z1/z2 convention)  (phi1="<<phi1<<", phi2="<<phi2<<")"<<std::endl;
      p4Z1 = thep4Z2; p4M11 = thep4M21; p4M12 = thep4M22;
      p4Z2 = thep4Z1; p4M21 = thep4M11; p4M22 = thep4M12;               
      costhetastar = theZ2X_p3.CosTheta();
      swappedZ1=true;
  }
  else{
      p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
      p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
      costhetastar = theZ1X_p3.CosTheta();
      swappedZ1=false;
  }

        //std::cout << "phi1: " << phi1 << ", phi2: " << phi2 << std::endl;
     // now helicity angles................................
  // ...................................................
  TVector3 boostZ1 = -(p4Z1.BoostVector());
  TLorentzVector p4Z2Z1(p4Z2);
  p4Z2Z1.Boost(boostZ1);
  //find the decay axis
  /////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
  TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
  norm = 1/(unitx_1.Mag());
  unitx_1*=norm;
  //boost daughters of z2
  TLorentzVector p4M21Z1(p4M21);
  TLorentzVector p4M22Z1(p4M22);
  p4M21Z1.Boost(boostZ1);
  p4M22Z1.Boost(boostZ1);
  //create z and y axes
  /////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
  TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
  TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
  TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
  norm = 1/(unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1);
  //caculate theta1
  TLorentzVector p4M11Z1(p4M11);
  p4M11Z1.Boost(boostZ1);
  TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
  TVector3 unitM11 = p3M11.Unit();
  double x_m11 = unitM11.Dot(unitx_1); 
  double y_m11 = unitM11.Dot(unity_1); 
  double z_m11 = unitM11.Dot(unitz_1);
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
  costheta1 = M11_Z1frame.CosTheta();
  //std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
  //////-----------------------old way of calculating phi---------------/////////
  phi = M11_Z1frame.Phi();
  //set axes for other system
  TVector3 boostZ2 = -(p4Z2.BoostVector());
  TLorentzVector p4Z1Z2(p4Z1);
  p4Z1Z2.Boost(boostZ2);
  TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
  norm = 1/(unitx_2.Mag());
  unitx_2*=norm;
  //boost daughters of z2
  TLorentzVector p4M11Z2(p4M11);
  TLorentzVector p4M12Z2(p4M12);
  p4M11Z2.Boost(boostZ2);
  p4M12Z2.Boost(boostZ2);
  TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
  TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
  TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
  norm = 1/(unitz_2.Mag());
  unitz_2*=norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2);
  //calcuate theta2
  TLorentzVector p4M21Z2(p4M21);
  p4M21Z2.Boost(boostZ2);
  TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
  TVector3 unitM21 = p3M21.Unit();
  double x_m21 = unitM21.Dot(unitx_2); 
  double y_m21 = unitM21.Dot(unity_2); 
  double z_m21 = unitM21.Dot(unitz_2);
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
  costheta2 = M21_Z2frame.CosTheta();
     // calculate phi
  //calculating phi_n
  TLorentzVector n_p4Z1inXFrame( p4Z1 );
  TLorentzVector n_p4M11inXFrame( p4M11 );
  n_p4Z1inXFrame.Boost( boostX );
  n_p4M11inXFrame.Boost( boostX );           
  TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
  TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();     
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  TVector3 n_unitz_1( n_p4Z1inXFrame_unit );//z' -> Z1 fly direction 
  TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross( n_unitz_1 );//vector orthog to neg lept from Z1 and z' (i.e., Z1)-> vect normal to Z1 decay plane
  TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );//in the Z1 decay plane, points opposite to neg lepton  
     ///////-----------------new way of calculating phi-----------------///////
  TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
  TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );
  //transform the beam axis to the coordinate system defined by unit_x, unit_y e unit_z as above
  //take the phi of it:phistar1
  phistar1 = n_p4PartoninXFrame_unitprime.Phi();
  // and the calculate phistar2
  TLorentzVector n_p4M21inXFrame( p4M21 );
  n_p4M21inXFrame.Boost( boostX );
  TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();//3-vector normalized parallel to pos charge lepton from Z1
  //   rotate into other plane
  TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );
  std::cout<<"n_p4M21inXFrame_unitprime = ("<< n_p4M21inXFrame_unitprime.X()<<", "<<n_p4M21inXFrame_unitprime.Y() <<", "<<n_p4M21inXFrame_unitprime.Z() <<" )"<<std::endl; 
  TLorentzVector n_p4Z2inXFrame( p4Z2 );
  n_p4Z2inXFrame.Boost( boostX );
  TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
  TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross( n_unitz_2 );//lepton #1 from Z2
  TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );
  TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
  phistar2 = n_p4PartoninZ2PlaneFrame_unitprime.Phi();
  double phistar12_0 = phistar1 + phistar2;
  if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
  else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
  else phistar12 = phistar12_0;

} //end calculateangles
*/
