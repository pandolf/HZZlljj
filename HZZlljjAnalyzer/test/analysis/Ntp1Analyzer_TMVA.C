#include "Ntp1Analyzer_TMVA.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"

#include "AnalysisElectron.h"
#include "AnalysisMuon.h"

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

  //btags:
  float trackCountingHighEffBJetTag;

};




HelicityLikelihoodDiscriminant::HelicityAngles computeHelicityAngles(TLorentzVector leptMinus, TLorentzVector leptPlus, TLorentzVector jet1, TLorentzVector jet2 );





Ntp1Analyzer_TMVA::Ntp1Analyzer_TMVA( const std::string& dataset, const std::string& jetChoice, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "TMVA", dataset, flags, tree ) {


  h1_nCounter_Zee_ = new TH1D("nCounter_Zee", "", 1, 0., 1.);
  h1_nCounter_Zmumu_ = new TH1D("nCounter_Zmumu", "", 1, 0., 1.);

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
  reducedTree_->Branch("leptTypeMC",&leptTypeMC_,"leptTypeMC_/I");
  reducedTree_->Branch("leptType",&leptType_,"leptType_/I");

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

  reducedTree_->Branch("kinfit_chiSquare",  &kinfit_chiSquare_,  "kinfit_chiSquare_/F");
  reducedTree_->Branch("kinfit_chiSquareProb",  &kinfit_chiSquareProb_,  "kinfit_chiSquareProb_/F");

  reducedTree_->Branch("pfMet",  &epfMet_,  "epfMet_/F");

  reducedTree_->Branch("helicityLD", &helicityLD_, "helicityLD_/F");
  reducedTree_->Branch("helicityLD_kinFit", &helicityLD_kinFit_, "helicityLD_kinFit_/F");

  reducedTree_->Branch("nBTags", &nBTags_, "nBTags_/I");
  reducedTree_->Branch("nBTagsLoose", &nBTagsLoose_, "nBTagsLoose_/I");

  h1_mZjj = new TH1F("mZjj", "", 50, 0., 300.);
  h1_mZjj_matched = new TH1F("mZjj_matched", "", 50, 0., 300.);

  h1_mZjj_bestZ = new TH1F("mZjj_bestZ", "", 50, 0., 300.);
  h1_mZjj_bestH = new TH1F("mZjj_bestH", "", 50, 0., 300.);
  h1_mZjj_closestPair = new TH1F("mZjj_closestPair", "", 50, 0., 300.);

} 



Ntp1Analyzer_TMVA::~Ntp1Analyzer_TMVA() {

  outfile_->cd();
  h1_nCounter_Zee_->Write();
  h1_nCounter_Zmumu_->Write();
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

   // to fix Z->ee BR bug, compute number of events having Z->ee and Z->mumu
   int nCounterZee = (isMC_) ? fChain->GetEntries("statusMc==3 && idMc==11") : 0;
   int nCounterZmumu = (isMC_) ? fChain->GetEntries("statusMc==3 && idMc==13") : 0;
   h1_nCounter_Zee_->SetBinContent( 1, nCounterZee );
   h1_nCounter_Zmumu_->SetBinContent( 1, nCounterZmumu );

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

     std::vector<TLorentzVector> electronsMC;
     std::vector<TLorentzVector> muonsMC;

     for( unsigned iMc=0; iMc<nMc; ++iMc ) {

       // partons only
       if( statusMc[iMc] != 3 ) continue;

       TLorentzVector* thisParticle = new TLorentzVector();
       thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

       if( fabs(idMc[iMc])==11 && idMc[mothMc[iMc]]==23 ) electronsMC.push_back( *thisParticle );
       if( fabs(idMc[iMc])==13 && idMc[mothMc[iMc]]==23 ) muonsMC.push_back( *thisParticle );

       delete thisParticle;
       thisParticle = 0;

     }


     if( muonsMC.size() > 0 ) leptTypeMC_ = 0;
     else if( electronsMC.size() > 0 ) leptTypeMC_ = 1;






     // -----------------------------
     //      FROM NOW ON RECO
     // -----------------------------


     epfMet_ = energyPFMet[0];
     //phipfMet_ = phiPFMet[0];


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

       thisMuon.pixelHits = numberOfValidPixelBarrelHitsTrack[trackIndexMuon[iMuon]]+numberOfValidPixelEndcapHitsTrack[trackIndexMuon[iMuon]];
       thisMuon.trackerHits = trackValidHitsTrack[trackIndexMuon[iMuon]];


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
       thisEle.sigmaIetaIeta = (superClusterIndexEle[iEle]>=0) ? sqrt(covIEtaIEtaSC[superClusterIndexEle[iEle]]) : sqrt(covIEtaIEtaSC[PFsuperClusterIndexEle[iEle]]);
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

       thisJet.trackCountingHighEffBJetTag =  trackCountingHighEffBJetTagsAK5PFPUcorrJet[iJet];

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
     QGLikelihoodJetRecoil_ = (jetRecoil.Pt()>15.) ? qglikeli->computeQGLikelihood( jetRecoil.Pt(), jetRecoil.nCharged, jetRecoil.nNeutral, jetRecoil.ptD, -1. ) : -1.;

     QGLikelihoodJet1Jet2_ = QGLikelihoodJet1_*QGLikelihoodJet2_;
     QGLikelihoodJet1Jet2Recoil_ = (QGLikelihoodJetRecoil_>=0.) ? QGLikelihoodJet1_*QGLikelihoodJet2_*QGLikelihoodJetRecoil_ : QGLikelihoodJet1_*QGLikelihoodJet2_;


     // fill mZjj before the kinfit:
     TLorentzVector diJet_pre = jet1 + jet2;
     if( diJet_pre.M() < 50. || diJet_pre.M() > 150. ) continue; //loose cut to avoid kinfit crashes
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


     // ------------------------
     //   KINEMATIC FIT: BEGIN
     // ------------------------
     
     DiJetKinFitter* fitter_jets = new DiJetKinFitter( "fitter_jets", "fitter_jets", Zmass );
     std::pair<TLorentzVector,TLorentzVector> jets_kinfit = fitter_jets->fit(jet1, jet2);
     TLorentzVector jet1_kinfit(jets_kinfit.first);
     TLorentzVector jet2_kinfit(jets_kinfit.second);

     kinfit_chiSquare_ = fitter_jets->getS()/fitter_jets->getNDF();
     kinfit_chiSquareProb_ = TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF());

  //   jet1.SetPtEtaPhiE( jet1_kinfit.Pt(), jet1_kinfit.Eta(), jet1_kinfit.Phi(), jet1_kinfit.E() );
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

     
     bool twotags_loose  = (jet1.trackCountingHighEffBJetTag>=4. && jet2.trackCountingHighEffBJetTag>1.85)
                        || (jet1.trackCountingHighEffBJetTag>1.85 && jet2.trackCountingHighEffBJetTag>1.85);
     bool onetag_loose  = ( !twotags_loose ) && ( jet1.trackCountingHighEffBJetTag>1.85 || jet2.trackCountingHighEffBJetTag>1.85 );
     bool zerotags_loose = ( !twotags_loose && !onetag_loose );

     if( twotags_loose ) nBTagsLoose_ = 2;
     else if( onetag_loose ) nBTagsLoose_ = 1;
     else if( zerotags_loose ) nBTagsLoose_ = 0;



     bool twotags = jet1.trackCountingHighEffBJetTag>=4. && jet2.trackCountingHighEffBJetTag>=4.;
     bool onetag  = (jet1.trackCountingHighEffBJetTag>=4. && jet2.trackCountingHighEffBJetTag<4.)
                 || (jet1.trackCountingHighEffBJetTag<4. && jet2.trackCountingHighEffBJetTag>=4.);
     bool zerotags = jet1.trackCountingHighEffBJetTag<4. && jet2.trackCountingHighEffBJetTag<4.;
     

     if( twotags ) nBTags_ = 2;
     else if( onetag ) nBTags_ = 1;
     else if( zerotags ) nBTags_ = 0;

     reducedTree_->Fill(); 

   } //for entries


// std::cout << "BEST Z: " << std::endl;
// std::cout << "Correct Pairs: " << nCorrectPairs_bestZ << " (" << (float)nCorrectPairs_bestZ/((float)nCorrectPairs_bestZ+(float)nIncorrectPairs_bestZ)*100. << "%)" << std::endl;
// std::cout << "Incorrect Pairs: " << nIncorrectPairs_bestZ << " (" << (float)nIncorrectPairs_bestZ/((float)nCorrectPairs_bestZ+(float)nIncorrectPairs_bestZ)*100. << "%)" << std::endl;

// std::cout << std::endl << "BEST H: " << std::endl;
// std::cout << "Correct Pairs: " << nCorrectPairs_bestH << " (" << (float)nCorrectPairs_bestH/((float)nCorrectPairs_bestH+(float)nIncorrectPairs_bestH)*100. << "%)" << std::endl;
// std::cout << "Incorrect Pairs: " << nIncorrectPairs_bestH << " (" << (float)nIncorrectPairs_bestH/((float)nCorrectPairs_bestH+(float)nIncorrectPairs_bestH)*100. << "%)" << std::endl;

// std::cout << std::endl << "CLOSEST PAIR:" << std::endl;
// std::cout << "Correct Pairs: " << nCorrectPairs_closestPair << " (" << (float)nCorrectPairs_closestPair/((float)nCorrectPairs_closestPair+(float)nIncorrectPairs_closestPair)*100. << "%)" << std::endl;
// std::cout << "Incorrect Pairs: " << nIncorrectPairs_closestPair << " (" << (float)nIncorrectPairs_closestPair/((float)nCorrectPairs_closestPair+(float)nIncorrectPairs_closestPair)*100. << "%)" << std::endl;

} //loop



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
