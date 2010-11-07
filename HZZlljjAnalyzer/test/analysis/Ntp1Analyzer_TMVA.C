#include "Ntp1Analyzer_TMVA.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"







Ntp1Analyzer_TMVA::Ntp1Analyzer_TMVA( const std::string& dataset, const std::string& jetChoice, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "TMVA", dataset, flags, tree ) {


  jetChoice_ = jetChoice;


} //constructor



void Ntp1Analyzer_TMVA::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");

  reducedTree_->Branch("ptLept1",  &ptLept1_,  "ptLept1_/F");
  reducedTree_->Branch("absetaLept1",  &absetaLept1_,  "absetaLept1_/F");

  reducedTree_->Branch("ptLept2",  &ptLept2_,  "ptLept2_/F");
  reducedTree_->Branch("absetaLept2",  &absetaLept2_,  "absetaLept2_/F");

  reducedTree_->Branch("mZll",  &mZll_,  "mZll_/F");
  reducedTree_->Branch("ptZll",  &ptZll_,  "ptZll_/F");
  reducedTree_->Branch("deltaRll",  &deltaRll_,  "deltaRll_/F");

  reducedTree_->Branch( "ptJet1",  &ptJet1_,  "ptJet1_/F");
  reducedTree_->Branch("absetaJet1", &absetaJet1_, "absetaJet1_/F");

  reducedTree_->Branch( "ptJet2",  &ptJet2_,  "ptJet2_/F");
  reducedTree_->Branch("absetaJet2", &absetaJet2_, "absetaJet2_/F");

  reducedTree_->Branch( "ptJetRecoil",  &ptJetRecoil_,  "ptJetRecoil_/F");
  reducedTree_->Branch("absetaJetRecoil", &absetaJetRecoil_, "absetaJetRecoil_/F");
  reducedTree_->Branch("deltaR_recoil_jet1", &deltaR_recoil_jet1_, "deltaR_recoil_jet1_/F");
  reducedTree_->Branch("deltaR_recoil_Zjj", &deltaR_recoil_Zjj_, "deltaR_recoil_Zjj_/F");
  reducedTree_->Branch("deltaR_recoil_Higgs", &deltaR_recoil_Higgs_, "deltaR_recoil_Higgs_/F");

  reducedTree_->Branch("mZjj",  &mZjj_,  "mZjj_/F");
  reducedTree_->Branch("ptZjj",  &ptZjj_,  "ptZjj_/F");
  reducedTree_->Branch("deltaRjj",  &deltaRjj_,  "deltaRjj_/F");

  reducedTree_->Branch("deltaRZZ",  &deltaRZZ_,  "deltaRZZ_/F");
  reducedTree_->Branch("deltaAbsEtaZZ",  &deltaAbsEtaZZ_,  "deltaAbsEtaZZ_/F");
  reducedTree_->Branch("absDeltaEtaZZ",  &absDeltaEtaZZ_,  "absDeltaEtaZZ_/F");
  reducedTree_->Branch("deltaPhiZZ",  &deltaPhiZZ_,  "deltaPhiZZ_/F");
  reducedTree_->Branch("ptZZ",  &ptZZ_,  "ptZZ_/F");
  reducedTree_->Branch("absEtaZZ",  &absEtaZZ_,  "absEtaZZ_/F");

  reducedTree_->Branch("pfMet",  &epfMet_,  "epfMet_/F");

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

     TRegexp re("JHUgen_HiggsSM");
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

     ptLept1_ = leptons[0].Pt();
     absetaLept1_ = fabs(leptons[0].Eta());
     
     ptLept2_ = leptons[1].Pt();
     absetaLept2_ = fabs(leptons[1].Eta());

     TLorentzVector diLepton = leptons[0] + leptons[1];

     mZll_ = diLepton.M();
     ptZll_ = diLepton.Pt();
     deltaRll_ = leptons[0].DeltaR(leptons[1]);




     // ------------------
     // JETS
     // ------------------

     float jetPt_thresh = 30.;


     // first save leading jets in event:
     std::vector<TLorentzVector> leadJets;

     for( unsigned int iJet=0; iJet<nAK5PFJet; ++iJet ) {

       TLorentzVector thisJet( pxAK5PFJet[iJet], pyAK5PFJet[iJet], pzAK5PFJet[iJet], energyAK5PFJet[iJet] );

       // save at least 3 lead jets (if event has them) and all jets with pt>thresh:
       if( leadJets.size()>=3 && thisJet.Pt()<jetPt_thresh ) break;

       // far away from leptons:
       if( thisJet.DeltaR( leptons[0] ) <= 0.5 ) continue;
       if( thisJet.DeltaR( leptons[1] ) <= 0.5 ) continue;

       leadJets.push_back(thisJet);

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


       for( unsigned int jJet=iJet+1; jJet<leadJets.size(); ++jJet ) {

         TLorentzVector otherJet = leadJets[jJet];

         // --------------
         // kinematics:
         // --------------
         if( otherJet.Pt() < jetPt_thresh ) continue;

         if( nPairs_>=50 ) {
        
           std::cout << "MORE than 50 jet pairs found. SKIPPING!!" << std::endl;

         } else {

   //      ptJet1_[nPairs_] = leadJets[iJet].Pt();
   //      absetaJet1_[nPairs_] = fabs(leadJets[iJet].Eta());
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
       std::cout << "leadjets size: " << leadJets.size() << " best_i: " << best_i_bestZ << " best_j: " << best_j_bestZ << std::endl;
       continue; //means that less than 2 jets were found
     }

     TLorentzVector jet1, jet2, jetRecoil;

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


     TLorentzVector diJet_bestZ = leadJets[best_i_bestZ] + leadJets[best_j_bestZ];
     TLorentzVector diJet_bestH = leadJets[best_i_bestH] + leadJets[best_j_bestH];
     TLorentzVector diJet_closestPair = leadJets[best_i_closestPair] + leadJets[best_j_closestPair];

     h1_mZjj_bestZ->Fill( diJet_bestZ.M() );
     h1_mZjj_bestH->Fill( diJet_bestH.M() );
     h1_mZjj_closestPair->Fill( diJet_closestPair.M() );

     ptJet1_ = jet1.Pt();
     absetaJet1_ = fabs(jet1.Eta());
     
     ptJet2_ = jet2.Pt();
     absetaJet2_ = fabs(jet2.Eta());


     ptJetRecoil_ = (jetRecoil.E()>0.) ? jetRecoil.Pt() : 0.;
     absetaJetRecoil_ = (jetRecoil.E()>0.) ? fabs(jetRecoil.Eta()) : -1.;
     deltaR_recoil_jet1_ = (jetRecoil.E()>0.) ? jetRecoil.DeltaR(jet1) : -1.;

     TLorentzVector diJet = jet1 + jet2;
     TLorentzVector ZZ = diJet + diLepton;

     mZjj_ = diJet.M();

     deltaR_recoil_Zjj_ = (jetRecoil.E()>0.) ? jetRecoil.DeltaR(diJet) : -1.;
     deltaR_recoil_Higgs_ = (jetRecoil.E()>0.) ? jetRecoil.DeltaR(ZZ) : -1.;

     deltaRZZ_ = diLepton.DeltaR(diJet);
     deltaAbsEtaZZ_ = fabs(diLepton.Eta()) - fabs(diJet.Eta());
     absDeltaEtaZZ_ = fabs(diLepton.Eta() - diJet.Eta());
     deltaPhiZZ_ = diJet.DeltaPhi(diLepton);
     ptZZ_ = ZZ.Pt();
     absEtaZZ_ = fabs(ZZ.Eta());

     
     if( ZZ.M()<450. && ZZ.M()>350. )
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



