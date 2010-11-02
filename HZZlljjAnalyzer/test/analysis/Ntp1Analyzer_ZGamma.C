#include "Ntp1Analyzer_ZGamma.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"





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

  reducedTree_->Branch("nJets_total",  &nJets_total_,  "nJets_total_/I");

  reducedTree_->Branch("nJet",  &nJet_,  "nJet_/I");
  reducedTree_->Branch("eJet",  &eJet_,  "eJet_[nJet_]/F");
  reducedTree_->Branch("ptJet",  &ptJet_,  "ptJet_[nJet_]/F");
  reducedTree_->Branch("etaJet",  &etaJet_,  "etaJet_[nJet_]/F");
  reducedTree_->Branch("phiJet",  &phiJet_,  "phiJet_[nJet_]/F");


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


     bool noLeptons = false;
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
       noLeptons = true;
     }


     if( noLeptons ) continue;

     TLorentzVector ZGamma;
     ZGamma.SetPtEtaPhiE( pMc[zIndexll]*sin(thetaMc[zIndexll]), etaMc[zIndexll], phiMc[zIndexll], energyMc[zIndexll] );

     ptZGamma_  = ZGamma.Pt();
     eZGamma_   = ZGamma.Energy();
     etaZGamma_ = ZGamma.Eta();
     phiZGamma_ = ZGamma.Phi();

     nJet_ = nJets_total_ = 0;

     std::vector<AnalysisJet> jets;

     for( unsigned int iJet=0; iJet<nAK5GenJet; ++iJet ) {

       AnalysisJet thisJet( pxAK5GenJet[iJet], pyAK5GenJet[iJet], pzAK5GenJet[iJet], energyAK5GenJet[iJet] );

       // far away from leptons:
       if( thisJet.DeltaR(lept1MC) <= 0.5 ) continue;
       if( thisJet.DeltaR(lept2MC) <= 0.5 ) continue;

       nJets_total_++;

       if( nJet_<10 ) {
         ptJet_[nJet_] = thisJet.Pt();
         etaJet_[nJet_] = thisJet.Eta();
         phiJet_[nJet_] = thisJet.Phi();
         eJet_[nJet_] = thisJet.Energy();
         nJet_++;
       }

     } //for jets


     reducedTree_->Fill(); 


   } //for entries

} //loop



