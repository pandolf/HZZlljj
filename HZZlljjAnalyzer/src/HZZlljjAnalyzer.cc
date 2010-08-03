// -*- C++ -*-
// //
// // Package:    HZZlljjAnalyzer
// // Class:      HZZlljjAnalyzer
// // 
// /**\class HZZlljjAnalyzer HZZlljjAnalyzer.cc HZZlljj/HZZlljjAnalyzer/src/HZZlljjAnalyzer.cc
//
//  Description: [one line class summary]
//
//   Implementation:
//        [Notes on implementation]
//        */
//        //
//        // Original Author:  Francesco Pandolfi,32 4-C03,+41227672087,
//        //         Created:  Wed Jul 28 14:59:09 CEST 2010
//        // $Id$
//        //
//        //
//

//
// constructors and destructor
//

#include "HZZlljj/HZZlljjAnalyzer/interface/HZZlljjAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "MyAnalysis/IsolationTools/interface/SuperClusterHitsEcalIsolation.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"

// HLT trigger
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#define MAXHLTBITS    200

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TLorentzVector.h"

#include <set>
#include <algorithm>

using namespace reco;

// Signed version of delta_phi
inline float HZZlljjAnalyzer::delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}

// Mirror-symmetric delta_eta
inline float HZZlljjAnalyzer::delta_eta(float eta1, float eta2) {

  return (eta2 >= 0 ? eta1 - eta2 : eta2 - eta1);
}

inline double HZZlljjAnalyzer::oplus(double a, double b) {
  return sqrt(a*a + b*b);
}


HZZlljjAnalyzer::HZZlljjAnalyzer(const edm::ParameterSet& iConfig)
{
  _debug = iConfig.getParameter<bool>("debug");
  MCTruthCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("MCTruthCollection");
  triggerTag_ = iConfig.getUntrackedParameter<edm::InputTag>("TriggerTag");
  trackTags_ = iConfig.getUntrackedParameter<edm::InputTag>("tracks");
  Vertexsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("vertices");
  Electronsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("Electronsrc");
  Muonsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("Muonsrc");
  JetAlgosrc_ = iConfig.getUntrackedParameter<std::string>("jetalgo");
  METsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("met");

  METGensrc_ = iConfig.getUntrackedParameter<edm::InputTag>("genMet");
//HBhitsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("hbhits");
  ecalHitsCollection_    = iConfig.getParameter<std::string>("ecalHitsCollection");
  ecalHitsProducer_      = iConfig.getParameter<std::string>("ecalHitsProducer");
//JetCorrector_ = iConfig.getParameter<std::string>("JetCorrectionService");
  genjetptthr_ = iConfig.getParameter<double>("genjetptthr");
  jetptthr_ = iConfig.getParameter<double>("jetptthr");
  jetnmin_ = iConfig.getParameter<int>("jetnmin");
 
}


HZZlljjAnalyzer::~HZZlljjAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HZZlljjAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   nMC = nJet = nJetGen = 0;
   nMu = nEle = 0;
   nSIM = nPF = 0;

//   using reco::TrackCollection;
  

   isMC = !iEvent.isRealData(); // separate MC processing
   store = iEvent.eventAuxiliary().storeNumber(); // study stability across a store
   lbn = iEvent.luminosityBlock(); // sum LBN lumi for normalization
   bx = iEvent.bunchCrossing(); // study effect of out-of-time pile-up
   orbit = iEvent.orbitNumber(); // study beam heating with time (longer bunches)
   run = iEvent.id().run(); // unique ID - part 1
   event = iEvent.id().event(); // unique ID - part 2


   Handle<GenEventInfoProduct> hEventInfo;
   if( isMC ) iEvent.getByLabel("generator", hEventInfo);


   // ------ MC INFORMATION:

   // get MC info from GenParticleCandidates 
   Handle<GenParticleCollection> genParticles;
   if( isMC ) iEvent.getByLabel("genParticles", genParticles);
   
   // get GEANT sim tracks and vertices (includes conversions)
   Handle<SimTrackContainer> simTracks_h;
   const SimTrackContainer* simTracks;
   if( isMC ) iEvent.getByLabel("g4SimHits", simTracks_h);
   simTracks = (simTracks_h.isValid()) ? simTracks_h.product() : 0;
   
   Handle<SimVertexContainer> simVert_h;
   const SimVertexContainer* simVertices;
   if( isMC ) iEvent.getByLabel("g4SimHits", simVert_h);
   simVertices = (simVert_h.isValid()) ? simVert_h.product() : 0;
   

// edm::Handle<DcsStatusCollection> dcsHandle;
// iEvent.getByLabel("scalersRawToDigi", dcsHandle);
// double evt_bField;
// if( !isMC ) {
//   float currentToBFieldScaleFactor = 2.09237036221512717e-04;
//   float current = (*dcsHandle)[0].magnetCurrent();
//   evt_bField = current*currentToBFieldScaleFactor;
// }  else {
//   ESHandle<MagneticField> magneticField;
//   iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
//   evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
// }



   // ------ RECO INFORMATION:

   // get tracks
   Handle<TrackCollection> tracks;
   iEvent.getByLabel(trackTags_,tracks);
   
   // get primary vertices
   //Handle<VertexCollection> VertexHandle;
   Handle<vector<Vertex> > VertexHandle;
   //iEvent.getByLabel("offlinePrimaryVertices", VertexHandle);
   iEvent.getByLabel(Vertexsrc_, VertexHandle);
   
// // get photons
// Handle<PhotonCollection>  PhotonHandle;
// iEvent.getByLabel(Photonsrc_, PhotonHandle);

   // get muons
   Handle<MuonCollection>  muons;
   iEvent.getByLabel(Muonsrc_, muons);

   // get electrons
   Handle<GsfElectronCollection>  electrons;
   iEvent.getByLabel(Electronsrc_, electrons);

   // get PFCandidates
   Handle<PFCandidateCollection>  PFCandidates;
   iEvent.getByLabel("particleFlow", PFCandidates);

   // get (PF) jets collection
   std::string jetcollection = JetAlgosrc_ + "PFJets";
   Handle<PFJetCollection> jets;
   iEvent.getByLabel((edm::InputTag)(jetcollection), jets);

   //get jet correctors
   std::string correctedjets = JetAlgosrc_ + "PFL2L3";
   const JetCorrector* jetcorrector = JetCorrector::getJetCorrector (correctedjets, iSetup);
 
   // get gen jet collection
   std::string genjetcollection = JetAlgosrc_ + "GenJets";
   Handle<GenJetCollection> jetsgen;
   if( isMC ) iEvent.getByLabel((edm::InputTag)(genjetcollection), jetsgen);

   // get pfMET
   Handle<PFMETCollection> pfmethandle;
   iEvent.getByLabel(METsrc_, pfmethandle);

   // get gen MET
   Handle<GenMETCollection> genmethandle;
   if( isMC ) iEvent.getByLabel(METGensrc_, genmethandle);

   Handle<GenMETCollection> genmethandle2;
   if( isMC ) iEvent.getByLabel("genMetCalo", genmethandle2);
  
// // get HCAL info
// Handle<HBHERecHitCollection> hbhe;
// iEvent.getByLabel(HBhitsrc_, hbhe);
// const HBHERecHitMetaCollection mhbhe(*hbhe);

   // get ECAL reco hits
   Handle<EBRecHitCollection> ecalhits;
   const EBRecHitCollection* rhits=0;
   iEvent.getByLabel(ecalHitsProducer_, ecalHitsCollection_, ecalhits);
   const EcalRecHitMetaCollection mecalhits(*ecalhits);    
   rhits = ecalhits.product(); // get a ptr to the product

   Handle<EERecHitCollection> ecalhitsee;
   const EERecHitCollection* rhitsee=0;
   iEvent.getByLabel(ecalHitsProducer_, "EcalRecHitsEE", ecalhitsee);
   rhitsee = ecalhits.product(); // get a ptr to the product

// // get geometry
// edm::ESHandle<CaloGeometry> geoHandle;
// //   iSetup.get<IdealGeometryRecord>().get(geoHandle);
// iSetup.get<CaloGeometryRecord>().get(geoHandle);
// const CaloGeometry* geometry = geoHandle.product();
// const CaloSubdetectorGeometry* geometry_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

   // get topology
// const CaloSubdetectorTopology *topology_p;
// edm::ESHandle<CaloTopology> topoHandle;
// iSetup.get<CaloTopologyRecord>().get(topoHandle);
// topology_p = topoHandle->getSubdetectorTopology(DetId::Ecal, EcalBarrel);

// edm::ESHandle<CaloTopology> pTopology;
// iSetup.get<CaloTopologyRecord>().get(pTopology);
// const CaloTopology *topology = pTopology.product();

   //   edm::ESHandle<TrackerGeometry> trackerHandle_;
   //edm::ESHandle<MagneticField> theMagField;
   //iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
   //   iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle_);
  
   //ClusterShapeAlgo algo;

//---------------------------HLT Trigger ---------------------------------------------------------------------------------------------
// You Can See HLT Name list ->  " JetMETCorrections/GammaJet/test/HLTList.txt " file 

   for (unsigned int iHLT = 0; iHLT != MAXHLTBITS; ++iHLT) 
   {
   	aHLTResults[iHLT] = false;
   }
   hltPass = false;

   strcpy(aHLTNames,"");
   hltNamesLen = 0;

   edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
   iEvent.getByLabel(triggerTag_, hltTriggerResultHandle);

   if (!hltTriggerResultHandle.isValid()) 
   {
   	std::cout << "invalid handle for HLT TriggerResults" << std::endl;
   } 
   else 
   {
   	edm::TriggerNames HLTNames;
//   	HLTNames.init(*hltTriggerResultHandle);
      HLTNames = iEvent.triggerNames(*hltTriggerResultHandle);

   	std::string tempnames;  
   	hltCount = hltTriggerResultHandle->size();
   	//std::cout << "hltTriggerResult->size(): " << hltCount << std::endl;

   	for (int i = 0 ; i != hltCount; ++i) {

	  tempnames += HLTNames.triggerName(i) + ":";
	  //aHLTResults[i] = hltTriggerResultHandle->accept(i);
	  //cout << i <<"....." << HLTNames.triggerName(i).c_str() << ".... : " << hltTriggerResultHandle->accept(i) << endl;

	  map<string, int>::const_iterator it 
	    = hltTriggers.find(HLTNames.triggerName(i));
	  if (it != hltTriggers.end()) {
	    int itrig = it->second;
	    aHLTResults[itrig] = hltTriggerResultHandle->accept(i);
	    hltPass |= aHLTResults[itrig];
	  }
   	} // for i

   	hltNamesLen = tempnames.length();
   	strcpy(aHLTNames,tempnames.c_str());

   } // HLT isValid
//--------------------------------------------------------------------------------------------------------------------------------------

   // Loop over MC truth

   genpt = 0.;

   if( isMC ) {
     
     //   genpt = *genEventScale;   
     if (hEventInfo->binningValues().size() > 0)
       genpt = hEventInfo->binningValues()[0];

     // Figure out the true vertex from the partons:
     // Each event has exactly eight partons (status==3);
     // 0-1: incoming protons (pdgId==2212),
     // 2-3: ISR partons (|pdgId|<=21, udsg only, no c/b),
     // 4-5: incoming partons (|pdgId|<=21, udscg, no b?),
     // 6-7: outgoing partons (|pdgId|<=22, udscg, photons)
     // Every parton 2-7 has the same production vertex (0 for 0-1)
     //assert(genParticles->at(2).status()==3);
     //assert(fabs(genParticles->at(2).pdgId())<100); //<25-><100 to include Z,W
     vxMC = genParticles->at(2).vx();
     vyMC = genParticles->at(2).vy();
     vzMC = genParticles->at(2).vz();

     // Momentum conservation in px (py, pt?): 
     // outgoing + ISR - incoming ~ 0 (couple of GeV)
     // Is the small discrepancy due to the virtual mass of the partons?
     set<int> mothers;
     map<const GenParticle*, int> mapMC;

     for (GenParticleCollection::const_iterator p = genParticles->begin();
          p != genParticles->end(); ++p) {
       
       if (nMC>=(nMaxMC-1)) {continue;}  // to reduce the root file size
       // need space for the mom, eventually

       // Select only a subset of particles to reduce size:
       // All the partons (8)
       // All the stable (status==1) particles within deltaR<1.0
       // around outgoing partons
       // Parents (status==2) of the stable photons and electrons in above
     
       double deltaR = 0;
       if (p->status() == 1) {
         double deltaR1 = oplus(delta_eta(p->eta(),etaMC[6]),
                          delta_phi(p->phi(),phiMC[6]));
         double deltaR2 = oplus(delta_eta(p->eta(),etaMC[7]),
                          delta_phi(p->phi(),phiMC[7]));
         deltaR = min(deltaR1, deltaR2);
       }

       // Neutral particles kept with >200 MeV (ECAL ZS threshold)
       // Charged particles kept with >75 MeV (tracking threshold)
     //if (p->status()==3 || (p->status()==1 && deltaR<1.0 &&
     //                      (p->pt()>0.200 ||
     //                      (p->charge()!=0 && p->pt()>0.075))) ) {

       if (p->status()==3 || (p->status()==1 && 
                             (p->pt()>0.200 ||
                             (p->charge()!=0 && p->pt()>0.075))) ) {

         pdgIdMC[nMC] = p->pdgId();
         statusMC[nMC] = p->status();
         //massMC[nMC] = p->mass();
         ptMC[nMC] = p->pt();
         eMC[nMC] = p->energy();	 
         etaMC[nMC] = p->eta();	 
         phiMC[nMC] = p->phi();	 
         
         if (p->numberOfMothers() > 0) { 
           const Candidate * mom = p->mother();
           for (size_t j = 0; j != genParticles->size(); ++j) {
             const Candidate * ref = &((*genParticles)[j]);
             //if (mom->px() == ref->px() && mom->py() == ref->py()
             //&& mom->pz() == ref->pz() && mom->status() == ref->status()
             //&& mom->pdgId()==ref->pdgId()) {
           
               //assert(mom==ref); // address of the candidate is the same?
               //above works in about 99.7% of events
           
             if (mom==ref) {
               motherIDMC[nMC] = j;
               //motherIDMC = j;
             } //if mom ref
           } //for
         } //if

         mapMC[&(*p)] = nMC;
         ++nMC; 

      // // if stable photon/electron, find parent
      // if (p->status() == 1 && motherIDMC[nMC] != -1
      //     && (p->pdgId() == kPhoton || p->pdgId() == kElectron)) {
      //     
      //   //const Candidate * mom = p->mother();
      //   const GenParticle *mom = (const GenParticle*)p->mother();
      //   if (mom->status() == 2
      //       && (mom->pdgId()<81 || mom->pdgId()>100) // remove MC internal 
      //       && mothers.find(motherIDMC[nMC]) == mothers.end()) {

      //   mothers.insert(motherIDMC[nMC]);
      //   
      //   if (nMC>=nMaxMC) {continue;}  // to reduce the root file size
      //   pdgIdMC[nMC] = mom->pdgId();
      //   statusMC[nMC] = mom->status();
      //   //massMC[nMC] = mom->mass();
      //   ptMC[nMC] = mom->pt();
      //   eMC[nMC] = mom->energy();
      //   etaMC[nMC] = mom->eta();
      //   phiMC[nMC] = mom->phi(); 

      //   mapMC[mom] = nMC;
      //   ++nMC; 
      // }
      // } // stable photon has parent
      
       } // if keep particle

     } // loop particles

   
     //const double genjetptthr = 5.; // already implicit in GenJet reco
     //const int genjetnmin = 4;
     // Loop over gen Jets

     for (GenJetCollection::const_iterator it = jetsgen->begin(); 
       it != jetsgen->end(); ++it) {
     
       if (nJetGen>=(nMaxJet-1)) {cout << "number of gen jets is larger than nMaxJet = " << nMaxJet << ". Skipping" << endl; continue;}
       ptJetGen[nJetGen] = it->pt();	 
       eJetGen[nJetGen] = it->energy();	 
       etaJetGen[nJetGen] = it->eta();	 
       phiJetGen[nJetGen] = it->phi();	      
       
       nJetGen++;
     }


  } //if(isMC)



  //
  //  RECO
  //

   // Get the primary vertex coordinates
   nvertex = VertexHandle->size();

   for( unsigned i=0; i< (unsigned)(abs(nvertex)); ++i ) {

     if( i>9 ) { std::cout << "Number of vertexes is larger than 10. Skipping." << std::endl; continue;}

     vx[i] = (VertexHandle->at(i).isValid()) ? VertexHandle->at(i).x() : 999.;
     vy[i] = (VertexHandle->at(i).isValid()) ? VertexHandle->at(i).y() : 999.;
     vz[i] = (VertexHandle->at(i).isValid()) ? VertexHandle->at(i).z() : 999.;

     vntracks[i] = (VertexHandle->at(i).isValid()) ? VertexHandle->at(i).tracksSize() : 0;
     vchi2[i]    = (VertexHandle->at(i).isValid()) ? VertexHandle->at(i).normalizedChi2() : 100.;
     vndof[i]    = (VertexHandle->at(i).isValid()) ? VertexHandle->at(i).ndof() : 0.;

   } //for vertex


   //
   // loop over muons:
   //
   for (MuonCollection::const_iterator it = muons->begin();  it != muons->end(); ++it) {

     if( !(it->isGlobalMuon()) && !(it->isTrackerMuon()) ) continue;
     if (nMu>=(nMaxMu-1)) {cout << "number of reco muons is larger than nMaxMu = " << nMaxMu << ". Skipping" << endl; continue;}

       ptMu[nMu] = it->pt();
       eMu[nMu] = it->energy();	 
       etaMu[nMu] = it->eta();	 
       phiMu[nMu] = it->phi();	      
       
       isGlobalMu[nMu] = it->isGlobalMuon();
       isGoodMu[nMu] = muon::isGoodMuon((*it),muon::GlobalMuonPromptTight);
       isTrackerMu[nMu] = it->isTrackerMuon();
       nMatchesMu[nMu] = it->numberOfMatches();
       nValidHitsTrackMu[nMu] = ( &(*(it->innerTrack()))!=0 ) ? it->innerTrack()->numberOfValidHits() : 0;
       nPixelLayerMu[nMu] = ( &(*(it->innerTrack()))!=0 ) ? it->innerTrack()->hitPattern().pixelLayersWithMeasurement() : 0;
       dxyMu[nMu] = ( &(*(it->innerTrack()))!=0 ) ? it->innerTrack()->dxy() : 0.;

       nMu++;

   } //muons
     

   //
   // loop over electrons:
   //
   for (GsfElectronCollection::const_iterator it = electrons->begin();  it != electrons->end(); ++it) {

     if (nEle>=nMaxEle) {cout << "number of reco electrons is larger than nMaxEle = " << nMaxEle << ". Skipping" << endl; continue;}

       ptEle[nEle] = it->pt();
       eEle[nEle] = it->energy();	 
       etaEle[nEle] = it->eta();	 
       phiEle[nEle] = it->phi();	      
       
       trkIso03Ele[nEle] = it->dr03TkSumPt();
       ecalIso03Ele[nEle] = it->dr03EcalRecHitSumEt();
       // ecal isolation with SC seed rechits removal
       SuperClusterHitsEcalIsolation scBasedIsolation(rhits,rhitsee);
       scBasedIsolation.setExtRadius(0.3);
       scBasedIsolation.excludeHalo(false);
       ecalIsoGT03Ele[nEle] = scBasedIsolation.getSum(iEvent,iSetup,&(*(it->superCluster())));;
       hcalIso03Ele[nEle] = it->dr03HcalTowerSumEt();

       hOverEmEle[nEle] = it->hadronicOverEm();
       deltaPhiSCTrkEle[nEle] = it->deltaPhiSuperClusterTrackAtVtx();
       deltaEtaSCTrkEle[nEle] = it->deltaEtaSuperClusterTrackAtVtx();
       sigmaIetaIeta[nEle] = it->sigmaIetaIeta();

       const reco::Track *el_track = (const reco::Track*)(it->gsfTrack().get());  
       const reco::HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
       nInnerHitsEle[nEle] = p_inner.numberOfHits();

//     ConversionFinder convFinder;
//     ConversionInfo convInfo = convFinder.getConversionInfo(*it, tracks, evt_bField);
//    
//     convDistEle[nEle] = convInfo.dist();
//     convDcotThetaEle[nEle] = convInfo.dcot();

       nEle++;

   } //electrons
     

   //
   // loop over PFCandidates::
   //
   for (PFCandidateCollection::const_iterator it = PFCandidates->begin();  it != PFCandidates->end(); ++it) {

     if (nPF>=nMaxPF) {cout << "number of reco PFCandidates is larger than nMaxPF = " << nMaxPF << ". Skipping" << endl; continue;}

       ptPF[nPF] = it->pt();
       ePF[nPF] = it->energy();	 
       etaPF[nPF] = it->eta();	 
       phiPF[nPF] = it->phi();	      
       pTypePF[nPF] = it->particleId();	      
       
       ++nPF;

   } //PFCandidates
     

   //
   // loop over jets:
   //
   for (PFJetCollection::const_iterator it = jets->begin();  it != jets->end(); ++it) {
    
     if (nJet>=(nMaxJet-1)) {cout << "number of reco jets is larger than nMaxJet = " << nMaxJet << ". Skipping" << endl; continue;}

     if (nJet < jetnmin_ || it->pt() > jetptthr_) {

       ptJet[nJet] = it->pt();
       eJet[nJet] = it->energy();	 
       etaJet[nJet] = it->eta();	 
       phiJet[nJet] = it->phi();	      
       
       // Jet Energy Scale Corrections on-the-fly     
       PFJet  correctedJet = *it;
       double scale = jetcorrector->correction(it->p4());
       correctedJet.scaleEnergy(scale);
       ptCorrJet[nJet] = correctedJet.pt();

       // jet composition variables:
       nChargedHadronsJet[nJet] =  it->chargedHadronMultiplicity();
       eChargedHadronsJet[nJet] = it->chargedHadronEnergyFraction();
       nPhotonsJet[nJet] =  it->photonMultiplicity();
       ePhotonsJet[nJet] =  it->photonEnergyFraction();
       nNeutralHadronsJet[nJet] =  it->neutralHadronMultiplicity();
       eNeutralHadronsJet[nJet] =  it->neutralHadronEnergyFraction();
       nElectronsJet[nJet] =  it->electronMultiplicity();
       eElectronsJet[nJet] =  it->electronEnergyFraction();
       nMuonsJet[nJet] =  it->muonMultiplicity();
       eMuonsJet[nJet] =  it->muonEnergyFraction();
       nHFHadronsJet[nJet] =  it->HFHadronMultiplicity();
       eHFHadronsJet[nJet] =  it->HFHadronEnergyFraction();
       nHFEMJet[nJet] =  it->HFEMMultiplicity();
       eHFEMJet[nJet] =  it->HFEMEnergyFraction();

       ++nJet;

     } // if >jetptthr     
   } // jets


     
   // Fill pfMET
   const PFMETCollection *pfmetcol = pfmethandle.product();
   PFMET const& pfmet = pfmetcol->front();

   spfMet = pfmet.sumEt();
   epfMet = pfmet.pt();
   phipfMet = pfmet.phi();

   sMetGen = 0.;
   eMetGen = 0.;
   phiMetGen = 0.;

   sMetGen2 = 0.;
   eMetGen2 = 0.;
   phiMetGen2 = 0.;

   if( isMC ) {
     // Fill gen MET

     const GenMETCollection *genmetcol = genmethandle.product();
     GenMET const& genmet = genmetcol->front();

     sMetGen = genmet.sumEt();
     eMetGen = genmet.energy();
     phiMetGen = genmet.phi();

     const GenMETCollection *genmetcol2 = genmethandle2.product();
     GenMET const& genmet2 = genmetcol2->front();

     sMetGen2 = genmet2.sumEt();
     eMetGen2 = genmet2.energy();
     phiMetGen2 = genmet2.phi();

   } //if is MC
   
   m_tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
HZZlljjAnalyzer::beginJob()
{

  //m_tree = fs_->make<TTree>("pippo","Analysis tree");
  outfile = TFile::Open("output.root", "RECREATE");
  outfile->mkdir("myanalysis");
  outfile->cd("myanalysis");

  m_tree = new TTree ("pippo","Analysis tree") ;
  //  m_tree->SetAutoSave (10000000) ;
  m_tree->Branch("genpt",&genpt,"genpt/F");

  m_tree->Branch("isMC",&isMC,"isMC/O");
  m_tree->Branch("store",&store,"store/I");
  m_tree->Branch("lbn",&lbn,"lbn/I");
  m_tree->Branch("bx",&bx,"bx/I");
  m_tree->Branch("orbit",&orbit,"orbit/I");
  m_tree->Branch("run",&run,"run/I");
  m_tree->Branch("event",&event,"event/I");

  // Problem: nMC==100 always, and sometimes last particle has very high pT
  // => could be losing interesting particles, even quarks/gluons (status==2)
  // hmmm, status==1 particles have pT less than ~4.5 GeV for pThat>500 => ok
  m_tree->Branch("nMC",&nMC,"nMC/I");
  m_tree->Branch("pdgIdMC",&pdgIdMC,"pdgIdMC[nMC]/I");
  m_tree->Branch("statusMC",&statusMC,"statusMC[nMC]/I");
  m_tree->Branch("motherIDMC",&motherIDMC,"motherIDMC[nMC]/I");
  // Most MC particles have mass, but why do photons (status=1 and 3) have mass?
  //m_tree->Branch("massMC ",&massMC ,"massMC[nMC]/F");
  //to a good approximation, m = sqrt(e^2 - (pt*cosh(eta))^2), when m>1e-6 GeV
  m_tree->Branch("ptMC ",&ptMC ,"ptMC[nMC]/F");
  m_tree->Branch("eMC  ",&eMC  ,"eMC[nMC]/F");
  m_tree->Branch("etaMC",&etaMC,"etaMC[nMC]/F");
  m_tree->Branch("phiMC",&phiMC,"phiMC[nMC]/F");

  m_tree->Branch("nSIM",&nSIM,"nSIM/I");
  m_tree->Branch("pdgIdSIM",&pdgIdSIM,"pdgIdSIM[nSIM]/I");
  m_tree->Branch("statusSIM",&statusSIM,"statusSIM[nSIM]/I");
  m_tree->Branch("ptSIM ",&ptSIM ,"ptSIM[nSIM]/F");
  m_tree->Branch("eSIM  ",&eSIM  ,"eSIM[nSIM]/F");
  m_tree->Branch("etaSIM",&etaSIM,"etaSIM[nSIM]/F");
  m_tree->Branch("phiSIM",&phiSIM,"phiSIM[nSIM]/F");
  m_tree->Branch("rSIM",&rSIM,"rSIM[nSIM]/F");
  m_tree->Branch("zSIM",&zSIM,"zSIM[nSIM]/F");

  m_tree->Branch("nPF",&nPF,"nPF/I");
  m_tree->Branch("pdgIdPF",&pdgIdPF,"pdgIdPF[nPF]/I");
  m_tree->Branch("ptPF ",&ptPF ,"ptPF[nPF]/F");
  m_tree->Branch("ePF  ",&ePF  ,"ePF[nPF]/F");
  m_tree->Branch("etaPF",&etaPF,"etaPF[nPF]/F");
  m_tree->Branch("phiPF",&phiPF,"phiPF[nPF]/F");
  m_tree->Branch("pTypePF",&pTypePF,"pTypePF[nPF]/F");

  m_tree->Branch("nMu",&nMu,"nMu/I");
  m_tree->Branch("ptMu ",&ptMu ,"ptMu[nMu]/F");
  m_tree->Branch("eMu  ",&eMu  ,"eMu[nMu]/F");
  m_tree->Branch("etaMu",&etaMu,"etaMu[nMu]/F");
  m_tree->Branch("phiMu",&phiMu,"phiMu[nMu]/F");

  m_tree->Branch("isGlobalMu", &isGlobalMu, "isGlobalMu[nMu]/O");
  m_tree->Branch("isGoodMu", &isGoodMu, "isGoodMu[nMu]/O");
  m_tree->Branch("isTrackerMu", &isTrackerMu, "isTrackerMu[nMu]/O");
  m_tree->Branch("nValidHitsTrackMu", &nValidHitsTrackMu, "nValidHitsTrackMu[nMu]/I");
  m_tree->Branch("nPixelLayerMu", &nPixelLayerMu, "nPixelLayerMu[nMu]/I");
  m_tree->Branch("nMatchesMu", &nMatchesMu, "nMatchesMu[nMu]/I");
  m_tree->Branch("dxyMu", &dxyMu, "dxyMu[nMu]/F");

  m_tree->Branch("nEle",&nEle,"nEle/I");
  m_tree->Branch("ptEle ",&ptEle ,"ptEle[nEle]/F");
  m_tree->Branch("eEle  ",&eEle  ,"eEle[nEle]/F");
  m_tree->Branch("etaEle",&etaEle,"etaEle[nEle]/F");
  m_tree->Branch("phiEle",&phiEle,"phiEle[nEle]/F");
  m_tree->Branch("trkIso03Ele", &trkIso03Ele, "trkIso03Ele[nEle]/F");
  m_tree->Branch("ecalIso03Ele", &ecalIso03Ele, "ecalIso03Ele[nEle]/F");
  m_tree->Branch("ecalIsoGT03Ele", &ecalIsoGT03Ele, "ecalIsoGT03Ele[nEle]/F");
  m_tree->Branch("hcalIso03Ele", &hcalIso03Ele, "hcalIso03Ele[nEle]/F");
  m_tree->Branch("hOverEmEle", &hOverEmEle, "hOverEmEle[nEle]/F");
  m_tree->Branch("deltaPhiSCTrkEle", &hOverEmEle, "hOverEmEle[nEle]/F");
  m_tree->Branch("deltaEtaSCTrkEle", &deltaEtaSCTrkEle, "deltaEtaSCTrkEle[nEle]/F");
  m_tree->Branch("sigmaIetaIeta", &sigmaIetaIeta, "sigmaIetaIeta[nEle]/F");
  m_tree->Branch("nInnerHitsEle", &nInnerHitsEle, "nInnerHitsEle[nEle]/F");
  m_tree->Branch("convDistEle", &convDistEle, "convDistEle[nEle]/F");
  m_tree->Branch("convDcotThetaEle", &convDcotThetaEle, "convDcotThetaEle[nEle]/F");


  m_tree->Branch("nJet",&nJet,"nJet/I");
  m_tree->Branch("ptJet ",&ptJet ,"ptJet[nJet]/F");
  m_tree->Branch("ptCorrJet ",&ptCorrJet ,"ptCorrJet[nJet]/F");
  m_tree->Branch("eJet  ",&eJet  ,"eJet[nJet]/F");
  m_tree->Branch("etaJet",&etaJet,"etaJet[nJet]/F");
  m_tree->Branch("phiJet",&phiJet,"phiJet[nJet]/F");

  m_tree->Branch("nChargedHadronsJet",nChargedHadronsJet,"nChargedHadronsJet[nJet]/I");
  m_tree->Branch("nPhotonsJet",       nPhotonsJet,       "nPhotonsJet[nJet]/I");
  m_tree->Branch("nMuonsJet",         nMuonsJet,         "nMuonsJet[nJet]/I");
  m_tree->Branch("nElectronsJet",     nElectronsJet,     "nElectronsJet[nJet]/I");
  m_tree->Branch("nNeutralHadronsJet",nNeutralHadronsJet,"nNeutralHadronsJet[nJet]/I");
  m_tree->Branch("nHFHadronsJet",     nHFHadronsJet,     "nHFHadronsJet[nJet]/I");
  m_tree->Branch("nHFEMJet",     nHFEMJet,     "nHFEMJet[nJet]/I");

  m_tree->Branch("eChargedHadronsJet",eChargedHadronsJet,"eChargedHadronsJet[nJet]/F");
  m_tree->Branch("ePhotonsJet",ePhotonsJet,"ePhotonsJet[nJet]/F");
  m_tree->Branch("eMuonsJet",eMuonsJet,"eMuonsJet[nJet]/F");
  m_tree->Branch("eElectronsJet",eElectronsJet,"eElectronsJet[nJet]/F");
  m_tree->Branch("eNeutralHadronsJet",eNeutralHadronsJet,"eNeutralHadronsJet[nJet]/F");
  m_tree->Branch("eHFHadronsJet",eHFHadronsJet,"eHFHadronsJet[nJet]/F");
  m_tree->Branch("eHFEMJet",eHFEMJet,"eHFEMJet[nJet]/F");


  m_tree->Branch("nJetGen",&nJetGen,"nJetGen/I");
  m_tree->Branch("ptJetGen",&ptJetGen, "ptJetGen[nJetGen]/F");
  m_tree->Branch("eJetGen",&eJetGen, "eJetGen[nJetGen]/F");
  m_tree->Branch("etaJetGen",&etaJetGen, "etaJetGen[nJetGen]/F");
  m_tree->Branch("phiJetGen",&phiJetGen, "phiJetGen[nJetGen]/F");


  m_tree->Branch("spfMet  ",&spfMet  ,"spfMet/F");
  m_tree->Branch("epfMet  ",&epfMet  ,"epfMet/F");
  m_tree->Branch("phipfMet",&phipfMet,"phipfMet/F");

  m_tree->Branch("sMetGen  ",&sMetGen  ,"sMetGen/F");
  m_tree->Branch("eMetGen  ",&eMetGen  ,"eMetGen/F");
  m_tree->Branch("phiMetGen",&phiMetGen,"phiMetGen/F");

  m_tree->Branch("sMetGen2  ",&sMetGen2  ,"sMetGen2/F");
  m_tree->Branch("eMetGen2  ",&eMetGen2  ,"eMetGen2/F");
  m_tree->Branch("phiMetGen2",&phiMetGen2,"phiMetGen2/F");


  //vertex info
  m_tree->Branch("nvertex",&nvertex,"nvertex/I");

  m_tree->Branch("vxMC",&vxMC,"vxMC/F");
  m_tree->Branch("vyMC",&vyMC,"vyMC/F");
  m_tree->Branch("vzMC",&vzMC,"vzMC/F");

  m_tree->Branch("vx",&vx,"vx/F");
  m_tree->Branch("vy",&vy,"vy/F");
  m_tree->Branch("vz",&vz,"vz/F");
  m_tree->Branch("vntracks",&vntracks,"vntracks/F");
  m_tree->Branch("vchi2",&vchi2,"vchi2/F");
  m_tree->Branch("vndof",&vndof,"vndof/F");

  // Set trigger bits of interest
  nHLT = 10;

  // photon triggers:
  hltTriggers["HLT_Photon10_L1R"] = 0;
  hltTriggers["HLT_Photon15_L1R"] = 1;
  hltTriggers["HLT_Photon25_L1R"] = 2;
  hltTriggers["HLT_Photon30_L1R_1E31"] = 3;
  hltTriggers["HLT_Photon70_L1R"] = 4;
  hltTriggers["HLT_DoublePhoton5_L1R"] = 5;
  hltTriggers["HLT_DoublePhoton10_L1R"] = 6;

  // electron triggers:
  hltTriggers["HLT_Ele15_LW_L1R"] = 7;
  hltTriggers["HLT_Ele15_LW_L1R'"] = 8;
  hltTriggers["HLT_Ele15_SW_EleId_L1R'"] = 9;
  hltTriggers["HLT_Ele15_SW_L1R'"] = 10;
  hltTriggers["HLT_Ele15_SiStrip_L1R'"] = 11;
  hltTriggers["HLT_Ele20_SW_L1R'"] = 12;
  hltTriggers["HLT_Ele20_SiStrip_L1R"] = 13;
  hltTriggers["HLT_DoubleEle5_SW_L1R"] = 14;

  // muon triggers:
  hltTriggers["HLT_L1DoubleMuOpen_Tight"] = 15;
  hltTriggers["HLT_Mu3"] = 16;
  hltTriggers["HLT_Mu5"] = 17;
  hltTriggers["HLT_Mu7"] = 18;
  hltTriggers["HLT_DoubleMu0"] = 19;
  hltTriggers["HLT_DoubleMu3"] = 20;


  m_tree->Branch("hltPass",&hltPass,"hltPass/O");
  //m_tree->Branch("hltCount",&hltCount,"hltCount/I");
  m_tree->Branch("nHLT",&nHLT,"nHLT/I");
  m_tree->Branch("hltNamesLen",&hltNamesLen,"hltNamesLen/I");
  m_tree->Branch("HLTNames",&aHLTNames,"HLTNames[hltNamesLen]/C,6000");
  //m_tree->Branch("HLTResults",&aHLTResults,"HLTResults[hltCount]/O");
  m_tree->Branch("HLTResults",&aHLTResults,"HLTResults[nHLT]/O");

  event = 0;  
}

// ------------ method called once each job just after ending the event loop  ------------
void HZZlljjAnalyzer::endJob() {

  
  outfile->cd("myanalysis");
  m_tree->Write();
  outfile->Close();
  //outfile->Delete();
  
//   //avoid writing the tree second time (automatically)
//  m_tree->Delete();

}


//define this as a plug-in
DEFINE_FWK_MODULE(HZZlljjAnalyzer);
