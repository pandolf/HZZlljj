// system include files
#include <memory>
#include <iostream>
#include <string>

// user include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

#include <map>
#include <set>
class SimTrack;
class SimVertex;

#define MAXHLTBITS    200

using namespace edm;

//
// class declaration
//

class HZZlljjAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HZZlljjAnalyzer(const edm::ParameterSet&);
      ~HZZlljjAnalyzer();


   private:
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      //  calculate phi1-phi2 keeping value between 0 and pi
      inline float delta_phi(float phi1, float phi2);
      // calculate eta1-eta2 keeping eta2 positive
      inline float delta_eta(float eta1, float eta2);
      // calculate sum in quadrature
      inline double oplus(double a, double b);

      // Constants
      static const int kParton = 3;
      static const int kPhoton = 22;
      static const int kElectron = 11;
      static const int kMuon = 13;


      // ----------member data ---------------------------
      bool _debug;

      edm::InputTag MCTruthCollection_; 
      edm::InputTag triggerTag_;
      edm::InputTag Vertexsrc_;
      edm::InputTag Electronsrc_; 
      edm::InputTag Muonsrc_; 
      std::string JetAlgosrc_; 
      edm::InputTag METsrc_; 
      edm::InputTag METGensrc_; 
      edm::InputTag trackTags_; 
      edm::InputTag HBhitsrc_; 
      std::string ecalHitsCollection_; 
      std::string ecalHitsProducer_; 
      std::string JetCorrector_; 
      double genjetptthr_;
      double jetptthr_;
      int jetnmin_;

//      edm::Service<TFileService> fs_;
      TFile* outfile;

      // Tree with multiple info
      TTree * m_tree ;

      // Auxiliary event info will help to study correction stability
      // for different stores, as a function of instantaneous lumi (pile-up),
      // bunch-crossing (out-of-time pile-up), orbit number (beam heating) etc.
      Bool_t isMC;
      Int_t store;
      Int_t lbn;
      Int_t bx;
      Int_t orbit;
      Int_t run;
      Int_t event;

      // Vertex distribution
      Int_t nvertex;
      Float_t vx[10];
      Float_t vy[10];
      Float_t vz[10];
      Float_t vntracks[10];
      Float_t vchi2[10];
      Float_t vndof[10];

 
      // Vertex distribution at MC truth level (only primary vertex)
      Float_t vxMC;
      Float_t vyMC;
      Float_t vzMC;

      // MC particles help to reconstruct original jet parton,
      // particle jet and jet properties (fake photon from pi0, rho0?)
      static const int nMaxMC = 150;//100;
      Int_t nMC;
      Int_t pdgIdMC[nMaxMC];
      Int_t statusMC[nMaxMC];
      //Float_t massMC[nMaxMC];
      Int_t motherIDMC[nMaxMC];
      Float_t ptMC[nMaxMC];
      Float_t eMC[nMaxMC];
      Float_t etaMC[nMaxMC];
      Float_t phiMC[nMaxMC];

      // SIM particles (those not already in MC particles list)
      // help to study in-flight decays of Kshort, Lambda etc.
      // These are also useful to study photon conversions and
      // performance of PFlow
      // NB: go back and mark decayed MC particles with status = -1 
      static const int nMaxSIM = 150;
      Int_t nSIM;
      Int_t pdgIdSIM[nMaxSIM];
      Int_t statusSIM[nMaxSIM];
      Float_t ptSIM[nMaxSIM];
      Float_t eSIM[nMaxSIM];
      Float_t etaSIM[nMaxSIM];
      Float_t phiSIM[nMaxSIM];
      Float_t rSIM[nMaxSIM];
      Float_t zSIM[nMaxSIM];

      // PFlow candindates from the two leading jets
      static const int nMaxPF = 1000;
      Int_t nPF;
      Int_t pdgIdPF[nMaxPF];
      Float_t ptPF[nMaxPF];
      Float_t ePF[nMaxPF];
      Float_t etaPF[nMaxPF];
      Float_t phiPF[nMaxPF];
      Float_t pTypePF[nMaxPF];

      Float_t genpt;

      static const int nMaxMu = 20;

      Int_t nMu;
      Float_t ptMu[nMaxMu];
      Float_t eMu[nMaxMu];
      Float_t etaMu[nMaxMu];
      Float_t phiMu[nMaxMu];
      
      Bool_t isGlobalMu[nMaxMu];
      Bool_t isGoodMu[nMaxMu];
      Bool_t isTrackerMu[nMaxMu];
      Int_t  nValidHitsTrackMu[nMaxMu];
      Int_t  nPixelLayerMu[nMaxMu];
      Int_t  nMatchesMu[nMaxMu];
      Float_t dxyMu[nMaxMu];

      static const int nMaxEle = 30;

      Int_t nEle;
      Float_t ptEle[nMaxEle];
      Float_t eEle[nMaxEle];
      Float_t etaEle[nMaxEle];
      Float_t phiEle[nMaxEle];
      
      Float_t trkIso03Ele[nMaxEle];
      Float_t ecalIso03Ele[nMaxEle];
      Float_t ecalIsoGT03Ele[nMaxEle];
      Float_t hcalIso03Ele[nMaxEle];

      Float_t hOverEmEle[nMaxEle];
      Float_t deltaPhiSCTrkEle[nMaxEle];
      Float_t deltaEtaSCTrkEle[nMaxEle];
      Float_t sigmaIetaIeta[nMaxEle];

      Int_t nInnerHitsEle[nMaxEle];
      Float_t convDistEle[nMaxEle];
      Float_t convDcotThetaEle[nMaxEle];

      static const int nMaxJet = 30;

      Int_t nJet;
      Float_t ptJet[nMaxJet];
      Float_t ptCorrJet[nMaxJet];
      Float_t eJet[nMaxJet];
      Float_t etaJet[nMaxJet];
      Float_t phiJet[nMaxJet];
      Float_t emfJet[nMaxJet];

      Int_t nChargedHadronsJet[nMaxJet];
      Int_t nPhotonsJet[nMaxJet];
      Int_t nElectronsJet[nMaxJet];
      Int_t nMuonsJet[nMaxJet];
      Int_t nNeutralHadronsJet[nMaxJet];
      Int_t nHFHadronsJet[nMaxJet];
      Int_t nHFEMJet[nMaxJet];

      Float_t eChargedHadronsJet[nMaxJet];
      Float_t ePhotonsJet[nMaxJet];
      Float_t eElectronsJet[nMaxJet];
      Float_t eMuonsJet[nMaxJet];
      Float_t eNeutralHadronsJet[nMaxJet];
      Float_t eHFHadronsJet[nMaxJet];
      Float_t eHFEMJet[nMaxJet];


      Int_t nJetGen;
      Float_t ptJetGen[nMaxJet];
      Float_t eJetGen[nMaxJet];
      Float_t etaJetGen[nMaxJet];
      Float_t phiJetGen[nMaxJet];


      Float_t spfMet;
      Float_t epfMet;
      Float_t phipfMet;

      Float_t sMetGen;
      Float_t eMetGen;
      Float_t phiMetGen;

      Float_t sMetGen2;
      Float_t eMetGen2;
      Float_t phiMetGen2;

      Bool_t   hltPass;
      char     aHLTNames[6000];
      Int_t    hltNamesLen;
      Int_t    hltCount;
      bool     aHLTResults[MAXHLTBITS];

      int nHLT;
      std::map<std::string, int> hltTriggers;
};

