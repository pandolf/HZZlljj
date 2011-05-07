//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 12 18:22:37 2009 by ROOT version 5.22/00
// from TChain pippo/
//////////////////////////////////////////////////////////


//------------------------------------------------------------------
//
//    Base Ntp1Analyzer class. Abstract.
//    Reads output of GammaJetAnalyzer, and produces subtrees. 
//    Loop() and BookStuff() methods have to be implemented.
//    Additional data members should be assigned.
//
//------------------------------------------------------------------



#ifndef Ntp1Analyzer_h
#define Ntp1Analyzer_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>


struct GenEventParameters{

  Float_t crossSection;
  Float_t ptHatMin;
  Float_t ptHatMax;

};


class Ntp1Analyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nl1Technical;
   Int_t           l1Technical[3];   //[nl1Technical]
   Int_t           nl1Global;
   Int_t           l1Global[5];   //[nl1Global]
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiBlock;
   Int_t           bunchCrossing;
   Int_t           orbitNumber;
   Int_t           nMc;
   Float_t         pMc[1000];   //[nMc]
   Float_t         thetaMc[1000];   //[nMc]
   Float_t         etaMc[1000];   //[nMc]
   Float_t         phiMc[1000];   //[nMc]
   Float_t         energyMc[1000];   //[nMc]
   Int_t           idMc[1000];   //[nMc]
   Int_t           mothMc[1000];   //[nMc]
   Int_t           statusMc[1000];   //[nMc]
   Int_t           nTrg;
   Int_t           firedTrg[4];   //[nTrg]
   Int_t           nHLT;
   Int_t           indexHLT[3000];   //[nHLT]
   std::vector<std::string>  *nameHLT;
   Int_t           nEle;
   Int_t           chargeEle[150];   //[nEle]
   Float_t         energyEle[150];   //[nEle]
   Float_t         thetaEle[150];   //[nEle]
   Float_t         etaEle[150];   //[nEle]
   Float_t         phiEle[150];   //[nEle]
   Float_t         pxEle[150];   //[nEle]
   Float_t         pyEle[150];   //[nEle]
   Float_t         pzEle[150];   //[nEle]
   Float_t         vertexXEle[150];   //[nEle]
   Float_t         vertexYEle[150];   //[nEle]
   Float_t         vertexZEle[150];   //[nEle]
   Int_t           fiducialFlagsEle[150];   //[nEle]
   Int_t           recoFlagsEle[150];   //[nEle]
   Int_t           energyCorrectionsEle[150];   //[nEle]
   Float_t         esEnergyEle[150];   //[nEle]
   Int_t           superClusterIndexEle[150];   //[nEle]
   Int_t           PFsuperClusterIndexEle[150];   //[nEle]
   Int_t           trackIndexEle[150];   //[nEle]
   Int_t           gsfTrackIndexEle[150];   //[nEle]
   Float_t         convDistEle[150];   //[nEle]
   Float_t         convDcotEle[150];   //[nEle]
   Float_t         convRadiusEle[150];   //[nEle]
   Int_t           convTrackIndexEle[150];   //[nEle]
   Float_t         convXEle[150];   //[nEle]
   Float_t         convYEle[150];   //[nEle]
   Float_t         convZEle[150];   //[nEle]
   Float_t         convChi2ProbEle[150];   //[nEle]
   Int_t           scPixChargeEle[150];   //[nEle]
   Int_t           classificationEle[150];   //[nEle]
   Int_t           standardClassificationEle[150];   //[nEle]
   Float_t         fbremEle[150];   //[nEle]
   Int_t           nbremsEle[150];   //[nEle]
   Float_t         hOverEEle[150];   //[nEle]
   Float_t         eSuperClusterOverPEle[150];   //[nEle]
   Float_t         eSeedOverPoutEle[150];   //[nEle]
   Float_t         deltaEtaAtVtxEle[150];   //[nEle]
   Float_t         deltaPhiAtVtxEle[150];   //[nEle]
   Float_t         deltaEtaAtCaloEle[150];   //[nEle]
   Float_t         deltaPhiAtCaloEle[150];   //[nEle]
   Float_t         tipEle[150];   //[nEle]
   Float_t         dr03TkSumPtEle[150];   //[nEle]
   Float_t         dr03EcalRecHitSumEtEle[150];   //[nEle]
   Float_t         dr03HcalTowerSumEtEle[150];   //[nEle]
   Float_t         dr04TkSumPtEle[150];   //[nEle]
   Float_t         dr04EcalRecHitSumEtEle[150];   //[nEle]
   Float_t         dr04HcalTowerSumEtEle[150];   //[nEle]
   Float_t         scBasedEcalSum03Ele[150];   //[nEle]
   Float_t         scBasedEcalSum04Ele[150];   //[nEle]
   Int_t           eleIdCutsEle[150];   //[nEle]
   Float_t         eleIdLikelihoodEle[150];   //[nEle]
   Float_t         pflowMVAEle[150];   //[nEle]
   Int_t           nSC;
   Int_t           nBCSC[200];   //[nSC]
   Int_t           nCrystalsSC[200];   //[nSC]
   Float_t         rawEnergySC[200];   //[nSC]
   Float_t         energySC[200];   //[nSC]
   Float_t         etaSC[200];   //[nSC]
   Float_t         thetaSC[200];   //[nSC]
   Float_t         phiSC[200];   //[nSC]
   Float_t         phiWidthSC[200];   //[nSC]
   Float_t         etaWidthSC[200];   //[nSC]
   Float_t         e3x3SC[200];   //[nSC]
   Float_t         e5x5SC[200];   //[nSC]
   Float_t         eMaxSC[200];   //[nSC]
   Float_t         e2x2SC[200];   //[nSC]
   Float_t         e2ndSC[200];   //[nSC]
   Float_t         e1x5SC[200];   //[nSC]
   Float_t         e2x5MaxSC[200];   //[nSC]
   Float_t         e4SwissCrossSC[200];   //[nSC]
   Float_t         covIEtaIEtaSC[200];   //[nSC]
   Float_t         covIEtaIPhiSC[200];   //[nSC]
   Float_t         covIPhiIPhiSC[200];   //[nSC]
   Float_t         hOverESC[200];   //[nSC]
   Int_t           recoFlagSC[200];   //[nSC]
   Int_t           channelStatusSC[200];   //[nSC]
   Float_t         timeSC[200];   //[nSC]
   Float_t         chi2SC[200];   //[nSC]
   Float_t         seedEnergySC[200];   //[nSC]
   Int_t           idClosProblSC[200];   //[nSC]
   Int_t           sevClosProblSC[200];   //[nSC]
   Float_t         fracClosProblSC[200];   //[nSC]
   Float_t         scBasedEcalSum03SC[200];   //[nSC]
   Float_t         scBasedEcalSum04SC[200];   //[nSC]
   Float_t         ecalRecHitSumEtConeDR03SC[200];   //[nSC]
   Float_t         hcalTowerSumEtConeDR03SC[200];   //[nSC]
   Float_t         trkSumPtSolidConeDR03SC[200];   //[nSC]
   Float_t         ecalRecHitSumEtConeDR04SC[200];   //[nSC]
   Float_t         hcalTowerSumEtConeDR04SC[200];   //[nSC]
   Float_t         trkSumPtSolidConeDR04SC[200];   //[nSC]
   Float_t         sMajSC[200];   //[nSC]
   Float_t         sMinSC[200];   //[nSC]
   Float_t         alphaSC[200];   //[nSC]
   Int_t           nPFSC;
   Int_t           nBCPFSC[200];   //[nPFSC]
   Int_t           nCrystalsPFSC[200];   //[nPFSC]
   Float_t         rawEnergyPFSC[200];   //[nPFSC]
   Float_t         energyPFSC[200];   //[nPFSC]
   Float_t         etaPFSC[200];   //[nPFSC]
   Float_t         thetaPFSC[200];   //[nPFSC]
   Float_t         phiPFSC[200];   //[nPFSC]
   Float_t         phiWidthPFSC[200];   //[nPFSC]
   Float_t         etaWidthPFSC[200];   //[nPFSC]
   Float_t         e3x3PFSC[200];   //[nPFSC]
   Float_t         e5x5PFSC[200];   //[nPFSC]
   Float_t         eMaxPFSC[200];   //[nPFSC]
   Float_t         e2x2PFSC[200];   //[nPFSC]
   Float_t         e2ndPFSC[200];   //[nPFSC]
   Float_t         e1x5PFSC[200];   //[nPFSC]
   Float_t         e2x5MaxPFSC[200];   //[nPFSC]
   Float_t         e4SwissCrossPFSC[200];   //[nPFSC]
   Float_t         covIEtaIEtaPFSC[200];   //[nPFSC]
   Float_t         covIEtaIPhiPFSC[200];   //[nPFSC]
   Float_t         covIPhiIPhiPFSC[200];   //[nPFSC]
   Float_t         hOverEPFSC[200];   //[nPFSC]
   Int_t           recoFlagPFSC[200];   //[nPFSC]
   Int_t           channelStatusPFSC[200];   //[nPFSC]
   Float_t         timePFSC[200];   //[nPFSC]
   Float_t         chi2PFSC[200];   //[nPFSC]
   Float_t         seedEnergyPFSC[200];   //[nPFSC]
   Int_t           idClosProblPFSC[200];   //[nPFSC]
   Int_t           sevClosProblPFSC[200];   //[nPFSC]
   Float_t         fracClosProblPFSC[200];   //[nPFSC]
   Float_t         scBasedEcalSum03PFSC[200];   //[nPFSC]
   Float_t         scBasedEcalSum04PFSC[200];   //[nPFSC]
   Int_t           nBC;
   Int_t           nCrystalsBC[200];   //[nBC]
   Float_t         energyBC[200];   //[nBC]
   Float_t         etaBC[200];   //[nBC]
   Float_t         thetaBC[200];   //[nBC]
   Float_t         phiBC[200];   //[nBC]
   Float_t         e3x3BC[200];   //[nBC]
   Float_t         e5x5BC[200];   //[nBC]
   Float_t         eMaxBC[200];   //[nBC]
   Float_t         e2x2BC[200];   //[nBC]
   Float_t         e2ndBC[200];   //[nBC]
   Float_t         covIEtaIEtaBC[200];   //[nBC]
   Float_t         covIEtaIPhiBC[200];   //[nBC]
   Float_t         covIPhiIPhiBC[200];   //[nBC]
   Int_t           recoFlagBC[200];   //[nBC]
   Float_t         timeBC[200];   //[nBC]
   Float_t         chi2BC[200];   //[nBC]
   Float_t         seedEnergyBC[200];   //[nBC]
   Int_t           idClosProblBC[200];   //[nBC]
   Int_t           sevClosProblBC[200];   //[nBC]
   Float_t         fracClosProblBC[200];   //[nBC]
   Int_t           indexSCBC[200];   //[nBC]
   Int_t           nTrack;
   Float_t         pxTrack[500];   //[nTrack]
   Float_t         pyTrack[500];   //[nTrack]
   Float_t         pzTrack[500];   //[nTrack]
   Int_t           vtxIndexTrack[500];   //[nTrack]
   Float_t         vtxWeightTrack[500];   //[nTrack]
   Float_t         chargeTrack[500];   //[nTrack]
   Float_t         ptErrorTrack[500];   //[nTrack]
   Float_t         trackValidHitsTrack[500];   //[nTrack]
   Float_t         trackLostHitsTrack[500];   //[nTrack]
   Float_t         trackNormalizedChi2Track[500];   //[nTrack]
   Int_t           qualityMaskTrack[500];   //[nTrack]
   Float_t         impactPar3DTrack[500];   //[nTrack]
   Float_t         impactPar3DErrorTrack[500];   //[nTrack]
   Float_t         transvImpactParTrack[500];   //[nTrack]
   Float_t         transvImpactParErrorTrack[500];   //[nTrack]
   Float_t         trackVxTrack[500];   //[nTrack]
   Float_t         trackVyTrack[500];   //[nTrack]
   Float_t         trackVzTrack[500];   //[nTrack]
   Float_t         pxAtOuterTrack[500];   //[nTrack]
   Float_t         pyAtOuterTrack[500];   //[nTrack]
   Float_t         pzAtOuterTrack[500];   //[nTrack]
   Float_t         xAtOuterTrack[500];   //[nTrack]
   Float_t         yAtOuterTrack[500];   //[nTrack]
   Float_t         zAtOuterTrack[500];   //[nTrack]
   Float_t         pxAtInnerTrack[500];   //[nTrack]
   Float_t         pyAtInnerTrack[500];   //[nTrack]
   Float_t         pzAtInnerTrack[500];   //[nTrack]
   Float_t         xAtInnerTrack[500];   //[nTrack]
   Float_t         yAtInnerTrack[500];   //[nTrack]
   Float_t         zAtInnerTrack[500];   //[nTrack]
   Float_t         recHitsSizeTrack[500];   //[nTrack]
   Int_t           pixelHitsTrack[500];   //[nTrack]
   Int_t           expInnerLayersTrack[500];   //[nTrack]
   Int_t           numberOfValidPixelBarrelHitsTrack[500];   //[nTrack]
   Int_t           numberOfValidPixelEndcapHitsTrack[500];   //[nTrack]
   Int_t           numberOfValidStripTIBHitsTrack[500];   //[nTrack]
   Int_t           numberOfValidStripTIDHitsTrack[500];   //[nTrack]
   Int_t           numberOfValidStripTOBHitsTrack[500];   //[nTrack]
   Int_t           numberOfValidStripTECHitsTrack[500];   //[nTrack]
   Int_t           nGsfTrack;
   Float_t         pxGsfTrack[150];   //[nGsfTrack]
   Float_t         pyGsfTrack[150];   //[nGsfTrack]
   Float_t         pzGsfTrack[150];   //[nGsfTrack]
   Int_t           vtxIndexGsfTrack[150];   //[nGsfTrack]
   Float_t         vtxWeightGsfTrack[150];   //[nGsfTrack]
   Float_t         chargeGsfTrack[150];   //[nGsfTrack]
   Float_t         ptErrorGsfTrack[150];   //[nGsfTrack]
   Float_t         trackValidHitsGsfTrack[150];   //[nGsfTrack]
   Float_t         trackLostHitsGsfTrack[150];   //[nGsfTrack]
   Float_t         trackNormalizedChi2GsfTrack[150];   //[nGsfTrack]
   Int_t           qualityMaskGsfTrack[150];   //[nGsfTrack]
   Float_t         impactPar3DGsfTrack[150];   //[nGsfTrack]
   Float_t         impactPar3DErrorGsfTrack[150];   //[nGsfTrack]
   Float_t         transvImpactParGsfTrack[150];   //[nGsfTrack]
   Float_t         transvImpactParErrorGsfTrack[150];   //[nGsfTrack]
   Float_t         trackVxGsfTrack[150];   //[nGsfTrack]
   Float_t         trackVyGsfTrack[150];   //[nGsfTrack]
   Float_t         trackVzGsfTrack[150];   //[nGsfTrack]
   Float_t         pxAtOuterGsfTrack[150];   //[nGsfTrack]
   Float_t         pyAtOuterGsfTrack[150];   //[nGsfTrack]
   Float_t         pzAtOuterGsfTrack[150];   //[nGsfTrack]
   Float_t         xAtOuterGsfTrack[150];   //[nGsfTrack]
   Float_t         yAtOuterGsfTrack[150];   //[nGsfTrack]
   Float_t         zAtOuterGsfTrack[150];   //[nGsfTrack]
   Float_t         pxAtInnerGsfTrack[150];   //[nGsfTrack]
   Float_t         pyAtInnerGsfTrack[150];   //[nGsfTrack]
   Float_t         pzAtInnerGsfTrack[150];   //[nGsfTrack]
   Float_t         xAtInnerGsfTrack[150];   //[nGsfTrack]
   Float_t         yAtInnerGsfTrack[150];   //[nGsfTrack]
   Float_t         zAtInnerGsfTrack[150];   //[nGsfTrack]
   Float_t         recHitsSizeGsfTrack[150];   //[nGsfTrack]
   Int_t           pixelHitsGsfTrack[150];   //[nGsfTrack]
   Int_t           expInnerLayersGsfTrack[150];   //[nGsfTrack]
   Int_t           numberOfValidPixelBarrelHitsGsfTrack[150];   //[nGsfTrack]
   Int_t           numberOfValidPixelEndcapHitsGsfTrack[150];   //[nGsfTrack]
   Int_t           numberOfValidStripTIBHitsGsfTrack[150];   //[nGsfTrack]
   Int_t           numberOfValidStripTIDHitsGsfTrack[150];   //[nGsfTrack]
   Int_t           numberOfValidStripTOBHitsGsfTrack[150];   //[nGsfTrack]
   Int_t           numberOfValidStripTECHitsGsfTrack[150];   //[nGsfTrack]
   Int_t           chargeModeGsfTrack[150];   //[nGsfTrack]
   Float_t         pxModeGsfTrack[150];   //[nGsfTrack]
   Float_t         pyModeGsfTrack[150];   //[nGsfTrack]
   Float_t         pzModeGsfTrack[150];   //[nGsfTrack]
   Int_t           recoFlagsGsfTrack[150];   //[nGsfTrack]
   Int_t           nGlobalMuonTrack;
   Float_t         pxGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         pyGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         pzGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           vtxIndexGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         vtxWeightGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         chargeGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         ptErrorGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         trackValidHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         trackLostHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         trackNormalizedChi2GlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           qualityMaskGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         impactPar3DGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         impactPar3DErrorGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         transvImpactParGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         transvImpactParErrorGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         trackVxGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         trackVyGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         trackVzGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         pxAtOuterGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         pyAtOuterGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         pzAtOuterGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         xAtOuterGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         yAtOuterGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         zAtOuterGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         pxAtInnerGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         pyAtInnerGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         pzAtInnerGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         xAtInnerGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         yAtInnerGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         zAtInnerGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Float_t         recHitsSizeGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           pixelHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           expInnerLayersGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           numberOfValidPixelBarrelHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           numberOfValidPixelEndcapHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTIBHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTIDHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTOBHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTECHitsGlobalMuonTrack[20];   //[nGlobalMuonTrack]
   Int_t           nSTAMuonTrack;
   Float_t         pxSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         pySTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         pzSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           vtxIndexSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         vtxWeightSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         chargeSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         ptErrorSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         trackValidHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         trackLostHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         trackNormalizedChi2STAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           qualityMaskSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         impactPar3DSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         impactPar3DErrorSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         transvImpactParSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         transvImpactParErrorSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         trackVxSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         trackVySTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         trackVzSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         pxAtOuterSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         pyAtOuterSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         pzAtOuterSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         xAtOuterSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         yAtOuterSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         zAtOuterSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         pxAtInnerSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         pyAtInnerSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         pzAtInnerSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         xAtInnerSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         yAtInnerSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         zAtInnerSTAMuonTrack[20];   //[nSTAMuonTrack]
   Float_t         recHitsSizeSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           pixelHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           expInnerLayersSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           numberOfValidPixelBarrelHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           numberOfValidPixelEndcapHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTIBHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTIDHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTOBHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTECHitsSTAMuonTrack[20];   //[nSTAMuonTrack]
   Int_t           nPV;
   Float_t         PVxPV[10];   //[nPV]
   Float_t         PVyPV[10];   //[nPV]
   Float_t         PVzPV[10];   //[nPV]
   Float_t         PVErrxPV[10];   //[nPV]
   Float_t         PVErryPV[10];   //[nPV]
   Float_t         PVErrzPV[10];   //[nPV]
   Float_t         SumPtPV[10];   //[nPV]
   Float_t         ndofPV[10];   //[nPV]
   Float_t         chi2PV[10];   //[nPV]
   Int_t           nMuon;
   Int_t           chargeMuon[150];   //[nMuon]
   Float_t         energyMuon[150];   //[nMuon]
   Float_t         thetaMuon[150];   //[nMuon]
   Float_t         etaMuon[150];   //[nMuon]
   Float_t         phiMuon[150];   //[nMuon]
   Float_t         pxMuon[150];   //[nMuon]
   Float_t         pyMuon[150];   //[nMuon]
   Float_t         pzMuon[150];   //[nMuon]
   Float_t         vertexXMuon[150];   //[nMuon]
   Float_t         vertexYMuon[150];   //[nMuon]
   Float_t         vertexZMuon[150];   //[nMuon]
   Int_t           trackIndexMuon[150];   //[nMuon]
   Int_t           standAloneTrackIndexMuon[150];   //[nMuon]
   Int_t           combinedTrackIndexMuon[150];   //[nMuon]
   Int_t           muonIdMuon[150];   //[nMuon]
   Float_t         sumPt03Muon[150];   //[nMuon]
   Float_t         emEt03Muon[150];   //[nMuon]
   Float_t         hadEt03Muon[150];   //[nMuon]
   Float_t         hoEt03Muon[150];   //[nMuon]
   Float_t         nTrk03Muon[150];   //[nMuon]
   Float_t         nJets03Muon[150];   //[nMuon]
   Float_t         sumPt05Muon[150];   //[nMuon]
   Float_t         emEt05Muon[150];   //[nMuon]
   Float_t         hadEt05Muon[150];   //[nMuon]
   Float_t         hoEt05Muon[150];   //[nMuon]
   Float_t         nTrk05Muon[150];   //[nMuon]
   Float_t         nJets05Muon[150];   //[nMuon]
   Float_t         EcalExpDepoMuon[150];   //[nMuon]
   Float_t         HcalExpDepoMuon[150];   //[nMuon]
   Float_t         HoExpDepoMuon[150];   //[nMuon]
   Float_t         emS9Muon[150];   //[nMuon]
   Float_t         hadS9Muon[150];   //[nMuon]
   Float_t         hoS9Muon[150];   //[nMuon]
   Float_t         CaloCompMuon[150];   //[nMuon]
   Int_t           nMet;
   Int_t           chargeMet[1];   //[nMet]
   Float_t         energyMet[1];   //[nMet]
   Float_t         thetaMet[1];   //[nMet]
   Float_t         etaMet[1];   //[nMet]
   Float_t         phiMet[1];   //[nMet]
   Float_t         pxMet[1];   //[nMet]
   Float_t         pyMet[1];   //[nMet]
   Float_t         pzMet[1];   //[nMet]
   Float_t         vertexXMet[1];   //[nMet]
   Float_t         vertexYMet[1];   //[nMet]
   Float_t         vertexZMet[1];   //[nMet]
   Int_t           nTCMet;
   Int_t           chargeTCMet[1];   //[nTCMet]
   Float_t         energyTCMet[1];   //[nTCMet]
   Float_t         thetaTCMet[1];   //[nTCMet]
   Float_t         etaTCMet[1];   //[nTCMet]
   Float_t         phiTCMet[1];   //[nTCMet]
   Float_t         pxTCMet[1];   //[nTCMet]
   Float_t         pyTCMet[1];   //[nTCMet]
   Float_t         pzTCMet[1];   //[nTCMet]
   Float_t         vertexXTCMet[1];   //[nTCMet]
   Float_t         vertexYTCMet[1];   //[nTCMet]
   Float_t         vertexZTCMet[1];   //[nTCMet]
   Int_t           nPFMet;
   Int_t           chargePFMet[1];   //[nPFMet]
   Float_t         energyPFMet[1];   //[nPFMet]
   Float_t         thetaPFMet[1];   //[nPFMet]
   Float_t         etaPFMet[1];   //[nPFMet]
   Float_t         phiPFMet[1];   //[nPFMet]
   Float_t         pxPFMet[1];   //[nPFMet]
   Float_t         pyPFMet[1];   //[nPFMet]
   Float_t         pzPFMet[1];   //[nPFMet]
   Float_t         vertexXPFMet[1];   //[nPFMet]
   Float_t         vertexYPFMet[1];   //[nPFMet]
   Float_t         vertexZPFMet[1];   //[nPFMet]
   Int_t           nGenMet;
   Int_t           chargeGenMet[1];   //[nGenMet]
   Float_t         energyGenMet[1];   //[nGenMet]
   Float_t         thetaGenMet[1];   //[nGenMet]
   Float_t         etaGenMet[1];   //[nGenMet]
   Float_t         phiGenMet[1];   //[nGenMet]
   Float_t         pxGenMet[1];   //[nGenMet]
   Float_t         pyGenMet[1];   //[nGenMet]
   Float_t         pzGenMet[1];   //[nGenMet]
   Float_t         vertexXGenMet[1];   //[nGenMet]
   Float_t         vertexYGenMet[1];   //[nGenMet]
   Float_t         vertexZGenMet[1];   //[nGenMet]
   Int_t           nPFCand;
   Int_t           chargePFCand[600];   //[nPFCand]
   Float_t         energyPFCand[600];   //[nPFCand]
   Float_t         thetaPFCand[600];   //[nPFCand]
   Float_t         etaPFCand[600];   //[nPFCand]
   Float_t         phiPFCand[600];   //[nPFCand]
   Float_t         pxPFCand[600];   //[nPFCand]
   Float_t         pyPFCand[600];   //[nPFCand]
   Float_t         pzPFCand[600];   //[nPFCand]
   Float_t         vertexXPFCand[600];   //[nPFCand]
   Float_t         vertexYPFCand[600];   //[nPFCand]
   Float_t         vertexZPFCand[600];   //[nPFCand]
   Int_t           particleTypePFCand[600];   //[nPFCand]
   Int_t           iPFJetPFCand[600];   //[nPFCand]
   Int_t           nAK5Jet;
   Int_t           chargeAK5Jet[300];   //[nAK5Jet]
   Float_t         energyAK5Jet[300];   //[nAK5Jet]
   Float_t         thetaAK5Jet[300];   //[nAK5Jet]
   Float_t         etaAK5Jet[300];   //[nAK5Jet]
   Float_t         phiAK5Jet[300];   //[nAK5Jet]
   Float_t         pxAK5Jet[300];   //[nAK5Jet]
   Float_t         pyAK5Jet[300];   //[nAK5Jet]
   Float_t         pzAK5Jet[300];   //[nAK5Jet]
   Float_t         vertexXAK5Jet[300];   //[nAK5Jet]
   Float_t         vertexYAK5Jet[300];   //[nAK5Jet]
   Float_t         vertexZAK5Jet[300];   //[nAK5Jet]
   Float_t         emFracAK5Jet[300];   //[nAK5Jet]
   Float_t         hadFracAK5Jet[300];   //[nAK5Jet]
   Float_t         combinedSecondaryVertexBJetTagsAK5Jet[300];   //[nAK5Jet]
   Float_t         combinedSecondaryVertexMVABJetTagsAK5Jet[300];   //[nAK5Jet]
   Float_t         jetBProbabilityBJetTagsAK5Jet[300];   //[nAK5Jet]
   Float_t         jetProbabilityBJetTagsAK5Jet[300];   //[nAK5Jet]
   Float_t         simpleSecondaryVertexBJetTagsAK5Jet[300];   //[nAK5Jet]
   Float_t         softMuonBJetTagsAK5Jet[300];   //[nAK5Jet]
   Float_t         trackCountingHighPurBJetTagsAK5Jet[300];   //[nAK5Jet]
   Float_t         trackCountingHighEffBJetTagsAK5Jet[300];   //[nAK5Jet]
   Float_t         uncorrEnergyAK5Jet[300];   //[nAK5Jet]
   Int_t           nAK5PFPUcorrJet;
   Int_t           chargeAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         energyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         thetaAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         etaAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         phiAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         pxAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         pyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         pzAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         vertexXAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         vertexYAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         vertexZAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         chargedHadronEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         neutralHadronEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         photonEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         electronEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         muonEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         HFHadronEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         HFEMEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Int_t           chargedHadronMultiplicityAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Int_t           neutralHadronMultiplicityAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Int_t           photonMultiplicityAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Int_t           electronMultiplicityAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Int_t           muonMultiplicityAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Int_t           HFHadronMultiplicityAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Int_t           HFEMMultiplicityAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         chargedEmEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         neutralEmEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         jetBProbabilityBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         jetProbabilityBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         softMuonBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         softMuonByIP3dBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         softMuonByPtBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         softElectronBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         softElectronByIP3dBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         softElectronByPtBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         trackCountingHighPurBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         trackCountingHighEffBJetTagsAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         uncorrEnergyAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         ptDAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
   Float_t         rmsCandAK5PFPUcorrJet[300];   //[nAK5PFPUcorrJet]
 //Int_t           nAK5JPTJet;
 //Int_t           chargeAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         energyAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         thetaAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         etaAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         phiAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         pxAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         pyAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         pzAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         vertexXAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         vertexYAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         vertexZAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         emFracAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         hadFracAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         combinedSecondaryVertexBJetTagsAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         combinedSecondaryVertexMVABJetTagsAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         jetBProbabilityBJetTagsAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         jetProbabilityBJetTagsAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         simpleSecondaryVertexBJetTagsAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         softMuonBJetTagsAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         trackCountingHighPurBJetTagsAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         trackCountingHighEffBJetTagsAK5JPTJet[300];   //[nAK5JPTJet]
 //Float_t         uncorrEnergyAK5JPTJet[300];   //[nAK5JPTJet]
   Int_t           nAK5GenJet;
   Int_t           chargeAK5GenJet[300];   //[nAK5GenJet]
   Float_t         energyAK5GenJet[300];   //[nAK5GenJet]
   Float_t         thetaAK5GenJet[300];   //[nAK5GenJet]
   Float_t         etaAK5GenJet[300];   //[nAK5GenJet]
   Float_t         phiAK5GenJet[300];   //[nAK5GenJet]
   Float_t         pxAK5GenJet[300];   //[nAK5GenJet]
   Float_t         pyAK5GenJet[300];   //[nAK5GenJet]
   Float_t         pzAK5GenJet[300];   //[nAK5GenJet]
   Float_t         vertexXAK5GenJet[300];   //[nAK5GenJet]
   Float_t         vertexYAK5GenJet[300];   //[nAK5GenJet]
   Float_t         vertexZAK5GenJet[300];   //[nAK5GenJet]
   Double_t        genPtHat;
   Double_t        genProcessId;
   Double_t        genWeight;
   Double_t        genAlphaQCD;
   Double_t        genAlphaQED;


   std::string dataset_;
   bool isMC_;

   std::string recoType_;
   std::string jetAlgo_;
   std::string algoType_;

   std::string analyzerType_;
   std::string flags_;

   typedef std::vector< std::pair<int,int> > GoodLSVector;
   typedef std::map< int, GoodLSVector  >    LSRange ;
   LSRange goodLS_;
   LSRange::const_iterator goodLSCache_; // ptr to list of good LS for last run
   bool filterGoodRuns_;

   std::vector<std::string> requiredTriggers_;
   std::vector<int> index_requiredTriggers_;

   typedef std::pair< int, int > RunLumiPair;
   typedef std::map< RunLumiPair, double > LSLumi;
   LSLumi LSLumimap_;
   double totalIntLumi_;

   TFile* outfile_;
   TTree* reducedTree_;
   // this histogram saves total number of analyzed events in file:
   // (needed for event weight determination in last step)
   TH1F* h1_nCounter_;

   Int_t run_;
   Int_t oldrun_;
   Int_t cachedLS_;
   Int_t cachedRun_;
   Int_t LS_;
   Int_t nvertex_;
   Int_t event_;
   Float_t ptHat_;
   Float_t ptHatMin_;
   Float_t ptHatMax_;
   Float_t eventWeight_;
  

   bool DEBUG_;


   // List of branches
   TBranch        *b_nl1Technical;   //!
   TBranch        *b_l1Technical;   //!
   TBranch        *b_nl1Global;   //!
   TBranch        *b_l1Global;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
   TBranch        *b_nMc;   //!
   TBranch        *b_pMc;   //!
   TBranch        *b_thetaMc;   //!
   TBranch        *b_etaMc;   //!
   TBranch        *b_phiMc;   //!
   TBranch        *b_energyMc;   //!
   TBranch        *b_idMc;   //!
   TBranch        *b_mothMc;   //!
   TBranch        *b_statusMc;   //!
   TBranch        *b_nTrg;   //!
   TBranch        *b_firedTrg;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_nameHLT;   //!
   TBranch        *b_indexHLT;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_energyEle;   //!
   TBranch        *b_thetaEle;   //!
   TBranch        *b_etaEle;   //!
   TBranch        *b_phiEle;   //!
   TBranch        *b_pxEle;   //!
   TBranch        *b_pyEle;   //!
   TBranch        *b_pzEle;   //!
   TBranch        *b_vertexXEle;   //!
   TBranch        *b_vertexYEle;   //!
   TBranch        *b_vertexZEle;   //!
   TBranch        *b_fiducialFlagsEle;   //!
   TBranch        *b_recoFlagsEle;   //!
   TBranch        *b_energyCorrectionsEle;   //!
   TBranch        *b_esEnergyEle;   //!
   TBranch        *b_superClusterIndexEle;   //!
   TBranch        *b_PFsuperClusterIndexEle;   //!
   TBranch        *b_trackIndexEle;   //!
   TBranch        *b_gsfTrackIndexEle;   //!
   TBranch        *b_convDistEle;   //!
   TBranch        *b_convDcotEle;   //!
   TBranch        *b_convRadiusEle;   //!
   TBranch        *b_convTrackIndexEle;   //!
   TBranch        *b_convXEle;   //!
   TBranch        *b_convYEle;   //!
   TBranch        *b_convZEle;   //!
   TBranch        *b_convChi2ProbEle;   //!
   TBranch        *b_scPixChargeEle;   //!
   TBranch        *b_classificationEle;   //!
   TBranch        *b_standardClassificationEle;   //!
   TBranch        *b_fbremEle;   //!
   TBranch        *b_nbremsEle;   //!
   TBranch        *b_hOverEEle;   //!
   TBranch        *b_eSuperClusterOverPEle;   //!
   TBranch        *b_eSeedOverPoutEle;   //!
   TBranch        *b_deltaEtaAtVtxEle;   //!
   TBranch        *b_deltaPhiAtVtxEle;   //!
   TBranch        *b_deltaEtaAtCaloEle;   //!
   TBranch        *b_deltaPhiAtCaloEle;   //!
   TBranch        *b_tipEle;   //!
   TBranch        *b_dr03TkSumPtEle;   //!
   TBranch        *b_dr03EcalRecHitSumEtEle;   //!
   TBranch        *b_dr03HcalTowerSumEtEle;   //!
   TBranch        *b_dr04TkSumPtEle;   //!
   TBranch        *b_dr04EcalRecHitSumEtEle;   //!
   TBranch        *b_dr04HcalTowerSumEtEle;   //!
   TBranch        *b_scBasedEcalSum03Ele;   //!
   TBranch        *b_scBasedEcalSum04Ele;   //!
   TBranch        *b_eleIdCutsEle;   //!
   TBranch        *b_eleIdLikelihoodEle;   //!
   TBranch        *b_pflowMVAEle;   //!
   TBranch        *b_nPFEle;   //!
   TBranch        *b_chargePFEle;   //!
   TBranch        *b_energyPFEle;   //!
   TBranch        *b_thetaPFEle;   //!
   TBranch        *b_etaPFEle;   //!
   TBranch        *b_phiPFEle;   //!
   TBranch        *b_pxPFEle;   //!
   TBranch        *b_pyPFEle;   //!
   TBranch        *b_pzPFEle;   //!
   TBranch        *b_vertexXPFEle;   //!
   TBranch        *b_vertexYPFEle;   //!
   TBranch        *b_vertexZPFEle;   //!
   TBranch        *b_MvaOutputPFEle;   //!
   TBranch        *b_PS1EnergyPFEle;   //!
   TBranch        *b_PS2EnergyPFEle;   //!
   TBranch        *b_EcalEnergyPFEle;   //!
   TBranch        *b_HcalEnergyPFEle;   //!
   TBranch        *b_RawEcalEnergyPFEle;   //!
   TBranch        *b_RawHcalEnergyPFEle;   //!
   TBranch        *b_PositionAtEcalXPFEle;   //!
   TBranch        *b_PositionAtEcalYPFEle;   //!
   TBranch        *b_PositionAtEcalZPFEle;   //!
   TBranch        *b_gsfTrackIndexPFEle;   //!
   TBranch        *b_trackIndexPFEle;   //!
   TBranch        *b_chIso03vetoPFEle;   //!
   TBranch        *b_chIso04vetoPFEle;   //!
   TBranch        *b_chIso05vetoPFEle;   //!
   TBranch        *b_chIso03noVetoPFEle;   //!
   TBranch        *b_chIso04noVetoPFEle;   //!
   TBranch        *b_chIso05noVetoPFEle;   //!
   TBranch        *b_nhIso03vetoPFEle;   //!
   TBranch        *b_nhIso04vetoPFEle;   //!
   TBranch        *b_nhIso05vetoPFEle;   //!
   TBranch        *b_nhIso03noVetoPFEle;   //!
   TBranch        *b_nhIso04noVetoPFEle;   //!
   TBranch        *b_nhIso05noVetoPFEle;   //!
   TBranch        *b_phIso03vetoPFEle;   //!
   TBranch        *b_phIso04vetoPFEle;   //!
   TBranch        *b_phIso05vetoPFEle;   //!
   TBranch        *b_phIso03noVetoPFEle;   //!
   TBranch        *b_phIso04noVetoPFEle;   //!
   TBranch        *b_phIso05noVetoPFEle;   //!
   TBranch        *b_nSC;   //!
   TBranch        *b_nBCSC;   //!
   TBranch        *b_nCrystalsSC;   //!
   TBranch        *b_rawEnergySC;   //!
   TBranch        *b_energySC;   //!
   TBranch        *b_etaSC;   //!
   TBranch        *b_thetaSC;   //!
   TBranch        *b_phiSC;   //!
   TBranch        *b_phiWidthSC;   //!
   TBranch        *b_etaWidthSC;   //!
   TBranch        *b_e3x3SC;   //!
   TBranch        *b_e5x5SC;   //!
   TBranch        *b_eMaxSC;   //!
   TBranch        *b_e2x2SC;   //!
   TBranch        *b_e2ndSC;   //!
   TBranch        *b_e1x5SC;   //!
   TBranch        *b_e2x5MaxSC;   //!
   TBranch        *b_e4SwissCrossSC;   //!
   TBranch        *b_covIEtaIEtaSC;   //!
   TBranch        *b_covIEtaIPhiSC;   //!
   TBranch        *b_covIPhiIPhiSC;   //!
   TBranch        *b_hOverESC;   //!
   TBranch        *b_recoFlagSC;   //!
   TBranch        *b_channelStatusSC;   //!
   TBranch        *b_timeSC;   //!
   TBranch        *b_chi2SC;   //!
   TBranch        *b_seedEnergySC;   //!
   TBranch        *b_idClosProblSC;   //!
   TBranch        *b_sevClosProblSC;   //!
   TBranch        *b_fracClosProblSC;   //!
   TBranch        *b_scBasedEcalSum03SC;   //!
   TBranch        *b_scBasedEcalSum04SC;   //!
   TBranch        *b_ecalRecHitSumEtConeDR03SC;   //!
   TBranch        *b_hcalTowerSumEtConeDR03SC;   //!
   TBranch        *b_trkSumPtSolidConeDR03SC;   //!
   TBranch        *b_ecalRecHitSumEtConeDR04SC;   //!
   TBranch        *b_hcalTowerSumEtConeDR04SC;   //!
   TBranch        *b_trkSumPtSolidConeDR04SC;   //!
   TBranch        *b_sMajSC;   //!
   TBranch        *b_sMinSC;   //!
   TBranch        *b_alphaSC;   //!
   TBranch        *b_nPFSC;   //!
   TBranch        *b_nBCPFSC;   //!
   TBranch        *b_nCrystalsPFSC;   //!
   TBranch        *b_rawEnergyPFSC;   //!
   TBranch        *b_energyPFSC;   //!
   TBranch        *b_etaPFSC;   //!
   TBranch        *b_thetaPFSC;   //!
   TBranch        *b_phiPFSC;   //!
   TBranch        *b_phiWidthPFSC;   //!
   TBranch        *b_etaWidthPFSC;   //!
   TBranch        *b_e3x3PFSC;   //!
   TBranch        *b_e5x5PFSC;   //!
   TBranch        *b_eMaxPFSC;   //!
   TBranch        *b_e2x2PFSC;   //!
   TBranch        *b_e2ndPFSC;   //!
   TBranch        *b_e1x5PFSC;   //!
   TBranch        *b_e2x5MaxPFSC;   //!
   TBranch        *b_e4SwissCrossPFSC;   //!
   TBranch        *b_covIEtaIEtaPFSC;   //!
   TBranch        *b_covIEtaIPhiPFSC;   //!
   TBranch        *b_covIPhiIPhiPFSC;   //!
   TBranch        *b_hOverEPFSC;   //!
   TBranch        *b_recoFlagPFSC;   //!
   TBranch        *b_channelStatusPFSC;   //!
   TBranch        *b_timePFSC;   //!
   TBranch        *b_chi2PFSC;   //!
   TBranch        *b_seedEnergyPFSC;   //!
   TBranch        *b_idClosProblPFSC;   //!
   TBranch        *b_sevClosProblPFSC;   //!
   TBranch        *b_fracClosProblPFSC;   //!
   TBranch        *b_scBasedEcalSum03PFSC;   //!
   TBranch        *b_scBasedEcalSum04PFSC;   //!
   TBranch        *b_nBC;   //!
   TBranch        *b_nCrystalsBC;   //!
   TBranch        *b_energyBC;   //!
   TBranch        *b_etaBC;   //!
   TBranch        *b_thetaBC;   //!
   TBranch        *b_phiBC;   //!
   TBranch        *b_e3x3BC;   //!
   TBranch        *b_e5x5BC;   //!
   TBranch        *b_eMaxBC;   //!
   TBranch        *b_e2x2BC;   //!
   TBranch        *b_e2ndBC;   //!
   TBranch        *b_covIEtaIEtaBC;   //!
   TBranch        *b_covIEtaIPhiBC;   //!
   TBranch        *b_covIPhiIPhiBC;   //!
   TBranch        *b_recoFlagBC;   //!
   TBranch        *b_timeBC;   //!
   TBranch        *b_chi2BC;   //!
   TBranch        *b_seedEnergyBC;   //!
   TBranch        *b_idClosProblBC;   //!
   TBranch        *b_sevClosProblBC;   //!
   TBranch        *b_fracClosProblBC;   //!
   TBranch        *b_indexSCBC;   //!
   TBranch        *b_nTrack;   //!
   TBranch        *b_pxTrack;   //!
   TBranch        *b_pyTrack;   //!
   TBranch        *b_pzTrack;   //!
   TBranch        *b_vtxIndexTrack;   //!
   TBranch        *b_vtxWeightTrack;   //!
   TBranch        *b_chargeTrack;   //!
   TBranch        *b_ptErrorTrack;   //!
   TBranch        *b_trackValidHitsTrack;   //!
   TBranch        *b_trackLostHitsTrack;   //!
   TBranch        *b_trackNormalizedChi2Track;   //!
   TBranch        *b_qualityMaskTrack;   //!
   TBranch        *b_impactPar3DTrack;   //!
   TBranch        *b_impactPar3DErrorTrack;   //!
   TBranch        *b_transvImpactParTrack;   //!
   TBranch        *b_transvImpactParErrorTrack;   //!
   TBranch        *b_trackVxTrack;   //!
   TBranch        *b_trackVyTrack;   //!
   TBranch        *b_trackVzTrack;   //!
   TBranch        *b_pxAtOuterTrack;   //!
   TBranch        *b_pyAtOuterTrack;   //!
   TBranch        *b_pzAtOuterTrack;   //!
   TBranch        *b_xAtOuterTrack;   //!
   TBranch        *b_yAtOuterTrack;   //!
   TBranch        *b_zAtOuterTrack;   //!
   TBranch        *b_pxAtInnerTrack;   //!
   TBranch        *b_pyAtInnerTrack;   //!
   TBranch        *b_pzAtInnerTrack;   //!
   TBranch        *b_xAtInnerTrack;   //!
   TBranch        *b_yAtInnerTrack;   //!
   TBranch        *b_zAtInnerTrack;   //!
   TBranch        *b_recHitsSizeTrack;   //!
   TBranch        *b_pixelHitsTrack;   //!
   TBranch        *b_expInnerLayersTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsTrack;   //!
   TBranch        *b_nGsfTrack;   //!
   TBranch        *b_pxGsfTrack;   //!
   TBranch        *b_pyGsfTrack;   //!
   TBranch        *b_pzGsfTrack;   //!
   TBranch        *b_vtxIndexGsfTrack;   //!
   TBranch        *b_vtxWeightGsfTrack;   //!
   TBranch        *b_chargeGsfTrack;   //!
   TBranch        *b_ptErrorGsfTrack;   //!
   TBranch        *b_trackValidHitsGsfTrack;   //!
   TBranch        *b_trackLostHitsGsfTrack;   //!
   TBranch        *b_trackNormalizedChi2GsfTrack;   //!
   TBranch        *b_qualityMaskGsfTrack;   //!
   TBranch        *b_impactPar3DGsfTrack;   //!
   TBranch        *b_impactPar3DErrorGsfTrack;   //!
   TBranch        *b_transvImpactParGsfTrack;   //!
   TBranch        *b_transvImpactParErrorGsfTrack;   //!
   TBranch        *b_trackVxGsfTrack;   //!
   TBranch        *b_trackVyGsfTrack;   //!
   TBranch        *b_trackVzGsfTrack;   //!
   TBranch        *b_pxAtOuterGsfTrack;   //!
   TBranch        *b_pyAtOuterGsfTrack;   //!
   TBranch        *b_pzAtOuterGsfTrack;   //!
   TBranch        *b_xAtOuterGsfTrack;   //!
   TBranch        *b_yAtOuterGsfTrack;   //!
   TBranch        *b_zAtOuterGsfTrack;   //!
   TBranch        *b_pxAtInnerGsfTrack;   //!
   TBranch        *b_pyAtInnerGsfTrack;   //!
   TBranch        *b_pzAtInnerGsfTrack;   //!
   TBranch        *b_xAtInnerGsfTrack;   //!
   TBranch        *b_yAtInnerGsfTrack;   //!
   TBranch        *b_zAtInnerGsfTrack;   //!
   TBranch        *b_recHitsSizeGsfTrack;   //!
   TBranch        *b_pixelHitsGsfTrack;   //!
   TBranch        *b_expInnerLayersGsfTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsGsfTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsGsfTrack;   //!
   TBranch        *b_chargeModeGsfTrack;   //!
   TBranch        *b_pxModeGsfTrack;   //!
   TBranch        *b_pyModeGsfTrack;   //!
   TBranch        *b_pzModeGsfTrack;   //!
   TBranch        *b_recoFlagsGsfTrack;   //!
   TBranch        *b_nGlobalMuonTrack;   //!
   TBranch        *b_pxGlobalMuonTrack;   //!
   TBranch        *b_pyGlobalMuonTrack;   //!
   TBranch        *b_pzGlobalMuonTrack;   //!
   TBranch        *b_vtxIndexGlobalMuonTrack;   //!
   TBranch        *b_vtxWeightGlobalMuonTrack;   //!
   TBranch        *b_chargeGlobalMuonTrack;   //!
   TBranch        *b_ptErrorGlobalMuonTrack;   //!
   TBranch        *b_trackValidHitsGlobalMuonTrack;   //!
   TBranch        *b_trackLostHitsGlobalMuonTrack;   //!
   TBranch        *b_trackNormalizedChi2GlobalMuonTrack;   //!
   TBranch        *b_qualityMaskGlobalMuonTrack;   //!
   TBranch        *b_impactPar3DGlobalMuonTrack;   //!
   TBranch        *b_impactPar3DErrorGlobalMuonTrack;   //!
   TBranch        *b_transvImpactParGlobalMuonTrack;   //!
   TBranch        *b_transvImpactParErrorGlobalMuonTrack;   //!
   TBranch        *b_trackVxGlobalMuonTrack;   //!
   TBranch        *b_trackVyGlobalMuonTrack;   //!
   TBranch        *b_trackVzGlobalMuonTrack;   //!
   TBranch        *b_pxAtOuterGlobalMuonTrack;   //!
   TBranch        *b_pyAtOuterGlobalMuonTrack;   //!
   TBranch        *b_pzAtOuterGlobalMuonTrack;   //!
   TBranch        *b_xAtOuterGlobalMuonTrack;   //!
   TBranch        *b_yAtOuterGlobalMuonTrack;   //!
   TBranch        *b_zAtOuterGlobalMuonTrack;   //!
   TBranch        *b_pxAtInnerGlobalMuonTrack;   //!
   TBranch        *b_pyAtInnerGlobalMuonTrack;   //!
   TBranch        *b_pzAtInnerGlobalMuonTrack;   //!
   TBranch        *b_xAtInnerGlobalMuonTrack;   //!
   TBranch        *b_yAtInnerGlobalMuonTrack;   //!
   TBranch        *b_zAtInnerGlobalMuonTrack;   //!
   TBranch        *b_recHitsSizeGlobalMuonTrack;   //!
   TBranch        *b_pixelHitsGlobalMuonTrack;   //!
   TBranch        *b_expInnerLayersGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsGlobalMuonTrack;   //!
   TBranch        *b_nSTAMuonTrack;   //!
   TBranch        *b_pxSTAMuonTrack;   //!
   TBranch        *b_pySTAMuonTrack;   //!
   TBranch        *b_pzSTAMuonTrack;   //!
   TBranch        *b_vtxIndexSTAMuonTrack;   //!
   TBranch        *b_vtxWeightSTAMuonTrack;   //!
   TBranch        *b_chargeSTAMuonTrack;   //!
   TBranch        *b_ptErrorSTAMuonTrack;   //!
   TBranch        *b_trackValidHitsSTAMuonTrack;   //!
   TBranch        *b_trackLostHitsSTAMuonTrack;   //!
   TBranch        *b_trackNormalizedChi2STAMuonTrack;   //!
   TBranch        *b_qualityMaskSTAMuonTrack;   //!
   TBranch        *b_impactPar3DSTAMuonTrack;   //!
   TBranch        *b_impactPar3DErrorSTAMuonTrack;   //!
   TBranch        *b_transvImpactParSTAMuonTrack;   //!
   TBranch        *b_transvImpactParErrorSTAMuonTrack;   //!
   TBranch        *b_trackVxSTAMuonTrack;   //!
   TBranch        *b_trackVySTAMuonTrack;   //!
   TBranch        *b_trackVzSTAMuonTrack;   //!
   TBranch        *b_pxAtOuterSTAMuonTrack;   //!
   TBranch        *b_pyAtOuterSTAMuonTrack;   //!
   TBranch        *b_pzAtOuterSTAMuonTrack;   //!
   TBranch        *b_xAtOuterSTAMuonTrack;   //!
   TBranch        *b_yAtOuterSTAMuonTrack;   //!
   TBranch        *b_zAtOuterSTAMuonTrack;   //!
   TBranch        *b_pxAtInnerSTAMuonTrack;   //!
   TBranch        *b_pyAtInnerSTAMuonTrack;   //!
   TBranch        *b_pzAtInnerSTAMuonTrack;   //!
   TBranch        *b_xAtInnerSTAMuonTrack;   //!
   TBranch        *b_yAtInnerSTAMuonTrack;   //!
   TBranch        *b_zAtInnerSTAMuonTrack;   //!
   TBranch        *b_recHitsSizeSTAMuonTrack;   //!
   TBranch        *b_pixelHitsSTAMuonTrack;   //!
   TBranch        *b_expInnerLayersSTAMuonTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsSTAMuonTrack;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVxPV;   //!
   TBranch        *b_PVyPV;   //!
   TBranch        *b_PVzPV;   //!
   TBranch        *b_PVErrxPV;   //!
   TBranch        *b_PVErryPV;   //!
   TBranch        *b_PVErrzPV;   //!
   TBranch        *b_SumPtPV;   //!
   TBranch        *b_ndofPV;   //!
   TBranch        *b_chi2PV;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_chargeMuon;   //!
   TBranch        *b_energyMuon;   //!
   TBranch        *b_thetaMuon;   //!
   TBranch        *b_etaMuon;   //!
   TBranch        *b_phiMuon;   //!
   TBranch        *b_pxMuon;   //!
   TBranch        *b_pyMuon;   //!
   TBranch        *b_pzMuon;   //!
   TBranch        *b_vertexXMuon;   //!
   TBranch        *b_vertexYMuon;   //!
   TBranch        *b_vertexZMuon;   //!
   TBranch        *b_trackIndexMuon;   //!
   TBranch        *b_standAloneTrackIndexMuon;   //!
   TBranch        *b_combinedTrackIndexMuon;   //!
   TBranch        *b_muonIdMuon;   //!
   TBranch        *b_sumPt03Muon;   //!
   TBranch        *b_emEt03Muon;   //!
   TBranch        *b_hadEt03Muon;   //!
   TBranch        *b_hoEt03Muon;   //!
   TBranch        *b_nTrk03Muon;   //!
   TBranch        *b_nJets03Muon;   //!
   TBranch        *b_sumPt05Muon;   //!
   TBranch        *b_emEt05Muon;   //!
   TBranch        *b_hadEt05Muon;   //!
   TBranch        *b_hoEt05Muon;   //!
   TBranch        *b_nTrk05Muon;   //!
   TBranch        *b_nJets05Muon;   //!
   TBranch        *b_EcalExpDepoMuon;   //!
   TBranch        *b_HcalExpDepoMuon;   //!
   TBranch        *b_HoExpDepoMuon;   //!
   TBranch        *b_emS9Muon;   //!
   TBranch        *b_hadS9Muon;   //!
   TBranch        *b_hoS9Muon;   //!
   TBranch        *b_CaloCompMuon;   //!
   TBranch        *b_nMet;   //!
   TBranch        *b_chargeMet;   //!
   TBranch        *b_energyMet;   //!
   TBranch        *b_thetaMet;   //!
   TBranch        *b_etaMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_pxMet;   //!
   TBranch        *b_pyMet;   //!
   TBranch        *b_pzMet;   //!
   TBranch        *b_vertexXMet;   //!
   TBranch        *b_vertexYMet;   //!
   TBranch        *b_vertexZMet;   //!
   TBranch        *b_nTCMet;   //!
   TBranch        *b_chargeTCMet;   //!
   TBranch        *b_energyTCMet;   //!
   TBranch        *b_thetaTCMet;   //!
   TBranch        *b_etaTCMet;   //!
   TBranch        *b_phiTCMet;   //!
   TBranch        *b_pxTCMet;   //!
   TBranch        *b_pyTCMet;   //!
   TBranch        *b_pzTCMet;   //!
   TBranch        *b_vertexXTCMet;   //!
   TBranch        *b_vertexYTCMet;   //!
   TBranch        *b_vertexZTCMet;   //!
   TBranch        *b_nPFMet;   //!
   TBranch        *b_chargePFMet;   //!
   TBranch        *b_energyPFMet;   //!
   TBranch        *b_thetaPFMet;   //!
   TBranch        *b_etaPFMet;   //!
   TBranch        *b_phiPFMet;   //!
   TBranch        *b_pxPFMet;   //!
   TBranch        *b_pyPFMet;   //!
   TBranch        *b_pzPFMet;   //!
   TBranch        *b_vertexXPFMet;   //!
   TBranch        *b_vertexYPFMet;   //!
   TBranch        *b_vertexZPFMet;   //!
   TBranch        *b_nGenMet;   //!
   TBranch        *b_chargeGenMet;   //!
   TBranch        *b_energyGenMet;   //!
   TBranch        *b_thetaGenMet;   //!
   TBranch        *b_etaGenMet;   //!
   TBranch        *b_phiGenMet;   //!
   TBranch        *b_pxGenMet;   //!
   TBranch        *b_pyGenMet;   //!
   TBranch        *b_pzGenMet;   //!
   TBranch        *b_vertexXGenMet;   //!
   TBranch        *b_vertexYGenMet;   //!
   TBranch        *b_vertexZGenMet;   //!
   TBranch        *b_nPFCand;   //!
   TBranch        *b_chargePFCand;   //!
   TBranch        *b_energyPFCand;   //!
   TBranch        *b_thetaPFCand;   //!
   TBranch        *b_etaPFCand;   //!
   TBranch        *b_phiPFCand;   //!
   TBranch        *b_pxPFCand;   //!
   TBranch        *b_pyPFCand;   //!
   TBranch        *b_pzPFCand;   //!
   TBranch        *b_vertexXPFCand;   //!
   TBranch        *b_vertexYPFCand;   //!
   TBranch        *b_vertexZPFCand;   //!
   TBranch        *b_particleTypePFCand;   //!
   TBranch        *b_iPFJetPFCand;   //!
   TBranch        *b_nAK5Jet;   //!
   TBranch        *b_chargeAK5Jet;   //!
   TBranch        *b_energyAK5Jet;   //!
   TBranch        *b_thetaAK5Jet;   //!
   TBranch        *b_etaAK5Jet;   //!
   TBranch        *b_phiAK5Jet;   //!
   TBranch        *b_pxAK5Jet;   //!
   TBranch        *b_pyAK5Jet;   //!
   TBranch        *b_pzAK5Jet;   //!
   TBranch        *b_vertexXAK5Jet;   //!
   TBranch        *b_vertexYAK5Jet;   //!
   TBranch        *b_vertexZAK5Jet;   //!
   TBranch        *b_emFracAK5Jet;   //!
   TBranch        *b_hadFracAK5Jet;   //!
   TBranch        *b_combinedSecondaryVertexBJetTagsAK5Jet;   //!
   TBranch        *b_combinedSecondaryVertexMVABJetTagsAK5Jet;   //!
   TBranch        *b_jetBProbabilityBJetTagsAK5Jet;   //!
   TBranch        *b_jetProbabilityBJetTagsAK5Jet;   //!
   TBranch        *b_simpleSecondaryVertexBJetTagsAK5Jet;   //!
   TBranch        *b_softMuonBJetTagsAK5Jet;   //!
   TBranch        *b_trackCountingHighPurBJetTagsAK5Jet;   //!
   TBranch        *b_trackCountingHighEffBJetTagsAK5Jet;   //!
   TBranch        *b_uncorrEnergyAK5Jet;   //!
   TBranch        *b_nAK5PFPUcorrJet;   //!
   TBranch        *b_chargeAK5PFPUcorrJet;   //!
   TBranch        *b_energyAK5PFPUcorrJet;   //!
   TBranch        *b_thetaAK5PFPUcorrJet;   //!
   TBranch        *b_etaAK5PFPUcorrJet;   //!
   TBranch        *b_phiAK5PFPUcorrJet;   //!
   TBranch        *b_pxAK5PFPUcorrJet;   //!
   TBranch        *b_pyAK5PFPUcorrJet;   //!
   TBranch        *b_pzAK5PFPUcorrJet;   //!
   TBranch        *b_vertexXAK5PFPUcorrJet;   //!
   TBranch        *b_vertexYAK5PFPUcorrJet;   //!
   TBranch        *b_vertexZAK5PFPUcorrJet;   //!
   TBranch        *b_chargedHadronEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_neutralHadronEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_photonEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_electronEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_muonEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_HFHadronEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_HFEMEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_chargedHadronMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_neutralHadronMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_photonMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_electronMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_muonMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_HFHadronMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_HFEMMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_chargedEmEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_neutralEmEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_combinedSecondaryVertexBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_jetBProbabilityBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_jetProbabilityBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softMuonBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softMuonByIP3dBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softMuonByPtBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softElectronBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softElectronByIP3dBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softElectronByPtBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_trackCountingHighPurBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_trackCountingHighEffBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_uncorrEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_ptDAK5PFPUcorrJet;   //!
   TBranch        *b_rmsCandAK5PFPUcorrJet;   //!
 //TBranch        *b_nAK5JPTJet;   //!
 //TBranch        *b_chargeAK5JPTJet;   //!
 //TBranch        *b_energyAK5JPTJet;   //!
 //TBranch        *b_thetaAK5JPTJet;   //!
 //TBranch        *b_etaAK5JPTJet;   //!
 //TBranch        *b_phiAK5JPTJet;   //!
 //TBranch        *b_pxAK5JPTJet;   //!
 //TBranch        *b_pyAK5JPTJet;   //!
 //TBranch        *b_pzAK5JPTJet;   //!
 //TBranch        *b_vertexXAK5JPTJet;   //!
 //TBranch        *b_vertexYAK5JPTJet;   //!
 //TBranch        *b_vertexZAK5JPTJet;   //!
 //TBranch        *b_emFracAK5JPTJet;   //!
 //TBranch        *b_hadFracAK5JPTJet;   //!
 //TBranch        *b_combinedSecondaryVertexBJetTagsAK5JPTJet;   //!
 //TBranch        *b_combinedSecondaryVertexMVABJetTagsAK5JPTJet;   //!
 //TBranch        *b_jetBProbabilityBJetTagsAK5JPTJet;   //!
 //TBranch        *b_jetProbabilityBJetTagsAK5JPTJet;   //!
 //TBranch        *b_simpleSecondaryVertexBJetTagsAK5JPTJet;   //!
 //TBranch        *b_softMuonBJetTagsAK5JPTJet;   //!
 //TBranch        *b_trackCountingHighPurBJetTagsAK5JPTJet;   //!
 //TBranch        *b_trackCountingHighEffBJetTagsAK5JPTJet;   //!
 //TBranch        *b_uncorrEnergyAK5JPTJet;   //!
   TBranch        *b_nAK5GenJet;   //!
   TBranch        *b_chargeAK5GenJet;   //!
   TBranch        *b_energyAK5GenJet;   //!
   TBranch        *b_thetaAK5GenJet;   //!
   TBranch        *b_etaAK5GenJet;   //!
   TBranch        *b_phiAK5GenJet;   //!
   TBranch        *b_pxAK5GenJet;   //!
   TBranch        *b_pyAK5GenJet;   //!
   TBranch        *b_pzAK5GenJet;   //!
   TBranch        *b_vertexXAK5GenJet;   //!
   TBranch        *b_vertexYAK5GenJet;   //!
   TBranch        *b_vertexZAK5GenJet;   //!
   TBranch        *b_genPtHat;   //!
   TBranch        *b_genProcessId;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genAlphaQCD;   //!
   TBranch        *b_genAlphaQED;   //!

   Ntp1Analyzer( const std::string& analyzerType, const std::string& dataset, const std::string& flags="", TTree *tree=0);
   virtual ~Ntp1Analyzer();

   virtual void SetFlags( const std::string& flags ) { flags_ = flags; };
   virtual void SetRequiredTriggers( const std::vector<std::string>& reqTrigz ) { requiredTriggers_ = reqTrigz; };
   virtual void AddRequiredTrigger( const std::string& trigger ) { requiredTriggers_.push_back(trigger); };
   virtual bool PassedHLT( const std::string& HLTName="" );

   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     LoadInput();
   virtual void     LoadInputFromFile( const std::string& fileName );
   virtual void     LoadTrigger( TFile* condFile=0 );
   virtual void     CreateOutputFile();
   virtual void     Init(TTree *tree);
   virtual void     Loop()=0;
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     ReadJSONFile(const std::string& json);
   virtual void     ReadCSVFile(const std::string& csv);
   virtual void     UpdateCache();
   virtual bool     isGoodEvent();
   virtual GenEventParameters     getGenEventParameters ();
   double           trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);


   //virtual void BookStuff()=0;

};

#endif
