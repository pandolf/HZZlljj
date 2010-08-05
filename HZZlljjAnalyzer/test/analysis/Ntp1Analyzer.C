#include "Ntp1Analyzer.h"
#include "TH1F.h"
#include <cstdlib>
#include <fstream>



Ntp1Analyzer::Ntp1Analyzer(const std::string& analyzerType, const std::string& dataset, const std::string& flags, TTree* tree)
{

   dataset_ = dataset;

   DEBUG_ = false;
   filterGoodRuns_ = false; //default: do not filter
   totalIntLumi_ = 0.;

   analyzerType_ = analyzerType;

   flags_ = flags;

}



// destructor
Ntp1Analyzer::~Ntp1Analyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   outfile_->cd();
   h1_nCounter_->Write();
   reducedTree_->Write();
   outfile_->Write();
   outfile_->Close();
   

}


void Ntp1Analyzer::LoadInput() {

   std::cout << "-> Loading input... (could take a while)" << std::endl;

   char treePath[400];
   TChain * chain = new TChain("ntp1","");
   if( dataset_=="Wenu_Summer10_START37_V5_S09_v1" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/Wenu_Summer10_START37_V5_S09_v1/default_*.root/ntp1");
   } else {
     sprintf(treePath, "%s/default_*.root/ntp1", dataset_.c_str());
   }

   int addInt = chain->Add(treePath);

   if( addInt==0 ) {
     std::cout << "Didn't find files to add for dataset: '" << dataset_ << "'. Looking for a list..." << std::endl;
     std::string fileName = "files_HZZlljj_" + dataset_ + ".txt";
     this->LoadInputFromFile(fileName);
   } else {
     TTree* tree = chain;
     std::cout << "-> Tree has " << tree->GetEntries() << " entries." << std::endl;
     Init(tree);
     this->CreateOutputFile();
   }

}


void Ntp1Analyzer::LoadInputFromFile( const std::string& fileName ) {

   FILE* iff = fopen(fileName.c_str(),"r");
   if(iff == 0) {
     std::cout << "cannot open input file " << fileName << " ... now exiting." << std::endl;
     exit(-1);
   }

   TChain * chain = new TChain("ntp1","");

   char singleLine[500];
 
   while( fscanf(iff, "%s", singleLine) !=EOF ) {
   
     std::string singleLine_str(singleLine);
     singleLine_str = singleLine_str + "/ntp1";
     std::cout << "-> Adding " << singleLine_str << std::endl;
     chain->Add(singleLine_str.c_str());

   }
   fclose(iff);

   TTree* tree = chain;
   std::cout << "-> Tree has " << tree->GetEntries() << " entries." << std::endl;
   Init(tree);
   this->CreateOutputFile();

}





void Ntp1Analyzer::CreateOutputFile() {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova2ndLevel_"+dataset_;
   else {
    if(dataset_!="") outfileName = analyzerType_ + "_2ndLevelTree_"+dataset_;
    else outfileName = analyzerType_ + "_2ndLevelTree";
   }


   if( flags_=="" )
     outfileName = outfileName + ".root";
   else 
     outfileName = outfileName + "_" + flags_ + ".root";

   outfile_ = TFile::Open(outfileName.c_str(), "RECREATE");
   
   outfile_->cd();

   reducedTree_ = new TTree("reducedTree", "Reduced Tree");
   reducedTree_->SetMaxTreeSize(100000000000ULL); //setting max tree size to 100 GB

   h1_nCounter_ = new TH1F("nCounter", "", 1, 0., 1.);

}


Int_t Ntp1Analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Ntp1Analyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Ntp1Analyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   if (!tree) return;

   GenEventParameters genPars = this->getGenEventParameters();
   ptHatMin_ = genPars.ptHatMin;
   ptHatMax_ = genPars.ptHatMax;

   //will cut on pt_hat, so have to divide only by correct number of events:
   char cutOnPtHat[70];
   sprintf( cutOnPtHat, "genPtHat>%lf && genPtHat<%lf", (Double_t)ptHatMin_, (Double_t)ptHatMax_);
   Int_t nEntries_cut = tree->GetEntries(cutOnPtHat);
   h1_nCounter_->SetBinContent( 1, nEntries_cut );


   // Set branch addresses and branch pointers
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   std::string branchName;

   fChain->SetBranchAddress("nl1Technical", &nl1Technical, &b_nl1Technical);
   fChain->SetBranchAddress("l1Technical", l1Technical, &b_l1Technical);
   fChain->SetBranchAddress("nl1Global", &nl1Global, &b_nl1Global);
   fChain->SetBranchAddress("l1Global", l1Global, &b_l1Global);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("nMc", &nMc, &b_nMc);
   fChain->SetBranchAddress("pMc", pMc, &b_pMc);
   fChain->SetBranchAddress("thetaMc", thetaMc, &b_thetaMc);
   fChain->SetBranchAddress("etaMc", etaMc, &b_etaMc);
   fChain->SetBranchAddress("phiMc", phiMc, &b_phiMc);
   fChain->SetBranchAddress("energyMc", energyMc, &b_energyMc);
   fChain->SetBranchAddress("idMc", idMc, &b_idMc);
   fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
   fChain->SetBranchAddress("statusMc", statusMc, &b_statusMc);
   fChain->SetBranchAddress("nTrg", &nTrg, &b_nTrg);
   fChain->SetBranchAddress("firedTrg", firedTrg, &b_firedTrg);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("energyEle", energyEle, &b_energyEle);
   fChain->SetBranchAddress("thetaEle", thetaEle, &b_thetaEle);
   fChain->SetBranchAddress("etaEle", etaEle, &b_etaEle);
   fChain->SetBranchAddress("phiEle", phiEle, &b_phiEle);
   fChain->SetBranchAddress("pxEle", pxEle, &b_pxEle);
   fChain->SetBranchAddress("pyEle", pyEle, &b_pyEle);
   fChain->SetBranchAddress("pzEle", pzEle, &b_pzEle);
   fChain->SetBranchAddress("vertexXEle", vertexXEle, &b_vertexXEle);
   fChain->SetBranchAddress("vertexYEle", vertexYEle, &b_vertexYEle);
   fChain->SetBranchAddress("vertexZEle", vertexZEle, &b_vertexZEle);
   fChain->SetBranchAddress("fiducialFlagsEle", fiducialFlagsEle, &b_fiducialFlagsEle);
   fChain->SetBranchAddress("recoFlagsEle", recoFlagsEle, &b_recoFlagsEle);
   fChain->SetBranchAddress("energyCorrectionsEle", energyCorrectionsEle, &b_energyCorrectionsEle);
   fChain->SetBranchAddress("esEnergyEle", esEnergyEle, &b_esEnergyEle);
   fChain->SetBranchAddress("superClusterIndexEle", superClusterIndexEle, &b_superClusterIndexEle);
   fChain->SetBranchAddress("PFsuperClusterIndexEle", PFsuperClusterIndexEle, &b_PFsuperClusterIndexEle);
   fChain->SetBranchAddress("trackIndexEle", trackIndexEle, &b_trackIndexEle);
   fChain->SetBranchAddress("gsfTrackIndexEle", gsfTrackIndexEle, &b_gsfTrackIndexEle);
   fChain->SetBranchAddress("convDistEle", convDistEle, &b_convDistEle);
   fChain->SetBranchAddress("convDcotEle", convDcotEle, &b_convDcotEle);
   fChain->SetBranchAddress("convRadiusEle", convRadiusEle, &b_convRadiusEle);
   fChain->SetBranchAddress("convTrackIndexEle", convTrackIndexEle, &b_convTrackIndexEle);
   fChain->SetBranchAddress("convXEle", convXEle, &b_convXEle);
   fChain->SetBranchAddress("convYEle", convYEle, &b_convYEle);
   fChain->SetBranchAddress("convZEle", convZEle, &b_convZEle);
   fChain->SetBranchAddress("convChi2ProbEle", convChi2ProbEle, &b_convChi2ProbEle);
   fChain->SetBranchAddress("scPixChargeEle", scPixChargeEle, &b_scPixChargeEle);
   fChain->SetBranchAddress("classificationEle", classificationEle, &b_classificationEle);
   fChain->SetBranchAddress("standardClassificationEle", standardClassificationEle, &b_standardClassificationEle);
   fChain->SetBranchAddress("fbremEle", fbremEle, &b_fbremEle);
   fChain->SetBranchAddress("nbremsEle", nbremsEle, &b_nbremsEle);
   fChain->SetBranchAddress("hOverEEle", hOverEEle, &b_hOverEEle);
   fChain->SetBranchAddress("eSuperClusterOverPEle", eSuperClusterOverPEle, &b_eSuperClusterOverPEle);
   fChain->SetBranchAddress("eSeedOverPoutEle", eSeedOverPoutEle, &b_eSeedOverPoutEle);
   fChain->SetBranchAddress("deltaEtaAtVtxEle", deltaEtaAtVtxEle, &b_deltaEtaAtVtxEle);
   fChain->SetBranchAddress("deltaPhiAtVtxEle", deltaPhiAtVtxEle, &b_deltaPhiAtVtxEle);
   fChain->SetBranchAddress("deltaEtaAtCaloEle", deltaEtaAtCaloEle, &b_deltaEtaAtCaloEle);
   fChain->SetBranchAddress("deltaPhiAtCaloEle", deltaPhiAtCaloEle, &b_deltaPhiAtCaloEle);
   fChain->SetBranchAddress("tipEle", tipEle, &b_tipEle);
   fChain->SetBranchAddress("dr03TkSumPtEle", dr03TkSumPtEle, &b_dr03TkSumPtEle);
   fChain->SetBranchAddress("dr03EcalRecHitSumEtEle", dr03EcalRecHitSumEtEle, &b_dr03EcalRecHitSumEtEle);
   fChain->SetBranchAddress("dr03HcalTowerSumEtEle", dr03HcalTowerSumEtEle, &b_dr03HcalTowerSumEtEle);
   fChain->SetBranchAddress("dr04TkSumPtEle", dr04TkSumPtEle, &b_dr04TkSumPtEle);
   fChain->SetBranchAddress("dr04EcalRecHitSumEtEle", dr04EcalRecHitSumEtEle, &b_dr04EcalRecHitSumEtEle);
   fChain->SetBranchAddress("dr04HcalTowerSumEtEle", dr04HcalTowerSumEtEle, &b_dr04HcalTowerSumEtEle);
   fChain->SetBranchAddress("scBasedEcalSum03Ele", scBasedEcalSum03Ele, &b_scBasedEcalSum03Ele);
   fChain->SetBranchAddress("scBasedEcalSum04Ele", scBasedEcalSum04Ele, &b_scBasedEcalSum04Ele);
   fChain->SetBranchAddress("eleIdCutsEle", eleIdCutsEle, &b_eleIdCutsEle);
   fChain->SetBranchAddress("eleIdLikelihoodEle", eleIdLikelihoodEle, &b_eleIdLikelihoodEle);
   fChain->SetBranchAddress("pflowMVAEle", pflowMVAEle, &b_pflowMVAEle);
   fChain->SetBranchAddress("nPFEle", &nPFEle, &b_nPFEle);
   fChain->SetBranchAddress("chargePFEle", chargePFEle, &b_chargePFEle);
   fChain->SetBranchAddress("energyPFEle", energyPFEle, &b_energyPFEle);
   fChain->SetBranchAddress("thetaPFEle", thetaPFEle, &b_thetaPFEle);
   fChain->SetBranchAddress("etaPFEle", etaPFEle, &b_etaPFEle);
   fChain->SetBranchAddress("phiPFEle", phiPFEle, &b_phiPFEle);
   fChain->SetBranchAddress("pxPFEle", pxPFEle, &b_pxPFEle);
   fChain->SetBranchAddress("pyPFEle", pyPFEle, &b_pyPFEle);
   fChain->SetBranchAddress("pzPFEle", pzPFEle, &b_pzPFEle);
   fChain->SetBranchAddress("vertexXPFEle", vertexXPFEle, &b_vertexXPFEle);
   fChain->SetBranchAddress("vertexYPFEle", vertexYPFEle, &b_vertexYPFEle);
   fChain->SetBranchAddress("vertexZPFEle", vertexZPFEle, &b_vertexZPFEle);
   fChain->SetBranchAddress("MvaOutputPFEle", MvaOutputPFEle, &b_MvaOutputPFEle);
   fChain->SetBranchAddress("PS1EnergyPFEle", PS1EnergyPFEle, &b_PS1EnergyPFEle);
   fChain->SetBranchAddress("PS2EnergyPFEle", PS2EnergyPFEle, &b_PS2EnergyPFEle);
   fChain->SetBranchAddress("EcalEnergyPFEle", EcalEnergyPFEle, &b_EcalEnergyPFEle);
   fChain->SetBranchAddress("HcalEnergyPFEle", HcalEnergyPFEle, &b_HcalEnergyPFEle);
   fChain->SetBranchAddress("RawEcalEnergyPFEle", RawEcalEnergyPFEle, &b_RawEcalEnergyPFEle);
   fChain->SetBranchAddress("RawHcalEnergyPFEle", RawHcalEnergyPFEle, &b_RawHcalEnergyPFEle);
   fChain->SetBranchAddress("PositionAtEcalXPFEle", PositionAtEcalXPFEle, &b_PositionAtEcalXPFEle);
   fChain->SetBranchAddress("PositionAtEcalYPFEle", PositionAtEcalYPFEle, &b_PositionAtEcalYPFEle);
   fChain->SetBranchAddress("PositionAtEcalZPFEle", PositionAtEcalZPFEle, &b_PositionAtEcalZPFEle);
   fChain->SetBranchAddress("gsfTrackIndexPFEle", gsfTrackIndexPFEle, &b_gsfTrackIndexPFEle);
   fChain->SetBranchAddress("trackIndexPFEle", trackIndexPFEle, &b_trackIndexPFEle);
   fChain->SetBranchAddress("chIso03vetoPFEle", chIso03vetoPFEle, &b_chIso03vetoPFEle);
   fChain->SetBranchAddress("chIso04vetoPFEle", chIso04vetoPFEle, &b_chIso04vetoPFEle);
   fChain->SetBranchAddress("chIso05vetoPFEle", chIso05vetoPFEle, &b_chIso05vetoPFEle);
   fChain->SetBranchAddress("chIso03noVetoPFEle", chIso03noVetoPFEle, &b_chIso03noVetoPFEle);
   fChain->SetBranchAddress("chIso04noVetoPFEle", chIso04noVetoPFEle, &b_chIso04noVetoPFEle);
   fChain->SetBranchAddress("chIso05noVetoPFEle", chIso05noVetoPFEle, &b_chIso05noVetoPFEle);
   fChain->SetBranchAddress("nhIso03vetoPFEle", nhIso03vetoPFEle, &b_nhIso03vetoPFEle);
   fChain->SetBranchAddress("nhIso04vetoPFEle", nhIso04vetoPFEle, &b_nhIso04vetoPFEle);
   fChain->SetBranchAddress("nhIso05vetoPFEle", nhIso05vetoPFEle, &b_nhIso05vetoPFEle);
   fChain->SetBranchAddress("nhIso03noVetoPFEle", nhIso03noVetoPFEle, &b_nhIso03noVetoPFEle);
   fChain->SetBranchAddress("nhIso04noVetoPFEle", nhIso04noVetoPFEle, &b_nhIso04noVetoPFEle);
   fChain->SetBranchAddress("nhIso05noVetoPFEle", nhIso05noVetoPFEle, &b_nhIso05noVetoPFEle);
   fChain->SetBranchAddress("phIso03vetoPFEle", phIso03vetoPFEle, &b_phIso03vetoPFEle);
   fChain->SetBranchAddress("phIso04vetoPFEle", phIso04vetoPFEle, &b_phIso04vetoPFEle);
   fChain->SetBranchAddress("phIso05vetoPFEle", phIso05vetoPFEle, &b_phIso05vetoPFEle);
   fChain->SetBranchAddress("phIso03noVetoPFEle", phIso03noVetoPFEle, &b_phIso03noVetoPFEle);
   fChain->SetBranchAddress("phIso04noVetoPFEle", phIso04noVetoPFEle, &b_phIso04noVetoPFEle);
   fChain->SetBranchAddress("phIso05noVetoPFEle", phIso05noVetoPFEle, &b_phIso05noVetoPFEle);
   fChain->SetBranchAddress("nSC", &nSC, &b_nSC);
   fChain->SetBranchAddress("nBCSC", nBCSC, &b_nBCSC);
   fChain->SetBranchAddress("nCrystalsSC", nCrystalsSC, &b_nCrystalsSC);
   fChain->SetBranchAddress("rawEnergySC", rawEnergySC, &b_rawEnergySC);
   fChain->SetBranchAddress("energySC", energySC, &b_energySC);
   fChain->SetBranchAddress("etaSC", etaSC, &b_etaSC);
   fChain->SetBranchAddress("thetaSC", thetaSC, &b_thetaSC);
   fChain->SetBranchAddress("phiSC", phiSC, &b_phiSC);
   fChain->SetBranchAddress("phiWidthSC", phiWidthSC, &b_phiWidthSC);
   fChain->SetBranchAddress("etaWidthSC", etaWidthSC, &b_etaWidthSC);
   fChain->SetBranchAddress("e3x3SC", e3x3SC, &b_e3x3SC);
   fChain->SetBranchAddress("e5x5SC", e5x5SC, &b_e5x5SC);
   fChain->SetBranchAddress("eMaxSC", eMaxSC, &b_eMaxSC);
   fChain->SetBranchAddress("e2x2SC", e2x2SC, &b_e2x2SC);
   fChain->SetBranchAddress("e2ndSC", e2ndSC, &b_e2ndSC);
   fChain->SetBranchAddress("e1x5SC", e1x5SC, &b_e1x5SC);
   fChain->SetBranchAddress("e2x5MaxSC", e2x5MaxSC, &b_e2x5MaxSC);
   fChain->SetBranchAddress("e4SwissCrossSC", e4SwissCrossSC, &b_e4SwissCrossSC);
   fChain->SetBranchAddress("covIEtaIEtaSC", covIEtaIEtaSC, &b_covIEtaIEtaSC);
   fChain->SetBranchAddress("covIEtaIPhiSC", covIEtaIPhiSC, &b_covIEtaIPhiSC);
   fChain->SetBranchAddress("covIPhiIPhiSC", covIPhiIPhiSC, &b_covIPhiIPhiSC);
   fChain->SetBranchAddress("hOverESC", hOverESC, &b_hOverESC);
   fChain->SetBranchAddress("recoFlagSC", recoFlagSC, &b_recoFlagSC);
   fChain->SetBranchAddress("channelStatusSC", channelStatusSC, &b_channelStatusSC);
   fChain->SetBranchAddress("timeSC", timeSC, &b_timeSC);
   fChain->SetBranchAddress("chi2SC", chi2SC, &b_chi2SC);
   fChain->SetBranchAddress("seedEnergySC", seedEnergySC, &b_seedEnergySC);
   fChain->SetBranchAddress("idClosProblSC", idClosProblSC, &b_idClosProblSC);
   fChain->SetBranchAddress("sevClosProblSC", sevClosProblSC, &b_sevClosProblSC);
   fChain->SetBranchAddress("fracClosProblSC", fracClosProblSC, &b_fracClosProblSC);
   fChain->SetBranchAddress("scBasedEcalSum03SC", scBasedEcalSum03SC, &b_scBasedEcalSum03SC);
   fChain->SetBranchAddress("scBasedEcalSum04SC", scBasedEcalSum04SC, &b_scBasedEcalSum04SC);
   fChain->SetBranchAddress("ecalRecHitSumEtConeDR03SC", ecalRecHitSumEtConeDR03SC, &b_ecalRecHitSumEtConeDR03SC);
   fChain->SetBranchAddress("hcalTowerSumEtConeDR03SC", hcalTowerSumEtConeDR03SC, &b_hcalTowerSumEtConeDR03SC);
   fChain->SetBranchAddress("trkSumPtSolidConeDR03SC", trkSumPtSolidConeDR03SC, &b_trkSumPtSolidConeDR03SC);
   fChain->SetBranchAddress("ecalRecHitSumEtConeDR04SC", ecalRecHitSumEtConeDR04SC, &b_ecalRecHitSumEtConeDR04SC);
   fChain->SetBranchAddress("hcalTowerSumEtConeDR04SC", hcalTowerSumEtConeDR04SC, &b_hcalTowerSumEtConeDR04SC);
   fChain->SetBranchAddress("trkSumPtSolidConeDR04SC", trkSumPtSolidConeDR04SC, &b_trkSumPtSolidConeDR04SC);
   fChain->SetBranchAddress("nPFSC", &nPFSC, &b_nPFSC);
   fChain->SetBranchAddress("nBCPFSC", nBCPFSC, &b_nBCPFSC);
   fChain->SetBranchAddress("nCrystalsPFSC", nCrystalsPFSC, &b_nCrystalsPFSC);
   fChain->SetBranchAddress("rawEnergyPFSC", rawEnergyPFSC, &b_rawEnergyPFSC);
   fChain->SetBranchAddress("energyPFSC", energyPFSC, &b_energyPFSC);
   fChain->SetBranchAddress("etaPFSC", etaPFSC, &b_etaPFSC);
   fChain->SetBranchAddress("thetaPFSC", thetaPFSC, &b_thetaPFSC);
   fChain->SetBranchAddress("phiPFSC", phiPFSC, &b_phiPFSC);
   fChain->SetBranchAddress("phiWidthPFSC", phiWidthPFSC, &b_phiWidthPFSC);
   fChain->SetBranchAddress("etaWidthPFSC", etaWidthPFSC, &b_etaWidthPFSC);
   fChain->SetBranchAddress("e3x3PFSC", e3x3PFSC, &b_e3x3PFSC);
   fChain->SetBranchAddress("e5x5PFSC", e5x5PFSC, &b_e5x5PFSC);
   fChain->SetBranchAddress("eMaxPFSC", eMaxPFSC, &b_eMaxPFSC);
   fChain->SetBranchAddress("e2x2PFSC", e2x2PFSC, &b_e2x2PFSC);
   fChain->SetBranchAddress("e2ndPFSC", e2ndPFSC, &b_e2ndPFSC);
   fChain->SetBranchAddress("e1x5PFSC", e1x5PFSC, &b_e1x5PFSC);
   fChain->SetBranchAddress("e2x5MaxPFSC", e2x5MaxPFSC, &b_e2x5MaxPFSC);
   fChain->SetBranchAddress("e4SwissCrossPFSC", e4SwissCrossPFSC, &b_e4SwissCrossPFSC);
   fChain->SetBranchAddress("covIEtaIEtaPFSC", covIEtaIEtaPFSC, &b_covIEtaIEtaPFSC);
   fChain->SetBranchAddress("covIEtaIPhiPFSC", covIEtaIPhiPFSC, &b_covIEtaIPhiPFSC);
   fChain->SetBranchAddress("covIPhiIPhiPFSC", covIPhiIPhiPFSC, &b_covIPhiIPhiPFSC);
   fChain->SetBranchAddress("hOverEPFSC", hOverEPFSC, &b_hOverEPFSC);
   fChain->SetBranchAddress("recoFlagPFSC", recoFlagPFSC, &b_recoFlagPFSC);
   fChain->SetBranchAddress("channelStatusPFSC", channelStatusPFSC, &b_channelStatusPFSC);
   fChain->SetBranchAddress("timePFSC", timePFSC, &b_timePFSC);
   fChain->SetBranchAddress("chi2PFSC", chi2PFSC, &b_chi2PFSC);
   fChain->SetBranchAddress("seedEnergyPFSC", seedEnergyPFSC, &b_seedEnergyPFSC);
   fChain->SetBranchAddress("idClosProblPFSC", idClosProblPFSC, &b_idClosProblPFSC);
   fChain->SetBranchAddress("sevClosProblPFSC", sevClosProblPFSC, &b_sevClosProblPFSC);
   fChain->SetBranchAddress("fracClosProblPFSC", fracClosProblPFSC, &b_fracClosProblPFSC);
   fChain->SetBranchAddress("scBasedEcalSum03PFSC", scBasedEcalSum03PFSC, &b_scBasedEcalSum03PFSC);
   fChain->SetBranchAddress("scBasedEcalSum04PFSC", scBasedEcalSum04PFSC, &b_scBasedEcalSum04PFSC);
   fChain->SetBranchAddress("nBC", &nBC, &b_nBC);
   fChain->SetBranchAddress("nCrystalsBC", nCrystalsBC, &b_nCrystalsBC);
   fChain->SetBranchAddress("energyBC", energyBC, &b_energyBC);
   fChain->SetBranchAddress("etaBC", etaBC, &b_etaBC);
   fChain->SetBranchAddress("thetaBC", thetaBC, &b_thetaBC);
   fChain->SetBranchAddress("phiBC", phiBC, &b_phiBC);
   fChain->SetBranchAddress("e3x3BC", e3x3BC, &b_e3x3BC);
   fChain->SetBranchAddress("e5x5BC", e5x5BC, &b_e5x5BC);
   fChain->SetBranchAddress("eMaxBC", eMaxBC, &b_eMaxBC);
   fChain->SetBranchAddress("e2x2BC", e2x2BC, &b_e2x2BC);
   fChain->SetBranchAddress("e2ndBC", e2ndBC, &b_e2ndBC);
   fChain->SetBranchAddress("covIEtaIEtaBC", covIEtaIEtaBC, &b_covIEtaIEtaBC);
   fChain->SetBranchAddress("covIEtaIPhiBC", covIEtaIPhiBC, &b_covIEtaIPhiBC);
   fChain->SetBranchAddress("covIPhiIPhiBC", covIPhiIPhiBC, &b_covIPhiIPhiBC);
   fChain->SetBranchAddress("recoFlagBC", recoFlagBC, &b_recoFlagBC);
   fChain->SetBranchAddress("timeBC", timeBC, &b_timeBC);
   fChain->SetBranchAddress("chi2BC", chi2BC, &b_chi2BC);
   fChain->SetBranchAddress("seedEnergyBC", seedEnergyBC, &b_seedEnergyBC);
   fChain->SetBranchAddress("idClosProblBC", idClosProblBC, &b_idClosProblBC);
   fChain->SetBranchAddress("sevClosProblBC", sevClosProblBC, &b_sevClosProblBC);
   fChain->SetBranchAddress("fracClosProblBC", fracClosProblBC, &b_fracClosProblBC);
   fChain->SetBranchAddress("indexSCBC", indexSCBC, &b_indexSCBC);
   fChain->SetBranchAddress("nTrack", &nTrack, &b_nTrack);
   fChain->SetBranchAddress("pxTrack", pxTrack, &b_pxTrack);
   fChain->SetBranchAddress("pyTrack", pyTrack, &b_pyTrack);
   fChain->SetBranchAddress("pzTrack", pzTrack, &b_pzTrack);
   fChain->SetBranchAddress("vtxIndexTrack", vtxIndexTrack, &b_vtxIndexTrack);
   fChain->SetBranchAddress("vtxWeightTrack", vtxWeightTrack, &b_vtxWeightTrack);
   fChain->SetBranchAddress("chargeTrack", chargeTrack, &b_chargeTrack);
   fChain->SetBranchAddress("ptErrorTrack", ptErrorTrack, &b_ptErrorTrack);
   fChain->SetBranchAddress("trackValidHitsTrack", trackValidHitsTrack, &b_trackValidHitsTrack);
   fChain->SetBranchAddress("trackLostHitsTrack", trackLostHitsTrack, &b_trackLostHitsTrack);
   fChain->SetBranchAddress("trackNormalizedChi2Track", trackNormalizedChi2Track, &b_trackNormalizedChi2Track);
   fChain->SetBranchAddress("qualityMaskTrack", qualityMaskTrack, &b_qualityMaskTrack);
   fChain->SetBranchAddress("impactPar3DTrack", impactPar3DTrack, &b_impactPar3DTrack);
   fChain->SetBranchAddress("impactPar3DErrorTrack", impactPar3DErrorTrack, &b_impactPar3DErrorTrack);
   fChain->SetBranchAddress("transvImpactParTrack", transvImpactParTrack, &b_transvImpactParTrack);
   fChain->SetBranchAddress("transvImpactParErrorTrack", transvImpactParErrorTrack, &b_transvImpactParErrorTrack);
   fChain->SetBranchAddress("trackVxTrack", trackVxTrack, &b_trackVxTrack);
   fChain->SetBranchAddress("trackVyTrack", trackVyTrack, &b_trackVyTrack);
   fChain->SetBranchAddress("trackVzTrack", trackVzTrack, &b_trackVzTrack);
   fChain->SetBranchAddress("pxAtOuterTrack", pxAtOuterTrack, &b_pxAtOuterTrack);
   fChain->SetBranchAddress("pyAtOuterTrack", pyAtOuterTrack, &b_pyAtOuterTrack);
   fChain->SetBranchAddress("pzAtOuterTrack", pzAtOuterTrack, &b_pzAtOuterTrack);
   fChain->SetBranchAddress("xAtOuterTrack", xAtOuterTrack, &b_xAtOuterTrack);
   fChain->SetBranchAddress("yAtOuterTrack", yAtOuterTrack, &b_yAtOuterTrack);
   fChain->SetBranchAddress("zAtOuterTrack", zAtOuterTrack, &b_zAtOuterTrack);
   fChain->SetBranchAddress("pxAtInnerTrack", pxAtInnerTrack, &b_pxAtInnerTrack);
   fChain->SetBranchAddress("pyAtInnerTrack", pyAtInnerTrack, &b_pyAtInnerTrack);
   fChain->SetBranchAddress("pzAtInnerTrack", pzAtInnerTrack, &b_pzAtInnerTrack);
   fChain->SetBranchAddress("xAtInnerTrack", xAtInnerTrack, &b_xAtInnerTrack);
   fChain->SetBranchAddress("yAtInnerTrack", yAtInnerTrack, &b_yAtInnerTrack);
   fChain->SetBranchAddress("zAtInnerTrack", zAtInnerTrack, &b_zAtInnerTrack);
   fChain->SetBranchAddress("recHitsSizeTrack", recHitsSizeTrack, &b_recHitsSizeTrack);
   fChain->SetBranchAddress("pixelHitsTrack", pixelHitsTrack, &b_pixelHitsTrack);
   fChain->SetBranchAddress("expInnerLayersTrack", expInnerLayersTrack, &b_expInnerLayersTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsTrack", numberOfValidPixelBarrelHitsTrack, &b_numberOfValidPixelBarrelHitsTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsTrack", numberOfValidPixelEndcapHitsTrack, &b_numberOfValidPixelEndcapHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsTrack", numberOfValidStripTIBHitsTrack, &b_numberOfValidStripTIBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsTrack", numberOfValidStripTIDHitsTrack, &b_numberOfValidStripTIDHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsTrack", numberOfValidStripTOBHitsTrack, &b_numberOfValidStripTOBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsTrack", numberOfValidStripTECHitsTrack, &b_numberOfValidStripTECHitsTrack);
   fChain->SetBranchAddress("truncatedDeDxTrack", truncatedDeDxTrack, &b_truncatedDeDxTrack);
   fChain->SetBranchAddress("truncatedDeDxErrorTrack", truncatedDeDxErrorTrack, &b_truncatedDeDxErrorTrack);
   fChain->SetBranchAddress("truncatedDeDxNoMTrack", truncatedDeDxNoMTrack, &b_truncatedDeDxNoMTrack);
   fChain->SetBranchAddress("medianDeDxTrack", medianDeDxTrack, &b_medianDeDxTrack);
   fChain->SetBranchAddress("medianDeDxErrorTrack", medianDeDxErrorTrack, &b_medianDeDxErrorTrack);
   fChain->SetBranchAddress("medianDeDxNoMTrack", medianDeDxNoMTrack, &b_medianDeDxNoMTrack);
   fChain->SetBranchAddress("harmonic2DeDxTrack", harmonic2DeDxTrack, &b_harmonic2DeDxTrack);
   fChain->SetBranchAddress("harmonic2DeDxErrorTrack", harmonic2DeDxErrorTrack, &b_harmonic2DeDxErrorTrack);
   fChain->SetBranchAddress("harmonic2DeDxNoMTrack", harmonic2DeDxNoMTrack, &b_harmonic2DeDxNoMTrack);
   fChain->SetBranchAddress("nGsfTrack", &nGsfTrack, &b_nGsfTrack);
   fChain->SetBranchAddress("pxGsfTrack", pxGsfTrack, &b_pxGsfTrack);
   fChain->SetBranchAddress("pyGsfTrack", pyGsfTrack, &b_pyGsfTrack);
   fChain->SetBranchAddress("pzGsfTrack", pzGsfTrack, &b_pzGsfTrack);
   fChain->SetBranchAddress("vtxIndexGsfTrack", vtxIndexGsfTrack, &b_vtxIndexGsfTrack);
   fChain->SetBranchAddress("vtxWeightGsfTrack", vtxWeightGsfTrack, &b_vtxWeightGsfTrack);
   fChain->SetBranchAddress("chargeGsfTrack", chargeGsfTrack, &b_chargeGsfTrack);
   fChain->SetBranchAddress("ptErrorGsfTrack", ptErrorGsfTrack, &b_ptErrorGsfTrack);
   fChain->SetBranchAddress("trackValidHitsGsfTrack", trackValidHitsGsfTrack, &b_trackValidHitsGsfTrack);
   fChain->SetBranchAddress("trackLostHitsGsfTrack", trackLostHitsGsfTrack, &b_trackLostHitsGsfTrack);
   fChain->SetBranchAddress("trackNormalizedChi2GsfTrack", trackNormalizedChi2GsfTrack, &b_trackNormalizedChi2GsfTrack);
   fChain->SetBranchAddress("qualityMaskGsfTrack", qualityMaskGsfTrack, &b_qualityMaskGsfTrack);
   fChain->SetBranchAddress("impactPar3DGsfTrack", impactPar3DGsfTrack, &b_impactPar3DGsfTrack);
   fChain->SetBranchAddress("impactPar3DErrorGsfTrack", impactPar3DErrorGsfTrack, &b_impactPar3DErrorGsfTrack);
   fChain->SetBranchAddress("transvImpactParGsfTrack", transvImpactParGsfTrack, &b_transvImpactParGsfTrack);
   fChain->SetBranchAddress("transvImpactParErrorGsfTrack", transvImpactParErrorGsfTrack, &b_transvImpactParErrorGsfTrack);
   fChain->SetBranchAddress("trackVxGsfTrack", trackVxGsfTrack, &b_trackVxGsfTrack);
   fChain->SetBranchAddress("trackVyGsfTrack", trackVyGsfTrack, &b_trackVyGsfTrack);
   fChain->SetBranchAddress("trackVzGsfTrack", trackVzGsfTrack, &b_trackVzGsfTrack);
   fChain->SetBranchAddress("pxAtOuterGsfTrack", pxAtOuterGsfTrack, &b_pxAtOuterGsfTrack);
   fChain->SetBranchAddress("pyAtOuterGsfTrack", pyAtOuterGsfTrack, &b_pyAtOuterGsfTrack);
   fChain->SetBranchAddress("pzAtOuterGsfTrack", pzAtOuterGsfTrack, &b_pzAtOuterGsfTrack);
   fChain->SetBranchAddress("xAtOuterGsfTrack", xAtOuterGsfTrack, &b_xAtOuterGsfTrack);
   fChain->SetBranchAddress("yAtOuterGsfTrack", yAtOuterGsfTrack, &b_yAtOuterGsfTrack);
   fChain->SetBranchAddress("zAtOuterGsfTrack", zAtOuterGsfTrack, &b_zAtOuterGsfTrack);
   fChain->SetBranchAddress("pxAtInnerGsfTrack", pxAtInnerGsfTrack, &b_pxAtInnerGsfTrack);
   fChain->SetBranchAddress("pyAtInnerGsfTrack", pyAtInnerGsfTrack, &b_pyAtInnerGsfTrack);
   fChain->SetBranchAddress("pzAtInnerGsfTrack", pzAtInnerGsfTrack, &b_pzAtInnerGsfTrack);
   fChain->SetBranchAddress("xAtInnerGsfTrack", xAtInnerGsfTrack, &b_xAtInnerGsfTrack);
   fChain->SetBranchAddress("yAtInnerGsfTrack", yAtInnerGsfTrack, &b_yAtInnerGsfTrack);
   fChain->SetBranchAddress("zAtInnerGsfTrack", zAtInnerGsfTrack, &b_zAtInnerGsfTrack);
   fChain->SetBranchAddress("recHitsSizeGsfTrack", recHitsSizeGsfTrack, &b_recHitsSizeGsfTrack);
   fChain->SetBranchAddress("pixelHitsGsfTrack", pixelHitsGsfTrack, &b_pixelHitsGsfTrack);
   fChain->SetBranchAddress("expInnerLayersGsfTrack", expInnerLayersGsfTrack, &b_expInnerLayersGsfTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsGsfTrack", numberOfValidPixelBarrelHitsGsfTrack, &b_numberOfValidPixelBarrelHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsGsfTrack", numberOfValidPixelEndcapHitsGsfTrack, &b_numberOfValidPixelEndcapHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsGsfTrack", numberOfValidStripTIBHitsGsfTrack, &b_numberOfValidStripTIBHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsGsfTrack", numberOfValidStripTIDHitsGsfTrack, &b_numberOfValidStripTIDHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsGsfTrack", numberOfValidStripTOBHitsGsfTrack, &b_numberOfValidStripTOBHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsGsfTrack", numberOfValidStripTECHitsGsfTrack, &b_numberOfValidStripTECHitsGsfTrack);
   fChain->SetBranchAddress("chargeModeGsfTrack", chargeModeGsfTrack, &b_chargeModeGsfTrack);
   fChain->SetBranchAddress("pxModeGsfTrack", pxModeGsfTrack, &b_pxModeGsfTrack);
   fChain->SetBranchAddress("pyModeGsfTrack", pyModeGsfTrack, &b_pyModeGsfTrack);
   fChain->SetBranchAddress("pzModeGsfTrack", pzModeGsfTrack, &b_pzModeGsfTrack);
   fChain->SetBranchAddress("recoFlagsGsfTrack", recoFlagsGsfTrack, &b_recoFlagsGsfTrack);
   fChain->SetBranchAddress("nGlobalMuonTrack", &nGlobalMuonTrack, &b_nGlobalMuonTrack);
   fChain->SetBranchAddress("pxGlobalMuonTrack", pxGlobalMuonTrack, &b_pxGlobalMuonTrack);
   fChain->SetBranchAddress("pyGlobalMuonTrack", pyGlobalMuonTrack, &b_pyGlobalMuonTrack);
   fChain->SetBranchAddress("pzGlobalMuonTrack", pzGlobalMuonTrack, &b_pzGlobalMuonTrack);
   fChain->SetBranchAddress("vtxIndexGlobalMuonTrack", vtxIndexGlobalMuonTrack, &b_vtxIndexGlobalMuonTrack);
   fChain->SetBranchAddress("vtxWeightGlobalMuonTrack", vtxWeightGlobalMuonTrack, &b_vtxWeightGlobalMuonTrack);
   fChain->SetBranchAddress("chargeGlobalMuonTrack", chargeGlobalMuonTrack, &b_chargeGlobalMuonTrack);
   fChain->SetBranchAddress("ptErrorGlobalMuonTrack", ptErrorGlobalMuonTrack, &b_ptErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackValidHitsGlobalMuonTrack", trackValidHitsGlobalMuonTrack, &b_trackValidHitsGlobalMuonTrack);
   fChain->SetBranchAddress("trackLostHitsGlobalMuonTrack", trackLostHitsGlobalMuonTrack, &b_trackLostHitsGlobalMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2GlobalMuonTrack", trackNormalizedChi2GlobalMuonTrack, &b_trackNormalizedChi2GlobalMuonTrack);
   fChain->SetBranchAddress("qualityMaskGlobalMuonTrack", qualityMaskGlobalMuonTrack, &b_qualityMaskGlobalMuonTrack);
   fChain->SetBranchAddress("impactPar3DGlobalMuonTrack", impactPar3DGlobalMuonTrack, &b_impactPar3DGlobalMuonTrack);
   fChain->SetBranchAddress("impactPar3DErrorGlobalMuonTrack", impactPar3DErrorGlobalMuonTrack, &b_impactPar3DErrorGlobalMuonTrack);
   fChain->SetBranchAddress("transvImpactParGlobalMuonTrack", transvImpactParGlobalMuonTrack, &b_transvImpactParGlobalMuonTrack);
   fChain->SetBranchAddress("transvImpactParErrorGlobalMuonTrack", transvImpactParErrorGlobalMuonTrack, &b_transvImpactParErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackVxGlobalMuonTrack", trackVxGlobalMuonTrack, &b_trackVxGlobalMuonTrack);
   fChain->SetBranchAddress("trackVyGlobalMuonTrack", trackVyGlobalMuonTrack, &b_trackVyGlobalMuonTrack);
   fChain->SetBranchAddress("trackVzGlobalMuonTrack", trackVzGlobalMuonTrack, &b_trackVzGlobalMuonTrack);
   fChain->SetBranchAddress("pxAtOuterGlobalMuonTrack", pxAtOuterGlobalMuonTrack, &b_pxAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("pyAtOuterGlobalMuonTrack", pyAtOuterGlobalMuonTrack, &b_pyAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("pzAtOuterGlobalMuonTrack", pzAtOuterGlobalMuonTrack, &b_pzAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("xAtOuterGlobalMuonTrack", xAtOuterGlobalMuonTrack, &b_xAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("yAtOuterGlobalMuonTrack", yAtOuterGlobalMuonTrack, &b_yAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("zAtOuterGlobalMuonTrack", zAtOuterGlobalMuonTrack, &b_zAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("pxAtInnerGlobalMuonTrack", pxAtInnerGlobalMuonTrack, &b_pxAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("pyAtInnerGlobalMuonTrack", pyAtInnerGlobalMuonTrack, &b_pyAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("pzAtInnerGlobalMuonTrack", pzAtInnerGlobalMuonTrack, &b_pzAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("xAtInnerGlobalMuonTrack", xAtInnerGlobalMuonTrack, &b_xAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("yAtInnerGlobalMuonTrack", yAtInnerGlobalMuonTrack, &b_yAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("zAtInnerGlobalMuonTrack", zAtInnerGlobalMuonTrack, &b_zAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("recHitsSizeGlobalMuonTrack", recHitsSizeGlobalMuonTrack, &b_recHitsSizeGlobalMuonTrack);
   fChain->SetBranchAddress("pixelHitsGlobalMuonTrack", pixelHitsGlobalMuonTrack, &b_pixelHitsGlobalMuonTrack);
   fChain->SetBranchAddress("expInnerLayersGlobalMuonTrack", expInnerLayersGlobalMuonTrack, &b_expInnerLayersGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsGlobalMuonTrack", numberOfValidPixelBarrelHitsGlobalMuonTrack, &b_numberOfValidPixelBarrelHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsGlobalMuonTrack", numberOfValidPixelEndcapHitsGlobalMuonTrack, &b_numberOfValidPixelEndcapHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsGlobalMuonTrack", numberOfValidStripTIBHitsGlobalMuonTrack, &b_numberOfValidStripTIBHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsGlobalMuonTrack", numberOfValidStripTIDHitsGlobalMuonTrack, &b_numberOfValidStripTIDHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsGlobalMuonTrack", numberOfValidStripTOBHitsGlobalMuonTrack, &b_numberOfValidStripTOBHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsGlobalMuonTrack", numberOfValidStripTECHitsGlobalMuonTrack, &b_numberOfValidStripTECHitsGlobalMuonTrack);
   fChain->SetBranchAddress("nSTAMuonTrack", &nSTAMuonTrack, &b_nSTAMuonTrack);
   fChain->SetBranchAddress("pxSTAMuonTrack", pxSTAMuonTrack, &b_pxSTAMuonTrack);
   fChain->SetBranchAddress("pySTAMuonTrack", pySTAMuonTrack, &b_pySTAMuonTrack);
   fChain->SetBranchAddress("pzSTAMuonTrack", pzSTAMuonTrack, &b_pzSTAMuonTrack);
   fChain->SetBranchAddress("vtxIndexSTAMuonTrack", vtxIndexSTAMuonTrack, &b_vtxIndexSTAMuonTrack);
   fChain->SetBranchAddress("vtxWeightSTAMuonTrack", vtxWeightSTAMuonTrack, &b_vtxWeightSTAMuonTrack);
   fChain->SetBranchAddress("chargeSTAMuonTrack", chargeSTAMuonTrack, &b_chargeSTAMuonTrack);
   fChain->SetBranchAddress("ptErrorSTAMuonTrack", ptErrorSTAMuonTrack, &b_ptErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackValidHitsSTAMuonTrack", trackValidHitsSTAMuonTrack, &b_trackValidHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackLostHitsSTAMuonTrack", trackLostHitsSTAMuonTrack, &b_trackLostHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2STAMuonTrack", trackNormalizedChi2STAMuonTrack, &b_trackNormalizedChi2STAMuonTrack);
   fChain->SetBranchAddress("qualityMaskSTAMuonTrack", qualityMaskSTAMuonTrack, &b_qualityMaskSTAMuonTrack);
   fChain->SetBranchAddress("impactPar3DSTAMuonTrack", impactPar3DSTAMuonTrack, &b_impactPar3DSTAMuonTrack);
   fChain->SetBranchAddress("impactPar3DErrorSTAMuonTrack", impactPar3DErrorSTAMuonTrack, &b_impactPar3DErrorSTAMuonTrack);
   fChain->SetBranchAddress("transvImpactParSTAMuonTrack", transvImpactParSTAMuonTrack, &b_transvImpactParSTAMuonTrack);
   fChain->SetBranchAddress("transvImpactParErrorSTAMuonTrack", transvImpactParErrorSTAMuonTrack, &b_transvImpactParErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackVxSTAMuonTrack", trackVxSTAMuonTrack, &b_trackVxSTAMuonTrack);
   fChain->SetBranchAddress("trackVySTAMuonTrack", trackVySTAMuonTrack, &b_trackVySTAMuonTrack);
   fChain->SetBranchAddress("trackVzSTAMuonTrack", trackVzSTAMuonTrack, &b_trackVzSTAMuonTrack);
   fChain->SetBranchAddress("pxAtOuterSTAMuonTrack", pxAtOuterSTAMuonTrack, &b_pxAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("pyAtOuterSTAMuonTrack", pyAtOuterSTAMuonTrack, &b_pyAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("pzAtOuterSTAMuonTrack", pzAtOuterSTAMuonTrack, &b_pzAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("xAtOuterSTAMuonTrack", xAtOuterSTAMuonTrack, &b_xAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("yAtOuterSTAMuonTrack", yAtOuterSTAMuonTrack, &b_yAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("zAtOuterSTAMuonTrack", zAtOuterSTAMuonTrack, &b_zAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("pxAtInnerSTAMuonTrack", pxAtInnerSTAMuonTrack, &b_pxAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("pyAtInnerSTAMuonTrack", pyAtInnerSTAMuonTrack, &b_pyAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("pzAtInnerSTAMuonTrack", pzAtInnerSTAMuonTrack, &b_pzAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("xAtInnerSTAMuonTrack", xAtInnerSTAMuonTrack, &b_xAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("yAtInnerSTAMuonTrack", yAtInnerSTAMuonTrack, &b_yAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("zAtInnerSTAMuonTrack", zAtInnerSTAMuonTrack, &b_zAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("recHitsSizeSTAMuonTrack", recHitsSizeSTAMuonTrack, &b_recHitsSizeSTAMuonTrack);
   fChain->SetBranchAddress("pixelHitsSTAMuonTrack", pixelHitsSTAMuonTrack, &b_pixelHitsSTAMuonTrack);
   fChain->SetBranchAddress("expInnerLayersSTAMuonTrack", expInnerLayersSTAMuonTrack, &b_expInnerLayersSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsSTAMuonTrack", numberOfValidPixelBarrelHitsSTAMuonTrack, &b_numberOfValidPixelBarrelHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsSTAMuonTrack", numberOfValidPixelEndcapHitsSTAMuonTrack, &b_numberOfValidPixelEndcapHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsSTAMuonTrack", numberOfValidStripTIBHitsSTAMuonTrack, &b_numberOfValidStripTIBHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsSTAMuonTrack", numberOfValidStripTIDHitsSTAMuonTrack, &b_numberOfValidStripTIDHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsSTAMuonTrack", numberOfValidStripTOBHitsSTAMuonTrack, &b_numberOfValidStripTOBHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsSTAMuonTrack", numberOfValidStripTECHitsSTAMuonTrack, &b_numberOfValidStripTECHitsSTAMuonTrack);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVxPV", PVxPV, &b_PVxPV);
   fChain->SetBranchAddress("PVyPV", PVyPV, &b_PVyPV);
   fChain->SetBranchAddress("PVzPV", PVzPV, &b_PVzPV);
   fChain->SetBranchAddress("PVErrxPV", PVErrxPV, &b_PVErrxPV);
   fChain->SetBranchAddress("PVErryPV", PVErryPV, &b_PVErryPV);
   fChain->SetBranchAddress("PVErrzPV", PVErrzPV, &b_PVErrzPV);
   fChain->SetBranchAddress("SumPtPV", SumPtPV, &b_SumPtPV);
   fChain->SetBranchAddress("ndofPV", ndofPV, &b_ndofPV);
   fChain->SetBranchAddress("chi2PV", chi2PV, &b_chi2PV);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("chargeMuon", chargeMuon, &b_chargeMuon);
   fChain->SetBranchAddress("energyMuon", energyMuon, &b_energyMuon);
   fChain->SetBranchAddress("thetaMuon", thetaMuon, &b_thetaMuon);
   fChain->SetBranchAddress("etaMuon", etaMuon, &b_etaMuon);
   fChain->SetBranchAddress("phiMuon", phiMuon, &b_phiMuon);
   fChain->SetBranchAddress("pxMuon", pxMuon, &b_pxMuon);
   fChain->SetBranchAddress("pyMuon", pyMuon, &b_pyMuon);
   fChain->SetBranchAddress("pzMuon", pzMuon, &b_pzMuon);
   fChain->SetBranchAddress("vertexXMuon", vertexXMuon, &b_vertexXMuon);
   fChain->SetBranchAddress("vertexYMuon", vertexYMuon, &b_vertexYMuon);
   fChain->SetBranchAddress("vertexZMuon", vertexZMuon, &b_vertexZMuon);
   fChain->SetBranchAddress("trackIndexMuon", trackIndexMuon, &b_trackIndexMuon);
   fChain->SetBranchAddress("standAloneTrackIndexMuon", standAloneTrackIndexMuon, &b_standAloneTrackIndexMuon);
   fChain->SetBranchAddress("combinedTrackIndexMuon", combinedTrackIndexMuon, &b_combinedTrackIndexMuon);
   fChain->SetBranchAddress("muonIdMuon", muonIdMuon, &b_muonIdMuon);
   fChain->SetBranchAddress("sumPt03Muon", sumPt03Muon, &b_sumPt03Muon);
   fChain->SetBranchAddress("emEt03Muon", emEt03Muon, &b_emEt03Muon);
   fChain->SetBranchAddress("hadEt03Muon", hadEt03Muon, &b_hadEt03Muon);
   fChain->SetBranchAddress("hoEt03Muon", hoEt03Muon, &b_hoEt03Muon);
   fChain->SetBranchAddress("nTrk03Muon", nTrk03Muon, &b_nTrk03Muon);
   fChain->SetBranchAddress("nJets03Muon", nJets03Muon, &b_nJets03Muon);
   fChain->SetBranchAddress("sumPt05Muon", sumPt05Muon, &b_sumPt05Muon);
   fChain->SetBranchAddress("emEt05Muon", emEt05Muon, &b_emEt05Muon);
   fChain->SetBranchAddress("hadEt05Muon", hadEt05Muon, &b_hadEt05Muon);
   fChain->SetBranchAddress("hoEt05Muon", hoEt05Muon, &b_hoEt05Muon);
   fChain->SetBranchAddress("nTrk05Muon", nTrk05Muon, &b_nTrk05Muon);
   fChain->SetBranchAddress("nJets05Muon", nJets05Muon, &b_nJets05Muon);
   fChain->SetBranchAddress("EcalExpDepoMuon", EcalExpDepoMuon, &b_EcalExpDepoMuon);
   fChain->SetBranchAddress("HcalExpDepoMuon", HcalExpDepoMuon, &b_HcalExpDepoMuon);
   fChain->SetBranchAddress("HoExpDepoMuon", HoExpDepoMuon, &b_HoExpDepoMuon);
   fChain->SetBranchAddress("emS9Muon", emS9Muon, &b_emS9Muon);
   fChain->SetBranchAddress("hadS9Muon", hadS9Muon, &b_hadS9Muon);
   fChain->SetBranchAddress("hoS9Muon", hoS9Muon, &b_hoS9Muon);
   fChain->SetBranchAddress("CaloCompMuon", CaloCompMuon, &b_CaloCompMuon);
   fChain->SetBranchAddress("nMet", &nMet, &b_nMet);
   fChain->SetBranchAddress("chargeMet", chargeMet, &b_chargeMet);
   fChain->SetBranchAddress("energyMet", energyMet, &b_energyMet);
   fChain->SetBranchAddress("thetaMet", thetaMet, &b_thetaMet);
   fChain->SetBranchAddress("etaMet", etaMet, &b_etaMet);
   fChain->SetBranchAddress("phiMet", phiMet, &b_phiMet);
   fChain->SetBranchAddress("pxMet", pxMet, &b_pxMet);
   fChain->SetBranchAddress("pyMet", pyMet, &b_pyMet);
   fChain->SetBranchAddress("pzMet", pzMet, &b_pzMet);
   fChain->SetBranchAddress("vertexXMet", vertexXMet, &b_vertexXMet);
   fChain->SetBranchAddress("vertexYMet", vertexYMet, &b_vertexYMet);
   fChain->SetBranchAddress("vertexZMet", vertexZMet, &b_vertexZMet);
   fChain->SetBranchAddress("nTCMet", &nTCMet, &b_nTCMet);
   fChain->SetBranchAddress("chargeTCMet", chargeTCMet, &b_chargeTCMet);
   fChain->SetBranchAddress("energyTCMet", energyTCMet, &b_energyTCMet);
   fChain->SetBranchAddress("thetaTCMet", thetaTCMet, &b_thetaTCMet);
   fChain->SetBranchAddress("etaTCMet", etaTCMet, &b_etaTCMet);
   fChain->SetBranchAddress("phiTCMet", phiTCMet, &b_phiTCMet);
   fChain->SetBranchAddress("pxTCMet", pxTCMet, &b_pxTCMet);
   fChain->SetBranchAddress("pyTCMet", pyTCMet, &b_pyTCMet);
   fChain->SetBranchAddress("pzTCMet", pzTCMet, &b_pzTCMet);
   fChain->SetBranchAddress("vertexXTCMet", vertexXTCMet, &b_vertexXTCMet);
   fChain->SetBranchAddress("vertexYTCMet", vertexYTCMet, &b_vertexYTCMet);
   fChain->SetBranchAddress("vertexZTCMet", vertexZTCMet, &b_vertexZTCMet);
   fChain->SetBranchAddress("nPFMet", &nPFMet, &b_nPFMet);
   fChain->SetBranchAddress("chargePFMet", chargePFMet, &b_chargePFMet);
   fChain->SetBranchAddress("energyPFMet", energyPFMet, &b_energyPFMet);
   fChain->SetBranchAddress("thetaPFMet", thetaPFMet, &b_thetaPFMet);
   fChain->SetBranchAddress("etaPFMet", etaPFMet, &b_etaPFMet);
   fChain->SetBranchAddress("phiPFMet", phiPFMet, &b_phiPFMet);
   fChain->SetBranchAddress("pxPFMet", pxPFMet, &b_pxPFMet);
   fChain->SetBranchAddress("pyPFMet", pyPFMet, &b_pyPFMet);
   fChain->SetBranchAddress("pzPFMet", pzPFMet, &b_pzPFMet);
   fChain->SetBranchAddress("vertexXPFMet", vertexXPFMet, &b_vertexXPFMet);
   fChain->SetBranchAddress("vertexYPFMet", vertexYPFMet, &b_vertexYPFMet);
   fChain->SetBranchAddress("vertexZPFMet", vertexZPFMet, &b_vertexZPFMet);
   fChain->SetBranchAddress("nGenMet", &nGenMet, &b_nGenMet);
   fChain->SetBranchAddress("chargeGenMet", chargeGenMet, &b_chargeGenMet);
   fChain->SetBranchAddress("energyGenMet", energyGenMet, &b_energyGenMet);
   fChain->SetBranchAddress("thetaGenMet", thetaGenMet, &b_thetaGenMet);
   fChain->SetBranchAddress("etaGenMet", etaGenMet, &b_etaGenMet);
   fChain->SetBranchAddress("phiGenMet", phiGenMet, &b_phiGenMet);
   fChain->SetBranchAddress("pxGenMet", pxGenMet, &b_pxGenMet);
   fChain->SetBranchAddress("pyGenMet", pyGenMet, &b_pyGenMet);
   fChain->SetBranchAddress("pzGenMet", pzGenMet, &b_pzGenMet);
   fChain->SetBranchAddress("vertexXGenMet", vertexXGenMet, &b_vertexXGenMet);
   fChain->SetBranchAddress("vertexYGenMet", vertexYGenMet, &b_vertexYGenMet);
   fChain->SetBranchAddress("vertexZGenMet", vertexZGenMet, &b_vertexZGenMet);
   fChain->SetBranchAddress("nAK5Jet", &nAK5Jet, &b_nAK5Jet);
   fChain->SetBranchAddress("chargeAK5Jet", chargeAK5Jet, &b_chargeAK5Jet);
   fChain->SetBranchAddress("energyAK5Jet", energyAK5Jet, &b_energyAK5Jet);
   fChain->SetBranchAddress("thetaAK5Jet", thetaAK5Jet, &b_thetaAK5Jet);
   fChain->SetBranchAddress("etaAK5Jet", etaAK5Jet, &b_etaAK5Jet);
   fChain->SetBranchAddress("phiAK5Jet", phiAK5Jet, &b_phiAK5Jet);
   fChain->SetBranchAddress("pxAK5Jet", pxAK5Jet, &b_pxAK5Jet);
   fChain->SetBranchAddress("pyAK5Jet", pyAK5Jet, &b_pyAK5Jet);
   fChain->SetBranchAddress("pzAK5Jet", pzAK5Jet, &b_pzAK5Jet);
   fChain->SetBranchAddress("vertexXAK5Jet", vertexXAK5Jet, &b_vertexXAK5Jet);
   fChain->SetBranchAddress("vertexYAK5Jet", vertexYAK5Jet, &b_vertexYAK5Jet);
   fChain->SetBranchAddress("vertexZAK5Jet", vertexZAK5Jet, &b_vertexZAK5Jet);
   fChain->SetBranchAddress("emFracAK5Jet", emFracAK5Jet, &b_emFracAK5Jet);
   fChain->SetBranchAddress("hadFracAK5Jet", hadFracAK5Jet, &b_hadFracAK5Jet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5Jet", combinedSecondaryVertexBJetTagsAK5Jet, &b_combinedSecondaryVertexBJetTagsAK5Jet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5Jet", combinedSecondaryVertexMVABJetTagsAK5Jet, &b_combinedSecondaryVertexMVABJetTagsAK5Jet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5Jet", jetBProbabilityBJetTagsAK5Jet, &b_jetBProbabilityBJetTagsAK5Jet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5Jet", jetProbabilityBJetTagsAK5Jet, &b_jetProbabilityBJetTagsAK5Jet);
   fChain->SetBranchAddress("simpleSecondaryVertexBJetTagsAK5Jet", simpleSecondaryVertexBJetTagsAK5Jet, &b_simpleSecondaryVertexBJetTagsAK5Jet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5Jet", softMuonBJetTagsAK5Jet, &b_softMuonBJetTagsAK5Jet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5Jet", trackCountingHighPurBJetTagsAK5Jet, &b_trackCountingHighPurBJetTagsAK5Jet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5Jet", trackCountingHighEffBJetTagsAK5Jet, &b_trackCountingHighEffBJetTagsAK5Jet);
   fChain->SetBranchAddress("uncorrEnergyAK5Jet", uncorrEnergyAK5Jet, &b_uncorrEnergyAK5Jet);
   fChain->SetBranchAddress("nAK5PFJet", &nAK5PFJet, &b_nAK5PFJet);
   fChain->SetBranchAddress("chargeAK5PFJet", chargeAK5PFJet, &b_chargeAK5PFJet);
   fChain->SetBranchAddress("energyAK5PFJet", energyAK5PFJet, &b_energyAK5PFJet);
   fChain->SetBranchAddress("thetaAK5PFJet", thetaAK5PFJet, &b_thetaAK5PFJet);
   fChain->SetBranchAddress("etaAK5PFJet", etaAK5PFJet, &b_etaAK5PFJet);
   fChain->SetBranchAddress("phiAK5PFJet", phiAK5PFJet, &b_phiAK5PFJet);
   fChain->SetBranchAddress("pxAK5PFJet", pxAK5PFJet, &b_pxAK5PFJet);
   fChain->SetBranchAddress("pyAK5PFJet", pyAK5PFJet, &b_pyAK5PFJet);
   fChain->SetBranchAddress("pzAK5PFJet", pzAK5PFJet, &b_pzAK5PFJet);
   fChain->SetBranchAddress("vertexXAK5PFJet", vertexXAK5PFJet, &b_vertexXAK5PFJet);
   fChain->SetBranchAddress("vertexYAK5PFJet", vertexYAK5PFJet, &b_vertexYAK5PFJet);
   fChain->SetBranchAddress("vertexZAK5PFJet", vertexZAK5PFJet, &b_vertexZAK5PFJet);
   fChain->SetBranchAddress("chargedHadronEnergyAK5PFJet", chargedHadronEnergyAK5PFJet, &b_chargedHadronEnergyAK5PFJet);
   fChain->SetBranchAddress("neutralHadronEnergyAK5PFJet", neutralHadronEnergyAK5PFJet, &b_neutralHadronEnergyAK5PFJet);
   fChain->SetBranchAddress("chargedEmEnergyAK5PFJet", chargedEmEnergyAK5PFJet, &b_chargedEmEnergyAK5PFJet);
   fChain->SetBranchAddress("neutralEmEnergyAK5PFJet", neutralEmEnergyAK5PFJet, &b_neutralEmEnergyAK5PFJet);
   fChain->SetBranchAddress("neutralMultiplicityAK5PFJet", neutralMultiplicityAK5PFJet, &b_neutralMultiplicityAK5PFJet);
   fChain->SetBranchAddress("chargedMultiplicityAK5PFJet", chargedMultiplicityAK5PFJet, &b_chargedMultiplicityAK5PFJet);
   fChain->SetBranchAddress("muonMultiplicityAK5PFJet", muonMultiplicityAK5PFJet, &b_muonMultiplicityAK5PFJet);
   fChain->SetBranchAddress("uncorrEnergyAK5PFJet", uncorrEnergyAK5PFJet, &b_uncorrEnergyAK5PFJet);
   fChain->SetBranchAddress("nAK5JPTJet", &nAK5JPTJet, &b_nAK5JPTJet);
   fChain->SetBranchAddress("chargeAK5JPTJet", chargeAK5JPTJet, &b_chargeAK5JPTJet);
   fChain->SetBranchAddress("energyAK5JPTJet", energyAK5JPTJet, &b_energyAK5JPTJet);
   fChain->SetBranchAddress("thetaAK5JPTJet", thetaAK5JPTJet, &b_thetaAK5JPTJet);
   fChain->SetBranchAddress("etaAK5JPTJet", etaAK5JPTJet, &b_etaAK5JPTJet);
   fChain->SetBranchAddress("phiAK5JPTJet", phiAK5JPTJet, &b_phiAK5JPTJet);
   fChain->SetBranchAddress("pxAK5JPTJet", pxAK5JPTJet, &b_pxAK5JPTJet);
   fChain->SetBranchAddress("pyAK5JPTJet", pyAK5JPTJet, &b_pyAK5JPTJet);
   fChain->SetBranchAddress("pzAK5JPTJet", pzAK5JPTJet, &b_pzAK5JPTJet);
   fChain->SetBranchAddress("vertexXAK5JPTJet", vertexXAK5JPTJet, &b_vertexXAK5JPTJet);
   fChain->SetBranchAddress("vertexYAK5JPTJet", vertexYAK5JPTJet, &b_vertexYAK5JPTJet);
   fChain->SetBranchAddress("vertexZAK5JPTJet", vertexZAK5JPTJet, &b_vertexZAK5JPTJet);
   fChain->SetBranchAddress("emFracAK5JPTJet", emFracAK5JPTJet, &b_emFracAK5JPTJet);
   fChain->SetBranchAddress("hadFracAK5JPTJet", hadFracAK5JPTJet, &b_hadFracAK5JPTJet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5JPTJet", combinedSecondaryVertexBJetTagsAK5JPTJet, &b_combinedSecondaryVertexBJetTagsAK5JPTJet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5JPTJet", combinedSecondaryVertexMVABJetTagsAK5JPTJet, &b_combinedSecondaryVertexMVABJetTagsAK5JPTJet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5JPTJet", jetBProbabilityBJetTagsAK5JPTJet, &b_jetBProbabilityBJetTagsAK5JPTJet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5JPTJet", jetProbabilityBJetTagsAK5JPTJet, &b_jetProbabilityBJetTagsAK5JPTJet);
   fChain->SetBranchAddress("simpleSecondaryVertexBJetTagsAK5JPTJet", simpleSecondaryVertexBJetTagsAK5JPTJet, &b_simpleSecondaryVertexBJetTagsAK5JPTJet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5JPTJet", softMuonBJetTagsAK5JPTJet, &b_softMuonBJetTagsAK5JPTJet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5JPTJet", trackCountingHighPurBJetTagsAK5JPTJet, &b_trackCountingHighPurBJetTagsAK5JPTJet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5JPTJet", trackCountingHighEffBJetTagsAK5JPTJet, &b_trackCountingHighEffBJetTagsAK5JPTJet);
   fChain->SetBranchAddress("uncorrEnergyAK5JPTJet", uncorrEnergyAK5JPTJet, &b_uncorrEnergyAK5JPTJet);
   fChain->SetBranchAddress("nAK5GenJet", &nAK5GenJet, &b_nAK5GenJet);
   fChain->SetBranchAddress("chargeAK5GenJet", chargeAK5GenJet, &b_chargeAK5GenJet);
   fChain->SetBranchAddress("energyAK5GenJet", energyAK5GenJet, &b_energyAK5GenJet);
   fChain->SetBranchAddress("thetaAK5GenJet", thetaAK5GenJet, &b_thetaAK5GenJet);
   fChain->SetBranchAddress("etaAK5GenJet", etaAK5GenJet, &b_etaAK5GenJet);
   fChain->SetBranchAddress("phiAK5GenJet", phiAK5GenJet, &b_phiAK5GenJet);
   fChain->SetBranchAddress("pxAK5GenJet", pxAK5GenJet, &b_pxAK5GenJet);
   fChain->SetBranchAddress("pyAK5GenJet", pyAK5GenJet, &b_pyAK5GenJet);
   fChain->SetBranchAddress("pzAK5GenJet", pzAK5GenJet, &b_pzAK5GenJet);
   fChain->SetBranchAddress("vertexXAK5GenJet", vertexXAK5GenJet, &b_vertexXAK5GenJet);
   fChain->SetBranchAddress("vertexYAK5GenJet", vertexYAK5GenJet, &b_vertexYAK5GenJet);
   fChain->SetBranchAddress("vertexZAK5GenJet", vertexZAK5GenJet, &b_vertexZAK5GenJet);
   fChain->SetBranchAddress("genPtHat", &genPtHat, &b_genPtHat);
   fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genAlphaQCD", &genAlphaQCD, &b_genAlphaQCD);
   fChain->SetBranchAddress("genAlphaQED", &genAlphaQED, &b_genAlphaQED);

   Notify();
}

Bool_t Ntp1Analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Ntp1Analyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Ntp1Analyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void Ntp1Analyzer::UpdateCache() {
     // cache current run
     if( oldrun_ != runNumber ) {
       oldrun_ = runNumber;
       goodLSCache_ = goodLS_.find( runNumber );
     }
}


void Ntp1Analyzer::ReadCSVFile(const std::string& csv) {

  std::cout << "Reading CSV file of LS luminosities " << csv << std::endl;
  FILE* iff = fopen(csv.c_str(),"r");

  if(iff == 0) {
    std::cout << "cannot open CSV file " << csv << " ... now exiting." << std::endl;
    exit(-2);
  }

  int irun, iLS, iHFcounts, iVtxcounts;
  double iHFlumi, iVtxlumi;

  while( fscanf(iff,"%*[ \t\n]%d,%*[ \t\n]%d,%*[ \t\n]%d,%*[ \t\n]%d,%*[ \t\n]%lf,%*[ \t\n]%lf", &irun, &iLS, &iHFcounts, &iVtxcounts, &iHFlumi, &iVtxlumi) != EOF ) {

    RunLumiPair rlpair( std::pair<int, int>(irun, iLS) );
    double istLumi = (iVtxlumi>0.) ? iVtxlumi : iHFlumi; //using Vtx lumi as default
    double intLumi = 23.*istLumi/1.0e30; //in microb-1 (1 LS = 23 seconds)
    LSLumimap_[rlpair] = intLumi; 

  }
  fclose(iff);

}


void Ntp1Analyzer::ReadJSONFile(const std::string& json) {


  std::cout << "Reading JSON file of good runs " << json << std::endl;
  FILE* iff = fopen(json.c_str(),"r");

  if(iff == 0) {
    std::cout << "cannot open JSON file " << json << " ... now exiting." << std::endl;
    exit(-1);
  }

  char c1, c2, c3;
  int run1, run2, LS1, LS2;

  std::cout << "Following LS will be used" << std::endl;
  std::cout << "-------------------------" << std::endl;
  while( fscanf(iff,"%*[ \t\n]%c%d:%d-%d:%d%c%c",&c1,&run1,&LS1,&run2,&LS2,&c2,&c3 ) != EOF ) {
      std::cout << "run: " << run1 << "  LS range: " << LS1
         << " --> " << LS2 << std::endl;
      goodLS_[run1].push_back(  std::pair<int,int>(LS1,LS2) );
  }
  fclose(iff);
  filterGoodRuns_ = true; // will run only on good runs/LS

}


bool Ntp1Analyzer::isGoodLS() {

     bool returnBool = false;

     if(!filterGoodRuns_) returnBool = true; // if filtered not requested all events are good

     this->UpdateCache();


     // check whether this run is part of the good runs. else retrun false
     if( goodLSCache_ != goodLS_.end() ) {

        // get list of LS intervals
        const GoodLSVector& lsvector =   goodLSCache_->second; 
        // loop over good LS intervals and return as soon as one interval contains this event
        for(GoodLSVector::const_iterator iLS = lsvector.begin(); (iLS != lsvector.end())&&(returnBool==false); iLS++) {
     
           if(lumiBlock >= iLS->first && lumiBlock <= iLS->second ) {
    
             returnBool = true;
          } // check current LS being in the interval
        } // loop over good LS for this run
     }

     if( returnBool==true ) {
   
       if( currentLS_ != lumiBlock )  {
          currentLS_ = lumiBlock;
          RunLumiPair rlpair = (std::pair<int, int>(runNumber,lumiBlock));
          totalIntLumi_ += LSLumimap_[rlpair];
       }
     }

     return returnBool;

}


GenEventParameters Ntp1Analyzer::getGenEventParameters() {

   GenEventParameters returnGenPars;
   returnGenPars.ptHatMin = -1.;

   if( dataset_=="PhotonJet_Summer09_Pt15" ) {
     returnGenPars.crossSection = 288813. - 32203.8;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="PhotonJet_Summer09_Pt30" ) {
     returnGenPars.crossSection = 32203.8 - 1012.08;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="PhotonJet_Summer09_Pt80" ) {
     returnGenPars.crossSection = 1012.08 - 51.36;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="PhotonJet_Summer09_Pt170" ) {
     returnGenPars.crossSection = 51.36 - 4.193;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="PhotonJet_Summer09_Pt300" ) {
     returnGenPars.crossSection =  4.193 - 0.45125;
     returnGenPars.ptHatMax = 470.;
   } else if( dataset_=="PhotonJet_Summer09_Pt470" ) {
     returnGenPars.crossSection =  0.45125 - 0.02;
     returnGenPars.ptHatMax = 800.;
   } else if( dataset_=="PhotonJet_Summer09_Pt800" ) {
     returnGenPars.crossSection = 0.02 - 0.000268;
     returnGenPars.ptHatMax = 1400.;
   } else if( dataset_=="PhotonJet_Summer09_Pt1400" ) {
     returnGenPars.crossSection = 0.000268;
     returnGenPars.ptHatMax = 14000.;
   } else if( dataset_=="QCD_Summer09_Pt15" ) {
     returnGenPars.crossSection = 1458126879.8;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="QCD_Summer09_Pt30" ) {
     returnGenPars.crossSection = 109005537.31;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="QCD_Summer09_Pt80" ) {
     returnGenPars.crossSection = 1936120.4893;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="QCD_Summer09_Pt170" ) {
     returnGenPars.crossSection = 62508.776856;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="QCD_Summer09_Pt300" ) {
     returnGenPars.crossSection = 3669.4197667;
     returnGenPars.ptHatMax = 470.;
   } else if( dataset_=="QCD_Summer09_Pt470" ) {
     returnGenPars.crossSection = 315.32221016;
     returnGenPars.ptHatMax = 800.;
   } else if( dataset_=="QCD_Summer09_Pt800" ) {
     returnGenPars.crossSection = 11.94070485;
     returnGenPars.ptHatMax = 1400.;
   } else if( dataset_=="QCD_Summer09_Pt1400" ) {
     returnGenPars.crossSection = 0.17207350709;
     returnGenPars.ptHatMax = 2200.;
   } else if( dataset_=="MC_PhotonJet_Summer09_Pt0to15" ) {
     returnGenPars.crossSection = 84460000.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="MC_PhotonJet_Summer09_Pt15to20" ) {
     returnGenPars.crossSection = 114700.;
     returnGenPars.ptHatMax = 20.;
   } else if( dataset_=="MC_PhotonJet_Summer09_Pt20to30" ) {
     returnGenPars.crossSection = 57180.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="MinBias_Spring10-START3X_V26A_356ReReco-v1"||dataset_=="MinBias_Spring10-START3X_V26A_357ReReco-v3"||dataset_=="MinBias_357ReReco_v3"||dataset_=="MinBias_357ReReco_v3_Pt0to15" ) {
     returnGenPars.crossSection = 71260000000.;
     returnGenPars.ptHatMax = ( dataset_=="MinBias_357ReReco_v3_Pt0to15" ) ? 15. : 7000.;
   } else if( dataset_=="PhotonJet_Spring10_Pt0to15" || dataset_ == "PhotonJet_Summer1036X_Pt0to15" ) {
     returnGenPars.crossSection = 84460000.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="PhotonJet_Spring10_Pt5to15" || dataset_=="PhotonJet_Summer1036X_Pt5to15" ) {
     returnGenPars.crossSection = 4030000.;
     returnGenPars.ptHatMax = 5.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="PhotonJet_Spring10_Pt15" ) {
     returnGenPars.crossSection = 192200.-20070.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="PhotonJet_Spring10_Pt30" ) {
     returnGenPars.crossSection = 20070.-556.5;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="PhotonJet_Spring10_Pt80" ) {
     returnGenPars.crossSection = 556.5-24.37;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="PhotonJet_Spring10_Pt170" ) {
     returnGenPars.crossSection = 24.37-1.636;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="PhotonJet_Spring10_Pt300" ) {
     returnGenPars.crossSection = 1.636-0.136;
     returnGenPars.ptHatMax = 470.;
   } else if( dataset_=="PhotonJet_Spring10_Pt470" ) {
     returnGenPars.crossSection = 0.136-0.003477;
     returnGenPars.ptHatMax = 800.;
   } else if( dataset_=="PhotonJet_Spring10_Pt800" ) {
     returnGenPars.crossSection = 0.003477-0.00001286;
     returnGenPars.ptHatMax = 1400.;
   } else if( dataset_=="PhotonJet_Spring10_Pt1400" ) {
     returnGenPars.crossSection = 0.00001286;
     returnGenPars.ptHatMax = 7000.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt15to20" ) {
     returnGenPars.crossSection = 114700.;
     returnGenPars.ptHatMin = 15.;
     returnGenPars.ptHatMax = 20.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt20to30" ) {
     returnGenPars.crossSection = 57180.;
     returnGenPars.ptHatMin = 20.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt30to50" ) {
     returnGenPars.crossSection = 16520.;
     returnGenPars.ptHatMin = 30.;
     returnGenPars.ptHatMax = 50.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt50to80" ) {
     returnGenPars.crossSection = 2723.;
     returnGenPars.ptHatMin = 50.;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt80to120" ) {
     returnGenPars.crossSection = 446.2;
     returnGenPars.ptHatMin = 80.;
     returnGenPars.ptHatMax = 120.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt120to170" ) {
     returnGenPars.crossSection = 84.43;
     returnGenPars.ptHatMin = 120.;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt170to300" ) {
     returnGenPars.crossSection = 22.55;
     returnGenPars.ptHatMin = 170.;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="QCD_Spring10_Pt0to15" ) {
     returnGenPars.crossSection = 48445000000.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="QCD_Spring10_Pt5to15" ) {
     returnGenPars.crossSection = 36640000000.;
     returnGenPars.ptHatMin = 5.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="QCD_Spring10_Pt15" ) {
     returnGenPars.crossSection = 876215000.-60411000.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="QCD_Spring10_Pt15to20" ) {
     returnGenPars.crossSection = 579411000.;
     returnGenPars.ptHatMin = 15.;
     returnGenPars.ptHatMax = 20.;
   } else if( dataset_=="QCD_Spring10_Pt20to30" ) {
     returnGenPars.crossSection = 236051000.;
     returnGenPars.ptHatMin = 20.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="QCD_Spring10_Pt30" ) {
     returnGenPars.crossSection = 60411000.-923821.;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="QCD_Spring10_Pt30to50" ) {
     returnGenPars.crossSection = 53114800.;
     returnGenPars.ptHatMin = 30.;
     returnGenPars.ptHatMax = 50.;
   } else if( dataset_=="QCD_Spring10_Pt50to80" ) {
     returnGenPars.crossSection = 6358210.;
     returnGenPars.ptHatMin = 50.;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="QCD_Spring10_Pt80" ) {
     returnGenPars.crossSection = 923821.-25474.9;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="QCD_Spring10_Pt170" ) {
     returnGenPars.crossSection = 25474.9-1255.87;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="QCD_Spring10_Pt300" ) {
     returnGenPars.crossSection = 1255.87-87.9799;
     returnGenPars.ptHatMax = 470.;
   } else if( dataset_=="QCD_Spring10_Pt470" ) {
     returnGenPars.crossSection = 87.9799-2.18608;
     returnGenPars.ptHatMax = 800.;
   } else if( dataset_=="QCD_Spring10_Pt800" ) {
     returnGenPars.crossSection = 2.18608-0.0112233;
     returnGenPars.ptHatMax = 1400.;
   } else if( dataset_=="QCD_Spring10_Pt1400" ) {
     returnGenPars.crossSection = 0.0112233;
     returnGenPars.ptHatMax = 10000.;
   } else if( dataset_=="Wenu_Summer10_START37_V5_S09_v1" ) {
     returnGenPars.crossSection = 7899.;
     returnGenPars.ptHatMax = 10000.;
   } else {
     std::cout << "-> (No ptHat cuts introduced.)" << std::endl;
     returnGenPars.crossSection = -1.;
     returnGenPars.ptHatMin = 0.;
     returnGenPars.ptHatMax = 10000.;
   }

//   if( returnGenPars.crossSection != -1 ) 
//     std::cout << "-> Dataset was in database. Cross-section correctly set." << std::endl;

   return returnGenPars;

}
