#include "Ntp1Analyzer.h"
#include "TH1F.h"
#include "TRandom.h"
#include <cstdlib>
#include <fstream>
#include <cmath>



Ntp1Analyzer::Ntp1Analyzer(const std::string& analyzerType, const std::string& dataset, const std::string& flags, TTree* tree)
{

   dataset_ = dataset;

   DEBUG_ = false;
   filterGoodRuns_ = false; //default: do not filter
   totalIntLumi_ = 0.;

   analyzerType_ = analyzerType;

   flags_ = flags;

   ptHatMin_ = 0.;
   ptHatMax_ = 10000.;

   cachedLS_ = 0;
   cachedRun_ = 0;

   rand_ = new TRandom();

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

   std::string treeDir;
   char treePath[400];
   TChain * chain = new TChain("ntp1","");
// if( dataset_=="Wenu_Summer10_START37_V5_S09_v1" ) {
//   treeDir = "/cmsrm/pc21_2/pandolf/MC/Wenu_Summer10_START37_V5_S09_v1";
// } else if( dataset_=="HZZ_qqll_gluonfusion_M200" ) {
//   treeDir = "/cmsrm/pc21_2/pandolf/MC/HZZ_qqll_gluonfusion_M200";
// } else if( dataset_=="HZZ_qqll_gluonfusion_M300" ) {
//   treeDir = "/cmsrm/pc21_2/pandolf/MC/HZZ_qqll_gluonfusion_M300";
// } else if( dataset_=="HZZ_qqll_gluonfusion_M400" ) {
//   treeDir = "/cmsrm/pc21_2/pandolf/MC/HZZ_qqll_gluonfusion_M400";
// } else {
// }


   sprintf(treePath, "%s/default_*.root/ntp1", treeDir.c_str());

   if( dataset_=="HZZ_qqll_gluonfusion_M300_CANDS") {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/HZZ_qqll_gluonfusion_M300/default_CANDS_1000ev.root");
   } else if( dataset_=="HZZ_qqll_gluonfusion_M400_CANDS") {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/HZZ_qqll_gluonfusion_M400/default_CANDS_1000ev.root");
   } 


   int addInt = chain->Add(treePath);

   if( addInt==0 ) {
     std::cout << "Didn't find files to add for dataset: '" << dataset_ << "'. Looking for a list..." << std::endl;
     std::string fileName = "files_" + dataset_ + ".txt";
     this->LoadInputFromFile(fileName);
   } 
   /*else {
     TTree* tree = chain;
     std::cout << "-> Tree has " << tree->GetEntries() << " entries." << std::endl;
     this->CreateOutputFile();
     Init(tree);
     //load trigger mask:
     std::string firstFileName = treeDir + "/default_1.root";
     if( dataset_=="HZZ_qqll_gluonfusion_M300_CANDS" || dataset_=="HZZ_qqll_gluonfusion_M400_CANDS" )
       firstFileName = treePath;
     TFile* firstFile = TFile::Open( firstFileName.c_str(), "read" );
     this->LoadTrigger(firstFile);
   }*/

}


void Ntp1Analyzer::LoadInputFromFile( const std::string& fileName ) {

   FILE* iff = fopen(fileName.c_str(),"r");
   if(iff == 0) {
     std::cout << "cannot open input file " << fileName << " ... now exiting." << std::endl;
     exit(-1);
   }

   TChain * chain = new TChain("ntp1","");

   char singleLine[500];
   bool isFirstFile=true;

   TFile* firstFile = 0;

   while( fscanf(iff, "%s", singleLine) !=EOF ) {
   
     std::string singleLine_str(singleLine);
     std::string treeName_str = singleLine_str + "/ntp1";
     std::cout << "-> Adding " << treeName_str << std::endl;
     chain->Add(treeName_str.c_str());
     if( isFirstFile ) {
       firstFile = TFile::Open(singleLine_str.c_str(), "read");
       isFirstFile=false;
     }

   }
   fclose(iff);

   TTree* tree = chain;
   std::cout << "-> Tree has " << tree->GetEntries() << " entries." << std::endl;
   this->CreateOutputFile();
   Init(tree);
   //this->LoadTrigger(firstFile);

}




void Ntp1Analyzer::LoadTrigger( int iEntry, bool verbose, TFile* condFile ) {

  
  TTree* treeCond = (condFile==0) ? 0 : (TTree*)(condFile->Get("Conditions"));

  std::vector<std::string> foundTriggers;
  std::vector<std::string> foundTriggersNOT;

  //new version: trigger loaded directly from ntp1 tree:
  if( treeCond==0 ) { 

    fChain->GetEntry(iEntry);
    //fChain->GetEntry(0);

    // required triggers:
    std::vector<int> triggerMask_required;
    for (std::vector< std::string >::const_iterator fIter=requiredTriggers_.begin();fIter!=requiredTriggers_.end();++fIter)
      {
        bool foundThisTrigger = false;
//std::cout << "looking for: " << (*fIter) << std::endl;
        for(unsigned int i=0; i<nameHLT->size(); i++) 
          {
//std::cout << std::endl << indexHLT[i] << " " << nameHLT->at(i);
            TString nameHLT_tstr(nameHLT->at(i));
            if( nameHLT_tstr.Contains((*fIter)) )
            //if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) ) 
              {
//std::cout << " <----- HERE IT IS!" << std::endl;
                foundThisTrigger = true;
                triggerMask_required.push_back( indexHLT[i] ) ;
                foundTriggers.push_back( nameHLT->at(i) ) ;
                break;
              }
          }
          if( !foundThisTrigger && verbose ) std::cout << "-> WARNING!! Didn't find HLT path: " << (*fIter).c_str() << ". Ignoring it." << std::endl;
      }
    index_requiredTriggers_ = triggerMask_required;

    // NOT triggers
    std::vector<int> triggerMask_NOT;
    for (std::vector< std::string >::const_iterator fIter=notTriggers_.begin();fIter!=notTriggers_.end();++fIter)
      {
        bool foundThisTrigger = false;
        for(unsigned int i=0; i<nameHLT->size(); i++) 
          {
//std::cout << std::endl << nameHLT->at(i);
            TString nameHLT_tstr(nameHLT->at(i));
            if( nameHLT_tstr.Contains((*fIter)) )
            //if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) ) 
              {
//std::cout << " <----- HERE IT IS!" << std::endl;
                foundThisTrigger = true;
                triggerMask_NOT.push_back( indexHLT[i] ) ;
                foundTriggersNOT.push_back( nameHLT->at(i) ) ;
                break;
              }
          }
          if( !foundThisTrigger && verbose ) std::cout << "-> WARNING!! Didn't find HLT path: " << (*fIter).c_str() << ". Ignoring it." << std::endl;
      }
    index_notTriggers_ = triggerMask_NOT;

  } else { //old version: Conditions

    int           nHLT_;
    std::vector<std::string>  *nameHLT_;
    std::vector<unsigned int> *indexHLT_;

    //To get the pointers for the vectors
    nameHLT_=0;
    indexHLT_=0;

    treeCond->SetBranchAddress("nHLT", &nHLT_);
    treeCond->SetBranchAddress("nameHLT", &nameHLT_);
    treeCond->SetBranchAddress("indexHLT", &indexHLT_);
    treeCond->GetEntry(0);

    std::vector<int> triggerMask;
    for (std::vector< std::string >::const_iterator fIter=requiredTriggers_.begin();fIter!=requiredTriggers_.end();++fIter)
      {
        bool foundThisTrigger = false;
        for(unsigned int i=0; i<nameHLT_->size(); i++) 
          {
            TString nameHLT_tstr(nameHLT_->at(i));
            //if( !strcmp ((*fIter).c_str(), nameHLT_->at(i).c_str() ) ) 
            if( nameHLT_tstr.Contains((*fIter).c_str()) ) 
              {
                foundThisTrigger = true;
                triggerMask.push_back( indexHLT_->at(i) ) ;
                break;
              }
          }
          if( !foundThisTrigger && verbose ) std::cout << "-> WARNING!! Didn't find HLT path: " << (*fIter).c_str() << ". Ignoring it." << std::endl;
      }
    index_requiredTriggers_ = triggerMask;

  }


  if( requiredTriggers_.size()==0 && notTriggers_.size()==0 && verbose )
    std::cout << "-> No trigger selection required." << std::endl;

  for (int i=0;i<index_requiredTriggers_.size();++i)
    if( verbose ) std::cout << "[ReloadTriggerMask]::Requiring bit " << index_requiredTriggers_[i] << " " << foundTriggers[i] << std::endl;

  for (int i=0;i<index_notTriggers_.size();++i)
    if( verbose ) std::cout << "[ReloadTriggerMask]::Vetoing bit " << index_notTriggers_[i] << " " << foundTriggersNOT[i] << std::endl;


} // LoadTrigger




bool Ntp1Analyzer::PassedHLT( int iEntry, const std::string& HLTName ) { //default is OR of all required triggers (HLTName=="")


  if ( index_requiredTriggers_.size()==0 && index_notTriggers_.size()==0 && HLTName=="" ) return true;


  bool rememberToReset = false;
  std::vector<int> index_requiredTriggers_tmp = index_requiredTriggers_;
  std::vector<std::string> requiredTriggers_tmp = requiredTriggers_;
  if( HLTName!="" ) {
    index_requiredTriggers_.clear();
    requiredTriggers_.clear();
    requiredTriggers_.push_back(HLTName);
    this->LoadTrigger(iEntry, false);
    rememberToReset = true;
  }



  // first NOT triggers:
  for( int i=0; i<index_notTriggers_.size(); i++ ) {

    int block_veto =  index_notTriggers_[i]/30;
    int pos_veto = index_notTriggers_[i]%30;
    int word_veto = firedTrg[block_veto];
    
    if( (word_veto >> pos_veto)%2 ) {
      if( rememberToReset ) {
        index_requiredTriggers_ = index_requiredTriggers_tmp;
        requiredTriggers_ = requiredTriggers_tmp;
      }
      return false;
    }

  } // for not triggers
  

  // now required triggers:
  for( int i=0; i<index_requiredTriggers_.size(); i++ ) {

    if( HLTName=="" || requiredTriggers_[i]==HLTName ) {

      int block_required =  index_requiredTriggers_[i]/30;
      int pos_required = index_requiredTriggers_[i]%30;
      int word_required = firedTrg[block_required];

      if( (word_required >> pos_required)%2 ) {
        if( rememberToReset ) {
          index_requiredTriggers_ = index_requiredTriggers_tmp;
          requiredTriggers_ = requiredTriggers_tmp;
        }
        return true;
      }

    } // if required

  } // required trigger loop


  if( rememberToReset ) {
    index_requiredTriggers_ = index_requiredTriggers_tmp;
    requiredTriggers_ = requiredTriggers_tmp;
  }


  return false;

}




void Ntp1Analyzer::CreateOutputFile() {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova2ndLevel_"+dataset_;
   else {
    if(dataset_!="") outfileName = analyzerType_ + "_2ndLevelTree_"+dataset_;
    else outfileName = analyzerType_ + "_2ndLevelTree";
   }


   if( flags_!="" )
     outfileName = outfileName + "_" + flags_;
   outfileName = outfileName + ".root";

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

   // Set object pointer
   nameHLT = 0;
   // Set branch addresses and branch pointers
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->GetEntry(0);
   isMC_ = (runNumber < 10);

   GenEventParameters genPars = this->getGenEventParameters();
   ptHatMin_ = genPars.ptHatMin;
   ptHatMax_ = genPars.ptHatMax;

   //will cut on pt_hat, so have to divide only by correct number of events:
   char cutOnPtHat[70];
   sprintf( cutOnPtHat, "genPtHat>%lf && genPtHat<%lf", (Double_t)ptHatMin_, (Double_t)ptHatMax_);
   Int_t nEntries_cut = (isMC_) ? fChain->GetEntries(cutOnPtHat) : fChain->GetEntries();
   h1_nCounter_->SetBinContent( 1, nEntries_cut );


   std::string branchName;

   fChain->SetBranchAddress("nl1Technical", &nl1Technical, &b_nl1Technical);
   fChain->SetBranchAddress("l1Technical", l1Technical, &b_l1Technical);
   fChain->SetBranchAddress("nl1Global", &nl1Global, &b_nl1Global);
   fChain->SetBranchAddress("l1Global", l1Global, &b_l1Global);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("rhoFastjet", &rhoFastjet, &b_rhoFastjet);
   if( isMC_ ) {
     fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
     fChain->SetBranchAddress("nPU", nPU, &b_nPU);
     fChain->SetBranchAddress("bxPU", bxPU, &b_bxPU);
     fChain->SetBranchAddress("nMc", &nMc, &b_nMc);
     fChain->SetBranchAddress("pMc", pMc, &b_pMc);
     fChain->SetBranchAddress("thetaMc", thetaMc, &b_thetaMc);
     fChain->SetBranchAddress("etaMc", etaMc, &b_etaMc);
     fChain->SetBranchAddress("phiMc", phiMc, &b_phiMc);
     fChain->SetBranchAddress("energyMc", energyMc, &b_energyMc);
     fChain->SetBranchAddress("idMc", idMc, &b_idMc);
     fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
     fChain->SetBranchAddress("statusMc", statusMc, &b_statusMc);
   }
   fChain->SetBranchAddress("nTrg", &nTrg, &b_nTrg);
   fChain->SetBranchAddress("firedTrg", firedTrg, &b_firedTrg);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("nameHLT", &nameHLT, &b_nameHLT);
   fChain->SetBranchAddress("indexHLT", indexHLT, &b_indexHLT);
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
   //fChain->SetBranchAddress("esEnergyEle", esEnergyEle, &b_esEnergyEle);
   fChain->SetBranchAddress("superClusterIndexEle", superClusterIndexEle, &b_superClusterIndexEle);
   fChain->SetBranchAddress("PFsuperClusterIndexEle", PFsuperClusterIndexEle, &b_PFsuperClusterIndexEle);
   fChain->SetBranchAddress("trackIndexEle", trackIndexEle, &b_trackIndexEle);
   fChain->SetBranchAddress("gsfTrackIndexEle", gsfTrackIndexEle, &b_gsfTrackIndexEle);
   fChain->SetBranchAddress("convDistEle", convDistEle, &b_convDistEle);
   fChain->SetBranchAddress("convDcotEle", convDcotEle, &b_convDcotEle);
   fChain->SetBranchAddress("convRadiusEle", convRadiusEle, &b_convRadiusEle);
   fChain->SetBranchAddress("convTrackIndexEle", convTrackIndexEle, &b_convTrackIndexEle);
   //fChain->SetBranchAddress("convXEle", convXEle, &b_convXEle);
   //fChain->SetBranchAddress("convYEle", convYEle, &b_convYEle);
   //fChain->SetBranchAddress("convZEle", convZEle, &b_convZEle);
   //fChain->SetBranchAddress("convChi2ProbEle", convChi2ProbEle, &b_convChi2ProbEle);
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
   //fChain->SetBranchAddress("tipEle", tipEle, &b_tipEle);
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
   fChain->SetBranchAddress("sMajSC", sMajSC, &b_sMajSC);
   fChain->SetBranchAddress("sMinSC", sMinSC, &b_sMinSC);
   fChain->SetBranchAddress("alphaSC", alphaSC, &b_alphaSC);
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
   //fChain->SetBranchAddress("nBC", &nBC, &b_nBC);
   //fChain->SetBranchAddress("nCrystalsBC", nCrystalsBC, &b_nCrystalsBC);
   //fChain->SetBranchAddress("energyBC", energyBC, &b_energyBC);
   //fChain->SetBranchAddress("etaBC", etaBC, &b_etaBC);
   //fChain->SetBranchAddress("thetaBC", thetaBC, &b_thetaBC);
   //fChain->SetBranchAddress("phiBC", phiBC, &b_phiBC);
   //fChain->SetBranchAddress("e3x3BC", e3x3BC, &b_e3x3BC);
   //fChain->SetBranchAddress("e5x5BC", e5x5BC, &b_e5x5BC);
   //fChain->SetBranchAddress("eMaxBC", eMaxBC, &b_eMaxBC);
   //fChain->SetBranchAddress("e2x2BC", e2x2BC, &b_e2x2BC);
   //fChain->SetBranchAddress("e2ndBC", e2ndBC, &b_e2ndBC);
   //fChain->SetBranchAddress("covIEtaIEtaBC", covIEtaIEtaBC, &b_covIEtaIEtaBC);
   //fChain->SetBranchAddress("covIEtaIPhiBC", covIEtaIPhiBC, &b_covIEtaIPhiBC);
   //fChain->SetBranchAddress("covIPhiIPhiBC", covIPhiIPhiBC, &b_covIPhiIPhiBC);
   //fChain->SetBranchAddress("recoFlagBC", recoFlagBC, &b_recoFlagBC);
   //fChain->SetBranchAddress("timeBC", timeBC, &b_timeBC);
   //fChain->SetBranchAddress("chi2BC", chi2BC, &b_chi2BC);
   //fChain->SetBranchAddress("seedEnergyBC", seedEnergyBC, &b_seedEnergyBC);
   //fChain->SetBranchAddress("idClosProblBC", idClosProblBC, &b_idClosProblBC);
   //fChain->SetBranchAddress("sevClosProblBC", sevClosProblBC, &b_sevClosProblBC);
   //fChain->SetBranchAddress("fracClosProblBC", fracClosProblBC, &b_fracClosProblBC);
   //fChain->SetBranchAddress("indexSCBC", indexSCBC, &b_indexSCBC);
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
   //fChain->SetBranchAddress("pxAtOuterTrack", pxAtOuterTrack, &b_pxAtOuterTrack);
   //fChain->SetBranchAddress("pyAtOuterTrack", pyAtOuterTrack, &b_pyAtOuterTrack);
   //fChain->SetBranchAddress("pzAtOuterTrack", pzAtOuterTrack, &b_pzAtOuterTrack);
   //fChain->SetBranchAddress("xAtOuterTrack", xAtOuterTrack, &b_xAtOuterTrack);
   //fChain->SetBranchAddress("yAtOuterTrack", yAtOuterTrack, &b_yAtOuterTrack);
   //fChain->SetBranchAddress("zAtOuterTrack", zAtOuterTrack, &b_zAtOuterTrack);
   //fChain->SetBranchAddress("pxAtInnerTrack", pxAtInnerTrack, &b_pxAtInnerTrack);
   //fChain->SetBranchAddress("pyAtInnerTrack", pyAtInnerTrack, &b_pyAtInnerTrack);
   //fChain->SetBranchAddress("pzAtInnerTrack", pzAtInnerTrack, &b_pzAtInnerTrack);
   //fChain->SetBranchAddress("xAtInnerTrack", xAtInnerTrack, &b_xAtInnerTrack);
   //fChain->SetBranchAddress("yAtInnerTrack", yAtInnerTrack, &b_yAtInnerTrack);
   //fChain->SetBranchAddress("zAtInnerTrack", zAtInnerTrack, &b_zAtInnerTrack);
   //fChain->SetBranchAddress("recHitsSizeTrack", recHitsSizeTrack, &b_recHitsSizeTrack);
   fChain->SetBranchAddress("pixelHitsTrack", pixelHitsTrack, &b_pixelHitsTrack);
   fChain->SetBranchAddress("expInnerLayersTrack", expInnerLayersTrack, &b_expInnerLayersTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsTrack", numberOfValidPixelBarrelHitsTrack, &b_numberOfValidPixelBarrelHitsTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsTrack", numberOfValidPixelEndcapHitsTrack, &b_numberOfValidPixelEndcapHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsTrack", numberOfValidStripTIBHitsTrack, &b_numberOfValidStripTIBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsTrack", numberOfValidStripTIDHitsTrack, &b_numberOfValidStripTIDHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsTrack", numberOfValidStripTOBHitsTrack, &b_numberOfValidStripTOBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsTrack", numberOfValidStripTECHitsTrack, &b_numberOfValidStripTECHitsTrack);
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
   //fChain->SetBranchAddress("pxAtOuterGsfTrack", pxAtOuterGsfTrack, &b_pxAtOuterGsfTrack);
   //fChain->SetBranchAddress("pyAtOuterGsfTrack", pyAtOuterGsfTrack, &b_pyAtOuterGsfTrack);
   //fChain->SetBranchAddress("pzAtOuterGsfTrack", pzAtOuterGsfTrack, &b_pzAtOuterGsfTrack);
   //fChain->SetBranchAddress("xAtOuterGsfTrack", xAtOuterGsfTrack, &b_xAtOuterGsfTrack);
   //fChain->SetBranchAddress("yAtOuterGsfTrack", yAtOuterGsfTrack, &b_yAtOuterGsfTrack);
   //fChain->SetBranchAddress("zAtOuterGsfTrack", zAtOuterGsfTrack, &b_zAtOuterGsfTrack);
   //fChain->SetBranchAddress("pxAtInnerGsfTrack", pxAtInnerGsfTrack, &b_pxAtInnerGsfTrack);
   //fChain->SetBranchAddress("pyAtInnerGsfTrack", pyAtInnerGsfTrack, &b_pyAtInnerGsfTrack);
   //fChain->SetBranchAddress("pzAtInnerGsfTrack", pzAtInnerGsfTrack, &b_pzAtInnerGsfTrack);
   //fChain->SetBranchAddress("xAtInnerGsfTrack", xAtInnerGsfTrack, &b_xAtInnerGsfTrack);
   //fChain->SetBranchAddress("yAtInnerGsfTrack", yAtInnerGsfTrack, &b_yAtInnerGsfTrack);
   //fChain->SetBranchAddress("zAtInnerGsfTrack", zAtInnerGsfTrack, &b_zAtInnerGsfTrack);
   //fChain->SetBranchAddress("recHitsSizeGsfTrack", recHitsSizeGsfTrack, &b_recHitsSizeGsfTrack);
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
   //fChain->SetBranchAddress("recoFlagsGsfTrack", recoFlagsGsfTrack, &b_recoFlagsGsfTrack);
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
   //fChain->SetBranchAddress("pxAtOuterGlobalMuonTrack", pxAtOuterGlobalMuonTrack, &b_pxAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("pyAtOuterGlobalMuonTrack", pyAtOuterGlobalMuonTrack, &b_pyAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("pzAtOuterGlobalMuonTrack", pzAtOuterGlobalMuonTrack, &b_pzAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("xAtOuterGlobalMuonTrack", xAtOuterGlobalMuonTrack, &b_xAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("yAtOuterGlobalMuonTrack", yAtOuterGlobalMuonTrack, &b_yAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("zAtOuterGlobalMuonTrack", zAtOuterGlobalMuonTrack, &b_zAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("pxAtInnerGlobalMuonTrack", pxAtInnerGlobalMuonTrack, &b_pxAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("pyAtInnerGlobalMuonTrack", pyAtInnerGlobalMuonTrack, &b_pyAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("pzAtInnerGlobalMuonTrack", pzAtInnerGlobalMuonTrack, &b_pzAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("xAtInnerGlobalMuonTrack", xAtInnerGlobalMuonTrack, &b_xAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("yAtInnerGlobalMuonTrack", yAtInnerGlobalMuonTrack, &b_yAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("zAtInnerGlobalMuonTrack", zAtInnerGlobalMuonTrack, &b_zAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("recHitsSizeGlobalMuonTrack", recHitsSizeGlobalMuonTrack, &b_recHitsSizeGlobalMuonTrack);
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
   //fChain->SetBranchAddress("vtxIndexSTAMuonTrack", vtxIndexSTAMuonTrack, &b_vtxIndexSTAMuonTrack);
   //fChain->SetBranchAddress("vtxWeightSTAMuonTrack", vtxWeightSTAMuonTrack, &b_vtxWeightSTAMuonTrack);
   fChain->SetBranchAddress("chargeSTAMuonTrack", chargeSTAMuonTrack, &b_chargeSTAMuonTrack);
   fChain->SetBranchAddress("ptErrorSTAMuonTrack", ptErrorSTAMuonTrack, &b_ptErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackValidHitsSTAMuonTrack", trackValidHitsSTAMuonTrack, &b_trackValidHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackLostHitsSTAMuonTrack", trackLostHitsSTAMuonTrack, &b_trackLostHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2STAMuonTrack", trackNormalizedChi2STAMuonTrack, &b_trackNormalizedChi2STAMuonTrack);
   fChain->SetBranchAddress("qualityMaskSTAMuonTrack", qualityMaskSTAMuonTrack, &b_qualityMaskSTAMuonTrack);
   //fChain->SetBranchAddress("impactPar3DSTAMuonTrack", impactPar3DSTAMuonTrack, &b_impactPar3DSTAMuonTrack);
   //fChain->SetBranchAddress("impactPar3DErrorSTAMuonTrack", impactPar3DErrorSTAMuonTrack, &b_impactPar3DErrorSTAMuonTrack);
   //fChain->SetBranchAddress("transvImpactParSTAMuonTrack", transvImpactParSTAMuonTrack, &b_transvImpactParSTAMuonTrack);
   //fChain->SetBranchAddress("transvImpactParErrorSTAMuonTrack", transvImpactParErrorSTAMuonTrack, &b_transvImpactParErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackVxSTAMuonTrack", trackVxSTAMuonTrack, &b_trackVxSTAMuonTrack);
   fChain->SetBranchAddress("trackVySTAMuonTrack", trackVySTAMuonTrack, &b_trackVySTAMuonTrack);
   fChain->SetBranchAddress("trackVzSTAMuonTrack", trackVzSTAMuonTrack, &b_trackVzSTAMuonTrack);
   //fChain->SetBranchAddress("pxAtOuterSTAMuonTrack", pxAtOuterSTAMuonTrack, &b_pxAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("pyAtOuterSTAMuonTrack", pyAtOuterSTAMuonTrack, &b_pyAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("pzAtOuterSTAMuonTrack", pzAtOuterSTAMuonTrack, &b_pzAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("xAtOuterSTAMuonTrack", xAtOuterSTAMuonTrack, &b_xAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("yAtOuterSTAMuonTrack", yAtOuterSTAMuonTrack, &b_yAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("zAtOuterSTAMuonTrack", zAtOuterSTAMuonTrack, &b_zAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("pxAtInnerSTAMuonTrack", pxAtInnerSTAMuonTrack, &b_pxAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("pyAtInnerSTAMuonTrack", pyAtInnerSTAMuonTrack, &b_pyAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("pzAtInnerSTAMuonTrack", pzAtInnerSTAMuonTrack, &b_pzAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("xAtInnerSTAMuonTrack", xAtInnerSTAMuonTrack, &b_xAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("yAtInnerSTAMuonTrack", yAtInnerSTAMuonTrack, &b_yAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("zAtInnerSTAMuonTrack", zAtInnerSTAMuonTrack, &b_zAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("recHitsSizeSTAMuonTrack", recHitsSizeSTAMuonTrack, &b_recHitsSizeSTAMuonTrack);
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
   fChain->SetBranchAddress("numberOfMatchesMuon", numberOfMatchesMuon, &b_numberOfMatchesMuon);
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
   fChain->SetBranchAddress("significancePFMet", significancePFMet, &b_significancePFMet);
   fChain->SetBranchAddress("mEtSigPFMet", mEtSigPFMet, &b_mEtSigPFMet);
   fChain->SetBranchAddress("sumEtPFMet", sumEtPFMet, &b_sumEtPFMet);
   if( isMC_ ) {
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
   }
 //fChain->SetBranchAddress("nPFCand", &nPFCand, &b_nPFCand);
 //fChain->SetBranchAddress("chargePFCand", chargePFCand, &b_chargePFCand);
 //fChain->SetBranchAddress("energyPFCand", energyPFCand, &b_energyPFCand);
 //fChain->SetBranchAddress("thetaPFCand", thetaPFCand, &b_thetaPFCand);
 //fChain->SetBranchAddress("etaPFCand", etaPFCand, &b_etaPFCand);
 //fChain->SetBranchAddress("phiPFCand", phiPFCand, &b_phiPFCand);
 //fChain->SetBranchAddress("pxPFCand", pxPFCand, &b_pxPFCand);
 //fChain->SetBranchAddress("pyPFCand", pyPFCand, &b_pyPFCand);
 //fChain->SetBranchAddress("pzPFCand", pzPFCand, &b_pzPFCand);
 //fChain->SetBranchAddress("vertexXPFCand", vertexXPFCand, &b_vertexXPFCand);
 //fChain->SetBranchAddress("vertexYPFCand", vertexYPFCand, &b_vertexYPFCand);
 //fChain->SetBranchAddress("vertexZPFCand", vertexZPFCand, &b_vertexZPFCand);
 //fChain->SetBranchAddress("particleTypePFCand", particleTypePFCand, &b_particleTypePFCand);
 //fChain->SetBranchAddress("iPFJetPFCand", iPFJetPFCand, &b_iPFJetPFCand);
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
   //fChain->SetBranchAddress("simpleSecondaryVertexBJetTagsAK5Jet", simpleSecondaryVertexBJetTagsAK5Jet, &b_simpleSecondaryVertexBJetTagsAK5Jet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5Jet", softMuonBJetTagsAK5Jet, &b_softMuonBJetTagsAK5Jet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5Jet", trackCountingHighPurBJetTagsAK5Jet, &b_trackCountingHighPurBJetTagsAK5Jet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5Jet", trackCountingHighEffBJetTagsAK5Jet, &b_trackCountingHighEffBJetTagsAK5Jet);
   fChain->SetBranchAddress("uncorrEnergyAK5Jet", uncorrEnergyAK5Jet, &b_uncorrEnergyAK5Jet);
   fChain->SetBranchAddress("nAK5PFPUcorrJet", &nAK5PFPUcorrJet, &b_nAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargeAK5PFPUcorrJet", chargeAK5PFPUcorrJet, &b_chargeAK5PFPUcorrJet);
   fChain->SetBranchAddress("energyAK5PFPUcorrJet", energyAK5PFPUcorrJet, &b_energyAK5PFPUcorrJet);
   fChain->SetBranchAddress("thetaAK5PFPUcorrJet", thetaAK5PFPUcorrJet, &b_thetaAK5PFPUcorrJet);
   fChain->SetBranchAddress("etaAK5PFPUcorrJet", etaAK5PFPUcorrJet, &b_etaAK5PFPUcorrJet);
   fChain->SetBranchAddress("phiAK5PFPUcorrJet", phiAK5PFPUcorrJet, &b_phiAK5PFPUcorrJet);
   fChain->SetBranchAddress("pxAK5PFPUcorrJet", pxAK5PFPUcorrJet, &b_pxAK5PFPUcorrJet);
   fChain->SetBranchAddress("pyAK5PFPUcorrJet", pyAK5PFPUcorrJet, &b_pyAK5PFPUcorrJet);
   fChain->SetBranchAddress("pzAK5PFPUcorrJet", pzAK5PFPUcorrJet, &b_pzAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexXAK5PFPUcorrJet", vertexXAK5PFPUcorrJet, &b_vertexXAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexYAK5PFPUcorrJet", vertexYAK5PFPUcorrJet, &b_vertexYAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexZAK5PFPUcorrJet", vertexZAK5PFPUcorrJet, &b_vertexZAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedHadronEnergyAK5PFPUcorrJet", chargedHadronEnergyAK5PFPUcorrJet, &b_chargedHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralHadronEnergyAK5PFPUcorrJet", neutralHadronEnergyAK5PFPUcorrJet, &b_neutralHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("photonEnergyAK5PFPUcorrJet", photonEnergyAK5PFPUcorrJet, &b_photonEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("electronEnergyAK5PFPUcorrJet", electronEnergyAK5PFPUcorrJet, &b_electronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("muonEnergyAK5PFPUcorrJet", muonEnergyAK5PFPUcorrJet, &b_muonEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFHadronEnergyAK5PFPUcorrJet", HFHadronEnergyAK5PFPUcorrJet, &b_HFHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFEMEnergyAK5PFPUcorrJet", HFEMEnergyAK5PFPUcorrJet, &b_HFEMEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedHadronMultiplicityAK5PFPUcorrJet", chargedHadronMultiplicityAK5PFPUcorrJet, &b_chargedHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralHadronMultiplicityAK5PFPUcorrJet", neutralHadronMultiplicityAK5PFPUcorrJet, &b_neutralHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("photonMultiplicityAK5PFPUcorrJet", photonMultiplicityAK5PFPUcorrJet, &b_photonMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("electronMultiplicityAK5PFPUcorrJet", electronMultiplicityAK5PFPUcorrJet, &b_electronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("muonMultiplicityAK5PFPUcorrJet", muonMultiplicityAK5PFPUcorrJet, &b_muonMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFHadronMultiplicityAK5PFPUcorrJet", HFHadronMultiplicityAK5PFPUcorrJet, &b_HFHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFEMMultiplicityAK5PFPUcorrJet", HFEMMultiplicityAK5PFPUcorrJet, &b_HFEMMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedEmEnergyAK5PFPUcorrJet", chargedEmEnergyAK5PFPUcorrJet, &b_chargedEmEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralEmEnergyAK5PFPUcorrJet", neutralEmEnergyAK5PFPUcorrJet, &b_neutralEmEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5PFPUcorrJet", combinedSecondaryVertexBJetTagsAK5PFPUcorrJet, &b_combinedSecondaryVertexBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet", combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet, &b_combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5PFPUcorrJet", jetBProbabilityBJetTagsAK5PFPUcorrJet, &b_jetBProbabilityBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5PFPUcorrJet", jetProbabilityBJetTagsAK5PFPUcorrJet, &b_jetProbabilityBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet", simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet, &b_simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet", simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet, &b_simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5PFPUcorrJet", softMuonBJetTagsAK5PFPUcorrJet, &b_softMuonBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonByIP3dBJetTagsAK5PFPUcorrJet", softMuonByIP3dBJetTagsAK5PFPUcorrJet, &b_softMuonByIP3dBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonByPtBJetTagsAK5PFPUcorrJet", softMuonByPtBJetTagsAK5PFPUcorrJet, &b_softMuonByPtBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronBJetTagsAK5PFPUcorrJet", softElectronBJetTagsAK5PFPUcorrJet, &b_softElectronBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronByIP3dBJetTagsAK5PFPUcorrJet", softElectronByIP3dBJetTagsAK5PFPUcorrJet, &b_softElectronByIP3dBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronByPtBJetTagsAK5PFPUcorrJet", softElectronByPtBJetTagsAK5PFPUcorrJet, &b_softElectronByPtBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5PFPUcorrJet", trackCountingHighPurBJetTagsAK5PFPUcorrJet, &b_trackCountingHighPurBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5PFPUcorrJet", trackCountingHighEffBJetTagsAK5PFPUcorrJet, &b_trackCountingHighEffBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("uncorrEnergyAK5PFPUcorrJet", uncorrEnergyAK5PFPUcorrJet, &b_uncorrEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("ptDAK5PFPUcorrJet", ptDAK5PFPUcorrJet, &b_ptDAK5PFPUcorrJet);
   fChain->SetBranchAddress("rmsCandAK5PFPUcorrJet", rmsCandAK5PFPUcorrJet, &b_rmsCandAK5PFPUcorrJet);
   if( isMC_ ) {
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
   }

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
   // cache current run for lumi measurement:
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


bool Ntp1Analyzer::isGoodEvent( int iEntry ) {

     bool okForJSON = false;
     bool okForHLT = false;

     // FIRST STEP: check if LS is a good one in the JSON file:

     if(!filterGoodRuns_) okForJSON = true; // if filtered not requested all events are good

     this->UpdateCache();


     // check whether this run is part of the good runs. else retrun false
     if( goodLSCache_ != goodLS_.end() ) {

        // get list of LS intervals
        const GoodLSVector& lsvector =   goodLSCache_->second; 
        // loop over good LS intervals and return as soon as one interval contains this event
        for(GoodLSVector::const_iterator iLS = lsvector.begin(); (iLS != lsvector.end())&&(okForJSON==false); iLS++) {
     
           if(lumiBlock >= iLS->first && lumiBlock <= iLS->second ) {
    
             okForJSON = true;
          } // check current LS being in the interval
        } // loop over good LS for this run
     }

     // SECOND STEP: if it's ok in the JSON, check if it has passed the trigger
     // once per run, reload trigger mask:
     if( cachedRun_!=runNumber ) {
       if( cachedRun_==0 )
         std::cout << "-> Loading Trigger Mask for run " << runNumber << "." << std::endl;
       else 
         std::cout << "-> Passing from run " << cachedRun_ << " to run " << runNumber << ". Reloading Trigger Mask." << std::endl;
       cachedRun_ = runNumber;
       this->LoadTrigger( iEntry );
     }

     okForHLT = this->PassedHLT(iEntry);

  // if( okForJSON && okForHLT ) { //will take lumi, so (once per LS) increment luminosity
  
  //   if( cachedLS_ != lumiBlock )  {
  //      cachedLS_ = lumiBlock;
  //      RunLumiPair rlpair = (std::pair<int, int>(runNumber,lumiBlock));
  //      totalIntLumi_ += LSLumimap_[rlpair];
  //   }
  // }

     bool returnBool = ( okForJSON && okForHLT );
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





double Ntp1Analyzer::trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}
