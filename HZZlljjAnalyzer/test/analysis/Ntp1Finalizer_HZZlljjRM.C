#include "Ntp1Finalizer_HZZlljjRM.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "CommonTools/fitTools.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"
#include "KinematicFit/DiJetKinFitter.h"

#include "CommonTools/PUWeight.h"




bool ANALYZE_SIDEBANDS_=true;
bool USE_MC_MASS=false;

int DEBUG_EVENTNUMBER = 98901397;





int getNJets( int nPairs );

float getMuonHLTSF_DoubleTrigger( float pt, float eta, const std::string& runPeriod );
float getMuonHLTSF_SingleTrigger( float pt, float eta, const std::string& runPeriod );
float getEventHLTSF( float effSingle1, float effSingle2, float effDouble1, float effDouble2 );

std::vector<TH1D*> getHistoVector(int nPtBins, Double_t *ptBins, std::string histoName, int nBins, float xMin, float xMax );


// constructor:

Ntp1Finalizer_HZZlljjRM::Ntp1Finalizer_HZZlljjRM( const std::string& dataset, const std::string& selectionType, const std::string& PUType, const std::string& leptType ) : Ntp1Finalizer( "HZZlljjRM", dataset, leptType ) {

  if( leptType!="ALL" && leptType!="MU" && leptType!="ELE" ) {
    std::cout << "Lept type '" << leptType << "' currently not supported. Exiting." << std::endl;
    exit(9176);
  }

  leptType_ = leptType;
  PUType_ = PUType;

  setSelectionType(selectionType);

}




void Ntp1Finalizer_HZZlljjRM::finalize() {

  //if( outFile_==0 ) this->createOutputFile();
  
  Int_t run;
  tree_->SetBranchAddress("run", &run);
  tree_->GetEntry(0);
  bool isMC = (run < 10);
  std::string fullFlags = selectionType_;
  if( isMC ) fullFlags = fullFlags + "_PU" + PUType_;
  fullFlags = fullFlags + "_" + leptType_;
  this->set_flags(fullFlags); //this is for the outfile name
  this->createOutputFile();


  TString dataset_tstr(dataset_);


  TTree* tree_passedEvents = new TTree("tree_passedEvents", "Unbinned data for statistical treatment");

  TH1D* h1_nCounter = new TH1D("nCounter", "", 1, 0., 1.);
  h1_nCounter->Sumw2();
  TH1D* h1_nCounterW = new TH1D("nCounterW", "", 1, 0., 1.);
  h1_nCounterW->Sumw2();
  TH1D* h1_nCounterPU = new TH1D("nCounterPU", "", 1, 0., 1.);
  h1_nCounterPU->Sumw2();

  TH1D* h1_nEventsCategories_presel = new TH1D("nEventsCategories_presel", "", 4, -1.5, 2.5);
  h1_nEventsCategories_presel->Sumw2();
  h1_nEventsCategories_presel->GetXaxis()->SetLabelSize(0.1);
  h1_nEventsCategories_presel->GetXaxis()->SetBinLabel(1, "Gluon-tag"); 
  h1_nEventsCategories_presel->GetXaxis()->SetBinLabel(2, "0 b-tag"); 
  h1_nEventsCategories_presel->GetXaxis()->SetBinLabel(3, "1 b-tag"); 
  h1_nEventsCategories_presel->GetXaxis()->SetBinLabel(4, "2 b-tag"); 

  TH1D* h1_nEvents_partFlavor_beforeQG = new TH1D("nEvents_partFlavor_beforeQG", "", 2, -0.5, 1.5);
  h1_nEvents_partFlavor_beforeQG->Sumw2();
  h1_nEvents_partFlavor_beforeQG->GetXaxis()->SetLabelSize(0.1);
  h1_nEvents_partFlavor_beforeQG->GetXaxis()->SetBinLabel(1, "2 Quarks");
  h1_nEvents_partFlavor_beforeQG->GetXaxis()->SetBinLabel(1, "#geq 1 Gluon");

  TH1D* h1_nEvents_partFlavor_afterQG = new TH1D("nEvents_partFlavor_afterQG", "", 2, -0.5, 1.5);
  h1_nEvents_partFlavor_afterQG->Sumw2();
  h1_nEvents_partFlavor_afterQG->GetXaxis()->SetLabelSize(0.1);
  h1_nEvents_partFlavor_afterQG->GetXaxis()->SetBinLabel(1, "2 Quarks");
  h1_nEvents_partFlavor_afterQG->GetXaxis()->SetBinLabel(1, "#geq 1 Gluon");


  TH1D* h1_nEventsMW400_partFlavor_beforeQG = new TH1D("nEventsMW400_partFlavor_beforeQG", "", 2, -0.5, 1.5);
  h1_nEventsMW400_partFlavor_beforeQG->Sumw2();
  h1_nEventsMW400_partFlavor_beforeQG->GetXaxis()->SetLabelSize(0.1);
  h1_nEventsMW400_partFlavor_beforeQG->GetXaxis()->SetBinLabel(1, "2 Quarks");
  h1_nEventsMW400_partFlavor_beforeQG->GetXaxis()->SetBinLabel(1, "#geq 1 Gluon");

  TH1D* h1_nEventsMW400_partFlavor_afterQG = new TH1D("nEventsMW400_partFlavor_afterQG", "", 2, -0.5, 1.5);
  h1_nEventsMW400_partFlavor_afterQG->Sumw2();
  h1_nEventsMW400_partFlavor_afterQG->GetXaxis()->SetLabelSize(0.1);
  h1_nEventsMW400_partFlavor_afterQG->GetXaxis()->SetBinLabel(1, "2 Quarks");
  h1_nEventsMW400_partFlavor_afterQG->GetXaxis()->SetBinLabel(1, "#geq 1 Gluon");


  TH1D* h1_nEventsMW500_partFlavor_beforeQG = new TH1D("nEventsMW500_partFlavor_beforeQG", "", 2, -0.5, 1.5);
  h1_nEventsMW500_partFlavor_beforeQG->Sumw2();
  h1_nEventsMW500_partFlavor_beforeQG->GetXaxis()->SetLabelSize(0.1);
  h1_nEventsMW500_partFlavor_beforeQG->GetXaxis()->SetBinLabel(1, "2 Quarks");
  h1_nEventsMW500_partFlavor_beforeQG->GetXaxis()->SetBinLabel(1, "#geq 1 Gluon");

  TH1D* h1_nEventsMW500_partFlavor_afterQG = new TH1D("nEventsMW500_partFlavor_afterQG", "", 2, -0.5, 1.5);
  h1_nEventsMW500_partFlavor_afterQG->Sumw2();
  h1_nEventsMW500_partFlavor_afterQG->GetXaxis()->SetLabelSize(0.1);
  h1_nEventsMW500_partFlavor_afterQG->GetXaxis()->SetBinLabel(1, "2 Quarks");
  h1_nEventsMW500_partFlavor_afterQG->GetXaxis()->SetBinLabel(1, "#geq 1 Gluon");


  // these histograms will save the final yields and efficiencies:

  TH1D* h1_nEvents_fb_gluetag_250 = new TH1D("nEvents_fb_gluetag_250", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_250 = new TH1D("nEvents_fb_0btag_250", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_250 = new TH1D("nEvents_fb_1btag_250", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_250 = new TH1D("nEvents_fb_2btag_250", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_250_ELE = new TH1D("nEvents_fb_gluetag_250_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_250_ELE = new TH1D("nEvents_fb_0btag_250_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_250_ELE = new TH1D("nEvents_fb_1btag_250_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_250_ELE = new TH1D("nEvents_fb_2btag_250_ELE", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_250_MU = new TH1D("nEvents_fb_gluetag_250_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_250_MU = new TH1D("nEvents_fb_0btag_250_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_250_MU = new TH1D("nEvents_fb_1btag_250_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_250_MU = new TH1D("nEvents_fb_2btag_250_MU", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_250 = new TH1D("eff_gluetag_250", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_250 = new TH1D("eff_0btag_250", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_250 = new TH1D("eff_1btag_250", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_250 = new TH1D("eff_2btag_250", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_250_ELE = new TH1D("eff_gluetag_250_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_250_ELE = new TH1D("eff_0btag_250_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_250_ELE = new TH1D("eff_1btag_250_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_250_ELE = new TH1D("eff_2btag_250_ELE", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_250_MU = new TH1D("eff_gluetag_250_MU", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_250_MU = new TH1D("eff_0btag_250_MU", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_250_MU = new TH1D("eff_1btag_250_MU", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_250_MU = new TH1D("eff_2btag_250_MU", "", 1, 0., 1.);


  TH1D* h1_nEvents_fb_gluetag_300 = new TH1D("nEvents_fb_gluetag_300", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_300 = new TH1D("nEvents_fb_0btag_300", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_300 = new TH1D("nEvents_fb_1btag_300", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_300 = new TH1D("nEvents_fb_2btag_300", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_300_ELE = new TH1D("nEvents_fb_gluetag_300_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_300_ELE = new TH1D("nEvents_fb_0btag_300_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_300_ELE = new TH1D("nEvents_fb_1btag_300_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_300_ELE = new TH1D("nEvents_fb_2btag_300_ELE", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_300_MU = new TH1D("nEvents_fb_gluetag_300_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_300_MU = new TH1D("nEvents_fb_0btag_300_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_300_MU = new TH1D("nEvents_fb_1btag_300_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_300_MU = new TH1D("nEvents_fb_2btag_300_MU", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_300 = new TH1D("eff_gluetag_300", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_300 = new TH1D("eff_0btag_300", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_300 = new TH1D("eff_1btag_300", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_300 = new TH1D("eff_2btag_300", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_300_ELE = new TH1D("eff_gluetag_300_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_300_ELE = new TH1D("eff_0btag_300_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_300_ELE = new TH1D("eff_1btag_300_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_300_ELE = new TH1D("eff_2btag_300_ELE", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_300_MU = new TH1D("eff_gluetag_300_MU", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_300_MU = new TH1D("eff_0btag_300_MU", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_300_MU = new TH1D("eff_1btag_300_MU", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_300_MU = new TH1D("eff_2btag_300_MU", "", 1, 0., 1.);


  TH1D* h1_nEvents_fb_gluetag_350 = new TH1D("nEvents_fb_gluetag_350", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_0btag_350 = new TH1D("nEvents_fb_0btag_350", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_1btag_350 = new TH1D("nEvents_fb_1btag_350", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_2btag_350 = new TH1D("nEvents_fb_2btag_350", "", 1, 0., 1.0001);

  TH1D* h1_nEvents_fb_gluetag_350_ELE = new TH1D("nEvents_fb_gluetag_350_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_350_ELE = new TH1D("nEvents_fb_0btag_350_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_350_ELE = new TH1D("nEvents_fb_1btag_350_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_350_ELE = new TH1D("nEvents_fb_2btag_350_ELE", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_350_MU = new TH1D("nEvents_fb_gluetag_350_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_350_MU = new TH1D("nEvents_fb_0btag_350_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_350_MU = new TH1D("nEvents_fb_1btag_350_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_350_MU = new TH1D("nEvents_fb_2btag_350_MU", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_350 = new TH1D("eff_gluetag_350", "", 1, 0., 1.0001);
  TH1D* h1_eff_0btag_350 = new TH1D("eff_0btag_350", "", 1, 0., 1.0001);
  TH1D* h1_eff_1btag_350 = new TH1D("eff_1btag_350", "", 1, 0., 1.0001);
  TH1D* h1_eff_2btag_350 = new TH1D("eff_2btag_350", "", 1, 0., 1.0001);

  TH1D* h1_eff_gluetag_350_ELE = new TH1D("eff_gluetag_350_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_350_ELE = new TH1D("eff_0btag_350_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_350_ELE = new TH1D("eff_1btag_350_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_350_ELE = new TH1D("eff_2btag_350_ELE", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_350_MU = new TH1D("eff_gluetag_350_MU", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_350_MU = new TH1D("eff_0btag_350_MU", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_350_MU = new TH1D("eff_1btag_350_MU", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_350_MU = new TH1D("eff_2btag_350_MU", "", 1, 0., 1.);


  TH1D* h1_nEvents_fb_gluetag_400 = new TH1D("nEvents_fb_gluetag_400", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_0btag_400 = new TH1D("nEvents_fb_0btag_400", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_1btag_400 = new TH1D("nEvents_fb_1btag_400", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_2btag_400 = new TH1D("nEvents_fb_2btag_400", "", 1, 0., 1.0001);

  TH1D* h1_nEvents_fb_gluetag_400_ELE = new TH1D("nEvents_fb_gluetag_400_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_400_ELE = new TH1D("nEvents_fb_0btag_400_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_400_ELE = new TH1D("nEvents_fb_1btag_400_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_400_ELE = new TH1D("nEvents_fb_2btag_400_ELE", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_400_MU = new TH1D("nEvents_fb_gluetag_400_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_400_MU = new TH1D("nEvents_fb_0btag_400_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_400_MU = new TH1D("nEvents_fb_1btag_400_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_400_MU = new TH1D("nEvents_fb_2btag_400_MU", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_400 = new TH1D("eff_gluetag_400", "", 1, 0., 1.0001);
  TH1D* h1_eff_0btag_400 = new TH1D("eff_0btag_400", "", 1, 0., 1.0001);
  TH1D* h1_eff_1btag_400 = new TH1D("eff_1btag_400", "", 1, 0., 1.0001);
  TH1D* h1_eff_2btag_400 = new TH1D("eff_2btag_400", "", 1, 0., 1.0001);

  TH1D* h1_eff_gluetag_400_ELE = new TH1D("eff_gluetag_400_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_400_ELE = new TH1D("eff_0btag_400_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_400_ELE = new TH1D("eff_1btag_400_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_400_ELE = new TH1D("eff_2btag_400_ELE", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_400_MU = new TH1D("eff_gluetag_400_MU", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_400_MU = new TH1D("eff_0btag_400_MU", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_400_MU = new TH1D("eff_1btag_400_MU", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_400_MU = new TH1D("eff_2btag_400_MU", "", 1, 0., 1.);


  TH1D* h1_nEvents_fb_gluetag_450 = new TH1D("nEvents_fb_gluetag_450", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_0btag_450 = new TH1D("nEvents_fb_0btag_450", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_1btag_450 = new TH1D("nEvents_fb_1btag_450", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_2btag_450 = new TH1D("nEvents_fb_2btag_450", "", 1, 0., 1.0001);

  TH1D* h1_nEvents_fb_gluetag_450_ELE = new TH1D("nEvents_fb_gluetag_450_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_450_ELE = new TH1D("nEvents_fb_0btag_450_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_450_ELE = new TH1D("nEvents_fb_1btag_450_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_450_ELE = new TH1D("nEvents_fb_2btag_450_ELE", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_450_MU = new TH1D("nEvents_fb_gluetag_450_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_450_MU = new TH1D("nEvents_fb_0btag_450_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_450_MU = new TH1D("nEvents_fb_1btag_450_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_450_MU = new TH1D("nEvents_fb_2btag_450_MU", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_450 = new TH1D("eff_gluetag_450", "", 1, 0., 1.0001);
  TH1D* h1_eff_0btag_450 = new TH1D("eff_0btag_450", "", 1, 0., 1.0001);
  TH1D* h1_eff_1btag_450 = new TH1D("eff_1btag_450", "", 1, 0., 1.0001);
  TH1D* h1_eff_2btag_450 = new TH1D("eff_2btag_450", "", 1, 0., 1.0001);

  TH1D* h1_eff_gluetag_450_ELE = new TH1D("eff_gluetag_450_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_450_ELE = new TH1D("eff_0btag_450_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_450_ELE = new TH1D("eff_1btag_450_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_450_ELE = new TH1D("eff_2btag_450_ELE", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_450_MU = new TH1D("eff_gluetag_450_MU", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_450_MU = new TH1D("eff_0btag_450_MU", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_450_MU = new TH1D("eff_1btag_450_MU", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_450_MU = new TH1D("eff_2btag_450_MU", "", 1, 0., 1.);


  TH1D* h1_nEvents_fb_gluetag_500 = new TH1D("nEvents_fb_gluetag_500", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_0btag_500 = new TH1D("nEvents_fb_0btag_500", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_1btag_500 = new TH1D("nEvents_fb_1btag_500", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_2btag_500 = new TH1D("nEvents_fb_2btag_500", "", 1, 0., 1.0001);

  TH1D* h1_nEvents_fb_gluetag_500_ELE = new TH1D("nEvents_fb_gluetag_500_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_500_ELE = new TH1D("nEvents_fb_0btag_500_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_500_ELE = new TH1D("nEvents_fb_1btag_500_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_500_ELE = new TH1D("nEvents_fb_2btag_500_ELE", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_500_MU = new TH1D("nEvents_fb_gluetag_500_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_500_MU = new TH1D("nEvents_fb_0btag_500_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_500_MU = new TH1D("nEvents_fb_1btag_500_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_500_MU = new TH1D("nEvents_fb_2btag_500_MU", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_500 = new TH1D("eff_gluetag_500", "", 1, 0., 1.0001);
  TH1D* h1_eff_0btag_500 = new TH1D("eff_0btag_500", "", 1, 0., 1.0001);
  TH1D* h1_eff_1btag_500 = new TH1D("eff_1btag_500", "", 1, 0., 1.0001);
  TH1D* h1_eff_2btag_500 = new TH1D("eff_2btag_500", "", 1, 0., 1.0001);

  TH1D* h1_eff_gluetag_500_ELE = new TH1D("eff_gluetag_500_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_500_ELE = new TH1D("eff_0btag_500_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_500_ELE = new TH1D("eff_1btag_500_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_500_ELE = new TH1D("eff_2btag_500_ELE", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_500_MU = new TH1D("eff_gluetag_500_MU", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_500_MU = new TH1D("eff_0btag_500_MU", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_500_MU = new TH1D("eff_1btag_500_MU", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_500_MU = new TH1D("eff_2btag_500_MU", "", 1, 0., 1.);


  TH1D* h1_nEvents_fb_gluetag_600 = new TH1D("nEvents_fb_gluetag_600", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_0btag_600 = new TH1D("nEvents_fb_0btag_600", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_1btag_600 = new TH1D("nEvents_fb_1btag_600", "", 1, 0., 1.0001);
  TH1D* h1_nEvents_fb_2btag_600 = new TH1D("nEvents_fb_2btag_600", "", 1, 0., 1.0001);

  TH1D* h1_nEvents_fb_gluetag_600_ELE = new TH1D("nEvents_fb_gluetag_600_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_600_ELE = new TH1D("nEvents_fb_0btag_600_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_600_ELE = new TH1D("nEvents_fb_1btag_600_ELE", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_600_ELE = new TH1D("nEvents_fb_2btag_600_ELE", "", 1, 0., 1.);

  TH1D* h1_nEvents_fb_gluetag_600_MU = new TH1D("nEvents_fb_gluetag_600_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_0btag_600_MU = new TH1D("nEvents_fb_0btag_600_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_1btag_600_MU = new TH1D("nEvents_fb_1btag_600_MU", "", 1, 0., 1.);
  TH1D* h1_nEvents_fb_2btag_600_MU = new TH1D("nEvents_fb_2btag_600_MU", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_600 = new TH1D("eff_gluetag_600", "", 1, 0., 1.0001);
  TH1D* h1_eff_0btag_600 = new TH1D("eff_0btag_600", "", 1, 0., 1.0001);
  TH1D* h1_eff_1btag_600 = new TH1D("eff_1btag_600", "", 1, 0., 1.0001);
  TH1D* h1_eff_2btag_600 = new TH1D("eff_2btag_600", "", 1, 0., 1.0001);

  TH1D* h1_eff_gluetag_600_ELE = new TH1D("eff_gluetag_600_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_600_ELE = new TH1D("eff_0btag_600_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_600_ELE = new TH1D("eff_1btag_600_ELE", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_600_ELE = new TH1D("eff_2btag_600_ELE", "", 1, 0., 1.);

  TH1D* h1_eff_gluetag_600_MU = new TH1D("eff_gluetag_600_MU", "", 1, 0., 1.);
  TH1D* h1_eff_0btag_600_MU = new TH1D("eff_0btag_600_MU", "", 1, 0., 1.);
  TH1D* h1_eff_1btag_600_MU = new TH1D("eff_1btag_600_MU", "", 1, 0., 1.);
  TH1D* h1_eff_2btag_600_MU = new TH1D("eff_2btag_600_MU", "", 1, 0., 1.);



  TH1F* h1_run = new TH1F("run", "", 15149, 132440, 147589);

  TH1D* h1_nvertex = new TH1D("nvertex", "", 36, -0.5, 35.5);
  h1_nvertex->Sumw2();
  TH1D* h1_nvertex_PUW = new TH1D("nvertex_PUW", "", 36, -0.5, 35.5);
  h1_nvertex_PUW->Sumw2();
  TH1D* h1_nvertex_PUW_ave = new TH1D("nvertex_PUW_ave", "", 36, -0.5, 35.5);
  h1_nvertex_PUW_ave->Sumw2();

  TH1D* h1_pfMet = new TH1D("pfMet", "", 60, 0., 120.);
  h1_pfMet->Sumw2();
  TH1D* h1_pfMet_2btag = new TH1D("pfMet_2btag", "", 60, 0., 120.);
  h1_pfMet_2btag->Sumw2();

  TH1D* h1_pfMetOverMZZ= new TH1D("pfMetOverMZZ", "", 100, 0., 0.6);
  h1_pfMetOverMZZ->Sumw2();
  TH1D* h1_pfMetOverMZZ_2btag = new TH1D("pfMetOverMZZ_2btag", "", 100, 0., 0.6);
  h1_pfMetOverMZZ_2btag->Sumw2();

  TH1D* h1_metSignificance= new TH1D("metSignificance", "", 80, 0., 40.);
  h1_metSignificance->Sumw2();
  TH1D* h1_metSignificance_2btag = new TH1D("metSignificance_2btag", "", 80, 0., 40.);
  h1_metSignificance_2btag->Sumw2();

  TH1D* h1_mEtSig= new TH1D("mEtSig", "", 60, 0., 15.);
  h1_mEtSig->Sumw2();
  TH1D* h1_mEtSig_2btag = new TH1D("mEtSig_2btag", "", 60, 0., 15.);
  h1_mEtSig_2btag->Sumw2();


  TH1D* h1_rhoPF_presel = new TH1D("rhoPF_presel", "", 50, 0., 20.);
  h1_rhoPF_presel->Sumw2();
  TH1D* h1_rhoPF = new TH1D("rhoPF", "", 50, 0., 20.);
  h1_rhoPF->Sumw2();

  TH1D* h1_nCandidates = new TH1D("nCandidates", "", 6, -0.5, 5.5 );
  h1_nCandidates->Sumw2();

  TH1D* h1_ptLept1_presel = new TH1D("ptLept1_presel", "", 400, 40., 440.);
  h1_ptLept1_presel->Sumw2();
  TH1D* h1_ptLept2_presel = new TH1D("ptLept2_presel", "", 200, 20., 220.);
  h1_ptLept2_presel->Sumw2();
  TH1D* h1_etaLept1_presel = new TH1D("etaLept1_presel", "", 25, -2.5, 2.5);
  h1_etaLept1_presel->Sumw2();
  TH1D* h1_etaLept2_presel = new TH1D("etaLept2_presel", "", 25, -2.5, 2.5);
  h1_etaLept2_presel->Sumw2();

  TH1D* h1_ptJet_all_presel = new TH1D("ptJet_all_presel", "", 200, 30., 430.);
  h1_ptJet_all_presel->Sumw2();
  TH1D* h1_etaJet_all_presel = new TH1D("etaJet_all_presel", "", 25, -5., 5.);
  h1_etaJet_all_presel->Sumw2();

  const int nPtBins = 20;
  Double_t ptBins[nPtBins+1];
  fitTools::getBins_int( nPtBins+1, ptBins, 15., 1000. );

//std::vector<TH1D*> vh1_ptDJet_all_presel = getHistoVector(nPtBins, ptBins, "ptDJet_all_presel", 60, 0., 1.);
//std::vector<TH1D*> vh1_rmsCandJet_all_presel = getHistoVector(nPtBins, ptBins, "rmsCandJet_all_presel", 60, 0., 0.07);
//std::vector<TH1D*> vh1_nChargedJet_all_presel = getHistoVector(nPtBins, ptBins, "nChargedJet_all_presel", 41, -0.5, 40.5);
//std::vector<TH1D*> vh1_nNeutralJet_all_presel = getHistoVector(nPtBins, ptBins, "nNeutralJet_all_presel", 41, -0.5, 40.5);

  TH1D* h1_nJets_presel = new TH1D("nJets_presel", "", 7, 1.5, 8.5);
  h1_nJets_presel->Sumw2();
  TH1D* h1_nPairs_presel = new TH1D("nPairs_presel", "", 21, 0.5, 21.5);
  h1_nPairs_presel->Sumw2();

  TH1D* h1_ptLept1= new TH1D("ptLept1", "", 400, 40., 440.);
  h1_ptLept1->Sumw2();
  TH1D* h1_ptLept2= new TH1D("ptLept2", "", 200, 20., 220.);
  h1_ptLept2->Sumw2();
  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 400, 30., 430.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 200, 30., 230.);
  h1_ptJet2->Sumw2();
  TH1D* h1_ptJet1_prekin = new TH1D("ptJet1_prekin", "", 200, 30., 430.);
  h1_ptJet1_prekin->Sumw2();
  TH1D* h1_ptJet2_prekin = new TH1D("ptJet2_prekin", "", 100, 30., 230.);
  h1_ptJet2_prekin->Sumw2();

  TH1D* h1_ptJetRecoil = new TH1D("ptJetRecoil", "", 150, 0., 150.);
  h1_ptJetRecoil->Sumw2();
  TH1D* h1_ptHiggs = new TH1D("ptHiggs", "", 150, 0., 150.);
  h1_ptHiggs->Sumw2();

  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 100, -2.4, 2.4);
  h1_etaJet1->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 100, -2.4, 2.4);
  h1_etaJet2->Sumw2();

  TH1D* h1_tcheJet = new TH1D("tcheJet", "", 100, 0., 20.);
  h1_tcheJet->Sumw2();
  TH1D* h1_tcheJet1 = new TH1D("tcheJet1", "", 100, 0., 20.);
  h1_tcheJet1->Sumw2();
  TH1D* h1_tcheJet2 = new TH1D("tcheJet2", "", 100, 0., 20.);
  h1_tcheJet2->Sumw2();

  TH1D* h1_eElectronsJet1 = new TH1D("eElectronsJet1", "", 60, 0., 1.0001);
  h1_eElectronsJet1->Sumw2();
  TH1D* h1_eElectronsJet2 = new TH1D("eElectronsJet2", "", 60, 0., 1.0001);
  h1_eElectronsJet2->Sumw2();
  TH1D* h1_eMuonsJet1 = new TH1D("eMuonsJet1", "", 60, 0., 1.0001);
  h1_eMuonsJet1->Sumw2();
  TH1D* h1_eMuonsJet2 = new TH1D("eMuonsJet2", "", 60, 0., 1.0001);
  h1_eMuonsJet2->Sumw2();

  TH1D* h1_ptResoJet1_beforeKin = new TH1D("ptResoJet1_beforeKin", "", 100, -1., 1.);
  h1_ptResoJet1_beforeKin->Sumw2();
  TH1D* h1_ptResoJet2_beforeKin = new TH1D("ptResoJet2_beforeKin", "", 100, -1., 1.);
  h1_ptResoJet2_beforeKin->Sumw2();

  TH1D* h1_ptResoJet1_afterKin = new TH1D("ptResoJet1_afterKin", "", 100, -1., 1.);
  h1_ptResoJet1_afterKin->Sumw2();
  TH1D* h1_ptResoJet2_afterKin = new TH1D("ptResoJet2_afterKin", "", 100, -1., 1.);
  h1_ptResoJet2_afterKin->Sumw2();

  TH1D* h1_ptZreso_beforeKin = new TH1D("ptZreso_beforeKin", "", 100, -1., 1.);
  h1_ptZreso_beforeKin->Sumw2();
  TH1D* h1_ptZreso_afterKin = new TH1D("ptZreso_afterKin", "", 100, -1., 1.);
  h1_ptZreso_afterKin->Sumw2();

  TH1D* h1_mHreso_beforeKin = new TH1D("mHreso_beforeKin", "", 100, -1., 1.);
  h1_mHreso_beforeKin->Sumw2();
  TH1D* h1_mHreso_afterKin = new TH1D("mHreso_afterKin", "", 100, -1., 1.);
  h1_mHreso_afterKin->Sumw2();


  TH1D* h1_mZjjMC = new TH1D("mZjjMC", "", 400, 30., 430.);
  h1_mZjjMC->Sumw2();

  TH1D* h1_mZjj_all_presel = new TH1D("mZjj_all_presel", "", 400, 30., 430.);
  h1_mZjj_all_presel->Sumw2();

  TH1D* h1_deltaRjj_all_presel = new TH1D("deltaRjj_all_presel", "", 500, 0.5, 5.);
  h1_deltaRjj_all_presel->Sumw2();
  TH1D* h1_deltaRll_presel = new TH1D("deltaRll_presel", "", 100, 0., 5.);
  h1_deltaRll_presel->Sumw2();

  TH1D* h1_mZll = new TH1D("mZll", "", 120, 60., 120.);
  h1_mZll->Sumw2();
  TH1D* h1_mZll_presel = new TH1D("mZll_presel", "", 120, 60., 120.);
  h1_mZll_presel->Sumw2();
  TH1D* h1_mZll_presel_0jets = new TH1D("mZll_presel_0jets", "", 120, 60., 120.);
  h1_mZll_presel_0jets->Sumw2();

  TH1D* h1_mZmumu = new TH1D("mZmumu", "", 120, 60., 120.);
  h1_mZmumu->Sumw2();
  TH1D* h1_mZmumu_presel = new TH1D("mZmumu_presel", "", 120, 60., 120.);
  h1_mZmumu_presel->Sumw2();
  TH1D* h1_mZmumu_presel_0jets = new TH1D("mZmumu_presel_0jets", "", 120, 60., 120.);
  h1_mZmumu_presel_0jets->Sumw2();
  TH1D* h1_mZee = new TH1D("mZee", "", 120, 60., 120.);
  h1_mZee->Sumw2();
  TH1D* h1_mZee_presel = new TH1D("mZee_presel", "", 120, 60., 120.);
  h1_mZee_presel->Sumw2();
  TH1D* h1_mZee_presel_0jets = new TH1D("mZee_presel_0jets", "", 120, 60., 120.);
  h1_mZee_presel_0jets->Sumw2();


  TH1D* h1_deltaR_part1 = new TH1D("deltaR_part1", "", 50, 0., 0.8);
  h1_deltaR_part1->Sumw2();
  //TH1D* h1_partFlavorJet1= new TH1D("partFlavorJet1", "", 30, -7.5, 22.5);
  TH1D* h1_partFlavorJet1= new TH1D("partFlavorJet1", "", 38, -15.5, 22.5);
  h1_partFlavorJet1->Sumw2();
  TH1D* h1_partFlavorJet1_MW400= new TH1D("partFlavorJet1_MW400", "", 38, -15.5, 22.5);
  h1_partFlavorJet1_MW400->Sumw2();
  TH1D* h1_partFlavorJet1_MW500= new TH1D("partFlavorJet1_MW500", "", 38, -15.5, 22.5);
  h1_partFlavorJet1_MW500->Sumw2();
  TH1D* h1_partFlavorJet1_matched= new TH1D("partFlavorJet1_matched", "", 38, -15.5, 22.5);
  h1_partFlavorJet1_matched->Sumw2();
  TH1D* h1_partFlavorJet1_notmatched= new TH1D("partFlavorJet1_notmatched", "", 38, -15.5, 22.5);
  h1_partFlavorJet1_notmatched->Sumw2();
  TH1D* h1_partFlavorJet1_0btag= new TH1D("partFlavorJet1_0btag", "", 38, -15.5, 22.5);
  h1_partFlavorJet1_0btag->Sumw2();
  TH1D* h1_partFlavorJet1_1btag= new TH1D("partFlavorJet1_1btag", "", 38, -15.5, 22.5);
  h1_partFlavorJet1_1btag->Sumw2();
  TH1D* h1_partFlavorJet1_2btag= new TH1D("partFlavorJet1_2btag", "", 38, -15.5, 22.5);
  h1_partFlavorJet1_2btag->Sumw2();
  TH1D* h1_partFlavorJet1_gluetag= new TH1D("partFlavorJet1_gluetag", "", 38, -15.5, 22.5);
  h1_partFlavorJet1_gluetag->Sumw2();
//std::vector<TH1D*> vh1_ptDJet1 = getHistoVector(nPtBins, ptBins, "ptDJet1", 60, 0., 1.);
//std::vector<TH1D*> vh1_rmsCandJet1 = getHistoVector(nPtBins, ptBins, "rmsCandJet1", 60, 0., 0.1);
//std::vector<TH1D*> vh1_nChargedJet1 = getHistoVector(nPtBins, ptBins, "nChargedJet1", 51, -0.5, 50.5);
//std::vector<TH1D*> vh1_nNeutralJet1 = getHistoVector(nPtBins, ptBins, "nNeutralJet1", 51, -0.5, 50.5);
//std::vector<TH1D*> vh1_QGLikelihoodJet1 = getHistoVector(nPtBins, ptBins, "QGLikelihoodJet1", 60, 0., 1.);
  TH1D* h1_nChargedJet1 = new TH1D("nChargedJet1", "", 31, -0.5, 30.5);
  h1_nChargedJet1->Sumw2();
  TH1D* h1_nNeutralJet1 = new TH1D("nNeutralJet1", "", 31, -0.5, 30.5);
  h1_nNeutralJet1->Sumw2();
  TH1D* h1_ptDJet1 = new TH1D("ptDJet1", "", 50, 0., 1.0001);
  h1_ptDJet1->Sumw2();
  TH1D* h1_QGLikelihoodNoPUJet1 = new TH1D("QGLikelihoodNoPUJet1", "", 60, 0., 1.0001);
  h1_QGLikelihoodNoPUJet1->Sumw2();
  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet1->Sumw2();
  TH1D* h1_QGLikelihoodJet1_MW300 = new TH1D("QGLikelihoodJet1_MW300", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet1_MW300->Sumw2();
  TH1D* h1_QGLikelihoodJet1_MW400 = new TH1D("QGLikelihoodJet1_MW400", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet1_MW400->Sumw2();
  TH1D* h1_QGLikelihoodJet1_MW500 = new TH1D("QGLikelihoodJet1_MW500", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet1_MW500->Sumw2();
  TH1D* h1_QGLikelihoodJet1_gluonMatched = new TH1D("QGLikelihoodJet1_gluonMatched", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet1_gluonMatched->Sumw2();
  TH1D* h1_QGLikelihoodJet1_quarkMatched = new TH1D("QGLikelihoodJet1_quarkMatched", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet1_quarkMatched->Sumw2();
  TH1D* h1_QGLikelihoodJetRecoil = new TH1D("QGLikelihoodJetRecoil", "", 60, 0., 1.0001);
  h1_QGLikelihoodJetRecoil->Sumw2();
  TH1D* h1_QGLikelihoodJetRecoil_MW400 = new TH1D("QGLikelihoodJetRecoil_MW400", "", 60, 0., 1.0001);
  h1_QGLikelihoodJetRecoil_MW400->Sumw2();
  TH1D* h1_QGLikelihoodProdRecoil = new TH1D("QGLikelihoodProdRecoil", "", 60, 0., 1.0001);
  h1_QGLikelihoodProdRecoil->Sumw2();
  TH1D* h1_QGLikelihoodProdRecoil_MW400 = new TH1D("QGLikelihoodProdRecoil_MW400", "", 60, 0., 1.0001);
  h1_QGLikelihoodProdRecoil_MW400->Sumw2();

  TH1D* h1_QGLikelihood_100_123 = new TH1D("QGLikelihood_100_123", "", 50, 0., 1.0001);
  h1_QGLikelihood_100_123->Sumw2();
  TH1D* h1_QGLikelihood_66_81 = new TH1D("QGLikelihood_66_81", "", 50, 0., 1.0001);
  h1_QGLikelihood_66_81->Sumw2();


  TH1D* h1_deltaR_part2 = new TH1D("deltaR_part2", "", 60, 0., 0.8);
  h1_deltaR_part2->Sumw2();
  TH1D* h1_partFlavorJet2= new TH1D("partFlavorJet2", "", 38, -15.5, 22.5);
  h1_partFlavorJet2->Sumw2();
  TH1D* h1_partFlavorJet2_MW400= new TH1D("partFlavorJet2_MW400", "", 38, -15.5, 22.5);
  h1_partFlavorJet2_MW400->Sumw2();
  TH1D* h1_partFlavorJet2_MW500= new TH1D("partFlavorJet2_MW500", "", 38, -15.5, 22.5);
  h1_partFlavorJet2_MW500->Sumw2();
  TH1D* h1_partFlavorJet2_matched= new TH1D("partFlavorJet2_matched", "", 38, -15.5, 22.5);
  h1_partFlavorJet2_matched->Sumw2();
  TH1D* h1_partFlavorJet2_notmatched= new TH1D("partFlavorJet2_notmatched", "", 38, -15.5, 22.5);
  h1_partFlavorJet2_notmatched->Sumw2();
  TH1D* h1_partFlavorJet2_0btag= new TH1D("partFlavorJet2_0btag", "", 38, -15.5, 22.5);
  h1_partFlavorJet2_0btag->Sumw2();
  TH1D* h1_partFlavorJet2_1btag= new TH1D("partFlavorJet2_1btag", "", 38, -15.5, 22.5);
  h1_partFlavorJet2_1btag->Sumw2();
  TH1D* h1_partFlavorJet2_2btag= new TH1D("partFlavorJet2_2btag", "", 38, -15.5, 22.5);
  h1_partFlavorJet2_2btag->Sumw2();
  TH1D* h1_partFlavorJet2_gluetag= new TH1D("partFlavorJet2_gluetag", "", 38, -15.5, 22.5);
  h1_partFlavorJet2_gluetag->Sumw2();

  //std::vector<TH1D*> vh1_ptDJet2 = getHistoVector(nPtBins, ptBins, "ptDJet2", 60, 0., 1.);
  //std::vector<TH1D*> vh1_rmsCandJet2 = getHistoVector(nPtBins, ptBins, "rmsCandJet2", 60, 0., 0.1);
  //std::vector<TH1D*> vh1_nChargedJet2 = getHistoVector(nPtBins, ptBins, "nChargedJet2", 51, -0.5, 50.5);
  //std::vector<TH1D*> vh1_nNeutralJet2 = getHistoVector(nPtBins, ptBins, "nNeutralJet2", 51, -0.5, 50.5);
  //std::vector<TH1D*> vh1_QGLikelihoodJet2 = getHistoVector(nPtBins, ptBins, "QGLikelihoodJet2", 60, 0., 1.);
  TH1D* h1_nChargedJet2 = new TH1D("nChargedJet2", "", 31, -0.5, 30.5);
  h1_nChargedJet2->Sumw2();
  TH1D* h1_nNeutralJet2 = new TH1D("nNeutralJet2", "", 31, -0.5, 30.5);
  h1_nNeutralJet2->Sumw2();
  TH1D* h1_ptDJet2 = new TH1D("ptDJet2", "", 50, 0., 1.0001);
  h1_ptDJet2->Sumw2();
  TH1D* h1_QGLikelihoodNoPUJet2 = new TH1D("QGLikelihoodNoPUJet2", "", 60, 0., 1.0001);
  h1_QGLikelihoodNoPUJet2->Sumw2();
  TH1D* h1_QGLikelihoodJet2 = new TH1D("QGLikelihoodJet2", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet2->Sumw2();
  TH1D* h1_QGLikelihoodJet2_MW300 = new TH1D("QGLikelihoodJet2_MW300", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet2_MW300->Sumw2();
  TH1D* h1_QGLikelihoodJet2_MW400 = new TH1D("QGLikelihoodJet2_MW400", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet2_MW400->Sumw2();
  TH1D* h1_QGLikelihoodJet2_MW500 = new TH1D("QGLikelihoodJet2_MW500", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet2_MW500->Sumw2();
  TH1D* h1_QGLikelihoodJet2_gluonMatched = new TH1D("QGLikelihoodJet2_gluonMatched", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet2_gluonMatched->Sumw2();
  TH1D* h1_QGLikelihoodJet2_quarkMatched = new TH1D("QGLikelihoodJet2_quarkMatched", "", 60, 0., 1.0001);
  h1_QGLikelihoodJet2_quarkMatched->Sumw2();

  TH1D* h1_QGLikelihoodProd_oneGluon = new TH1D("QGLikelihoodProd_oneGluon", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_oneGluon->Sumw2();
  TH1D* h1_QGLikelihoodProd_twoGluon = new TH1D("QGLikelihoodProd_twoGluon", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_twoGluon->Sumw2();
  TH1D* h1_QGLikelihoodProd_twoQuark = new TH1D("QGLikelihoodProd_twoQuark", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_twoQuark->Sumw2();
  TH1D* h1_QGLikelihoodProd_twoQuark_MW400 = new TH1D("QGLikelihoodProd_twoQuark_MW400", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_twoQuark_MW400->Sumw2();
  TH1D* h1_QGLikelihoodProd_twoQuark_MW500 = new TH1D("QGLikelihoodProd_twoQuark_MW500", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_twoQuark_MW500->Sumw2();

  TH1D* h1_QGLikelihoodNoPUProd = new TH1D("QGLikelihoodNoPUProd", "", 60, 0., 1.0001);
  h1_QGLikelihoodNoPUProd->Sumw2();
  TH1D* h1_QGLikelihoodProd = new TH1D("QGLikelihoodProd", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd->Sumw2();
  TH1D* h1_QGLikelihoodProd_MW300 = new TH1D("QGLikelihoodProd_MW300", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_MW300->Sumw2();
  TH1D* h1_QGLikelihoodProd_MW400 = new TH1D("QGLikelihoodProd_MW400", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_MW400->Sumw2();
  TH1D* h1_QGLikelihoodProd_MW500 = new TH1D("QGLikelihoodProd_MW500", "", 60, 0., 1.0001);
  h1_QGLikelihoodProd_MW500->Sumw2();


  TH1D* h1_mZjj= new TH1D("mZjj", "", 400, 30., 430.);
  h1_mZjj->Sumw2();
  TH1D* h1_mZjj_MU = new TH1D("mZjj_MU", "", 400, 30., 430.);
  h1_mZjj_MU->Sumw2();
  TH1D* h1_mZjj_ELE= new TH1D("mZjj_ELE", "", 400, 30., 430.);
  h1_mZjj_ELE->Sumw2();
  TH1D* h1_mZjj_nogluetag = new TH1D("mZjj_nogluetag", "", 400, 30., 430.);
  h1_mZjj_nogluetag->Sumw2();
  TH1D* h1_mZjj_nogluetag_MU = new TH1D("mZjj_nogluetag_MU", "", 400, 30., 430.);
  h1_mZjj_nogluetag_MU->Sumw2();
  TH1D* h1_mZjj_nogluetag_ELE= new TH1D("mZjj_nogluetag_ELE", "", 400, 30., 430.);
  h1_mZjj_nogluetag_ELE->Sumw2();
  TH1D* h1_mZjj_0btag= new TH1D("mZjj_0btag", "", 400, 30., 430.);
  h1_mZjj_0btag->Sumw2();
  TH1D* h1_mZjj_1btag= new TH1D("mZjj_1btag", "", 400, 30., 430.);
  h1_mZjj_1btag->Sumw2();
  TH1D* h1_mZjj_2btag= new TH1D("mZjj_2btag", "", 400, 30., 430.);
  h1_mZjj_2btag->Sumw2();
  TH1D* h1_mZjj_loMass= new TH1D("mZjj_loMass", "", 400, 30., 430.);
  h1_mZjj_loMass->Sumw2();
  TH1D* h1_mZjj_medMass= new TH1D("mZjj_medMass", "", 400, 30., 430.);
  h1_mZjj_medMass->Sumw2();
  TH1D* h1_mZjj_hiMass= new TH1D("mZjj_hiMass", "", 400, 30., 430.);
  h1_mZjj_hiMass->Sumw2();
  TH1D* h1_mZjj_loChiSquareProb= new TH1D("mZjj_loChiSquareProb", "", 100, 30., 200.);
  h1_mZjj_loChiSquareProb->Sumw2();
  TH1D* h1_mZjj_hiChiSquareProb= new TH1D("mZjj_hiChiSquareProb", "", 100, 30., 200.);
  h1_mZjj_hiChiSquareProb->Sumw2();

  TH1D* h1_ptZll_presel = new TH1D("ptZll_presel", "", 400, 0., 400.);
  h1_ptZll_presel->Sumw2();
  TH1D* h1_ptZjj_all_presel = new TH1D("ptZjj_all_presel", "", 400, 0., 400.);
  h1_ptZjj_all_presel->Sumw2();

  TH1D* h1_ptZll = new TH1D("ptZll", "", 400, 0., 400.);
  h1_ptZll->Sumw2();
  TH1D* h1_ptZjj = new TH1D("ptZjj", "", 400, 0., 400.);
  h1_ptZjj->Sumw2();

  TH1D* h1_deltaRjj= new TH1D("deltaRjj", "", 450, 0.5, 5.);
  h1_deltaRjj->Sumw2();
  TH1D* h1_deltaRjj_prekin= new TH1D("deltaRjj_prekin", "", 50, 0.5, 5.);
  h1_deltaRjj_prekin->Sumw2();
  TH1D* h1_deltaRjj_sidebands= new TH1D("deltaRjj_sidebands", "", 50, 0.5, 5.);
  h1_deltaRjj_sidebands->Sumw2();
  TH1D* h1_deltaRjj_prekin_sidebands= new TH1D("deltaRjj_prekin_sidebands", "", 50, 0.5, 5.);
  h1_deltaRjj_prekin_sidebands->Sumw2();

  TH1D* h1_cosThetaStar = new TH1D("cosThetaStar", "", 45, -1.001, 1.001);
  h1_cosThetaStar->Sumw2();
  TH1D* h1_cosTheta1 = new TH1D("cosTheta1", "", 45, -1.001, 1.001);
  h1_cosTheta1->Sumw2();
  TH1D* h1_cosTheta2 = new TH1D("cosTheta2", "", 45, 0., 1.001);
  h1_cosTheta2->Sumw2();
  TH1D* h1_phi = new TH1D("phi", "", 45, -3.1416, 3.1416);
  h1_phi->Sumw2();
  TH1D* h1_phi1 = new TH1D("phi1", "", 45, -3.1416, 3.1416);
  h1_phi1->Sumw2();

  TH1D* h1_cosThetaStar_kinfit = new TH1D("cosThetaStar_kinfit", "", 45, -1.0001, 1.001);
  h1_cosThetaStar_kinfit->Sumw2();
  TH1D* h1_cosTheta1_kinfit = new TH1D("cosTheta1_kinfit", "", 45, -1.001, 1.001);
  h1_cosTheta1_kinfit->Sumw2();
  TH1D* h1_cosTheta2_kinfit = new TH1D("cosTheta2_kinfit", "", 45, 0., 1.001);
  h1_cosTheta2_kinfit->Sumw2();
  TH1D* h1_phi_kinfit = new TH1D("phi_kinfit", "", 45, -3.1416, 3.1416);
  h1_phi_kinfit->Sumw2();
  TH1D* h1_phi1_kinfit = new TH1D("phi1_kinfit", "", 45, -3.1416, 3.1416);
  h1_phi1_kinfit->Sumw2();

  
  TH1D* h1_kinfit_chiSquare = new TH1D("kinfit_chiSquare", "", 60, 0., 10.);
  h1_kinfit_chiSquare->Sumw2();
  TH1D* h1_kinfit_chiSquareProb = new TH1D("kinfit_chiSquareProb", "", 60, 0., 1.0001);
  h1_kinfit_chiSquareProb->Sumw2();

  TH1D* h1_helicityLD = new TH1D("helicityLD", "", 60, 0., 1.);
  h1_helicityLD->Sumw2();
  TH1D* h1_helicityLD_nogluetag = new TH1D("helicityLD_nogluetag", "", 60, 0., 1.);
  h1_helicityLD_nogluetag->Sumw2();
  TH1D* h1_helicityLD_nokinfit = new TH1D("helicityLD_nokinfit", "", 60, 0., 1.);
  h1_helicityLD_nokinfit->Sumw2();
  TH1D* h1_helicityLD_sidebands = new TH1D("helicityLD_sidebands", "", 60, 0., 1.);
  h1_helicityLD_sidebands->Sumw2();
  TH1D* h1_helicityLD_MW200 = new TH1D("helicityLD_MW200", "", 100, 0., 1.);
  h1_helicityLD_MW200->Sumw2();
  TH1D* h1_helicityLD_MW250 = new TH1D("helicityLD_MW250", "", 100, 0., 1.);
  h1_helicityLD_MW250->Sumw2();
  TH1D* h1_helicityLD_MW300 = new TH1D("helicityLD_MW300", "", 100, 0., 1.);
  h1_helicityLD_MW300->Sumw2();
  TH1D* h1_helicityLD_MW400 = new TH1D("helicityLD_MW400", "", 100, 0., 1.);
  h1_helicityLD_MW400->Sumw2();
  TH1D* h1_helicityLD_MW500 = new TH1D("helicityLD_MW500", "", 100, 0., 1.);
  h1_helicityLD_MW500->Sumw2();
  TH1D* h1_helicityLD_kinfit = new TH1D("helicityLD_kinfit", "", 60, 0., 1.);
  h1_helicityLD_kinfit->Sumw2();

  TH2D* h2_helicityLD_vs_mZZ = new TH2D("helicityLD_vs_mZZ", "", 50, 100., 600., 50, 0., 1.);
  h2_helicityLD_vs_mZZ->Sumw2();

  TH1D* h1_deltaRZZ= new TH1D("deltaRZZ", "", 60, 0., 6.);
  h1_deltaRZZ->Sumw2();

  TH1D* h1_mZZ_hiChiSquareProb = new TH1D("mZZ_hiChiSquareProb", "", 200, 100., 700.);
  h1_mZZ_hiChiSquareProb->Sumw2();
  TH1D* h1_mZZ_loChiSquareProb = new TH1D("mZZ_loChiSquareProb", "", 200, 100., 700.);
  h1_mZZ_loChiSquareProb->Sumw2();
  TH1D* h1_mZZ_mZjj_cut = new TH1D("mZZ_mZjj_cut", "", 200, 100., 700.);
  h1_mZZ_mZjj_cut->Sumw2();
  TH1D* h1_mZZ_mZjj_notcut = new TH1D("mZZ_mZjj_notcut", "", 200, 100., 700.);
  h1_mZZ_mZjj_notcut->Sumw2();


  TH1D* h1_mZZ_nokinfit_hiMass_all = new TH1D("mZZ_nokinfit_hiMass_all", "", 1500, 150., 1650.);
  h1_mZZ_nokinfit_hiMass_all->Sumw2();

  TH1D* h1_mZZ_kinfit_hiMass_all = new TH1D("mZZ_kinfit_hiMass_all", "",1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_all->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_hiPU = new TH1D("mZZ_kinfit_hiMass_hiPU", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_hiPU->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_loPU = new TH1D("mZZ_kinfit_hiMass_loPU", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_loPU->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_nogluetag = new TH1D("mZZ_kinfit_hiMass_nogluetag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_nogluetag->Sumw2();

  TH1D* h1_mZZ_kinfit_hiMass_oneGluon = new TH1D("mZZ_kinfit_hiMass_oneGluon", "",1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_oneGluon->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_twoGluon = new TH1D("mZZ_kinfit_hiMass_twoGluon", "",1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_twoGluon->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_twoQuark = new TH1D("mZZ_kinfit_hiMass_twoQuark", "",1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_twoQuark->Sumw2();

  TH1D* h1_mZZ_kinfit_hiMass_gluetag= new TH1D("mZZ_kinfit_hiMass_gluetag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_gluetag->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_0btag= new TH1D("mZZ_kinfit_hiMass_0btag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_0btag->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_1btag= new TH1D("mZZ_kinfit_hiMass_1btag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_1btag->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_2btag= new TH1D("mZZ_kinfit_hiMass_2btag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_2btag->Sumw2();

  TH1D* h1_mZZ_kinfit_hiMass_gluetag_ELE = new TH1D("mZZ_kinfit_hiMass_gluetag_ELE ", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_gluetag_ELE->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_0btag_ELE = new TH1D("mZZ_kinfit_hiMass_0btag_ELE ", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_0btag_ELE->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_1btag_ELE = new TH1D("mZZ_kinfit_hiMass_1btag_ELE ", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_1btag_ELE->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_2btag_ELE = new TH1D("mZZ_kinfit_hiMass_2btag_ELE ", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_2btag_ELE->Sumw2();

  TH1D* h1_mZZ_kinfit_hiMass_gluetag_MU = new TH1D("mZZ_kinfit_hiMass_gluetag_MU ", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_gluetag_MU->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_0btag_MU = new TH1D("mZZ_kinfit_hiMass_0btag_MU ", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_0btag_MU->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_1btag_MU = new TH1D("mZZ_kinfit_hiMass_1btag_MU ", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_1btag_MU->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_2btag_MU = new TH1D("mZZ_kinfit_hiMass_2btag_MU ", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_2btag_MU->Sumw2();

  TH1D* h1_mZZ_kinfit_hiMass_sidebands_gluetag= new TH1D("mZZ_kinfit_hiMass_sidebands_gluetag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_gluetag->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_0btag= new TH1D("mZZ_kinfit_hiMass_sidebands_0btag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_0btag->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_1btag= new TH1D("mZZ_kinfit_hiMass_sidebands_1btag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_1btag->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_2btag= new TH1D("mZZ_kinfit_hiMass_sidebands_2btag", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_2btag->Sumw2();

  TH1D* h1_mZZ_kinfit_hiMass_sidebands_gluetag_ELE= new TH1D("mZZ_kinfit_hiMass_sidebands_gluetag_ELE", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_gluetag_ELE->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_0btag_ELE= new TH1D("mZZ_kinfit_hiMass_sidebands_0btag_ELE", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_0btag_ELE->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_1btag_ELE= new TH1D("mZZ_kinfit_hiMass_sidebands_1btag_ELE", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_1btag_ELE->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_2btag_ELE= new TH1D("mZZ_kinfit_hiMass_sidebands_2btag_ELE", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_2btag_ELE->Sumw2();

  TH1D* h1_mZZ_kinfit_hiMass_sidebands_gluetag_MU= new TH1D("mZZ_kinfit_hiMass_sidebands_gluetag_MU", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_gluetag_MU->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_0btag_MU= new TH1D("mZZ_kinfit_hiMass_sidebands_0btag_MU", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_0btag_MU->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_1btag_MU= new TH1D("mZZ_kinfit_hiMass_sidebands_1btag_MU", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_1btag_MU->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_sidebands_2btag_MU= new TH1D("mZZ_kinfit_hiMass_sidebands_2btag_MU", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_sidebands_2btag_MU->Sumw2();

/*
  TH1D* h1_mZZ_nokinfit_hiMass_gluetag= new TH1D("mZZ_nokinfit_hiMass_gluetag", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_gluetag->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_0btag= new TH1D("mZZ_nokinfit_hiMass_0btag", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_0btag->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_1btag= new TH1D("mZZ_nokinfit_hiMass_1btag", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_1btag->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_2btag= new TH1D("mZZ_nokinfit_hiMass_2btag", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_2btag->Sumw2();

  TH1D* h1_mZZ_nokinfit_hiMass_gluetag_ELE = new TH1D("mZZ_nokinfit_hiMass_gluetag_ELE ", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_gluetag_ELE->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_0btag_ELE = new TH1D("mZZ_nokinfit_hiMass_0btag_ELE ", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_0btag_ELE->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_1btag_ELE = new TH1D("mZZ_nokinfit_hiMass_1btag_ELE ", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_1btag_ELE->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_2btag_ELE = new TH1D("mZZ_nokinfit_hiMass_2btag_ELE ", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_2btag_ELE->Sumw2();

  TH1D* h1_mZZ_nokinfit_hiMass_gluetag_MU = new TH1D("mZZ_nokinfit_hiMass_gluetag_MU ", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_gluetag_MU->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_0btag_MU = new TH1D("mZZ_nokinfit_hiMass_0btag_MU ", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_0btag_MU->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_1btag_MU = new TH1D("mZZ_nokinfit_hiMass_1btag_MU ", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_1btag_MU->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_2btag_MU = new TH1D("mZZ_nokinfit_hiMass_2btag_MU ", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_2btag_MU->Sumw2();

  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_gluetag= new TH1D("mZZ_nokinfit_hiMass_sidebands_gluetag", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_gluetag->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_0btag= new TH1D("mZZ_nokinfit_hiMass_sidebands_0btag", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_0btag->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_1btag= new TH1D("mZZ_nokinfit_hiMass_sidebands_1btag", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_1btag->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_2btag= new TH1D("mZZ_nokinfit_hiMass_sidebands_2btag", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_2btag->Sumw2();

  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_gluetag_ELE= new TH1D("mZZ_nokinfit_hiMass_sidebands_gluetag_ELE", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_gluetag_ELE->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_0btag_ELE= new TH1D("mZZ_nokinfit_hiMass_sidebands_0btag_ELE", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_0btag_ELE->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_1btag_ELE= new TH1D("mZZ_nokinfit_hiMass_sidebands_1btag_ELE", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_1btag_ELE->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_2btag_ELE= new TH1D("mZZ_nokinfit_hiMass_sidebands_2btag_ELE", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_2btag_ELE->Sumw2();

  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_gluetag_MU= new TH1D("mZZ_nokinfit_hiMass_sidebands_gluetag_MU", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_gluetag_MU->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_0btag_MU= new TH1D("mZZ_nokinfit_hiMass_sidebands_0btag_MU", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_0btag_MU->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_1btag_MU= new TH1D("mZZ_nokinfit_hiMass_sidebands_1btag_MU", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_1btag_MU->Sumw2();
  TH1D* h1_mZZ_nokinfit_hiMass_sidebands_2btag_MU= new TH1D("mZZ_nokinfit_hiMass_sidebands_2btag_MU", "", 1250, 150., 1400.);
  h1_mZZ_nokinfit_hiMass_sidebands_2btag_MU->Sumw2();
*/


  TH1D* h1_mZZ_kinfit_hiMass_hiQG= new TH1D("mZZ_kinfit_hiMass_hiQG", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_hiQG->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_loQG= new TH1D("mZZ_kinfit_hiMass_loQG", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_loQG->Sumw2();
  TH1D* h1_mZZ_ZjjMassConstr_hiMass  = new TH1D("mZZ_ZjjMassConstr_hiMass", "", 1500, 150., 1650.);
  h1_mZZ_ZjjMassConstr_hiMass->Sumw2();

  TH1D* h1_deltaRZmatching = new TH1D("deltaRZmatching", "", 50, 0., 1.0);
  h1_deltaRZmatching->Sumw2();
  TH1D* h1_mZZ_kinfit_hiMass_0btag_matched = new TH1D("mZZ_kinfit_hiMass_0btag_matched", "", 1500, 150., 1650.);
  h1_mZZ_kinfit_hiMass_0btag_matched->Sumw2();

  TH1D* h1_ptZZ  = new TH1D("ptZZ", "", 100, 0., 300.);
  h1_ptZZ->Sumw2();
  TH1D* h1_ptZZ_kinfit  = new TH1D("ptZZ_kinfit", "", 100, 0., 300.);
  h1_ptZZ_kinfit->Sumw2();
  TH1D* h1_etaZZ  = new TH1D("etaZZ", "", 100, -5.5, 5.5);
  h1_etaZZ->Sumw2();
  TH1D* h1_etaZZ_kinfit  = new TH1D("etaZZ_kinfit", "", 100, -5.5, 5.5);
  h1_etaZZ_kinfit->Sumw2();
  

//TH1D* h1_mZZ_MCassoc  = new TH1D("mZZ_MCassoc", "", 100, 200., 600.);
//h1_mZZ_MCassoc->Sumw2();
//TH1D* h1_mZZ_MCassoc_ZjjMassConstr  = new TH1D("mZZ_MCassoc_ZjjMassConstr", "", 100, 200., 600.);
//h1_mZZ_MCassoc_ZjjMassConstr->Sumw2();
//TH1D* h1_mZZ_MCassoc_kinfit  = new TH1D("mZZ_MCassoc_kinfit", "", 100, 200., 600.);
//h1_mZZ_MCassoc_kinfit->Sumw2();
//TH1D* h1_mZZ_MCassoc_kinfit_cands  = new TH1D("mZZ_MCassoc_kinfit_cands", "", 100, 200., 600.);
//h1_mZZ_MCassoc_kinfit_cands->Sumw2();

  TH1D* h1_partFlavor_tight = new TH1D("partFlavor_tight", "", 31, -8.5, 22.5);

  TH2D* h2_mZjj_vs_mZZ = new TH2D("mZjj_vs_mZZ", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_0btag = new TH2D("mZjj_vs_mZZ_0btag", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_0btag->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_1btag = new TH2D("mZjj_vs_mZZ_1btag", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_1btag->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_2btag = new TH2D("mZjj_vs_mZZ_2btag", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_2btag->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_gluetag = new TH2D("mZjj_vs_mZZ_gluetag", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_gluetag->Sumw2();

  TH2D* h2_mZjj_vs_mZZ_kinfit = new TH2D("mZjj_vs_mZZ_kinfit", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_kinfit->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_kinfit_0btag = new TH2D("mZjj_vs_mZZ_kinfit_0btag", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_kinfit_0btag->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_kinfit_1btag = new TH2D("mZjj_vs_mZZ_kinfit_1btag", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_kinfit_1btag->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_kinfit_2btag = new TH2D("mZjj_vs_mZZ_kinfit_2btag", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_kinfit_2btag->Sumw2();
  TH2D* h2_mZjj_vs_mZZ_kinfit_gluetag = new TH2D("mZjj_vs_mZZ_kinfit_gluetag", "", 600, 150., 750., 200, 40., 240.);
  h2_mZjj_vs_mZZ_kinfit_gluetag->Sumw2();


  TH1D* h1_deltaE_ch = new TH1D("deltaE_ch", "", 100, -1., 1.);
  TH1D* h1_deltaE_gamma = new TH1D("deltaE_gamma", "", 100, -1., 1.);
  TH1D* h1_deltaE_nh = new TH1D("deltaE_nh", "", 100, -1., 1.);

  TH1D* h1_deltaEta_ch    = new TH1D("deltaEta_ch", "", 100, -1., 1.);
  TH1D* h1_deltaEta_gamma = new TH1D("deltaEta_gamma", "", 100, -1., 1.);
  TH1D* h1_deltaEta_nh    = new TH1D("deltaEta_nh", "", 100, -1., 1.);

  TH1D* h1_deltaPhi_ch    = new TH1D("deltaPhi_ch", "", 100, -1., 1.);
  TH1D* h1_deltaPhi_gamma = new TH1D("deltaPhi_gamma", "", 100, -1., 1.);
  TH1D* h1_deltaPhi_nh    = new TH1D("deltaPhi_nh", "", 100, -1., 1.);

  TH1D* h1_deltaPt_ch = new TH1D("deltaPt_ch", "", 100, -1., 1.);
  TH1D* h1_deltaPt_gamma = new TH1D("deltaPt_gamma", "", 100, -1., 1.);
  TH1D* h1_deltaPt_nh = new TH1D("deltaPt_nh", "", 100, -1., 1.);





  Int_t nPU;
  tree_->SetBranchAddress("nPU", &nPU);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);
  Int_t LS;
  tree_->SetBranchAddress("LS", &LS);
  unsigned int event;
  tree_->SetBranchAddress("event", &event);
  Float_t genWeight=1.;
  if( dataset_tstr.Contains("sherpa") )
    tree_->SetBranchAddress("genWeight", &genWeight);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);
  Float_t eventWeightPU;
  tree_->SetBranchAddress("eventWeightPU", &eventWeightPU);
  Float_t eventWeightPU_ave;
  tree_->SetBranchAddress("eventWeightPU_ave", &eventWeightPU_ave);
  Float_t eventWeight_Zee;
  tree_->SetBranchAddress("eventWeight_Zee", &eventWeight_Zee);
  Float_t eventWeight_Zmm;
  tree_->SetBranchAddress("eventWeight_Zmm", &eventWeight_Zmm);

  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Float_t pfMet;
  tree_->SetBranchAddress("epfMet", &pfMet);
  Float_t metSignificance;
  tree_->SetBranchAddress("metSignificance", &metSignificance);
  Float_t mEtSig;
  tree_->SetBranchAddress("mEtSig", &mEtSig);
  Float_t phiMet;
  tree_->SetBranchAddress("phipfMet", &phiMet);


  int leptType;
  tree_->SetBranchAddress("leptType", &leptType);

  Float_t eLept1;
  tree_->SetBranchAddress("eLept1", &eLept1);
  Float_t ptLept1;
  tree_->SetBranchAddress("ptLept1", &ptLept1);
  Float_t etaLept1;
  tree_->SetBranchAddress("etaLept1", &etaLept1);
  Float_t phiLept1;
  tree_->SetBranchAddress("phiLept1", &phiLept1);
  Int_t chargeLept1;
  tree_->SetBranchAddress("chargeLept1", &chargeLept1);

  Float_t eLept2;
  tree_->SetBranchAddress("eLept2", &eLept2);
  Float_t ptLept2;
  tree_->SetBranchAddress("ptLept2", &ptLept2);
  Float_t etaLept2;
  tree_->SetBranchAddress("etaLept2", &etaLept2);
  Float_t phiLept2;
  tree_->SetBranchAddress("phiLept2", &phiLept2);
  Int_t chargeLept2;
  tree_->SetBranchAddress("chargeLept2", &chargeLept2);


  Int_t nPairs;
  tree_->SetBranchAddress("nPairs", &nPairs);

  Int_t iJet1[50];
  tree_->SetBranchAddress("iJet1", iJet1);
  Float_t eJet1[50];
  tree_->SetBranchAddress("eJet1", eJet1);
  Float_t ptJet1[50];
  tree_->SetBranchAddress("ptJet1", ptJet1);
  Float_t etaJet1[50];
  tree_->SetBranchAddress("etaJet1", etaJet1);
  Float_t phiJet1[50];
  tree_->SetBranchAddress("phiJet1", phiJet1);
  Float_t eJet1Gen[50];
  tree_->SetBranchAddress("eJet1Gen", eJet1Gen);
  Float_t ptJet1Gen[50];
  tree_->SetBranchAddress("ptJet1Gen", ptJet1Gen);
  Float_t etaJet1Gen[50];
  tree_->SetBranchAddress("etaJet1Gen", etaJet1Gen);
  Float_t phiJet1Gen[50];
  tree_->SetBranchAddress("phiJet1Gen", phiJet1Gen);
  Float_t eChargedHadronsJet1[50];
  tree_->SetBranchAddress("eChargedHadronsJet1", eChargedHadronsJet1);
  Float_t rmsCandJet1[50];
  tree_->SetBranchAddress("rmsCandJet1", rmsCandJet1);
  Float_t ptDJet1[50];
  tree_->SetBranchAddress("ptDJet1", ptDJet1);
  Int_t nChargedJet1[50];
  tree_->SetBranchAddress("nChargedJet1", nChargedJet1);
  Int_t nNeutralJet1[50];
  tree_->SetBranchAddress("nNeutralJet1", nNeutralJet1);
  Float_t eMuonsJet1[50];
  tree_->SetBranchAddress("eMuonsJet1", eMuonsJet1);
  Float_t eElectronsJet1[50];
  tree_->SetBranchAddress("eElectronsJet1", eElectronsJet1);
  Float_t trackCountingHighEffBJetTagJet1[50];
  tree_->SetBranchAddress("trackCountingHighEffBJetTagJet1", trackCountingHighEffBJetTagJet1);
  Float_t trackCountingHighPurBJetTagJet1[50];
  tree_->SetBranchAddress("trackCountingHighPurBJetTagJet1", trackCountingHighPurBJetTagJet1);
  Float_t simpleSecondaryVertexHighEffBJetTagJet1[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagJet1", simpleSecondaryVertexHighEffBJetTagJet1);
  Float_t simpleSecondaryVertexHighPurBJetTagJet1[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagJet1", simpleSecondaryVertexHighPurBJetTagJet1);
  Float_t jetBProbabilityBJetTagJet1[50];
  tree_->SetBranchAddress("jetBProbabilityBJetTagJet1", jetBProbabilityBJetTagJet1);
  Float_t jetProbabilityBJetTagJet1[50];
  tree_->SetBranchAddress("jetProbabilityBJetTagJet1", jetProbabilityBJetTagJet1);

  Int_t nPFCand1;
  tree_->SetBranchAddress("nPFCand1", &nPFCand1);
  Float_t ePFCand1[100];
  tree_->SetBranchAddress("ePFCand1", ePFCand1);
  Float_t ptPFCand1[100];
  tree_->SetBranchAddress("ptPFCand1", ptPFCand1);
  Float_t etaPFCand1[100];
  tree_->SetBranchAddress("etaPFCand1", etaPFCand1);
  Float_t phiPFCand1[100];
  tree_->SetBranchAddress("phiPFCand1", phiPFCand1);
  Int_t particleTypePFCand1[100];
  tree_->SetBranchAddress("particleTypePFCand1", particleTypePFCand1);

  Int_t iJet2[50];
  tree_->SetBranchAddress("iJet2", iJet2);
  Float_t eJet2[50];
  tree_->SetBranchAddress("eJet2", eJet2);
  Float_t ptJet2[50];
  tree_->SetBranchAddress("ptJet2", ptJet2);
  Float_t etaJet2[50];
  tree_->SetBranchAddress("etaJet2", etaJet2);
  Float_t phiJet2[50];
  tree_->SetBranchAddress("phiJet2", phiJet2);
  Float_t eJet2Gen[50];
  tree_->SetBranchAddress("eJet2Gen", eJet2Gen);
  Float_t ptJet2Gen[50];
  tree_->SetBranchAddress("ptJet2Gen", ptJet2Gen);
  Float_t etaJet2Gen[50];
  tree_->SetBranchAddress("etaJet2Gen", etaJet2Gen);
  Float_t phiJet2Gen[50];
  tree_->SetBranchAddress("phiJet2Gen", phiJet2Gen);
  Float_t eChargedHadronsJet2[50];
  tree_->SetBranchAddress("eChargedHadronsJet2", eChargedHadronsJet2);
  Float_t rmsCandJet2[50];
  tree_->SetBranchAddress("rmsCandJet2", rmsCandJet2);
  Float_t ptDJet2[50];
  tree_->SetBranchAddress("ptDJet2", ptDJet2);
  Int_t nChargedJet2[50];
  tree_->SetBranchAddress("nChargedJet2", nChargedJet2);
  Int_t nNeutralJet2[50];
  tree_->SetBranchAddress("nNeutralJet2", nNeutralJet2);
  Float_t eMuonsJet2[50];
  tree_->SetBranchAddress("eMuonsJet2", eMuonsJet2);
  Float_t eElectronsJet2[50];
  tree_->SetBranchAddress("eElectronsJet2", eElectronsJet2);
  Float_t trackCountingHighEffBJetTagJet2[50];
  tree_->SetBranchAddress("trackCountingHighEffBJetTagJet2", trackCountingHighEffBJetTagJet2);
  Float_t trackCountingHighPurBJetTagJet2[50];
  tree_->SetBranchAddress("trackCountingHighPurBJetTagJet2", trackCountingHighPurBJetTagJet2);
  Float_t simpleSecondaryVertexHighEffBJetTagJet2[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagJet2", simpleSecondaryVertexHighEffBJetTagJet2);
  Float_t simpleSecondaryVertexHighPurBJetTagJet2[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagJet2", simpleSecondaryVertexHighPurBJetTagJet2);
  Float_t jetBProbabilityBJetTagJet2[50];
  tree_->SetBranchAddress("jetBProbabilityBJetTagJet2", jetBProbabilityBJetTagJet2);
  Float_t jetProbabilityBJetTagJet2[50];
  tree_->SetBranchAddress("jetProbabilityBJetTagJet2", jetProbabilityBJetTagJet2);

  Float_t ptJetRecoil[50];
  tree_->SetBranchAddress("ptJetRecoil", ptJetRecoil);
  Float_t etaJetRecoil[50];
  tree_->SetBranchAddress("etaJetRecoil", etaJetRecoil);
  Float_t phiJetRecoil[50];
  tree_->SetBranchAddress("phiJetRecoil", phiJetRecoil);
  Float_t eJetRecoil[50];
  tree_->SetBranchAddress("eJetRecoil", eJetRecoil);
  Float_t ptDJetRecoil[50];
  tree_->SetBranchAddress("ptDJetRecoil", ptDJetRecoil);
  Int_t nChargedJetRecoil[50];
  tree_->SetBranchAddress("nChargedJetRecoil", nChargedJetRecoil);
  Int_t nNeutralJetRecoil[50];
  tree_->SetBranchAddress("nNeutralJetRecoil", nNeutralJetRecoil);

  Int_t nPFCand2;
  tree_->SetBranchAddress("nPFCand2", &nPFCand2);
  Float_t ePFCand2[100];
  tree_->SetBranchAddress("ePFCand2", ePFCand2);
  Float_t ptPFCand2[100];
  tree_->SetBranchAddress("ptPFCand2", ptPFCand2);
  Float_t etaPFCand2[100];
  tree_->SetBranchAddress("etaPFCand2", etaPFCand2);
  Float_t phiPFCand2[100];
  tree_->SetBranchAddress("phiPFCand2", phiPFCand2);
  Int_t particleTypePFCand2[100];
  tree_->SetBranchAddress("particleTypePFCand2", particleTypePFCand2);

  Int_t nPart;
  tree_->SetBranchAddress("nPart", &nPart);
  Float_t ePart[20];
  tree_->SetBranchAddress("ePart", ePart);
  Float_t ptPart[20];
  tree_->SetBranchAddress("ptPart", ptPart);
  Float_t etaPart[20];
  tree_->SetBranchAddress("etaPart", etaPart);
  Float_t phiPart[20];
  tree_->SetBranchAddress("phiPart", phiPart);
  Int_t pdgIdPart[20];
  tree_->SetBranchAddress("pdgIdPart", pdgIdPart);


  // HLT:
  Bool_t passed_HLT_DoubleMu6;
  tree_->SetBranchAddress("passed_HLT_DoubleMu6", &passed_HLT_DoubleMu6);
  Bool_t passed_HLT_DoubleMu7;
  tree_->SetBranchAddress("passed_HLT_DoubleMu7", &passed_HLT_DoubleMu7);
  Bool_t passed_HLT_Mu13_Mu8;
  tree_->SetBranchAddress("passed_HLT_Mu13_Mu8", &passed_HLT_Mu13_Mu8);
  Bool_t passed_HLT_IsoMu17;
  tree_->SetBranchAddress("passed_HLT_IsoMu17", &passed_HLT_IsoMu17);
  Bool_t passed_HLT_IsoMu24;
  tree_->SetBranchAddress("passed_HLT_IsoMu24", &passed_HLT_IsoMu24);
  Bool_t passed_HLT_Mu8_Jet40;
  tree_->SetBranchAddress("passed_HLT_Mu8_Jet40", &passed_HLT_Mu8_Jet40);
  Bool_t passed_HLT_L2DoubleMu23_NoVertex;
  tree_->SetBranchAddress("passed_HLT_L2DoubleMu23_NoVertex", &passed_HLT_L2DoubleMu23_NoVertex);
  Bool_t passed_HLT_L2DoubleMu30_NoVertex;
  tree_->SetBranchAddress("passed_HLT_L2DoubleMu30_NoVertex", &passed_HLT_L2DoubleMu30_NoVertex);
  Bool_t passed_HLT_TripleMu5;
  tree_->SetBranchAddress("passed_HLT_TripleMu5", &passed_HLT_TripleMu5);

  Bool_t passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
  tree_->SetBranchAddress("passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL", &passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
  Bool_t passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
  tree_->SetBranchAddress("passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
  


  float nEventsPassed_fb_kinfit=0.;
  int nEventsPassed_kinfit=0;

  float nEventsPassed_fb_gluetag_250=0.;
  int nEventsPassed_gluetag_250=0;
  float nEventsPassed_fb_0btag_250=0.;
  int nEventsPassed_0btag_250=0;
  float nEventsPassed_fb_1btag_250=0.;
  int nEventsPassed_1btag_250=0;
  float nEventsPassed_fb_2btag_250=0.;
  int nEventsPassed_2btag_250=0;

  float nEventsPassed_fb_gluetag_250_ELE=0.;
  int nEventsPassed_gluetag_250_ELE=0;
  float nEventsPassed_fb_0btag_250_ELE=0.;
  int nEventsPassed_0btag_250_ELE=0;
  float nEventsPassed_fb_1btag_250_ELE=0.;
  int nEventsPassed_1btag_250_ELE=0;
  float nEventsPassed_fb_2btag_250_ELE=0.;
  int nEventsPassed_2btag_250_ELE=0;

  float nEventsPassed_fb_gluetag_250_MU=0.;
  int nEventsPassed_gluetag_250_MU=0;
  float nEventsPassed_fb_0btag_250_MU=0.;
  int nEventsPassed_0btag_250_MU=0;
  float nEventsPassed_fb_1btag_250_MU=0.;
  int nEventsPassed_1btag_250_MU=0;
  float nEventsPassed_fb_2btag_250_MU=0.;
  int nEventsPassed_2btag_250_MU=0;

  float nEventsPassed_fb_gluetag_300=0.;
  int nEventsPassed_gluetag_300=0;
  float nEventsPassed_fb_0btag_300=0.;
  int nEventsPassed_0btag_300=0;
  float nEventsPassed_fb_1btag_300=0.;
  int nEventsPassed_1btag_300=0;
  float nEventsPassed_fb_2btag_300=0.;
  int nEventsPassed_2btag_300=0;

  float nEventsPassed_fb_gluetag_300_ELE=0.;
  int nEventsPassed_gluetag_300_ELE=0;
  float nEventsPassed_fb_0btag_300_ELE=0.;
  int nEventsPassed_0btag_300_ELE=0;
  float nEventsPassed_fb_1btag_300_ELE=0.;
  int nEventsPassed_1btag_300_ELE=0;
  float nEventsPassed_fb_2btag_300_ELE=0.;
  int nEventsPassed_2btag_300_ELE=0;

  float nEventsPassed_fb_gluetag_300_MU=0.;
  int nEventsPassed_gluetag_300_MU=0;
  float nEventsPassed_fb_0btag_300_MU=0.;
  int nEventsPassed_0btag_300_MU=0;
  float nEventsPassed_fb_1btag_300_MU=0.;
  int nEventsPassed_1btag_300_MU=0;
  float nEventsPassed_fb_2btag_300_MU=0.;
  int nEventsPassed_2btag_300_MU=0;

  float nEventsPassed_fb_gluetag_350=0.;
  int nEventsPassed_gluetag_350=0;
  float nEventsPassed_fb_0btag_350=0.;
  int nEventsPassed_0btag_350=0;
  float nEventsPassed_fb_1btag_350=0.;
  int nEventsPassed_1btag_350=0;
  float nEventsPassed_fb_2btag_350=0.;
  int nEventsPassed_2btag_350=0;

  float nEventsPassed_fb_gluetag_350_ELE=0.;
  int nEventsPassed_gluetag_350_ELE=0;
  float nEventsPassed_fb_0btag_350_ELE=0.;
  int nEventsPassed_0btag_350_ELE=0;
  float nEventsPassed_fb_1btag_350_ELE=0.;
  int nEventsPassed_1btag_350_ELE=0;
  float nEventsPassed_fb_2btag_350_ELE=0.;
  int nEventsPassed_2btag_350_ELE=0;

  float nEventsPassed_fb_gluetag_350_MU=0.;
  int nEventsPassed_gluetag_350_MU=0;
  float nEventsPassed_fb_0btag_350_MU=0.;
  int nEventsPassed_0btag_350_MU=0;
  float nEventsPassed_fb_1btag_350_MU=0.;
  int nEventsPassed_1btag_350_MU=0;
  float nEventsPassed_fb_2btag_350_MU=0.;
  int nEventsPassed_2btag_350_MU=0;

  float nEventsPassed_fb_gluetag_400=0.;
  int nEventsPassed_gluetag_400=0;
  float nEventsPassed_fb_0btag_400=0.;
  int nEventsPassed_0btag_400=0;
  float nEventsPassed_fb_1btag_400=0.;
  int nEventsPassed_1btag_400=0;
  float nEventsPassed_fb_2btag_400=0.;
  int nEventsPassed_2btag_400=0;

  float nEventsPassed_fb_gluetag_400_ELE=0.;
  int nEventsPassed_gluetag_400_ELE=0;
  float nEventsPassed_fb_0btag_400_ELE=0.;
  int nEventsPassed_0btag_400_ELE=0;
  float nEventsPassed_fb_1btag_400_ELE=0.;
  int nEventsPassed_1btag_400_ELE=0;
  float nEventsPassed_fb_2btag_400_ELE=0.;
  int nEventsPassed_2btag_400_ELE=0;

  float nEventsPassed_fb_gluetag_400_MU=0.;
  int nEventsPassed_gluetag_400_MU=0;
  float nEventsPassed_fb_0btag_400_MU=0.;
  int nEventsPassed_0btag_400_MU=0;
  float nEventsPassed_fb_1btag_400_MU=0.;
  int nEventsPassed_1btag_400_MU=0;
  float nEventsPassed_fb_2btag_400_MU=0.;
  int nEventsPassed_2btag_400_MU=0;

  float nEventsPassed_fb_gluetag_450=0.;
  int nEventsPassed_gluetag_450=0;
  float nEventsPassed_fb_0btag_450=0.;
  int nEventsPassed_0btag_450=0;
  float nEventsPassed_fb_1btag_450=0.;
  int nEventsPassed_1btag_450=0;
  float nEventsPassed_fb_2btag_450=0.;
  int nEventsPassed_2btag_450=0;

  float nEventsPassed_fb_gluetag_450_ELE=0.;
  int nEventsPassed_gluetag_450_ELE=0;
  float nEventsPassed_fb_0btag_450_ELE=0.;
  int nEventsPassed_0btag_450_ELE=0;
  float nEventsPassed_fb_1btag_450_ELE=0.;
  int nEventsPassed_1btag_450_ELE=0;
  float nEventsPassed_fb_2btag_450_ELE=0.;
  int nEventsPassed_2btag_450_ELE=0;

  float nEventsPassed_fb_gluetag_450_MU=0.;
  int nEventsPassed_gluetag_450_MU=0;
  float nEventsPassed_fb_0btag_450_MU=0.;
  int nEventsPassed_0btag_450_MU=0;
  float nEventsPassed_fb_1btag_450_MU=0.;
  int nEventsPassed_1btag_450_MU=0;
  float nEventsPassed_fb_2btag_450_MU=0.;
  int nEventsPassed_2btag_450_MU=0;

  float nEventsPassed_fb_gluetag_500=0.;
  int nEventsPassed_gluetag_500=0;
  float nEventsPassed_fb_0btag_500=0.;
  int nEventsPassed_0btag_500=0;
  float nEventsPassed_fb_1btag_500=0.;
  int nEventsPassed_1btag_500=0;
  float nEventsPassed_fb_2btag_500=0.;
  int nEventsPassed_2btag_500=0;

  float nEventsPassed_fb_gluetag_500_ELE=0.;
  int nEventsPassed_gluetag_500_ELE=0;
  float nEventsPassed_fb_0btag_500_ELE=0.;
  int nEventsPassed_0btag_500_ELE=0;
  float nEventsPassed_fb_1btag_500_ELE=0.;
  int nEventsPassed_1btag_500_ELE=0;
  float nEventsPassed_fb_2btag_500_ELE=0.;
  int nEventsPassed_2btag_500_ELE=0;

  float nEventsPassed_fb_gluetag_500_MU=0.;
  int nEventsPassed_gluetag_500_MU=0;
  float nEventsPassed_fb_0btag_500_MU=0.;
  int nEventsPassed_0btag_500_MU=0;
  float nEventsPassed_fb_1btag_500_MU=0.;
  int nEventsPassed_1btag_500_MU=0;
  float nEventsPassed_fb_2btag_500_MU=0.;
  int nEventsPassed_2btag_500_MU=0;

  float nEventsPassed_fb_gluetag_600=0.;
  int nEventsPassed_gluetag_600=0;
  float nEventsPassed_fb_0btag_600=0.;
  int nEventsPassed_0btag_600=0;
  float nEventsPassed_fb_1btag_600=0.;
  int nEventsPassed_1btag_600=0;
  float nEventsPassed_fb_2btag_600=0.;
  int nEventsPassed_2btag_600=0;

  float nEventsPassed_fb_gluetag_600_ELE=0.;
  int nEventsPassed_gluetag_600_ELE=0;
  float nEventsPassed_fb_0btag_600_ELE=0.;
  int nEventsPassed_0btag_600_ELE=0;
  float nEventsPassed_fb_1btag_600_ELE=0.;
  int nEventsPassed_1btag_600_ELE=0;
  float nEventsPassed_fb_2btag_600_ELE=0.;
  int nEventsPassed_2btag_600_ELE=0;

  float nEventsPassed_fb_gluetag_600_MU=0.;
  int nEventsPassed_gluetag_600_MU=0;
  float nEventsPassed_fb_0btag_600_MU=0.;
  int nEventsPassed_0btag_600_MU=0;
  float nEventsPassed_fb_1btag_600_MU=0.;
  int nEventsPassed_1btag_600_MU=0;
  float nEventsPassed_fb_2btag_600_MU=0.;
  int nEventsPassed_2btag_600_MU=0;



  int nEntries = tree_->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_1_3/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root");
  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_1_3/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1.root");
  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1.root");
  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");
  float Zmass = 91.1876;
  DiJetKinFitter* fitter_jets = new DiJetKinFitter( "fitter_jets", "fitter_jets", Zmass );
  HelicityLikelihoodDiscriminant *LD = new HelicityLikelihoodDiscriminant();
  float helicityLD_selected = -1.;
  float helicityLD_nokinfit_selected = -1.;
  HelicityLikelihoodDiscriminant::HelicityAngles hangles_selected;

  BTagSFUtil* btsfutil = new BTagSFUtil(13);
  TFile* file_loose = TFile::Open("BTagPayloads_TCHEL.root");
  TFile* file_medium = TFile::Open("BTagPayloads_TCHEM.root");
  btsfutil->set_fileLoose(file_loose);
  btsfutil->set_fileMedium(file_medium);

  std::string puType = "Spring11_Flat10";
  std::string puType_ave = "Spring11_Flat10";
  if( dataset_tstr.Contains("Summer11") && dataset_tstr.Contains("PU_S4") ) {
    puType = "Summer11_S4";
    puType_ave = "Summer11_S4_ave";
  } else if( dataset_tstr.Contains("Fall11") ) {
    puType = "Fall11";
  }
  PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
  PUWeight* fPUWeight_ave = new PUWeight(-1, "2011A", puType_ave);
  std::string puFileName;
  if( PUType_=="HR11" || PUType_=="HR11_v3")
    puFileName = "Pileup_DATA_up_to_178479_SiXie.root";
    //puFileName = "Pileup_DATA_up_to_178078.root";
  else if( PUType_=="Run2011A" )
    puFileName = "Pileup_DATA_up_to_173692.root";
  else if( PUType_=="Run2011A_73pb" )
    puFileName = "all2011A.pileup_v2_73mb.root";
  else if( PUType_=="Run2011B" )
    puFileName = "Pileup_DATA_Run2011B.root";
    //puFileName = "Pileup_DATA_173692_to_178078.root";
  else if( PUType_=="Run2011B_73pb" )
    puFileName = "all2011B.pileup_v2_73mb.root";
  else if( PUType_=="HR11_73pb" || PUType_=="HR11_73pb_DY" )
    puFileName = "all2011AB.pileup_v2_73mb.root";
  else if( PUType_!="HR11_v2" ) {
    std::cout << "-> Unknown PU Type: '" << PUType_ << "'. Will use HR11 default." << std::endl;
    puFileName = "Pileup_DATA_up_to_178078.root";
  }


  if( PUType_!="HR11_v2" ) {
    std::cout << std::endl << "-> Using data pileup file: " << puFileName << std::endl;
    TFile* filePU = TFile::Open(puFileName.c_str());
    TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
    fPUWeight->SetDataHistogram(h1_nPU_data);
    fPUWeight_ave->SetDataHistogram(h1_nPU_data);
  } else {  // HR11_v2: 4.6fb-1 = 2.1 (A) + 2.5 (B)
    TFile* filePU_RunA = TFile::Open("all2011A.pileup_v2_73mb.root");
    TFile* filePU_RunB = TFile::Open("all2011B.pileup_v2_73mb.root");
    TH1F* h1_PURunA = (TH1F*)filePU_RunA->Get("pileup");
    TH1F* h1_PURunB = (TH1F*)filePU_RunB->Get("pileup");
    h1_PURunA->Scale(2.1/h1_PURunA->Integral());
    h1_PURunB->Scale(2.5/h1_PURunB->Integral());
    TH1F* h1_PU_weightedAverage = new TH1F(*h1_PURunA);
    h1_PU_weightedAverage->Add(h1_PURunB);
    fPUWeight->SetDataHistogram(h1_PU_weightedAverage);
    fPUWeight_ave->SetDataHistogram(h1_PU_weightedAverage);
  }
    
     


  if( PUType_=="HR11_73pb_DY" ) {
    TFile* filePUMC = TFile::Open("generatedpileup_Zjets_MADGRAPH_AOD423.root");
    TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("GenLevelInfoModule/npileup");
    std::cout << "-> Switching MC PU file to: generatedpileup_Zjets_MADGRAPH_AOD423.root" << std::endl;
    fPUWeight->SetMCHistogram(h1_nPU_mc);
  } else if( dataset_tstr.Contains("Summer11") && dataset_tstr.Contains("PU_S4") && PUType_!="HR11_v3" ) {
    TFile* filePUMC = TFile::Open("Pileup_MC_Summer11_S4.root");
    TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("hNPU");
    std::cout << "-> Switching MC PU file to: Pileup_MC_Summer11_S4.root" << std::endl;
    fPUWeight->SetMCHistogram(h1_nPU_mc);
  } else if( dataset_tstr.Contains("Fall11") ) {
    TFile* filePUMC = TFile::Open("s6MCPileUp.root");
    TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("pileup");
    std::cout << "-> Switching MC PU file to: s6MCPileUp.root" << std::endl;
    fPUWeight->SetMCHistogram(h1_nPU_mc);
  }

  int maxBTag_found = -1;
  float mZZ, mZZ_nokinfit, mZjj, mZll;
  float ptLept1_t, ptLept2_t, etaLept1_t, etaLept2_t;
  float ptJet1_t, ptJet2_t, etaJet1_t, etaJet2_t;
  float QGLikelihoodJet1_t, QGLikelihoodJet2_t;
  float HLTSF;
  bool isSidebands=false;
  bool isSignalRegion=false;
  bool foundSignalRegionMjj=false;
  bool atLeastOneInSignalRegionMjj=false;

  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "LS", &LS, "LS/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "leptType", &leptType, "leptType/I" );
  tree_passedEvents->Branch( "helicityLD", &helicityLD_selected, "helicityLD_selected/F" );
  tree_passedEvents->Branch( "ptLept1", &ptLept1_t, "ptLept1_t/F" );
  tree_passedEvents->Branch( "ptLept2", &ptLept2_t, "ptLept2_t/F" );
  tree_passedEvents->Branch( "etaLept1", &etaLept1_t, "etaLept1_t/F" );
  tree_passedEvents->Branch( "etaLept2", &etaLept2_t, "etaLept2_t/F" );
  tree_passedEvents->Branch( "ptJet1", &ptJet1_t, "ptJet1_t/F" );
  tree_passedEvents->Branch( "ptJet2", &ptJet2_t, "ptJet2_t/F" );
  tree_passedEvents->Branch( "QGLikelihoodJet1", &QGLikelihoodJet1_t, "QGLikelihoodJet1_t/F" );
  tree_passedEvents->Branch( "QGLikelihoodJet2", &QGLikelihoodJet2_t, "QGLikelihoodJet2_t/F" );
  tree_passedEvents->Branch( "etaJet1", &etaJet1_t, "etaJet1_t/F" );
  tree_passedEvents->Branch( "etaJet2", &etaJet2_t, "etaJet2_t/F" );
  tree_passedEvents->Branch( "mZjj", &mZjj, "mZjj/F" );
  tree_passedEvents->Branch( "mZll", &mZll, "mZll/F" );
  tree_passedEvents->Branch( "CMS_hzz2l2q_mZZ", &mZZ, "mZZ/F" );
  tree_passedEvents->Branch( "mZZ_nokinfit", &mZZ_nokinfit, "mZZ_nokinfit/F" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "HLTSF", &HLTSF, "HLTSF/F" );
  tree_passedEvents->Branch( "PUWeight", &eventWeightPU, "eventWeightPU/F" );
  tree_passedEvents->Branch( "nBTags", &maxBTag_found, "maxBTag_found/I" );
  tree_passedEvents->Branch( "isSidebands", &isSidebands, "isSidebands/O" );
  tree_passedEvents->Branch( "isSignalRegion", &isSignalRegion, "isSignalRegion/O" );




ofstream ofs("run_event.txt");


  
  float nEvents_presel=0.;
  float nEvents_presel_mZll=0.;
  float nEvents_presel_mZll_mZjj=0.;
  float nEvents_presel_mZll_mZjj_0btag=0.;
  float nEvents_presel_mZll_mZjj_1btag=0.;
  float nEvents_presel_mZll_mZjj_2btag=0.;


  std::cout << std::endl << std::endl;
  std::cout << "+++ BEGINNING ANALYSIS LOOP" << std::endl;
  std::cout << "----> DATASET: " << dataset_ << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << std::endl;
  if( isMC ) std::cout << "----> PU REWEIGHING: " << PUType_ << std::endl;
  std::cout << std::endl << std::endl;



  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 20000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    eventWeight *= genWeight; //=1 for every dataset except sherpa (for now)


    if( leptType_!="ALL" ) {
      if( leptType_=="ELE" && leptType==0 ) continue;
      if( leptType_=="MU" && leptType==1 ) continue;
    }


  //// BUG FIX in Z->ll BR in Spring11 Alpgen Z+jets:
  //if( dataset_=="ZJets_alpgen_TuneZ2_Spring11_v2" ) {
  //  if( leptType==0 )      eventWeight = eventWeight_Zmm;
  //  else if( leptType==1 ) eventWeight = eventWeight_Zee;
  //}




    h1_nvertex->Fill(nvertex, eventWeight);

    if( isMC ) {

      // scale factor for double mu triggers:
      if( leptType==0 ) {

        float effDouble1_Run2011A = getMuonHLTSF_DoubleTrigger( ptLept1, etaLept1, "Run2011A" );
        float effDouble2_Run2011A = getMuonHLTSF_DoubleTrigger( ptLept2, etaLept2, "Run2011A" );

        float effDouble1_Run2011B = getMuonHLTSF_DoubleTrigger( ptLept1, etaLept1, "Run2011B" );
        float effDouble2_Run2011B = getMuonHLTSF_DoubleTrigger( ptLept2, etaLept2, "Run2011B" );

        float effSingle1_Run2011A1 = getMuonHLTSF_SingleTrigger( ptLept1, etaLept1, "Run2011A1");
        float effSingle2_Run2011A1 = getMuonHLTSF_SingleTrigger( ptLept2, etaLept2, "Run2011A1");

        float effSingle1_Run2011A2 = getMuonHLTSF_SingleTrigger( ptLept1, etaLept1, "Run2011A2");
        float effSingle2_Run2011A2 = getMuonHLTSF_SingleTrigger( ptLept2, etaLept2, "Run2011A2");

        float effSingle1_Run2011A3 = getMuonHLTSF_SingleTrigger( ptLept1, etaLept1, "Run2011A3");
        float effSingle2_Run2011A3 = getMuonHLTSF_SingleTrigger( ptLept2, etaLept2, "Run2011A3");

        float effSingle1_Run2011B = getMuonHLTSF_SingleTrigger( ptLept1, etaLept1, "Run2011B");
        float effSingle2_Run2011B = getMuonHLTSF_SingleTrigger( ptLept2, etaLept2, "Run2011B");


        float HLTSF_Run2011A1 = getEventHLTSF( effSingle1_Run2011A1, effSingle2_Run2011A1, effDouble1_Run2011A, effDouble2_Run2011A );
        float HLTSF_Run2011A2 = getEventHLTSF( effSingle1_Run2011A2, effSingle2_Run2011A2, effDouble1_Run2011A, effDouble2_Run2011A );
        float HLTSF_Run2011A3 = getEventHLTSF( effSingle1_Run2011A3, effSingle2_Run2011A3, effDouble1_Run2011A, effDouble2_Run2011A );
        float HLTSF_Run2011B  = getEventHLTSF( effSingle1_Run2011B, effSingle2_Run2011B, effDouble1_Run2011B, effDouble2_Run2011B );


        // weighted average over full run (weighted with lumi):
        // LP11:
        //HLTSF = (217.*HLTSF_Run2011A1 + 920.*HLTSF_Run2011A2 + 478.*HLTSF_Run2011A3)/(217.+920.+478.);
        if( PUType_=="Run2011A" || PUType_=="Run2011A_73pb" )
          HLTSF = (217.*HLTSF_Run2011A1 + 920.*HLTSF_Run2011A2 + 1000.*HLTSF_Run2011A3)/(217.+920.+1000.);
        else if( PUType_=="HR11" ) 
          HLTSF = (217.*HLTSF_Run2011A1 + 920.*HLTSF_Run2011A2 + 1000.*HLTSF_Run2011A3 + 2100.*HLTSF_Run2011B)/(217.+920.+1000.+2100.);
        else if( PUType_=="HR11_v2" || PUType_=="HR11_73pb" )
          HLTSF = (217.*HLTSF_Run2011A1 + 920.*HLTSF_Run2011A2 + 1000.*HLTSF_Run2011A3 + 2500.*HLTSF_Run2011B)/(217.+920.+1000.+2500.);

        eventWeight *= HLTSF;

      } else { //electrons

        HLTSF = 1.;

      }

      // PU reweighting:
      h1_nvertex_PUW_ave->Fill(nvertex, eventWeight*eventWeightPU_ave);
      //eventWeight *= eventWeightPU;
     
//for( unsigned ii=0; ii<36; ++ii )
//  std::cout << ii << " " << fPUWeight->GetWeight(ii) << std::endl;
//exit(1);
      eventWeight *= fPUWeight->GetWeight(nPU);
      //std::cout << event << " " << nPU << " " << nvertex << " " << fPUWeight->GetWeight(nPU) << std::endl;

    } // if is MC



    h1_nvertex_PUW->Fill(nvertex, eventWeight);
    if( !isMC ) 
      h1_nvertex_PUW_ave->Fill(nvertex, eventWeight);


    if( !isMC ) { 

      // remove duplicate events:

      std::map<int, std::map<int, std::vector<int> > >::iterator it;

      it = run_lumi_ev_map.find(run);


      if( it==run_lumi_ev_map.end() ) {

        std::vector<int> events;
        events.push_back(event);
        std::map<int, std::vector<int> > lumi_ev_map;
        lumi_ev_map.insert( std::pair<int,std::vector<int> >(LS, events));
        run_lumi_ev_map.insert( std::pair<int, std::map<int, std::vector<int> > > (run, lumi_ev_map) );

      } else { //run exists, look for LS


        std::map<int, std::vector<int> >::iterator it_LS;
        it_LS = it->second.find( LS );

        if( it_LS==(it->second.end())  ) {

          std::vector<int> events;
          events.push_back(event);
          it->second.insert( std::pair<int, std::vector<int> > (LS, events) );

        } else { //LS exists, look for event

          std::vector<int>::iterator ev;
          for( ev=it_LS->second.begin(); ev!=it_LS->second.end(); ++ev )
            if( *ev==event ) break;


          if( ev==it_LS->second.end() ) {

            it_LS->second.push_back(event);

          } else {

            std::cout << "DISCARDING DUPLICATE EVENT!! Run: " << run << " LS: " << LS << " event: " << event << std::endl;

            continue;

          }
        }
      }


      h1_run->Fill( run, eventWeight );


      // HLT requirement:
  
      if( dataset_tstr.BeginsWith("SingleMu") ) {

        bool passedHLT = passed_HLT_IsoMu24
                      && !passed_HLT_DoubleMu7 && !passed_HLT_Mu13_Mu8
                      && !passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL && !passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;

        if( !passedHLT ) continue;

      }

    //if( dataset_tstr.BeginsWith("DoubleMu") ) {

    //  bool passedHLT = !passed_HLT_IsoMu24
    //                && passed_HLT_DoubleMu7 && passed_HLT_Mu13_Mu8
    //                && !passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL && !passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;

    //  if( !passedHLT ) continue;

    //}

    
    } //if is not mc


//std::cout << std::endl << "new event" << std::endl;

    h1_rhoPF_presel->Fill( rhoPF, eventWeight);


    TLorentzVector lept1, lept2;
    lept1.SetPtEtaPhiE( ptLept1, etaLept1, phiLept1, eLept1 );
    lept2.SetPtEtaPhiE( ptLept2, etaLept2, phiLept2, eLept2 );

    TLorentzVector diLepton = lept1+lept2;



    nEvents_presel += eventWeight;

    std::vector< std::pair< AnalysisJet, AnalysisJet > > jetPairs_selected;

    if( nPairs>0 ) {

      h1_mZll_presel->Fill( diLepton.M(), eventWeight );
      if( leptType==0 )
        h1_mZmumu_presel->Fill( diLepton.M(), eventWeight );
      else
        h1_mZee_presel->Fill( diLepton.M(), eventWeight );
      h1_ptZll_presel->Fill( diLepton.Pt(), eventWeight );

      h1_deltaRll_presel->Fill(lept1.DeltaR(lept2), eventWeight );

      h1_ptLept1_presel->Fill( ptLept1, eventWeight );
      h1_ptLept2_presel->Fill( ptLept2, eventWeight );

      h1_etaLept1_presel->Fill( etaLept1, eventWeight );
      h1_etaLept2_presel->Fill( etaLept2, eventWeight );

    } else {

      h1_mZll_presel_0jets->Fill( diLepton.M(), eventWeight );
      if( leptType==0 )
        h1_mZmumu_presel_0jets->Fill( diLepton.M(), eventWeight );
      else
        h1_mZee_presel_0jets->Fill( diLepton.M(), eventWeight );

    }




    if( event==DEBUG_EVENTNUMBER ) {
      std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
      std::cout << "** LOG FOR RUN: " << run << "   EVENT: " << DEBUG_EVENTNUMBER << std::endl << std::endl;
      std::cout << "leptType: " << leptType << std::endl; 
      std::cout << "lept1.Pt(): " << lept1.Pt() << " lept1.Eta(): " << lept1.Eta() << std::endl;
      std::cout << "lept2.Pt(): " << lept2.Pt() << " lept2.Eta(): " << lept2.Eta() << std::endl;
      std::cout << "diLepton.M(): " << diLepton.M() << std::endl;
    }


    // ----------------------------
    // KINEMATIC SELECTION: LEPTONS
    // ----------------------------

    if( lept1.Pt() < ptLept1_thresh_ ) continue;
    if( lept2.Pt() < ptLept2_thresh_ ) continue;
    if( fabs(lept1.Eta()) > etaLept1_thresh_ ) continue;
    if( fabs(lept2.Eta()) > etaLept2_thresh_ ) continue;
    if( invert_mZll_ ) {
      if( diLepton.M() > mZll_threshLo_ && diLepton.M() < mZll_threshHi_ ) continue;
      if( diLepton.M() < 50. ) continue; //quick fix for now: madgraph DY does not go under 50 gev
    } else {
      if( diLepton.M() < mZll_threshLo_ || diLepton.M() > mZll_threshHi_ ) continue;
    }

    nEvents_presel_mZll += eventWeight;


    int nJets = getNJets(nPairs);
    h1_nJets_presel->Fill( nJets , eventWeight );
    h1_nPairs_presel->Fill( nPairs , eventWeight );


    float cached_jetpt = 0.;
    AnalysisJet jet1_selected, jet2_selected, jetRecoil_selected;
    float bestMass = 0.;
    int  foundJets = 0;
    int  foundJets_signalRegion = 0;
    bool foundJets_ZZmass = false;
    maxBTag_found = -1;
    foundSignalRegionMjj = false;
    atLeastOneInSignalRegionMjj = false;

    for( unsigned iJetPair=0; iJetPair<nPairs; ++iJetPair ) {

      bool isSignalRegionMjj = false;

      AnalysisJet jet1, jet2;
      jet1.SetPtEtaPhiE( ptJet1[iJetPair], etaJet1[iJetPair], phiJet1[iJetPair], eJet1[iJetPair]);
      jet2.SetPtEtaPhiE( ptJet2[iJetPair], etaJet2[iJetPair], phiJet2[iJetPair], eJet2[iJetPair]);

      jet1.rmsCand = rmsCandJet1[iJetPair];
      jet1.ptD = ptDJet1[iJetPair];
      jet1.nCharged = nChargedJet1[iJetPair];
      jet1.nNeutral = nNeutralJet1[iJetPair];
      jet1.muonEnergyFraction = eMuonsJet1[iJetPair]/jet1.Energy();
      jet1.electronEnergyFraction = eElectronsJet1[iJetPair]/jet1.Energy();

      jet1.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet1[iJetPair];
      jet1.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet1[iJetPair];
      jet1.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet1[iJetPair];
      jet1.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet1[iJetPair];
      jet1.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet1[iJetPair];
      jet1.jetProbabilityBJetTag               = jetProbabilityBJetTagJet1[iJetPair];

      jet1.ptGen = ptJet1Gen[iJetPair];
      jet1.etaGen = etaJet1Gen[iJetPair];
      jet1.phiGen = phiJet1Gen[iJetPair];
      jet1.eGen = eJet1Gen[iJetPair];


      jet2.rmsCand = rmsCandJet2[iJetPair];
      jet2.ptD = ptDJet2[iJetPair];
      jet2.nCharged = nChargedJet2[iJetPair];
      jet2.nNeutral = nNeutralJet2[iJetPair];
      jet2.muonEnergyFraction = eMuonsJet2[iJetPair]/jet2.Energy();
      jet2.electronEnergyFraction = eElectronsJet2[iJetPair]/jet2.Energy();

      jet2.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet2[iJetPair];
      jet2.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet2[iJetPair];
      jet2.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet2[iJetPair];
      jet2.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet2[iJetPair];
      jet2.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet2[iJetPair];
      jet2.jetProbabilityBJetTag               = jetProbabilityBJetTagJet2[iJetPair];

      jet2.ptGen =   ptJet2Gen[iJetPair];
      jet2.etaGen = etaJet2Gen[iJetPair];
      jet2.phiGen = phiJet2Gen[iJetPair];
      jet2.eGen =     eJet2Gen[iJetPair];

      AnalysisJet jetRecoil;
      jetRecoil.SetPtEtaPhiE( ptJetRecoil[iJetPair], etaJetRecoil[iJetPair], phiJetRecoil[iJetPair], eJetRecoil[iJetPair]);
      jetRecoil.ptD = ptDJetRecoil[iJetPair];
      jetRecoil.nCharged = nChargedJetRecoil[iJetPair];
      jetRecoil.nNeutral = nNeutralJetRecoil[iJetPair];


      TLorentzVector diJet = jet1 + jet2;

      if( event==DEBUG_EVENTNUMBER ) {
        std::cout << std::endl << "jet pair N.: " << iJetPair << std::endl; 
        std::cout << "jet1.Pt(): " << jet1.Pt() << " jet1.Eta(): " << jet1.Eta() << std::endl;
        std::cout << "jet2.Pt(): " << jet2.Pt() << " jet2.Eta(): " << jet2.Eta() << std::endl;
        std::cout << "diJet.M(): " << diJet.M() << std::endl;
      }


      // fill histos before selection
      h1_mZjj_all_presel->Fill( diJet.M(), eventWeight );
      h1_ptZjj_all_presel->Fill( diJet.Pt(), eventWeight );
      h1_deltaRjj_all_presel->Fill(jet1.DeltaR(jet2), eventWeight );

      if( jet1.Pt()!=cached_jetpt ) {
        h1_ptJet_all_presel->Fill( jet1.Pt(), eventWeight );
        h1_etaJet_all_presel->Fill( jet1.Eta(), eventWeight );
        cached_jetpt = jet1.Pt();
      }



      //match to parton:
      int partFlavor1=0;
      float deltaRmin1=999.;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = jet1.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin1 ) {
          partFlavor1 = pdgIdPart[iPart];
          deltaRmin1 = thisDeltaR;
        }
      }
      jet1.pdgIdPart = partFlavor1;

      int partFlavor2=0;
      float deltaRmin2=999.;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = jet2.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin2 ) {
          partFlavor2 = pdgIdPart[iPart];
          deltaRmin2 = thisDeltaR;
        }
      }
      jet2.pdgIdPart = partFlavor2;



      // -------------------------
      // KINEMATIC SELECTION: JETS
      // -------------------------

      if( jet1.Pt() < ptJet1_thresh_ ) continue;
      if( event==DEBUG_EVENTNUMBER ) std::cout << "first jet pt OK" << std::endl;
      if( jet2.Pt() < ptJet2_thresh_ ) continue;
      if( event==DEBUG_EVENTNUMBER ) std::cout << "second jet pt OK" << std::endl;
      if( fabs(jet1.Eta()) > etaJet1_thresh_ ) continue;
      if( event==DEBUG_EVENTNUMBER ) std::cout << "first jet eta OK" << std::endl;
      if( fabs(jet2.Eta()) > etaJet2_thresh_ ) continue;
      if( event==DEBUG_EVENTNUMBER ) std::cout << "second jet eta OK" << std::endl;


      // sideband logic:
      if( foundSignalRegionMjj && ( diJet.M() < mZjj_threshLo_ || diJet.M() > mZjj_threshHi_ ) ) continue;


      if( diJet.M() >= mZjj_threshLo_ && diJet.M() <= mZjj_threshHi_ ) {
        if( !atLeastOneInSignalRegionMjj ) nEvents_presel_mZll_mZjj += eventWeight;
        atLeastOneInSignalRegionMjj = true;
        isSignalRegionMjj = true;
      }




      // ----------
      // B-TAGGING:
      // ----------


      int nBTags = this->get_nBTags( jet1, jet2, btsfutil, use_looseBTags_ );
      if( event==DEBUG_EVENTNUMBER ) std::cout << "nBTags: " << nBTags << std::endl;

      if( (isSignalRegionMjj && foundSignalRegionMjj) ||
          (!isSignalRegionMjj && !foundSignalRegionMjj) ) {
        if( nBTags<maxBTag_found ) continue; //speed it up a little
      }






      // -------------------
      // Q-G DISCRIMINATION:
      // -------------------

      jet1.QGLikelihood = qglikeli->computeQGLikelihoodPU( jet1.Pt(), rhoPF, jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );
      jet2.QGLikelihood = qglikeli->computeQGLikelihoodPU( jet2.Pt(), rhoPF, jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );
      jet1.QGLikelihoodNoPU = qglikeli->computeQGLikelihood( jet1.Pt(), jet1.nCharged, jet1.nNeutral, jet1.ptD, -1. );
      jet2.QGLikelihoodNoPU = qglikeli->computeQGLikelihood( jet2.Pt(), jet2.nCharged, jet2.nNeutral, jet2.ptD, -1. );
      float QGLikelihoodProd = jet1.QGLikelihood*jet2.QGLikelihood;
      if( event==DEBUG_EVENTNUMBER ) std::cout << "rhoPF: " << rhoPF << std::endl;
      if( event==DEBUG_EVENTNUMBER ) std::cout << "jet1: nCharged: " << jet1.nCharged << " nNeutral: " << jet1.nNeutral << " ptD: " << jet1.ptD << " QGlikelihood: " << jet1.QGLikelihood << std::endl;
      if( event==DEBUG_EVENTNUMBER ) std::cout << "jet2: nCharged: " << jet2.nCharged << " nNeutral: " << jet2.nNeutral << " ptD: " << jet2.ptD << " QGlikelihood: " << jet2.QGLikelihood << std::endl;
      if( event==DEBUG_EVENTNUMBER ) std::cout << "QGLikelihoodProd: " << QGLikelihoodProd << std::endl;
      if( nBTags==0 ) {
        //if( QGLikelihoodProd < QGLikelihoodProd_thresh_ ) continue;
        if( QGLikelihoodProd < QGLikelihoodProd_thresh_ ) nBTags=-1; //glue-tag category
      }





      // --------------
      // KINEMATIC FIT:
      // --------------

  
      if( USE_MC_MASS ) {
        TLorentzVector ZjjMC;
        TLorentzVector ZllMC;
        int iZllMC=-1;
        // first look for ZllMC
        float bestDeltaRToZll=9999.;
        for( unsigned iPart=0; iPart<nPart; ++iPart ) {
          if( pdgIdPart[iPart]!=23 ) continue;
          TLorentzVector thisZMC;
          thisZMC.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
          if( thisZMC.DeltaR(diLepton) < bestDeltaRToZll ) {
            ZllMC=thisZMC;
            iZllMC = iPart;
            bestDeltaRToZll=thisZMC.DeltaR(diLepton);
          }
        }

        if( iZllMC>=0 ) {

          bool foundZjjMC=false;
          for( unsigned iPart=0; iPart<nPart && !foundZjjMC; ++iPart ) {
            if( pdgIdPart[iPart]==23 && iPart!=iZllMC ) {
              ZjjMC.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
              foundZjjMC=true;
            }
          }

          if( foundZjjMC ) {
            h1_mZjjMC->Fill( ZjjMC.M(), eventWeight );
            fitter_jets->set_mass(ZjjMC.M());
          }

        } // if found ZllMC

      } //if use mc z mass

      std::pair<TLorentzVector,TLorentzVector> jets_kinfit = fitter_jets->fit(jet1, jet2);
      TLorentzVector jet1_kinfit(jets_kinfit.first);
      TLorentzVector jet2_kinfit(jets_kinfit.second);

      jet1.pt_preKinFit  = jet1.Pt();
      jet1.eta_preKinFit = jet1.Eta();
      jet1.phi_preKinFit = jet1.Phi();
      jet1.e_preKinFit   = jet1.Energy();

      jet2.pt_preKinFit  = jet2.Pt();
      jet2.eta_preKinFit = jet2.Eta();
      jet2.phi_preKinFit = jet2.Phi();
      jet2.e_preKinFit   = jet2.Energy();

      // update 4-vector to kinfit results:
      jet1.SetPtEtaPhiE( jet1_kinfit.Pt(), jet1_kinfit.Eta(), jet1_kinfit.Phi(), jet1_kinfit.Energy() );
      jet2.SetPtEtaPhiE( jet2_kinfit.Pt(), jet2_kinfit.Eta(), jet2_kinfit.Phi(), jet2_kinfit.Energy() );

      TLorentzVector jet1_nokinfit, jet2_nokinfit;
      jet1_nokinfit.SetPtEtaPhiE( jet1.pt_preKinFit, jet1.eta_preKinFit, jet1.phi_preKinFit, jet1.e_preKinFit );
      jet2_nokinfit.SetPtEtaPhiE( jet2.pt_preKinFit, jet2.eta_preKinFit, jet2.phi_preKinFit, jet2.e_preKinFit );

      TLorentzVector diJet_kinfit = jet1_kinfit + jet2_kinfit;
      TLorentzVector ZZ_kinfit_tmp = diJet_kinfit + diLepton;

      if( event==DEBUG_EVENTNUMBER ) std::cout << "ZZ_kinfit_tmp.M(): " << ZZ_kinfit_tmp.M() << std::endl;
      if( ZZ_kinfit_tmp.M()<150. ) continue; //speed it up a little


      // --------------------
      // FULL EVENT VARIABLES
      // --------------------
   
      if( event==DEBUG_EVENTNUMBER ) std::cout << "metSignificance: " << metSignificance << std::endl;
      if( nBTags==2 )  {
        if( metSignificance > metSignificance_thresh_ ) continue;
      }



      // ------------
      // HELICITY LD:
      // ------------

      HelicityLikelihoodDiscriminant::HelicityAngles hangles;
      if( chargeLept1<0 ) hangles = LD->computeHelicityAngles(lept1, lept2, jet1, jet2);
      else                hangles = LD->computeHelicityAngles(lept2, lept1, jet1, jet2);
    
      LD->setMeasurables(hangles);
      double sProb = LD->getSignalProbability();
      double bProb = LD->getBkgdProbability();
      double helicityLD = sProb/(sProb+bProb);
    
      HelicityLikelihoodDiscriminant::HelicityAngles hangles_nokinfit;
      if( chargeLept1<0 ) hangles_nokinfit = LD->computeHelicityAngles(lept1, lept2, jet1_nokinfit, jet2_nokinfit);
      else                hangles_nokinfit = LD->computeHelicityAngles(lept2, lept1, jet1_nokinfit, jet2_nokinfit);
    
      LD->setMeasurables(hangles_nokinfit);
      double sProb_nokinfit=LD->getSignalProbability();
      double bProb_nokinfit=LD->getBkgdProbability();
      double helicityLD_nokinfit=sProb_nokinfit/(sProb_nokinfit+bProb_nokinfit);
    
      float helicityLD_thresh = (nBTags>=0) ? this->get_helicityLD_thresh(ZZ_kinfit_tmp.M(), nBTags) : this->get_helicityLD_thresh(ZZ_kinfit_tmp.M(), 0);

      if( event==DEBUG_EVENTNUMBER ) std::cout << "helicityLD: " << helicityLD << " (threshold: " << helicityLD_thresh << ")" << std::endl;
      if( helicityLD < helicityLD_thresh ) continue;



      // this pair has passed selections
      // if it is a signal region pair, dont consider any more sideband events:
      if( isSignalRegionMjj && !foundSignalRegionMjj) {
        if( event==DEBUG_EVENTNUMBER ) std::cout << "updating maxBTag_found: " << nBTags << std::endl;
        maxBTag_found = nBTags; //reset maxbtag found
        foundSignalRegionMjj = true;
      }
      if( foundSignalRegionMjj && !isSignalRegionMjj ) continue;



      // logic: choose jet pair with highest number of btags which passes selection
      // if more than one pair found, choose the one with the invariant mass closest to Z mass.
      // keep also best-Zmass dijet pair which passes all cuts except mZjj for sideband study

      float invMass = diJet.M();

      if( foundJets==0 || (foundJets_signalRegion==0 && isSignalRegionMjj) ) {

        if( event==DEBUG_EVENTNUMBER ) {
          if(foundJets_signalRegion==0 && isSignalRegionMjj)   std::cout << "first **SIGNAL REGION** jet pair passing selections" << std::endl;
          else std::cout << "first jet pair passing selections" << std::endl;
        }
        bestMass = invMass;
        jet1_selected = jet1;
        jet2_selected = jet2; 
        jetRecoil_selected = jetRecoil; 
        hangles_selected = hangles;
        helicityLD_selected = helicityLD;
        helicityLD_nokinfit_selected = helicityLD_nokinfit;
        foundJets += 1;
        if( isSignalRegionMjj ) foundJets_signalRegion += 1;
        maxBTag_found = nBTags;

      } else { 

        foundJets += 1;

        if( event==DEBUG_EVENTNUMBER ) std::cout << "jet pair passed selections, going to fight for first place:" << std::endl;
        if( event==DEBUG_EVENTNUMBER ) std::cout << "nBTags: " << nBTags << " > " << maxBTag_found << " ?" << std::endl;
        if( event==DEBUG_EVENTNUMBER ) std::cout << "invMass: " << invMass << " better than " << bestMass << " ?" << std::endl;

        if( ( nBTags == maxBTag_found && (fabs(invMass-Zmass) < fabs(bestMass-Zmass)))
         || (nBTags > maxBTag_found)  ) {
          if( event==DEBUG_EVENTNUMBER ) std::cout << "---> YES! this pair is now in pole position." << std::endl;
          bestMass = invMass;
          jet1_selected = jet1;
          jet2_selected = jet2;
          jetRecoil_selected = jetRecoil; 
          hangles_selected = hangles;
          helicityLD_selected = helicityLD;
          helicityLD_nokinfit_selected = helicityLD_nokinfit;
          maxBTag_found = nBTags;
        } else {
          if( event==DEBUG_EVENTNUMBER ) std::cout << "NOPE. keeping old pair (" << bestMass << ")" << std::endl;
        }

      }

      h1_nEventsCategories_presel->Fill( nBTags, eventWeight );

    } //for on jet pairs


    if( event==DEBUG_EVENTNUMBER ) std::cout << "----------------------------------" << std::endl;
    if( event==DEBUG_EVENTNUMBER ) std::cout << "foundJets:" << foundJets << std::endl;
    if( event==DEBUG_EVENTNUMBER ) std::cout << "maxBTag_found:" << maxBTag_found << std::endl;
    if( event==DEBUG_EVENTNUMBER ) {
      if( foundSignalRegionMjj )  std::cout << "foundSignalRegionMjj is TRUE" << std::endl;
      else std::cout << "foundSignalRegionMjj is FALSE" << std::endl;
    }



    h1_nCandidates->Fill( foundJets, eventWeight );


    if( foundJets==0 ) continue;


    // event has passed selection (except mZjj cut). define all objects:

    TLorentzVector Zjj_kinfit = jet1_selected + jet2_selected;

    TLorentzVector jet1_nokinfit, jet2_nokinfit;
    jet1_nokinfit.SetPtEtaPhiE( jet1_selected.pt_preKinFit, jet1_selected.eta_preKinFit, jet1_selected.phi_preKinFit, jet1_selected.e_preKinFit );
    jet2_nokinfit.SetPtEtaPhiE( jet2_selected.pt_preKinFit, jet2_selected.eta_preKinFit, jet2_selected.phi_preKinFit, jet2_selected.e_preKinFit );

    TLorentzVector Zjj_nokinfit = jet1_nokinfit + jet2_nokinfit;

    TLorentzVector ZZ_nokinfit = Zjj_nokinfit + diLepton;
    TLorentzVector ZZ_kinfit = diLepton + Zjj_kinfit;

    ptLept1_t = ptLept1;
    ptLept2_t = ptLept2;
    etaLept1_t = etaLept1;
    etaLept2_t = etaLept2;
    ptJet1_t = jet1_selected.Pt();
    ptJet2_t = jet2_selected.Pt();
    QGLikelihoodJet1_t = jet1_selected.QGLikelihood;
    QGLikelihoodJet2_t = jet2_selected.QGLikelihood;
    etaJet1_t = jet1_selected.Pt();
    etaJet2_t = jet2_selected.Eta();

    mZjj = Zjj_nokinfit.M();
    mZll = diLepton.M();
    mZZ = ZZ_kinfit.M();
    mZZ_nokinfit = ZZ_nokinfit.M();

    isSignalRegion = (Zjj_nokinfit.M()>=75. && Zjj_nokinfit.M()<=105.);
    isSidebands = !isSignalRegion &&  Zjj_nokinfit.M()>60. && Zjj_nokinfit.M()<130.;

    // and fill tree:
    tree_passedEvents->Fill();



    // fill mZjj plots for all events:

    h1_mZjj->Fill( Zjj_nokinfit.M(), eventWeight);
 
    if( leptType==0 ) {
      h1_mZjj_MU->Fill( Zjj_nokinfit.M(), eventWeight);
      if( maxBTag_found>=0 ) h1_mZjj_nogluetag_MU->Fill( Zjj_nokinfit.M(), eventWeight);
    }
    if( leptType==1 ) {
      h1_mZjj_ELE->Fill( Zjj_nokinfit.M(), eventWeight);
      if( maxBTag_found>=0 ) h1_mZjj_nogluetag_ELE->Fill( Zjj_nokinfit.M(), eventWeight);
    }
    if( mZZ>150. && mZZ<250. ) h1_mZjj_loMass->Fill( Zjj_nokinfit.M(), eventWeight);
    else if( mZZ>250. && mZZ<400. ) h1_mZjj_medMass->Fill( Zjj_nokinfit.M(), eventWeight);
    else if( mZZ>400. ) h1_mZjj_hiMass->Fill( Zjj_nokinfit.M(), eventWeight);
    if( maxBTag_found==0 ) h1_mZjj_0btag->Fill( Zjj_nokinfit.M(), eventWeight);
    else if( maxBTag_found==1 ) h1_mZjj_1btag->Fill( Zjj_nokinfit.M(), eventWeight);
    else if( maxBTag_found==2 ) h1_mZjj_2btag->Fill( Zjj_nokinfit.M(), eventWeight);

    h2_mZjj_vs_mZZ->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    if( maxBTag_found==0 ) h2_mZjj_vs_mZZ_0btag->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==1 ) h2_mZjj_vs_mZZ_1btag->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==2 ) h2_mZjj_vs_mZZ_2btag->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==-1 ) h2_mZjj_vs_mZZ_gluetag->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );

    h2_mZjj_vs_mZZ_kinfit->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M() );
    if( maxBTag_found==0 )       h2_mZjj_vs_mZZ_kinfit_0btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==1 )  h2_mZjj_vs_mZZ_kinfit_1btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==2 )  h2_mZjj_vs_mZZ_kinfit_2btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==-1 ) h2_mZjj_vs_mZZ_kinfit_gluetag->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M() );


    //if( Zjj_nokinfit.M()<60. || Zjj_nokinfit.M()>130. ) continue;


    
    if( isSidebands ) {
 
      //fill sideband plots:

      h2_mZjj_vs_mZZ->Fill( mZZ, Zjj_nokinfit.M(), eventWeight );
      h2_mZjj_vs_mZZ_kinfit->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M(), eventWeight );

      h1_helicityLD_sidebands->Fill( helicityLD_selected, eventWeight );

      if( maxBTag_found==0 ) {
        h2_mZjj_vs_mZZ_0btag->Fill( mZZ, Zjj_nokinfit.M() );
        h2_mZjj_vs_mZZ_kinfit_0btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
        h1_mZZ_kinfit_hiMass_sidebands_0btag->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==0 ) h1_mZZ_kinfit_hiMass_sidebands_0btag_MU->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==1 ) h1_mZZ_kinfit_hiMass_sidebands_0btag_ELE->Fill( ZZ_kinfit.M(), eventWeight );
      } else if( maxBTag_found==1 ) {
        h2_mZjj_vs_mZZ_1btag->Fill( mZZ, Zjj_nokinfit.M() );
        h2_mZjj_vs_mZZ_kinfit_1btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
        h1_mZZ_kinfit_hiMass_sidebands_1btag->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==0 ) h1_mZZ_kinfit_hiMass_sidebands_1btag_MU->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==1 ) h1_mZZ_kinfit_hiMass_sidebands_1btag_ELE->Fill( ZZ_kinfit.M(), eventWeight );
      } else if( maxBTag_found==2 ) {
        h2_mZjj_vs_mZZ_2btag->Fill( mZZ, Zjj_nokinfit.M() );
        h2_mZjj_vs_mZZ_kinfit_2btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
        h1_mZZ_kinfit_hiMass_sidebands_2btag->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==0 ) h1_mZZ_kinfit_hiMass_sidebands_2btag_MU->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==1 ) h1_mZZ_kinfit_hiMass_sidebands_2btag_ELE->Fill( ZZ_kinfit.M(), eventWeight );
      } else if( maxBTag_found==-1 ) {
        h2_mZjj_vs_mZZ_gluetag->Fill( mZZ, Zjj_nokinfit.M() );
        h2_mZjj_vs_mZZ_kinfit_gluetag->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M() );
        h1_mZZ_kinfit_hiMass_sidebands_gluetag->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==0 ) h1_mZZ_kinfit_hiMass_sidebands_gluetag_MU->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==1 ) h1_mZZ_kinfit_hiMass_sidebands_gluetag_ELE->Fill( ZZ_kinfit.M(), eventWeight );
      }

      

    } else if( isSignalRegion ) {
    


      if( maxBTag_found==0 ) nEvents_presel_mZll_mZjj_0btag += eventWeight;
      if( maxBTag_found==1 ) nEvents_presel_mZll_mZjj_1btag += eventWeight;
      if( maxBTag_found==2 ) nEvents_presel_mZll_mZjj_2btag += eventWeight;
  
  
  
      if( helicityLD_selected < 0. ) 
        std::cout << "helicityLD_selected is less than 0!!! THIS IS NOT POSSIBLE!!" << std::endl;
  
  
        if( maxBTag_found>=0 ) ofs << run << " " << LS << " " << event << std::endl;
  
  
  
  
  
      //if( Zjj_nokinfit.M() < mZjj_threshLo_ || Zjj_nokinfit.M() > mZjj_threshHi_ ) continue;
  
    
  
      // percentages which define the cut and count windows:
      float mZZ_minPerc = 0.94;  //settled on -6/+10%
      float mZZ_maxPerc = 1.1;
  
      if( ZZ_kinfit.M() > 250.*mZZ_minPerc && ZZ_kinfit.M() < 250.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_250  += eventWeight;
          nEventsPassed_0btag_250++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_250_MU  += eventWeight;
            nEventsPassed_0btag_250_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_250_ELE  += eventWeight;
            nEventsPassed_0btag_250_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_250 += eventWeight;
          nEventsPassed_1btag_250++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_250_MU  += eventWeight;
            nEventsPassed_1btag_250_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_250_ELE  += eventWeight;
            nEventsPassed_1btag_250_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_250 += eventWeight;
          nEventsPassed_2btag_250++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_250_MU  += eventWeight;
            nEventsPassed_2btag_250_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_250_ELE  += eventWeight;
            nEventsPassed_2btag_250_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_250 += eventWeight;
          nEventsPassed_gluetag_250++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_250_MU  += eventWeight;
            nEventsPassed_gluetag_250_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_250_ELE  += eventWeight;
            nEventsPassed_gluetag_250_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 300.*mZZ_minPerc && ZZ_kinfit.M() < 300.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_300  += eventWeight;
          nEventsPassed_0btag_300++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_300_MU  += eventWeight;
            nEventsPassed_0btag_300_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_300_ELE  += eventWeight;
            nEventsPassed_0btag_300_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_300 += eventWeight;
          nEventsPassed_1btag_300++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_300_MU  += eventWeight;
            nEventsPassed_1btag_300_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_300_ELE  += eventWeight;
            nEventsPassed_1btag_300_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_300 += eventWeight;
          nEventsPassed_2btag_300++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_300_MU  += eventWeight;
            nEventsPassed_2btag_300_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_300_ELE  += eventWeight;
            nEventsPassed_2btag_300_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_300 += eventWeight;
          nEventsPassed_gluetag_300++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_300_MU  += eventWeight;
            nEventsPassed_gluetag_300_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_300_ELE  += eventWeight;
            nEventsPassed_gluetag_300_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 350.*mZZ_minPerc && ZZ_kinfit.M() < 350.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_350  += eventWeight;
          nEventsPassed_0btag_350++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_350_MU  += eventWeight;
            nEventsPassed_0btag_350_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_350_ELE  += eventWeight;
            nEventsPassed_0btag_350_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_350 += eventWeight;
          nEventsPassed_1btag_350++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_350_MU  += eventWeight;
            nEventsPassed_1btag_350_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_350_ELE  += eventWeight;
            nEventsPassed_1btag_350_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_350 += eventWeight;
          nEventsPassed_2btag_350++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_350_MU  += eventWeight;
            nEventsPassed_2btag_350_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_350_ELE  += eventWeight;
            nEventsPassed_2btag_350_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_350 += eventWeight;
          nEventsPassed_gluetag_350++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_350_MU  += eventWeight;
            nEventsPassed_gluetag_350_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_350_ELE  += eventWeight;
            nEventsPassed_gluetag_350_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 400.*mZZ_minPerc && ZZ_kinfit.M() < 400.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_400  += eventWeight;
          nEventsPassed_0btag_400++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_400_MU  += eventWeight;
            nEventsPassed_0btag_400_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_400_ELE  += eventWeight;
            nEventsPassed_0btag_400_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_400 += eventWeight;
          nEventsPassed_1btag_400++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_400_MU  += eventWeight;
            nEventsPassed_1btag_400_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_400_ELE  += eventWeight;
            nEventsPassed_1btag_400_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_400 += eventWeight;
          nEventsPassed_2btag_400++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_400_MU  += eventWeight;
            nEventsPassed_2btag_400_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_400_ELE  += eventWeight;
            nEventsPassed_2btag_400_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_400 += eventWeight;
          nEventsPassed_gluetag_400++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_400_MU  += eventWeight;
            nEventsPassed_gluetag_400_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_400_ELE  += eventWeight;
            nEventsPassed_gluetag_400_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 450.*mZZ_minPerc && ZZ_kinfit.M() < 450.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_450  += eventWeight;
          nEventsPassed_0btag_450++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_450_MU  += eventWeight;
            nEventsPassed_0btag_450_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_450_ELE  += eventWeight;
            nEventsPassed_0btag_450_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_450 += eventWeight;
          nEventsPassed_1btag_450++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_450_MU  += eventWeight;
            nEventsPassed_1btag_450_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_450_ELE  += eventWeight;
            nEventsPassed_1btag_450_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_450 += eventWeight;
          nEventsPassed_2btag_450++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_450_MU  += eventWeight;
            nEventsPassed_2btag_450_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_450_ELE  += eventWeight;
            nEventsPassed_2btag_450_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_450 += eventWeight;
          nEventsPassed_gluetag_450++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_450_MU  += eventWeight;
            nEventsPassed_gluetag_450_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_450_ELE  += eventWeight;
            nEventsPassed_gluetag_450_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 500.*mZZ_minPerc && ZZ_kinfit.M() < 500.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_500  += eventWeight;
          nEventsPassed_0btag_500++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_500_MU  += eventWeight;
            nEventsPassed_0btag_500_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_500_ELE  += eventWeight;
            nEventsPassed_0btag_500_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_500 += eventWeight;
          nEventsPassed_1btag_500++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_500_MU  += eventWeight;
            nEventsPassed_1btag_500_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_500_ELE  += eventWeight;
            nEventsPassed_1btag_500_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_500 += eventWeight;
          nEventsPassed_2btag_500++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_500_MU  += eventWeight;
            nEventsPassed_2btag_500_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_500_ELE  += eventWeight;
            nEventsPassed_2btag_500_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_500 += eventWeight;
          nEventsPassed_gluetag_500++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_500_MU  += eventWeight;
            nEventsPassed_gluetag_500_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_500_ELE  += eventWeight;
            nEventsPassed_gluetag_500_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 600.*mZZ_minPerc && ZZ_kinfit.M() < 600.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_600  += eventWeight;
          nEventsPassed_0btag_600++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_600_MU  += eventWeight;
            nEventsPassed_0btag_600_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_600_ELE  += eventWeight;
            nEventsPassed_0btag_600_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_600 += eventWeight;
          nEventsPassed_1btag_600++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_600_MU  += eventWeight;
            nEventsPassed_1btag_600_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_600_ELE  += eventWeight;
            nEventsPassed_1btag_600_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_600 += eventWeight;
          nEventsPassed_2btag_600++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_600_MU  += eventWeight;
            nEventsPassed_2btag_600_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_600_ELE  += eventWeight;
            nEventsPassed_2btag_600_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_600 += eventWeight;
          nEventsPassed_gluetag_600++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_600_MU  += eventWeight;
            nEventsPassed_gluetag_600_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_600_ELE  += eventWeight;
            nEventsPassed_gluetag_600_ELE++;
          }
        }
      } 
  
  
  
      // match to partons:
      TLorentzVector matchedPart1, matchedPart2;
      float bestDeltaRPart1=999.;
      float bestDeltaRPart2=999.;
      for( unsigned iPart=0; iPart<nPart; ++iPart ) {
        if( abs(pdgIdPart[iPart])>6 ) continue;
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        if( jet1_selected.DeltaR(thisPart) < bestDeltaRPart1 ) {
          bestDeltaRPart1 = jet1_selected.DeltaR(thisPart);
          matchedPart1 = thisPart;
        }
        if( jet2_selected.DeltaR(thisPart) < bestDeltaRPart2 ) {
          bestDeltaRPart2 = jet2_selected.DeltaR(thisPart);
          matchedPart2 = thisPart;
        }
      }
      float ptReso1_before = (isMC) ? ( jet1_nokinfit.Pt()-matchedPart1.Pt() )/matchedPart1.Pt() : 0.;
      float ptReso2_before = (isMC) ? ( jet2_nokinfit.Pt()-matchedPart2.Pt() )/matchedPart2.Pt() : 0.;
      h1_ptResoJet1_beforeKin->Fill( ptReso1_before, eventWeight );
      h1_ptResoJet2_beforeKin->Fill( ptReso2_before, eventWeight );
  
      float ptReso1_after = (isMC) ? ( jet1_selected.Pt()-matchedPart1.Pt() )/matchedPart1.Pt() : 0.;
      float ptReso2_after = (isMC) ? ( jet2_selected.Pt()-matchedPart2.Pt() )/matchedPart2.Pt() : 0.;
      h1_ptResoJet1_afterKin->Fill( ptReso1_after, eventWeight );
      h1_ptResoJet2_afterKin->Fill( ptReso2_after, eventWeight );
  
  
  
      TLorentzVector matchedZ;
      float bestDeltaRZ=999.;
      for( unsigned iPart=0; iPart<nPart; ++iPart ) {
        if( pdgIdPart[iPart]!=23 ) continue;
        TLorentzVector thisZ;
        thisZ.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        if( Zjj_kinfit.DeltaR(thisZ) < bestDeltaRZ ) {
          bestDeltaRZ = Zjj_kinfit.DeltaR(thisZ);
          matchedZ = thisZ;
        }
      }
  
      bool eventIsMatched = bestDeltaRZ<0.2;
  
      float ptZreso_before = (isMC) ? (Zjj_nokinfit.Pt()-matchedZ.Pt())/matchedZ.Pt() : 0.;
      float ptZreso_after  = (isMC) ? (Zjj_kinfit.Pt()-matchedZ.Pt())/matchedZ.Pt() : 0.;
      h1_ptZreso_beforeKin->Fill( ptZreso_before, eventWeight);
      h1_ptZreso_afterKin->Fill( ptZreso_after, eventWeight);
  
  
      TLorentzVector matchedH;
      bool foundH=false;
      for( unsigned iPart=0; iPart<nPart && !foundH; ++iPart ) {
        if( pdgIdPart[iPart]!=25 && pdgIdPart[iPart]!=39 ) continue;
        matchedH.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        foundH=true;
      }
  
      if( foundH ) {
        float mHreso_beforeKin = (ZZ_nokinfit.M()-matchedH.M())/matchedH.M();
        float mHreso_afterKin = (ZZ_kinfit.M()-matchedH.M())/matchedH.M();
        h1_mHreso_beforeKin->Fill(mHreso_beforeKin, eventWeight);
        h1_mHreso_afterKin->Fill(mHreso_afterKin, eventWeight);
      }
  
      //compare kinfit to Zmass constraint
      TLorentzVector Zjj_constr;
      Zjj_constr.SetXYZM( Zjj_nokinfit.Px(), Zjj_nokinfit.Py(), Zjj_nokinfit.Pz(), Zmass);
  
      TLorentzVector ZZ_constr = diLepton + Zjj_constr;
  
  
      float chiSquareProb = TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF());
      h1_kinfit_chiSquare->Fill( fitter_jets->getS()/fitter_jets->getNDF(), eventWeight ); 
      h1_kinfit_chiSquareProb->Fill( chiSquareProb, eventWeight ); 
  
  
  
  
  
  
      h1_pfMet->Fill( pfMet, eventWeight );
      h1_pfMetOverMZZ->Fill( pfMet/ZZ_kinfit.M(), eventWeight );
      h1_metSignificance->Fill( metSignificance, eventWeight );
      h1_mEtSig->Fill( mEtSig, eventWeight );
      if( maxBTag_found==2 ) {
        h1_pfMet_2btag->Fill( pfMet, eventWeight );
        h1_pfMetOverMZZ_2btag->Fill( pfMet/ZZ_kinfit.M(), eventWeight );
        h1_metSignificance_2btag->Fill( metSignificance, eventWeight );
        h1_mEtSig_2btag->Fill( mEtSig, eventWeight );
      }
  
      h1_rhoPF->Fill( rhoPF, eventWeight );
  
      h2_helicityLD_vs_mZZ->Fill( ZZ_kinfit.M(), helicityLD_selected, eventWeight );

      h1_tcheJet->Fill( jet1_selected.trackCountingHighEffBJetTag, eventWeight );
      h1_tcheJet->Fill( jet2_selected.trackCountingHighEffBJetTag, eventWeight );

      h1_ptJetRecoil->Fill( jetRecoil_selected.Pt(), eventWeight );
      h1_ptHiggs->Fill( ZZ_kinfit.Pt(), eventWeight );

      if( jetRecoil_selected.Pt()>0. && fabs(jetRecoil_selected.Eta())<2.4 ) {
        float QGLikelihoodJetRecoil = qglikeli->computeQGLikelihoodPU( jetRecoil_selected.Pt(), rhoPF, jetRecoil_selected.nCharged, jetRecoil_selected.nNeutral, jetRecoil_selected.ptD, -1. );
        h1_QGLikelihoodJetRecoil->Fill( QGLikelihoodJetRecoil, eventWeight);
        h1_QGLikelihoodProdRecoil->Fill( QGLikelihoodJetRecoil*jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight);
        if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
          h1_QGLikelihoodJetRecoil_MW400->Fill( QGLikelihoodJetRecoil, eventWeight);
          h1_QGLikelihoodProdRecoil_MW400->Fill( QGLikelihoodJetRecoil*jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight);
        }
      } else {
        h1_QGLikelihoodJetRecoil->Fill( -1., eventWeight);
        h1_QGLikelihoodProdRecoil->Fill( -1., eventWeight );
        if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
          h1_QGLikelihoodJetRecoil_MW400->Fill( -1., eventWeight);
          h1_QGLikelihoodProdRecoil_MW400->Fill( -1., eventWeight);
        }
      }
     
  
      if( jet1_selected.Pt()>jet2_selected.Pt() ) {
        h1_ptJet1->Fill( jet1_selected.Pt(), eventWeight );
        h1_ptJet2->Fill( jet2_selected.Pt(), eventWeight );
        h1_ptJet1_prekin->Fill( jet1_nokinfit.Pt(), eventWeight );
        h1_ptJet2_prekin->Fill( jet2_nokinfit.Pt(), eventWeight );
        h1_etaJet1->Fill( jet1_selected.Eta(), eventWeight );
        h1_etaJet2->Fill( jet2_selected.Eta(), eventWeight );
        h1_tcheJet1->Fill( jet1_selected.trackCountingHighEffBJetTag, eventWeight );
        h1_tcheJet2->Fill( jet2_selected.trackCountingHighEffBJetTag, eventWeight );
      } else {
        h1_ptJet1->Fill( jet2_selected.Pt(), eventWeight );
        h1_ptJet2->Fill( jet1_selected.Pt(), eventWeight );
        h1_ptJet1_prekin->Fill( jet2_nokinfit.Pt(), eventWeight );
        h1_ptJet2_prekin->Fill( jet1_nokinfit.Pt(), eventWeight );
        h1_etaJet1->Fill( jet2_selected.Eta(), eventWeight );
        h1_etaJet2->Fill( jet1_selected.Eta(), eventWeight );
        h1_tcheJet1->Fill( jet2_selected.trackCountingHighEffBJetTag, eventWeight );
        h1_tcheJet2->Fill( jet1_selected.trackCountingHighEffBJetTag, eventWeight );
      }
      h1_eMuonsJet1->Fill( jet1_selected.muonEnergyFraction, eventWeight );
      h1_eMuonsJet2->Fill( jet2_selected.muonEnergyFraction, eventWeight );
      h1_eElectronsJet1->Fill( jet1_selected.electronEnergyFraction, eventWeight );
      h1_eElectronsJet2->Fill( jet2_selected.electronEnergyFraction, eventWeight );
  
      h1_ptLept1->Fill( lept1.Pt(), eventWeight );
      h1_ptLept2->Fill( lept2.Pt(), eventWeight );
      h1_deltaRjj->Fill( jet1_selected.DeltaR(jet2_selected), eventWeight);
      h1_deltaRjj_prekin->Fill( jet1_nokinfit.DeltaR(jet2_nokinfit), eventWeight);
      h1_ptZll->Fill( diLepton.Pt(), eventWeight);
      h1_ptZjj->Fill( Zjj_kinfit.Pt(), eventWeight);
      if( leptType==0 )
        h1_mZmumu->Fill( diLepton.M(), eventWeight );
      else
        h1_mZee->Fill( diLepton.M(), eventWeight );
      h1_mZll->Fill( diLepton.M(), eventWeight);
  
  
      // fill QG plots only for the 0- and glue-tag category:
      if( maxBTag_found<=0 ) {
  
        h1_nChargedJet1->Fill(jet1_selected.nCharged, eventWeight);
        h1_nNeutralJet1->Fill(jet1_selected.nNeutral, eventWeight);
        h1_ptDJet1->Fill(jet1_selected.ptD, eventWeight);
      
        h1_nChargedJet2->Fill(jet2_selected.nCharged, eventWeight);
        h1_nNeutralJet2->Fill(jet2_selected.nNeutral, eventWeight);
        h1_ptDJet2->Fill(jet2_selected.ptD, eventWeight);
      
        h1_QGLikelihoodJet1->Fill( jet1_selected.QGLikelihood, eventWeight );
        h1_QGLikelihoodJet2->Fill( jet2_selected.QGLikelihood, eventWeight );
        h1_QGLikelihoodProd->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
  
        h1_QGLikelihoodNoPUJet1->Fill( jet1_selected.QGLikelihoodNoPU, eventWeight );
        h1_QGLikelihoodNoPUJet2->Fill( jet2_selected.QGLikelihoodNoPU, eventWeight );
        h1_QGLikelihoodNoPUProd->Fill( jet1_selected.QGLikelihoodNoPU*jet2_selected.QGLikelihoodNoPU, eventWeight );
  
        if( jet1_selected.Pt()>100. && jet1_selected.Pt()<123. ) h1_QGLikelihood_100_123->Fill( jet1_selected.QGLikelihood, eventWeight );
        if( jet2_selected.Pt()>100. && jet2_selected.Pt()<123. ) h1_QGLikelihood_100_123->Fill( jet2_selected.QGLikelihood, eventWeight );
        if( jet1_selected.Pt()>66. && jet1_selected.Pt()<81. ) h1_QGLikelihood_66_81->Fill( jet1_selected.QGLikelihood, eventWeight );
        if( jet2_selected.Pt()>66. && jet2_selected.Pt()<81. ) h1_QGLikelihood_66_81->Fill( jet2_selected.QGLikelihood, eventWeight );
  
  
        if( ZZ_kinfit.M()>0.94*300. && ZZ_kinfit.M()<1.1*300. ) {
          h1_QGLikelihoodJet1_MW300->Fill( jet1_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodJet2_MW300->Fill( jet2_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodProd_MW300->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
        }
        if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
          h1_QGLikelihoodJet1_MW400->Fill( jet1_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodJet2_MW400->Fill( jet2_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodProd_MW400->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
        }
        if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) {
          h1_QGLikelihoodJet1_MW500->Fill( jet1_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodJet2_MW500->Fill( jet2_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodProd_MW500->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
        }
  
  
        if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood < 0.1 ) {
          h1_mZZ_kinfit_hiMass_loQG->Fill(ZZ_kinfit.M(), eventWeight);
        } else {
          h1_mZZ_kinfit_hiMass_hiQG->Fill(ZZ_kinfit.M(), eventWeight);
        }
  
      } // if 0/glue tags
  
      h1_mZZ_nokinfit_hiMass_all->Fill( ZZ_nokinfit.M(), eventWeight);
      h1_mZZ_ZjjMassConstr_hiMass->Fill(ZZ_constr.M(), eventWeight);
      h1_mZZ_kinfit_hiMass_all->Fill( ZZ_kinfit.M(), eventWeight);
      if( nvertex>5 ) h1_mZZ_kinfit_hiMass_hiPU->Fill( ZZ_kinfit.M(), eventWeight);
      else h1_mZZ_kinfit_hiMass_loPU->Fill( ZZ_kinfit.M(), eventWeight);
      if( maxBTag_found>=0 ) h1_mZZ_kinfit_hiMass_nogluetag->Fill( ZZ_kinfit.M(), eventWeight);
      if( maxBTag_found==0 ) {
        h1_mZZ_kinfit_hiMass_0btag->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==0 )  h1_mZZ_kinfit_hiMass_0btag_MU->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==1 )  h1_mZZ_kinfit_hiMass_0btag_ELE->Fill( ZZ_kinfit.M(), eventWeight);
      } else if( maxBTag_found==1 ) {
        h1_mZZ_kinfit_hiMass_1btag->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==0 )  h1_mZZ_kinfit_hiMass_1btag_MU->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==1 )  h1_mZZ_kinfit_hiMass_1btag_ELE->Fill( ZZ_kinfit.M(), eventWeight);
      } else if( maxBTag_found==2 ) {
        h1_mZZ_kinfit_hiMass_2btag->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==0 )  h1_mZZ_kinfit_hiMass_2btag_MU->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==1 )  h1_mZZ_kinfit_hiMass_2btag_ELE->Fill( ZZ_kinfit.M(), eventWeight);
      } else if( maxBTag_found==-1 ) {
        h1_mZZ_kinfit_hiMass_gluetag->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==0 )  h1_mZZ_kinfit_hiMass_gluetag_MU->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==1 )  h1_mZZ_kinfit_hiMass_gluetag_ELE->Fill( ZZ_kinfit.M(), eventWeight);
      }
  
      h1_deltaRZmatching->Fill( bestDeltaRZ, eventWeight );
      if( maxBTag_found==0 && eventIsMatched ) h1_mZZ_kinfit_hiMass_0btag_matched->Fill( ZZ_kinfit.M(), eventWeight);
  
  
      h1_deltaRZZ->Fill(Zjj_nokinfit.DeltaR(diLepton), eventWeight);
  
      h1_ptZZ->Fill( ZZ_nokinfit.Pt(), eventWeight );
      h1_ptZZ_kinfit->Fill( ZZ_kinfit.Pt(), eventWeight );
      h1_etaZZ->Fill( ZZ_nokinfit.Eta(), eventWeight );
      h1_etaZZ_kinfit->Fill( ZZ_kinfit.Eta(), eventWeight );
  
      h1_helicityLD->Fill( helicityLD_selected, eventWeight );
      if( maxBTag_found>=0 ) h1_helicityLD_nogluetag->Fill( helicityLD_selected, eventWeight );
      h1_helicityLD_nokinfit->Fill( helicityLD_nokinfit_selected, eventWeight );
  
      h1_cosThetaStar->Fill(hangles_selected.helCosThetaStar, eventWeight);
      h1_cosTheta1->Fill(hangles_selected.helCosTheta1, eventWeight);
      h1_cosTheta2->Fill(hangles_selected.helCosTheta2, eventWeight);
      h1_phi->Fill(hangles_selected.helPhi, eventWeight);
      h1_phi1->Fill(hangles_selected.helPhi1, eventWeight);
  
  
  
      int partFlavor1=0;
      float deltaRmin1=999.;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        if( abs(pdgIdPart[iPart])>6 && pdgIdPart[iPart]!=21 ) continue;
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = jet1_selected.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin1 ) {
          partFlavor1 = pdgIdPart[iPart];
          deltaRmin1 = thisDeltaR;
        }
      }
      h1_deltaR_part1->Fill(deltaRmin1, eventWeight);
      h1_partFlavorJet1->Fill( partFlavor1, eventWeight );
      if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) h1_partFlavorJet1_MW400->Fill( partFlavor1, eventWeight );
      if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) h1_partFlavorJet1_MW500->Fill( partFlavor1, eventWeight );
      if( deltaRmin1<0.5 ) h1_partFlavorJet1_matched->Fill( partFlavor1, eventWeight );
      else h1_partFlavorJet1_notmatched->Fill( partFlavor1, eventWeight );
  
      if( maxBTag_found==0 ) h1_partFlavorJet1_0btag->Fill( partFlavor1, eventWeight );
      else if( maxBTag_found==1 ) h1_partFlavorJet1_1btag->Fill( partFlavor1, eventWeight );
      else if( maxBTag_found==2 ) h1_partFlavorJet1_2btag->Fill( partFlavor1, eventWeight );
      else if( maxBTag_found==-1 ) h1_partFlavorJet1_gluetag->Fill( partFlavor1, eventWeight );
      jet1_selected.pdgIdPart = partFlavor1;
      if( partFlavor1==21 ) h1_QGLikelihoodJet1_gluonMatched->Fill( jet1_selected.QGLikelihood, eventWeight );
      else if( abs(partFlavor1)<5 ) h1_QGLikelihoodJet1_quarkMatched->Fill( jet1_selected.QGLikelihood, eventWeight );
  
      float deltaRmin2=999.;
      int partFlavor2=0;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        if( abs(pdgIdPart[iPart])>6 && pdgIdPart[iPart]!=21 ) continue;
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = jet2_selected.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin2 ) {
          partFlavor2 = pdgIdPart[iPart];
          deltaRmin2 = thisDeltaR;
        }
      }
      h1_deltaR_part2->Fill(deltaRmin2, eventWeight);
      h1_partFlavorJet2->Fill( partFlavor2, eventWeight );
      if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) h1_partFlavorJet2_MW400->Fill( partFlavor2, eventWeight );
      if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) h1_partFlavorJet2_MW500->Fill( partFlavor2, eventWeight );
      if( deltaRmin2<0.5 ) h1_partFlavorJet2_matched->Fill( partFlavor2, eventWeight );
      else h1_partFlavorJet2_notmatched->Fill( partFlavor2, eventWeight );
  
      if( maxBTag_found==0 ) h1_partFlavorJet2_0btag->Fill( partFlavor2, eventWeight );
      else if( maxBTag_found==1 ) h1_partFlavorJet2_1btag->Fill( partFlavor2, eventWeight );
      else if( maxBTag_found==2 ) h1_partFlavorJet2_2btag->Fill( partFlavor2, eventWeight );
      else if( maxBTag_found==-1 ) h1_partFlavorJet2_gluetag->Fill( partFlavor2, eventWeight );
      jet2_selected.pdgIdPart = partFlavor2;
      if( partFlavor2==21 ) h1_QGLikelihoodJet2_gluonMatched->Fill( jet2_selected.QGLikelihood, eventWeight );
      else if( abs(partFlavor2)<5 ) h1_QGLikelihoodJet2_quarkMatched->Fill( jet2_selected.QGLikelihood, eventWeight );
  
  
      if( maxBTag_found<=0 ) {

        if( partFlavor1==21 || partFlavor2==21 ) {
          h1_QGLikelihoodProd_oneGluon->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          h1_mZZ_kinfit_hiMass_oneGluon->Fill( mZZ, eventWeight );
        }
        if( partFlavor1==21 && partFlavor2==21 ) {
          h1_QGLikelihoodProd_twoGluon->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          h1_mZZ_kinfit_hiMass_twoGluon->Fill( mZZ, eventWeight );
        }
        if( abs(partFlavor1)<5 && abs(partFlavor2)<5 ) {
          h1_QGLikelihoodProd_twoQuark->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          h1_mZZ_kinfit_hiMass_twoQuark->Fill( mZZ, eventWeight );
          if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
            h1_QGLikelihoodProd_twoQuark_MW400->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          }
          if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) {
            h1_QGLikelihoodProd_twoQuark_MW500->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          }
          if( maxBTag_found<=0 ) {
            h1_nEvents_partFlavor_beforeQG->Fill( 0., eventWeight );
            if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEvents_partFlavor_afterQG->Fill( 0., eventWeight );
            if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
              h1_nEventsMW400_partFlavor_beforeQG->Fill( 0., eventWeight );
              if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEventsMW400_partFlavor_afterQG->Fill( 0., eventWeight );
            }
            if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) {
              h1_nEventsMW500_partFlavor_beforeQG->Fill( 0., eventWeight );
              if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEventsMW500_partFlavor_afterQG->Fill( 0., eventWeight );
            }
          }
        } else { 
          if( maxBTag_found<=0 ) {
            h1_nEvents_partFlavor_beforeQG->Fill( 1., eventWeight );
            if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEvents_partFlavor_afterQG->Fill( 1., eventWeight );
            if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
              h1_nEventsMW400_partFlavor_beforeQG->Fill( 1., eventWeight );
              if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEventsMW400_partFlavor_afterQG->Fill( 1., eventWeight );
            }
            if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) {
              h1_nEventsMW500_partFlavor_beforeQG->Fill( 1., eventWeight );
              if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEventsMW500_partFlavor_afterQG->Fill( 1., eventWeight );
            }
          }
        }

      } //if 0/glue tag
  

    } //if signal region
  
  } //for entries


  std::cout << "nEvents_presel: " << nEvents_presel << std::endl;
  std::cout << "nEvents_presel_mZll: " << nEvents_presel_mZll << std::endl;
  std::cout << "nEvents_presel_mZll_mZjj: " << nEvents_presel_mZll_mZjj << std::endl;
  std::cout << "nEvents_presel_mZll_mZjj_0btag: " << nEvents_presel_mZll_mZjj_0btag << std::endl;
  std::cout << "nEvents_presel_mZll_mZjj_1btag: " << nEvents_presel_mZll_mZjj_1btag << std::endl;
  std::cout << "nEvents_presel_mZll_mZjj_2btag: " << nEvents_presel_mZll_mZjj_2btag << std::endl;

  h1_nCounter->SetBinContent(1, nCounter_);
  h1_nCounterW->SetBinContent(1, nCounterW_);
  h1_nCounterPU->SetBinContent(1, nCounterPU_);

  float eff_gluetag_250 = nEventsPassed_fb_gluetag_250/nCounterW_;
  float eff_0btag_250 = nEventsPassed_fb_0btag_250/nCounterW_;
  float eff_1btag_250 = nEventsPassed_fb_1btag_250/nCounterW_;
  float eff_2btag_250 = nEventsPassed_fb_2btag_250/nCounterW_;

  float eff_gluetag_250_ELE = nEventsPassed_fb_gluetag_250_ELE/nCounterW_;
  float eff_0btag_250_ELE = nEventsPassed_fb_0btag_250_ELE/nCounterW_;
  float eff_1btag_250_ELE = nEventsPassed_fb_1btag_250_ELE/nCounterW_;
  float eff_2btag_250_ELE = nEventsPassed_fb_2btag_250_ELE/nCounterW_;

  float eff_gluetag_250_MU = nEventsPassed_fb_gluetag_250_MU/nCounterW_;
  float eff_0btag_250_MU = nEventsPassed_fb_0btag_250_MU/nCounterW_;
  float eff_1btag_250_MU = nEventsPassed_fb_1btag_250_MU/nCounterW_;
  float eff_2btag_250_MU = nEventsPassed_fb_2btag_250_MU/nCounterW_;

  float eff_gluetag_300 = nEventsPassed_fb_gluetag_300/nCounterW_;
  float eff_0btag_300 = nEventsPassed_fb_0btag_300/nCounterW_;
  float eff_1btag_300 = nEventsPassed_fb_1btag_300/nCounterW_;
  float eff_2btag_300 = nEventsPassed_fb_2btag_300/nCounterW_;

  float eff_gluetag_300_ELE = nEventsPassed_fb_gluetag_300_ELE/nCounterW_;
  float eff_0btag_300_ELE = nEventsPassed_fb_0btag_300_ELE/nCounterW_;
  float eff_1btag_300_ELE = nEventsPassed_fb_1btag_300_ELE/nCounterW_;
  float eff_2btag_300_ELE = nEventsPassed_fb_2btag_300_ELE/nCounterW_;

  float eff_gluetag_300_MU = nEventsPassed_fb_gluetag_300_MU/nCounterW_;
  float eff_0btag_300_MU = nEventsPassed_fb_0btag_300_MU/nCounterW_;
  float eff_1btag_300_MU = nEventsPassed_fb_1btag_300_MU/nCounterW_;
  float eff_2btag_300_MU = nEventsPassed_fb_2btag_300_MU/nCounterW_;

  float eff_gluetag_350 = nEventsPassed_fb_gluetag_350/nCounterW_;
  float eff_0btag_350 = nEventsPassed_fb_0btag_350/nCounterW_;
  float eff_1btag_350 = nEventsPassed_fb_1btag_350/nCounterW_;
  float eff_2btag_350 = nEventsPassed_fb_2btag_350/nCounterW_;

  float eff_gluetag_350_ELE = nEventsPassed_fb_gluetag_350_ELE/nCounterW_;
  float eff_0btag_350_ELE = nEventsPassed_fb_0btag_350_ELE/nCounterW_;
  float eff_1btag_350_ELE = nEventsPassed_fb_1btag_350_ELE/nCounterW_;
  float eff_2btag_350_ELE = nEventsPassed_fb_2btag_350_ELE/nCounterW_;

  float eff_gluetag_350_MU = nEventsPassed_fb_gluetag_350_MU/nCounterW_;
  float eff_0btag_350_MU = nEventsPassed_fb_0btag_350_MU/nCounterW_;
  float eff_1btag_350_MU = nEventsPassed_fb_1btag_350_MU/nCounterW_;
  float eff_2btag_350_MU = nEventsPassed_fb_2btag_350_MU/nCounterW_;

  float eff_gluetag_400 = nEventsPassed_fb_gluetag_400/nCounterW_;
  float eff_0btag_400 = nEventsPassed_fb_0btag_400/nCounterW_;
  float eff_1btag_400 = nEventsPassed_fb_1btag_400/nCounterW_;
  float eff_2btag_400 = nEventsPassed_fb_2btag_400/nCounterW_;

  float eff_gluetag_400_ELE = nEventsPassed_fb_gluetag_400_ELE/nCounterW_;
  float eff_0btag_400_ELE = nEventsPassed_fb_0btag_400_ELE/nCounterW_;
  float eff_1btag_400_ELE = nEventsPassed_fb_1btag_400_ELE/nCounterW_;
  float eff_2btag_400_ELE = nEventsPassed_fb_2btag_400_ELE/nCounterW_;

  float eff_gluetag_400_MU = nEventsPassed_fb_gluetag_400_MU/nCounterW_;
  float eff_0btag_400_MU = nEventsPassed_fb_0btag_400_MU/nCounterW_;
  float eff_1btag_400_MU = nEventsPassed_fb_1btag_400_MU/nCounterW_;
  float eff_2btag_400_MU = nEventsPassed_fb_2btag_400_MU/nCounterW_;

  float eff_gluetag_450 = nEventsPassed_fb_gluetag_450/nCounterW_;
  float eff_0btag_450 = nEventsPassed_fb_0btag_450/nCounterW_;
  float eff_1btag_450 = nEventsPassed_fb_1btag_450/nCounterW_;
  float eff_2btag_450 = nEventsPassed_fb_2btag_450/nCounterW_;

  float eff_gluetag_450_ELE = nEventsPassed_fb_gluetag_450_ELE/nCounterW_;
  float eff_0btag_450_ELE = nEventsPassed_fb_0btag_450_ELE/nCounterW_;
  float eff_1btag_450_ELE = nEventsPassed_fb_1btag_450_ELE/nCounterW_;
  float eff_2btag_450_ELE = nEventsPassed_fb_2btag_450_ELE/nCounterW_;

  float eff_gluetag_450_MU = nEventsPassed_fb_gluetag_450_MU/nCounterW_;
  float eff_0btag_450_MU = nEventsPassed_fb_0btag_450_MU/nCounterW_;
  float eff_1btag_450_MU = nEventsPassed_fb_1btag_450_MU/nCounterW_;
  float eff_2btag_450_MU = nEventsPassed_fb_2btag_450_MU/nCounterW_;

  float eff_gluetag_500 = nEventsPassed_fb_gluetag_500/nCounterW_;
  float eff_0btag_500 = nEventsPassed_fb_0btag_500/nCounterW_;
  float eff_1btag_500 = nEventsPassed_fb_1btag_500/nCounterW_;
  float eff_2btag_500 = nEventsPassed_fb_2btag_500/nCounterW_;

  float eff_gluetag_500_ELE = nEventsPassed_fb_gluetag_500_ELE/nCounterW_;
  float eff_0btag_500_ELE = nEventsPassed_fb_0btag_500_ELE/nCounterW_;
  float eff_1btag_500_ELE = nEventsPassed_fb_1btag_500_ELE/nCounterW_;
  float eff_2btag_500_ELE = nEventsPassed_fb_2btag_500_ELE/nCounterW_;

  float eff_gluetag_500_MU = nEventsPassed_fb_gluetag_500_MU/nCounterW_;
  float eff_0btag_500_MU = nEventsPassed_fb_0btag_500_MU/nCounterW_;
  float eff_1btag_500_MU = nEventsPassed_fb_1btag_500_MU/nCounterW_;
  float eff_2btag_500_MU = nEventsPassed_fb_2btag_500_MU/nCounterW_;

  float eff_gluetag_600 = nEventsPassed_fb_gluetag_600/nCounterW_;
  float eff_0btag_600 = nEventsPassed_fb_0btag_600/nCounterW_;
  float eff_1btag_600 = nEventsPassed_fb_1btag_600/nCounterW_;
  float eff_2btag_600 = nEventsPassed_fb_2btag_600/nCounterW_;

  float eff_gluetag_600_ELE = nEventsPassed_fb_gluetag_600_ELE/nCounterW_;
  float eff_0btag_600_ELE = nEventsPassed_fb_0btag_600_ELE/nCounterW_;
  float eff_1btag_600_ELE = nEventsPassed_fb_1btag_600_ELE/nCounterW_;
  float eff_2btag_600_ELE = nEventsPassed_fb_2btag_600_ELE/nCounterW_;

  float eff_gluetag_600_MU = nEventsPassed_fb_gluetag_600_MU/nCounterW_;
  float eff_0btag_600_MU = nEventsPassed_fb_0btag_600_MU/nCounterW_;
  float eff_1btag_600_MU = nEventsPassed_fb_1btag_600_MU/nCounterW_;
  float eff_2btag_600_MU = nEventsPassed_fb_2btag_600_MU/nCounterW_;


  h1_nEvents_fb_gluetag_250->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_250);
  h1_nEvents_fb_0btag_250->SetBinContent(1,1000.*nEventsPassed_fb_0btag_250);
  h1_nEvents_fb_1btag_250->SetBinContent(1,1000.*nEventsPassed_fb_1btag_250);
  h1_nEvents_fb_2btag_250->SetBinContent(1,1000.*nEventsPassed_fb_2btag_250);

  h1_nEvents_fb_gluetag_250_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_250_ELE);
  h1_nEvents_fb_0btag_250_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_250_ELE);
  h1_nEvents_fb_1btag_250_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_250_ELE);
  h1_nEvents_fb_2btag_250_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_250_ELE);

  h1_nEvents_fb_gluetag_250_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_250_MU);
  h1_nEvents_fb_0btag_250_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_250_MU);
  h1_nEvents_fb_1btag_250_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_250_MU);
  h1_nEvents_fb_2btag_250_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_250_MU);

  h1_eff_gluetag_250->SetBinContent(1,eff_gluetag_250);
  h1_eff_0btag_250->SetBinContent(1,eff_0btag_250);
  h1_eff_1btag_250->SetBinContent(1,eff_1btag_250);
  h1_eff_2btag_250->SetBinContent(1,eff_2btag_250);

  h1_eff_gluetag_250_ELE->SetBinContent(1,eff_gluetag_250_ELE);
  h1_eff_0btag_250_ELE->SetBinContent(1,eff_0btag_250_ELE);
  h1_eff_1btag_250_ELE->SetBinContent(1,eff_1btag_250_ELE);
  h1_eff_2btag_250_ELE->SetBinContent(1,eff_2btag_250_ELE);

  h1_eff_gluetag_250_MU->SetBinContent(1,eff_gluetag_250_MU);
  h1_eff_0btag_250_MU->SetBinContent(1,eff_0btag_250_MU);
  h1_eff_1btag_250_MU->SetBinContent(1,eff_1btag_250_MU);
  h1_eff_2btag_250_MU->SetBinContent(1,eff_2btag_250_MU);


  h1_nEvents_fb_gluetag_300->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_300);
  h1_nEvents_fb_0btag_300->SetBinContent(1,1000.*nEventsPassed_fb_0btag_300);
  h1_nEvents_fb_1btag_300->SetBinContent(1,1000.*nEventsPassed_fb_1btag_300);
  h1_nEvents_fb_2btag_300->SetBinContent(1,1000.*nEventsPassed_fb_2btag_300);

  h1_nEvents_fb_gluetag_300_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_300_ELE);
  h1_nEvents_fb_0btag_300_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_300_ELE);
  h1_nEvents_fb_1btag_300_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_300_ELE);
  h1_nEvents_fb_2btag_300_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_300_ELE);

  h1_nEvents_fb_gluetag_300_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_300_MU);
  h1_nEvents_fb_0btag_300_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_300_MU);
  h1_nEvents_fb_1btag_300_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_300_MU);
  h1_nEvents_fb_2btag_300_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_300_MU);

  h1_eff_gluetag_300->SetBinContent(1,eff_gluetag_300);
  h1_eff_0btag_300->SetBinContent(1,eff_0btag_300);
  h1_eff_1btag_300->SetBinContent(1,eff_1btag_300);
  h1_eff_2btag_300->SetBinContent(1,eff_2btag_300);

  h1_eff_gluetag_300_ELE->SetBinContent(1,eff_gluetag_300_ELE);
  h1_eff_0btag_300_ELE->SetBinContent(1,eff_0btag_300_ELE);
  h1_eff_1btag_300_ELE->SetBinContent(1,eff_1btag_300_ELE);
  h1_eff_2btag_300_ELE->SetBinContent(1,eff_2btag_300_ELE);

  h1_eff_gluetag_300_MU->SetBinContent(1,eff_gluetag_300_MU);
  h1_eff_0btag_300_MU->SetBinContent(1,eff_0btag_300_MU);
  h1_eff_1btag_300_MU->SetBinContent(1,eff_1btag_300_MU);
  h1_eff_2btag_300_MU->SetBinContent(1,eff_2btag_300_MU);


  h1_nEvents_fb_gluetag_350->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_350);
  h1_nEvents_fb_0btag_350->SetBinContent(1,1000.*nEventsPassed_fb_0btag_350);
  h1_nEvents_fb_1btag_350->SetBinContent(1,1000.*nEventsPassed_fb_1btag_350);
  h1_nEvents_fb_2btag_350->SetBinContent(1,1000.*nEventsPassed_fb_2btag_350);

  h1_nEvents_fb_gluetag_350_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_350_ELE);
  h1_nEvents_fb_0btag_350_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_350_ELE);
  h1_nEvents_fb_1btag_350_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_350_ELE);
  h1_nEvents_fb_2btag_350_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_350_ELE);

  h1_nEvents_fb_gluetag_350_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_350_MU);
  h1_nEvents_fb_0btag_350_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_350_MU);
  h1_nEvents_fb_1btag_350_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_350_MU);
  h1_nEvents_fb_2btag_350_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_350_MU);

  h1_eff_gluetag_350->SetBinContent(1,eff_gluetag_350);
  h1_eff_0btag_350->SetBinContent(1,eff_0btag_350);
  h1_eff_1btag_350->SetBinContent(1,eff_1btag_350);
  h1_eff_2btag_350->SetBinContent(1,eff_2btag_350);

  h1_eff_gluetag_350_ELE->SetBinContent(1,eff_gluetag_350_ELE);
  h1_eff_0btag_350_ELE->SetBinContent(1,eff_0btag_350_ELE);
  h1_eff_1btag_350_ELE->SetBinContent(1,eff_1btag_350_ELE);
  h1_eff_2btag_350_ELE->SetBinContent(1,eff_2btag_350_ELE);

  h1_eff_gluetag_350_MU->SetBinContent(1,eff_gluetag_350_MU);
  h1_eff_0btag_350_MU->SetBinContent(1,eff_0btag_350_MU);
  h1_eff_1btag_350_MU->SetBinContent(1,eff_1btag_350_MU);
  h1_eff_2btag_350_MU->SetBinContent(1,eff_2btag_350_MU);


  h1_nEvents_fb_gluetag_400->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_400);
  h1_nEvents_fb_0btag_400->SetBinContent(1,1000.*nEventsPassed_fb_0btag_400);
  h1_nEvents_fb_1btag_400->SetBinContent(1,1000.*nEventsPassed_fb_1btag_400);
  h1_nEvents_fb_2btag_400->SetBinContent(1,1000.*nEventsPassed_fb_2btag_400);

  h1_nEvents_fb_gluetag_400_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_400_ELE);
  h1_nEvents_fb_0btag_400_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_400_ELE);
  h1_nEvents_fb_1btag_400_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_400_ELE);
  h1_nEvents_fb_2btag_400_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_400_ELE);

  h1_nEvents_fb_gluetag_400_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_400_MU);
  h1_nEvents_fb_0btag_400_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_400_MU);
  h1_nEvents_fb_1btag_400_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_400_MU);
  h1_nEvents_fb_2btag_400_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_400_MU);

  h1_eff_gluetag_400->SetBinContent(1,eff_gluetag_400);
  h1_eff_0btag_400->SetBinContent(1,eff_0btag_400);
  h1_eff_1btag_400->SetBinContent(1,eff_1btag_400);
  h1_eff_2btag_400->SetBinContent(1,eff_2btag_400);

  h1_eff_gluetag_400_ELE->SetBinContent(1,eff_gluetag_400_ELE);
  h1_eff_0btag_400_ELE->SetBinContent(1,eff_0btag_400_ELE);
  h1_eff_1btag_400_ELE->SetBinContent(1,eff_1btag_400_ELE);
  h1_eff_2btag_400_ELE->SetBinContent(1,eff_2btag_400_ELE);

  h1_eff_gluetag_400_MU->SetBinContent(1,eff_gluetag_400_MU);
  h1_eff_0btag_400_MU->SetBinContent(1,eff_0btag_400_MU);
  h1_eff_1btag_400_MU->SetBinContent(1,eff_1btag_400_MU);
  h1_eff_2btag_400_MU->SetBinContent(1,eff_2btag_400_MU);


  h1_nEvents_fb_gluetag_450->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_450);
  h1_nEvents_fb_0btag_450->SetBinContent(1,1000.*nEventsPassed_fb_0btag_450);
  h1_nEvents_fb_1btag_450->SetBinContent(1,1000.*nEventsPassed_fb_1btag_450);
  h1_nEvents_fb_2btag_450->SetBinContent(1,1000.*nEventsPassed_fb_2btag_450);

  h1_nEvents_fb_gluetag_450_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_450_ELE);
  h1_nEvents_fb_0btag_450_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_450_ELE);
  h1_nEvents_fb_1btag_450_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_450_ELE);
  h1_nEvents_fb_2btag_450_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_450_ELE);

  h1_nEvents_fb_gluetag_450_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_450_MU);
  h1_nEvents_fb_0btag_450_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_450_MU);
  h1_nEvents_fb_1btag_450_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_450_MU);
  h1_nEvents_fb_2btag_450_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_450_MU);

  h1_eff_gluetag_450->SetBinContent(1,eff_gluetag_450);
  h1_eff_0btag_450->SetBinContent(1,eff_0btag_450);
  h1_eff_1btag_450->SetBinContent(1,eff_1btag_450);
  h1_eff_2btag_450->SetBinContent(1,eff_2btag_450);

  h1_eff_gluetag_450_ELE->SetBinContent(1,eff_gluetag_450_ELE);
  h1_eff_0btag_450_ELE->SetBinContent(1,eff_0btag_450_ELE);
  h1_eff_1btag_450_ELE->SetBinContent(1,eff_1btag_450_ELE);
  h1_eff_2btag_450_ELE->SetBinContent(1,eff_2btag_450_ELE);

  h1_eff_gluetag_450_MU->SetBinContent(1,eff_gluetag_450_MU);
  h1_eff_0btag_450_MU->SetBinContent(1,eff_0btag_450_MU);
  h1_eff_1btag_450_MU->SetBinContent(1,eff_1btag_450_MU);
  h1_eff_2btag_450_MU->SetBinContent(1,eff_2btag_450_MU);


  h1_nEvents_fb_gluetag_500->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_500);
  h1_nEvents_fb_0btag_500->SetBinContent(1,1000.*nEventsPassed_fb_0btag_500);
  h1_nEvents_fb_1btag_500->SetBinContent(1,1000.*nEventsPassed_fb_1btag_500);
  h1_nEvents_fb_2btag_500->SetBinContent(1,1000.*nEventsPassed_fb_2btag_500);

  h1_nEvents_fb_gluetag_500_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_500_ELE);
  h1_nEvents_fb_0btag_500_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_500_ELE);
  h1_nEvents_fb_1btag_500_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_500_ELE);
  h1_nEvents_fb_2btag_500_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_500_ELE);

  h1_nEvents_fb_gluetag_500_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_500_MU);
  h1_nEvents_fb_0btag_500_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_500_MU);
  h1_nEvents_fb_1btag_500_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_500_MU);
  h1_nEvents_fb_2btag_500_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_500_MU);

  h1_eff_gluetag_500->SetBinContent(1,eff_gluetag_500);
  h1_eff_0btag_500->SetBinContent(1,eff_0btag_500);
  h1_eff_1btag_500->SetBinContent(1,eff_1btag_500);
  h1_eff_2btag_500->SetBinContent(1,eff_2btag_500);

  h1_eff_gluetag_500_ELE->SetBinContent(1,eff_gluetag_500_ELE);
  h1_eff_0btag_500_ELE->SetBinContent(1,eff_0btag_500_ELE);
  h1_eff_1btag_500_ELE->SetBinContent(1,eff_1btag_500_ELE);
  h1_eff_2btag_500_ELE->SetBinContent(1,eff_2btag_500_ELE);

  h1_eff_gluetag_500_MU->SetBinContent(1,eff_gluetag_500_MU);
  h1_eff_0btag_500_MU->SetBinContent(1,eff_0btag_500_MU);
  h1_eff_1btag_500_MU->SetBinContent(1,eff_1btag_500_MU);
  h1_eff_2btag_500_MU->SetBinContent(1,eff_2btag_500_MU);


  h1_nEvents_fb_gluetag_600->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_600);
  h1_nEvents_fb_0btag_600->SetBinContent(1,1000.*nEventsPassed_fb_0btag_600);
  h1_nEvents_fb_1btag_600->SetBinContent(1,1000.*nEventsPassed_fb_1btag_600);
  h1_nEvents_fb_2btag_600->SetBinContent(1,1000.*nEventsPassed_fb_2btag_600);

  h1_nEvents_fb_gluetag_600_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_600_ELE);
  h1_nEvents_fb_0btag_600_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_600_ELE);
  h1_nEvents_fb_1btag_600_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_600_ELE);
  h1_nEvents_fb_2btag_600_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_600_ELE);

  h1_nEvents_fb_gluetag_600_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_600_MU);
  h1_nEvents_fb_0btag_600_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_600_MU);
  h1_nEvents_fb_1btag_600_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_600_MU);
  h1_nEvents_fb_2btag_600_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_600_MU);

  h1_eff_gluetag_600->SetBinContent(1,eff_gluetag_600);
  h1_eff_0btag_600->SetBinContent(1,eff_0btag_600);
  h1_eff_1btag_600->SetBinContent(1,eff_1btag_600);
  h1_eff_2btag_600->SetBinContent(1,eff_2btag_600);

  h1_eff_gluetag_600_ELE->SetBinContent(1,eff_gluetag_600_ELE);
  h1_eff_0btag_600_ELE->SetBinContent(1,eff_0btag_600_ELE);
  h1_eff_1btag_600_ELE->SetBinContent(1,eff_1btag_600_ELE);
  h1_eff_2btag_600_ELE->SetBinContent(1,eff_2btag_600_ELE);

  h1_eff_gluetag_600_MU->SetBinContent(1,eff_gluetag_600_MU);
  h1_eff_0btag_600_MU->SetBinContent(1,eff_0btag_600_MU);
  h1_eff_1btag_600_MU->SetBinContent(1,eff_1btag_600_MU);
  h1_eff_2btag_600_MU->SetBinContent(1,eff_2btag_600_MU);



  std::cout << std::endl << std::endl;
  std::cout << "----> DATASET: " << dataset_ << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << std::endl;
  if( isMC ) std::cout << "----> PU REWEIGHING: " << PUType_ << std::endl;
  std::cout << std::endl << std::endl;
  std::cout << "----> 250 GeV (235-275): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_250 << " ev/fb-1  (" << nEventsPassed_0btag_250 << " events)" << " Efficiency: " << 100.*eff_0btag_250 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_250 << " ev/fb-1  (" << nEventsPassed_1btag_250 << " events)" << " Efficiency: " << 100.*eff_1btag_250 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_250 << " ev/fb-1  (" << nEventsPassed_2btag_250 << " events)" << " Efficiency: " << 100.*eff_2btag_250 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 300 GeV (282-330): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_300 << " ev/fb-1  (" << nEventsPassed_0btag_300 << " events)" << " Efficiency: " << 100.*eff_0btag_300 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_300 << " ev/fb-1  (" << nEventsPassed_1btag_300 << " events)" << " Efficiency: " << 100.*eff_1btag_300 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_300 << " ev/fb-1  (" << nEventsPassed_2btag_300 << " events)" << " Efficiency: " << 100.*eff_2btag_300 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 350 GeV (329-385): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_350 << " ev/fb-1  (" << nEventsPassed_0btag_350 << " events)" << " Efficiency: " << 100.*eff_0btag_350 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_350 << " ev/fb-1  (" << nEventsPassed_1btag_350 << " events)" << " Efficiency: " << 100.*eff_1btag_350 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_350 << " ev/fb-1  (" << nEventsPassed_2btag_350 << " events)" << " Efficiency: " << 100.*eff_2btag_350 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 400 GeV (376-440): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_400 << " ev/fb-1  (" << nEventsPassed_0btag_400 << " events)" << " Efficiency: " << 100.*eff_0btag_400 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_400 << " ev/fb-1  (" << nEventsPassed_1btag_400 << " events)" << " Efficiency: " << 100.*eff_1btag_400 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_400 << " ev/fb-1  (" << nEventsPassed_2btag_400 << " events)" << " Efficiency: " << 100.*eff_2btag_400 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 450 GeV (423-495): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_450 << " ev/fb-1  (" << nEventsPassed_0btag_450 << " events)" << " Efficiency: " << 100.*eff_0btag_450 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_450 << " ev/fb-1  (" << nEventsPassed_1btag_450 << " events)" << " Efficiency: " << 100.*eff_1btag_450 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_450 << " ev/fb-1  (" << nEventsPassed_2btag_450 << " events)" << " Efficiency: " << 100.*eff_2btag_450 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 500 GeV (470-550): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_500 << " ev/fb-1  (" << nEventsPassed_0btag_500 << " events)" << " Efficiency: " << 100.*eff_0btag_500 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_500 << " ev/fb-1  (" << nEventsPassed_1btag_500 << " events)" << " Efficiency: " << 100.*eff_1btag_500 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_500 << " ev/fb-1  (" << nEventsPassed_2btag_500 << " events)" << " Efficiency: " << 100.*eff_2btag_500 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 600 GeV (564-660): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_600 << " ev/fb-1  (" << nEventsPassed_0btag_600 << " events)" << " Efficiency: " << 100.*eff_0btag_600 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_600 << " ev/fb-1  (" << nEventsPassed_1btag_600 << " events)" << " Efficiency: " << 100.*eff_1btag_600 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_600 << " ev/fb-1  (" << nEventsPassed_2btag_600 << " events)" << " Efficiency: " << 100.*eff_2btag_600 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << std::endl;



  outFile_->cd();

  tree_passedEvents->Write();

  h1_nCounter->Write();
  h1_nCounterW->Write();
  h1_nCounterPU->Write();

  h1_nEventsCategories_presel->Write();

  h1_nEvents_fb_gluetag_250->Write();
  h1_nEvents_fb_0btag_250->Write();
  h1_nEvents_fb_1btag_250->Write();
  h1_nEvents_fb_2btag_250->Write();

  h1_nEvents_fb_gluetag_250_ELE->Write();
  h1_nEvents_fb_0btag_250_ELE->Write();
  h1_nEvents_fb_1btag_250_ELE->Write();
  h1_nEvents_fb_2btag_250_ELE->Write();

  h1_nEvents_fb_gluetag_250_MU->Write();
  h1_nEvents_fb_0btag_250_MU->Write();
  h1_nEvents_fb_1btag_250_MU->Write();
  h1_nEvents_fb_2btag_250_MU->Write();

  h1_eff_gluetag_250->Write();
  h1_eff_0btag_250->Write();
  h1_eff_1btag_250->Write();
  h1_eff_2btag_250->Write();

  h1_eff_gluetag_250_ELE->Write();
  h1_eff_0btag_250_ELE->Write();
  h1_eff_1btag_250_ELE->Write();
  h1_eff_2btag_250_ELE->Write();

  h1_eff_gluetag_250_MU->Write();
  h1_eff_0btag_250_MU->Write();
  h1_eff_1btag_250_MU->Write();
  h1_eff_2btag_250_MU->Write();

  h1_nEvents_fb_gluetag_300->Write();
  h1_nEvents_fb_0btag_300->Write();
  h1_nEvents_fb_1btag_300->Write();
  h1_nEvents_fb_2btag_300->Write();

  h1_nEvents_fb_gluetag_300_ELE->Write();
  h1_nEvents_fb_0btag_300_ELE->Write();
  h1_nEvents_fb_1btag_300_ELE->Write();
  h1_nEvents_fb_2btag_300_ELE->Write();

  h1_nEvents_fb_gluetag_300_MU->Write();
  h1_nEvents_fb_0btag_300_MU->Write();
  h1_nEvents_fb_1btag_300_MU->Write();
  h1_nEvents_fb_2btag_300_MU->Write();

  h1_eff_gluetag_300->Write();
  h1_eff_0btag_300->Write();
  h1_eff_1btag_300->Write();
  h1_eff_2btag_300->Write();

  h1_eff_gluetag_300_ELE->Write();
  h1_eff_0btag_300_ELE->Write();
  h1_eff_1btag_300_ELE->Write();
  h1_eff_2btag_300_ELE->Write();

  h1_eff_gluetag_300_MU->Write();
  h1_eff_0btag_300_MU->Write();
  h1_eff_1btag_300_MU->Write();
  h1_eff_2btag_300_MU->Write();

  h1_nEvents_fb_gluetag_350->Write();
  h1_nEvents_fb_0btag_350->Write();
  h1_nEvents_fb_1btag_350->Write();
  h1_nEvents_fb_2btag_350->Write();

  h1_nEvents_fb_gluetag_350_ELE->Write();
  h1_nEvents_fb_0btag_350_ELE->Write();
  h1_nEvents_fb_1btag_350_ELE->Write();
  h1_nEvents_fb_2btag_350_ELE->Write();

  h1_nEvents_fb_gluetag_350_MU->Write();
  h1_nEvents_fb_0btag_350_MU->Write();
  h1_nEvents_fb_1btag_350_MU->Write();
  h1_nEvents_fb_2btag_350_MU->Write();

  h1_eff_gluetag_350->Write();
  h1_eff_0btag_350->Write();
  h1_eff_1btag_350->Write();
  h1_eff_2btag_350->Write();

  h1_eff_gluetag_350_ELE->Write();
  h1_eff_0btag_350_ELE->Write();
  h1_eff_1btag_350_ELE->Write();
  h1_eff_2btag_350_ELE->Write();

  h1_eff_gluetag_350_MU->Write();
  h1_eff_0btag_350_MU->Write();
  h1_eff_1btag_350_MU->Write();
  h1_eff_2btag_350_MU->Write();

  h1_nEvents_fb_gluetag_400->Write();
  h1_nEvents_fb_0btag_400->Write();
  h1_nEvents_fb_1btag_400->Write();
  h1_nEvents_fb_2btag_400->Write();

  h1_nEvents_fb_gluetag_400_ELE->Write();
  h1_nEvents_fb_0btag_400_ELE->Write();
  h1_nEvents_fb_1btag_400_ELE->Write();
  h1_nEvents_fb_2btag_400_ELE->Write();

  h1_nEvents_fb_gluetag_400_MU->Write();
  h1_nEvents_fb_0btag_400_MU->Write();
  h1_nEvents_fb_1btag_400_MU->Write();
  h1_nEvents_fb_2btag_400_MU->Write();

  h1_eff_gluetag_400->Write();
  h1_eff_0btag_400->Write();
  h1_eff_1btag_400->Write();
  h1_eff_2btag_400->Write();

  h1_eff_gluetag_400_ELE->Write();
  h1_eff_0btag_400_ELE->Write();
  h1_eff_1btag_400_ELE->Write();
  h1_eff_2btag_400_ELE->Write();

  h1_eff_gluetag_400_MU->Write();
  h1_eff_0btag_400_MU->Write();
  h1_eff_1btag_400_MU->Write();
  h1_eff_2btag_400_MU->Write();

  h1_nEvents_fb_gluetag_450->Write();
  h1_nEvents_fb_0btag_450->Write();
  h1_nEvents_fb_1btag_450->Write();
  h1_nEvents_fb_2btag_450->Write();

  h1_nEvents_fb_gluetag_450_ELE->Write();
  h1_nEvents_fb_0btag_450_ELE->Write();
  h1_nEvents_fb_1btag_450_ELE->Write();
  h1_nEvents_fb_2btag_450_ELE->Write();

  h1_nEvents_fb_gluetag_450_MU->Write();
  h1_nEvents_fb_0btag_450_MU->Write();
  h1_nEvents_fb_1btag_450_MU->Write();
  h1_nEvents_fb_2btag_450_MU->Write();

  h1_eff_gluetag_450->Write();
  h1_eff_0btag_450->Write();
  h1_eff_1btag_450->Write();
  h1_eff_2btag_450->Write();

  h1_eff_gluetag_450_ELE->Write();
  h1_eff_0btag_450_ELE->Write();
  h1_eff_1btag_450_ELE->Write();
  h1_eff_2btag_450_ELE->Write();

  h1_eff_gluetag_450_MU->Write();
  h1_eff_0btag_450_MU->Write();
  h1_eff_1btag_450_MU->Write();
  h1_eff_2btag_450_MU->Write();

  h1_nEvents_fb_gluetag_500->Write();
  h1_nEvents_fb_0btag_500->Write();
  h1_nEvents_fb_1btag_500->Write();
  h1_nEvents_fb_2btag_500->Write();

  h1_nEvents_fb_gluetag_500_ELE->Write();
  h1_nEvents_fb_0btag_500_ELE->Write();
  h1_nEvents_fb_1btag_500_ELE->Write();
  h1_nEvents_fb_2btag_500_ELE->Write();

  h1_nEvents_fb_gluetag_500_MU->Write();
  h1_nEvents_fb_0btag_500_MU->Write();
  h1_nEvents_fb_1btag_500_MU->Write();
  h1_nEvents_fb_2btag_500_MU->Write();

  h1_eff_gluetag_500->Write();
  h1_eff_0btag_500->Write();
  h1_eff_1btag_500->Write();
  h1_eff_2btag_500->Write();

  h1_eff_gluetag_500_ELE->Write();
  h1_eff_0btag_500_ELE->Write();
  h1_eff_1btag_500_ELE->Write();
  h1_eff_2btag_500_ELE->Write();

  h1_eff_gluetag_500_MU->Write();
  h1_eff_0btag_500_MU->Write();
  h1_eff_1btag_500_MU->Write();
  h1_eff_2btag_500_MU->Write();

  h1_nEvents_fb_gluetag_600->Write();
  h1_nEvents_fb_0btag_600->Write();
  h1_nEvents_fb_1btag_600->Write();
  h1_nEvents_fb_2btag_600->Write();

  h1_nEvents_fb_gluetag_600_ELE->Write();
  h1_nEvents_fb_0btag_600_ELE->Write();
  h1_nEvents_fb_1btag_600_ELE->Write();
  h1_nEvents_fb_2btag_600_ELE->Write();

  h1_nEvents_fb_gluetag_600_MU->Write();
  h1_nEvents_fb_0btag_600_MU->Write();
  h1_nEvents_fb_1btag_600_MU->Write();
  h1_nEvents_fb_2btag_600_MU->Write();

  h1_eff_gluetag_600->Write();
  h1_eff_0btag_600->Write();
  h1_eff_1btag_600->Write();
  h1_eff_2btag_600->Write();

  h1_eff_gluetag_600_ELE->Write();
  h1_eff_0btag_600_ELE->Write();
  h1_eff_1btag_600_ELE->Write();
  h1_eff_2btag_600_ELE->Write();

  h1_eff_gluetag_600_MU->Write();
  h1_eff_0btag_600_MU->Write();
  h1_eff_1btag_600_MU->Write();
  h1_eff_2btag_600_MU->Write();


  h1_run->Write();

  h1_nvertex->Write();
  h1_nvertex_PUW->Write();
  h1_nvertex_PUW_ave->Write();

  h1_rhoPF_presel->Write();
  h1_rhoPF->Write();
  
  h1_pfMet->Write();
  h1_metSignificance->Write();
  h1_metSignificance_2btag->Write();
  h1_mEtSig->Write();
  h1_mEtSig_2btag->Write();
  h1_pfMet_2btag->Write();
  h1_pfMetOverMZZ->Write();
  h1_pfMetOverMZZ_2btag->Write();

  h1_nCandidates->Write();

  h1_ptJet_all_presel->Write();
  h1_etaJet_all_presel->Write();
  h1_nJets_presel->Write();
  h1_nPairs_presel->Write();

  h1_ptResoJet1_beforeKin->Write();
  h1_ptResoJet2_beforeKin->Write();
  h1_ptResoJet1_afterKin->Write();
  h1_ptResoJet2_afterKin->Write();

  h1_ptZreso_beforeKin->Write();
  h1_ptZreso_afterKin->Write();

  h1_mHreso_beforeKin->Write();
  h1_mHreso_afterKin->Write();

  h1_deltaRll_presel->Write();
  h1_deltaRjj_all_presel->Write();

  h1_ptLept1->Write();
  h1_ptLept2->Write();

  h1_ptLept1_presel->Write();
  h1_ptLept2_presel->Write();
  h1_etaLept1_presel->Write();
  h1_etaLept2_presel->Write();

  h1_mZll->Write();
  h1_mZll_presel->Write();
  h1_mZmumu->Write();
  h1_mZmumu_presel->Write();
  h1_mZmumu_presel_0jets->Write();
  h1_mZee->Write();
  h1_mZee_presel->Write();
  h1_mZee_presel_0jets->Write();

  h1_mZjjMC->Write();

  h1_mZjj->Write();
  h1_mZjj_0btag->Write();
  h1_mZjj_1btag->Write();
  h1_mZjj_2btag->Write();
  h1_mZjj_loMass->Write();
  h1_mZjj_medMass->Write();
  h1_mZjj_hiMass->Write();
  h1_mZjj_MU->Write();
  h1_mZjj_ELE->Write();
  h1_mZjj_nogluetag->Write();
  h1_mZjj_nogluetag_MU->Write();
  h1_mZjj_nogluetag_ELE->Write();
  h1_mZjj_loChiSquareProb->Write();
  h1_mZjj_hiChiSquareProb->Write();
  h1_mZjj_all_presel->Write();

  h1_ptZll_presel->Write();
  h1_ptZjj_all_presel->Write();

  h1_deltaRjj_sidebands->Write();
  h1_deltaRjj_prekin_sidebands->Write();
  h1_deltaRjj->Write();
  h1_deltaRjj_prekin->Write();

  h1_ptZjj->Write();
  h1_ptZll->Write();

  h1_cosThetaStar->Write();
  h1_cosTheta1->Write();
  h1_cosTheta2->Write();
  h1_phi->Write();
  h1_phi1->Write();

  h1_cosThetaStar_kinfit->Write();
  h1_cosTheta1_kinfit->Write();
  h1_cosTheta2_kinfit->Write();
  h1_phi_kinfit->Write();
  h1_phi1_kinfit->Write();

  h1_kinfit_chiSquare->Write();
  h1_kinfit_chiSquareProb->Write();
  
  h1_helicityLD->Write();
  h1_helicityLD_nogluetag->Write();
  h1_helicityLD_nokinfit->Write();
  h1_helicityLD_sidebands->Write();
  h1_helicityLD_MW200->Write();
  h1_helicityLD_MW250->Write();
  h1_helicityLD_MW300->Write();
  h1_helicityLD_MW400->Write();
  h1_helicityLD_MW500->Write();
  h1_helicityLD_kinfit->Write();

  h2_helicityLD_vs_mZZ->Write();

  h1_mZZ_hiChiSquareProb->Write();
  h1_mZZ_loChiSquareProb->Write();
  h1_mZZ_mZjj_cut->Write();
  h1_mZZ_mZjj_notcut->Write();
  h1_mZZ_ZjjMassConstr_hiMass->Write();
  h1_mZZ_nokinfit_hiMass_all->Write();
  h1_mZZ_kinfit_hiMass_all->Write();
  h1_mZZ_kinfit_hiMass_hiPU->Write();
  h1_mZZ_kinfit_hiMass_loPU->Write();
  h1_mZZ_kinfit_hiMass_nogluetag->Write();
  h1_mZZ_kinfit_hiMass_gluetag->Write();
  h1_mZZ_kinfit_hiMass_gluetag_ELE->Write();
  h1_mZZ_kinfit_hiMass_gluetag_MU->Write();
  h1_mZZ_kinfit_hiMass_0btag->Write();
  h1_mZZ_kinfit_hiMass_0btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_0btag_MU->Write();
//h1_mZZ_nokinfit_hiMass_0btag->Write();
//h1_mZZ_nokinfit_hiMass_0btag_ELE->Write();
//h1_mZZ_nokinfit_hiMass_0btag_MU->Write();
  h1_mZZ_kinfit_hiMass_1btag->Write();
  h1_mZZ_kinfit_hiMass_1btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_1btag_MU->Write();
//h1_mZZ_nokinfit_hiMass_1btag->Write();
//h1_mZZ_nokinfit_hiMass_1btag_ELE->Write();
//h1_mZZ_nokinfit_hiMass_1btag_MU->Write();
  h1_mZZ_kinfit_hiMass_2btag->Write();
  h1_mZZ_kinfit_hiMass_2btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_2btag_MU->Write();
//h1_mZZ_nokinfit_hiMass_2btag->Write();
//h1_mZZ_nokinfit_hiMass_2btag_ELE->Write();
//h1_mZZ_nokinfit_hiMass_2btag_MU->Write();
  h1_mZZ_kinfit_hiMass_hiQG->Write();
  h1_mZZ_kinfit_hiMass_loQG->Write();

  h1_mZZ_kinfit_hiMass_sidebands_gluetag->Write();
  h1_mZZ_kinfit_hiMass_sidebands_gluetag_ELE->Write();
  h1_mZZ_kinfit_hiMass_sidebands_gluetag_MU->Write();
  h1_mZZ_kinfit_hiMass_sidebands_0btag->Write();
  h1_mZZ_kinfit_hiMass_sidebands_0btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_sidebands_0btag_MU->Write();
  h1_mZZ_kinfit_hiMass_sidebands_1btag->Write();
  h1_mZZ_kinfit_hiMass_sidebands_1btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_sidebands_1btag_MU->Write();
  h1_mZZ_kinfit_hiMass_sidebands_2btag->Write();
  h1_mZZ_kinfit_hiMass_sidebands_2btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_sidebands_2btag_MU->Write();

  h1_deltaRZmatching->Write();
  h1_mZZ_kinfit_hiMass_0btag_matched->Write();

  h1_ptZZ->Write();
  h1_ptZZ_kinfit->Write();
  h1_etaZZ->Write();
  h1_etaZZ_kinfit->Write();

  h1_ptJetRecoil->Write();
  h1_ptHiggs->Write();

  h1_deltaR_part1->Write();
  h1_ptJet1->Write();
  h1_ptJet1_prekin->Write();
  h1_etaJet1->Write();
  h1_eElectronsJet1->Write();
  h1_eMuonsJet1->Write();
  h1_partFlavorJet1->Write();
  h1_partFlavorJet1_MW400->Write();
  h1_partFlavorJet1_MW500->Write();
  h1_partFlavorJet1_matched->Write();
  h1_partFlavorJet1_notmatched->Write();
  h1_partFlavorJet1_0btag->Write();
  h1_partFlavorJet1_1btag->Write();
  h1_partFlavorJet1_2btag->Write();
  h1_partFlavorJet1_gluetag->Write();
  h1_QGLikelihoodJet1_gluonMatched->Write();
  h1_QGLikelihoodJet1_quarkMatched->Write();
  h1_tcheJet1->Write();

  h1_tcheJet->Write();

  h1_QGLikelihoodJetRecoil->Write();
  h1_QGLikelihoodJetRecoil_MW400->Write();
  h1_QGLikelihoodProdRecoil->Write();
  h1_QGLikelihoodProdRecoil_MW400->Write();

  h1_deltaR_part2->Write();
  h1_ptJet2->Write();
  h1_ptJet2_prekin->Write();
  h1_etaJet2->Write();
  h1_eElectronsJet2->Write();
  h1_eMuonsJet2->Write();
  h1_partFlavorJet2->Write();
  h1_partFlavorJet2_MW400->Write();
  h1_partFlavorJet2_MW500->Write();
  h1_partFlavorJet2_matched->Write();
  h1_partFlavorJet2_notmatched->Write();
  h1_partFlavorJet2_0btag->Write();
  h1_partFlavorJet2_1btag->Write();
  h1_partFlavorJet2_2btag->Write();
  h1_partFlavorJet2_gluetag->Write();
  h1_QGLikelihoodJet2_gluonMatched->Write();
  h1_QGLikelihoodJet2_quarkMatched->Write();
  h1_tcheJet2->Write();

  h1_QGLikelihoodProd_oneGluon->Write();
  h1_mZZ_kinfit_hiMass_oneGluon->Write();
  h1_QGLikelihoodProd_twoGluon->Write();
  h1_mZZ_kinfit_hiMass_twoGluon->Write();
  h1_QGLikelihoodProd_twoQuark->Write();
  h1_QGLikelihoodProd_twoQuark_MW400->Write();
  h1_QGLikelihoodProd_twoQuark_MW500->Write();
  h1_mZZ_kinfit_hiMass_twoQuark->Write();

  h1_nEvents_partFlavor_beforeQG->Write();
  h1_nEvents_partFlavor_afterQG->Write();

  h1_nEventsMW400_partFlavor_beforeQG->Write();
  h1_nEventsMW400_partFlavor_afterQG->Write();

  h1_nEventsMW500_partFlavor_beforeQG->Write();
  h1_nEventsMW500_partFlavor_afterQG->Write();


  h1_deltaRZZ->Write();

//h1_mZZ_MCassoc->Write();
//h1_mZZ_MCassoc_ZjjMassConstr->Write();
//h1_mZZ_MCassoc_kinfit->Write();
//h1_mZZ_MCassoc_kinfit_cands->Write();

  h2_mZjj_vs_mZZ->Write();
  h2_mZjj_vs_mZZ_0btag->Write();
  h2_mZjj_vs_mZZ_1btag->Write();
  h2_mZjj_vs_mZZ_2btag->Write();
  h2_mZjj_vs_mZZ_gluetag->Write();

  h2_mZjj_vs_mZZ_kinfit->Write();
  h2_mZjj_vs_mZZ_kinfit_0btag->Write();
  h2_mZjj_vs_mZZ_kinfit_1btag->Write();
  h2_mZjj_vs_mZZ_kinfit_2btag->Write();
  h2_mZjj_vs_mZZ_kinfit_gluetag->Write();


  h1_deltaE_ch->Write();
  h1_deltaE_gamma->Write();
  h1_deltaE_nh->Write();
  h1_deltaEta_ch->Write();
  h1_deltaEta_gamma->Write();
  h1_deltaEta_nh->Write();
  h1_deltaPhi_ch->Write();
  h1_deltaPhi_gamma->Write();
  h1_deltaPhi_nh->Write();
  h1_deltaPt_ch->Write();
  h1_deltaPt_gamma->Write();
  h1_deltaPt_nh->Write();

  h1_ptDJet1->Write();
  h1_ptDJet2->Write();

  h1_nChargedJet1->Write();
  h1_nNeutralJet1->Write();
  h1_ptDJet1->Write();

  h1_nChargedJet2->Write();
  h1_nNeutralJet2->Write();
  h1_ptDJet2->Write();

  h1_QGLikelihoodJet1->Write();
  h1_QGLikelihoodJet2->Write();
  h1_QGLikelihoodProd->Write();

  h1_QGLikelihoodNoPUJet1->Write();
  h1_QGLikelihoodNoPUJet2->Write();
  h1_QGLikelihoodNoPUProd->Write();

  h1_QGLikelihoodJet1_MW300->Write();
  h1_QGLikelihoodJet2_MW300->Write();
  h1_QGLikelihoodProd_MW300->Write();

  h1_QGLikelihoodJet1_MW400->Write();
  h1_QGLikelihoodJet2_MW400->Write();
  h1_QGLikelihoodProd_MW400->Write();

  h1_QGLikelihoodJet1_MW500->Write();
  h1_QGLikelihoodJet2_MW500->Write();
  h1_QGLikelihoodProd_MW500->Write();


  h1_QGLikelihood_100_123->Write();
  h1_QGLikelihood_66_81->Write();
//outFile_->mkdir("QGbins");
//outFile_->cd("QGbins");

//for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

//  vh1_rmsCandJet1[iPtBin]->Write();
//  vh1_ptDJet1[iPtBin]->Write();
//  vh1_nChargedJet1[iPtBin]->Write();
//  vh1_nNeutralJet1[iPtBin]->Write();
//  vh1_QGLikelihoodJet1[iPtBin]->Write();

//  vh1_rmsCandJet2[iPtBin]->Write();
//  vh1_ptDJet2[iPtBin]->Write();
//  vh1_nChargedJet2[iPtBin]->Write();
//  vh1_nNeutralJet2[iPtBin]->Write();
//  vh1_QGLikelihoodJet2[iPtBin]->Write();

//}


  outFile_->Close();


} // finalize()



void Ntp1Finalizer_HZZlljjRM::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;

  if( selectionType_=="presel" ) {

    ptLept1_thresh_ = 10.;
    ptLept2_thresh_ = 10.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.;
    helicityLD_intercept_1btags_ = 0.;
    helicityLD_intercept_2btags_ = 0.;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 9999999999.;
    use_looseBTags_ = true;

  } else if( selectionType_=="presel_invMZll" ) {

    ptLept1_thresh_ = 10.;
    ptLept2_thresh_ = 10.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = true;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.;
    helicityLD_intercept_1btags_ = 0.;
    helicityLD_intercept_2btags_ = 0.;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 9999999999.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_looseBTags_v1" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00124;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.1433;
    helicityLD_intercept_1btags_ = 0.55;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 9999999999.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_looseBTags_v2" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00025;
    helicityLD_slope_1btags_ = 0.000656;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.55;
    helicityLD_intercept_1btags_ = 0.302;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_noBTagCat" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00025;
    helicityLD_slope_1btags_ = 0.00025;;
    helicityLD_slope_2btags_ = 0.00025;
    helicityLD_intercept_0btags_ = 0.55;
    helicityLD_intercept_1btags_ = 0.55;
    helicityLD_intercept_2btags_ = 0.55;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_mediumBTags_v1" ) { //"option 6"

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00054;
    helicityLD_slope_1btags_ = 0.00098;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.428;
    helicityLD_intercept_1btags_ = 0.156;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.5;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = false;

  } else if( selectionType_=="optLD_looseBTags_noQG" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00124;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.1433;
    helicityLD_intercept_1btags_ = 0.55;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_looseBTags_v2_noQG" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00025;
    helicityLD_slope_1btags_ = 0.000656;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.55;
    helicityLD_intercept_1btags_ = 0.302;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_looseBTags_fix400OLD" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    //helicityLD_intercept_0btags_ = 0.639;
    helicityLD_intercept_0btags_ = 0.65;
    helicityLD_intercept_1btags_ = 0.55;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="noCutLD" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.;
    helicityLD_intercept_1btags_ = 0.;
    helicityLD_intercept_2btags_ = 0.;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType



int getNJets( int nPairs ) {

  // fast and easy
  if( nPairs==0 ) return 0;


  int nJets=1;
  while(nJets<50) {
    int trial = nJets*(nJets-1)/2;
    if(trial==nPairs) return nJets;
    nJets++;
  }

  return 0;

}

//else if( nPairs==1 ) return 2;
//else if( nPairs==3 ) return 3;
//else if( nPairs==6 ) return 4;
//else if( nPairs==10 ) return 5;
//else if( nPairs==15 ) return 6;
//else if( nPairs==21 ) return 7;
//else if( nPairs==28 ) return 8;
//else if( nPairs==28 ) return 8;
//else if( nPairs>21 ) 
//  std::cout << "Careful!! Impossible number of pairs!!" << std::endl;

//return 0;

//int nJets = 0;

//int i=1;
//bool found=false;
//while( !found) {
//  nJets += i;
//  i++;
//  if( nJets== nPairs ) found=true;
//}

//return i;

//}


std::vector<TH1D*> getHistoVector(int nPtBins, Double_t *ptBins, std::string histoName, int nBins, float xMin, float xMax ) {

  std::vector<TH1D*> returnVector;

  for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

    Double_t ptMin = ptBins[iPtBin];
    Double_t ptMax = ptBins[iPtBin+1];

    char histoName_thisBin[200];
    sprintf( histoName_thisBin, "%s_pt_%.0lf_%.0lf", histoName.c_str(), ptMin, ptMax);

    TH1D* newHisto = new TH1D(histoName_thisBin, "", nBins, xMin, xMax);
    newHisto->Sumw2();

    returnVector.push_back(newHisto);

  } //for pt bins

  return returnVector;

}



/*
int Ntp1Finalizer_HZZlljjRM::get_nBTags( const AnalysisJet& jet1, const AnalysisJet& jet2, BTagSFUtil* btsfutil, bool loosebtags ) {

  int nBTags;

  bool jet1_tagged_medium = jet1.btag_medium();
  bool jet1_tagged_loose  = jet1.btag_loose();
  bool jet2_tagged_medium = jet2.btag_medium();
  bool jet2_tagged_loose  = jet2.btag_loose();

  btsfutil->modifyBTagsWithSF( jet1_tagged_loose, jet1_tagged_medium, jet1.Pt(), jet1.Eta(), jet1.pdgIdPart );
  btsfutil->modifyBTagsWithSF( jet2_tagged_loose, jet2_tagged_medium, jet2.Pt(), jet2.Eta(), jet2.pdgIdPart );

  if( loosebtags ) {

    bool twoBTags  = ( jet1_tagged_medium && jet2_tagged_loose  )
                  || ( jet1_tagged_loose  && jet2_tagged_medium );
    bool oneBTag   = (!twoBTags) && ( jet1_tagged_loose || jet2_tagged_loose );

    if( twoBTags ) nBTags=2;
    else if( oneBTag ) nBTags=1;
    else nBTags=0;

  } else {

    bool twoBTags  = ( jet1_tagged_medium && jet2_tagged_medium );
    bool oneBTag   = (!twoBTags) && ( jet1_tagged_medium || jet2_tagged_medium );

    if( twoBTags ) nBTags=2;
    else if( oneBTag ) nBTags=1;
    else nBTags=0;

  }

  return nBTags;

}
*/


float Ntp1Finalizer_HZZlljjRM::get_helicityLD_thresh(float mass, int nBTags) {

  float helicityLD_thresh;

  if( nBTags==0 ) {

    helicityLD_thresh = helicityLD_slope_0btags_*mass + helicityLD_intercept_0btags_;
    if( helicityLD_thresh<helicityLD_minThresh_0btags_ ) helicityLD_thresh=helicityLD_minThresh_0btags_;
    if( helicityLD_thresh>helicityLD_maxThresh_0btags_ ) helicityLD_thresh=helicityLD_maxThresh_0btags_;

  } else if( nBTags==1 ) {

    helicityLD_thresh = helicityLD_slope_1btags_*mass + helicityLD_intercept_1btags_;
    if( helicityLD_thresh<helicityLD_minThresh_2btags_ ) helicityLD_thresh=helicityLD_minThresh_2btags_;
    if( helicityLD_thresh>helicityLD_maxThresh_2btags_ ) helicityLD_thresh=helicityLD_maxThresh_2btags_;

  } else if( nBTags==2 ) {

    helicityLD_thresh = helicityLD_slope_2btags_*mass + helicityLD_intercept_2btags_;
    if( helicityLD_thresh<helicityLD_minThresh_2btags_ ) helicityLD_thresh=helicityLD_minThresh_2btags_;
    if( helicityLD_thresh>helicityLD_maxThresh_2btags_ ) helicityLD_thresh=helicityLD_maxThresh_2btags_;

  } else {

    std::cout << "Unexpected number of btags (" << nBTags << "). Returning 0." << std::endl;
    helicityLD_thresh = 0.;
  
  }

  return helicityLD_thresh;

}








float getMuonHLTSF_DoubleTrigger( float pt, float eta, const std::string& runPeriod ) {

  float hltsf = 0.;

  if( runPeriod=="Run2011A" ) {

    if( fabs(eta)<0.8 ) 
      hltsf = 0.975;
    else if( fabs(eta)<2.1 )
      hltsf = 0.955;
    else 
      hltsf = 0.910;

  } else if( runPeriod=="Run2011B" ) {

    if( fabs(eta)<0.8 ) 
      hltsf = 0.972;
    else if( fabs(eta)<2.1 )
      hltsf = 0.945;
    else { // eta 2.1 -> 2.4
      if( pt<40. ) {
        hltsf = 0.85;
      } else {
        hltsf = 0.87;
      }
    }

  } else {

    std::cout << "WARNING! Unknown run period: " << runPeriod << "! Returning HLTSF=0." << std::endl;

  }

  return hltsf;

}




float getMuonHLTSF_SingleTrigger( float pt, float eta, const std::string& runPeriod ) {

  if( pt<25. ) return 0.;

  float hltsf = 0.;

  if( runPeriod=="Run2011A1" ) { //up to may10 technical stop

    if( fabs(eta)<0.8 ) 
      hltsf = 0.896;
    else if( fabs(eta)<2.1 )
      hltsf = 0.807;
    else 
      hltsf = 0.608;

  } else if( runPeriod=="Run2011A2" ) { //from may10 to EPS

    if( fabs(eta)<0.8 ) 
      hltsf = 0.895;
    else if( fabs(eta)<2.1 )
      hltsf = 0.838;
    else 
      hltsf = 0.738;

  } else if( runPeriod=="Run2011A3" ) { //from EPS to end of Run2011A

    if( fabs(eta)<0.8 ) 
      hltsf = 0.890;
    else if( fabs(eta)<2.1 )
      hltsf = 0.809;
    else 
      hltsf = 0.493;

  } else if( runPeriod=="Run2011B" ) {

    if( fabs(eta)<0.8 ) 
      hltsf = 0.87;
    else if( fabs(eta)<2.1 )
      hltsf = 0.79;
    else  //using HLT_IsoMu24_eta2p1
      hltsf = 0.;

  } else {

    std::cout << "WARNING! Unknown run period: " << runPeriod << "! Returning HLTSF=0." << std::endl;

  }

  return hltsf;

}


float getEventHLTSF( float effSingle1, float effSingle2, float effDouble1, float effDouble2 ) {

  float HLTSF = effDouble1 * effDouble2 +
                effSingle2 * (1. - effDouble2 ) +
                effSingle1 * (1. - effDouble1 );

//float HLTSF = effDouble1 * effDouble2 +
//              effSingle2 * (1. - effDouble1*effDouble2 ) +
//              effSingle1 * (1. - effDouble1*effDouble2 ) * (1. - effSingle2);

  return HLTSF;

}
