#include "TFile.h"
#include "TH1D.h"



void computeBGyields() {


  TFile* fileDY = TFile::Open("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  TFile* fileVV = TFile::Open("HZZlljjRM_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_optLD_looseBTags_v2_ALL.root");
  TFile* fileTT = TFile::Open("HZZlljjRM_TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_optLD_looseBTags_v2_ALL.root");


  TTree* treeDY = (TTree*)fileDY->Get("tree_passedEvents");
  TTree* treeVV = (TTree*)fileVV->Get("tree_passedEvents");
  TTree* treeTT = (TTree*)fileTT->Get("tree_passedEvents");

  TH1D* h1_mZZ_DY_EE_0btag = new TH1D("mZZ_DY_EE_0btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_VV_EE_0btag = new TH1D("mZZ_VV_EE_0btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_TT_EE_0btag = new TH1D("mZZ_TT_EE_0btag", "", 600, 183., 800.);

  TH1D* h1_mZZ_DY_MM_0btag = new TH1D("mZZ_DY_MM_0btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_VV_MM_0btag = new TH1D("mZZ_VV_MM_0btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_TT_MM_0btag = new TH1D("mZZ_TT_MM_0btag", "", 600, 183., 800.);

  TH1D* h1_mZZ_DY_EE_1btag = new TH1D("mZZ_DY_EE_1btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_VV_EE_1btag = new TH1D("mZZ_VV_EE_1btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_TT_EE_1btag = new TH1D("mZZ_TT_EE_1btag", "", 600, 183., 800.);

  TH1D* h1_mZZ_DY_MM_1btag = new TH1D("mZZ_DY_MM_1btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_VV_MM_1btag = new TH1D("mZZ_VV_MM_1btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_TT_MM_1btag = new TH1D("mZZ_TT_MM_1btag", "", 600, 183., 800.);

  TH1D* h1_mZZ_DY_EE_2btag = new TH1D("mZZ_DY_EE_2btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_VV_EE_2btag = new TH1D("mZZ_VV_EE_2btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_TT_EE_2btag = new TH1D("mZZ_TT_EE_2btag", "", 600, 183., 800.);

  TH1D* h1_mZZ_DY_MM_2btag = new TH1D("mZZ_DY_MM_2btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_VV_MM_2btag = new TH1D("mZZ_VV_MM_2btag", "", 600, 183., 800.);
  TH1D* h1_mZZ_TT_MM_2btag = new TH1D("mZZ_TT_MM_2btag", "", 600, 183., 800.);


  treeDY->Project("mZZ_DY_EE_0btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==0 && !isSidebands && leptType==1)");
  treeVV->Project("mZZ_VV_EE_0btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==0 && !isSidebands && leptType==1)");
  treeTT->Project("mZZ_TT_EE_0btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==0 && !isSidebands && leptType==1)");

  treeDY->Project("mZZ_DY_MM_0btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==0 && !isSidebands && leptType==0)");
  treeVV->Project("mZZ_VV_MM_0btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==0 && !isSidebands && leptType==0)");
  treeTT->Project("mZZ_TT_MM_0btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==0 && !isSidebands && leptType==0)");

  treeDY->Project("mZZ_DY_EE_1btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==1 && !isSidebands && leptType==1)");
  treeVV->Project("mZZ_VV_EE_1btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==1 && !isSidebands && leptType==1)");
  treeTT->Project("mZZ_TT_EE_1btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==1 && !isSidebands && leptType==1)");

  treeDY->Project("mZZ_DY_MM_1btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==1 && !isSidebands && leptType==0)");
  treeVV->Project("mZZ_VV_MM_1btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==1 && !isSidebands && leptType==0)");
  treeTT->Project("mZZ_TT_MM_1btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==1 && !isSidebands && leptType==0)");

  treeDY->Project("mZZ_DY_EE_2btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==2 && !isSidebands && leptType==1)");
  treeVV->Project("mZZ_VV_EE_2btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==2 && !isSidebands && leptType==1)");
  treeTT->Project("mZZ_TT_EE_2btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==2 && !isSidebands && leptType==1)");

  treeDY->Project("mZZ_DY_MM_2btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==2 && !isSidebands && leptType==0)");
  treeVV->Project("mZZ_VV_MM_2btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==2 && !isSidebands && leptType==0)");
  treeTT->Project("mZZ_TT_MM_2btag", "mZZ", "eventWeight*(mZZ>183. && mZZ<800. && nBTags==2 && !isSidebands && leptType==0)");


  float yield_EE_0btag = 1000.*(h1_mZZ_DY_EE_0btag->Integral() + h1_mZZ_VV_EE_0btag->Integral() + h1_mZZ_TT_EE_0btag->Integral());
  float yield_EE_1btag = 1000.*(h1_mZZ_DY_EE_1btag->Integral() + h1_mZZ_VV_EE_1btag->Integral() + h1_mZZ_TT_EE_1btag->Integral());
  float yield_EE_2btag = 1000.*(h1_mZZ_DY_EE_2btag->Integral() + h1_mZZ_VV_EE_2btag->Integral() + h1_mZZ_TT_EE_2btag->Integral());

  float yield_MM_0btag = 1000.*(h1_mZZ_DY_MM_0btag->Integral() + h1_mZZ_VV_MM_0btag->Integral() + h1_mZZ_TT_MM_0btag->Integral());
  float yield_MM_1btag = 1000.*(h1_mZZ_DY_MM_1btag->Integral() + h1_mZZ_VV_MM_1btag->Integral() + h1_mZZ_TT_MM_1btag->Integral());
  float yield_MM_2btag = 1000.*(h1_mZZ_DY_MM_2btag->Integral() + h1_mZZ_VV_MM_2btag->Integral() + h1_mZZ_TT_MM_2btag->Integral());

  float entries_EE_0btag = h1_mZZ_DY_EE_0btag->GetEntries();
  float entries_EE_1btag = h1_mZZ_DY_EE_1btag->GetEntries();
  float entries_EE_2btag = h1_mZZ_DY_EE_2btag->GetEntries();

  float entries_MM_0btag = h1_mZZ_DY_MM_0btag->GetEntries();
  float entries_MM_1btag = h1_mZZ_DY_MM_1btag->GetEntries();
  float entries_MM_2btag = h1_mZZ_DY_MM_2btag->GetEntries();

  float weight_EE_0btag = (yield_EE_0btag/entries_EE_0btag);
  float weight_EE_1btag = (yield_EE_1btag/entries_EE_1btag);
  float weight_EE_2btag = (yield_EE_2btag/entries_EE_2btag);

  float weight_MM_0btag = (yield_MM_0btag/entries_MM_0btag);
  float weight_MM_1btag = (yield_MM_1btag/entries_MM_1btag);
  float weight_MM_2btag = (yield_MM_2btag/entries_MM_2btag);

  float error_EE_0btag = weight_EE_0btag*sqrt(entries_EE_0btag);
  float error_EE_1btag = weight_EE_1btag*sqrt(entries_EE_1btag);
  float error_EE_2btag = weight_EE_2btag*sqrt(entries_EE_2btag);

  float error_MM_0btag = weight_MM_0btag*sqrt(entries_MM_0btag);
  float error_MM_1btag = weight_MM_1btag*sqrt(entries_MM_1btag);
  float error_MM_2btag = weight_MM_2btag*sqrt(entries_MM_2btag);
  

  std::cout << "ELECTRONS:" << std::endl;
  std::cout << yield_EE_0btag << "+-" << error_EE_0btag << "\t" << yield_EE_1btag << "+-" << error_EE_1btag << "\t" << yield_EE_2btag << "+-" << error_EE_2btag << std::endl;
  std::cout << "MUONS: " << std::endl;
  std::cout << yield_MM_0btag << "+-" << error_MM_0btag << "\t" << yield_MM_1btag << "+-" << error_MM_1btag << "\t" << yield_MM_2btag << "+-" << error_MM_2btag << std::endl;

}
