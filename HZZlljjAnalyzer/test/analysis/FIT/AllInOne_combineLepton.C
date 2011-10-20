#include "TH2D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "THStack.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooPolynomial.h"
#include "RooFFTConvPdf.h"
#include "TStyle.h"
#include "RooAddPdf.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "PDFs/RooCB.h"
#include "PDFs/RooRelBW.h"
#include "PDFs/RooFermi.h"
#include "PDFs/RooDoubleCB.h"

using namespace std;
using namespace RooFit;

void AllInOne_combineLepton(const std::string& dataset,int btag=0, double timesX=1, bool useNewShape=false){

  std::string useNewShapeText = (useNewShape) ? "_newshape" : "";

  char wsFileName[700];
  sprintf(wsFileName,"datacards/400/hzz2l2q_ee%db.input.root",btag);
  //sprintf(wsFileName,"PROVA/datacards_%s%s/hzz2l2q_ee%db.input.root",dataset.c_str(), useNewShapeText.c_str(),btag);

  gSystem->Load("libRooFit");
  gSystem->Load("libFFTW");

  string histoName[3];
  histoName[0]="mZZ_kinfit_hiMass_0btag";
  histoName[1]="mZZ_kinfit_hiMass_1btag";
  histoName[2]="mZZ_kinfit_hiMass_2btag";
  string btagName[3]={"0b","1b","2b"};

  double LumiScale=0.;
  if( dataset=="Run2011A_FULL" ) LumiScale = 2100.;
  else if( dataset=="LP11" ) LumiScale = 1600.;
  else {
    std::cout << "Unknown dataset '" << dataset << "'. Exiting." << std::endl;
  }

  RooDataSet *data_bkg;  
  RooDataSet *data_temp;
  TFile *file;

  string cutString[3];
  cutString[0]="nBTags==0 && (mZjj>75 && mZjj<105) && mZZ>183";
  cutString[1]="nBTags==1 && (mZjj>75 && mZjj<105) && mZZ>183";
  cutString[2]="nBTags==2 && (mZjj>75 && mZjj<105) && mZZ>183";

  int binWidth=20;
  int highBin=750;
  int lowBin=150;

//double muonEff[3]={.576,.548,.489};

//double expSig[3];
//expSig[0]=5.65*LumiScale ;
//expSig[1]=4.89*LumiScale ;
//expSig[2]=2.37*LumiScale ;

  char alphaFileName[200];
  sprintf( alphaFileName, "alphaFile_%s_%dbtag_ALL.root", dataset.c_str(), btag);
  TFile* alphaFile = TFile::Open(alphaFileName);
  TTree* treeSidebandsDATA_alphaCorr = (TTree*)alphaFile->Get("sidebandsDATA_alpha");
  TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, 150., 800.);
  char sidebandsCut_alpha[500];
  sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d)", btag);
  treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "mZZ", sidebandsCut_alpha);
  float expBkg = h1_mZZ_sidebands_alpha->Integral();
  std::cout <<  "++++ expBkg: " << expBkg << std::endl;



  stringstream convert;

  // --------------------- measurable (ZZ invariant mass) ----------------
  string temp;
  temp="m_{ZZ}";
  RooRealVar mZZ("mZZ",temp.c_str(),lowBin,highBin);
  RooRealVar nBTags("nBTags","nBTags",-1.,3.);
  RooRealVar leptType("leptType","leptType",0,1);
  RooRealVar mZjj("mZjj","mZjj",0,200.);

  // ----------------- get parameters from data cards! -----------------
  
  TFile *wsFile = new TFile(wsFileName);
  RooWorkspace *ws = (RooWorkspace*) wsFile->Get("w");

  // ==================== defining bkg PDF ==========================
  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",ws->var("cutOff_BKG")->getVal());
  cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",ws->var("beta_BKG")->getVal());
  beta.setConstant(kTRUE);
	     		       
  RooFermi fermi("fermi","fermi function",mZZ,cutOff,beta);
  // -------------------- double gauss ---------------------------
  temp="CMS_hzz2l2q_bkg"+btagName[btag]+"p1";
  RooRealVar m("m","m",ws->var(temp.c_str())->getVal());
  m.setConstant(kTRUE);
  temp="CMS_hzz2l2q_bkg"+btagName[btag]+"p2";
  RooRealVar wdth("wdth","wdth",ws->var(temp.c_str())->getVal());
  wdth.setConstant(kTRUE);
  temp="CMS_hzz2l2q_bkg"+btagName[btag]+"p3";
  RooRealVar n("n","n",ws->var(temp.c_str())->getVal());
  n.setConstant(kTRUE);
  temp="CMS_hzz2l2q_bkg"+btagName[btag]+"p4";
  RooRealVar alpha("alpha","alpha",ws->var(temp.c_str())->getVal());
  alpha.setConstant(kTRUE);
  temp="CMS_hzz2l2q_bkg"+btagName[btag]+"p5";
  RooRealVar theta("theta","theta",ws->var(temp.c_str())->getVal());
  theta.setConstant(kTRUE);

  RooCB CB("CB","Crystal ball",mZZ,m,wdth,alpha,n,theta);
  
  RooProdPdf background("background","background",RooArgSet(fermi,CB));


  // ------------------ get data --------------------------
  // for reading sideband extrapolated data...
  std::string fileName = "HZZlljjRM_DATA_" + dataset + "_optLD_looseBTags_v2_ALL.root";
  file = new TFile(fileName.c_str());
  //file = new TFile("../HZZlljjRM_DATA_Run2011A_FULL_optLD_looseBTags_v2_ALL.root");
  TTree* t=(TTree*)file->Get("tree_passedEvents");
  data_bkg=new RooDataSet("data_bkg","data_bkg",t,
			  RooArgSet(mZZ,leptType,nBTags,mZjj),
			  cutString[btag].c_str());

  // --------- draw MC data -------------------
  RooPlot *plot_MCbkg = mZZ.frame(lowBin,highBin,(int)(highBin-lowBin)/binWidth);

  //-----------------------------------------------------------------------
  TCanvas *c2 = new TCanvas("c2","c2",600,600);

  TPaveText* cmslabel = new TPaveText( 0.145, 0.953, 0.6, 0.975, "brNDC" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextAlign(11);
  cmslabel->SetTextFont(62);
  cmslabel->SetBorderSize(0);
  char lumilabel[500];
  sprintf( lumilabel, "CMS Preliminary 2011, %.1f fb^{-1}", LumiScale/1000.);
  cmslabel->AddText(lumilabel);

  TPaveText* label_sqrt = new TPaveText(0.7,0.953,0.96,0.975, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetTextAlign(31); // align right
  label_sqrt->SetBorderSize(0);
  label_sqrt->AddText("#sqrt{s} = 7 TeV");
  //-----------------------------------------------------------------------

  background.plotOn(plot_MCbkg,Normalization(expBkg));
  data_bkg->plotOn(plot_MCbkg,Binning((int)(highBin-lowBin)/binWidth));
  // -------------------- get histograms -----------------
  
  TFile *ZjetsFile = new TFile("HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  TFile *TTFile = new TFile("HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  TFile *VVFile = new TFile("HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");
  TFile *H400File = new TFile("HZZlljjRM_GluGluToHToZZTo2L2Q_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root");

  TH1F *hZjets =(TH1F*)ZjetsFile->Get(histoName[btag].c_str());
  hZjets->SetName("hZjets");
  hZjets->Scale(LumiScale);
  hZjets->Rebin(binWidth);
  hZjets->SetFillColor(30);

  TH1F *hTT =(TH1F*)TTFile->Get(histoName[btag].c_str());
  hTT->SetName("hTT");
  convert << binWidth;
  temp=";m_{ZZ} [GeV]; Events / "+convert.str()+" GeV";
  hTT->SetTitle(temp.c_str());
  hTT->Scale(LumiScale);
  hTT->Rebin(binWidth);
  hTT->SetFillColor(39);
  gStyle->SetOptStat(0);

  //if(btag==0)
  //  hTT->GetYaxis()->SetRangeUser(0.0001,245);
  //if(btag==1)
  //  hTT->GetYaxis()->SetRangeUser(0.0001,245);
  //if(btag==2)
  //  hTT->GetYaxis()->SetRangeUser(0.0001,25);

  //hTT->Draw();


  TH1F *hVV =(TH1F*)VVFile->Get(histoName[btag].c_str());
  hVV->SetName("hVV");
  hVV->Scale(LumiScale);
  hVV->Rebin(binWidth);
  hVV->SetFillColor(38);

  TH1F* hH400 = (TH1F*)H400File->Get(histoName[btag].c_str());
  hH400->SetName("hH400");
  hH400->Scale(LumiScale*timesX);
  cout << "SIGNAL NORMALIZATION: " << hH400->Integral() << endl;
  hH400->Rebin(binWidth);
  hH400->SetFillColor(kYellow);//kRed+3);

  temp = ";m_{ZZ} [GeV];Events / "+convert.str()+" GeV";
  THStack *hBkg = new THStack("hBkg",temp.c_str());
  convert.str("");

  hBkg->Add(hVV);
  hBkg->Add(hTT);
  hBkg->Add(hZjets);
  hBkg->Add(hH400);

  float xMin = 150.;
  float xMax = 750.;
  float yMin = 0.;
  float yMax = 1.3*hBkg->GetMaximum();

  char yAxisName[500];
  sprintf( yAxisName, "Events / %.0f GeV", hVV->GetXaxis()->GetBinWidth(1) );
  //std::string yAxisName = "Events / "+convert.str()+" GeV";
  convert.str("");

//char xAxisName[400];
//sprintf( xAxisName, "m_{%s%sjj} [GeV]", leptType_forLegend.c_str(), leptType_forLegend.c_str() );

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax);
  h2_axes->SetXTitle("m_{lljj} [GeV]");
  h2_axes->SetYTitle(yAxisName);
  h2_axes->Draw();

  hBkg->Draw("SAMEHIST");
  plot_MCbkg->Draw("SAME");
  cmslabel->Draw();
  label_sqrt->Draw();
  // ---------------legend ---------------------------

  TLegend *leg = new TLegend(.4,.5,.8,.9);
  leg->SetTextSize(0.036);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  convert.str("");
  convert << LumiScale/1000;
  temp="CMS Preliminary #sqrt{s}=7 TeV "+convert.str()+" fb^{-1}";
  //leg->SetHeader(temp.c_str());
  leg->AddEntry("background_Norm[mZZ]","Sideband Extrapolated Fit","l");
  //leg->AddEntry("model_Norm[mZZ]","Background+10#timesSignal","l");

  if(btag==0)
    temp="0 b-tag 2l2q data";
  if(btag==1)
    temp="1 b-tag 2l2q data";
  if(btag==2)
    temp="2 b-tag 2l2q data";

  leg->AddEntry("h_data_bkg",temp.c_str(),"p");
  leg->AddEntry("hZjets","Z + Jets","f");
  leg->AddEntry("hTT","tt/tW","f");
  leg->AddEntry("hVV","ZZ/WZ/WW","f");
  convert.str("");
  convert << timesX;
  temp="400 GeV SM Higgs#times"+convert.str();
  leg->AddEntry("hH400",temp.c_str(),"f");
  leg->Draw();

  string saveFileName="AllInOne_"+dataset + "_"+btagName[btag]+"tag_ll_LIMIT" + useNewShapeText + ".eps";

  c2->SaveAs(saveFileName.c_str()); 

}

