#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "TFile.h"
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

using namespace std;
using namespace RooFit;

void AllInOne(
	      bool isMuon=false,
	      int btag=0,
	      double timesX=1){


  gSystem->Load("libRooFit");
  gROOT->ProcessLine(".L PDFs/RooRelBW.cc+");//gSystem->Load("PDFs/RooRelBW_cc.so");
  gROOT->ProcessLine(".L PDFs/RooFermi.cc+");//gSystem->Load("PDFs/RooFermi_cc.so");
  gROOT->ProcessLine(".L PDFs/RooDoubleCB.cc+");//gSystem->Load("PDFs/RooDoubleCB_cc.so");
  gROOT->ProcessLine(".L PDFs/RooCB.cc+");//gSystem->Load("PDFs/RooCB_cc.so");
  gSystem->Load("libFFTW");
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  string histoName[3];
  histoName[0]=isMuon?"mZZ_kinfit_hiMass_0btag_MU":"mZZ_kinfit_hiMass_0btag_ELE";
  histoName[1]=isMuon?"mZZ_kinfit_hiMass_1btag_MU":"mZZ_kinfit_hiMass_1btag_ELE";
  histoName[2]=isMuon?"mZZ_kinfit_hiMass_2btag_MU":"mZZ_kinfit_hiMass_2btag_ELE";
  string btagName[3]={"0b","1b","2b"};

  double LumiScale=isMuon?1615.:1556.;

  RooDataSet *data_bkg;  
  RooDataSet *data_temp;
  TFile *file;

  string cutString[3];
  cutString[0]=isMuon?"nBTags==0 && (mZjj>75 && mZjj<105) && leptType==0 && mZZ>183":"nBTags==0 && (mZjj>75 && mZjj<105) && leptType==1 && mZZ>183";
  cutString[1]=isMuon?"nBTags==1 && (mZjj>75 && mZjj<105) && leptType==0 && mZZ>183":"nBTags==1 && (mZjj>75 && mZjj<105) && leptType==1 && mZZ>183";
  cutString[2]=isMuon?"nBTags==2 && (mZjj>75 && mZjj<105) && leptType==0 && mZZ>183":"nBTags==2 && (mZjj>75 && mZjj<105) && leptType==1 && mZZ>183";

  int binWidth=20;
  int highBin=750;
  int lowBin=150;

  double muonEff[3]={.576,.548,.489};

  double expSig[3];
  expSig[0]=isMuon ? muonEff[0]*5.65*LumiScale : (1-muonEff[0])*5.65*LumiScale ;
  expSig[1]=isMuon ? muonEff[1]*4.89*LumiScale : (1-muonEff[1])*4.89*LumiScale ;
  expSig[2]=isMuon ? muonEff[2]*2.37*LumiScale : (1-muonEff[2])*2.37*LumiScale ;
  double expBkg[3];
  expBkg[0]=isMuon ? 575.85: 490.54;
  expBkg[1]=isMuon ? 606.78: 526.86;
  expBkg[2]=isMuon ? 41.15 : 35.37 ;

  stringstream convert;

  // --------------------- measurable (ZZ invariant mass) ----------------
  string temp;
  if(isMuon)
    temp="m_{#mu#mujj}";  else
    temp="m_{eejj}";
  RooRealVar mZZ("mZZ",temp.c_str(),lowBin,highBin);
  RooRealVar nBTags("nBTags","nBTags",-1.,3.);
  RooRealVar leptType("leptType","leptType",0,1);
  RooRealVar mZjj("mZjj","mZjj",0,200.);

  // ----------------- get parameters from data cards! -----------------

  std::string leptTag = (isMuon) ? "mm" : "ee";
  char wsFileName[400];
  sprintf( wsFileName, "hzz2l2q_%s%db.input.root", leptTag.c_str(), btag );
  
  TFile *wsFile = new TFile(wsFileName);
  RooWorkspace *ws = (RooWorkspace*) wsFile->Get("w");

  // ==================== defining bkg PDF ==========================
  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",ws->var("cutOff_BKG")->getVal());
  cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",ws->var("beta_BKG")->getVal());
  beta.setConstant(kTRUE);
	     		       
  RooFermi fermi("fermi","fermi function",mZZ,cutOff,beta);
  // -------------------- crystal ball ---------------------------
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

  // ------------------ signal PDF ------------------------
  // -------------------- fermi ------------------------
  
  RooRealVar cutOffSig("cutOffSig","cutOff",ws->var("cutOff_SIG")->getVal());
  RooRealVar gSig("gSig","g",ws->var("g_SIG")->getVal());

  RooFermi fermi2("fermi2","fermi2",mZZ,cutOffSig,gSig);

  RooRealVar cutOffSig2("cutOffSig2","cutOff2",ws->var("cutOff2_SIG")->getVal());
  RooRealVar gSig2("gSig2","g2",ws->var("g2_SIG")->getVal());

  RooFermi fermi3("fermi2","fermi2",mZZ,cutOffSig2,gSig2);

  // ------------------- Relativistic BW --------------------------------
  //
 
  RooRealVar BW_mean("BW_mean", "mean",ws->var("BW_mean")->getVal());
  BW_mean.setConstant(kTRUE);
  RooRealVar BW_sigma("BW_sigma", "sigma",ws->var("BW_sigma")->getVal());
  BW_sigma.setConstant(kTRUE);
  RooRealVar BW_n("BW_n","n",0.,0.,1.);

  RooRelBW BW("BW","Relativistic B-W",mZZ,BW_mean,BW_sigma,BW_n);

  // ------------------- Crystal Ball -------------------------------
  temp="CMS_hzz2l2q_sig"+btagName[btag]+"p1";
  RooRealVar CB_mean("CB_mean","param 1 of CB",ws->var(temp.c_str())->getVal());
  temp="CMS_hzz2l2q_sig"+btagName[btag]+"p2";
  RooRealVar CB_sigma("CB_sigma","param 2 of CB",ws->var(temp.c_str())->getVal());
  temp="CB_alpha1";
  RooRealVar CB_alpha1("CB_alpha1","param 3 of CB",ws->var(temp.c_str())->getVal());
  temp="CB_n1";
  RooRealVar CB_n1("CB_n1","param 4 of CB",ws->var(temp.c_str())->getVal());
  temp="CB_alpha2";
  RooRealVar CB_alpha2("CB_alpha2","param 4 of CB",ws->var(temp.c_str())->getVal());
  temp="CB_n2";
  RooRealVar CB_n2("CB_n2","param 5 of CB",ws->var(temp.c_str())->getVal());

  RooDoubleCB CBSig("CBSig","Crystal Ball",mZZ,CB_mean,CB_sigma,CB_alpha1,CB_n1,CB_alpha2,CB_n2);
  
  //------------------------ convolution -------------------------
  //Set #bins to be used for FFT sampling to ...
  mZZ.setBins(10000,"fft");

  //RooFFTConvPdf signal("signal","Rel-BW (X) CB",mZZ,BW,CB);
  //signal.setBufferFraction(1.0);
  RooFFTConvPdf sig("sig","Rel-BW (X) CB",mZZ,BW,CBSig);
  sig.setBufferFraction(1.0);
  
  RooProdPdf signal("signal","signal",RooArgSet(fermi2,fermi3,sig));

  // ------------------ bkg+sig PDF ----------------------
  
  RooRealVar nSig("nSig","nSig",expSig[btag]*10/(expSig[btag]*10+expBkg[btag]));

  RooAddPdf model("model","model",signal,background,nSig);

  // ------------------ get data --------------------------
  // for reading sideband extrapolated data...
  file = new TFile("HZZlljjRM_DATA_LP11_optLD_looseBTags_v2_ALL.root");
  TTree* t=(TTree*)file->Get("tree_passedEvents");
  data_bkg=new RooDataSet("data_bkg","data_bkg",t,
			  RooArgSet(mZZ,leptType,nBTags,mZjj),
			  cutString[btag].c_str());

  // --------- draw MC data -------------------
  RooPlot *plot_MCbkg = mZZ.frame(lowBin,highBin,(int)(highBin-lowBin)/binWidth);

  //-----------------------------------------------------------------------
  TCanvas *c2 = new TCanvas("c2","c2",600,600);

  TPaveText* cmslabel = new TPaveText( 0.147, 0.953, 0.64, 0.975, "brNDC" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextAlign(11);
  cmslabel->SetTextFont(62);
  cmslabel->SetBorderSize(0);

  if(isMuon)
    cmslabel->AddText("CMS Preliminary 2011, 1.615 fb^{-1}");
  else
    cmslabel->AddText("CMS Preliminary 2011, 1.556 fb^{-1}");

  TPaveText* label_sqrt = new TPaveText(0.7,0.953,.96,0.975, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetTextAlign(31); // align right
  label_sqrt->SetBorderSize(0);
  label_sqrt->AddText("#sqrt{s} = 7 TeV");
  //-----------------------------------------------------------------------

  background.plotOn(plot_MCbkg,Normalization(expBkg[btag]));
  //model.plotOn(plot_MCbkg,Normalization(expSig[btag]*10+expBkg[btag]),LineStyle(2));
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
  if(isMuon)
    temp=";m_{#mu#mujj} [GeV]; Events / "+convert.str()+" GeV";
  else
    temp=";m_{eejj} [GeV]; Events / "+convert.str()+" GeV";
  hTT->SetTitle(temp.c_str());
  hTT->Scale(LumiScale);
  hTT->Rebin(binWidth);
  hTT->SetFillColor(39);
  gStyle->SetOptStat(0);

  if(btag==0)
    hTT->GetYaxis()->SetRangeUser(0.0001,isMuon?130:120);
  if(btag==1)
    hTT->GetYaxis()->SetRangeUser(0.0001,isMuon?150:130);
  if(btag==2)
    hTT->GetYaxis()->SetRangeUser(0.0001,isMuon?25:25);

  hTT->GetXaxis()->SetRangeUser(183., 800.);
  hTT->Draw();

  TH1F *hVV =(TH1F*)VVFile->Get(histoName[btag].c_str());
  hVV->SetName("hVV");
  hVV->Scale(LumiScale);
  hVV->Rebin(binWidth);
  hVV->SetFillColor(38);

  TH1F* hH400 = (TH1F*)H400File->Get(histoName[btag].c_str());
  hH400->SetName("hH400");
  hH400->Scale(LumiScale*timesX);
  hH400->Rebin(binWidth);
  hH400->SetFillColor(kYellow);//kRed+3);

  temp = ";m_{lljj} [GeV];Events / "+convert.str()+" GeV";
  THStack *hBkg = new THStack("hBkg",temp.c_str());
  convert.str("");

  hBkg->Add(hVV);
  hBkg->Add(hTT);
  hBkg->Add(hZjets);
  hBkg->Add(hH400);

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
  convert << (double)LumiScale/1000.;
  if(isMuon)
    temp="CMS Preliminary #sqrt{s}=7 TeV "+convert.str()+" fb^{-1}";
  else
    temp="CMS Preliminary #sqrt{s}=7 TeV "+convert.str()+" fb^{-1}";
  //leg->SetHeader(temp.c_str());
  leg->AddEntry("background_Norm[mZZ]","Sideband Extrapolated Fit","l");
  //leg->AddEntry("model_Norm[mZZ]","Background+10#timesSignal","l");

  if(btag==0)
    temp=isMuon?"0 b-tag 2#mu2q data":"0 b-tag 2e2q data";
  if(btag==1)
    temp=isMuon?"1 b-tag 2#mu2q data":"1 b-tag 2e2q data";
  if(btag==2)
    temp=isMuon?"2 b-tag 2#mu2q data":"2 b-tag 2e2q data";

  leg->AddEntry("h_data_bkg",temp.c_str(),"p");
  leg->AddEntry("hZjets","Z + Jets","f");
  leg->AddEntry("hTT","tt/tW","f");
  leg->AddEntry("hVV","ZZ/WZ/WW","f");
  convert.str("");
  convert << timesX;
  temp="400 GeV SM Higgs #times"+convert.str();
  leg->AddEntry("hH400",temp.c_str(),"f");
  leg->Draw();

  string saveFileName="AllInOne_"+btagName[btag]+"tag"+(isMuon?"_mm.eps":"_ee.eps");

  c2->SaveAs(saveFileName.c_str()); 

}

