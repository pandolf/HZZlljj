#include "TCanvas.h"
#include "RooRealVar.h"
#include "PDFs/RooCB.h"
#include "RooProdPdf.h"
#include "TFile.h"
#include "TTree.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "PDFs/RooFermi.h"
#include "RooEllipse.h"
#include "RooFitResult.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include "TSystem.h"

using namespace std;
using namespace RooFit;

RooFitResult* fitBkgZZinvMass_CB_DATA(double rot=-0.627836893, int btag=0){

  //gSystem->Load("libRooFit");
  //gROOT->ProcessLine(".L PDFs/RooFermi_cc.so");
  //gROOT->ProcessLine(".L PDFs/RooCB_cc.so");
  gSystem->Load("libFFTW");
  //gROOT->ProcessLine(".L ~/tdrstyle.C");
  //setTDRStyle();


  // --------------------- initial values -----------------
  string plotName;
  vector<string> cutString;
  cutString.push_back("nBTags==0 && ((mZjj>60 && mZjj<75) || (mZjj>105 && mZjj<130))");
  cutString.push_back("nBTags==1 && ((mZjj>60 && mZjj<75) || (mZjj>105 && mZjj<130))");
  cutString.push_back("nBTags==2 && ((mZjj>60 && mZjj<75) || (mZjj>105 && mZjj<130))");
  vector<double> paramVal;
  if(btag==0){//-1.5545 rot for Run2011A
    paramVal.push_back(186.41); //cutOff   
    paramVal.push_back(5.257);	//beta	    
    paramVal.push_back(222.72);	//CB_mean  
    //paramVal.push_back(53.876);	//CB_wdth   54.686 --- 1.0 /fb fits
    paramVal.push_back(27.6372);	//CB_wdth   2.0 /fb fits    
    paramVal.push_back(25.775);	//CB_n	    13.39 
    //paramVal.push_back(-.8312);	//CB_alpha  -.879 
    paramVal.push_back(-0.427402);	//CB_alpha  -.879 
  }
  if(btag==1){//-1.552 rot for Run2011A
    paramVal.push_back(184.04);  //cutOff   
    paramVal.push_back(3.225);	 //beta	    
    paramVal.push_back(166.6);	 //CB_mean  
    //paramVal.push_back(94.7190 );//CB_wdth   87.446    --- 1.0/fb fits
    paramVal.push_back(83.3205);//CB_wdth   87.446    --- 2.0/fb fits
    paramVal.push_back(6.5670  );//CB_n	     21.499  
    //paramVal.push_back(-1.47023);//CB_alpha  -1.16716
    paramVal.push_back(-1.31797);//CB_alpha  -1.16716
  }
  if(btag==2){//-1.5574 rot for Run2011A
    paramVal.push_back(182.32);  //cutOff   
    paramVal.push_back(2.8492);	 //beta	    
    paramVal.push_back(226.23);	 //CB_mean  
    //paramVal.push_back(39.4034); //CB_wdth  72.11  --- 1.0/fb fits
    paramVal.push_back(46.4583); //CB_wdth  72.11  --- 2.0/fb fits
    paramVal.push_back(6.1580 ); //CB_n	    2.7583
    //paramVal.push_back(-.74438); //CB_alpha -1.654
    paramVal.push_back(-0.622234); //CB_alpha -1.654
  }
  cout << "parameters set..."<< endl;
  // --------------------- measurable (ZZ invariant mass) ----------------

  RooRealVar mZZ("mZZ", "zz inv mass", 160.,800.);
  RooRealVar nBTags("nBTags","nBTags",-1.,3.);
  //RooRealVar isSB("isSB","isSB",-1.,1.);
  RooRealVar eventWeight_alpha("eventWeight_alpha","eventWeight_alpha",0,100.);
  RooRealVar mZjj("mZjj","mZjj",0,150.);
  RooRealVar leptType("leptType","leptType",0,2);

  // ==================== defining bkg PDF ==========================
  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",paramVal.at(0),0,1000);
  cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",paramVal.at(1),0,50);
  beta.setConstant(kTRUE);
	     		       
  RooFermi fermi("fermi","fermi function",mZZ,cutOff,beta);
  // -------------------- double gauss ---------------------------

  double a=cos(rot)*paramVal.at(5) + sin(rot)*paramVal.at(3);
  double w=-sin(rot)*paramVal.at(5) + cos(rot)*paramVal.at(3);
  cout << "alpha true: " << a << endl;
  cout << "width true: " << w << endl;
  
  RooRealVar m("m","m",paramVal.at(2),200,1000);
  //m.setConstant(kTRUE);
  RooRealVar wdth("wdth","#beta",w,-1000,1000);
  //wdth.setConstant(kTRUE);
  RooRealVar n("n","n",paramVal.at(4),0,100);
  n.setConstant(kTRUE);
  RooRealVar alpha("alpha","#alpha",a,-1000,1000); 
  //alpha.setConstant(kTRUE);
  RooRealVar theta("theta","theta",rot,-3.1415,3.1415); 
  theta.setConstant(kTRUE);
  
  RooCB CB("CB","Crystal ball",mZZ,m,wdth,alpha,n,theta);
  
  RooProdPdf background("background","background",RooArgSet(fermi,CB));

  // ------------------ get data --------------------------
  // for reading sideband extrapolated data...
  TFile *file = new TFile("HZZlljjRM_DATA_LP11_optLD_looseBTags_v2_ALL.root");
  RooDataSet* data_bkg=new RooDataSet("data_bkg","data_bkg",(TTree*)file->Get("tree_passedEvents"),RooArgSet(mZZ,nBTags,mZjj,eventWeight_alpha,leptType),cutString.at(btag).c_str(),"eventWeight_alpha");

  //RooDataSet* data_bkg=new RooDataSet("data_bkg","data_bkg",(TTree*)file->Get("selectedEvents"),RooArgSet(mZZ,nBTags,mZjj,wght,lep),cutString.at(btag).c_str());//,"wght");

  cout << "check" << endl;

  // --------- draw MC data -------------------

  //RooFitResult *r = background.fitTo(data_bkg,SumW2Error(kTRUE),InitialHesse(kTRUE),Minos(kTRUE),Save());  
  RooFitResult *r = background.fitTo(*data_bkg,SumW2Error(kTRUE),Save());  

  RooPlot *plot_MCbkg = mZZ.frame(160,800,64);
  data_bkg->plotOn(plot_MCbkg,DataError(RooAbsData::SumW2));
  background.plotOn(plot_MCbkg,VisualizeError(*r,2.0,kFALSE),FillColor(kYellow));
  background.plotOn(plot_MCbkg,VisualizeError(*r,1.0,kFALSE),FillColor(kGreen));  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  background.plotOn(plot_MCbkg);
  data_bkg->plotOn(plot_MCbkg,DataError(RooAbsData::SumW2));
  
  plot_MCbkg->Draw();
  
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  RooRealVar* alpha_0 = (RooRealVar*)(r->floatParsFinal().find("wdth"));
  RooRealVar* beta_0 = (RooRealVar*) (r->floatParsFinal().find("alpha"));

  Double_t x1= alpha_0->getVal();
  Double_t x2= beta_0->getVal();
  Double_t s1= alpha_0->getError();// magic number goes here
  Double_t s2= beta_0->getError();// magic number goes here
  Double_t rho= r->correlation("alpha", "wdth");
  r->Print("V");
  r->correlationMatrix().Print("v");
  r->covarianceMatrix().Print("v");

  RooEllipse *contour= new RooEllipse("contour",x1,x2,s1,s2,rho,500);
  contour->SetLineWidth(2) ;
  

  RooPlot *plot = new RooPlot(*alpha_0,*beta_0,wdth.getVal()-2*wdth.getError(),wdth.getVal()+2*wdth.getError(),
  			      alpha.getVal()-2*alpha.getError(),alpha.getVal()+2*alpha.getError());
  //RooPlot *plot = new RooPlot(*alpha_0,*beta_0,40,100,-1.5,.4);
  plot->addPlotable(contour);

  r->plotOn(plot,*alpha_0,*beta_0,"ME12");

  plot->Draw();

  a=cos(-rot)* x2+ sin(-rot)*x1;
  w=-sin(-rot)*x2 + cos(-rot)*x1;
  cout << "alpha true: " << a << endl;
  cout << "width true: " << w << endl;


  return r;
}

// btag0 -1.5545
//     Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
//   --------------------  ------------  --------------------------  --------
//                  alpha   -2.7640e+01   -2.7643e+01 +/-  3.16e+00  <none>
//                   wdth    2.3020e-02    2.3036e-02 +/-  1.15e-02  <none>

// btag1 -1.552 
//     Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
//   --------------------  ------------  --------------------------  --------
//                  alpha   -8.3331e+01   -8.3298e+01 +/-  3.55e+00  <none>
//                   wdth    2.4829e-01    2.4819e-01 +/-  4.09e-02  <none>

// btag1  -1.5574
//     Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
//   --------------------  ------------  --------------------------  --------
//                  alpha   -4.6462e+01   -4.6463e+01 +/-  1.44e+01  <none>
//                   wdth    1.7379e-04    1.3547e-04 +/-  7.90e-02  <none>


