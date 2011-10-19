#include "TH1F.h"
#include "TCanvas.h"
#include <Riostream.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>


#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
#include "RooExponential.h"

#include "PDFs/RooRodenbach.h"
#include "PDFs/RooFermi.h"
#include "PDFs/RooDoubleCB.h"
#include "PDFs/RooCB.h"
#include "PDFs/RooRelBW.h"

#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"

using namespace std;
//using namespace ROOT::Math;
using namespace RooFit;


//inputs: btag category, observed bkg yield (-> expected one for MC limit calc)
//        mass of Higgs, sigma of Higgs, 4 parameters of CB (depend on mass)

void make_roofitfiles(const std::string& dataset, int btag, int chan, double massH, double sigmaH, double &obs_yield, double &exp_yield, vector<double> cb_pars){

  gSystem->Load("libRooFit");
  cout<<"Trying to load custom PDFS ..."<<flush<<endl;
  gSystem->Load("libFFTW");

  string str_btag="dummyb";
  if(btag==0)str_btag="0b";
  else if(btag==1)str_btag="1b";
  else if(btag==2)str_btag="2b";
  else cout<<"Unrecognized number of btags: "<<btag<<endl;

  string str_chan="dummychan";
  if(chan==0)str_chan="ee";
  else if(chan==1)str_chan="mm";
  else cout<<"Unrecognized number of channels: "<<chan<<endl;


  //integration window
  double effWidth=sqrt(sigmaH*sigmaH+100);  // effective width used for defining window
  double fitRangeLow=99.;
  double fitRangeHigh=101.0; 

  ///expo functions/////////////////////////////////////////////////////
  /*if(massH-10*effWidth<230.0)  fitRangeLow=230.0;
  else fitRangeLow=massH-10*effWidth;

  if(massH+10*effWidth>800.0)  fitRangeHigh=800.;
  else fitRangeHigh=massH+10*effWidth;
  cout<<"----- FIT RANGE : "<<fitRangeLow<<" - "<< fitRangeHigh<<endl;*/
  /////////////////////////////////////////////////////////////////////

  //////////////////////////eps functions /////////////////////////////
  if(massH-10*effWidth<183.0)  fitRangeLow=183.0;
  else fitRangeLow=massH-10*effWidth;

  if(massH+10*effWidth>800.0)  fitRangeHigh=800.;
  else fitRangeHigh=massH+10*effWidth;

  cout<<"----- FIT RANGE : "<<fitRangeLow<<" - "<< fitRangeHigh<<endl;

  ////////////////////////////////////////////////////////////////

  std::ostringstream ossm;
  ossm<<massH; 
  string str_massH=ossm.str() ;

  // --------------------- initial values -----------------
  RooRealVar CMS_hzz2l2q_mZZ("CMS_hzz2l2q_mZZ", "zz inv mass",fitRangeLow,fitRangeHigh);

  // ==================== defining bkg PDF ==========================


  //exponential ////////////////////////////////////////////////////////////
  //par0 will be the overall bkgd normalization, used in the datacard
  /*string bkgp1name="CMS_hzz2l2q_bkg"+str_btag+"p1";
 
  double slope_Val[3]={-0.0150692, -0.0124851, -0.0127107};

  RooRealVar slope(bkgp1name.c_str(),"slope",slope_Val[btag],-1.0,0.0);
  slope.setConstant(kTRUE);

  RooExponential background("background","Exponential background",CMS_hzz2l2q_mZZ,slope);
  */

  //new shapes from Andrew///////////////////////
  /*vector<double> BKGparam;
  if(btag==0){
    BKGparam.push_back(186.405);   //fermi: cutOff
    BKGparam.push_back(5.25694);   //fermi: beta
    BKGparam.push_back(226.78);   //m
    BKGparam.push_back(44.1607);   //width
    BKGparam.push_back(254.615);   //alpha
  }else if(btag==1){
    BKGparam.push_back(184.04);   //fermi: cutOff
    BKGparam.push_back(3.225);   //fermi: beta
    BKGparam.push_back(200);   //m
    BKGparam.push_back(59.9061);   //width
    BKGparam.push_back(244.851);   //alpha
  }else if(btag==2){
    BKGparam.push_back(182.32);   //fermi: cutOff
    BKGparam.push_back(2.8492);   //fermi: beta
    BKGparam.push_back(200.002);   //m
    BKGparam.push_back(103.69);   //width
    BKGparam.push_back(646.467);   //alpha
  }
  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",BKGparam.at(0),0,1000);
  cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",BKGparam.at(1),0,50);
  beta.setConstant(kTRUE);
	     		       
  RooFermi fermi("fermi","fermi function",CMS_hzz2l2q_mZZ,cutOff,beta);
  // -------------------- double gauss ---------------------------
  string bkgp1name="CMS_hzz2l2q_bkg"+str_btag+"p1";
  string bkgp2name="CMS_hzz2l2q_bkg"+str_btag+"p2";
  string bkgp3name="CMS_hzz2l2q_bkg"+str_btag+"p3";
  RooRealVar m(bkgp1name.c_str(),bkgp1name.c_str(),BKGparam.at(2),200,1000);
  m.setConstant(kTRUE);
  RooRealVar wdth(bkgp2name.c_str(),bkgp2name.c_str(),BKGparam.at(3),0,1000);
  wdth.setConstant(kTRUE);
  RooRealVar alpha(bkgp3name.c_str(),bkgp3name.c_str(),BKGparam.at(4),200,1000); 
  alpha.setConstant(kTRUE);

  RooRodenbach Rod("Rod","Rod",CMS_hzz2l2q_mZZ,m,wdth,alpha);

  RooProdPdf background("background","background",RooArgSet(fermi,Rod));
  */
  //CB for EPS/////////////////////////////////////////////////////////////

  // ------------------------ fermi ------------------------------            
  vector<double> BKGparam;
////new parameters from matthias (17-10-2011)
// if(btag==0){//-1.5545 rot for Run2011A
//   BKGparam.push_back(186.41); //cutOff
//   BKGparam.push_back(5.257);	//beta	
//   BKGparam.push_back(222.72);	//CB_mean
//   BKGparam.push_back(27.6372);	//CB_wdth   2.0 /fb fits
//   BKGparam.push_back(25.775);	//CB_n	    13.39
//   BKGparam.push_back(-0.427402);	//CB_alpha  -.879
//   BKGparam.push_back(0.); //theta
//   //BKGparam.push_back(-1.5545); //theta
// }
// if(btag==1){//-1.552 rot for Run2011A
//   BKGparam.push_back(184.04);  //cutOff
//   BKGparam.push_back(3.225);	 //beta	
//   BKGparam.push_back(166.6);	 //CB_mean
//   BKGparam.push_back(83.3205);//CB_wdth   87.446    --- 2.0/fb fits
//   BKGparam.push_back(6.5670  );//CB_n	     21.499
//   BKGparam.push_back(-1.31797);//CB_alpha  -1.16716
//   BKGparam.push_back(0.); //theta
//   //BKGparam.push_back(-1.552); //theta
// }
// if(btag==2){//-1.5574 rot for Run2011A
//   BKGparam.push_back(182.32);  //cutOff
//   BKGparam.push_back(2.8492);	 //beta	
//   BKGparam.push_back(226.23);	 //CB_mean
//   BKGparam.push_back(46.4583); //CB_wdth  72.11  --- 2.0/fb fits
//   BKGparam.push_back(6.1580 ); //CB_n	    2.7583
//   BKGparam.push_back(-0.622234); //CB_alpha -1.654
//   BKGparam.push_back(0.); //theta
//   //BKGparam.push_back(-1.5574); //theta
// }


  // eps ones:
  if(btag==0){
     BKGparam.push_back(186.41); //cutoff
     BKGparam.push_back(5.257); //beta
     BKGparam.push_back(222.72); //mean
     BKGparam.push_back(-.116428); //width
     BKGparam.push_back(13.39); //n
     BKGparam.push_back(54.704); //alpha
     BKGparam.push_back(1.589); //theta
     }else if(btag==1){
     BKGparam.push_back(184.04);
     BKGparam.push_back(3.225);
     BKGparam.push_back(166.6);
     BKGparam.push_back(87.467);
     BKGparam.push_back(21.499);
     BKGparam.push_back(.24946);
     BKGparam.push_back(0.0162); //theta
     }else if(btag==2){
     BKGparam.push_back(182.32);
     BKGparam.push_back(2.84922);
     BKGparam.push_back(226.23);
     BKGparam.push_back(72.124);
     BKGparam.push_back(2.75833);
     BKGparam.push_back(-.45718);
     BKGparam.push_back(0.0166); //theta
     }

  //LP ones
  /*
  if(btag==0){
    BKGparam.push_back(186.41);  //cutOff 
    BKGparam.push_back(5.257);   //beta  
    BKGparam.push_back(222.72);  //CB_mean   
    BKGparam.push_back(53.876);  //CB_wdth  ---- uncorrelated  
    BKGparam.push_back(25.775);   //CB_n              
    BKGparam.push_back(.11276);   //CB_alpha ---- uncorrelated          
    BKGparam.push_back(.01752);   //theta                       
    }
  if(btag==1){
    BKGparam.push_back(184.04);  //cutOff    
    BKGparam.push_back(3.225);   //beta   
    BKGparam.push_back(166.6);   //CB_mean                                       
    BKGparam.push_back(94.7425);  //CB_wdth   --- uncorrelated     
    BKGparam.push_back(6.5670);  //CB_n                                                                                           
    BKGparam.push_back(.26330 );  //CB_alpha   --- uncorrelated  
    BKGparam.push_back(.0183  );  //theta                                         
    }
  if(btag==2){
    BKGparam.push_back(182.32);  //cutOff   
    BKGparam.push_back(2.8492);  //beta      
    BKGparam.push_back(226.23);  //CB_mean     
    BKGparam.push_back(39.41);    //CB_wdth   --- uncorrelated  
    BKGparam.push_back(6.1580);  //CB_n   
    BKGparam.push_back(-.2518);   //CB_alpha   --- uncorrelated   
    BKGparam.push_back(.0125);    //theta                                        
    }
  */

  // -------------------- fermi ---------------------------
  RooRealVar cutOff_BKG("cutOff_BKG","position of fermi",BKGparam.at(0),0,1000);
  cutOff_BKG.setConstant(kTRUE);
  RooRealVar beta_BKG("beta_BKG","width of fermi",BKGparam.at(1),0,50);
  beta_BKG.setConstant(kTRUE);
	     		       
  RooFermi fermi_BKG("fermi_BKG","fermi function",CMS_hzz2l2q_mZZ,cutOff_BKG,beta_BKG);
 // -------------------- double gauss ---------------------------
  //par0 will be the overall bkgd normalization, used in the datacard
  string bkgp1name="CMS_hzz2l2q_bkg"+str_btag+"p1"; //m
  string bkgp2name="CMS_hzz2l2q_bkg"+str_btag+"p2"; //width
  string bkgp3name="CMS_hzz2l2q_bkg"+str_btag+"p3"; //n
  string bkgp4name="CMS_hzz2l2q_bkg"+str_btag+"p4"; //alpha
  string bkgp5name="CMS_hzz2l2q_bkg"+str_btag+"p5"; //theta (rotation)

  RooRealVar m(bkgp1name.c_str(),bkgp1name.c_str(),BKGparam.at(2),100.,1000.);
  m.setConstant(kTRUE);
  RooRealVar wdth(bkgp2name.c_str(),bkgp2name.c_str(),BKGparam.at(3),0,1000);
  wdth.setConstant(kTRUE);
  RooRealVar n(bkgp3name.c_str(),bkgp3name.c_str(),BKGparam.at(4),0.,1001.);//2.75833,0,1000);
  n.setConstant(kTRUE);
  RooRealVar alpha(bkgp4name.c_str(),bkgp4name.c_str(),BKGparam.at(5),-100,100);  //0,100);  //,-100,0);
  alpha.setConstant(kTRUE);
  RooRealVar theta(bkgp5name.c_str(),bkgp5name.c_str(),BKGparam.at(6),-3.14159,3.14159); 
  theta.setConstant(kTRUE);

  RooCB CB_BKG("CB_BKG","Crystal ball",CMS_hzz2l2q_mZZ,m,wdth,alpha,n, theta);
  RooProdPdf background("background","background",RooArgSet(fermi_BKG,CB_BKG));
  
  ///////////////////////////////////////////////////////////////////////////////////

  ////Fill dataset with REAL DATA 

 RooRealVar nBTags("nBTags","nBTags",-1.,3.);
  RooRealVar eventWeight("eventWeight","eventWeight",0,100.);
  RooRealVar mZjj("mZjj","mZjj",0,150.);
  RooRealVar leptType("leptType","lepton type",-1,2);

  string btag_sel="dummy";
  if(btag==0)btag_sel="nBTags==0.0";
  else if(btag==1)btag_sel="nBTags==1.0";
  else if(btag==2)btag_sel="nBTags==2.0";
  else btag_sel="DUMMYnBTags==99.0";
  string lept_sel= chan==0 ? "leptType==1.0" :"leptType==0.0" ;//opposite convention btw Francesco and me
  string tree_sel= btag_sel+" && mZjj>75.0 && mZjj<105.0 && "+lept_sel;
  stringstream ossmzz1;
  ossmzz1 << float(fitRangeLow);
  string mzzcut="CMS_hzz2l2q_mZZ>"+ossmzz1.str(); 
  stringstream ossmzz2;
  ossmzz2 << float(fitRangeHigh);
  mzzcut+="&&CMS_hzz2l2q_mZZ<"+ossmzz2.str();
  cout<<"$$$$$$ TEMP SEL:  "<<mzzcut.c_str()<<"  $$$$$$$$$$$$$$$$$$$$$$ "<<fitRangeLow<<" - "<< fitRangeHigh<<endl;
  tree_sel+=" && "+mzzcut;
 
 
  /* TFile *dfile = new TFile("fileout-999invpb.root");
  //RooArgList arg1(mZZ);
  RooFormulaVar cut1("mycut1",tree_sel.c_str(),RooArgList(mZZ,nBTags,mZjj,leptType));

  RooDataSet *data_b=new RooDataSet("data_bkg","data_bkg",(TTree*)dfile->Get("tree_passedEvents"),
				    RooArgSet(mZZ,nBTags,mZjj,leptType),cut1,"eventWeight");
  obs_yield=double(data_b->numEntries());
  //RooDataSet *data_b = background.generate(x,int(obs_yield));
  cout<<"\nBTAG "<<btag<<"   OBS_YIELDS: "<<obs_yield<<" ->   "<<int(obs_yield)<<endl<<endl;*/

  std::string fileName = "convertedTree_"+dataset+".root";
  TFile* file = new TFile(fileName.c_str());
  RooFormulaVar cut1("mycut1",tree_sel.c_str(),RooArgList(CMS_hzz2l2q_mZZ,nBTags,mZjj,leptType));
  RooDataSet *dataset_obs_orig=new RooDataSet("dataset_obs_orig","dataset_obs_orig",(TTree*)file->Get("tree_passedEvents"),
					      RooArgSet(CMS_hzz2l2q_mZZ,nBTags,mZjj,leptType),
					      cut1,"eventWeight");

  obs_yield=double(dataset_obs_orig->numEntries());

  RooArgSet *newMZZargset= new RooArgSet(CMS_hzz2l2q_mZZ);
  RooDataSet *dataset_obs=(RooDataSet*) dataset_obs_orig->reduce(*newMZZargset);
  dataset_obs->SetName("dataset_obs");
  cout<<"Dataset entries: ORIG "<< dataset_obs_orig->sumEntries()<< "   NEW "<<dataset_obs->sumEntries()<<endl;
  // ----------------------------------------------

  // ====================== defining signal PDF =========================

  vector<double> param;
  if(btag==0){
    param.push_back(70.6146-.697703*massH+0.00212559*massH*massH-0.00000180624*massH*massH*massH);
    param.push_back(-5.967+0.05885*massH-0.00006977*massH*massH);
    param.push_back(1.0);
    param.push_back(3.38183-0.00421732*massH);
    param.push_back(1.0);
    param.push_back(-1.37066+0.0190719*massH-0.0000250673*massH*massH);
  }else if(btag==1){
    param.push_back(50.6113-.536745*massH+0.00174203*massH*massH-.00000152642*massH*massH*massH);	
    param.push_back(-4.08947+0.0385981*massH);											
    param.push_back(1.0);															
    param.push_back(.824239+0.00236893*massH);											
    param.push_back(1.0);															
    param.push_back(.444549+0.00495338*massH);                                                                                     
  }else if(btag==2){
    param.push_back(37.2265-0.391693*massH+0.00128062*massH*massH-.00000111444*massH*massH*massH); 
    param.push_back(-2.46367+0.022368*massH);											 
    param.push_back(1.0);															 
    param.push_back(1.61113+0.0015772*massH);											 
    param.push_back(1.0);															 
    param.push_back(1.95681+.00090888*massH);                                                                                       
  }

    for(int i=0; i<param.size(); i++){
    cout << "param[" << i << "]: " << param.at(i) << endl;
  }

  // -------------------- fermi ------------------------
  
  RooRealVar cutOff_SIG("cutOff_SIG","cutOff",190-32.5+65*massH/400,0,1000); 
  cutOff_SIG.setConstant(kTRUE);  
  RooRealVar g_SIG("g_SIG","g",5-12.5+25*massH/400,0,100); 
  g_SIG.setConstant(kTRUE);

  RooFermi fermi_SIG("fermi_SIG","fermi",CMS_hzz2l2q_mZZ,cutOff_SIG,g_SIG);

  // ------------------- fermi for high mass cutoff --------------

  RooRealVar cutOff2_SIG("cutOff2_SIG","cutOff2",700,0,1000);
  cutOff2_SIG.setConstant(kTRUE);
  RooRealVar g2_SIG("g2_SIG","g2",-70.0,-100.0,0.0);
  g2_SIG.setConstant(kTRUE);

  RooFermi fermi2_SIG("fermi2_SIG","fermi2",CMS_hzz2l2q_mZZ,cutOff2_SIG,g2_SIG);

  // ------------------- Relativistic BW --------------------------------
  //
 
  RooRealVar BW_mean("BW_mean", "mean",massH,0,1000);
  BW_mean.setConstant(kTRUE);
  RooRealVar BW_sigma("BW_sigma", "sigma",sigmaH,0,200);
  BW_sigma.setConstant(kTRUE);
  RooRealVar BW_n("BW_n","n",0.,0.,1.);
  BW_n.setConstant(kTRUE);

  RooRelBW BW("BW","Relativistic B-W",CMS_hzz2l2q_mZZ,BW_mean,BW_sigma,BW_n);

  // ------------------- Crystal Ball -------------------------------
  string sigp1name="CMS_hzz2l2q_sig"+str_btag+"p1"; //m
  string sigp2name="CMS_hzz2l2q_sig"+str_btag+"p2"; //width
   RooRealVar CB_mean(sigp1name.c_str(),sigp1name.c_str(),param[0],0.,100.);
  CB_mean.setConstant(kTRUE);
  RooRealVar CB_sigma(sigp2name.c_str(),sigp2name.c_str(),param[1],0.,100.);
  CB_sigma.setConstant(kTRUE);
  RooRealVar CB_alpha1("CB_alpha1","param 3 of CB",param[2],0.,100.);
  CB_alpha1.setConstant(kTRUE);
  RooRealVar CB_n1("CB_n1","param 4 of CB",param[3],0.,100.);
  CB_n1.setConstant(kTRUE);
  RooRealVar CB_alpha2("CB_alpha2","param 3 of CB",param[4],0.,100.);
  CB_alpha2.setConstant(kTRUE);
  RooRealVar CB_n2("CB_n2","param 4 of CB",param[5],0.,100.);
  CB_n2.setConstant(kTRUE);

  RooDoubleCB CB_SIG("CB_SIG","Crystal Ball",CMS_hzz2l2q_mZZ,CB_mean,CB_sigma,CB_alpha1,CB_n1,CB_alpha2,CB_n2);
  //------------------------ convolution -------------------------
  CMS_hzz2l2q_mZZ.setBins(10000,"fft");

  RooFFTConvPdf sig("sig","Rel-BW (X) CB",CMS_hzz2l2q_mZZ,BW,CB_SIG);
  sig.setBufferFraction(1.0);
  
  RooProdPdf signal("signal","signal",RooArgSet(sig,fermi_SIG,fermi2_SIG));

  //RooProdPdf signal_ggH(signal, "ggH_prodPDF");
  //RooProdPdf signal_VBF(signal, "qqH_prodPDF");
 

  //--- write everything into the workspace -------  

  RooWorkspace* w = new RooWorkspace("w","w");
  w->addClassDeclImportDir("/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms3/include/");
  // w->addClassDeclImportDir("/afs/cern.ch/user/w/whitbeck/scratch0/HiggsStats/newHiggsStats/CMSSW_4_1_3/src/HiggsAnalysis/CombinedLimit/data/PDFs/");
  //w->addClassDeclImportDir("/afs/cern.ch/user/b/bonato/scratch0/PhysAnalysis/CMSSW_4_2_3_patch5/src/ZJetsAnalysis/ZJetsAnalysisV1/test/statistical_tools/PDFs/");
  //w->addClassDeclImportDir("/afs/cern.ch/user/s/sbologne/scratch0/CMSSW/CMSSW_4_2_4/src/HiggsAnalysis/CombinedLimit/test/rotatedEPSForLP/PDFs/");
  w->addClassDeclImportDir("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/HZZlljj/HZZlljjAnalyzer/test/analysis/FIT/PDFs");

  w->importClassCode(RooFermi::Class(),kTRUE);
  w->importClassCode("RooFermi",kTRUE);
  w->importClassCode(RooRelBW::Class(),kTRUE);
  w->importClassCode("RooRelBW",kTRUE);
  w->importClassCode(RooDoubleCB::Class(),kTRUE);
  w->importClassCode("RooDoubleCB",kTRUE);
  w->importClassCode(RooCB::Class(),kTRUE);
  w->importClassCode("RooCB",kTRUE);
  //fro roorodenbach!!!!!!!!!!!!!!!!!!!!!!
  //w->importClassCode(RooRodenbach::Class(),kTRUE);
  //w->importClassCode("RooRodenbach",kTRUE);
  w->import(background);
  w->import(signal);
  // w->import(signal_ggH);
  // w->import(signal_VBF);
  w->import(*dataset_obs);

  //string outFileName="datacards_20110803_epsrotatedRange/"+str_massH+"/hzz2l2q_"+str_chan+str_btag+".input.root";
  string outFileName="datacards/"+str_massH+"/hzz2l2q_"+str_chan+str_btag+".input.root";
  
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");
  w->Write();
  outFile->Close();

  //calculate expected bkg events
  //eps functions ////////////////////////////
  RooRealVar CMS_hzz2l2q_mZZfull("CMS_hzz2l2q_mZZfull", "zz inv mass",183.0 ,800.0);
  //expo functions
  //RooRealVar CMS_hzz2l2q_mZZfull("CMS_hzz2l2q_mZZfull", "zz inv mass",230.0 ,800.0);

  //expo////////////////////////////////////////////////////////////
  //RooExponential backgroundFull("backgroundFull","Exponential background over Full range",CMS_hzz2l2q_mZZfull,slope);
  
  //eps ////////////////////////////////////////////////////////////////////////////////////////
  RooGenericPdf fermiFull("fermiFull","fermi function","1/(1+exp((@1-@0)/@2))",RooArgList(CMS_hzz2l2q_mZZfull,cutOff_BKG,beta_BKG));
  RooCB CBbkgFull("CBbkgFull","Crystal ball for background",CMS_hzz2l2q_mZZfull,m,wdth,alpha,n, theta);
  RooProdPdf backgroundFull("backgroundFull","backgroundFull",RooArgSet(fermiFull,CBbkgFull));

  //new shape andrew////////////////////////////////////////////////////////////////////////////////
  /*RooRodenbach RodFull("RodFull","Rod",CMS_hzz2l2q_mZZfull,m,wdth,alpha);
  RooGenericPdf fermiFull("fermiFull","fermi function","1/(1+exp((@1-@0)/@2))",RooArgList(CMS_hzz2l2q_mZZfull,cutOff,beta));
  RooProdPdf backgroundFull("backgroundFull","backgroundFull",RooArgSet(fermiFull,RodFull));
  */
  vector<double> EvtNorm;
  //first muon then electrons ///for expo /////////////////////////////////
  
  /*EvtNorm.push_back(chan==1? 228.10 : 200.12 );  // 0btag 
  EvtNorm.push_back(chan==1? 230.80 : 195.80 );  // 1btag 
  EvtNorm.push_back(chan==1?  16.82 :  13.79 );  // 2btag
  */
  //for eps///////////////////////////
  /*EvtNorm.push_back(chan==1? 345.7 : 286.4 );  // 0btag 
  EvtNorm.push_back(chan==1? 376.4 : 334.7 );  // 1btag 
  EvtNorm.push_back(chan==1? 24.3 : 20.3 );  // 2btag*/

  //for LP
  EvtNorm.push_back(chan==1? 575.85 : 490.54 );  // 0btag 
  EvtNorm.push_back(chan==1? 606.78 : 526.86 );  // 1btag 
  EvtNorm.push_back(chan==1? 41.14 : 35.37 );  // 2btag

  string mzzcut2="CMS_hzz2l2q_mZZfull>"+ossmzz1.str(); 
  mzzcut2+="&&CMS_hzz2l2q_mZZfull<"+ossmzz2.str();
  RooDataHist *BkgHisto = backgroundFull.generateBinned(CMS_hzz2l2q_mZZfull,EvtNorm.at(btag),kTRUE,kFALSE);  
  exp_yield=float( BkgHisto->sumEntries(mzzcut2.c_str() )  );
  cout<<"MH"<<massH<<"  With this cut: "<<mzzcut2.c_str()<<"  ===> "<<exp_yield<<"   TOT "<< BkgHisto->sumEntries( )<<"  mZZ<300 "<< BkgHisto->sumEntries("CMS_hzz2l2q_mZZfull<300.0" )<<endl;

  delete file;
  
}


