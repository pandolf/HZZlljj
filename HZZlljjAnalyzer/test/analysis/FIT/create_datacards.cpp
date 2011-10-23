#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooCB.h"
#include "RooDoubleCB.h"
#include "RooRelBW.h"
#include "RooFermi.h"
#include "RooProdPdf.h"
#include "RooFFTConvPdf.h"
#include "RooWorkspace.h"



struct HiggsParameters {

  float mH;
  float CSgg;
  float CSgg_p;
  float CSgg_m;
  float CSpdfgg_p;
  float CSpdfgg_m;
  float CSvbf;
  float CSvbf_p;
  float CSvbf_m;
  float CSpdfvbf_p;
  float CSpdfvbf_m;
  float Gamma;
  float BRHZZ;
  float BRZZ2l2q;

};


struct CBParameters {

  float m;
  float wdth;
  float alpha;
  float n;
  float theta;

};


double sign( double x ) {

  double returnSign = 0.;

  if( x>=0. ) returnSign =  1.;
  else returnSign =  -1.;

  return returnSign;

}


int convert_leptType( const std::string& leptType );
std::string leptType_datacards( const std::string& leptType_str );

void create_singleDatacard( const std::string& dataset, float mass, float lumi, const std::string& leptType_str, int nbtags, TF1* f1_eff_vs_mass );
HiggsParameters get_higgsParameters( float mass );
double linear_interp( double x, double x_old, double mass, double mH, double mH_old );
std::pair<double,double> get_massWindow( HiggsParameters hp );
TF1* get_eff_vs_mass( const std::string& leptType_str, int nbtags );
double get_observedYield( const std::string& dataset, const std::string& leptType_str, int nbtags, HiggsParameters hp, RooWorkspace* w );
double get_expectedYield_background( const std::string& dataset, const std::string& leptType_str, int nbtags, HiggsParameters hp, RooWorkspace* w );
void import_signalShape( int nbtags, HiggsParameters hp, RooWorkspace* w );

std::string systString( std::pair<double,double> systPair, double maxDiff=0.01 );
std::pair<double,double> getTheorSyst( double errMinus, double errPlus, double addMinus=0., double addPlus=0. );

std::pair<double,double> leptTriggerSyst( const std::string& leptType_str);
std::pair<double,double> leptEffSyst( const std::string& leptType_str);
std::pair<double,double> leptScaleSyst( const std::string& leptType_str);

std::pair<double,double> jetScaleSyst( double mass );
std::pair<double,double> bTagEffSyst( const std::string& leptType_str, int nbtags, double mass );









int main( int argc, char* argv[] ) {

  if( argc!=2 ) {
    std::cout << "USAGE: ./create_datacards [dataset]" << std::endl;
    exit(113);
  }

  
  std::string dataset(argv[1]);

  float lumi_ELE;
  float lumi_MU;
  if( dataset=="LP11") {
    lumi_ELE=1556.; //pb^-1
    lumi_MU =1615.; //pb^-1
  } else if( dataset=="Run2011A_FULL" ) {
    lumi_ELE=2100.; //pb^-1
    lumi_MU =2100.; //pb^-1
  } else {
    std::cout << "Unknown dataset '" << dataset << "'. Exiting." << std::endl;
    exit(333);
  }


  //first loop over available signal MC files to fit efficiency:
  TF1* f1_eff_vs_mass_MU_0btag = get_eff_vs_mass("MU", 0);
  TF1* f1_eff_vs_mass_MU_1btag = get_eff_vs_mass("MU", 1);
  TF1* f1_eff_vs_mass_MU_2btag = get_eff_vs_mass("MU", 2);

  TF1* f1_eff_vs_mass_ELE_0btag = get_eff_vs_mass("ELE", 0);
  TF1* f1_eff_vs_mass_ELE_1btag = get_eff_vs_mass("ELE", 1);
  TF1* f1_eff_vs_mass_ELE_2btag = get_eff_vs_mass("ELE", 2);



  std::ifstream ifs("masses.txt");
  
  while( ifs.good() ) {
    
    float mass;
    ifs >> mass;

    std::cout << std::endl;
    std::cout << "++++++++++++++++++++++" << std::endl;
    std::cout << "+++ MASS: " << mass << std::endl;
    std::cout << "++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    char mkdir_command[100];
    sprintf( mkdir_command, "mkdir -p datacardsPROVA/%.0f", mass);
    system(mkdir_command);

    create_singleDatacard( dataset, mass, lumi_ELE, "ELE", 0, f1_eff_vs_mass_ELE_0btag);
    create_singleDatacard( dataset, mass, lumi_ELE, "ELE", 1, f1_eff_vs_mass_ELE_1btag);
    create_singleDatacard( dataset, mass, lumi_ELE, "ELE", 2, f1_eff_vs_mass_ELE_2btag);
    create_singleDatacard( dataset, mass, lumi_MU,   "MU", 0, f1_eff_vs_mass_MU_0btag);
    create_singleDatacard( dataset, mass, lumi_MU,   "MU", 1, f1_eff_vs_mass_MU_1btag);
    create_singleDatacard( dataset, mass, lumi_MU,   "MU", 2, f1_eff_vs_mass_MU_2btag);

  } //while masses

  return 0;

}



void create_singleDatacard( const std::string& dataset, float mass, float lumi, const std::string& leptType_str, int nbtags, TF1* f1_eff_vs_mass ) {

  if( leptType_str!="ELE" && leptType_str!="MU" ) {
    std::cout << "Unkown Lept Type '" << leptType_str << "'. Exiting." << std::endl;
    exit(12333);
  }


  HiggsParameters hp = get_higgsParameters(mass);

  
  std::string leptType_forDatacard = (leptType_str=="ELE") ? "ee" : "mm";
  char suffix[100];
  sprintf( suffix, "%s%db", leptType_forDatacard.c_str(), nbtags);
  std::string suffix_str(suffix);

  char datacardName[400];
  sprintf( datacardName, "datacardsPROVA/%.0f/hzz2l2q_%s.%.0f.txt", mass, suffix, mass);


  std::ofstream ofs(datacardName);
  ofs << "# Simple counting experiment, with one signal and one background process" << std::endl;
  ofs << "#imax 1  number of channels" << std::endl;
  ofs << "#jmax 1  number of backgrounds" << std::endl;
  ofs << "#kmax *  number of nuisance parameters (sources of systematical uncertainties)" << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "shapes ggH CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root  w:signal" << std::endl;
  ofs << "shapes VBF CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root  w:signal " << std::endl;
  ofs << "shapes background CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root w:background " << std::endl;
  ofs << "shapes data_obs   CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root w:dataset_obs" << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "bin         CMS_hzz2l2q_" << suffix << std::endl;



  // now will compute the rates. in the meantime define the workspace, 
  // so that when computing signal/BG rate, write the signal/BG PDFs in there

  RooWorkspace* w = new RooWorkspace("w","w");
  w->addClassDeclImportDir("/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms3/include/");
  w->addClassDeclImportDir("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/HZZlljj/HZZlljjAnalyzer/test/analysis/FIT/PDFs");

  w->importClassCode(RooFermi::Class(),kTRUE);
  w->importClassCode("RooFermi",kTRUE);
  w->importClassCode(RooRelBW::Class(),kTRUE);
  w->importClassCode("RooRelBW",kTRUE);
  w->importClassCode(RooDoubleCB::Class(),kTRUE);
  w->importClassCode("RooDoubleCB",kTRUE);
  w->importClassCode(RooCB::Class(),kTRUE);
  w->importClassCode("RooCB",kTRUE);

  // immport signal shape: 
  // (should be a way of doing this automatically)
  // (in any case new shapes from matthias will change stuff here)
  import_signalShape( nbtags, hp, w );


  double observed = get_observedYield( dataset, leptType_str, nbtags, hp, w );
  ofs << "observation   " << observed << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "bin                CMS_hzz2l2q_" << suffix << "\tCMS_hzz2l2q_" << suffix << "\tCMS_hzz2l2q_" << suffix << std::endl;
  ofs << "process            ggH\t\t\tVBF\t\t\tbackground" << std::endl;
  ofs << "process            -1\t\t\t0\t\t\t1" << std::endl;




  float eff = f1_eff_vs_mass->Eval(hp.mH);
  float rate_gg   = eff*hp.CSgg *hp.BRHZZ*hp.BRZZ2l2q*lumi*0.5; //xsect has both ee and mm
  float rate_vbf  = eff*hp.CSvbf*hp.BRHZZ*hp.BRZZ2l2q*lumi*0.5; //xsect has both ee and mm

  double rate_background = get_expectedYield_background( dataset, leptType_str, nbtags, hp, w );

  ofs << "rate               " << rate_gg << "\t\t" << rate_vbf << "\t\t" << rate_background << std::endl;
  ofs << "------------ " << std::endl;


  // and now systematics:

  ofs << "lumi\t\t\tlnN\t1.045\t\t\t1.045\t\t\t1.0" << std::endl;

  std::pair<double,double> pdf_gg  = getTheorSyst( hp.CSpdfgg_m, hp.CSpdfgg_p, 0.04, 0.015 );
  ofs << "pdf_gg\t\tlnN\t" << systString(pdf_gg) << "\t1.0\t\t\t1.0" << std::endl;

  std::pair<double,double> pdf_vbf = getTheorSyst( hp.CSpdfvbf_m, hp.CSpdfvbf_p, 0.04, 0.015 );
  ofs << "pdf_qqbar\t\tlnN\t1.0\t\t\t" << systString(pdf_vbf) << "\t1.0" << std::endl;

  std::pair<double,double> QCDscale_ggH = getTheorSyst( hp.CSgg_m, hp.CSgg_p);
  ofs << "QCDscale_ggH\tlnN\t" << systString(QCDscale_ggH) << "\t1.0\t\t\t1.0" << std::endl;

  std::pair<double,double> QCDscale_qqH = getTheorSyst( hp.CSvbf_m, hp.CSvbf_p);
  ofs << "QCDscale_qqH\tlnN\t1.0\t\t\t" << systString(QCDscale_qqH) << "\t1.0" << std::endl;


  ofs << "CMS_trigger_" << leptType_datacards(leptType_str) << "\tlnN\t" << systString(leptTriggerSyst(leptType_str)) << "\t" << systString(leptTriggerSyst(leptType_str)) << "\t1.0" << std::endl;

  ofs << "CMS_eff_" << leptType_datacards(leptType_str) << "\t\tlnN\t" << systString(leptEffSyst(leptType_str)) << "\t" << systString(leptEffSyst(leptType_str)) << "\t1.0" << std::endl;

  ofs << "CMS_scale_" << leptType_datacards(leptType_str) << "\t\tlnN\t" << systString(leptScaleSyst(leptType_str)) << "\t" << systString(leptScaleSyst(leptType_str)) << "\t1.0" << std::endl;

  ofs << "CMS_scale_j\t\tlnN\t" << systString(jetScaleSyst(hp.mH)) << "\t" << systString(jetScaleSyst(hp.mH)) << "\t1.0" << std::endl;

  ofs << "CMS_eff_b\t\tlnN\t" << systString(bTagEffSyst(leptType_str, nbtags, hp.mH)) << "\t" << systString(bTagEffSyst(leptType_str, nbtags, hp.mH)) << "\t1.0" << std::endl;

  ofs << "CMS_hzz2l2q_pu\t\tlnN\t1.02\t\t\t1.02\t\t\t1.0" << std::endl;


  ofs.close();


  // datacard is done. now write rootfile:

}










TF1* get_eff_vs_mass( const std::string& leptType_str, int nbtags ) {

  int leptType_int = convert_leptType(leptType_str);

  ifstream ifsMC("massesMC.txt"); //the points at which we have MC samples

  TGraph* gr_eff_vs_mass = new TGraph(0);

  int iPoint=0;
  
  while( ifsMC.good() ) {

    double mass;
    ifsMC >> mass;

    HiggsParameters hp = get_higgsParameters(mass);

    std::pair<double,double> massWindow = get_massWindow(hp);
    double fitRangeLow = massWindow.first;
    double fitRangeHigh = massWindow.second;


    char signalfileName[800];
    sprintf( signalfileName, "HZZlljjRM_GluGluToHToZZTo2L2Q_M-%.0f_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_ALL.root", hp.mH );

    TFile* signalFile = TFile::Open(signalfileName);
    TTree* signalTree = (TTree*)signalFile->Get("tree_passedEvents");

    char signalCut_MW[500];
    sprintf( signalCut_MW, "HLTSF*PUWeight*( mZjj>75. && mZjj<105. && mZZ>%f && mZZ<%f && nBTags==%d && leptType==%d)", fitRangeLow, fitRangeHigh, nbtags, leptType_int);
    TH1D* h1_mZZ_signal = new TH1D("mZZ_signal", "", 65, 150., 800.);
    h1_mZZ_signal->Sumw2();
    signalTree->Project( "mZZ_signal", "mZZ", signalCut_MW );

    float signalYield = h1_mZZ_signal->Integral();

    TH1F* h1_generatedEventsPU = (TH1F*)signalFile->Get("nCounterPU");
    float generatedYield = h1_generatedEventsPU->Integral();

    float efficiency = signalYield/(generatedYield/3.); //three lept types in powheg

    gr_eff_vs_mass->SetPoint( iPoint++, mass, efficiency );

  } //while masses


  char functName[200];
  sprintf( functName, "eff_vs_mass_%s_%dbtag", leptType_str.c_str(), nbtags );
  TF1* f1_eff_vs_mass = new TF1(functName, "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 150., 700.);
  gr_eff_vs_mass->Fit(f1_eff_vs_mass, "RQ");

  system("mkdir -p EfficiencyFits");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  gr_eff_vs_mass->SetMarkerStyle(20);
  gr_eff_vs_mass->SetMarkerSize(1.3);
  gr_eff_vs_mass->Draw("APE");

  char canvasName[500];
  sprintf( canvasName, "EfficiencyFits/effFit_%s_%dbtag.eps", leptType_str.c_str(), nbtags);
  c1->SaveAs(canvasName);

  delete c1;

  return f1_eff_vs_mass;

}




int convert_leptType( const std::string& leptType ) {

  if( leptType!="ELE" && leptType!="MU" ) {
    std::cout << "WARNING!!! LeptType '" << leptType << "' is NOT supported!!! Returning -1." << std::endl;
    return -1;
  }

  int leptType_int = (leptType=="MU" ) ? 0 : 1;

  return leptType_int;

}



std::string leptType_datacards( const std::string& leptType_str ) {
 
  std::string returnString="";

  if( leptType_str=="MU" ) returnString = "m";
  if( leptType_str=="ELE" ) returnString = "e";

  return returnString;

}


HiggsParameters get_higgsParameters( float mass ) {

  std::string nameXsecFile = "xsect_higgs_173points_new.txt";
  std::ifstream xsect_file(nameXsecFile.c_str());

  if (! xsect_file.is_open()) { 
    std::cout << "Failed to open file with xsections"<<endl;
    exit(13111);
  }

  xsect_file.clear();
  xsect_file.seekg(0);
 
  float mH, CSgg, CSgg_p, CSgg_m, CSpdfgg_p,CSpdfgg_m,CSvbf, CSvbf_p, CSvbf_m,CSpdfvbf_p,CSpdfvbf_m, 
        Gamma, BRHZZ, BRZZ2l2q;

  float mH_old, CSgg_old, CSgg_p_old, CSgg_m_old, CSpdfgg_p_old,CSpdfgg_m_old,CSvbf_old, CSvbf_p_old, CSvbf_m_old,CSpdfvbf_p_old,CSpdfvbf_m_old, 
        Gamma_old, BRHZZ_old, BRZZ2l2q_old;

  HiggsParameters hp;

  while(xsect_file.good()) {

    mH_old = mH;
    CSgg_old = CSgg;
    CSgg_p_old = CSgg_p;
    CSgg_m_old = CSgg_m;
    CSpdfgg_p_old = CSpdfgg_p;
    CSpdfgg_m_old = CSpdfgg_m;
    CSvbf_old = CSvbf;
    CSvbf_p_old = CSvbf_p;
    CSvbf_m_old = CSvbf_m;
    CSpdfvbf_p_old = CSpdfvbf_p;
    CSpdfvbf_m_old = CSpdfvbf_m;
    Gamma_old = Gamma;
    BRHZZ_old = BRHZZ;
    BRZZ2l2q_old = BRZZ2l2q;

    xsect_file >> mH >> CSgg>> CSgg_p >> CSgg_m >> CSpdfgg_p >> CSpdfgg_m >> CSvbf >> CSvbf_p >> CSvbf_m >> CSpdfvbf_p >> CSpdfvbf_m >> Gamma >> BRHZZ >> BRZZ2l2q;

    if( mH == mass ) {

      hp.mH = mH;
      hp.CSgg = CSgg;
      hp.CSgg_p = CSgg_p;
      hp.CSgg_m = CSgg_m;
      hp.CSpdfgg_p = CSpdfgg_p;
      hp.CSpdfgg_m = CSpdfgg_m;
      hp.CSvbf = CSvbf;
      hp.CSvbf_p = CSvbf_p;
      hp.CSvbf_m = CSvbf_m;
      hp.CSpdfvbf_p = CSpdfvbf_p;
      hp.CSpdfvbf_m = CSpdfvbf_m;
      hp.Gamma = Gamma;
      hp.BRHZZ = BRHZZ;
      hp.BRZZ2l2q = BRZZ2l2q;
    
      break;

    } if( mass > mH_old && mass < mH ) {

      hp.mH = mass;
      hp.CSgg       = linear_interp( CSgg, CSgg_old, mass, mH, mH_old );
      hp.CSgg_p     = linear_interp( CSgg_p, CSgg_p_old, mass, mH, mH_old );
      hp.CSgg_m     = linear_interp( CSgg_m, CSgg_m_old, mass, mH, mH_old );
      hp.CSpdfgg_p  = linear_interp( CSpdfgg_p, CSpdfgg_p_old, mass, mH, mH_old );
      hp.CSpdfgg_m  = linear_interp( CSpdfgg_m, CSpdfgg_m_old, mass, mH, mH_old );
      hp.CSvbf      = linear_interp( CSvbf, CSvbf_old, mass, mH, mH_old );
      hp.CSvbf_p    = linear_interp( CSvbf_p, CSvbf_p_old, mass, mH, mH_old );
      hp.CSvbf_m    = linear_interp( CSvbf_m, CSvbf_m_old, mass, mH, mH_old );
      hp.CSpdfvbf_p = linear_interp( CSpdfvbf_p, CSpdfvbf_p_old, mass, mH, mH_old );
      hp.CSpdfvbf_m = linear_interp( CSpdfvbf_m, CSpdfvbf_m_old, mass, mH, mH_old );
      hp.Gamma      = linear_interp( Gamma, Gamma_old, mass, mH, mH_old );
      hp.BRHZZ      = linear_interp( BRHZZ, BRHZZ_old, mass, mH, mH_old );
      hp.BRZZ2l2q   = linear_interp( BRZZ2l2q, BRZZ2l2q_old, mass, mH, mH_old );
    
      break;

    } // if

  } //while ifs

  return hp;

}



std::pair<double,double> get_massWindow( HiggsParameters hp ) {

  double effWidth = sqrt(hp.Gamma*hp.Gamma+100.);
  double fitRangeLow  = (hp.mH-10.*effWidth<183.) ? 183. : hp.mH-10.*effWidth;
  double fitRangeHigh = (hp.mH+10.*effWidth>800.) ? 800. : hp.mH+10.*effWidth;

  std::pair<double,double> returnPair;
  returnPair.first = fitRangeLow;
  returnPair.second = fitRangeHigh;

  return returnPair;

}



double linear_interp( double x, double x_old, double mass, double mH, double mH_old ) {

  return (x_old + ( x-x_old ) * ( mass-mH_old ) / ( mH-mH_old ));

}



double get_observedYield( const std::string& dataset, const std::string& leptType_str, int nbtags, HiggsParameters hp, RooWorkspace* w ) {

  int leptType_int = convert_leptType(leptType_str);

  std::string dataFileName = "HZZlljjRM_DATA_" + dataset + "_optLD_looseBTags_v2_ALL.root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  TTree* tree_data = (TTree*)dataFile->Get("tree_passedEvents");
  tree_data->GetBranch("mZZ")->SetName("CMS_hzz2l2q_mZZ"); //needed for combination


  //integration window
  std::pair<double,double> massWindow = get_massWindow(hp);
  double fitRangeLow = massWindow.first;
  double fitRangeHigh = massWindow.second;
  
  RooRealVar CMS_hzz2l2q_mZZ("CMS_hzz2l2q_mZZ", "zz inv mass", fitRangeLow, fitRangeHigh);

  RooRealVar nBTags("nBTags","nBTags",-1.,3.);
  RooRealVar eventWeight("eventWeight","eventWeight",0,100.);
  RooRealVar mZjj("mZjj","mZjj",0,150.);
  RooRealVar leptType("leptType","lepton type",-1,2);

  char selection_MW[900];
  sprintf( selection_MW, "mZjj>75. && mZjj<105. && nBTags==%d && leptType==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f", nbtags, leptType_int, fitRangeLow, fitRangeHigh);


  RooFormulaVar massWindowSelection("massWindowSelection", selection_MW, RooArgList(CMS_hzz2l2q_mZZ,nBTags,mZjj,leptType));
  RooDataSet *dataset_obs = new RooDataSet("dataset_obs_orig", "dataset_obs_orig", tree_data,
                                     RooArgSet(CMS_hzz2l2q_mZZ, nBTags, mZjj, leptType, eventWeight),
                                     selection_MW, "eventWeight");


  // import it in the workspace:
  w->import(*dataset_obs);


  return double(dataset_obs->numEntries());

}



double get_expectedYield_background( const std::string& dataset, const std::string& leptType_str, int nbtags, HiggsParameters hp, RooWorkspace* w ) {

  int leptType_int = convert_leptType(leptType_str);




  // read background parametrizations from fit results file
  char fitResultsFile[900];
  sprintf( fitResultsFile, "FitSidebands_%s/fitresultsDATA_%dbtag.txt", dataset.c_str(), nbtags);
  
  ifstream ifs(fitResultsFile);
  CBParameters cbp;
  float fermi_cutoff0;
  float fermi_beta0;

  while( ifs.good() ) {

    std::string varName;
    float value, error;
    ifs >> varName >> value >> error;

    if( varName=="beta" ) {
      fermi_beta0 = value;
    }
    if( varName=="cutOff" ) {
      fermi_cutoff0 = value;
    }
    if( varName=="m" ) {
      cbp.m = value;
    }
    if( varName=="n" ) {
      cbp.n = value;
    }
    if( varName=="alpha_rot" ) {
    //if( varName=="alpha" ) {
      cbp.alpha = value;
    }
    if( varName=="wdth_rot" ) {
    //if( varName=="wdth" ) {
      cbp.wdth = value;
    }
    if( varName=="theta_best" ) {
    //if( varName=="theta" ) {
      cbp.theta = value;
    }

  } // while ifs fitresults


  //integration window
  std::pair<double,double> massWindow = get_massWindow(hp);
  double fitRangeLow = massWindow.first;
  double fitRangeHigh = massWindow.second;


  // define mZZ variable
  RooRealVar CMS_hzz2l2q_mZZ("CMS_hzz2l2q_mZZ", "zz inv mass",fitRangeLow,fitRangeHigh);
  RooRealVar CMS_hzz2l2q_mZZfull("CMS_hzz2l2q_mZZfull", "zz inv mass", 183., 800.);


  // define background PDF:

  RooRealVar fermi_cutoff("fermi_cutoff", "position of fermi", fermi_cutoff0, 0., 1000.);
  fermi_cutoff.setConstant(kTRUE);
  RooRealVar fermi_beta("fermi_beta", "width of fermi", fermi_beta0, 0., 50.);
  fermi_beta.setConstant(kTRUE);

  RooFermi fermi_BKG("fermi_BKG", "fermi function", CMS_hzz2l2q_mZZ, fermi_cutoff, fermi_beta);


  char backgroundParName[100];
  sprintf( backgroundParName, "CMS_hzz2l2q_bkg%dp1", nbtags);
  RooRealVar m(backgroundParName, backgroundParName, cbp.m, 100., 1000.);
  m.setConstant(kTRUE);
  sprintf( backgroundParName, "CMS_hzz2l2q_bkg%dp2", nbtags);
  RooRealVar wdth(backgroundParName, backgroundParName, cbp.wdth, 0., 1000.);
  wdth.setConstant(kTRUE);
  sprintf( backgroundParName, "CMS_hzz2l2q_bkg%dp3", nbtags);
  RooRealVar n(backgroundParName, backgroundParName, cbp.n, 0., 1001.);
  n.setConstant(kTRUE);
  sprintf( backgroundParName, "CMS_hzz2l2q_bkg%dp4", nbtags);
  RooRealVar alpha(backgroundParName, backgroundParName, cbp.alpha, -100., 100.);
  alpha.setConstant(kTRUE);
  sprintf( backgroundParName, "CMS_hzz2l2q_bkg%dp5", nbtags);
  RooRealVar theta(backgroundParName, backgroundParName, cbp.theta, -3.14159, 3.14159); 
  theta.setConstant(kTRUE);
  

  RooCB CB_BKG("CB_BKG", "Crystal ball", CMS_hzz2l2q_mZZ, m, wdth, alpha, n, theta);
  RooProdPdf background("background", "background", RooArgSet(fermi_BKG,CB_BKG));

  // import in the workspace:
  w->import(background);




  // now get the expected yield:

  char alphaFileName[300];
  sprintf(alphaFileName, "alphaFile_%s_%dbtag_ALL.root", dataset.c_str(), nbtags);
  TFile* alphaFile = TFile::Open(alphaFileName);
  TTree* treeSidebandsDATA_alphaCorr = (TTree*)alphaFile->Get("sidebandsDATA_alpha");
  //treeSidebandsDATA_alphaCorr->GetBranch("mZZ")->SetName("CMS_hzz2l2q_mZZfull"); //makes things easier
  TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, 150., 800.);
  h1_mZZ_sidebands_alpha->Sumw2();
  char sidebandsCut_alpha[500];
  sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d && leptType==%d)", nbtags, leptType_int);
  treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "mZZ", sidebandsCut_alpha);
  float EvtNorm = h1_mZZ_sidebands_alpha->Integral();

  //RooGenericPdf fermiFull("fermiFull","fermi function","1/(1+exp((@1-@0)/@2))",RooArgList(CMS_hzz2l2q_mZZfull, fermi_cutoff, fermi_beta));
  RooFermi fermi_BKG_FULL("fermi_BKG_FULL","fermi function", CMS_hzz2l2q_mZZfull, fermi_cutoff, fermi_beta);
  RooCB CB_BKG_FULL("CBbkgFull","Crystal ball for background", CMS_hzz2l2q_mZZfull, m, wdth, alpha, n, theta);
  RooProdPdf backgroundFull("backgroundFull","backgroundFull",RooArgSet(fermi_BKG_FULL,CB_BKG_FULL));

  RooDataHist *BkgHisto = backgroundFull.generateBinned(CMS_hzz2l2q_mZZfull,EvtNorm,kTRUE,kFALSE);  
  //RooDataHist *BkgHisto = background.generateBinned(CMS_hzz2l2q_mZZ,EvtNorm,kTRUE,kFALSE);  


/*   CANT MAKE THIS WORK, GETTING MAD

  // check that everything is correct with a plot:
  RooRealVar mZZ("mZZ","mZZ",fitRangeLow,fitRangeHigh);
  RooRealVar nBTags("nBTags","nBTags",-1.,3.);
  RooRealVar leptType("leptType","leptType",0,1);
  RooRealVar mZjj("mZjj","mZjj",0,200.);
  TH1D* h1_mZZ_signal = new TH1D("mZZ_signal", "", 30, 150., 750.);
  char signalCut[500];
  //sprintf(signalCut, "eventWeight_alpha*(mZjj>75. && mZjj<105. && nBTags==%d && leptType==%d)", nbtags, leptType_int);
  sprintf(signalCut, "(mZjj>75. && mZjj<105. && nBTags==%d && leptType==%d)", nbtags, leptType_int);
  //treeSidebandsDATA_alphaCorr->Project("mZZ_signal", "mZZ", signalCut);
  RooDataSet *data_signalRegion = new RooDataSet("data_signalRegion", "data_signalRegion", treeSidebandsDATA_alphaCorr,
                                                 RooArgSet(CMS_hzz2l2q_mZZfull,leptType,nBTags,mZjj), signalCut);

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  RooPlot *plot_MCbkg = CMS_hzz2l2q_mZZfull.frame();
  data_signalRegion->plotOn(plot_MCbkg, RooFit::Binning(65));
  //backgroundFull.plotOn(plot_MCbkg, RooFit::Normalization(EvtNorm));
  BkgHisto->plotOn(plot_MCbkg);
  //h1_mZZ_signal->SetMarkerStyle(20);
  //h1_mZZ_signal->SetMarkerSize(1.4);
  //h1_mZZ_signal->SetMarkerColor(kBlack);
  //h1_mZZ_signal->Draw("e");
  plot_MCbkg->Draw("");

  char canvasName[500];
  sprintf( canvasName, "FitSidebands_%s/crossCheck_%dbtag_%s.eps", dataset.c_str(), nbtags, leptType_str.c_str() );
  c1->SaveAs(canvasName);
*/



  char selection_massWindow[300];
  sprintf( selection_massWindow, "CMS_hzz2l2q_mZZfull>%f && CMS_hzz2l2q_mZZfull<%f", fitRangeLow, fitRangeHigh);
  //sprintf( selection_massWindow, "CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f", fitRangeLow, fitRangeHigh);
  double expectedYield_background = double( BkgHisto->sumEntries(selection_massWindow)  );

//// put it back as it was:
//treeSidebandsDATA_alphaCorr->GetBranch("CMS_hzz2l2q_mZZfull")->SetName("mZZ");

  return expectedYield_background;

}


//float CSggxs =1000.0*CSgg  *BRHZZ *BRZZ2l2q;
//float CSvbfxs=1000.0*CSvbf*BRHZZ  *BRZZ2l2q;
////   cout<<"THEOR xsect: "<<CSgg<<" , Vbf= "<<CSvbf<<"  BRHZZ="<<BRHZZ<<"  BRZZ2L2q="<<BRZZ2l2q<<endl;
//vector<float> myxsect;
//myxsect.push_back(CSggxs);
//myxsect.push_back(CSvbfxs);





void import_signalShape( int nbtags, HiggsParameters hp, RooWorkspace* w ) {

  //integration window
  std::pair<double,double> massWindow = get_massWindow(hp);
  double fitRangeLow = massWindow.first;
  double fitRangeHigh = massWindow.second;

  RooRealVar CMS_hzz2l2q_mZZ("CMS_hzz2l2q_mZZ", "zz inv mass", fitRangeLow, fitRangeHigh);


  // ====================== defining signal PDF =========================

  vector<double> param;
  if(nbtags==0){
    param.push_back(70.6146-.697703*hp.mH+0.00212559*hp.mH*hp.mH-0.00000180624*hp.mH*hp.mH*hp.mH);
    param.push_back(-5.967+0.05885*hp.mH-0.00006977*hp.mH*hp.mH);
    param.push_back(1.0);
    param.push_back(3.38183-0.00421732*hp.mH);
    param.push_back(1.0);
    param.push_back(-1.37066+0.0190719*hp.mH-0.0000250673*hp.mH*hp.mH);
  }else if(nbtags==1){
    param.push_back(50.6113-.536745*hp.mH+0.00174203*hp.mH*hp.mH-.00000152642*hp.mH*hp.mH*hp.mH);	
    param.push_back(-4.08947+0.0385981*hp.mH);
    param.push_back(1.0);
    param.push_back(.824239+0.00236893*hp.mH);
    param.push_back(1.0);
    param.push_back(.444549+0.00495338*hp.mH);                                                                                     
  }else if(nbtags==2){
    param.push_back(37.2265-0.391693*hp.mH+0.00128062*hp.mH*hp.mH-.00000111444*hp.mH*hp.mH*hp.mH); 
    param.push_back(-2.46367+0.022368*hp.mH);
    param.push_back(1.0);
    param.push_back(1.61113+0.0015772*hp.mH);
    param.push_back(1.0);
    param.push_back(1.95681+.00090888*hp.mH);                                                                                       
  }

    for(int i=0; i<param.size(); i++){
    cout << "param[" << i << "]: " << param.at(i) << endl;
  }

  // -------------------- fermi ------------------------
  
  RooRealVar cutOff_SIG("cutOff_SIG","cutOff",190-32.5+65*hp.mH/400,0,1000); 
  cutOff_SIG.setConstant(kTRUE);  
  RooRealVar g_SIG("g_SIG","g",5-12.5+25*hp.mH/400,0,100); 
  g_SIG.setConstant(kTRUE);

  RooFermi fermi_SIG("fermi_SIG","fermi",CMS_hzz2l2q_mZZ,cutOff_SIG,g_SIG);

  // ------------------- fermi for high mass cutoff --------------

  RooRealVar cutOff2_SIG("cutOff2_SIG","cutOff2",700,0,1000);
  cutOff2_SIG.setConstant(kTRUE);
  RooRealVar g2_SIG("g2_SIG","g2",-70.0,-100.0,0.0);
  g2_SIG.setConstant(kTRUE);

  RooFermi fermi2_SIG("fermi2_SIG","fermi2",CMS_hzz2l2q_mZZ,cutOff2_SIG,g2_SIG);

  // ------------------- Relativistic BW --------------------------------
 
  RooRealVar BW_mean("BW_mean", "mean", hp.mH, 0., 1000.);
  BW_mean.setConstant(kTRUE);
  RooRealVar BW_sigma("BW_sigma", "sigma", hp.Gamma, 0., 200.);
  BW_sigma.setConstant(kTRUE);
  RooRealVar BW_n("BW_n","n",0.,0.,1.);
  BW_n.setConstant(kTRUE);

  RooRelBW BW("BW","Relativistic B-W",CMS_hzz2l2q_mZZ,BW_mean,BW_sigma,BW_n);

  // ------------------- Crystal Ball -------------------------------
  char signalParName[100];
  sprintf( signalParName, "CMS_hzz2l2q_sig%dp1", nbtags);
  RooRealVar CB_mean(signalParName,signalParName,param[0],0.,100.);
  CB_mean.setConstant(kTRUE);
  sprintf( signalParName, "CMS_hzz2l2q_sig%dp2", nbtags);
  RooRealVar CB_sigma(signalParName,signalParName,param[1],0.,100.);
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

  // import the hell out of it
  w->import(signal);

}




std::string systString( std::pair<double,double> systPair, double maxDiff ) {

  double syst_ave = 1. + 0.5*(fabs(systPair.first-1.) + fabs(systPair.second-1.));
  
  char syst_char[100];
  if( fabs(syst_ave-systPair.second)/syst_ave < maxDiff )
    sprintf( syst_char, "%f    ", syst_ave );
  else
    sprintf( syst_char, "%f/%f", systPair.first, systPair.second );

  std::string syst_str(syst_char);

  return syst_str;

}
 

std::pair<double,double> getTheorSyst( double errMinus, double errPlus, double addMinus, double addPlus ) {

  float systPlus  = sign(errPlus) *sqrt(errPlus*errPlus   + addPlus*addPlus);
  float systMinus = sign(errMinus)*sqrt(errMinus*errMinus + addMinus*addMinus);

  systPlus  += 1.;
  systMinus += 1.;

  std::pair<double,double> returnPair;
  returnPair.first = systMinus;
  returnPair.second = systPlus;

  return returnPair;

}


std::pair<double,double> leptTriggerSyst( const std::string& leptType_str) {

  double syst;

  if( leptType_str=="MU" )  syst = 1.02;
  if( leptType_str=="ELE" ) syst = 1.01;

  std::pair<double,double> returnPair;
  returnPair.first  = syst; //symmetrical for now
  returnPair.second = syst;

  return returnPair;

}

std::pair<double,double> leptEffSyst( const std::string& leptType_str) {

  double syst;

  if( leptType_str=="MU" )  syst = 1.008;
  if( leptType_str=="ELE" ) syst = 1.034;

  std::pair<double,double> returnPair;
  returnPair.first  = syst; //symmetrical for now
  returnPair.second = syst;

  return returnPair;

}

std::pair<double,double> leptScaleSyst( const std::string& leptType_str) {

  double syst;

  if( leptType_str=="MU" )  syst = 1.01;
  if( leptType_str=="ELE" ) syst = 1.03;

  std::pair<double,double> returnPair;
  returnPair.first  = syst; //symmetrical for now
  returnPair.second = syst;

  return returnPair;

}


std::pair<double,double> jetScaleSyst( double mass ) {

  float p0= 8.3  , p1=-0.0215 ;
  float m0=-8.6, m1=0.02 ;

  std::pair<double,double> returnPair;
  returnPair.first  = 1.0 + 0.01*(m0+m1*mass);
  returnPair.second = 1.0 + 0.01*(p0+p1*mass);

  return returnPair;

}


std::pair<double,double> bTagEffSyst( const std::string& leptType_str, int nbtags, double mass ) {

  float p0=0.0, p1=0.0;
  float m0=0.0, m1=0.0;

  if( leptType_str=="ELE" ) {
    if(nbtags==0){
      p0=0.983256647923;
      p1=-0.0000883532570978;
      m0=1.02907356239;
      m1=0.0000713061639147;
    }
    else if(nbtags==1){
      p0=1.04446110956;
      p1=-0.0000195508160829;
      m0=0.940063743731;
      m1=0.0000737044467898;
    }
    else if(nbtags==2){
      p0=1.13365470372;
      p1=0.00000584572717646;
      m0=0.82161771535;
      m1=-0.0000161054152592;
    }
  }
  else if(leptType_str=="MU" ) {
    if(nbtags==0){
      p0=0.984636818312;
      p1=-0.0000898705296203;
      m0=1.02836905579;
      m1= 0.0000726807344479;
    }
    else if(nbtags==1){
      p0=1.04385580002;
      p1=0.0000206096278947;
      m0=0.942713582987;
      m1=0.0000719882385098;
    }
    else if(nbtags==2){
      p0=1.1333366687;
      p1=0.00000462542786413;
      m0=0.813316607701;
      m1=-0.00000205840248842;
    }
  }

  std::pair<double,double> returnPair;
  returnPair.first  = m1*mass+m0;
  returnPair.second = p1*mass+p0;

  return returnPair;

}
