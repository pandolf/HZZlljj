#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooCB.h"
#include "RooFermi.h"
//#include "RooGenericPdf.h"
#include "RooProdPdf.h"



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


int convert_leptType( const std::string& leptType );
void create_singleDatacard( const std::string& dataset, float mass, float lumi, const std::string& leptType_str, int nbtags );
HiggsParameters get_higgsParameters( float mass );
double get_observedYield( const std::string& dataset, const std::string& leptType_str, int nbtags, float mass, float gammaHiggs );
double get_expectedYield_background( const std::string& dataset, const std::string& leptType_str, int nbtags, float mass, float gammaHiggs );





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


  std::ifstream ifs("masses.txt");
  
  while( ifs.good() ) {
    
    float mass;
    ifs >> mass;

    std::cout << "+++ MASS: " << mass << std::endl;

    char mkdir_command[100];
    sprintf( mkdir_command, "mkdir -p datacardsPROVA/%.0f", mass);
    system(mkdir_command);

    create_singleDatacard( dataset, mass, lumi_ELE, "ELE", 0);
    create_singleDatacard( dataset, mass, lumi_ELE, "ELE", 1);
    create_singleDatacard( dataset, mass, lumi_ELE, "ELE", 2);
    create_singleDatacard( dataset, mass, lumi_MU,  "MU", 0);
    create_singleDatacard( dataset, mass, lumi_MU,  "MU", 1);
    create_singleDatacard( dataset, mass, lumi_MU,  "MU", 2);

  } //while masses

  return 0;

}



void create_singleDatacard( const std::string& dataset, float mass, float lumi, const std::string& leptType_str, int nbtags ) {

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

  double observed = get_observedYield( dataset, leptType_str, nbtags, mass, hp.Gamma );
  ofs << "observation   " << observed << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "bin                CMS_hzz2l2q_" << suffix << "\tCMS_hzz2l2q_" << suffix << "\tCMS_hzz2l2q_" << suffix << std::endl;
  ofs << "process            ggH\t\t\tVBF\t\t\tbackground" << std::endl;
  ofs << "process            -1\t\t\t0\t\t\t1" << std::endl;

  double rate_gg(0.), rate_VBF(0.);

  double rate_background = get_expectedYield_background( dataset, leptType_str, nbtags, mass, hp.Gamma );

  ofs << "rate               " << rate_gg << "\t\t\t" << rate_VBF << "\t\t\t" << rate_background << std::endl;

  ofs.close();

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
 

  HiggsParameters hp;

  while(xsect_file.good()) {
    xsect_file >> hp.mH >> hp.CSgg>> hp.CSgg_p >> hp.CSgg_m >> hp.CSpdfgg_p >> hp.CSpdfgg_m >> hp.CSvbf >> hp.CSvbf_p >> hp.CSvbf_m >> hp.CSpdfvbf_p >> hp.CSpdfvbf_m >> hp.Gamma >> hp.BRHZZ >> hp.BRZZ2l2q;
    if( hp.mH == mass ) break;
  }

  return hp;

}


double get_observedYield( const std::string& dataset, const std::string& leptType_str, int nbtags, float mass, float gammaHiggs ) {

  int leptType_int = convert_leptType(leptType_str);

  std::string dataFileName = "HZZlljjRM_DATA_" + dataset + "_optLD_looseBTags_v2_ALL.root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  TTree* tree_data = (TTree*)dataFile->Get("tree_passedEvents");
  tree_data->GetBranch("mZZ")->SetName("CMS_hzz2l2q_mZZ"); //needed for combination


  //integration window
  double effWidth = sqrt(gammaHiggs*gammaHiggs+100.);
  double fitRangeLow  = (mass-10.*effWidth<183.) ? 183. : mass-10.*effWidth;
  double fitRangeHigh = (mass+10.*effWidth>800.) ? 800. : mass+10.*effWidth;
  
  RooRealVar CMS_hzz2l2q_mZZ("CMS_hzz2l2q_mZZ", "zz inv mass", fitRangeLow, fitRangeHigh);

  RooRealVar nBTags("nBTags","nBTags",-1.,3.);
  RooRealVar eventWeight("eventWeight","eventWeight",0,100.);
  RooRealVar mZjj("mZjj","mZjj",0,150.);
  RooRealVar leptType("leptType","lepton type",-1,2);

  char selection_MW[900];
  sprintf( selection_MW, "mZjj>75. && mZjj<105. && nBTags==%d && leptType==%d && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f", nbtags, leptType_int, fitRangeLow, fitRangeHigh);


  RooFormulaVar massWindowSelection("massWindowSelection", selection_MW, RooArgList(CMS_hzz2l2q_mZZ,nBTags,mZjj,leptType));
  RooDataSet *dataset_obs_orig = new RooDataSet("dataset_obs_orig", "dataset_obs_orig", tree_data,
                                     RooArgSet(CMS_hzz2l2q_mZZ,nBTags,mZjj,leptType, eventWeight),
                                     selection_MW, "eventWeight");


  return double(dataset_obs_orig->numEntries());

}



double get_expectedYield_background( const std::string& dataset, const std::string& leptType_str, int nbtags, float mass, float gammaHiggs ) {

  int leptType_int = convert_leptType(leptType_str);

  char alphaFileName[300];
  sprintf(alphaFileName, "alphaFile_%s_%dbtag_ALL.root", dataset.c_str(), nbtags);
  TFile* alphaFile = TFile::Open(alphaFileName);
  TTree* treeSidebandsDATA_alphaCorr = (TTree*)alphaFile->Get("sidebandsDATA_alpha");
  TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, 150., 800.);
  char sidebandsCut_alpha[500];
  sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d && leptType==%d)", nbtags, leptType_int);
  treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "mZZ", sidebandsCut_alpha);
  float EvtNorm = h1_mZZ_sidebands_alpha->Integral();



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
      cbp.alpha = value;
    }
    if( varName=="wdth_rot" ) {
      cbp.wdth = value;
    }
    if( varName=="theta_best" ) {
      cbp.theta = value;
    }

  } // while ifs fitresults



  //integration window
  double effWidth = sqrt(gammaHiggs*gammaHiggs+100.);
  double fitRangeLow  = (mass-10.*effWidth<183.) ? 183. : mass-10.*effWidth;
  double fitRangeHigh = (mass+10.*effWidth>800.) ? 800. : mass+10.*effWidth;

  RooRealVar CMS_hzz2l2q_mZZ("CMS_hzz2l2q_mZZ", "zz inv mass",fitRangeLow,fitRangeHigh);
  RooRealVar CMS_hzz2l2q_mZZfull("CMS_hzz2l2q_mZZfull", "zz inv mass", 183., 800.);

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

  //RooGenericPdf fermiFull("fermiFull","fermi function","1/(1+exp((@1-@0)/@2))",RooArgList(CMS_hzz2l2q_mZZfull, fermi_cutoff, fermi_beta));
  RooFermi fermi_BKG_FULL("fermi_BKG_FULL","fermi function", CMS_hzz2l2q_mZZfull, fermi_cutoff, fermi_beta);
  RooCB CB_BKG_FULL("CBbkgFull","Crystal ball for background", CMS_hzz2l2q_mZZfull, m, wdth, alpha, n, theta);
  RooProdPdf backgroundFull("backgroundFull","backgroundFull",RooArgSet(fermi_BKG_FULL,CB_BKG_FULL));

  RooDataHist *BkgHisto = backgroundFull.generateBinned(CMS_hzz2l2q_mZZfull,EvtNorm,kTRUE,kFALSE);  
  //RooDataHist *BkgHisto = background.generateBinned(CMS_hzz2l2q_mZZ,EvtNorm,kTRUE,kFALSE);  

  char selection_massWindow[300];
  sprintf( selection_massWindow, "CMS_hzz2l2q_mZZfull>%f && CMS_hzz2l2q_mZZfull<%f", fitRangeLow, fitRangeHigh);
  //sprintf( selection_massWindow, "CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<%f", fitRangeLow, fitRangeHigh);
  double expectedYield_background = double( BkgHisto->sumEntries(selection_massWindow)  );

  return expectedYield_background;

}


//float CSggxs =1000.0*CSgg  *BRHZZ *BRZZ2l2q;
//float CSvbfxs=1000.0*CSvbf*BRHZZ  *BRZZ2l2q;
////   cout<<"THEOR xsect: "<<CSgg<<" , Vbf= "<<CSvbf<<"  BRHZZ="<<BRHZZ<<"  BRZZ2L2q="<<BRZZ2l2q<<endl;
//vector<float> myxsect;
//myxsect.push_back(CSggxs);
//myxsect.push_back(CSvbfxs);


int convert_leptType( const std::string& leptType ) {

  if( leptType!="ELE" && leptType!="MU" ) {
    std::cout << "WARNING!!! LeptType '" << leptType << "' is NOT supported!!! Returning -1." << std::endl;
    return -1;
  }

  int leptType_int = (leptType=="MU" ) ? 0 : 1;

  return leptType_int;

}
