#ifndef SidebandFitter_h
#define SidebandFitter_h

#include <string>
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"


struct FitResults {

  float fermi_beta;
  float fermi_cutoff;
  float CB_m;
  float CB_wdth;
  float CB_alpha;
  float CB_n;
  float CB_theta;

  float fermi_beta_err;
  float fermi_cutoff_err;
  float CB_m_err;
  float CB_wdth_err;
  float CB_alpha_err;
  float CB_n_err;
  float CB_theta_err;

  float CB_alpha_rot;
  float CB_wdth_rot;
  float CB_theta_best;

};




class SidebandFitter {

 public:

  SidebandFitter( const std::string& dataset );
  ~SidebandFitter();

  TH1D* getAlphaHisto( int btagCategory, const std::string leptType_str, TTree* treeMC );
  
  FitResults fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed=-1 );

  std::string get_fitResultsName( int nbtags, const std::string& data_mc="DATA" );

  std::string get_outdir();
 
  TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name );

  TH1D* shuffle( TH1D* inhist, TRandom3* random, char *histName );

  void modifyFitResultError( const std::string& thisVar, double thisVarError, int nbtags );


 private:

  std::string dataset_;

};


#endif
