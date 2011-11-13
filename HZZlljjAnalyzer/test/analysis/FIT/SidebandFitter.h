#ifndef SidebandFitter_h
#define SidebandFitter_h

#include <string>
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"

#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"





class SidebandFitter {

 public:

  SidebandFitter( const std::string& dataset, const std::string PUType );

  ~SidebandFitter() { 
    if( workspace_!=0 ) delete workspace_;
  }


  TH1D* getAlphaHisto( int btagCategory, const std::string leptType_str, TTree* treeMC );
  TH1D* getAlphaHistoSmooth( int btagCategory, const std::string leptType_str, TTree* treeMC );
  
  //RooFitResult* fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed=-1 );
  
  RooFitResult* fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed=-1 , std::string init="MCSignal");
  void generateFixedPars(TTree* treeMC,int btagCategory, const std::string& leptType, TH1D* h1_alpha);
  void fitPseudo( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed , std::string init="MCSignal");

  std::string get_fitResultsName( int nbtags, const std::string& data_mc="DATA" );

  std::string get_outdir();
 
  TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name );

  TH1D* shuffle( TH1D* inhist, TRandom3* random, char *histName );

  void pseudoMassge(int ntoys, int btagCategory , const std::string& leptType, std::string init, RooFitResult* r_nominal);

  // this method returns only rate:
  Double_t get_backgroundNormalization( int nbtags, const std::string& leptType, const std::string& data_mc="DATA", float mZZmin=-1., float mZZmax=-1. );
  // this one returns both rate (first) and error on rate (second):
  std::pair<Double_t, Double_t> get_backgroundNormalizationAndError( int nbtags, const std::string& leptType, const std::string& data_mc="DATA", float mZZmin=-1., float mZZmax=-1.);

  RooDataSet* get_observedDataset( RooRealVar* CMS_hzz2l2q_mZZ, const std::string& leptType_str, int nbtags );

  static int convert_leptType( const std::string& leptType );

  RooPlot* ContourPlot(std::string var1, std::string var2, RooFitResult* r);

 private:

  std::string dataset_;
  std::string PUType_;

  float mZZmin_;
  float mZZmax_;

  RooWorkspace* workspace_;

};


#endif
