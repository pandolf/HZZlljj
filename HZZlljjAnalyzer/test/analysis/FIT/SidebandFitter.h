#ifndef SidebandFitter_h
#define SidebandFitter_h

#include <string>
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"

#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooCBShape.h"

#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"




class SidebandFitter {

 public:

  SidebandFitter( const std::string& dataset, const std::string& PUType, const std::string& init );
  ~SidebandFitter() {};

  TH1D* getAlphaHisto( int btagCategory, const std::string& leptType_str, TTree* treeMC );
  
  RooFitResult* fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha );

  std::string get_fitResultsName( int nbtags, const std::string& init );
  std::string get_fitResultsRootFileName( int btagCategory, const std::string& leptType );

  std::string get_outdir();
 
  TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name );

  TH1D* shuffle( TH1D* inhist, TRandom3* random, char *histName );
  RooPlot* ContourPlot(RooRealVar* var1,RooRealVar* var2, RooFitResult* r);

  // this method return only rate:
  Double_t get_backgroundNormalization( int nbtags, const std::string& leptType, const std::string& data_mc, float mZZmin=-1., float mZZmax=-1. );
  // this one return both rate (first) and error on rate (second):
  std::pair<Double_t, Double_t> get_backgroundNormalizationAndError( int nbtags, const std::string& leptType, const std::string& data_mc, float mZZmin=-1., float mZZmax=-1.);

  RooDataSet* get_observedDataset( RooRealVar* CMS_hzz2l2q_mZZ, const std::string& leptType_str, int nbtags );

  static int convert_leptType( const std::string& leptType );

  void fitPseudo( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed );
  void pseudoMassge(int ntoys, int btagCategory , const std::string& leptType, RooFitResult* r_nominal);


 private:

  std::string dataset_;
  std::string PUType_;
  std::string init_;

  float mZZmin_;
  float mZZmax_;

  RooRealVar* CMS_hzz2l2q_mZZ_;
  
  RooRealVar* beta_;
  RooRealVar* cutOff_;
  RooRealVar* m_;
  RooRealVar* n_;
  RooRealVar* wdth_;
  RooRealVar* alpha_;

  RooFermi* fermi_;
  RooCBShape* CBShape_;
 

  RooProdPdf* background_;
   
  RooAbsPdf* background_decorr_;

  RooWorkspace* fitWorkspace;
};


#endif
