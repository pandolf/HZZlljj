#ifndef SidebandFitter_h
#define SidebandFitter_h

#include <string>
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"

#include "RooFitResult.h"





class SidebandFitter {

 public:

  SidebandFitter( const std::string& dataset, const std::string PUType );
  ~SidebandFitter();

  TH1D* getAlphaHisto( int btagCategory, const std::string leptType_str, TTree* treeMC );
  
  RooFitResult* fitSidebands( TTree* treeMC, TTree* treeDATA, int btagCategory, const std::string& leptType, TH1D* h1_alpha, int seed=-1 );

  std::string get_fitResultsName( int nbtags, const std::string& data_mc="DATA" );

  std::string get_outdir();
 
  TTree* correctTreeWithAlpha( TTree* tree, TH1D* h1_alpha, int btagCategory, const std::string& name );

  TH1D* shuffle( TH1D* inhist, TRandom3* random, char *histName );

  static float get_backgroundNormalization( const std::string& data_dataset, const std::string& PUType, int nbtags, const std::string& leptType );

  static int convert_leptType( const std::string& leptType );


 private:

  std::string dataset_;
  std::string PUType_;

  float mZZmin_;
  float mZZmax_;

};


#endif
