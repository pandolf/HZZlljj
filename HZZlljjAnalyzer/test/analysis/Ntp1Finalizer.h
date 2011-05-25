// ------------------------------------------------------------
//  
//    Ntp1Finalizer - Base class for finalizing the analysis
//    Used in the Ntp1 tree workflow, after a Ntp1Analyzer
//    class has produced 2nd level trees 
//    (and merge_and_setWeights has done its job)
//
//    It is an abstract class: the method 'finalize' has to be
//    implemented.
//
// ------------------------------------------------------------


#include <vector>
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "AnalysisJet.h"
#include "BTagSFUtil/BTagSFUtil.h"



class Ntp1Finalizer {


 public:

  Ntp1Finalizer( const std::string& analyzerType, const std::string& dataset, const std::string& flags="" );
  virtual ~Ntp1Finalizer();

  void createOutputFile( const std::string& additionalFlags="" );
  virtual void addFile(const std::string& dataset);

  std::vector<TH1F*> getResponseHistos(const std::string& name, unsigned binArraySize, Double_t* ptBins);
  void writeResponseHistos( TFile* file, std::vector<TH1F*> h1_response, std::string dirName );

  TChain* get_tree() { return tree_; };
  TFile* get_outFile() { return outFile_; };
  bool get_DEBUG() { return DEBUG_; };
  int get_nBTags( const AnalysisJet& jet1, const AnalysisJet& jet2, BTagSFUtil* btsfutil, bool losebtags=true );

  void set_outFile( const std::string& fileName="", const std::string& suffix="" );
  void set_inputAnalyzerType( const std::string& analyzerType ) { inputAnalyzerType_ = analyzerType; };
  void set_flags( const std::string& flags ) { flags_ = flags; };
  void set_DEBUG( bool DEBUG ) { DEBUG_ = DEBUG; };

  virtual void finalize() = 0;


  TChain* tree_;

  float nCounter_;
  float nCounterW_;

  std::string analyzerType_;
  std::string inputAnalyzerType_;
  std::string dataset_;
  std::string flags_;

  TFile* outFile_;

  bool DEBUG_;

 private:

};
