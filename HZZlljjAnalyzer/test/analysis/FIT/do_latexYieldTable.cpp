#include "RooDataSet.h"
#include "RooRealVar.h"

#include "SidebandFitter.h"


//horrible but dont have time now:
double LUMI_ = 4600.;



struct EventYields {

  int observed;
  Double_t expectedDATA;
  Double_t expectedMC;
  Double_t expectedDATA_error;
  Double_t expectedMC_error;

};



EventYields getYields( const std::string& dataset, const std::string& PUType, int nbtags );



int main( int argc, char* argv[] ) {


  if( argc != 2 ) {
    std::cout << "USAGE: ./do_latexYieldTable [dataset]" << std::endl;
    exit(11);
  }


  std::string dataset(argv[1]);

  std::string PUType = "Run2011A";
  if( dataset=="HR11" )
    PUType = "HR11";
  if( dataset=="HR11_v2" )
    PUType = "HR11_73pb";

  
  EventYields yield_0tag = getYields( dataset, PUType, 0 );
  EventYields yield_1tag = getYields( dataset, PUType, 1 );
  EventYields yield_2tag = getYields( dataset, PUType, 2 );

  std::cout << "\\multicolumn{2}{|c|}{observed yield}     & " << yield_0tag.observed << "  & " << yield_1tag.observed << " & " << yield_2tag.observed << "\\\\" << std::endl;

  std::cout << "\\multicolumn{2}{|c|}{exp. background (data)}   & $" << yield_0tag.expectedDATA << "\\pm" << yield_0tag.expectedDATA_error;
  std::cout << "$  & $" << yield_1tag.expectedDATA << "\\pm" << yield_1tag.expectedDATA_error;
  std::cout << "$  & $" << yield_2tag.expectedDATA << "\\pm" << yield_2tag.expectedDATA_error;
  std::cout << "$  \\\\" << std::endl;

  std::cout << "\\multicolumn{2}{|c|}{exp. background (mc)}   & $" << yield_0tag.expectedMC << "\\pm" << yield_0tag.expectedMC_error;
  std::cout << "$  & $" << yield_1tag.expectedMC << "\\pm" << yield_1tag.expectedMC_error;
  std::cout << "$  & $" << yield_2tag.expectedMC << "\\pm" << yield_2tag.expectedMC_error;
  std::cout << "$  \\\\" << std::endl;

  return 0;

}



EventYields getYields( const std::string& dataset, const std::string& PUType, int nbtags ) {

  EventYields returnYields;
  
  SidebandFitter* sf = new SidebandFitter( dataset, PUType );


  std::pair<Double_t,Double_t> rate_backgroundDATA = sf->get_backgroundNormalizationAndError( nbtags, "ALL", "DATA" );
  returnYields.expectedDATA = rate_backgroundDATA.first;
  returnYields.expectedDATA_error = rate_backgroundDATA.second;

  std::pair<Double_t,Double_t> rate_backgroundMC = sf->get_backgroundNormalizationAndError( nbtags, "ALL", "MC" );
  returnYields.expectedMC = rate_backgroundMC.first;
  returnYields.expectedMC_error = rate_backgroundMC.second;

  returnYields.expectedMC *= LUMI_;
  returnYields.expectedMC_error *= LUMI_;


  RooRealVar* CMS_hzz2l2q_mZZ = new RooRealVar("CMS_hzz2l2q_mZZ", "m_{lljj}", 160., 800., "GeV");
  RooDataSet* observed = sf->get_observedDataset( CMS_hzz2l2q_mZZ, "ALL", nbtags );
  returnYields.observed = observed->sumEntries();

  return returnYields;

}

